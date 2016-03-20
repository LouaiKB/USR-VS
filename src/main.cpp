#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <vector>
#include <array>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <thread>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <mongo/client/dbclient.h>
#include "io_service_pool.hpp"
#include "safe_counter.hpp"
using namespace std;
using namespace std::chrono;
using namespace RDKit;
using namespace RDGeom;
using namespace boost::filesystem;
using namespace boost::gregorian;
using namespace boost::posix_time;
using namespace mongo;
using namespace bson;

inline static auto local_time()
{
	return to_simple_string(microsec_clock::local_time()) + " ";
}

inline static auto milliseconds_since_epoch()
{
	return Date_t(duration_cast<chrono::milliseconds>(system_clock::now().time_since_epoch()).count());
}

template <typename T>
inline vector<T> read(const path src)
{
	boost::filesystem::ifstream ifs(src, ios::binary | ios::ate);
	const size_t num_bytes = ifs.tellg();
	cout << local_time() << "Reading " << src << " of " << num_bytes << " bytes" << endl;
	vector<T> buf;
	buf.resize(num_bytes / sizeof(T));
	ifs.seekg(0);
	ifs.read(reinterpret_cast<char*>(buf.data()), num_bytes);
	return buf;
}

template <typename size_type>
class header_array
{
public:
	explicit header_array(path src)
	{
		src.replace_extension(".ftr");
		boost::filesystem::ifstream ifs(src, ios::binary | ios::ate);
		const size_t num_bytes = ifs.tellg();
		cout << local_time() << "Reading " << src << " of " << num_bytes << " bytes" << endl;
		hdr.resize(1 + num_bytes / sizeof(size_type));
		hdr.front() = 0;
		ifs.seekg(0);
		ifs.read(reinterpret_cast<char*>(hdr.data() + 1), num_bytes);
	}

	size_t size() const
	{
		return hdr.size() - 1;
	}

protected:
	vector<size_type> hdr;
};

template <typename size_type>
class string_array : public header_array<size_type>
{
public:
	explicit string_array(const path src) : header_array<size_type>(src)
	{
		boost::filesystem::ifstream ifs(src, ios::binary | ios::ate);
		const size_t num_bytes = ifs.tellg();
		cout << local_time() << "Reading " << src << " of " << num_bytes << " bytes" << endl;
		buf.resize(num_bytes);
		ifs.seekg(0);
		ifs.read(const_cast<char*>(buf.data()), num_bytes);
	}

	string operator[](const size_t index) const
	{
		const auto pos = this->hdr[index];
		const auto len = this->hdr[index + 1] - pos;
		return buf.substr(pos, len);
	}

protected:
	string buf;
};

template <typename size_type>
class stream_array : public header_array<size_type>
{
public:
	explicit stream_array(const path src) : header_array<size_type>(src), ifs(src, ios::binary)
	{
	}

	string operator[](const size_t index)
	{
		const auto pos = this->hdr[index];
		const auto len = this->hdr[index + 1] - pos;
		string buf;
		buf.resize(len);
		ifs.seekg(pos);
		ifs.read(const_cast<char*>(buf.data()), len);
		return buf;
	}

protected:
	boost::filesystem::ifstream ifs;
};

template<typename T>
double dist2(const T& p0, const T& p1)
{
	const auto d0 = p0[0] - p1[0];
	const auto d1 = p0[1] - p1[1];
	const auto d2 = p0[2] - p1[2];
	return d0 * d0 + d1 * d1 + d2 * d2;
}

array<Point3D, 4> calcRefPoints(const ROMol& mol, const vector<int>& heavyAtoms)
{
	const auto num_points = heavyAtoms.size();
	assert(num_points == mol.getNumHeavyAtoms());
	const auto& conf = mol.getConformer();
	array<Point3D, 4> refPoints;
	for (auto& ref : refPoints)
	{
		assert(ref[0] == 0);
		assert(ref[1] == 0);
		assert(ref[2] == 0);
	}
	auto& ctd = refPoints[0];
	auto& cst = refPoints[1];
	auto& fct = refPoints[2];
	auto& ftf = refPoints[3];
	for (const auto i : heavyAtoms)
	{
		const auto& a = conf.getAtomPos(i);
		ctd += a;
	}
	ctd /= num_points;
	double cst_dist = numeric_limits<double>::max();
	double fct_dist = numeric_limits<double>::lowest();
	double ftf_dist = numeric_limits<double>::lowest();
	for (const auto i : heavyAtoms)
	{
		const auto& a = conf.getAtomPos(i);
		const auto this_dist = dist2(a, ctd);
		if (this_dist < cst_dist)
		{
			cst = a;
			cst_dist = this_dist;
		}
		if (this_dist > fct_dist)
		{
			fct = a;
			fct_dist = this_dist;
		}
	}
	for (const auto i : heavyAtoms)
	{
		const auto& a = conf.getAtomPos(i);
		const auto this_dist = dist2(a, fct);
		if (this_dist > ftf_dist)
		{
			ftf = a;
			ftf_dist = this_dist;
		}
	}
	return refPoints;
}

int main(int argc, char* argv[])
{
	// Check the required number of command line arguments.
	if (argc != 5)
	{
		cout << "usr host user pwd jobs_path" << endl;
		return 0;
	}

	// Fetch command line arguments.
	const auto host = argv[1];
	const auto user = argv[2];
	const auto pwd = argv[3];
	const path jobs_path = argv[4];

	DBClientConnection conn;
	{
		// Connect to host and authenticate user.
		cout << local_time() << "Connecting to " << host << " and authenticating " << user << endl;
		string errmsg;
		if ((!conn.connect(host, errmsg)) || (!conn.auth("istar", user, pwd, errmsg)))
		{
			cerr << local_time() << errmsg << endl;
			return 1;
		}
	}

	// Initialize constants.
	cout << local_time() << "Initializing" << endl;
	const auto collection = "istar.usr2";
	const size_t num_usrs = 2;
	const array<string, 2> usr_names{{ "USR", "USRCAT" }};
	constexpr array<size_t, num_usrs> qn{{ 12, 60 }};
	constexpr array<double, num_usrs> qv{{ 1.0 / qn[0], 1.0 / qn[1] }};
	const size_t num_refPoints = 4;
	const size_t num_subsets = 5;
	const array<string, num_subsets> SubsetSMARTS
	{{
		"[!#1]", // heavy
		"[#6+0!$(*~[#7,#8,F]),SH0+0v2,s+0,S^3,Cl+0,Br+0,I+0]", // hydrophobic
		"[a]", // aromatic
		"[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N&v3;H1,H2]-[!$(*=[O,N,P,S])]),$([N;v3;H0]),$([n,o,s;+0]),F]", // acceptor
		"[N!H0v3,N!H0+v4,OH+0,SH+0,nH+0]", // donor
	}};
	const size_t num_hits = 100;

	// Wrap SMARTS strings to RWMol objects.
	array<unique_ptr<ROMol>, num_subsets> SubsetMols;
	for (size_t k = 0; k < num_subsets; ++k)
	{
		SubsetMols[k].reset(reinterpret_cast<ROMol*>(SmartsToMol(SubsetSMARTS[k])));
	}

	// Read ZINC ID file.
	const string_array<size_t> zincids("16/zincid.txt");
	const auto num_ligands = zincids.size();
	cout << local_time() << "Found " << num_ligands << " database molecules" << endl;

	// Read SMILES file.
	const string_array<size_t> smileses("16/smiles.txt");
	assert(smileses.size() == num_ligands);

	// Read supplier file.
	const string_array<size_t> suppliers("16/supplier.txt");
	assert(suppliers.size() == num_ligands);

	// Read property files of floating point types and integer types.
	const auto zfproperties = read<array<float, 4>>("16/zfprop.f32");
	assert(zfproperties.size() == num_ligands);
	const auto ziproperties = read<array<int16_t, 5>>("16/ziprop.i16");
	assert(ziproperties.size() == num_ligands);

	// Read cumulative number of conformers file.
	const auto mconfss = read<size_t>("16/mconfs.u64");
	const auto num_conformers = mconfss.back();
	assert(mconfss.size() == num_ligands);
	assert(num_conformers >= num_ligands);
	cout << local_time() << "Found " << num_conformers << " database conformers" << endl;

	// Read feature file.
	const auto features = read<array<double, qn.back()>>("16/usrcat.f64");
	assert(features.size() == num_conformers);

	// Read ligand footer file and open ligand SDF file for seeking and reading.
	stream_array<size_t> ligands("16/ligand.sdf");
	assert(ligands.size() == num_conformers);

	// Initialize variables.
	array<vector<int>, num_subsets> subsets;
	array<Point3D, num_refPoints> refPoints;
	array<vector<double>, num_refPoints> dista;
	alignas(32) array<double, qn.back()> q;

	// Initialize vectors to store compounds' primary score and their corresponding conformer.
	vector<double> scores(num_ligands); // Primary score of molecules.
	vector<size_t> cnfids(num_ligands); // ID of conformer with the best primary score.
	const auto compare = [&](const size_t val0, const size_t val1) // Sort by the primary score.
	{
		return scores[val0] < scores[val1];
	};

	// Initialize an io service pool and create worker threads for later use.
	const size_t num_threads = thread::hardware_concurrency();
	cout << local_time() << "Creating an io service pool of " << num_threads << " worker threads" << endl;
	io_service_pool io(num_threads);
	safe_counter<size_t> cnt;

	// Initialize the number of chunks and the number of molecules per chunk.
	const auto num_chunks = num_threads << 4;
	const auto chunk_size = 1 + (num_ligands - 1) / num_chunks;
	assert(chunk_size * num_chunks >= num_ligands);
	assert(chunk_size >= num_hits);
	cout << local_time() << "Using " << num_chunks << " chunks and a chunk size of " << chunk_size << endl;
	vector<size_t> scase(num_ligands);
	vector<size_t> zcase(num_hits * (num_chunks - 1) + min(num_hits, num_ligands - chunk_size * (num_chunks - 1))); // The last chunk might have fewer than num_hits records.

	// Enter event loop.
	cout << local_time() << "Entering event loop" << endl;
	cout.setf(ios::fixed, ios::floatfield);
	bool sleeping = false;
	while (true)
	{
		// Fetch an incompleted job in a first-come-first-served manner.
		if (!sleeping) cout << local_time() << "Fetching an incompleted job" << endl;
		BSONObj info;
		const auto started = milliseconds_since_epoch();
		conn.runCommand("istar", BSON("findandmodify" << "usr2" << "query" << BSON("started" << BSON("$exists" << false)) << "sort" << BSON("submitted" << 1) << "update" << BSON("$set" << BSON("started" << started))), info); // conn.findAndModify() is available since MongoDB C++ Driver legacy-1.0.0
		const auto value = info["value"];
		if (value.isNull())
		{
			// No incompleted jobs. Sleep for a while.
			if (!sleeping) cout << local_time() << "Sleeping" << endl;
			sleeping = true;
			this_thread::sleep_for(chrono::seconds(2));
			continue;
		}
		sleeping = false;
		const auto job = value.Obj();

		// Obtain job properties.
		const auto _id = job["_id"].OID();
		cout << local_time() << "Executing job " << _id.str() << endl;
		const auto job_path = jobs_path / _id.str();
		const size_t usr0 = job["usr"].Int(); // Specify the primary sorting score. 0: USR; 1: USRCAT.
		assert(usr0 == 0 || usr0 == 1);
		const auto usr1 = usr0 ^ 1;
		const auto qnu0 = qn[usr0];
		const auto qnu1 = qn[usr1];

		// Read and validate the user-supplied SDF file.
		cout << local_time() << "Reading and validating the query file" << endl;
		SDMolSupplier sup((job_path / "query.sdf").string(), true, false, true); // sanitize, removeHs, strictParsing
		if (!sup.length() || !sup.atEnd())
		{
			const auto error = 1;
			cout << local_time() << "Failed to parse the query file, error code = " << error << endl;
			conn.update(collection, BSON("_id" << _id), BSON("$set" << BSON("completed" << milliseconds_since_epoch() << "error" << error)));
			continue;
		}

		// Process each of the query molecules sequentially.
		const auto num_queries = 1; // Restrict the number of query molecules to 1. Setting num_queries = sup.length() to execute any number of query molecules.
		for (unsigned int query_number = 0; query_number < num_queries; ++query_number)
		{
			cout << local_time() << "Processing query molecule " << query_number << endl;
			const unique_ptr<ROMol> qry_ptr(sup.next()); // Calling next() may print "ERROR: Could not sanitize molecule on line XXXX" to stderr.
			auto& qryMol = *qry_ptr;

			// Get the number of points, excluding hydrogens.
			const auto num_points = qryMol.getNumHeavyAtoms();
			assert(num_points);
			cout << local_time() << "Found " << num_points << " heavy atoms" << endl;

			// Classify atoms to pharmacophoric subsets.
			cout << local_time() << "Classifying atoms into subsets" << endl;
			for (size_t k = 0; k < num_subsets; ++k)
			{
				vector<vector<pair<int, int>>> matchVect;
				SubstructMatch(qryMol, *SubsetMols[k], matchVect);
				const auto num_matches = matchVect.size();
				auto& subset = subsets[k];
				subset.resize(num_matches);
				for (size_t i = 0; i < num_matches; ++i)
				{
					subset[i] = matchVect[i].front().second;
				}
				cout << local_time() << "Found " << num_matches << " atoms for subset " << k << endl;
			}
			const auto& subset0 = subsets.front();
			assert(subset0.size() == num_points);

			// Calculate the four reference points.
			cout << local_time() << "Calculating " << num_refPoints << " reference points" << endl;
			const auto qryRefPoints = calcRefPoints(qryMol, subset0);

			// Precalculate the distances of heavy atoms to the reference points, given that subsets[1 to 4] are subsets of subsets[0].
			cout << local_time() << "Calculating " << num_points * num_refPoints << " pairwise distances" << endl;
			const auto& qryCnf = qryMol.getConformer();
			for (size_t k = 0; k < num_refPoints; ++k)
			{
				const auto& refPoint = refPoints[k];
				auto& distp = dista[k];
				distp.resize(num_points);
				for (size_t i = 0; i < num_points; ++i)
				{
					distp[subset0[i]] = sqrt(dist2(qryCnf.getAtomPos(subset0[i]), refPoint));
				}
			}

			// Loop over pharmacophoric subsets and reference points.
			cout << local_time() << "Calculating " << 3 * num_refPoints * num_subsets << " moments of USRCAT feature" << endl;
			size_t qo = 0;
			for (const auto& subset : subsets)
			{
				const auto n = subset.size();
				for (size_t k = 0; k < num_refPoints; ++k)
				{
					// Load distances from precalculated ones.
					const auto& distp = dista[k];
					vector<double> dists(n);
					for (size_t i = 0; i < n; ++i)
					{
						dists[i] = distp[subset[i]];
					}

					// Compute moments.
					array<double, 3> m{};
					if (n > 2)
					{
						const auto v = 1.0 / n;
						for (size_t i = 0; i < n; ++i)
						{
							const auto d = dists[i];
							m[0] += d;
						}
						m[0] *= v;
						for (size_t i = 0; i < n; ++i)
						{
							const auto d = dists[i] - m[0];
							m[1] += d * d;
						}
						m[1] = sqrt(m[1] * v);
						for (size_t i = 0; i < n; ++i)
						{
							const auto d = dists[i] - m[0];
							m[2] += d * d * d;
						}
						m[2] = cbrt(m[2] * v);
					}
					else if (n == 2)
					{
						m[0] = 0.5 *     (dists[0] + dists[1]);
						m[1] = 0.5 * fabs(dists[0] - dists[1]);
					}
					else if (n == 1)
					{
						m[0] = dists[0];
					}
					for (const auto e : m)
					{
						q[qo++] = e;
					}
				}
			}
			assert(qo == qn.back());

			// Compute USR and USRCAT scores.
			cout << local_time() << "Calculating " << num_ligands << " " << usr_names[usr0] << " scores" << endl;
			scores.assign(scores.size(), numeric_limits<double>::max());
			iota(scase.begin(), scase.end(), 0);
			cnt.init(num_chunks);
			for (size_t l = 0; l < num_chunks; ++l)
			{
				io.post([&,l]()
				{
					// Loop over molecules of the current chunk.
					const auto chunk_beg = chunk_size * l;
					const auto chunk_end = min(chunk_beg + chunk_size, num_ligands);
					for (size_t k = chunk_beg; k < chunk_end; ++k)
					{
						// Loop over conformers of the current molecule and calculate their primary score.
						auto& scorek = scores[k];
						size_t j = k ? mconfss[k - 1] : 0;
						for (const auto mconfs = mconfss[k]; j < mconfs; ++j)
						{
							const auto& d = features[j];
							double s = 0;
							for (size_t i = 0; i < qnu0; ++i)
							{
								s += abs(q[i] - d[i]);
								if (s >= scorek) break;
							}
							if (s < scorek)
							{
								scorek = s;
								cnfids[k] = j;
							}
						}
					}

					// Sort the scores of molecules of the current chunk.
					sort(scase.begin() + chunk_beg, scase.begin() + chunk_end, compare);

					// Copy the indexes of top hits of the current chunk to a global vector for final sorting.
					copy_n(scase.begin() + chunk_beg, min(num_hits, chunk_end - chunk_beg), zcase.begin() + num_hits * l);

					cnt.increment();
				});
			}
			cnt.wait();

			// Sort the top hits from chunks.
			cout << local_time() << "Sorting " << zcase.size() << " hits by " << usr_names[usr0] << " score" << endl;
			sort(zcase.begin(), zcase.end(), compare);

			// Create output directory and write output files.
			cout << local_time() << "Creating output directory and writing output files" << endl;
			const auto output_dir = job_path / to_string(query_number);
			create_directory(output_dir);
			boost::filesystem::ofstream hits_sdf(output_dir / "hits.sdf");
			boost::filesystem::ofstream hits_csv(output_dir / "hits.csv");
			hits_csv.setf(ios::fixed, ios::floatfield);
			hits_csv << "ZINC ID,USR score,USRCAT score,Molecular weight (g/mol),Partition coefficient xlogP,Apolar desolvation (kcal/mol),Polar desolvation (kcal/mol),Hydrogen bond donors,Hydrogen bond acceptors,Polar surface area tPSA (Å^2),Net charge,Rotatable bonds,SMILES,Suppliers and annotations\n";
			for (size_t l = 0; l < num_hits; ++l)
			{
				const auto k = zcase[l];

				// Calculate the secondary score of the saved conformer, which has the best primary score.
				const auto& d = features[cnfids[k]];
				double s = 0;
				for (size_t i = 0; i < qnu1; ++i)
				{
					s += abs(q[i] - d[i]);
				}

				const auto u0score = 1 / (1 + scores[k] * qv[usr0]); // Primary score of the current molecule.
				const auto u1score = 1 / (1 + s         * qv[usr1]); // Secondary score of the current molecule.
				const auto zincid = zincids[k].substr(0, 8); // Take another substr() to get rid of the trailing newline.
				const auto zfp = zfproperties[k];
				const auto zip = ziproperties[k];
				const auto smiles = smileses[k];    // A newline is already included in smileses[k].
				const auto supplier = suppliers[k]; // A newline is already included in suppliers[k].
				hits_csv
					<< zincid
					<< setprecision(8)
					<< ',' << (usr1 ? u0score : u1score)
					<< ',' << (usr1 ? u1score : u0score)
					<< setprecision(3)
					<< ',' << zfp[0]
					<< ',' << zfp[1]
					<< ',' << zfp[2]
					<< ',' << zfp[3]
					<< ',' << zip[0]
					<< ',' << zip[1]
					<< ',' << zip[2]
					<< ',' << zip[3]
					<< ',' << zip[4]
					<< ',' << smiles.substr(0, smiles.length() - 1)     // Get rid of the trailing newline.
					<< ',' << supplier.substr(0, supplier.length() - 1) // Get rid of the trailing newline.
					<< '\n'
				;
				const auto lig = ligands[cnfids[k]];
				hits_sdf.write(lig.data(), lig.size());
			}
		}

		// Update job status.
		cout << local_time() << "Setting completed time" << endl;
		const auto completed = milliseconds_since_epoch();
		conn.update(collection, BSON("_id" << _id), BSON("$set" << BSON("completed" << completed << "nqueries" << num_queries)));

		// Calculate runtime in seconds and screening speed in million conformers per second.
		const auto runtime = (completed - started) * 0.001;
		const auto speed = num_conformers * 0.000001 * num_queries / runtime;
		cout
			<< local_time() << "Completed " << num_queries << " " << (num_queries == 1 ? "query" : "queries") << " in " << setprecision(3) << runtime << " seconds" << endl
			<< local_time() << "Screening speed was " << setprecision(0) << speed << " M conformers per second" << endl
		;
	}
}
