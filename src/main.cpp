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
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
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

inline static string local_time()
{
	return to_simple_string(microsec_clock::local_time()) + " ";
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
	constexpr array<size_t, num_usrs> qn{{ 12, 60 }};
	constexpr array<double, num_usrs> qv{{ 1.0 / qn[0], 1.0 / qn[1] }};
	const size_t num_references = 4;
	const size_t num_subsets = 5;
	const array<string, num_subsets> SubsetSMARTS
	{{
		"[!#1]", // heavy
		"[#6+0!$(*~[#7,#8,F]),SH0+0v2,s+0,S^3,Cl+0,Br+0,I+0]", // hydrophobic
		"[a]", // aromatic
		"[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N&v3;H1,H2]-[!$(*=[O,N,P,S])]),$([N;v3;H0]),$([n,o,s;+0]),F]", // acceptor
		"[N!H0v3,N!H0+v4,OH+0,SH+0,nH+0]", // donor
	}};

	// Wrap SMARTS strings to RWMol objects.
	array<unique_ptr<ROMol>, num_subsets> SubsetMols;
	for (size_t k = 0; k < num_subsets; ++k)
	{
		SubsetMols[k].reset(reinterpret_cast<ROMol*>(SmartsToMol(SubsetSMARTS[k])));
	}

	// Read ZINC ID file.
	const string_array<size_t> zincids("16_zincid.txt");
	const auto num_ligands = zincids.size();

	// Read SMILES file.
	const string_array<size_t> smileses("16_smiles.txt");
	assert(smileses.size() == num_ligands);

	// Read supplier file.
	const string_array<size_t> suppliers("16_supplier.txt");
	assert(suppliers.size() == num_ligands);

	// Read property files of floating point types and integer types.
	const auto zfproperties = read<array<float, 4>>("16_zfprop.f32");
	assert(zfproperties.size() == num_ligands);
	const auto ziproperties = read<array<int16_t, 5>>("16_ziprop.i16");
	assert(ziproperties.size() == num_ligands);

	// Read cumulative number of conformers file.
	const auto mconfss = read<size_t>("16_mconfs.u64");
	const auto num_conformers = mconfss.back();
	assert(mconfss.size() == num_ligands);
	assert(num_conformers >= num_ligands);

	// Read feature file.
#ifdef PRELOAD_FEATURES
	const auto features = read<array<double, qn.back()>>("16_usrcat.f64");
	assert(features.size() == num_conformers);
#endif

	// Read ligand footer file and open ligand SDF file for seeking and reading.
	stream_array<size_t> ligands("16_ligand.sdf");
	assert(ligands.size() == num_conformers);

	// Initialize variables.
	array<vector<int>, num_subsets> subsets;
	array<Point3D, num_references> references;
	array<vector<double>, num_references> dista;
	alignas(32) array<double, qn.back()> q;

	// Initialize vectors to store compounds' USR and USRCAT scores and their corresponding conformer.
	array<vector<double>, 2> scores
	{{
		vector<double>(num_ligands),
		vector<double>(num_ligands)
	}};
	array<vector<size_t>, 2> cnfids
	{{
		vector<size_t>(num_ligands),
		vector<size_t>(num_ligands)
	}};
	vector<size_t> scase(num_ligands);

	// Initialize an io service pool and create worker threads for later use.
	const size_t num_threads = thread::hardware_concurrency();
	cout << local_time() << "Creating an io service pool of " << num_threads << " worker threads" << endl;
	io_service_pool io(num_threads);
	safe_counter<size_t> cnt;

	// Initialize the number of chunks and the number of molecules per chunk.
	const auto num_chunks = num_threads << 6;
	const auto chunk_size = 1 + (num_ligands - 1) / num_chunks;

	// Enter event loop.
	cout << local_time() << "Entering event loop" << endl;
	bool sleeping = false;
	while (true)
	{
		// Fetch an incompleted job in a first-come-first-served manner.
		if (!sleeping) cout << local_time() << "Fetching an incompleted job" << endl;
		BSONObj info;
		conn.runCommand("istar", BSON("findandmodify" << "usr2" << "query" << BSON("done" << BSON("$exists" << false) << "started" << BSON("$exists" << false)) << "sort" << BSON("submitted" << 1) << "update" << BSON("$set" << BSON("started" << Date_t(duration_cast<std::chrono::milliseconds>(system_clock::now().time_since_epoch()).count())))), info); // conn.findAndModify() is available since MongoDB C++ Driver legacy-1.0.0
		const auto value = info["value"];
		if (value.isNull())
		{
			// No incompleted jobs. Sleep for a while.
			if (!sleeping) cout << local_time() << "Sleeping" << endl;
			sleeping = true;
			this_thread::sleep_for(chrono::seconds(10));
			continue;
		}
		sleeping = false;
		const auto job = value.Obj();

		// Obtain job properties.
		const auto _id = job["_id"].OID();
		cout << local_time() << "Executing job " << _id.str() << endl;
		const auto job_path = jobs_path / _id.str();
//		const auto usr = job["usr"].Int();
		const auto usr = 1; // Specify the primary sorting score. 0: USR; 1: USRCAT.
		const auto& u0scores = scores[usr];   // Primary sorting score.
		const auto& u1scores = scores[usr^1]; // Secondary sorting score.

		// Parse the user-supplied SDF file, setting sanitize=true, removeHs=false, strictParsing=true.
		size_t query_number = 0;
		for (SDMolSupplier sup((job_path / "query.sdf").string(), true, false, true); !sup.atEnd();) // sanitize, removeHs, strictParsing
		{
			cout << local_time() << "Parsing query molecule " << ++query_number << endl;
			const unique_ptr<ROMol> mol_ptr(sup.next()); // Calling next() may print "ERROR: Could not sanitize molecule on line XXXX" to stderr.
			auto& mol = *mol_ptr;

			// Get the number of points, excluding hydrogens.
			const auto num_points = mol.getNumHeavyAtoms();
			cout << local_time() << "Found " << num_points << " heavy atoms" << endl;

			// Classify subset atoms.
			cout << local_time() << "Classifying atoms into subsets, whose sizes are:";
			for (size_t k = 0; k < num_subsets; ++k)
			{
				vector<vector<pair<int, int>>> matchVect;
				SubstructMatch(mol, *SubsetMols[k], matchVect);
				const auto num_matches = matchVect.size();
				auto& subset = subsets[k];
				subset.resize(num_matches);
				for (size_t i = 0; i < num_matches; ++i)
				{
					subset[i] = matchVect[i].front().second;
				}
				cout << ' ' << subset.size();
			}
			cout << endl;

			// Check user-provided ligand validity.
			const auto& subset0 = subsets.front();
			if (subset0.empty())
			{
				// Record job completion time stamp.
				// TODO: set an error code in the database.
				const auto millis_since_epoch = duration_cast<std::chrono::milliseconds>(system_clock::now().time_since_epoch()).count();
				conn.update(collection, BSON("_id" << _id), BSON("$set" << BSON("done" << Date_t(millis_since_epoch))));
				continue;
			}
			assert(subset0.size() == num_points);

			// Calculate the four reference points.
			cout << local_time() << "Calculating reference points" << endl;
			const auto& conf = mol.getConformer();
			for (auto& ref : references)
			{
				ref.x = ref.y = ref.z = 0;
			}
			auto& ctd = references[0];
			auto& cst = references[1];
			auto& fct = references[2];
			auto& ftf = references[3];
			for (const auto i : subset0)
			{
				const auto& a = conf.getAtomPos(i);
				ctd += a;
			}
			ctd /= num_points;
			double cst_dist = numeric_limits<double>::max();
			double fct_dist = numeric_limits<double>::lowest();
			double ftf_dist = numeric_limits<double>::lowest();
			for (const auto i : subset0)
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
			for (const auto i : subset0)
			{
				const auto& a = conf.getAtomPos(i);
				const auto this_dist = dist2(a, fct);
				if (this_dist > ftf_dist)
				{
					ftf = a;
					ftf_dist = this_dist;
				}
			}

			// Precalculate the distances of heavy atoms to the reference points, given that subsets[1 to 4] are subsets of subsets[0].
			cout << local_time() << "Calculating pairwise distances" << endl;
			for (size_t k = 0; k < num_references; ++k)
			{
				const auto& reference = references[k];
				auto& distp = dista[k];
				distp.resize(num_points);
				for (size_t i = 0; i < num_points; ++i)
				{
					distp[subset0[i]] = sqrt(dist2(conf.getAtomPos(subset0[i]), reference));
				}
			}

			// Loop over pharmacophoric subsets and reference points.
			cout << local_time() << "Calculating USRCAT moments" << endl;
			size_t qo = 0;
			for (const auto& subset : subsets)
			{
				const auto n = subset.size();
				for (size_t k = 0; k < num_references; ++k)
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
			cout << local_time() << "Calculating USR and USRCAT scores" << endl;
			for (auto& ss : scores)
			{
				ss.assign(ss.size(), numeric_limits<double>::max());
			}
			cnt.init(num_chunks);
			for (size_t l = 0; l < num_chunks; ++l)
			{
				io.post([&,l]()
				{
					for (size_t k = chunk_size * l, chunk_end = min<size_t>(k + chunk_size, num_ligands); k < chunk_end; ++k)
					{
						size_t j = k ? mconfss[k - 1] : 0;
						alignas(32) array<double, 4> a;
#ifndef PRELOAD_FEATURES
						alignas(32) array<double, qn.back()> l;
						std::ifstream usrcat_f64("16_usrcat.f64");
						usrcat_f64.seekg(sizeof(l) * j);
#endif
						for (const auto mconfs = mconfss[k]; j < mconfs; ++j)
						{
#ifdef PRELOAD_FEATURES
							const auto& l = features[j];
#else
							assert(usrcat_f64.tellg() == sizeof(l) * j);
							usrcat_f64.read(reinterpret_cast<char*>(l.data()), sizeof(l));
#endif
							double s = 0;
							#pragma unroll
							for (size_t i = 0, u = 0; u < num_usrs; ++u)
							{
								auto& scoreuk = scores[u][k];
								#pragma unroll
								for (const auto qnu = qn[u]; i < qnu; ++i)
								{
									s += abs(q[i] - l[i]);
									if (u == 1 && s >= scoreuk) break;
								}
								if (s < scoreuk)
								{
									scoreuk = s;
									cnfids[u][k] = j;
								}
							}
						}
					}
					cnt.increment();
				});
			}
			cnt.wait();

			// Sort ligands by USRCAT score, if equal then by USR score, if equal then by ZINC ID.
			cout << local_time() << "Sorting " << scase.size() << " scores" << endl;
			iota(scase.begin(), scase.end(), 0);
			sort(scase.begin(), scase.end(), [&](const size_t val0, const size_t val1)
			{
				const auto u0score0 = u0scores[val0];
				const auto u0score1 = u0scores[val1];
				if (u0score0 == u0score1)
				{
					const auto u1score0 = u1scores[val0];
					const auto u1score1 = u1scores[val1];
					if (u1score0 == u1score1)
					{
						return zincids[val0] < zincids[val1];
					}
					return u1score0 < u1score1;
				}
				return u0score0 < u0score1;
			});

			// Create output directory and write output files.
			cout << local_time() << "Creating output directory" << endl;
			const auto output_dir = job_path / to_string(query_number);
			create_directory(output_dir);
			cout << local_time() << "Writing output files" << endl;
			using namespace boost::iostreams;
			filtering_ostream log_csv_gz;
			log_csv_gz.push(gzip_compressor());
			log_csv_gz.push(file_sink((output_dir / "log.csv.gz").string()));
			log_csv_gz.setf(ios::fixed, ios::floatfield);
			log_csv_gz << "ZINC ID,USR score,USRCAT score,Molecular weight (g/mol),Partition coefficient xlogP,Apolar desolvation (kcal/mol),Polar desolvation (kcal/mol),Hydrogen bond donors,Hydrogen bond acceptors,Polar surface area tPSA (A^2),Net charge,Rotatable bonds,SMILES,Substance information,Suppliers and annotations\n";
			filtering_ostream hits_sdf_gz;
			hits_sdf_gz.push(gzip_compressor());
			hits_sdf_gz.push(file_sink((output_dir / "hits.sdf.gz").string()));
			for (size_t t = 0; t < 1000; ++t)
			{
				const auto k = scase[t];
				const auto zincid = zincids[k].substr(0, 8); // Take another substr() to get rid of the trailing newline.
				const auto u0score = 1 / (1 + scores[0][k] * qv[0]);
				const auto u1score = 1 / (1 + scores[1][k] * qv[1]);
				const auto zfp = zfproperties[k];
				const auto zip = ziproperties[k];
				const auto smiles = smileses[k];    // A newline is already included in smileses[k].
				const auto supplier = suppliers[k]; // A newline is already included in suppliers[k].
				log_csv_gz
					<< zincid
					<< setprecision(8)
					<< ',' << u0score
					<< ',' << u1score
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
					<< ",http://zinc.docking.org/substance/" << zincid
					<< ',' << supplier.substr(0, supplier.length() - 1) // Get rid of the trailing newline.
					<< '\n'
				;
				const auto lig = ligands[cnfids[usr][k]];
				hits_sdf_gz.write(lig.data(), lig.size());
			}
		}

		// Update progress.
		cout << local_time() << "Setting done time" << endl;
		const auto millis_since_epoch = duration_cast<std::chrono::milliseconds>(system_clock::now().time_since_epoch()).count();
		conn.update(collection, BSON("_id" << _id), BSON("$set" << BSON("done" << Date_t(millis_since_epoch))));
	}
}
