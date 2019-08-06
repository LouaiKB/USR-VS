#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <filesystem>
#include <vector>
#include <array>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <thread>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <Numerics/Alignment/AlignPoints.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <mongocxx/instance.hpp>
#include <mongocxx/pool.hpp>
#include <mongocxx/client.hpp>
#include <bsoncxx/json.hpp>
#include "io_service_pool.hpp"
#include "safe_counter.hpp"
using namespace std;
using namespace std::chrono;
using namespace std::filesystem;
using namespace RDKit;
using namespace RDKit::MolOps;
using namespace RDDepict;
using namespace RDKit::MorganFingerprints;
using namespace RDGeom;
using namespace RDNumeric::Alignments;
using namespace MolTransforms;
using namespace boost::posix_time;
using namespace mongocxx;
using bsoncxx::builder::basic::kvp;

inline static auto local_time()
{
	return to_simple_string(microsec_clock::local_time()) + " "; // TODO: update to std::chrono::format(std::chrono::system_clock::now()) when this c++20 feature is implemented in gcc or clang
}

template <typename T>
inline vector<T> read(const path src)
{
	ifstream ifs(src, ios::binary | ios::ate);
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
		ifstream ifs(src, ios::binary | ios::ate);
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
		ifstream ifs(src, ios::binary | ios::ate);
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

template<typename T>
auto dist2(const T& p0, const T& p1)
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
	auto cst_dist = numeric_limits<double>::max();
	auto fct_dist = numeric_limits<double>::lowest();
	auto ftf_dist = numeric_limits<double>::lowest();
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

	// Connect to host and authenticate user.
	cout << local_time() << "Connecting to " << host << " and authenticating " << user << endl;
	const instance inst;
	const uri uri("mongodb://localhost:27017/?minPoolSize=0&maxPoolSize=2"); // When connecting to a replica set, it is much more efficient to use a pool as opposed to manually constructing client objects.
	pool pool(uri);
	const auto client = pool.acquire(); // Return value of acquire() is an instance of entry. An entry is a handle on a client object acquired via the pool.
	const auto db = client->database("jstar");
	auto coll = db.collection("usr2");
	const auto jobid_filter = bsoncxx::from_json(R"({ "started" : { "$exists" : false }})");
	const auto jobid_foau_options = options::find_one_and_update().sort(bsoncxx::from_json(R"({ "submitted" : 1 })")).projection(bsoncxx::from_json(R"({ "_id" : 1, "usr": 1 })")); // By default, the original document is returned

	// Initialize constants.
	cout << local_time() << "Initializing" << endl;
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

	// Calculate number of conformers.
	const auto num_conformers = num_ligands << 2;
	cout << local_time() << "Found " << num_conformers << " database conformers" << endl;

	// Read feature file.
	const auto features = read<array<double, qn.back()>>("16/usrcat.f64");
	assert(features.size() == num_conformers);

	// Initialize variables.
	array<vector<int>, num_subsets> subsets;
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
		const auto started = system_clock::now();
		bsoncxx::builder::basic::document jobid_update_builder;
		jobid_update_builder.append(
			kvp("$set", [=](bsoncxx::builder::basic::sub_document set_subdoc) {
				set_subdoc.append(kvp("started", bsoncxx::types::b_date(started)));
			})
		);
		const auto jobid_document = coll.find_one_and_update(jobid_filter.view(), jobid_update_builder.extract(), jobid_foau_options);
		if (!jobid_document)
		{
			// No incompleted jobs. Sleep for a while.
			if (!sleeping) cout << local_time() << "Sleeping" << endl;
			sleeping = true;
			this_thread::sleep_for(std::chrono::seconds(2));
			continue;
		}
		sleeping = false;
		const auto jobid_view = jobid_document->view();

		// Obtain job properties.
		const auto _id = jobid_view["_id"].get_oid().value;
		cout << local_time() << "Executing job " << _id.to_string() << endl;
		const auto job_path = jobs_path / _id.to_string();
		const size_t usr0 = jobid_view["usr"].get_int32(); // Specify the primary sorting score. 0: USR; 1: USRCAT.
		assert(usr0 == 0 || usr0 == 1);
		const auto usr1 = usr0 ^ 1;
		const auto qnu0 = qn[usr0];
		const auto qnu1 = qn[usr1];

		// Read and validate the user-supplied SDF file.
		cout << local_time() << "Reading and validating the query file" << endl;
		SDMolSupplier sup((job_path / "query.sdf").string(), true, false, true); // sanitize, removeHs, strictParsing. Note: setting removeHs=true (which is the default setting) will lead to fewer hydrogen bond acceptors being matched.
		if (!sup.length() || !sup.atEnd())
		{
			const auto error = 1;
			cout << local_time() << "Failed to parse the query file, error code = " << error << endl;
			bsoncxx::builder::basic::document compt_update_builder;
			compt_update_builder.append(
				kvp("$set", [=](bsoncxx::builder::basic::sub_document set_subdoc) {
					set_subdoc.append(kvp("completed", bsoncxx::types::b_date(started)));
					set_subdoc.append(kvp("error", error));
				})
			);
			const auto compt_update = coll.update_one(bsoncxx::builder::basic::make_document(kvp("_id", _id)), compt_update_builder.extract(), options::update()); // stdx::optional<result::update>. options: write_concern
			assert(compt_update);
			assert(compt_update->matched_count() == 1);
			assert(compt_update->modified_count() == 1);
			continue;
		}

		// Process each of the query molecules sequentially.
		const auto num_queries = 1; // Restrict the number of query molecules to 1. Setting num_queries = sup.length() to execute any number of query molecules.
		for (unsigned int query_number = 0; query_number < num_queries; ++query_number)
		{
			cout << local_time() << "Parsing query molecule " << query_number << endl;
			const unique_ptr<ROMol> qry_ptr(sup.next()); // Calling next() may print "ERROR: Could not sanitize molecule on line XXXX" to stderr.
			auto& qryMol = *qry_ptr;

			// Get the number of atoms, including and excluding hydrogens.
			const auto num_atoms = qryMol.getNumAtoms();
			const auto num_heavy_atoms = qryMol.getNumHeavyAtoms();
			assert(num_heavy_atoms);
			cout << local_time() << "Found " << num_atoms << " atoms and " << num_heavy_atoms << " heavy atoms" << endl;

			// Create an output directory.
			cout << local_time() << "Creating output directory" << endl;
			const auto output_dir = job_path / to_string(query_number);
			create_directory(output_dir);

			// Draw a SVG.
			cout << local_time() << "Drawing a SVG" << endl;
			{
				const unique_ptr<ROMol> qrz_ptr(removeHs(qryMol));
				auto& qrzMol = *qrz_ptr;
				compute2DCoords(qrzMol);
				ofstream ofs(output_dir / "query.svg");
				MolDraw2DSVG drawer(600, 600, ofs); // width, height, output
				drawer.drawMolecule(qrzMol);
				drawer.finishDrawing();
			}

			// Calculate Morgan fingerprint.
			cout << local_time() << "Calculating Morgan fingerprint" << endl;
			const unique_ptr<SparseIntVect<uint32_t>> qryFp(getFingerprint(qryMol, 2));

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
			assert(subset0.size() == num_heavy_atoms);

			// Calculate the four reference points.
			cout << local_time() << "Calculating " << num_refPoints << " reference points" << endl;
			const auto qryRefPoints = calcRefPoints(qryMol, subset0);
			const Point3DConstPtrVect qryRefPointv
			{{
				&qryRefPoints[0],
				&qryRefPoints[1],
				&qryRefPoints[2],
				&qryRefPoints[3],
			}};

			// Precalculate the distances of heavy atoms to the reference points, given that subsets[1 to 4] are subsets of subsets[0].
			cout << local_time() << "Calculating " << num_heavy_atoms * num_refPoints << " pairwise distances" << endl;
			const auto& qryCnf = qryMol.getConformer();
			for (size_t k = 0; k < num_refPoints; ++k)
			{
				const auto& refPoint = qryRefPoints[k];
				auto& distp = dista[k];
				distp.resize(num_atoms);
				for (size_t i = 0; i < num_heavy_atoms; ++i)
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
						for (size_t j = k << 2; j < (k + 1) << 2; ++j)
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
			cout << local_time() << "Writing output files" << endl;
			SDWriter hits_sdf((output_dir / "hits.sdf").string());
			ofstream hits_csv(output_dir / "hits.csv");
			hits_csv.setf(ios::fixed, ios::floatfield);
			hits_csv << "ZINC ID,USR score,USRCAT score,2D Tanimoto score,Molecular weight (g/mol),Partition coefficient xlogP,Apolar desolvation (kcal/mol),Polar desolvation (kcal/mol),Hydrogen bond donors,Hydrogen bond acceptors,Polar surface area tPSA (Å^2),Net charge,Rotatable bonds,SMILES,Vendors and annotations\n";
			for (size_t l = 0; l < num_hits; ++l)
			{
				// Obtain indexes to the hit molecule and the hit conformer.
				const auto k = zcase[l];
				const auto j = cnfids[k];

				// Read SDF content of the hit conformer.
				// TODO: query the SDF conformer from MongoDB using _id, which is zincids[k]
				const auto lig = "";

				// Construct a RDKit ROMol object.
				istringstream iss(lig);
				SDMolSupplier sup(&iss, false, true, false, true);
				assert(sup.length() == 1);
				assert(sup.atEnd());
				const unique_ptr<ROMol> hit_ptr(sup.next());
				auto& hitMol = *hit_ptr;

				// Calculate Morgan fingerprint.
				const unique_ptr<SparseIntVect<uint32_t>> hitFp(getFingerprint(hitMol, 2));

				// Calculate Tanimoto similarity.
				const auto ts = TanimotoSimilarity(*qryFp, *hitFp);

				// Find heavy atoms.
				vector<vector<pair<int, int>>> matchVect;
				SubstructMatch(hitMol, *SubsetMols[0], matchVect);
				const auto num_matches = matchVect.size();
				assert(num_matches == hitMol.getNumHeavyAtoms());
				vector<int> hitHeavyAtoms(num_matches);
				for (size_t i = 0; i < num_matches; ++i)
				{
					hitHeavyAtoms[i] = matchVect[i].front().second;
					assert(hitHeavyAtoms[i] == i); // hitHeavyAtoms can be constructed using iota(hitHeavyAtoms.begin(), hitHeavyAtoms.end(), 0); because for RDKit-generated SDF molecules, heavy atom are always the first few atoms.
				}

				// Calculate the four reference points.
				const auto hitRefPoints = calcRefPoints(hitMol, hitHeavyAtoms);
				const Point3DConstPtrVect hitRefPointv
				{{
					&hitRefPoints[0],
					&hitRefPoints[1],
					&hitRefPoints[2],
					&hitRefPoints[3],
				}};

				// Calculate a 3D transform from the four reference points of the hit conformer to those of the query molecule.
				Transform3D trans;
				AlignPoints(qryRefPointv, hitRefPointv, trans);

				// Apply the 3D transform to all atoms of the hit conformer.
				auto& hitCnf = hitMol.getConformer();
				transformConformer(hitCnf, trans);

				// Write the aligned hit conformer.
				hits_sdf.write(hitMol);

				// Calculate the secondary score of the saved conformer, which has the best primary score.
				const auto& d = features[j];
				double s = 0;
				for (size_t i = 0; i < qnu1; ++i)
				{
					s += abs(q[i] - d[i]);
				}

				const auto u0score = 1 / (1 + scores[k] * qv[usr0]); // Primary score of the current molecule.
				const auto u1score = 1 / (1 + s         * qv[usr1]); // Secondary score of the current molecule.
				const auto zincid = zincids[k].substr(0, 8); // Take another substr() to get rid of the trailing newline.
				hits_csv
					<< zincid
					<< setprecision(8)
					<< ',' << (usr1 ? u0score : u1score)
					<< ',' << (usr1 ? u1score : u0score)
					<< ',' << ts
					<< '\n'
				;
			}
		}

		// Update job status.
		cout << local_time() << "Setting completed time" << endl;
		const auto completed = system_clock::now();
		bsoncxx::builder::basic::document compt_update_builder;
		compt_update_builder.append(
			kvp("$set", [=](bsoncxx::builder::basic::sub_document set_subdoc) {
				set_subdoc.append(kvp("completed", bsoncxx::types::b_date(completed)));
				set_subdoc.append(kvp("nqueries", num_queries));
			})
		);
		const auto compt_update = coll.update_one(bsoncxx::builder::basic::make_document(kvp("_id", _id)), compt_update_builder.extract(), options::update()); // stdx::optional<result::update>. options: write_concern
		assert(compt_update);
		assert(compt_update->matched_count() == 1);
		assert(compt_update->modified_count() == 1);

		// Calculate runtime in seconds and screening speed in million conformers per second.
		const auto runtime = (completed - started).count() * 1e-9; // in seconds
		const auto speed = num_conformers * 1e-6 * num_queries / runtime;
		cout
			<< local_time() << "Completed " << num_queries << " " << (num_queries == 1 ? "query" : "queries") << " in " << setprecision(3) << runtime << " seconds" << endl
			<< local_time() << "Screening speed was " << setprecision(0) << speed << " M conformers per second" << endl
		;
	}
}
