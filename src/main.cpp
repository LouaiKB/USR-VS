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
#include <immintrin.h>
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
#include <Poco/Net/MailMessage.h>
#include <Poco/Net/MailRecipient.h>
#include <Poco/Net/SMTPClientSession.h>
using namespace std;
using namespace std::chrono;
using namespace RDKit;
using namespace RDGeom;
using namespace boost::filesystem;
using namespace boost::gregorian;
using namespace boost::posix_time;
using namespace mongo;
using namespace bson;
using namespace Poco::Net;

inline static string local_time()
{
	return to_simple_string(microsec_clock::local_time()) + " ";
}

template <typename T>
inline vector<T> read(const path src)
{
	cout << local_time() << "Reading " << src;
	vector<T> buf;
	boost::filesystem::ifstream ifs(src, ios::binary | ios::ate);
	const size_t num_bytes = ifs.tellg();
	cout << ", " << num_bytes << " bytes" << endl;
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
		cout << local_time() << "Reading " << src;
		boost::filesystem::ifstream ifs(src, ios::binary | ios::ate);
		const size_t num_bytes = ifs.tellg();
		cout << ", " << num_bytes << " bytes" << endl;
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
		cout << local_time() << "Reading " << src;
		boost::filesystem::ifstream ifs(src, ios::binary | ios::ate);
		const size_t num_bytes = ifs.tellg();
		cout << ", " << num_bytes << " bytes" << endl;
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
	const auto epoch = date(1970, 1, 1);
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
	const auto m256s = _mm256_set1_pd(-0.); // -0. = 1 << 63

	// Wrap SMARTS strings to RWMol objects.
	array<unique_ptr<ROMol>, num_subsets> SubsetMols;
	for (size_t k = 0; k < num_subsets; ++k)
	{
		SubsetMols[k].reset(reinterpret_cast<ROMol*>(SmartsToMol(SubsetSMARTS[k])));
	}

	// Initialize variables.
	alignas(32) array<double, qn.back()> q;
	alignas(32) array<double, 4> a;

	// Read ZINC ID file.
	const string_array<size_t> zincids("16_zincid.txt");
	const auto num_ligands = zincids.size();

	// Read cumulative number of conformers file.
	const auto mconfss = read<size_t>("16_mconfs.u64");
	const auto num_conformers = mconfss.back();
	assert(mconfss.size() == num_ligands);
	assert(num_conformers >= num_ligands);

	// Read feature file.
	const auto features = read<array<double, qn.back()>>("16_usrcat.f64");
	assert(features.size() == num_conformers);

	// Read ligand header file.
	stream_array<size_t> ligands("16_ligand.sdf");
	assert(ligands.size() == num_conformers);

	array<vector<double>, 2> scores
	{{
		vector<double>(num_ligands, numeric_limits<double>::max()),
		vector<double>(num_ligands, numeric_limits<double>::max())
	}};
	array<vector<size_t>, 2> cnfids
	{{
		vector<size_t>(num_ligands),
		vector<size_t>(num_ligands)
	}};
	const auto& u0scores = scores[0];
	const auto& u1scores = scores[1];
	vector<size_t> scase(num_ligands);

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
		const auto email = job["email"].String();

		// Parse the user-supplied SDF file, setting sanitize=true, removeHs=false, strictParsing=true.
		cout << local_time() << "Parsing the query SDF file" << endl;
		SDMolSupplier sup((job_path / "query.sdf").string(), true, false, true);
		const unique_ptr<ROMol> mol_ptr(sup.next()); // Calling next() may print "ERROR: Could not sanitize molecule on line XXXX" to stderr.
		auto& mol = *mol_ptr;

		// Get the number of points, excluding hydrogens.
		const auto npt = mol.getNumHeavyAtoms();

		// Classify subset atoms.
		cout << local_time() << "Classifying atoms into subsets:";
		array<vector<int>, num_subsets> subsets;
		for (size_t k = 0; k < num_subsets; ++k)
		{
			auto& subset = subsets[k];
			subset.reserve(npt);
			vector<vector<pair<int, int>>> matchVect;
			SubstructMatch(mol, *SubsetMols[k], matchVect);
			for (const auto& v : matchVect)
			{
				subset.push_back(v.front().second);
			}
			cout << ' ' << subset.size();
		}
		cout << endl;
		const auto& subset0 = subsets.front();

		// Check user-provided ligand validity.
		if (subset0.empty())
		{
			// Record job completion time stamp.
			const auto millis_since_epoch = duration_cast<std::chrono::milliseconds>(system_clock::now().time_since_epoch()).count();
			conn.update(collection, BSON("_id" << _id), BSON("$set" << BSON("done" << Date_t(millis_since_epoch))));

			// Send error notification email.
			cout << local_time() << "Sending an error notification email to " << email << endl;
			MailMessage message;
			message.setSender("usr <noreply@cse.cuhk.edu.hk>");
			message.setSubject("Your usr job has failed");
			message.setContent("Description: " + job["description"].String() + "\nSubmitted: " + to_simple_string(ptime(epoch, boost::posix_time::milliseconds(job["submitted"].Date().millis))) + " UTC\nFailed: " + to_simple_string(ptime(epoch, boost::posix_time::milliseconds(millis_since_epoch))) + " UTC\nReason: failed to parse the provided ligand.");
			message.addRecipient(MailRecipient(MailRecipient::PRIMARY_RECIPIENT, email));
			SMTPClientSession session("137.189.91.190");
			session.login();
			session.sendMessage(message);
			session.close();
			continue;
		}

		// Calculate the four reference points.
		cout << local_time() << "Calculating reference points" << endl;
		const auto& conf = mol.getConformer();
		const auto n = subset0.size();
		assert(n == npt);
		const auto v = 1.0 / n;
		array<Point3D, num_references> references{};
		auto& ctd = references[0];
		auto& cst = references[1];
		auto& fct = references[2];
		auto& ftf = references[3];
		for (const auto i : subset0)
		{
			const auto& a = conf.getAtomPos(i);
			ctd += a;
		}
		ctd *= v;
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
		array<vector<double>, num_references> dista;
		for (size_t k = 0; k < num_references; ++k)
		{
			const auto& reference = references[k];
			auto& dists = dista[k];
			dists.resize(n);
			for (size_t i = 0; i < n; ++i)
			{
				dists[subset0[i]] = sqrt(dist2(conf.getAtomPos(subset0[i]), reference));
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
		cout << local_time() << "Calculating USRCAT scores" << endl;
		for (size_t j = 0, k = 0; k < num_ligands; ++k)
		{
			for (const auto mconfs = mconfss[k]; j < mconfs; ++j)
			{
				const auto& l = features[j];
				double s = 0;
				#pragma unroll
				for (size_t i = 0, u = 0; u < num_usrs; ++u)
				{
					#pragma unroll
					for (const auto qnu = qn[u]; i < qnu; i += 4)
					{
						const auto m256a = _mm256_andnot_pd(m256s, _mm256_sub_pd(_mm256_load_pd(&q[i]), _mm256_load_pd(&l[i])));
						_mm256_stream_pd(a.data(), _mm256_hadd_pd(m256a, m256a));
						s += a[0] + a[2];
					}
					if (s < scores[u][k])
					{
						scores[u][k] = s;
						cnfids[u][k] = j;
					}
				}
			}
		}

		// Sort ligands by USRCAT score, if equal then by USR score, if equal then by ZINC ID.
		cout << local_time() << "Sorting scores" << endl;
		iota(scase.begin(), scase.end(), 0);
		sort(scase.begin(), scase.end(), [&](const size_t val0, const size_t val1)
		{
			const auto u1score0 = u1scores[val0];
			const auto u1score1 = u1scores[val1];
			if (u1score0 == u1score1)
			{
				const auto u0score0 = u0scores[val0];
				const auto u0score1 = u0scores[val1];
				if (u0score0 == u0score1)
				{
					return zincids[val0] < zincids[val1];
				}
				return u0score0 < u0score1;
			}
			return u1score0 < u1score1;
		});

		// Write results.
		cout << local_time() << "Writing output files" << endl;
		using namespace boost::iostreams;
		filtering_ostream log_csv_gz;
		log_csv_gz.push(gzip_compressor());
		log_csv_gz.push(file_sink((job_path / "log.csv.gz").string()));
		log_csv_gz.setf(ios::fixed, ios::floatfield);
		log_csv_gz << "ZINC ID,USR score,USRCAT score\n" << setprecision(8);
		filtering_ostream hits_sdf_gz;
		hits_sdf_gz.push(gzip_compressor());
		hits_sdf_gz.push(file_sink((job_path / "hits.sdf.gz").string()));
		for (size_t t = 0; t < 100000; ++t)
		{
			const auto k = scase[t];
			const auto zincid = zincids[k].substr(0, 8); // Take another substr() to get rid of the trailing newline.
			const auto u0score = 1 / (1 + scores[0][k] * qv[0]);
			const auto u1score = 1 / (1 + scores[1][k] * qv[1]);
			log_csv_gz << zincid << ',' << u0score << ',' << u1score << '\n';

			// Only write conformations of the top ligands to ligands.pdbqt.gz.
			if (t >= 1000) continue;

			const auto lig = ligands[cnfids[1][k]];
			hits_sdf_gz.write(lig.data(), lig.size());
		}

		// Update progress.
		cout << local_time() << "Setting done time" << endl;
		const auto millis_since_epoch = duration_cast<std::chrono::milliseconds>(system_clock::now().time_since_epoch()).count();
		conn.update(collection, BSON("_id" << _id), BSON("$set" << BSON("done" << Date_t(millis_since_epoch))));

		// Send completion notification email.
		cout << local_time() << "Sending a completion notification email to " << email << endl;
		MailMessage message;
		message.setSender("istar <noreply@cse.cuhk.edu.hk>");
		message.setSubject("Your usr job has completed");
		message.setContent("Description: " + job["description"].String() + "\nSubmitted: " + to_simple_string(ptime(epoch, boost::posix_time::milliseconds(job["submitted"].Date().millis))) + " UTC\nCompleted: " + to_simple_string(ptime(epoch, boost::posix_time::milliseconds(millis_since_epoch))) + " UTC\nResult: http://istar.cse.cuhk.edu.hk/usr/iview/?" + _id.str());
		message.addRecipient(MailRecipient(MailRecipient::PRIMARY_RECIPIENT, email));
		SMTPClientSession session("137.189.91.190");
		session.login();
		session.sendMessage(message);
		session.close();
	}
}
