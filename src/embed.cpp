#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FragCatalog/FragFPGenerator.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
using namespace std;
using namespace RDKit;
using namespace RDKit::MolOps;
using namespace RDKit::DGeomHelpers;

int main(int argc, char* argv[])
{
	// Read a line from standard input.
	string str;
	getline(cin, str);

	// Pipe the line to the supplier.
	istringstream iss(str);
	SmilesMolSupplier sup(&iss, false, "\t", 0, -1, false);

	// Validate the input.
	if (!sup.length() || !sup.atEnd()) return 1;

	// Obtain a pointer to the current molecule with heavy atoms only.
	const unique_ptr<ROMol> smi_ptr(sup.next()); // Calling next() may print "ERROR: Could not sanitize molecule on line XXXX" to stderr.

	// Add hydrogens to the molecule. To get good conformations, it's almost always a good idea to add hydrogens to the molecule first. http://rdkit.org/docs/GettingStartedInPython.html#writing-molecules
	const unique_ptr<ROMol> mol_ptr(addHs(*smi_ptr));

	// Obtain a reference to the molecule to avoid writing *mol_ptr.
	auto& mol = *mol_ptr;

	// Generate conformers with knowledge.
	const auto confId = EmbedMolecule(mol, 0, -1, true, false, 2.0, true, 1, 0, 1e-3, false, true, true, true, false, 5.0);

	// Check if conformers are generated.
	if (confId == -1) return 1;

	// Create output streams.
	SDWriter writer(&cout);
	writer.write(mol, confId);
}
