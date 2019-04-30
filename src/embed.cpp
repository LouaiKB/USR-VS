#include <fstream>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FragCatalog/FragFPGenerator.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
using namespace std;
using namespace RDKit;
using namespace RDKit::MolOps;
using namespace RDKit::DGeomHelpers;
using namespace RDDepict;

int main(int argc, char* argv[])
{
	// Read a line from standard input.
	string str;
	getline(cin, str);

	// Pipe the line to the supplier.
	istringstream iss(str);
	SmilesMolSupplier sup(&iss, false, "\t", 0, 0, false);

	// Validate the input.
	if (!sup.length() || !sup.atEnd()) return 1;

	// Obtain a pointer to the current molecule with heavy atoms only.
	const unique_ptr<ROMol> smi_ptr(sup.next()); // Calling next() may print "ERROR: Could not sanitize molecule on line XXXX" to stderr.

	// Obtain a reference to the molecule to avoid writing *smi_ptr.
	auto& smi = *smi_ptr;

	// Add hydrogens to the molecule. To get good conformations, it's almost always a good idea to add hydrogens to the molecule first. http://rdkit.org/docs/GettingStartedInPython.html#writing-molecules
	const unique_ptr<ROMol> mol_ptr(addHs(smi));

	// Obtain a reference to the molecule to avoid writing *mol_ptr.
	auto& mol = *mol_ptr;

	// Generate a conformer with knowledge.
	EmbedParameters params(ETKDGv2);
	params.randomSeed = 209;
	const auto confId = EmbedMolecule(mol, params); // https://github.com/rdkit/rdkit/pull/1597

	// Check if a conformer was generated.
	if (confId == -1) return 1;

	// Create output streams.
	SDWriter writer(&cout);
	writer.write(mol, confId);

	// Draw SVG.
	if (argc >= 2)
	{
		compute2DCoords(smi);
		ofstream ofs(argv[1]);
		MolDraw2DSVG drawer(600, 600, ofs); // width, height, output
		drawer.drawMolecule(smi);
		drawer.finishDrawing();
	}
}
