#include <GraphMol/FileParsers/MolSupplier.h>
using namespace std;
using namespace RDKit;

int main(int argc, char* argv[])
{
	// Read contents from standard input.
	stringstream ss;
	for (string line; getline(cin, line); ss << line << endl);

	// Construct a molecule supplier.
	SDMolSupplier sup(&ss, false, true, false, true); // sanitize, removeHs, strictParsing

	// Ensure there is at least a molecule and the cursor is at the end.
	if (!sup.length() || !sup.atEnd()) return 1;

	// Try parsing the molecule.
	const unique_ptr<ROMol> qry_ptr(sup.next()); // Calling next() may print "ERROR: Could not sanitize molecule on line XXXX" to stderr.

	// Obtain a reference to the molecule.
	const auto& qryMol = *qry_ptr;

	// Get the number of heavy atoms, excluding hydrogens.
	const auto num_points = qryMol.getNumHeavyAtoms();

	// Ensure the molecule contains some heavy atoms.
	if (!num_points) return 1;
}
