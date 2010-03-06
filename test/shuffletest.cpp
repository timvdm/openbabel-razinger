#include "obtest.h"
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

#include <iostream>
#include <vector>
#include <algorithm>


//#define WITH_DEPICTION 1

using std::cout;
using std::endl;
using namespace OpenBabel;

static unsigned int failed = 0;
static unsigned int count = 0;

bool doShuffleTest(const std::string &smiles)
{
  cout << "Shuffling: " << smiles << endl;
  // read a smiles string
  OBMol mol;
  OBConversion conv;
  OB_REQUIRE( conv.SetInFormat("smi") );
  OB_REQUIRE( conv.SetOutFormat("can") );
  // read a smiles string
  OB_REQUIRE( conv.ReadString(&mol, smiles) );

  int N = 200;
  unsigned int errors = 0;
  count++;

  std::vector<OBAtom*> atoms;
  FOR_ATOMS_OF_MOL(atom, mol)
    atoms.push_back(&*atom);
      
  std::string ref = conv.WriteString(&mol);
  cout << "ref = " << ref;
  
  for (int i = 0; i < N; ++i) {
    // shuffle the atoms
    std::random_shuffle(atoms.begin(), atoms.end());
    mol.RenumberAtoms(atoms);
    // get can smiles
    std::string cansmi = conv.WriteString(&mol);
    // comapare with ref
    if (cansmi != ref) {
      cout << "ERROR " << cansmi;
      failed++;
      return false;
    }
  }

  return true;
}

int main(int argc, char **argv)
{
  OB_ASSERT( doShuffleTest("O[C@H]1CC[C@@H](O)CC1") );
  OB_ASSERT( doShuffleTest("O[C@H]1C[C@@H](O)C[C@H](O)C1") );
  OB_ASSERT( doShuffleTest("O[C@H]1C[C@@H](O)C[C@@H](O)C1") );
  OB_ASSERT( doShuffleTest("O[C@H]1[C@@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H]1O") );
  
  //OB_ASSERT( doShuffleTest("") );

  cout << "FAILED TESTS: " << failed << "/" << count << endl;

  return 0;
}

