#include "obtest.h"
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

#include <openbabel/graphsym.h>
#include <openbabel/stereo/tetrahedral.h>
#include <openbabel/stereo/stereoisomer.h>

#include <iostream>
#include <vector>
#include <algorithm>

//#define WITH_DEPICTION 1

using std::cout;
using std::endl;
using namespace OpenBabel;

std::string GetFilename(const std::string &filename)
{
  string path = TESTDATADIR + filename;
  return path;
}

static unsigned int failed = 0;
static unsigned int testCount = 0;



bool doShuffleTest(const std::string &smiles)
{
  cout << " True Shuffling: " << smiles << endl;
  // read a smiles string
  OBMol mol;
  OBConversion canConv, smiConv;
  OB_REQUIRE( canConv.SetInFormat("smi") );
  OB_REQUIRE( canConv.SetOutFormat("can") );
  OB_REQUIRE( smiConv.SetOutFormat("smi") );
  // read a smiles string
  OB_REQUIRE( canConv.ReadString(&mol, smiles) );

  int N = 20;
  testCount++;

  std::vector<OBAtom*> atoms;
  FOR_ATOMS_OF_MOL(atom, mol)
    atoms.push_back(&*atom);
  
  std::string ref = canConv.WriteString(&mol); // FIXME
  cout << "ref = " << ref;
 
  bool result = true;
  for (int i = 0; i < N; ++i) {
    // shuffle the atoms
    std::random_shuffle(atoms.begin(), atoms.end());
    mol.RenumberAtoms(atoms);
    // get can smiles
    std::string cansmi = canConv.WriteString(&mol);
    OB_ASSERT( cansmi == ref );
    // comapare with ref
    if (cansmi != ref) {
      cout << " " << cansmi;
      if (result)
        failed++;
      result = false;
    }
  }

  return result;
}

bool doShuffleTestFile(const std::string &filename)
{
  cout << " True Shuffling: " << filename << endl;
  std::string file = GetFilename(filename);
  // read a smiles string
  OBMol mol;
  OBConversion canConv, smiConv;
  OBFormat *format = canConv.FormatFromExt(file.c_str());
  OB_REQUIRE( format );
  OB_REQUIRE( canConv.SetInFormat(format) );
  OB_REQUIRE( canConv.ReadFile(&mol, file) );
  OB_REQUIRE( canConv.SetOutFormat("can") );
  OB_REQUIRE( smiConv.SetOutFormat("smi") );

  std::string smiles = canConv.WriteString(&mol);
  int N = 200;
  testCount++;

  std::vector<OBAtom*> atoms;
  FOR_ATOMS_OF_MOL(atom, mol)
    atoms.push_back(&*atom);
  
  std::string ref = canConv.WriteString(&mol); // FIXME
  cout << "ref = " << ref;
 
  bool result = true;
  for (int i = 0; i < N; ++i) {
    // shuffle the atoms
    std::random_shuffle(atoms.begin(), atoms.end());
    mol.RenumberAtoms(atoms);
    // get can smiles
    std::string cansmi = canConv.WriteString(&mol); // FIXME
    OB_ASSERT( cansmi == ref );
    // comapare with ref
    if (cansmi != ref) {
      cout << " " << cansmi;
      if (result)
        failed++;
      result = false;
    }
  }

  return result;
}

bool doShuffleTestMultiFile(const std::string &filename)
{
  cout << " True Shuffling: " << filename << endl;
  std::string file = GetFilename(filename);
  // read a smiles string
  OBMol mol;
  OBConversion canConv;
  OBFormat *format = canConv.FormatFromExt(file.c_str());
  OB_REQUIRE( format );
  OB_REQUIRE( canConv.SetInFormat(format) );
  OB_REQUIRE( canConv.SetOutFormat("can") );

  testCount++;

  std::ifstream ifs;
  ifs.open(file.c_str());
  OB_REQUIRE( ifs );
  OB_REQUIRE( canConv.Read(&mol, &ifs) );

  std::string ref = canConv.WriteString(&mol); // FIXME
  cout << "ref = " << ref;
 
  bool result = true;
  while (canConv.Read(&mol, &ifs)) {
    // get can smiles
    std::string cansmi = canConv.WriteString(&mol); // FIXME
    OB_ASSERT( cansmi == ref );
    // comapare with ref
    if (cansmi != ref) {
      cout << " " << cansmi;
      if (result)
        failed++;
      result = false;
    }
  }

  return result;
}

int main(int argc, char **argv)
{
  if (argc == 2) {
    OB_ASSERT( doShuffleTestFile(argv[1]) );
    return 0;
  }

  OB_ASSERT( doShuffleTestMultiFile("stereo/shuffle_multi1.smi") );
  OB_ASSERT( doShuffleTestMultiFile("stereo/shuffle_multi2.smi") );

  OB_ASSERT( doShuffleTest("O[C@H]1CC[C@@H](O)CC1") );
  OB_ASSERT( doShuffleTest("O[C@H]1C[C@@H](O)C[C@H](O)C1") );
  OB_ASSERT( doShuffleTest("O[C@H]1C[C@@H](O)C[C@@H](O)C1") );
  
  OB_ASSERT( doShuffleTest("[C@@H]1([C@H]([C@H]([C@H]1C)C)C)C") );
  OB_ASSERT( doShuffleTestFile("stereo/cyclobutane_D1.smi") );
  
  OB_ASSERT( doShuffleTestFile("stereo/razinger_fig7_5_spec.mol") );
  
  OB_ASSERT( doShuffleTestFile("stereo/cyclohexanediol_D1.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/cyclohexanediol_D2.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/cyclohexanetriol_D1.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/cyclohexanetriol_D2.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/cyclobutane_D1.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/cyclobutane_D2.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/cyclobutane_D3.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/cyclobutane_D4.mol") );
  OB_ASSERT( doShuffleTest("[C@@H]1([C@H]([C@H]([C@H]1C)C)C)C") );
  OB_ASSERT( doShuffleTest("[C@@H]1([C@@H]([C@H]([C@H]1C)C)C)C") );	
  OB_ASSERT( doShuffleTest("[C@@H]1([C@@H]([C@H]([C@@H]1C)C)C)C") );	
  OB_ASSERT( doShuffleTest("[C@@H]1([C@@H]([C@@H]([C@H]1C)C)C)C") );	

  OB_ASSERT( doShuffleTestFile("stereo/inositol_cis.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/inositol_epi.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/inositol_allo.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/inositol_myo.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/inositol_muco.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/inositol_neo.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/inositol_scyllo.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/inositol_chiroD.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/inositol_chiroL.mol") );
  
  OB_ASSERT( doShuffleTestFile("stereo/razinger_fig7_19_spec1.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/razinger_fig7_19_spec2.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/razinger_fig7_19_spec3.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/razinger_fig7_19_spec4.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/razinger_fig7_19_spec5.mol") );

  OB_ASSERT( doShuffleTestFile("stereo/razinger_fig7_20_spec1.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/razinger_fig7_26_spec1.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/razinger_fig7_59_spec1.mol") );
  
  OB_ASSERT( doShuffleTest("O[C@H]1[C@@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H]1O") );
  
  OB_ASSERT( doShuffleTestFile("stereo/razinger_fig7_30.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/razinger_fig7_34.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/razinger_fig7_35.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/razinger_fig7_36.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/razinger_fig7_37.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/razinger_fig7_38.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/razinger_fig7_39.mol") );

  OB_ASSERT( doShuffleTestFile("stereo/canon1.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/canon2.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/canon3.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/canon4.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/canon5.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/canon6.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/canon7.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/canon8.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/canon9.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/canon10.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/canon11.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/canon12.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/canon13.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/canon14.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/canon15.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/canon16.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/canon17.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/canon18.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/canon19.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/canon20.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/canon21.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/canon22.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/canon23.mol") );


  //OB_ASSERT( doShuffleTest("") );

  cout << "PASSED TESTS: " << testCount - failed << "/" << testCount << endl;

  return 0;
}

