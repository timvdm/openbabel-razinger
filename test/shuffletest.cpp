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

void IdsToSymClasses(OBMol *mol, OBTetrahedralStereo::Config &config, 
    const std::vector<unsigned int> &symClasses)
{
  OBAtom *atom;
  // center
  atom = mol->GetAtomById(config.center);
  if (atom) {
    if (atom->IsHydrogen())
      config.center = OBStereo::ImplicitRef;
    else
      config.center = symClasses.at(atom->GetIndex());
  }
  // from/towards
  atom = mol->GetAtomById(config.from);
  if (atom) {
    if (atom->IsHydrogen())
      config.from = OBStereo::ImplicitRef;
    else
      config.from = symClasses.at(atom->GetIndex());
  }
  // refs
  for (unsigned int i = 0; i < config.refs.size(); ++i) {
    atom = mol->GetAtomById(config.refs.at(i));
    if (atom) {
      if (atom->IsHydrogen())
        config.refs[i] = OBStereo::ImplicitRef;
      else
        config.refs[i] = symClasses.at(atom->GetIndex());
    }
  }
}


int configParity(const OBTetrahedralStereo::Config &config, OBMol *mol, 
    const std::vector<unsigned int> &symClasses = std::vector<unsigned int>())
{
  OBTetrahedralStereo::Config cfg = config;
  if (!symClasses.empty())
    IdsToSymClasses(mol, cfg, symClasses);
  //cout << "ERROR  " << cfg << endl;
  std::vector<unsigned long> refs = cfg.refs;
  refs.insert(refs.begin(), cfg.from); // doesn't matter if view is from or towards, will be sorted anyway
  //std::sort(refs.begin(), refs.end());

  bool p = (OBStereo::NumInversions(refs) % 2) ? true : false;
        
  if (p)
    return 1;
  else
    return -1;
}

void permutateConfig(std::vector<OBTetrahedralStereo::Config> &configs, unsigned int k)
{
  OBTetrahedralStereo::Config &config = configs[k];
  OBStereo::Permutate(config.refs, 0, 1);
}
        
void storeConfigs(OBMol *mol, const std::vector<OBTetrahedralStereo::Config> &configs, OBStereoFacade &stereoFacade)
{
  unsigned int idx = 0;
  FOR_ATOMS_OF_MOL (atom, mol) {
    if (stereoFacade.HasTetrahedralStereo(atom->GetId())) {
      OBTetrahedralStereo *ts = stereoFacade.GetTetrahedralStereo(atom->GetId());
      ts->SetConfig(configs.at(idx));
      ++idx;
    }
  }
}

std::string canonicalSmiles(const std::string &smiles)
{
  // read a smiles string
  OBMol mol;
  OBConversion conv, smiConv;
  OB_REQUIRE( conv.SetInFormat("smi") );
  OB_REQUIRE( conv.SetOutFormat("can") );
  // read a smiles string
  OB_REQUIRE( conv.ReadString(&mol, smiles) );

    OBStereoisomer isomers(&mol);
    OBStereoFacade stereoFacade(&mol);
    // print number of stereomers
    cout << "Enantiomer pairs: " << isomers.numEnantiomerPairs() << endl;
    cout << "Diastereomers: " << isomers.numDiastereomers() << endl;

    std::vector<unsigned int> symmetry_classes, canon_order;
    OBGraphSym gs(&mol);
    gs.GetSymmetry(symmetry_classes);
    //gs.GetSymmetry(symmetry_classes, false);
    gs.CanonicalLabels(canon_order);

    std::vector<unsigned long> canIds;

    // print the parities for the molecule
    cout << "ERROR XXX parities sym: ";
    OBStereoisomer::ParityVec parities, lastParities, canParities, idParities;
    std::vector<OBTetrahedralStereo::Config> configs;
    std::vector<unsigned long> atomIds;
    FOR_ATOMS_OF_MOL (atom, mol) {
      if (stereoFacade.HasTetrahedralStereo(atom->GetId())) {
        OBTetrahedralStereo::Config config = stereoFacade.GetTetrahedralStereo(atom->GetId())->GetConfig();
        atomIds.push_back(config.center);
        canIds.push_back(canon_order.at(atom->GetIndex()));
        int symParity = configParity(config, &mol, symmetry_classes);
        parities.push_back(symParity);
        cout << symParity << " ";
        int canParity = configParity(config, &mol, canon_order);
        canParities.push_back(canParity);
        int idParity = configParity(config, &mol);
        idParities.push_back(idParity);
        configs.push_back(config);
      }
    }
    cout << endl;
    cout << "ERROR XXX parities can = ";
    for (unsigned int i = 0; i < canParities.size(); ++i)
      cout << canParities.at(i) << " ";
    cout << endl;
    cout << "ERROR XXX parities id = ";
    for (unsigned int i = 0; i < idParities.size(); ++i)
      cout << idParities.at(i) << " ";
    cout << endl;
 
    cout << "ERROR XXX canIds = ";
    for (unsigned int i = 0; i < canIds.size(); ++i)
      cout << canIds.at(i) << " ";
    cout << endl;
 
    std::vector<std::string> canonicalCandidates;

    lastParities = canParities;
    // print all stereomers
    const std::vector<OBStereoisomer::Diastereomer> &diastereomers = isomers.diastereomers();
    for (unsigned int i = 0; i < diastereomers.size(); ++i) {
      std::vector<std::string> candidates;
      cout << "  diastereomer " << i+1 << endl;
      cout << "ERROR    parities: ";
      bool foundDiastereomer = false;
      for (unsigned int j = 0; j < diastereomers.at(i).parities.size(); ++j) {
        const OBStereoisomer::ParityVec &pv = diastereomers.at(i).parities.at(j);
        cout << "    ";
        for (unsigned int k = 0; k < pv.size(); ++k) {
          if (lastParities.at(k) != pv.at(k)) {
            permutateConfig(configs, k);
          }
          cout << pv.at(k) << " ";
        }
        cout << endl;
        if (pv == canParities) {
          cout << "ERROR  ----> found diastereomer matching input structure: " << smiConv.WriteString(&mol);
          foundDiastereomer = true;
        }
        lastParities = pv;
        storeConfigs(&mol, configs, stereoFacade);
        std::string candidate = conv.WriteString(&mol); 
        cout << candidate;
        candidates.push_back(candidate);
      }

      if (foundDiastereomer)   
        canonicalCandidates = candidates;
    }
    cout << "ERROR XXX: " << canonicalCandidates.size() << "   Inv = " << OBStereo::NumInversions(atomIds) << endl;

    cout << "Canonical candidates:" << endl;
    for (unsigned int i = 0; i < canonicalCandidates.size(); ++i) {
      cout << "ERROR  " << canonicalCandidates.at(i);
    }

    std::sort(canonicalCandidates.begin(), canonicalCandidates.end());
    cout << "ERROR  True canonical SMILES: ";
    cout << canonicalCandidates.front() << endl;


    return canonicalCandidates.front();
}

static unsigned int failed = 0;
static unsigned int testCount = 0;



bool doShuffleTest(const std::string &smiles)
{
  cout << "ERROR True Shuffling: " << smiles << endl;
  cerr << "ERROR True Shuffling: " << smiles << endl;
  // read a smiles string
  OBMol mol;
  OBConversion canConv, smiConv;
  OB_REQUIRE( canConv.SetInFormat("smi") );
  OB_REQUIRE( canConv.SetOutFormat("can") );
  OB_REQUIRE( smiConv.SetOutFormat("smi") );
  // read a smiles string
  OB_REQUIRE( canConv.ReadString(&mol, smiles) );

  int N = 200;
  testCount++;

  std::vector<OBAtom*> atoms;
  FOR_ATOMS_OF_MOL(atom, mol)
    atoms.push_back(&*atom);
      
  std::string ref = canonicalSmiles(smiles);
  cout << "ref = " << ref;
 
  bool result = true;
  for (int i = 0; i < N; ++i) {
    // shuffle the atoms
    std::random_shuffle(atoms.begin(), atoms.end());
    mol.RenumberAtoms(atoms);
    // get can smiles
    std::string smiles = smiConv.WriteString(&mol);
    cerr << smiles;
    std::string cansmi = canonicalSmiles(smiles);
    // comapare with ref
    if (cansmi != ref) {
      cout << "ERROR " << cansmi;
      failed++;
      result = false;
    }
  }

  return result;
}











/*
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
*/

int main(int argc, char **argv)
{
  OB_ASSERT( doShuffleTest("O[C@H]1CC[C@@H](O)CC1") );
  OB_ASSERT( doShuffleTest("O[C@H]1C[C@@H](O)C[C@H](O)C1") );
  OB_ASSERT( doShuffleTest("O[C@H]1C[C@@H](O)C[C@@H](O)C1") );
  //OB_ASSERT( doShuffleTest("O[C@H]1[C@@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H]1O") );
  
  //OB_ASSERT( doShuffleTest("") );

  cout << "FAILED TESTS: " << failed << "/" << testCount << endl;

  return 0;
}

