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


int configParity(const OBTetrahedralStereo::Config &config, OBMol *mol, const std::vector<unsigned int> &symClasses)
{
  OBTetrahedralStereo::Config cfg = config;
  //IdsToSymClasses(mol, cfg, symClasses);



  OBTetrahedralStereo::Config orderedConfig;    // constructor initializes view & winding --> canonical
  // copy specified flag and center Ref
  orderedConfig.center = cfg.center;
  orderedConfig.specified = cfg.specified;

  std::vector<unsigned long> refs = cfg.refs;
  refs.insert(refs.begin(), cfg.from); // doesn't matter if view is from or towards, will be sorted anyway
  std::sort(refs.begin(), refs.end());
  for (unsigned int i = 0; i < refs.size(); ++i)
    if (i)
      orderedConfig.refs.push_back(refs.at(i));
    else
      orderedConfig.from = refs.at(0);
  // store the match/mismatch result
  bool p = (cfg == orderedConfig) ? true : false;
/*

  //cout << "ERROR  " << cfg << endl;
  std::vector<unsigned long> refs = cfg.refs;
  refs.insert(refs.begin(), config.from); // doesn't matter if view is from or towards, will be sorted anyway
  //std::sort(refs.begin(), refs.end());

  bool p = (OBStereo::NumInversions(refs) % 2) ? true : false;
  */      
  if (p)
    return -1;
  else
    return 1;
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
  OBConversion conv;
  OB_REQUIRE( conv.SetInFormat("smi") );
  OB_REQUIRE( conv.SetOutFormat("can") );
  // read a smiles string
  OB_REQUIRE( conv.ReadString(&mol, smiles) );

 



  OBStereoisomer isomers(&mol);
    OBStereoFacade stereoFacade(&mol);
    // print number of stereomers
    cout << "Enantiomer pairs: " << isomers.numEnantiomerPairs() << endl;
    cout << "Diastereomers: " << isomers.numDiastereomers() << endl;

    std::vector<unsigned int> symmetry_classes;
    OBGraphSym gs(&mol);
    gs.GetSymmetry(symmetry_classes);
//    gs.CanonicalLabels(symmetry_classes);

    // print the parities for the molecule
    cout << "parities: ";
    OBStereoisomer::ParityVec parities, lastParities;
    std::vector<OBTetrahedralStereo::Config> configs;
    std::vector<unsigned long> atomIds;
    FOR_ATOMS_OF_MOL (atom, mol) {
      if (stereoFacade.HasTetrahedralStereo(atom->GetId())) {
        cout << "(id:" << atom->GetId() << ") ";
        OBTetrahedralStereo::Config config = stereoFacade.GetTetrahedralStereo(atom->GetId())->GetConfig();
        atomIds.push_back(config.center);
        int parity = configParity(config, &mol, symmetry_classes);
        parities.push_back(parity);
        cout << parity << " ";
        configs.push_back(config);
      }
    }
    cout << endl;
 
    std::vector<std::string> canonicalCandidates;

    lastParities = parities;
    // print all stereomers
    const std::vector<OBStereoisomer::Enantiomer> &enantiomers = isomers.enantiomers();
    for (unsigned int i = 0; i < enantiomers.size(); ++i) {
      std::vector<std::string> candidates;
      std::vector<std::string> invCandidates;
      cout << "  enantiomer " << i+1 << endl;
      cout << "ERROR    parities: ";
      bool foundEnantiomer = false;
      for (unsigned int j = 0; j < enantiomers.at(i).parities.size(); ++j) {
        const OBStereoisomer::ParityVec &pv = enantiomers.at(i).parities.at(j);
        cout << "    ";
        for (unsigned int k = 0; k < pv.size(); ++k) {
          if (lastParities.at(k) != pv.at(k)) {
            permutateConfig(configs, k);
          }
          cout << pv.at(k) << " ";
        }
        cout << endl;
        if (pv == parities) {
          cout << "ERROR  ----> found enantiomer matching input structure: " << conv.WriteString(&mol);
          foundEnantiomer = true;
        }
        lastParities = pv;
        storeConfigs(&mol, configs, stereoFacade);
        std::string candidate = conv.WriteString(&mol); 
        cout << candidate;
        candidates.push_back(candidate);
      }
      bool foundInvEnantiomer = false;
      for (unsigned int j = 0; j < enantiomers.at(i).inverseParities.size(); ++j) {
        const OBStereoisomer::ParityVec &pv = enantiomers.at(i).inverseParities.at(j);
        cout << "    ";
        for (unsigned int k = 0; k < pv.size(); ++k) {
          if (lastParities.at(k) != pv.at(k)) {
            permutateConfig(configs, k);
          }
          cout << pv.at(k) << " ";
        }
        cout << endl;
        if (pv == parities) {
          cout << "ERROR  ----> found enantiomer matching input structure: " << conv.WriteString(&mol);
          foundInvEnantiomer = true;
        }
        lastParities = pv;
        storeConfigs(&mol, configs, stereoFacade);
        std::string candidate = conv.WriteString(&mol); 
        cout << candidate;
        invCandidates.push_back(candidate);
      }

      if (foundEnantiomer)
        canonicalCandidates = candidates;
      if (foundInvEnantiomer)
        canonicalCandidates = invCandidates;

    }
 
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
        if (pv == parities) {
          cout << "ERROR  ----> found diastereomer matching input structure: " << conv.WriteString(&mol);
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
  
  for (int i = 0; i < N; ++i) {
    // shuffle the atoms
    std::random_shuffle(atoms.begin(), atoms.end());
    mol.RenumberAtoms(atoms);
    // get can smiles
    std::string smiles = smiConv.WriteString(&mol);
    std::string cansmi = canonicalSmiles(smiles);
    // comapare with ref
    if (cansmi != ref) {
      cout << "ERROR " << cansmi;
      failed++;
      return false;
    }
  }

  return true;
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

