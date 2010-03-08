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


int configParityPaper(OBTetrahedralStereo *ts, OBMol *mol, 
    const std::vector<unsigned int> &canon_order)
{
  OBTetrahedralStereo::Config cfg = ts->GetConfig();
  IdsToSymClasses(mol, cfg, canon_order); // FIXME

  // lowest priority = highest canonical
  // H = lowest priority
  unsigned long lowest = cfg.from;
  for (unsigned int i = 0; i < cfg.refs.size(); ++i)
    if (cfg.refs.at(i) > lowest)
      lowest = cfg.refs.at(i);

  // create -1 ordered config
  OBTetrahedralStereo::Config ordered = ts->GetConfig(lowest, OBStereo::Clockwise, OBStereo::ViewTowards);
  std::sort(ordered.refs.begin(), ordered.refs.end()); // increase clockwise -> -1
  //IdsToSymClasses(mol, ordered, canon_order); // FIXME

  if (cfg == ordered)
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

std::string canonicalSmiles(OBMol &mol_orig)
{
  OBMol mol;
  // read a smiles string
  OBConversion conv;
  conv.SetInFormat("smi");
  conv.SetOutFormat("can");

  std::vector<unsigned int> symmetry_classes, canon_order;


  // FIXME : why is this needed??
  std::string smiles = conv.WriteString(&mol_orig); 
  conv.ReadString(&mol, smiles);
  /*
  mol = mol_orig;
  std::vector<OBAtom*> atoms(mol.NumAtoms());
  FOR_ATOMS_OF_MOL(atom, mol)
    atoms[ canon_order.at(atom->GetIndex())-1 ] = &*atom;
  mol.RenumberAtoms(atoms);
  gs.GetSymmetry(symmetry_classes);
  //gs.CanonicalLabels(canon_order);
  */
  OBGraphSym gs(&mol);
  gs.GetSymmetry(symmetry_classes);
  gs.CanonicalLabels(canon_order);


  OBStereoisomer isomers(&mol);
  OBStereoFacade stereoFacade(&mol);
  // print number of stereomers
  //cout << "Enantiomer pairs: " << isomers.numEnantiomerPairs() << endl;
  //cout << "Diastereomers: " << isomers.numDiastereomers() << endl;

  // print the parities for the molecule
  OBStereoisomer::ParityVec parities, lastParities, canParities;
  std::vector<OBTetrahedralStereo::Config> configs;
  FOR_ATOMS_OF_MOL (atom, mol) {
    if (stereoFacade.HasTetrahedralStereo(atom->GetId())) {
      OBTetrahedralStereo *ts = stereoFacade.GetTetrahedralStereo(atom->GetId());
      OBTetrahedralStereo::Config config = ts->GetConfig();
      int canParity = configParityPaper(ts, &mol, canon_order);
      canParities.push_back(canParity);
      configs.push_back(config);
    }
  }

    
    
  std::vector<std::string> canonicalCandidates;
  lastParities = canParities;
  //
  // Handle enantiomers
  //
  const std::vector<OBStereoisomer::Enantiomer> &enantiomers = isomers.enantiomers();
  for (unsigned int i = 0; i < enantiomers.size(); ++i) {
    std::vector<std::string> candidates;
    //cout << "  enantiomer " << i+1 << endl;
    //cout << "    parities: ";
    bool foundEnantiomer = false;
    for (unsigned int j = 0; j < enantiomers.at(i).parities.size(); ++j) {
      const OBStereoisomer::ParityVec &pv = enantiomers.at(i).parities.at(j);
      //cout << "    ";
      for (unsigned int k = 0; k < pv.size(); ++k) {
        if (lastParities.at(k) != pv.at(k)) {
          permutateConfig(configs, k);
        }
        //cout << pv.at(k) << " ";
      }
      //cout << endl;
      if (pv == canParities) {
        //cout << "  ----> found enantiomer matching input structure: " << conv.WriteString(&mol);
        foundEnantiomer = true;
      }
      lastParities = pv;
      storeConfigs(&mol, configs, stereoFacade);
      std::string candidate = conv.WriteString(&mol); 
      //cout << "CANDIDATE: " << candidate;
      candidates.push_back(candidate);
    }

    if (foundEnantiomer)   
      canonicalCandidates = candidates;

    foundEnantiomer = false;
    candidates.clear();
    //cout << "    inverseParities: ";
    for (unsigned int j = 0; j < enantiomers.at(i).inverseParities.size(); ++j) {
      const OBStereoisomer::ParityVec &pv = enantiomers.at(i).inverseParities.at(j);
      //cout << "    ";
      for (unsigned int k = 0; k < pv.size(); ++k) {
        if (lastParities.at(k) != pv.at(k)) {
          permutateConfig(configs, k);
        }
        //cout << pv.at(k) << " ";
      }
      //cout << endl;
      if (pv == canParities) {
        //cout << "  ----> found enantiomer matching input structure: " << conv.WriteString(&mol);
        foundEnantiomer = true;
      }
      lastParities = pv;
      storeConfigs(&mol, configs, stereoFacade);
      std::string candidate = conv.WriteString(&mol); 
      //cout << "CANDIDATE: " << candidate;
      candidates.push_back(candidate);
    }

    if (foundEnantiomer)   
      canonicalCandidates = candidates;

  }

  //
  // Handle diastereomers
  //
  const std::vector<OBStereoisomer::Diastereomer> &diastereomers = isomers.diastereomers();
  for (unsigned int i = 0; i < diastereomers.size(); ++i) {
    std::vector<std::string> candidates;
    //cout << "  diastereomer " << i+1 << endl;
    //cout << "    parities: ";
    bool foundDiastereomer = false;
    for (unsigned int j = 0; j < diastereomers.at(i).parities.size(); ++j) {
      const OBStereoisomer::ParityVec &pv = diastereomers.at(i).parities.at(j);
      //cout << "    ";
      for (unsigned int k = 0; k < pv.size(); ++k) {
        if (lastParities.at(k) != pv.at(k)) {
          permutateConfig(configs, k);
        }
        //cout << pv.at(k) << " ";
      }
      //cout << endl;
      if (pv == canParities) {
        //cout << "  ----> found diastereomer matching input structure: " << conv.WriteString(&mol);
        foundDiastereomer = true;
      }
      lastParities = pv;
      storeConfigs(&mol, configs, stereoFacade);
      std::string candidate = conv.WriteString(&mol); 
      //cout << "CANDIDATE: " << candidate;
      candidates.push_back(candidate);
    }

    if (foundDiastereomer)   
      canonicalCandidates = candidates;
  }

  /*
  cout << "Canonical candidates:" << endl;
  for (unsigned int i = 0; i < canonicalCandidates.size(); ++i) {
    cout << "  " << canonicalCandidates.at(i);
  }
  */

  std::sort(canonicalCandidates.begin(), canonicalCandidates.end());
  //cout << "  True canonical SMILES: ";
  //cout << canonicalCandidates.front() << endl;

  return canonicalCandidates.front();
}




