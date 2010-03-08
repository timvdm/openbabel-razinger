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
  //cout << "  " << cfg << endl;
  std::vector<unsigned long> refs = cfg.refs;
  refs.insert(refs.begin(), cfg.from); // doesn't matter if view is from or towards, will be sorted anyway
  //std::sort(refs.begin(), refs.end());

  bool p = (OBStereo::NumInversions(refs) % 2) ? true : false;
        
  if (p)
    return 1;
  else
    return -1;
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

std::string canonicalSmiles(OBMol &mol_orig, std::vector<std::string> &out_candidates)
{
  OBMol mol;
  // read a smiles string
  OBConversion conv;
  OB_REQUIRE( conv.SetInFormat("smi") );
  OB_REQUIRE( conv.SetOutFormat("can") );

  std::vector<unsigned int> symmetry_classes, canon_order;

  // FIXME : why is this needed??
  std::string smiles = conv.WriteString(&mol_orig); 
  conv.ReadString(&mol, smiles);

  OBGraphSym gs(&mol);
  gs.GetSymmetry(symmetry_classes);
  gs.CanonicalLabels(canon_order);
 /*
  std::vector<OBAtom*> atoms(mol.NumAtoms());
  FOR_ATOMS_OF_MOL(atom, mol)
    atoms[canon_order.at(atom->GetIndex())-1] = &*atom;
    
  mol.RenumberAtoms(atoms);
   gs.GetSymmetry(symmetry_classes);
  gs.CanonicalLabels(canon_order);
*/ 
    /*
    cout << "SYMCLASSES: ";
    for (unsigned int i = 0; i < symmetry_classes.size(); ++i)
      cout << symmetry_classes.at(i) << " ";
    cout << endl;
    cout << "CANORDER: ";
    for (unsigned int i = 0; i < canon_order.size(); ++i)
      cout << canon_order.at(i) << " ";
    cout << endl;
    cout << "IDORDER: ";
    FOR_ATOMS_OF_MOL (atom, mol)
      cout << atom->GetId() << " ";
    cout << endl;
    */


  OBStereoisomer isomers(&mol);
  OBStereoFacade stereoFacade(&mol);
  // print number of stereomers
  cout << "Enantiomer pairs: " << isomers.numEnantiomerPairs() << endl;
  cout << "Diastereomers: " << isomers.numDiastereomers() << endl;


  std::vector<unsigned long> canIds;
 
    // print the parities for the molecule
    cout << "XXX parities sym: ";
    OBStereoisomer::ParityVec parities, lastParities, canParities, idParities;
    std::vector<OBTetrahedralStereo::Config> configs;
    std::vector<unsigned long> atomIds;
    FOR_ATOMS_OF_MOL (atom, mol) {
      if (stereoFacade.HasTetrahedralStereo(atom->GetId())) {
        OBTetrahedralStereo *ts = stereoFacade.GetTetrahedralStereo(atom->GetId());
        OBTetrahedralStereo::Config config = ts->GetConfig();
        atomIds.push_back(config.center);
        canIds.push_back(canon_order.at(atom->GetIndex()));
        int symParity = configParity(config, &mol, symmetry_classes);
        parities.push_back(symParity);
        cout << symParity << " ";
        //int canParity = configParity(config, &mol, canon_order);
        int canParity = configParityPaper(ts, &mol, canon_order);
        canParities.push_back(canParity);
        int idParity = configParity(config, &mol);
        idParities.push_back(idParity);
        configs.push_back(config);
      }
    }
    cout << endl;
    cout << "XXX parities can = ";
    for (unsigned int i = 0; i < canParities.size(); ++i)
      cout << canParities.at(i) << " ";
    cout << endl;
    cout << "XXX parities id = ";
    for (unsigned int i = 0; i < idParities.size(); ++i)
      cout << idParities.at(i) << " ";
    cout << endl;
 
    cout << "XXX canIds = ";
    for (unsigned int i = 0; i < canIds.size(); ++i)
      cout << canIds.at(i) << " ";
    cout << "    inversion = " << OBStereo::NumInversions(canIds) << endl;

    std::vector<std::string> canonicalCandidates;


    lastParities = canParities;
    //
    // Handle enantiomers
    //
    const std::vector<OBStereoisomer::Enantiomer> &enantiomers = isomers.enantiomers();
    for (unsigned int i = 0; i < enantiomers.size(); ++i) {
      std::vector<std::string> candidates;
      cout << "  enantiomer " << i+1 << endl;
      cout << "    parities: ";
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
        if (pv == canParities) {
          cout << "  ----> found enantiomer matching input structure: " << conv.WriteString(&mol);
          foundEnantiomer = true;
        }
        lastParities = pv;
        storeConfigs(&mol, configs, stereoFacade);
        std::string candidate = conv.WriteString(&mol); 
        cout << "CANDIDATE: " << candidate;
        candidates.push_back(candidate);
      }
        
      if (foundEnantiomer)   
        canonicalCandidates = candidates;

      foundEnantiomer = false;
      candidates.clear();
      cout << "    inverseParities: ";
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
        if (pv == canParities) {
          cout << "  ----> found enantiomer matching input structure: " << conv.WriteString(&mol);
          foundEnantiomer = true;
        }
        lastParities = pv;
        storeConfigs(&mol, configs, stereoFacade);
        std::string candidate = conv.WriteString(&mol); 
        cout << "CANDIDATE: " << candidate;
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
      cout << "  diastereomer " << i+1 << endl;
      cout << "    parities: ";
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
          cout << "  ----> found diastereomer matching input structure: " << conv.WriteString(&mol);
          foundDiastereomer = true;
        }
        lastParities = pv;
        storeConfigs(&mol, configs, stereoFacade);
        std::string candidate = conv.WriteString(&mol); 
        cout << "CANDIDATE: " << candidate;
        candidates.push_back(candidate);
      }
        
      if (foundDiastereomer)   
        canonicalCandidates = candidates;
    }
    cout << "XXX: " << canonicalCandidates.size() << "   Inv = " << OBStereo::NumInversions(canIds) << endl;


    cout << "Canonical candidates:" << endl;
    for (unsigned int i = 0; i < canonicalCandidates.size(); ++i) {
      cout << "  " << canonicalCandidates.at(i);
    }

    std::sort(canonicalCandidates.begin(), canonicalCandidates.end());
    cout << "  True canonical SMILES: ";
    cout << canonicalCandidates.front() << endl;

    out_candidates = canonicalCandidates;

    return canonicalCandidates.front();
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

  int N = 50;
  testCount++;

  std::vector<OBAtom*> atoms;
  FOR_ATOMS_OF_MOL(atom, mol)
    atoms.push_back(&*atom);
  
  std::vector< std::vector<std::string> > allCandidates; 

  std::vector<std::string> candidates;
  std::string ref = canonicalSmiles(mol, candidates); // FIXME
  cout << "ref = " << ref;
 
  bool result = true;
  for (int i = 0; i < N; ++i) {
    // shuffle the atoms
    std::random_shuffle(atoms.begin(), atoms.end());
    mol.RenumberAtoms(atoms);
    // get can smiles
    std::string cansmi = canonicalSmiles(mol, candidates); // FIXME
    allCandidates.push_back(candidates);
    OB_ASSERT( cansmi == ref );
    // comapare with ref
    if (cansmi != ref) {
      cout << " " << cansmi;
      if (result)
        failed++;
      result = false;
    }
  }

  unsigned int numCandidates = allCandidates.at(0).size();
  const std::vector<std::string> &candidates2 = allCandidates.at(0);
  for (unsigned int i = 0; i < allCandidates.size(); ++i) {
    cout << allCandidates.at(i).at(0);
    OB_ASSERT( allCandidates.at(i).size() == numCandidates );
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
  
  std::vector< std::vector<std::string> > allCandidates; 

  std::vector<std::string> candidates;
  std::string ref = canonicalSmiles(mol, candidates); // FIXME
  cout << "ref = " << ref;
 
  bool result = true;
  for (int i = 0; i < N; ++i) {
    // shuffle the atoms
    std::random_shuffle(atoms.begin(), atoms.end());
    mol.RenumberAtoms(atoms);
    // get can smiles
    std::string cansmi = canonicalSmiles(mol, candidates); // FIXME
    allCandidates.push_back(candidates);
    OB_ASSERT( cansmi == ref );
    // comapare with ref
    if (cansmi != ref) {
      cout << " " << cansmi;
      if (result)
        failed++;
      result = false;
    }
  }

  unsigned int numCandidates = allCandidates.at(0).size();
  const std::vector<std::string> &candidates2 = allCandidates.at(0);
  for (unsigned int i = 0; i < allCandidates.size(); ++i) {
    cout << allCandidates.at(i).at(0);
    OB_ASSERT( allCandidates.at(i).size() == numCandidates );
  }

  return result;
}

int main(int argc, char **argv)
{
  if (argc == 2) {
    OB_ASSERT( doShuffleTestFile(argv[1]) );
    return 0;
  }

  OB_ASSERT( doShuffleTest("O[C@H]1CC[C@@H](O)CC1") );
  OB_ASSERT( doShuffleTest("O[C@H]1C[C@@H](O)C[C@H](O)C1") );
  OB_ASSERT( doShuffleTest("O[C@H]1C[C@@H](O)C[C@@H](O)C1") );
  
  OB_ASSERT( doShuffleTest("[C@@H]1([C@H]([C@H]([C@H]1C)C)C)C") );
  OB_ASSERT( doShuffleTestFile("stereo/cyclobutane_D1.smi") );

  
  
  // 
  // Enantiomers only
  //
  OB_ASSERT( doShuffleTestFile("stereo/razinger_fig7_5_spec.mol") );
  
 
  //
  // Diastereomers only
  //
  OB_ASSERT( doShuffleTestFile("stereo/cyclohexanediol_D1.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/cyclohexanediol_D2.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/cyclohexanetriol_D1.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/cyclohexanetriol_D2.mol") );
  // These work for mol files, not for smiles????? See graphsymtest.cpp
  OB_ASSERT( doShuffleTestFile("stereo/cyclobutane_D1.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/cyclobutane_D2.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/cyclobutane_D3.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/cyclobutane_D4.mol") );
  //OB_ASSERT( doShuffleTest("[C@@H]1([C@H]([C@H]([C@H]1C)C)C)C") );
//  OB_ASSERT( doShuffleTest("[C@@H]1([C@@H]([C@H]([C@H]1C)C)C)C") );	
  //OB_ASSERT( doShuffleTest("[C@@H]1([C@@H]([C@H]([C@@H]1C)C)C)C") );	
  //OB_ASSERT( doShuffleTest("[C@@H]1([C@@H]([C@@H]([C@H]1C)C)C)C") );	


  // Mixed: enantiomers + diastereomers
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

  //OB_ASSERT( doShuffleTest("O[C@H]1[C@@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H]1O") );
  
  //OB_ASSERT( doShuffleTest("") );

  cout << "PASSED TESTS: " << testCount - failed << "/" << testCount << endl;

  return 0;
}

