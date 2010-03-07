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


  int getStereoIndex(OBAtom *atom, const Permutation &p, const std::vector<unsigned int> &symmetry_classes)
  {
    // construct map: symmetry class -> vector<id>
    std::map<unsigned int, std::vector<unsigned int> > symClass2index;
    OBBondIterator bi;
    for (OBAtom *nbr = atom->BeginNbrAtom(bi); nbr; nbr = atom->NextNbrAtom(bi)) {
      //std::cout << "    nbr symClass: " << symmetry_classes.at(nbr->GetIndex()) << std::endl;
      // create the vector if the symmetry class doesn't have one yet
      if (symClass2index.find(symmetry_classes.at(nbr->GetIndex())) == symClass2index.end())
        symClass2index[symmetry_classes.at(nbr->GetIndex())] = std::vector<unsigned int>();
      // add this nbr's id
      symClass2index[symmetry_classes.at(nbr->GetIndex())].push_back(nbr->GetIndex());
    }

    unsigned int numInversions = 0;
    std::map<unsigned int, std::vector<unsigned int> >::iterator ids;
    for (ids = symClass2index.begin(); ids != symClass2index.end(); ++ids) {
      // make sure the nbr ids are ordered as in the permutation
      std::vector<unsigned long> ordered;
      //cout << "ordered:";
      for (unsigned int pi = 0; pi < p.map.size(); ++pi) {
        if (std::find(ids->second.begin(), ids->second.end(), p.map.at(pi) - 1) != ids->second.end()) {
          ordered.push_back(p.map.at(pi));
          //cout << " " << p.map.at(pi);
        }
      }
      //cout << endl;

      // debug
      //for (unsigned int deb = 0; deb < ids->second.size(); ++deb)
      //  std::cout << "  " << ids->second[deb];
      //std::cout << std::endl;
      numInversions += OBStereo::NumInversions(ordered);
    }

    if (numInversions % 2)
      return -1; // odd # of permutations
    else
      return 1; // even # of permutations
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
  //IdsToSymClasses(mol, cfg, canon_order);

  // lowest priority = highest canonical
  // H = lowest priority
  unsigned long lowest = cfg.from;
  for (unsigned int i = 0; i < cfg.refs.size(); ++i)
    if (cfg.refs.at(i) > lowest)
      lowest = cfg.refs.at(i);

  // create -1 ordered config
  OBTetrahedralStereo::Config ordered = ts->GetConfig(lowest, OBStereo::Clockwise, OBStereo::ViewTowards);
  std::sort(ordered.refs.begin(), ordered.refs.end()); // increase clockwise -> -1
  //IdsToSymClasses(mol, ordered, canon_order);

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
      struct Holder {
        int parity;
        unsigned long canId;
      };
      struct compare {
        bool operator()(const Holder &a, const Holder &b) const { return a.canId > b.canId; }
      };
 
std::string canonicalSmiles(const std::string &smiles, std::vector<std::string> &out_candidates, 
    const std::string &ref = std::string())
{
  // read a smiles string
  OBMol mol;
  OBConversion conv, smiConv;
  OB_REQUIRE( conv.SetInFormat("smi") );
  OB_REQUIRE( conv.SetOutFormat("can") );
  // read a smiles string
  OB_REQUIRE( conv.ReadString(&mol, smiles) );

    std::vector<unsigned int> symmetry_classes, canon_order;

    std::string canOrdered = conv.WriteString(&mol);
    conv.ReadString(&mol, canOrdered);

    OBGraphSym gs2(&mol);
    gs2.GetSymmetry(symmetry_classes);
    gs2.CanonicalLabels(canon_order);
  
    cout << "SYMCLASSES: ";
    for (unsigned int i = 0; i < symmetry_classes.size(); ++i)
      cout << symmetry_classes.at(i) << " ";
    cout << endl;


   /* 
    std::vector<OBAtom*> atoms(mol.NumAtoms());
    for (unsigned int i = 0; i < atoms.size(); ++i)
      atoms[canon_order.at(i)-1] = mol.GetAtom(i+1);
    mol.RenumberAtoms(atoms);
    
    OBGraphSym gs(&mol);
    gs.GetSymmetry(symmetry_classes);
    gs.CanonicalLabels(canon_order);
 */
    cout << "CANORDER: ";
    for (unsigned int i = 0; i < canon_order.size(); ++i)
      cout << canon_order.at(i) << " ";
    cout << endl;

    cout << "IDORDER: ";
    FOR_ATOMS_OF_MOL (atom, mol)
      cout << atom->GetId() << " ";
    cout << endl;


    OBStereoisomer isomers(&mol);
    OBStereoFacade stereoFacade(&mol);
    // print number of stereomers
    cout << "Enantiomer pairs: " << isomers.numEnantiomerPairs() << endl;
    cout << "Diastereomers: " << isomers.numDiastereomers() << endl;


    std::vector<unsigned long> canIds;
    cout << "XXX parities2: ";
    OBStereoisomer::ParityVec parities2;
    FOR_ATOMS_OF_MOL (atom, mol) {
      if (stereoFacade.HasTetrahedralStereo(atom->GetId())) {
        int index = getStereoIndex(&*atom, canon_order, symmetry_classes);
        //int index = getStereoIndex(&*atom, canon_order, canon_order);
        cout << index << " ";
        parities2.push_back(index);
      }
    }
    cout << endl;
 
    // print the parities for the molecule
    cout << "XXX parities sym: ";
    OBStereoisomer::ParityVec parities, lastParities, canParities, idParities, parities3;
    std::vector<OBTetrahedralStereo::Config> configs;
    std::vector<unsigned long> atomIds;
    FOR_ATOMS_OF_MOL (atom, mol) {
      if (stereoFacade.HasTetrahedralStereo(atom->GetId())) {
        OBTetrahedralStereo *ts = stereoFacade.GetTetrahedralStereo(atom->GetId());
        parities3.push_back(configParityPaper(ts, &mol, canon_order));
        OBTetrahedralStereo::Config config = stereoFacade.GetTetrahedralStereo(atom->GetId())->GetConfig();
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

    /*
    if (OBStereo::NumInversions(canIds)) {
      std::vector<Holder> holders;
      for (unsigned int i = 0; i < canIds.size(); ++i) {
        Holder holder;
        holder.canId = canIds.at(i);
        holder.parity = canParities.at(i);
        holders.push_back(holder);
      }
      std::sort(holders.begin(), holders.end(), compare());
      canParities.clear();
      for (unsigned int i = 0; i < canIds.size(); ++i)
        canParities.push_back(holders.at(i).parity);
    }
*/
    std::vector<std::string> canonicalCandidates;


    lastParities = canParities;
    // print all stereomers
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


      /*
      if (!ref.empty())
        if (std::find(candidates.begin(), candidates.end(), ref) != candidates.end())
          canonicalCandidates = candidates;
      else
      */
        if (foundDiastereomer)   
          canonicalCandidates = candidates;
        /*
      if (std::find(candidates.begin(), candidates.end(), ref) != candidates.end()) {
        cout << "CANDIDATE FOUND : " << candidates.size() << endl;
               canonicalCandidates = candidates;
      }
      */
 
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

  int N = 200;
  testCount++;

  std::vector<OBAtom*> atoms;
  FOR_ATOMS_OF_MOL(atom, mol)
    atoms.push_back(&*atom);
  
  std::vector< std::vector<std::string> > allCandidates; 

  std::vector<std::string> candidates;
  std::string ref = canonicalSmiles(smiles, candidates);
  cout << "ref = " << ref;
 
  bool result = true;
  for (int i = 0; i < N; ++i) {
    // shuffle the atoms
    std::random_shuffle(atoms.begin(), atoms.end());
    mol.RenumberAtoms(atoms);
    // get can smiles
    std::string smiles = smiConv.WriteString(&mol);
    std::string cansmi = canonicalSmiles(smiles, candidates, ref);
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
    OB_ASSERT( allCandidates.at(i) == candidates2 );
  }

  return result;
}

bool doShuffleTestFile(const std::string &filename)
{
  cout << " True Shuffling: " << filename << endl;
  std::string file = GetFilename(filename);
  // read a smiles string
  OBMol mol;
  OBConversion canConv;
  OBFormat *format = canConv.FormatFromExt(file.c_str());
  OB_REQUIRE( format );
  OB_REQUIRE( canConv.SetInFormat(format) );
  OB_REQUIRE( canConv.ReadFile(&mol, file) );
  OB_REQUIRE( canConv.SetOutFormat("can") );

  std::string smiles = canConv.WriteString(&mol);
  int N = 200;
  testCount++;

  std::vector<OBAtom*> atoms;
  FOR_ATOMS_OF_MOL(atom, mol)
    atoms.push_back(&*atom);
  
  std::vector< std::vector<std::string> > allCandidates; 

  std::vector<std::string> candidates;
  std::string ref = canonicalSmiles(smiles, candidates);
  cout << "ref = " << ref;
 
  bool result = true;
  for (int i = 0; i < N; ++i) {
    // shuffle the atoms
    std::random_shuffle(atoms.begin(), atoms.end());
    mol.RenumberAtoms(atoms);
    // get can smiles
    std::string smiles = canConv.WriteString(&mol);
    std::string cansmi = canonicalSmiles(smiles, candidates, ref);
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
      cout << " " << cansmi;
      failed++;
      return false;
    }
  }

  return true;
}
*/

int main(int argc, char **argv)
{
  //OB_ASSERT( doShuffleTest("O[C@H]1CC[C@@H](O)CC1") );
  //OB_ASSERT( doShuffleTest("O[C@H]1C[C@@H](O)C[C@H](O)C1") );
  //OB_ASSERT( doShuffleTest("O[C@H]1C[C@@H](O)C[C@@H](O)C1") );


  OB_ASSERT( doShuffleTestFile("stereo/cyclohexanediol_D1.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/cyclohexanediol_D2.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/cyclohexanetriol_D1.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/cyclohexanetriol_D2.mol") );


  // These work for mol files, not for smiles?????
  OB_ASSERT( doShuffleTestFile("stereo/cyclobutane_D1.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/cyclobutane_D2.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/cyclobutane_D3.mol") );
  OB_ASSERT( doShuffleTestFile("stereo/cyclobutane_D4.mol") );
  //OB_ASSERT( doShuffleTest("[C@@H]1([C@H]([C@H]([C@H]1C)C)C)C") );
//  OB_ASSERT( doShuffleTest("[C@@H]1([C@@H]([C@H]([C@H]1C)C)C)C") );	
  //OB_ASSERT( doShuffleTest("[C@@H]1([C@@H]([C@H]([C@@H]1C)C)C)C") );	
  //OB_ASSERT( doShuffleTest("[C@@H]1([C@@H]([C@@H]([C@H]1C)C)C)C") );	



  //OB_ASSERT( doShuffleTest("O[C@H]1[C@@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H]1O") );
  
  //OB_ASSERT( doShuffleTest("") );

  cout << "PASSED TESTS: " << testCount - failed << "/" << testCount << endl;

  return 0;
}

