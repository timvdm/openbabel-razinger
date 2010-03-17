#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>

#include <openbabel/graphsym.h>
#include <openbabel/stereo/tetrahedral.h>
#include <openbabel/stereo/cistrans.h>
#include <openbabel/stereo/stereoisomer.h>

#include <iostream>
#include <vector>
#include <algorithm>


using namespace std;
namespace OpenBabel
{

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

void IdsToSymClasses(OBMol *mol, OBCisTransStereo::Config &config, 
    const std::vector<unsigned int> &symClasses)
{
  OBAtom *atom;
  // begin
  atom = mol->GetAtomById(config.begin);
  if (atom) {
    if (atom->IsHydrogen())
      config.begin = OBStereo::ImplicitRef;
    else
      config.begin = symClasses.at(atom->GetIndex());
  }
  //end 
  atom = mol->GetAtomById(config.end);
  if (atom) {
    if (atom->IsHydrogen())
      config.end = OBStereo::ImplicitRef;
    else
      config.end = symClasses.at(atom->GetIndex());
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


int configParity(OBTetrahedralStereo *ts, OBMol *mol, const std::vector<unsigned int> &canon_order)
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

int configParity(OBCisTransStereo *ct, OBMol *mol, const std::vector<unsigned int> &canon_order)
{
  OBCisTransStereo::Config cfg = ct->GetConfig();
  IdsToSymClasses(mol, cfg, canon_order); // FIXME

  // lowest priority = highest canonical
  // H = lowest priority
  unsigned long lowest = cfg.refs[0];
  for (unsigned int i = 0; i < cfg.refs.size(); ++i)
    if (cfg.refs.at(i) > lowest)
      lowest = cfg.refs.at(i);

  // create -1 ordered config
  OBCisTransStereo::Config ordered = ct->GetConfig(lowest, OBStereo::ShapeU);
  std::sort(ordered.refs.begin(), ordered.refs.end()); // increase clockwise -> -1
  //IdsToSymClasses(mol, ordered, canon_order); // FIXME

  if (cfg == ordered)
    return 1;
  else
    return -1;
}

void permutateConfig(std::vector<OBTetrahedralStereo::Config> &configs, unsigned int k)
{
  OBTetrahedralStereo::Config &config = configs[k];
  OBStereo::Permutate(config.refs, 0, 1);
}
 
void permutateConfig(std::vector<OBCisTransStereo::Config> &configs, unsigned int k)
{
  OBCisTransStereo::Config &config = configs[k];
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

void storeConfigs(OBMol *mol, const std::vector<OBCisTransStereo::Config> &configs, OBStereoFacade &stereoFacade)
{
  unsigned int idx = 0;
  FOR_BONDS_OF_MOL (bond, mol) {
    if (stereoFacade.HasCisTransStereo(bond->GetId())) {
      OBCisTransStereo *ct = stereoFacade.GetCisTransStereo(bond->GetId());
      ct->SetConfig(configs.at(idx));
      ++idx;
    }
  }
}

std::string removeNewLine(const std::string &str)
{
  return str.substr(0, str.size()-1);
}

std::string canonicalSmiles(OBMol &mol_orig, std::vector<std::string> &stereoisomers)
{
  OBMol mol;
  // read a smiles string
  OBConversion conv;
  conv.SetInFormat("smi");
  conv.SetOutFormat("can");

  std::vector<unsigned int> symmetry_classes, canon_order;

  std::string smiles = conv.WriteString(&mol_orig); 
  conv.ReadString(&mol, smiles);
  
  OBStereoFacade stereoFacade(&mol);
  OBStereoisomer isomers(&mol);
  
  OBGraphSym gs(&mol);
  gs.GetSymmetry(symmetry_classes);
  gs.CanonicalLabels(canon_order);


  // print number of stereomers
  //cout << "Enantiomer pairs: " << isomers.numEnantiomerPairs() << endl;
  //cout << "Diastereomers: " << isomers.numDiastereomers() << endl;

  // print the parities for the molecule
  OBStereoisomer::ParityVec parities, lastParities, canParities;
  std::vector<OBTetrahedralStereo::Config> tetrahedralConfigs;
  std::vector<OBCisTransStereo::Config> cistransConfigs;
  FOR_ATOMS_OF_MOL (atom, mol) {
    if (stereoFacade.HasTetrahedralStereo(atom->GetId())) {
      OBTetrahedralStereo *ts = stereoFacade.GetTetrahedralStereo(atom->GetId());
      OBTetrahedralStereo::Config config = ts->GetConfig();
      int canParity = configParity(ts, &mol, canon_order);
      canParities.push_back(canParity);
      tetrahedralConfigs.push_back(config);
    }
  }
  FOR_BONDS_OF_MOL (bond, mol) {
    if (stereoFacade.HasCisTransStereo(bond->GetId())) {
      OBCisTransStereo *ct = stereoFacade.GetCisTransStereo(bond->GetId());
      OBCisTransStereo::Config config = ct->GetConfig();
      int canParity = configParity(ct, &mol, canon_order);
      canParities.push_back(canParity);
      cistransConfigs.push_back(config);
    }
  }

  int numTetrahedral = tetrahedralConfigs.size();
  int numCisTrans = cistransConfigs.size();
 

    
    
  std::vector<std::string> canonicalCandidates;
  lastParities = canParities;
  //
  // Handle enantiomers
  //
  const std::vector<OBStereoisomer::Enantiomer> &enantiomers = isomers.GetEnantiomers();
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
          if (k < numTetrahedral)
            permutateConfig(tetrahedralConfigs, k);
          else
            permutateConfig(cistransConfigs, k - numTetrahedral);
        }
        //cout << pv.at(k) << " ";
      }
      //cout << endl;
      if (pv == canParities) {
        //cout << "  ----> found enantiomer matching input structure: " << conv.WriteString(&mol);
        foundEnantiomer = true;
      }
      lastParities = pv;
      storeConfigs(&mol, tetrahedralConfigs, stereoFacade);
      storeConfigs(&mol, cistransConfigs, stereoFacade);
      std::string candidate = conv.WriteString(&mol); 
      //cout << "CANDIDATE: " << candidate;
      candidates.push_back(candidate);
    }

    std::sort(candidates.begin(), candidates.end());

    std::stringstream ss;
    ss << removeNewLine(candidates.front());
    ss << " Enantiomer +" << i+1;
    stereoisomers.push_back(ss.str());

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
          if (k < numTetrahedral)
            permutateConfig(tetrahedralConfigs, k);
          else
            permutateConfig(cistransConfigs, k - numTetrahedral);
        }
        //cout << pv.at(k) << " ";
      }
      //cout << endl;
      if (pv == canParities) {
        //cout << "  ----> found enantiomer matching input structure: " << conv.WriteString(&mol);
        foundEnantiomer = true;
      }
      lastParities = pv;
      storeConfigs(&mol, tetrahedralConfigs, stereoFacade);
      storeConfigs(&mol, cistransConfigs, stereoFacade);
      std::string candidate = conv.WriteString(&mol); 
      //cout << "CANDIDATE: " << candidate;
      candidates.push_back(candidate);
    }

    std::sort(candidates.begin(), candidates.end());

    ss.str("");
    ss << removeNewLine(candidates.front());
    ss << " Enantiomer -" << i+1;
    stereoisomers.push_back(ss.str());
    
    if (foundEnantiomer)   
      canonicalCandidates = candidates;

  }

  //
  // Handle diastereomers
  //
  const std::vector<OBStereoisomer::Diastereomer> &diastereomers = isomers.GetDiastereomers();
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
          if (k < numTetrahedral)
            permutateConfig(tetrahedralConfigs, k);
          else
            permutateConfig(cistransConfigs, k - numTetrahedral);
        }
        //cout << pv.at(k) << " ";
      }
      //cout << endl;
      if (pv == canParities) {
        //cout << "  ----> found diastereomer matching input structure: " << conv.WriteString(&mol);
        foundDiastereomer = true;
      }
      lastParities = pv;
      storeConfigs(&mol, tetrahedralConfigs, stereoFacade);
      storeConfigs(&mol, cistransConfigs, stereoFacade);
      std::string candidate = conv.WriteString(&mol); 
      //cout << "CANDIDATE: " << candidate;
      candidates.push_back(candidate);
    }

    std::sort(candidates.begin(), candidates.end());

    std::stringstream ss;
    ss << removeNewLine(candidates.front());
    ss << " Diastereomer " << i+1;
    stereoisomers.push_back(ss.str());
    
    if (foundDiastereomer)   
      canonicalCandidates = candidates;
  }

  /*
  cout << "Canonical candidates:" << endl;
  for (unsigned int i = 0; i < canonicalCandidates.size(); ++i) {
    cout << "  " << canonicalCandidates.at(i);
  }
  */

  //cout << "  True canonical SMILES: ";
  //cout << canonicalCandidates.front() << endl;

  if (!canonicalCandidates.empty())
    return canonicalCandidates.front();
  return removeNewLine(conv.WriteString(&mol)); 
}








class Can2Format : public OBMoleculeFormat
{
  public:
    Can2Format()
    {
      OBConversion::RegisterFormat("Can2",this);
      OBConversion::RegisterOptionParam("f", this, 1);
      OBConversion::RegisterOptionParam("n", this);
      OBConversion::RegisterOptionParam("s", this, 0, OBConversion::INOPTIONS);
    }

    virtual const char* Description() //required
    {
      return "";
    }

    virtual const char* SpecificationURL()
    {
      return "";
    }

    virtual const char* GetMIMEType() 
    { return "chemical/x-xxx"; }


    virtual unsigned int Flags()
    {
      return NOTREADABLE;
    };

    virtual int SkipObjects(int n, OBConversion* pConv)
    {
      return 0;
    }

    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
  private:
};  

Can2Format theCan2Format;

////////////////////////////////////////////////////////////////

bool Can2Format::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
  OBMol* pmol = dynamic_cast<OBMol*>(pOb);
  if(pmol==NULL)
      return false;

  // Make a copy so the molecule doesn't change
  OBMol mol(*pmol);

  ostream& ofs = *pConv->GetOutStream();

  std::vector<std::string> stereoisomers;
  if (pConv->IsOption("s")) {
    OBStereoFacade stereoFacade(&mol);

    // Make sure all stereocenters have specified objects...
    OBGraphSym graphSym(&mol);
    std::vector<unsigned int> symmetry_classes;
    graphSym.GetSymmetry(symmetry_classes, false);

    vector<StereogenicUnit> stereoUnits = FindStereogenicUnits(&mol, symmetry_classes); // FIXME need to cache this

    for (unsigned int i = 0; i < stereoUnits.size(); ++i) {
      const StereogenicUnit &unit = stereoUnits.at(i);
      if (unit.type == OBStereo::Tetrahedral) {
        if (stereoFacade.HasTetrahedralStereo(unit.id)) {
          OBTetrahedralStereo *ts = stereoFacade.GetTetrahedralStereo(unit.id);
          OBTetrahedralStereo::Config config = ts->GetConfig();
          config.specified = true;
          ts->SetConfig(config);
        } else {
          OBAtom *atom = mol.GetAtomById(unit.id);
          std::vector<unsigned long> refs;
          FOR_NBORS_OF_ATOM (nbr, atom)
            refs.push_back(nbr->GetId());
          if (refs.size() < 4)
            refs.push_back(OBStereo::ImplicitRef);

          OBTetrahedralStereo *ts = new OBTetrahedralStereo(&mol);
          OBTetrahedralStereo::Config config;
          config.center = unit.id;
          config.from = refs.back();
          refs.pop_back();
          config.refs = refs;
          ts->SetConfig(config);
          mol.SetData(ts);
        }
      }
      if (unit.type == OBStereo::CisTrans) {
        if (stereoFacade.HasCisTransStereo(unit.id)) {
          OBCisTransStereo *ct = stereoFacade.GetCisTransStereo(unit.id);
          OBCisTransStereo::Config config = ct->GetConfig();
          config.specified = true;
          ct->SetConfig(config);
        } else {
          OBBond *bond = mol.GetBondById(unit.id);
          OBAtom *begin = bond->GetBeginAtom();
          OBAtom *end = bond->GetEndAtom();
          std::vector<unsigned long> refs;
          FOR_NBORS_OF_ATOM (nbr, begin)
            if (nbr->GetId() != end->GetId())
              refs.push_back(nbr->GetId());
          FOR_NBORS_OF_ATOM (nbr, end)
            if (nbr->GetId() != begin->GetId())
              refs.push_back(nbr->GetId());
          if (refs.size() < 4)
            refs.push_back(OBStereo::ImplicitRef);

          OBCisTransStereo *ct = new OBCisTransStereo(&mol);
          OBCisTransStereo::Config config;
          config.begin = begin->GetId();
          config.end = end->GetId();
          config.refs = refs;
          ct->SetConfig(config);
          mol.SetData(ct);
        }
      }
 
    }

    std::string cansmiles = canonicalSmiles(mol, stereoisomers);
    
    for (unsigned int i = 0; i < stereoisomers.size(); ++i)
      ofs << stereoisomers.at(i) << std::endl;
  } else {
    std::string cansmiles = canonicalSmiles(mol, stereoisomers);
    ofs << cansmiles;
  }

  return true; //or false to stop converting
}

} //namespace OpenBabel

