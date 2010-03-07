/**********************************************************************
obgraphgsym - Open Babel properties calculation

Copyright (C) 2010 Tim Vandermeersch
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/


// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif

#include <openbabel/babelconfig.h>

#include <openbabel/mol.h>
#include <openbabel/descriptor.h>
#include <openbabel/obconversion.h>
#include <openbabel/graphsym.h>
#include <openbabel/stereo/stereoisomer.h>
#include <openbabel/stereo/tetrahedral.h>
#ifndef _MSC_VER
  #include <unistd.h>
#endif

using namespace std;
using namespace OpenBabel;
void some_tests(OBMol &mol);
// PROTOTYPES /////////////////////////////////////////////////////////////////
int nrings(OBMol &mol);
string sequence(OBMol &mol);

///////////////////////////////////////////////////////////////////////////////
//! \brief Compute some properties easy to access from open babel
//
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



int main(int argc,char **argv)
{
  char *program_name= argv[0];
  int c;
  char *FileIn = NULL;

  if (argc != 2)
    {
      string err = "Usage: ";
      err += program_name;
      err += " <filename>\n";
      cout << err;
      exit(-1);
    }
  else
    {
      FileIn  = argv[1];
    }

  // Find Input filetype
  OBConversion conv, smileConv, smiConv;
  OBFormat *format = conv.FormatFromExt(FileIn);
  smileConv.SetOutFormat("can");
  smiConv.SetOutFormat("smi");
    
  if (!format || !conv.SetInFormat(format))
    {
      cerr << program_name << ": cannot read input format!" << endl;
      exit (-1);
    }

  ifstream ifs;

  // Read the file
  ifs.open(FileIn);
  if (!ifs)
    {
      cerr << program_name << ": cannot read input file!" << endl;
      exit (-1);
    }
  
  OBMol mol;
  OBFormat *canSMIFormat = conv.FindFormat("can");
  OBFormat *inchiFormat = conv.FindFormat("inchi");

  for (c = 1;; ++c) {
    mol.Clear();
    conv.Read(&mol, &ifs);
    if (mol.Empty())
      break;

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


    cout << "ERROR XXX parities2: ";
    OBStereoisomer::ParityVec parities2;
    FOR_ATOMS_OF_MOL (atom, mol) {
      if (stereoFacade.HasTetrahedralStereo(atom->GetId())) {
        int index = getStereoIndex(&*atom, canon_order, symmetry_classes);
        cout << index << " ";
        parities2.push_back(index);
      }
    }
    cout << endl;
 
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

    lastParities = parities2;
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
        if (pv == parities2) {
          cout << "ERROR  ----> found diastereomer matching input structure: " << smiConv.WriteString(&mol);
          foundDiastereomer = true;
        }
        lastParities = pv;
        storeConfigs(&mol, configs, stereoFacade);
        std::string candidate = smileConv.WriteString(&mol); 
        cout << candidate;
        candidates.push_back(candidate);
      }

      if (foundDiastereomer)   
        canonicalCandidates = candidates;
    }
    cout << "ERROR XXX: " << canonicalCandidates.size() << "   Inv = " << OBStereo::NumInversions(atomIds) << endl;

    cout << "Canonical candidates: " << canonicalCandidates.size() << endl;
    for (unsigned int i = 0; i < canonicalCandidates.size(); ++i) {
      cout << "ERROR  " << canonicalCandidates.at(i);
    }

    std::sort(canonicalCandidates.begin(), canonicalCandidates.end());
    cout << "ERROR  True canonical SMILES: ";
    cout << canonicalCandidates.front() << endl;

      
  } // end for loop
  
  return(1);
}


