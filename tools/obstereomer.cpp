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
  //cout << "ERROR  " << cfg << endl;
  std::vector<unsigned long> refs = cfg.refs;
  refs.insert(refs.begin(), config.from); // doesn't matter if view is from or towards, will be sorted anyway
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

    std::vector<unsigned int> symmetry_classes;
    OBGraphSym gs(&mol);
    gs.GetSymmetry(symmetry_classes, false);

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

    cout << "Canonical candidates:" << endl;
    for (unsigned int i = 0; i < canonicalCandidates.size(); ++i) {
      cout << "ERROR  " << canonicalCandidates.at(i);
    }

    std::sort(canonicalCandidates.begin(), canonicalCandidates.end());
    cout << "ERROR  True canonical SMILES: ";
    cout << canonicalCandidates.front() << endl;

      
  } // end for loop
  
  return(1);
}


