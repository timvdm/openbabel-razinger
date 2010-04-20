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
#include <openbabel/permutation.h>
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
  OBConversion conv;
  OBFormat *format = conv.FormatFromExt(FileIn);
    
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

    OBGraphSym sym(&mol);
    std::vector<unsigned int> symmetry_classes;
    sym.GetSymmetry(symmetry_classes, false);
  
    OBPermutationGroup G = FindAutomorphisms(&mol, symmetry_classes);


    for (unsigned int i = 0; i < G.Size(); ++i) {
      cout << "  p" << i+1 << " = [";
      for (unsigned int j = 0; j < G.At(i).map.size(); ++j) {
        if (j > 0)
          cout << ",";
        cout << G.At(i).map[j];
      }
      cout << "]" << endl;

    }
      
  } // end for loop
  
  return(1);
}


