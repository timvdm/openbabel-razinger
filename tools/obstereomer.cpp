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
  OBConversion conv, can2Conv;
  OBFormat *format = conv.FormatFromExt(FileIn);
  can2Conv.SetOutFormat("can2");
  can2Conv.AddOption("s");
    
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
  OBFormat *can2Format = conv.FindFormat("can2");
  OBFormat *inchiFormat = conv.FindFormat("inchi");

  for (c = 1;; ++c) {
    mol.Clear();
    conv.Read(&mol, &ifs);
    if (mol.Empty())
      break;

    cout << can2Conv.WriteString(&mol);
  } // end for loop
  
  return(1);
}


