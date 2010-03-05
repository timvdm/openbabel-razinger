/**********************************************************************
  stereo.h - OBStereo & OBStereoBase

  Copyright (C) 2009 by Tim Vandermeersch
 
  This file is part of the Open Babel project.
  For more information, see <http://openbabel.sourceforge.net/>

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
  02110-1301, USA.
 **********************************************************************/
#ifndef OB_STEREOMER_H
#define OB_STEREOMER_H

#include <openbabel/babelconfig.h> 

#include <vector> 

namespace OpenBabel {

  class OBMol;

  ///@addtogroup stereo Stereochemistry 
  ///@{

  class OBAPI OBStereoisomer 
  {
    public:
      OBStereoisomer(OBMol *mol);
//      OBStereoisomer(OBMol *mol, const std::vector<unsigned int> &symmetry_classes);

      unsigned int numEnantiomerPairs() const
      {
        return m_enantiomerPairs;
      }
      unsigned int numDiastereomers() const
      {
        return m_diastereomers;
      }

    protected:
      void FindStereoisomers(OBMol *mol, const std::vector<unsigned int> &symmetry_classes);
      unsigned int m_enantiomerPairs;
      unsigned int m_diastereomers;
  };
     
  ///@}  addtogroup
}

#endif
