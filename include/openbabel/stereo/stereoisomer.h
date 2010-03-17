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
#include <openbabel/permutation.h> 

#include <vector> 

namespace OpenBabel {

  class OBMol;

  ///@addtogroup stereo Stereochemistry 
  ///@{

  /**
   * @brief Class for finsing all stereoisomers for a molecule.
   *
   * The OBStereoisomer class finds all stereoisomers for a molecule.
   *
     @verbatim
     Reference:
     [1] M. Razinger, K. Balasubramanian, M. Perdih, M. E. Munk, Stereoisomer
     Generation in Computer-Enhanced Structure Elucidation, J. Chem. Inf.
     Comput. Sci. 1993, 33, 812-825
     @endverbatim
   */

  class OBAPI OBStereoisomer 
  {
    public:
      /**
       * A std::vector of stereo parities. The parities are sorted tetrahedral 
       * first, cis/trans second. Tetrahedral are sorted by the center ids. 
       * Cis/Trans bond parities are sorted by lowest begin or end id.
       */ 
      typedef std::vector<int> ParityVec; // e.g 1 -1 1
      /**
       * A single enantiomer. 
       * @see GetEnantiomers ParityVec
       */
      struct Enantiomer {
        std::vector<ParityVec> parities;
        std::vector<ParityVec> inverseParities;
      };
      /**
       * A single diastereomer. 
       * @see GetDiastereomers ParityVec
       */
      struct Diastereomer {
        std::vector<ParityVec> parities;
      };

      /**
       * Constructor.
       * @param mol The molecule.
       */
      OBStereoisomer(OBMol *mol);

      /**
       * Get the number of enantiomers found.
       */
      unsigned int NumEnantiomerPairs() const
      {
        return m_numEnantiomerPairs;
      }
      /**
       * Get the number of diastereomers found.
       */
      unsigned int NumDiastereomers() const
      {
        return m_numDiastereomers;
      }

      /**
       * Get the found enantiomers as a std::vector of Enantiomer structs 
       * containing the parities.
       */
      const std::vector<Enantiomer>& GetEnantiomers() const
      {
        return m_enantiomers;
      }
      /**
       * Get the found diastereomers as a std::vector of Diastereomer structs
       * containing the parities.
       */
      const std::vector<Diastereomer>& GetDiastereomers() const
      {
        return m_diastereomers;
      }

      const OBPermutationGroup& GetAutomorphisms() const
      {
        return m_automorphisms;
      }


    protected:
      void FindStereoisomers(OBMol *mol, const std::vector<unsigned int> &symmetry_classes);
      unsigned int m_numEnantiomerPairs;
      unsigned int m_numDiastereomers;

      std::vector<Enantiomer> m_enantiomers;
      std::vector<Diastereomer> m_diastereomers;
      OBPermutationGroup m_automorphisms;
  };
     
  ///@}  addtogroup
}

#endif
