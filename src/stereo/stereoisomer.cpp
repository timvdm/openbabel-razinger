/**********************************************************************
  stereoisomer.cpp - OBStereoisomer

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

#include <openbabel/stereo/stereoisomer.h>
#include <openbabel/stereo/tetrahedral.h>
#include <openbabel/stereo/cistrans.h>
#include <openbabel/mol.h>
#include <openbabel/graphsym.h>

#include <openbabel/permutation.h>
    
using namespace std;

#include "paritymatrix.h"
#include "temp.h"

//#define DEBUG_STEREOISOMER 1

namespace OpenBabel {

  OBStereoisomer::OBStereoisomer(OBMol *mol)
  {
    OBGraphSym graphSym(mol);
    std::vector<unsigned int> symmetry_classes;
    graphSym.GetSymmetry(symmetry_classes, false);

    FindStereoisomers(mol, symmetry_classes);
  }
  
  /*
  OBStereoisomer::OBStereoisomer(OBMol *mol, const std::vector<unsigned int> &symmetry_classes)
  {
    FindStereoisomers(mol, symmetry_classes);
  }
  */

  OBStereoisomer::ParityVec createParityVec(const Eigen::VectorXi &v)
  {
    OBStereoisomer::ParityVec pv;
    for (unsigned int i = 0; i < v.size(); ++i)
      pv.push_back(v[i]);
    return pv;
  }

  OBStereoisomer::Enantiomer createEnantiomer(const std::vector<OBStereoisomer::ParityVec> &parities, 
      const std::vector<OBStereoisomer::ParityVec> &inverseParities)
  {
    OBStereoisomer::Enantiomer e;
    e.parities = parities;
    e.inverseParities = inverseParities;
    return e;
  }

  OBStereoisomer::Diastereomer createDiastereomer(const std::vector<OBStereoisomer::ParityVec> &parities)
  {
    OBStereoisomer::Diastereomer d;
    d.parities = parities;
    return d;
  }

  bool contains(const std::vector<OBStereoisomer::ParityVec> &parities, const OBStereoisomer::ParityVec &p)
  {
    for (unsigned int i = 0; i < parities.size(); ++i) {
      if (parities.at(i) == p)
        return true;
    }

    return false;
  }

  void OBStereoisomer::FindStereoisomers(OBMol *mol, const std::vector<unsigned int> &symmetry_classes)
  {
    //
    // Find all automorphisms for the structure
    //
    PermutationGroup G = findAutomorphisms(mol, symmetry_classes);
    m_automorphisms = G;
#ifdef DEBUG_STEREOISOMER
    std::cout << "Automorphisms:" << std::endl;
    for (unsigned int g = 0; g < G.permutations.size(); ++g) {
      G.at(g).print();
    }
#endif

    //
    // Find all types of stereogenic units (i.e. tetrahedral, cis/trans, ...)
    //
    vector<StereogenicUnit> stereoUnits = FindStereogenicUnits(mol, symmetry_classes); // FIXME need to cache this

    unsigned int numTetrahedral = 0;
    unsigned int numCisTrans = 0;

    // make a list of stereo atom indexes
    std::vector<unsigned int> stereoAtoms;
    for (vector<StereogenicUnit>::iterator unit = stereoUnits.begin(); unit != stereoUnits.end(); ++unit) {
      if (unit->type == OBStereo::Tetrahedral) {
        OBAtom *atom = mol->GetAtomById(unit->id);
        stereoAtoms.push_back(atom->GetIndex()+1);
        numTetrahedral++;
      } else if (unit->type == OBStereo::CisTrans) {
        OBBond *bond = mol->GetBondById(unit->id);
        OBAtom *begin = bond->GetBeginAtom();
        OBAtom *end = bond->GetEndAtom();
        stereoAtoms.push_back(begin->GetIndex()+1);
        stereoAtoms.push_back(end->GetIndex()+1);
        numCisTrans++;
      }
    }
      
    unsigned int n = numTetrahedral + 2 * numCisTrans;

    if (!n) {
      m_numEnantiomerPairs = 0;
      m_numDiastereomers = 0;
    }

    //
    // Contract the permutations to only account for stereocenters
    // (i.e. remove non-stereogenic atoms)
    //
    PermutationGroup Gc; // contracted automorphism permutations
    contractPermutations(G, Gc, stereoAtoms);
    if (!Gc.size())
      return;

    //
    // Construct the stereoindex vectors, one for each automorphism:
    //
    // - true-stereo atoms: +1
    // - para-stereo atoms:
    //   - identical ligands not interchanged: +1 (or even number of permutation)
    //   - identical liganss interchanged: -1 (or odd number of permutations)
    //
    std::vector<Eigen::VectorXi> stereoIndexVectors;
    createStereoIndexVectors(stereoIndexVectors, Gc, G, stereoUnits, mol, symmetry_classes, n);

    n = numTetrahedral + numCisTrans;

    /*
    //
    // Reduce the stereoindex vectors (i.e. CisTrans 2 entries -> 1 entry)
    //
    if (numCisTrans) {
      for (unsigned int idx = 0; idx < stereoIndexVectors.size(); ++idx) {
        Eigen::VectorXi reducedIndex(n);
        // copy the tetrahedral entries
        for (unsigned int j = 0; j < numTetrahedral; ++j)
          reducedIndex[j] = stereoIndexVectors.at(idx)[j];
        // reduce the cis/trans entries
        unsigned int i = numTetrahedral;
        for (unsigned int j = numTetrahedral; j < numTetrahedral + 2 * numCisTrans; j += 2)
          reducedIndex[i++] = stereoIndexVectors.at(idx)[j] * stereoIndexVectors.at(idx)[j+1];

        stereoIndexVectors[idx] = reducedIndex;
      }
    }

    // debug: print out reduced stereoindex vectors
    cout << "Reduced StereoIndexVectors:" << endl;
    for (unsigned int idx = 0; idx < stereoIndexVectors.size(); ++idx)
      cout << "reduced stereoindex:" << endl << stereoIndexVectors.at(idx) << endl;
    */

    //
    // Compute (reduced) signed permutation matrices
    //
    std::vector<Eigen::MatrixXi> signedMatrices;
    signedPermutationMatrices(Gc, stereoIndexVectors, signedMatrices, numTetrahedral, numCisTrans, n);

    Eigen::MatrixXi parityMatrix = stereoParityMatrix(n);

    //
    // We now have:
    // - contracted permutation matrices (n x n)
    // - stereoindex vectors (n)
    // - the computer generated stereo parity matrix (2^n x n)
    // 

    // multiply each row of the stereoparity matrix with the signed 
    // permutation matrices and look for redundancies
    unsigned int N = pow(2, stereoUnits.size());
    unsigned int numEnantiomers = 0;
    unsigned int numDiastereomers = 0;
    std::vector<unsigned int> redundant;
    for (unsigned int i = 0; i < N; ++i) {
      if (std::find(redundant.begin(), redundant.end(), i) != redundant.end()) {
        continue;
      }
      
#ifdef DEBUG_STEREOISOMER
      cout << "isomer: " << i+1 << endl;
      cout << "r" << i+1 << ": ";
      for (unsigned int deb = 0; deb < n; deb++)
        cout << parityMatrix(i,deb) << " ";
      cout << endl;
#endif

      std::vector<ParityVec> parities;
      std::vector<ParityVec> inverseParities;

      std::vector<Eigen::VectorXi> products;
      for (unsigned int j = 0; j < signedMatrices.size(); ++j) {
        Eigen::VectorXi Pirj = multiplyMatrixByRow(signedMatrices.at(j), parityMatrix.row(i));

        OBStereoisomer::ParityVec pv = createParityVec(Pirj);
        if (!contains(parities, pv))
          parities.push_back(pv);
        products.push_back(Pirj);
#ifdef DEBUG_STEREOISOMER
        cout << "P" << j+1 << " * r" << i+1 << "^T: ";
        for (unsigned int debug = 0; debug < Pirj.size(); ++debug)
          cout << Pirj[debug] << " ";
        cout << endl;
#endif
        
        for (unsigned int k = 0; k < N; ++k) {
          // skip redundant rows
          if (std::find(redundant.begin(), redundant.end(), k) != redundant.end()) 
            continue;
 
          // compare the product with the row from the parityMatrix
          bool match = true;
          for (unsigned int l = 0; l < n; ++l)
            if (Pirj[l] != parityMatrix(k,l))
              match = false;

          if (match) {
            redundant.push_back(k);
#ifdef DEBUG_STEREOISOMER
            cout << "redundant: " << k+1 << " (found redundant row)" << endl;
#endif
          }
        }
      }
      
      bool foundEnantiomer = false;
      for (unsigned int j = 0; j < products.size(); ++j) {
        Eigen::VectorXi mirror = products.at(j);
        for (unsigned int l = 0; l < numTetrahedral; ++l)
          mirror[l] = -products.at(j)[l];
        
#ifdef DEBUG_STEREOISOMER
        cout << "r-" << i+1 << ": ";
        for (unsigned int deb = 0; deb < n; deb++)
          cout << mirror(deb) << " ";
        cout << endl;
#endif
 
        for (unsigned int k = 0; k < N; ++k) {
          // skip redundant rows
          if (std::find(redundant.begin(), redundant.end(), k) != redundant.end()) 
            continue;
          OBStereoisomer::ParityVec pv = createParityVec(mirror);
          if (!contains(inverseParities, pv))
            inverseParities.push_back(pv);
 
          // compare the product's mirror to the parityMatrix row
          bool match = true;
          for (unsigned int l = 0; l < n; ++l)
            if (mirror[l] != parityMatrix(k,l))
              match = false;

          if (match) {
            redundant.push_back(k);
#ifdef DEBUG_STEREOISOMER
            cout << "redundant: " << k+1 << " (found enantiomer)" << endl;
#endif
            foundEnantiomer = true;
            break;
          }
        }

        //if (foundEnantiomer)
        //  break;
      }
      
      if (foundEnantiomer) {
#ifdef DEBUG_STEREOISOMER
        cout << "--------------------> found enantiomer pair!!" << endl;
#endif
        numEnantiomers++;
        m_enantiomers.push_back(createEnantiomer(parities, inverseParities));
      } else {
#ifdef DEBUG_STEREOISOMER
        cout << "--------------------> found diastereoisomer!!" << endl;
#endif
        numDiastereomers++;
        m_diastereomers.push_back(createDiastereomer(parities));
      }

   
    }

    cout << "numEnantionmers: " << 2*numEnantiomers << endl;
    cout << "numDiastereomers: " << numDiastereomers << endl;

    m_numEnantiomerPairs = numEnantiomers;
    m_numDiastereomers = numDiastereomers;

 
  
  }


}

