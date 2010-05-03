/**********************************************************************
  perception.cpp - Stereochemistry perception

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
#include <openbabel/stereo/tetrahedral.h>
#include <openbabel/stereo/cistrans.h>
#include <openbabel/mol.h>
#include <openbabel/graphsym.h>
#include <openbabel/oberror.h>

#include <openbabel/depict/depict.h>
#include <openbabel/mcdlutil.h>

#include <limits>
#include <set>

namespace OpenBabel {
 
  ////////////////////////////////////////////////////////////////////////////
  //
  //  General
  //
  ////////////////////////////////////////////////////////////////////////////

  void PerceiveStereo(OBMol *mol, bool force)
  {
    switch (mol->GetDimension()) {
      case 3:
        StereoFrom3D(mol, force);
        break;
      case 2:
        StereoFrom2D(mol, 0, force);
        break;
      default:
        StereoFrom0D(mol);
        break;
    }
    
    obErrorLog.ThrowError(__FUNCTION__, "Ran OpenBabel::PerceiveStereo", obAuditMsg);
  }

  enum NeighborSymmetryClasses 
  { 
    // Tetrahedral
    T1234 = 1, // 4 different symmetry classes 
    T1123 = 2, // 3 different symmetry classes, 1 class duplicated (2 times)
    T1122 = 3, // 2 different symmetry classes, 1 class duplicated (3 times)
    T1112 = 4, // 2 different symmetry classes, each class duplicated (2 times)
    T1111 = 5, // 1 symmetry class, duplictaed 4 times
    // CisTrans
    C12 = 6, // 2 different symmetry classes
    C11 = 7 // the same symmetry class
  };

  bool mayHaveChiralCenter(OBMol *mol)
  {
    std::vector<OBAtom*>::iterator ia;
    for (OBAtom *atom = mol->BeginAtom(ia); atom; atom = mol->NextAtom(ia))
      if (atom->GetHyb() == 3 && atom->GetHvyValence() >= 3) {
        return true;
      }
    return false;
  }

  bool mayHaveCisTransBond(OBMol *mol)
  {
    std::vector<OBBond*>::iterator ib;
    for (OBBond *bond = mol->BeginBond(ib); bond; bond = mol->NextBond(ib))
      if (bond->GetBO() == 2 && !bond->IsInRing()) {
        return true;
      }
    return false;
  }

  bool isPotentialTetrahedral(OBAtom *atom)
  {
    // consider only potential steroecenters
    if (atom->GetHyb() != 3 || atom->GetHvyValence() < 3)
      return false;
    // skip non-chiral N
    if (atom->IsNitrogen()) {
      int nbrRingAtomCount = 0;
      FOR_NBORS_OF_ATOM (nbr, atom) {
        if (nbr->IsInRing())
          nbrRingAtomCount++;        
      }
      if (nbrRingAtomCount < 3)
        return false;
    }

    return true;
  }

  bool isPotentialCisTrans(OBBond *bond)
  {
    if (bond->GetBondOrder() != 2)
      return false;
    if (bond->IsInRing())
      return false;
    if (!bond->GetBeginAtom()->HasSingleBond() || !bond->GetEndAtom()->HasSingleBond()) 
      return false;
    return true;
  }

  int classifyTetrahedralNbrSymClasses(const std::vector<unsigned int> &symClasses, OBAtom *atom)
  {
    std::vector<unsigned int> nbrClasses, nbrClassesCopy, uniqueClasses;
    FOR_NBORS_OF_ATOM (nbr, atom)
      nbrClasses.push_back(symClasses.at(nbr->GetIndex()));
    if (nbrClasses.size() == 3)
      nbrClasses.push_back(OBStereo::ImplicitRef);

    nbrClassesCopy = nbrClasses; // keep copy for count below
    std::sort(nbrClasses.begin(), nbrClasses.end());
    std::vector<unsigned int>::iterator endLoc = std::unique(nbrClasses.begin(), nbrClasses.end());
    std::copy(nbrClasses.begin(), endLoc, std::back_inserter(uniqueClasses));
    
    switch (uniqueClasses.size()) {
      case 4:
        return T1234; // e.g. 1 2 3 4
      case 3:
        return T1123; // e.g. 1 1 2 3
      case 2:
        if (std::count(nbrClassesCopy.begin(), nbrClassesCopy.end(), uniqueClasses.at(0)) == 2)
          return T1122; // e.g. 1 1 2 2
        else
          return T1112; // e.g. 1 1 1 2
      case 1:
	  default:
        return T1111; // e.g. 1 1 1 1
    }
  }

  int classifyCisTransNbrSymClasses(const std::vector<unsigned int> &symClasses, OBBond *doubleBond, OBAtom *atom)
  {
    std::vector<unsigned int> nbrClasses, uniqueClasses;
    FOR_NBORS_OF_ATOM (nbr, atom) {
      if (nbr->GetIdx() != doubleBond->GetNbrAtom(atom)->GetIdx())
        nbrClasses.push_back(symClasses.at(nbr->GetIndex()));
    }

    if (nbrClasses.size() == 1)
      nbrClasses.push_back(OBStereo::ImplicitRef);

    if (nbrClasses.at(0) == nbrClasses.at(1))
      return C11; // e.g. 1 1
    else
      return C12; // e.g. 1 2
  }



  std::vector<OBBitVec> mergeRings(OBMol *mol)
  {
    std::vector<OBRing*> rings = mol->GetSSSR();

    std::vector<OBBitVec> result;
    for (unsigned int i = 0; i < rings.size(); ++i) {
      // check if ring shares atom with previously found ring
      bool found = false;
      for (unsigned int j = 0; j < result.size(); ++j) {
        for (unsigned int k = 0; k < rings[i]->_path.size(); ++k) {
          if (result[j].BitIsSet(rings[i]->_path[k])) {
            found = true;
            for (unsigned int l = 0; l < rings[i]->_path.size(); ++l)
              result[j].SetBitOn(rings[i]->_path[l]);
            break;
          }
        }

        if (found) 
          break;
      }

      if (!found) {
        OBBitVec r;
        for (unsigned int l = 0; l < rings[i]->_path.size(); ++l)
          r.SetBitOn(rings[i]->_path[l]);
        result.push_back(r);
      }
    }
      
    return result;      
  }

  bool isInSameMergedRing(const std::vector<OBBitVec> &mergedRings, unsigned int idx1, unsigned int idx2)
  {
    std::vector<OBBitVec>::const_iterator bits;
    for (bits = mergedRings.begin(); bits != mergedRings.end(); ++bits)
      if ((*bits).BitIsSet( idx1 ) && (*bits).BitIsSet( idx2 ))
        return true;
    return false;  
  }
  
  void addNbrs(OBBitVec &fragment, OBAtom *atom, OBAtom *skip)
  {
    FOR_NBORS_OF_ATOM (nbr, atom) {
      if (nbr->GetId() == skip->GetId())
        continue;
      if (!fragment.BitIsSet(nbr->GetId())) {
        fragment.SetBitOn(nbr->GetId());
        addNbrs(fragment, &*nbr, skip);
      }
    }
  }

  OBBitVec getFragment(OBAtom *atom, OBAtom *skip)
  { 
    OBBitVec fragment;

    fragment.SetBitOn(atom->GetId());
    addNbrs(fragment, atom, skip);
    
    return fragment; 
  }

  std::vector<std::vector<StereogenicUnit> > sortParaCentersByMergedRings(OBMol *mol,
      const std::vector<OBBitVec> &mergedRings, const std::vector<unsigned int> &paraAtoms,
      const std::vector<unsigned int> &paraBonds)
  {
    std::vector<std::vector<StereogenicUnit> > result;

    for (std::vector<OBBitVec>::const_iterator bv = mergedRings.begin(); bv != mergedRings.end(); ++bv) {
      std::vector<StereogenicUnit> subresult;

      for (std::vector<unsigned int>::const_iterator pa = paraAtoms.begin(); pa != paraAtoms.end(); ++pa)
        if (bv->BitIsOn(*pa))
          subresult.push_back(StereogenicUnit(OBStereo::Tetrahedral, mol->GetAtom(*pa)->GetId()));
      for (std::vector<unsigned int>::const_iterator pb = paraBonds.begin(); pb != paraBonds.end(); ++pb) {
        OBBond *bond  = mol->GetBond(*pb);
        if (bv->BitIsOn(bond->GetBeginAtomIdx()) || bv->BitIsOn(bond->GetEndAtomIdx()))
          subresult.push_back(StereogenicUnit(OBStereo::CisTrans, mol->GetBond(*pb)->GetId()));
      }

      // add subresult, even if empty (allows for ring identification by index)
      result.push_back(subresult);
    }

    return result;
  }



  std::vector<StereogenicUnit> FindStereogenicUnits(OBMol *mol, 
      const std::vector<unsigned int> &symClasses)
  {
    std::vector<StereogenicUnit> units;

    // do quick test to see if there are any possible stereogenic units
    if (!mayHaveChiralCenter(mol) && !mayHaveCisTransBond(mol))
      return units;

    // make sure we have symmetry classes for all atoms
    if (symClasses.size() != mol->NumAtoms())
      return units;

    // para-stereocenters canditates
    std::vector<unsigned int> paraAtoms; 
    std::vector<unsigned int> paraBonds; 

    /**
     * true Tetrahedral stereocenters:
     * - have four different symmetry classes for the ligands to the central atom
     */
    bool ischiral;
    std::vector<OBAtom*>::iterator ia;
    for (OBAtom *atom = mol->BeginAtom(ia); atom; atom = mol->NextAtom(ia)) {
      if (!isPotentialTetrahedral(atom))
        continue;
        
      // list containing neighbor symmetry classes
      std::vector<unsigned int> tlist; 
      ischiral = true;

      // check neighbors to see if this atom is stereogenic
      std::vector<OBBond*>::iterator j;
      for (OBAtom *nbr = atom->BeginNbrAtom(j); nbr; nbr = atom->NextNbrAtom(j)) {
        // check if we already have a neighbor with this symmetry class
        std::vector<unsigned int>::iterator k;
        for (k = tlist.begin(); k != tlist.end(); ++k)
          if (symClasses[nbr->GetIndex()] == *k) {
            ischiral = false;
            // if so, might still be a para-stereocenter
            paraAtoms.push_back(atom->GetIdx());
          }

        if (ischiral)
          // keep track of all neighbors, so we can detect duplicates
          tlist.push_back(symClasses[nbr->GetIndex()]);
        else
          break;
      }

      if (ischiral) {
        // true-stereocenter found
        units.push_back(StereogenicUnit(OBStereo::Tetrahedral, atom->GetId()));
      }
    }

    /**
     * true CisTrans stereocenters:
     * - each terminal has two different symmetry classes for it's ligands
     */
    bool isCisTrans;
    std::vector<OBBond*>::iterator ib;
    for (OBBond *bond = mol->BeginBond(ib); bond; bond = mol->NextBond(ib)) {
      if (bond->IsInRing())
        continue;

      if (bond->GetBO() == 2) {
        OBAtom *begin = bond->GetBeginAtom();
        OBAtom *end = bond->GetEndAtom();
        if (!begin || !end) 
          continue;

        // Needs to have at least one explicit single bond at either end
        // FIXME: timvdm: what about C=C=C=C
        if (!begin->HasSingleBond() || !end->HasSingleBond())
          continue;
          
        isCisTrans = true;
        std::vector<OBBond*>::iterator j;
         
        if (begin->GetValence() == 2) {
          // Begin atom has two explicit neighbors. One is the end atom. The other should 
          // be a heavy atom - this is what we test here.
          // (There is a third, implicit, neighbor which is either a hydrogen
          // or a lone pair.)
          if (begin->ExplicitHydrogenCount() == 1)
            isCisTrans = false;
        } else if (begin->GetValence() == 3) {
          std::vector<unsigned int> tlist;
          
          for (OBAtom *nbr = begin->BeginNbrAtom(j); nbr; nbr = begin->NextNbrAtom(j)) {
            // skip end atom
            if (nbr->GetId() == end->GetId())
              continue;
            // do we already have an atom with this symmetry class?
            if (tlist.size()) {
              // compare second with first
              if (symClasses[nbr->GetIndex()] == tlist.at(0)) {
                isCisTrans = false;
                // if same, might still be a para-stereocenter
                paraBonds.push_back(bond->GetIdx());
              }
              break;
            }
              
            // save first summetry class
            tlist.push_back(symClasses[nbr->GetIndex()]);
          }
        } else {
          // Valence is not 2 or 3, for example SR3=NR
          isCisTrans = false;
        }

        if (!isCisTrans) 
          continue;

        if (end->GetValence() == 2) {
          // see comment above for begin atom
          if (end->ExplicitHydrogenCount() == 1)
            isCisTrans = false;
        } else if (end->GetValence() == 3) { 
          std::vector<unsigned int> tlist;
          
          for (OBAtom *nbr = end->BeginNbrAtom(j); nbr; nbr = end->NextNbrAtom(j)) {
            // skip end atom
            if (nbr->GetId() == begin->GetId())
              continue;
            // do we already have an atom with this symmetry class?
            if (tlist.size()) {
              // compare second with first
              if (symClasses[nbr->GetIndex()] == tlist.at(0)) {
                // if same, might still be a para-stereocenter
                paraBonds.push_back(bond->GetIdx());
                isCisTrans = false;
              }
              break;
            }
                
            // save first summetry class
            tlist.push_back(symClasses[nbr->GetIndex()]);
          }
        } else {
          // Valence is not 2 or 3, for example SR3=NR
          isCisTrans = false;
        }

        if (isCisTrans)
          // true-stereocenter found
          units.push_back(StereogenicUnit(OBStereo::CisTrans, bond->GetId()));
      }
    }

    /**
     * Apply rule 1 from the Razinger paper recusively:
     * 
     * All rings are merged "mergedRings". A merged ring is simply a fragment consisting
     * of all atoms of a ring system (bridged, spiro, adjacent, ...). If two rings in the
     * SSSR set share an atom, they are merged.
     *
     * Each merged must at least have two para-stereocenters (or 1 true + 1 para) in order
     * for the para-stereocenter to be valid. This is repeated until no new stereocenters
     * are identified.
     *
     * rule 1 for double bonds: 
     * - bond atom in ring has two identical symmetry classes for it's neighbor atoms (-> para)
     * - other bond atom:
     *   - has two different symmetry classes for it's neighbours -> new stereocenter
     *   - has two identical symmetry classes, but the ligand contains at least 1 true or para stereocenter -> new stereocenter
     *
     * rule 1 for tetracoord atoms:
     * - at least two neighbour symmetry classes are the same (-> para)
     * - other pair:
     *   - has two different symmetry classes for it's neighbours -> new stereocenter
     *   - has two identical symmetry classes, but the ligand contains at least 1 true or para stereocenter -> new stereocenter
     *
     * NOTE: there must always be at least 2 new stereocenters (or one existing + 1 newly found) in order for them to be valid
     */
    std::vector<OBBitVec> mergedRings = mergeRings(mol);
    std::vector<std::vector<StereogenicUnit> > sortedParas = sortParaCentersByMergedRings(mol, mergedRings, paraAtoms, paraBonds);
    unsigned int lastSize = units.size();
    while (true) {
      // iterate over the merged rings
      for (unsigned int s = 0; s < sortedParas.size(); ++s) {
        int centersInRing = 0;
        std::vector<unsigned long> newAtoms;
        std::vector<unsigned long> newBonds;

        // check for true-stereocenters in the ring
        for (std::vector<StereogenicUnit>::iterator u = units.begin(); u != units.end(); ++u) {
          if ((*u).type == OBStereo::Tetrahedral) {
            OBAtom *atom = mol->GetAtomById((*u).id);
            if (mergedRings.at(s).BitIsOn(atom->GetIdx()))
              centersInRing++;     
          } 
        }

        // check for para-stereocenters in the ring
        for (std::vector<StereogenicUnit>::iterator u = sortedParas[s].begin(); u != sortedParas[s].end(); ++u) {
          if ((*u).type == OBStereo::Tetrahedral) {
            OBAtom *atom = mol->GetAtomById((*u).id);
            int classification = classifyTetrahedralNbrSymClasses(symClasses, atom);
            switch (classification) {
              case T1123:
                {
                  // find the two ligands with identical symmetry classes
                  OBAtom *ligandAtom1 = 0;
                  OBAtom *ligandAtom2 = 0;

                  std::vector<unsigned int> nbrClasses;
                  FOR_NBORS_OF_ATOM (nbr, atom)
                    nbrClasses.push_back(symClasses.at(nbr->GetIndex()));

                  unsigned int duplicatedSymmetryClass;
                  std::sort(nbrClasses.begin(), nbrClasses.end());
                  for (unsigned int i = 0; i < nbrClasses.size()-1; ++i)
                    if (nbrClasses.at(i) == nbrClasses.at(i+1))
                      duplicatedSymmetryClass = nbrClasses.at(i);
 
                  FOR_NBORS_OF_ATOM (nbr, atom) {
                    if (symClasses.at(nbr->GetIndex()) == duplicatedSymmetryClass) {
                      if (!ligandAtom1)
                        ligandAtom1 = &*nbr;
                      else
                        ligandAtom2 = &*nbr;
                    }
                  }

                  // if the two identical symmetry classes are in a ring, this center is stereogenic
                  if (mergedRings[s].BitIsOn(ligandAtom1->GetIdx()) && mergedRings[s].BitIsOn(ligandAtom1->GetIdx())) {
                    centersInRing++;
                    newAtoms.push_back((*u).id);
                  }
                }
                break;
              case T1122:
                {
                  // find the two different ligands
                  OBAtom *ligandAtom1 = 0;
                  OBAtom *ligandAtom2 = 0;
                  FOR_NBORS_OF_ATOM (nbr, atom) {
                    if (!ligandAtom1)
                      ligandAtom1 = &*nbr;
                    else {
                      if (symClasses.at(ligandAtom1->GetIndex()) != symClasses.at(nbr->GetIndex())) {
                        ligandAtom2 = &*nbr;
                      }
                    }
                  }

                  // check which ligand corresponds to the current merged ring
                  OBBitVec ligand;
                  bool spiroCarbon = false;
                  if (mergedRings[s].BitIsOn(ligandAtom1->GetIdx()) && mergedRings[s].BitIsOn(ligandAtom2->GetIdx())) {
                    ligand = getFragment(ligandAtom1, atom) | getFragment(ligandAtom2, atom);
                    spiroCarbon = true;
                  } else if (mergedRings[s].BitIsOn(ligandAtom1->GetIdx())) {
                    ligand = getFragment(ligandAtom2, atom);
                  } else {
                    ligand = getFragment(ligandAtom1, atom);
                  }

                  // make sure there is at least one steroecenter (true/para) in the ligand
                  bool foundStereoCenterInLigand = false;
                  for (std::vector<StereogenicUnit>::iterator u2 = units.begin(); u2 != units.end(); ++u2) {
                    if ((*u2).type == OBStereo::Tetrahedral) {
                      if (ligand.BitIsOn((*u2).id)) {
                        if (spiroCarbon)
                          centersInRing++;
                        foundStereoCenterInLigand = true;
                      }
                    } else if((*u2).type == OBStereo::CisTrans) {
                      OBBond *bond = mol->GetBondById((*u2).id);
                      OBAtom *begin = bond->GetBeginAtom();
                      OBAtom *end = bond->GetEndAtom();
                      if (ligand.BitIsOn(begin->GetId()) || ligand.BitIsOn(end->GetId())) {
                        if (spiroCarbon)
                          centersInRing++;
                        foundStereoCenterInLigand = true;
                      }
                    }
                  }

                  if (foundStereoCenterInLigand) {
                    centersInRing++;
                    newAtoms.push_back(atom->GetId());
                  }
                }
                break;
              case T1111:
                {
                  OBAtom *ligandAtom = 0;
                  FOR_NBORS_OF_ATOM (nbr, atom) {
                    ligandAtom = &*nbr;
                    break;
                  }
                  
                  OBBitVec ligand = getFragment(ligandAtom, atom);
                  bool foundStereoCenterInLigand = false;
                  for (std::vector<StereogenicUnit>::iterator u2 = units.begin(); u2 != units.end(); ++u2) {
                    if ((*u2).type == OBStereo::Tetrahedral) {
                      if (ligand.BitIsOn((*u2).id))
                        foundStereoCenterInLigand = true;
                    } else if((*u2).type == OBStereo::CisTrans) {
                      OBBond *bond = mol->GetBondById((*u2).id);
                      OBAtom *begin = bond->GetBeginAtom();
                      OBAtom *end = bond->GetEndAtom();
                      if (ligand.BitIsOn(begin->GetId()) || ligand.BitIsOn(end->GetId()))
                        foundStereoCenterInLigand = true;
                    }
                  }

                  if (foundStereoCenterInLigand) {
                    centersInRing++;
                    newAtoms.push_back(atom->GetId());
                  }
                }
                break;
              case T1112:
                break;
              default:
                break;
            }
          } else {
            OBBond *bond = mol->GetBondById((*u).id);
            OBAtom *atom = 0;
            if (mergedRings[s].BitIsOn(bond->GetBeginAtomIdx()))
              atom = bond->GetEndAtom();
            else
              atom = bond->GetBeginAtom();
            int classification = classifyCisTransNbrSymClasses(symClasses, bond, atom);
            switch (classification) {
              case C12:
                // again, easy :-)
                centersInRing++;
                newBonds.push_back((*u).id);
                break;
              case C11:
                {
                  // find the ligand
                  OBAtom *ligandAtom = 0;
                  FOR_NBORS_OF_ATOM (nbr, atom) {
                    if ((nbr->GetIdx() != bond->GetBeginAtomIdx()) && (nbr->GetIdx() != bond->GetEndAtomIdx())) {
                      ligandAtom = &*nbr;
                      break;
                    }
                  }

                  OBBitVec ligand = getFragment(ligandAtom, atom);
                  bool foundStereoCenterInLigand = false;
                  for (std::vector<StereogenicUnit>::iterator u2 = units.begin(); u2 != units.end(); ++u2) {
                    if ((*u2).type == OBStereo::Tetrahedral) {
                      if (ligand.BitIsOn((*u2).id))
                        foundStereoCenterInLigand = true;
                    } else if((*u2).type == OBStereo::CisTrans) {
                      OBBond *bond = mol->GetBondById((*u2).id);
                      OBAtom *begin = bond->GetBeginAtom();
                      OBAtom *end = bond->GetEndAtom();
                      if (ligand.BitIsOn(begin->GetId()) || ligand.BitIsOn(end->GetId()))
                        foundStereoCenterInLigand = true;
                    }
                  }

                  if (foundStereoCenterInLigand) {
                    centersInRing++;
                    newBonds.push_back((*u).id);
                  }
                }
                break;

              default:
                break;
            }
          }
        }

        if (centersInRing < 2)
          continue;
        
        // filter out the newly added para-stereocenters in the set of units sorted by ring (sortedParas)
        std::vector<StereogenicUnit> filtered;
        for (std::vector<StereogenicUnit>::iterator u = sortedParas[s].begin(); u != sortedParas[s].end(); ++u) {
          if ((*u).type == OBStereo::Tetrahedral) {
            if (std::find(newAtoms.begin(), newAtoms.end(), (*u).id) == newAtoms.end())
              filtered.push_back(*u);
            else {
              units.push_back(StereogenicUnit(OBStereo::Tetrahedral, (*u).id, true));
            }
          } else if ((*u).type == OBStereo::CisTrans) {
            if (std::find(newBonds.begin(), newBonds.end(), (*u).id) == newBonds.end())
              filtered.push_back(*u);
            else
              units.push_back(StereogenicUnit(OBStereo::CisTrans, (*u).id, true));
          }
        }
        sortedParas[s] = filtered;

      }


      if (units.size() == lastSize)
        break;
      lastSize = units.size();
    }

    /**
     * Apply rule 2a for tetracoordinate carbon:
     * - 1 or 2 pair identical substituents
     * - each pair contains at least 1 true-stereocenter or 2 para-stereocenters
     *
     * Apply rule 2b for tetracoordinate carbon:
     * - 3 or 4 identical substituents with at least
     *   - 2 true-stereocenters
     *   - 2 separate assemblies of para-stereocenters
     */
    for (std::vector<unsigned int>::iterator idx = paraAtoms.begin(); idx != paraAtoms.end(); ++idx) {
      OBAtom *atom = mol->GetAtom(*idx);
      // make sure we didn't add this atom already from rule 1
      bool alreadyAdded = false;
      for (std::vector<StereogenicUnit>::iterator u2 = units.begin(); u2 != units.end(); ++u2) {
        if ((*u2).type == OBStereo::Tetrahedral)
          if (atom->GetId() == (*u2).id) {
            alreadyAdded = true;
          }
      }
      if (alreadyAdded)
        continue;
          
 
      int classification = classifyTetrahedralNbrSymClasses(symClasses, atom);
      switch (classification) {
        case T1123:
          // rule 2a with 1 pair
          {
            // find the duplicated symmetry class
            unsigned int duplicatedSymClass = std::numeric_limits<unsigned int>::max();
            std::vector<unsigned int> nbrSymClasses;
            FOR_NBORS_OF_ATOM (nbr, atom) {
              nbrSymClasses.push_back(symClasses.at(nbr->GetIndex()));
            }
            for (unsigned int i = 0; i < nbrSymClasses.size(); ++i) {
              if (std::count(nbrSymClasses.begin(), nbrSymClasses.end(), nbrSymClasses.at(i)) == 2) {
                duplicatedSymClass = nbrSymClasses.at(i);
                break;
              }
            }
            if (duplicatedSymClass == std::numeric_limits<unsigned int>::max())
              continue;

            // find the ligand atom
            OBAtom *ligandAtom = 0;
            FOR_NBORS_OF_ATOM (nbr, atom)
              if (symClasses.at(nbr->GetIndex()) == duplicatedSymClass)
                ligandAtom = &*nbr;

            // check if ligand contains at least:
            // - 1 true-stereocenter
            // - 2 para-stereocenters
            OBBitVec ligand = getFragment(ligandAtom, atom);
            bool foundTrueStereoCenter = false;
            int paraStereoCenterCount = 0;
            for (std::vector<StereogenicUnit>::iterator u2 = units.begin(); u2 != units.end(); ++u2) {
              if ((*u2).type == OBStereo::Tetrahedral) {
                if (ligand.BitIsOn((*u2).id)) {
                  if ((*u2).para)
                    paraStereoCenterCount++;
                  else
                    foundTrueStereoCenter = true;
                }
              } else if((*u2).type == OBStereo::CisTrans) {
                OBBond *bond = mol->GetBondById((*u2).id);
                OBAtom *begin = bond->GetBeginAtom();
                OBAtom *end = bond->GetEndAtom();
                if (ligand.BitIsOn(begin->GetId()) || ligand.BitIsOn(end->GetId()))
                  if ((*u2).para)
                    paraStereoCenterCount++;
                  else
                    foundTrueStereoCenter = true;
              }
            }

            if (foundTrueStereoCenter || paraStereoCenterCount >= 2) {
              units.push_back(StereogenicUnit(OBStereo::Tetrahedral, atom->GetId(), true));
            }

          }
          break;
        case T1122:
          // rule 2a with 2 pairs
          {
            // find the two different ligands
            OBAtom *ligandAtom1 = 0;
            OBAtom *ligandAtom2 = 0;
            FOR_NBORS_OF_ATOM (nbr, atom) {
              if (!ligandAtom1)
                ligandAtom1 = &*nbr;
              else {
                if (symClasses.at(ligandAtom1->GetIndex()) != symClasses.at(nbr->GetIndex())) {
                  ligandAtom2 = &*nbr;
                }
              }
            }

            // check if ligand contains at least:
            // - 1 true-stereocenter
            // - 2 para-stereocenters
            OBBitVec ligand = getFragment(ligandAtom1, atom);
            bool foundTrueStereoCenter = false;
            int paraStereoCenterCount = 0;
            for (std::vector<StereogenicUnit>::iterator u2 = units.begin(); u2 != units.end(); ++u2) {
              if ((*u2).type == OBStereo::Tetrahedral) {
                if (ligand.BitIsOn((*u2).id)) {
                  if ((*u2).para)
                    paraStereoCenterCount++;
                  else
                    foundTrueStereoCenter = true;
                }
              } else if((*u2).type == OBStereo::CisTrans) {
                OBBond *bond = mol->GetBondById((*u2).id);
                OBAtom *begin = bond->GetBeginAtom();
                OBAtom *end = bond->GetEndAtom();
                if (ligand.BitIsOn(begin->GetId()) || ligand.BitIsOn(end->GetId()))
                  if ((*u2).para)
                    paraStereoCenterCount++;
                  else
                    foundTrueStereoCenter = true;
              }
            }

            if (!foundTrueStereoCenter && paraStereoCenterCount < 2)
              continue;

            ligand = getFragment(ligandAtom1, atom);
            foundTrueStereoCenter = false;
            paraStereoCenterCount = 0;
            for (std::vector<StereogenicUnit>::iterator u2 = units.begin(); u2 != units.end(); ++u2) {
              if ((*u2).type == OBStereo::Tetrahedral) {
                if (ligand.BitIsOn((*u2).id)) {
                  if ((*u2).para)
                    paraStereoCenterCount++;
                  else
                    foundTrueStereoCenter = true;
                }
              } else if((*u2).type == OBStereo::CisTrans) {
                OBBond *bond = mol->GetBondById((*u2).id);
                OBAtom *begin = bond->GetBeginAtom();
                OBAtom *end = bond->GetEndAtom();
                if (ligand.BitIsOn(begin->GetId()) || ligand.BitIsOn(end->GetId()))
                  if ((*u2).para)
                    paraStereoCenterCount++;
                  else
                    foundTrueStereoCenter = true;
              }
            }

            if (foundTrueStereoCenter || paraStereoCenterCount >= 2) {
              units.push_back(StereogenicUnit(OBStereo::Tetrahedral, atom->GetId(), true));
            }
          }
          break;
        case T1112:
          // rule 2b with 3 identical
          {
            // find the duplicated symmetry class
            unsigned int duplicatedSymClass = std::numeric_limits<unsigned int>::max();
            std::vector<unsigned int> nbrSymClasses;
            FOR_NBORS_OF_ATOM (nbr, atom) {
              nbrSymClasses.push_back(symClasses.at(nbr->GetIndex()));
            }
            for (unsigned int i = 0; i < nbrSymClasses.size(); ++i) {
              if (std::count(nbrSymClasses.begin(), nbrSymClasses.end(), nbrSymClasses.at(i)) == 3) {
                duplicatedSymClass = nbrSymClasses.at(i);
                break;
              }
            }
            if (duplicatedSymClass == std::numeric_limits<unsigned int>::max())
              continue;

            // find the ligand atom
            OBAtom *ligandAtom = 0;
            FOR_NBORS_OF_ATOM (nbr, atom)
              if (symClasses.at(nbr->GetIndex()) == duplicatedSymClass)
                ligandAtom = &*nbr;

            if (!ligandAtom)
              break;

            // check if ligand contains at least:
            // - 2 true-stereocenter
            // - 2 separate para-stereocenters assemblies
            OBBitVec ligand = getFragment(ligandAtom, atom);
            int trueStereoCenterCount = 0;
            std::vector<unsigned int> ringIndices;
            for (std::vector<StereogenicUnit>::iterator u2 = units.begin(); u2 != units.end(); ++u2) {
              if ((*u2).type == OBStereo::Tetrahedral) {
                if (ligand.BitIsOn((*u2).id)) {
                  if ((*u2).para) {
                    OBAtom *paraAtom = mol->GetAtomById((*u2).id);
                    for (unsigned int ringIdx = 0; ringIdx < mergedRings.size(); ++ringIdx) {
                      if (mergedRings.at(ringIdx).BitIsOn(paraAtom->GetIdx()))
                        if (std::find(ringIndices.begin(), ringIndices.end(), ringIdx) == ringIndices.end())
                          ringIndices.push_back(ringIdx);
                    }
                  } else
                    trueStereoCenterCount++;
                }
              } else if((*u2).type == OBStereo::CisTrans) {
                OBBond *bond = mol->GetBondById((*u2).id);
                OBAtom *begin = bond->GetBeginAtom();
                OBAtom *end = bond->GetEndAtom();
                if (ligand.BitIsOn(begin->GetId()) || ligand.BitIsOn(end->GetId()))
                  if ((*u2).para) {
                    for (unsigned int ringIdx = 0; ringIdx < mergedRings.size(); ++ringIdx) {
                      if (mergedRings.at(ringIdx).BitIsOn(begin->GetIdx()) || mergedRings.at(ringIdx).BitIsOn(end->GetIdx()))
                        if (std::find(ringIndices.begin(), ringIndices.end(), ringIdx) == ringIndices.end())
                          ringIndices.push_back(ringIdx);
                    }
                  } else
                    trueStereoCenterCount++;
              }
            }

            if (trueStereoCenterCount >= 2 || ringIndices.size() >= 2) {
              units.push_back(StereogenicUnit(OBStereo::Tetrahedral, atom->GetId(), true));
            }

          }
          break;
        case T1111:
          // rule 2b with 4 identical
          {
            // find the ligand atom
            OBAtom *ligandAtom = 0;
            FOR_NBORS_OF_ATOM (nbr, atom) {
              ligandAtom = &*nbr;
              break;
            }

            // check if ligand contains at least:
            // - 2 true-stereocenter
            // - 2 separate para-stereocenters assemblies
            OBBitVec ligand = getFragment(ligandAtom, atom);
            int trueStereoCenterCount = 0;
            std::vector<unsigned int> ringIndices;
            for (std::vector<StereogenicUnit>::iterator u2 = units.begin(); u2 != units.end(); ++u2) {
              if ((*u2).type == OBStereo::Tetrahedral) {
                if (ligand.BitIsOn((*u2).id)) {
                  if ((*u2).para) {
                    OBAtom *paraAtom = mol->GetAtomById((*u2).id);
                    for (unsigned int ringIdx = 0; ringIdx < mergedRings.size(); ++ringIdx) {
                      if (mergedRings.at(ringIdx).BitIsOn(paraAtom->GetIdx()))
                        if (std::find(ringIndices.begin(), ringIndices.end(), ringIdx) == ringIndices.end())
                          ringIndices.push_back(ringIdx);
                    }
                  } else
                    trueStereoCenterCount++;
                }
              } else if((*u2).type == OBStereo::CisTrans) {
                OBBond *bond = mol->GetBondById((*u2).id);
                OBAtom *begin = bond->GetBeginAtom();
                OBAtom *end = bond->GetEndAtom();
                if (ligand.BitIsOn(begin->GetId()) || ligand.BitIsOn(end->GetId()))
                  if ((*u2).para) {
                    for (unsigned int ringIdx = 0; ringIdx < mergedRings.size(); ++ringIdx) {
                      if (mergedRings.at(ringIdx).BitIsOn(begin->GetIdx()) || mergedRings.at(ringIdx).BitIsOn(end->GetIdx()))
                        if (std::find(ringIndices.begin(), ringIndices.end(), ringIdx) == ringIndices.end())
                          ringIndices.push_back(ringIdx);
                    }
                  } else
                    trueStereoCenterCount++;
              }
            }

            if (trueStereoCenterCount >= 2 || ringIndices.size() >= 2) {
              units.push_back(StereogenicUnit(OBStereo::Tetrahedral, atom->GetId(), true));
            }
          }
          break;
      
      }
    
    }

    /**
     * Apply rule 3 for double bonds.
     */
    for (std::vector<unsigned int>::iterator idx = paraBonds.begin(); idx != paraBonds.end(); ++idx) {
      OBBond *bond = mol->GetBond(*idx);

      // make sure we didn't add this atom already from rule 1
      bool alreadyAdded = false;
      for (std::vector<StereogenicUnit>::iterator u2 = units.begin(); u2 != units.end(); ++u2) {
        if ((*u2).type == OBStereo::CisTrans)
          if (bond->GetId() == (*u2).id) {
            alreadyAdded = true;
          }
      }
      if (alreadyAdded)
        continue;
       
      OBAtom *begin = bond->GetBeginAtom();
      OBAtom *end = bond->GetEndAtom();
          
      int beginClassification = classifyCisTransNbrSymClasses(symClasses, bond, bond->GetBeginAtom());
      bool beginValid = false;
      switch (beginClassification) {
        case C12:
          beginValid = true;
          break;
        case C11:
          {
            // find the ligand
            OBAtom *ligandAtom = 0;
            FOR_NBORS_OF_ATOM (nbr, begin) {
              if ((nbr->GetIdx() != bond->GetBeginAtomIdx()) && (nbr->GetIdx() != bond->GetEndAtomIdx())) {
                ligandAtom = &*nbr;
                break;
              }
            }

            OBBitVec ligand = getFragment(ligandAtom, begin);
            for (std::vector<StereogenicUnit>::iterator u2 = units.begin(); u2 != units.end(); ++u2) {
              if ((*u2).type == OBStereo::Tetrahedral) {
                if (ligand.BitIsOn((*u2).id))
                  beginValid = true;
              } else if((*u2).type == OBStereo::CisTrans) {
                OBBond *bond = mol->GetBondById((*u2).id);
                OBAtom *begin = bond->GetBeginAtom();
                OBAtom *end = bond->GetEndAtom();
                if (ligand.BitIsOn(begin->GetId()) || ligand.BitIsOn(end->GetId()))
                  beginValid = true;
              }
            }
          }
          break;
      }
          
      if (!beginValid)      
        continue;
 
      int endClassification = classifyCisTransNbrSymClasses(symClasses, bond, bond->GetEndAtom());
      bool endValid = false;
      switch (endClassification) {
        case C12:
          endValid = true;
          break;
        case C11:
          {
            // find the ligand
            OBAtom *ligandAtom = 0;
            FOR_NBORS_OF_ATOM (nbr, end) {
              if ((nbr->GetIdx() != bond->GetBeginAtomIdx()) && (nbr->GetIdx() != bond->GetEndAtomIdx())) {
                ligandAtom = &*nbr;
                break;
              }
            }

            OBBitVec ligand = getFragment(ligandAtom, end);
            for (std::vector<StereogenicUnit>::iterator u2 = units.begin(); u2 != units.end(); ++u2) {
              if ((*u2).type == OBStereo::Tetrahedral) {
                if (ligand.BitIsOn((*u2).id))
                  endValid = true;
              } else if((*u2).type == OBStereo::CisTrans) {
                OBBond *bond = mol->GetBondById((*u2).id);
                OBAtom *begin = bond->GetBeginAtom();
                OBAtom *end = bond->GetEndAtom();
                if (ligand.BitIsOn(begin->GetId()) || ligand.BitIsOn(end->GetId()))
                  endValid = true;
              }
            }
          }
          break;
      }

      if (endValid)
        units.push_back(StereogenicUnit(OBStereo::CisTrans, bond->GetId(), true));
    }
 
    /*
    cout << "Final True-Tetrahedral: ";
    for (std::vector<StereogenicUnit>::iterator u = units.begin(); u != units.end(); ++u)
      if ((*u).type == OBStereo::Tetrahedral)
        cout << (*u).id << " ";
    cout << endl;
    cout << "Final True-CisTrans: ";
    for (std::vector<StereogenicUnit>::iterator u = units.begin(); u != units.end(); ++u)
      if ((*u).type == OBStereo::CisTrans)
        cout << (*u).id << " ";
    cout << endl;
    */
 
    return units;
  }

  unsigned int findDuplicatedSymmetryClass(OBAtom *atom, const std::vector<unsigned int> &symClasses)
  { 
    // find the duplicated symmetry class
    unsigned int duplicatedSymClass = std::numeric_limits<unsigned int>::max();
    std::vector<unsigned int> nbrSymClasses;
    //cout << "findDuplicatedSymmetryClass(): ";
    FOR_NBORS_OF_ATOM (nbr, atom) {
      nbrSymClasses.push_back(symClasses.at(nbr->GetIndex()));
      //cout << nbrSymClasses.back() << " ";
    }
    //cout << endl;
    for (unsigned int i = 0; i < nbrSymClasses.size(); ++i) {
      if (std::count(nbrSymClasses.begin(), nbrSymClasses.end(), nbrSymClasses.at(i)) >= 2) {
        duplicatedSymClass = nbrSymClasses.at(i);
        break;
      }
    }
    return duplicatedSymClass;
  }

  void findDuplicatedSymmetryClasses(OBAtom *atom, const std::vector<unsigned int> &symClasses,
      unsigned int &duplicated1, unsigned int &duplicated2)
  { 
    std::vector<unsigned int> nbrSymClasses;
    FOR_NBORS_OF_ATOM (nbr, atom)
      nbrSymClasses.push_back(symClasses.at(nbr->GetIndex()));
    std::sort(nbrSymClasses.begin(), nbrSymClasses.end());
    duplicated1 = nbrSymClasses[0];
    duplicated2 = nbrSymClasses[2];
  }


  bool permutationInvertsTetrahedralCenter(const OBPermutation &p, OBAtom *center, 
      const std::vector<unsigned int> &symmetry_classes)
  {
    OBMol *mol = center->GetParent();
    std::vector<OBBitVec> mergedRings = mergeRings(mol);
    unsigned int duplicatedSymClass = findDuplicatedSymmetryClass(center, symmetry_classes);
  
    OBAtom *otherCenter = mol->GetAtom(p.map[center->GetIndex()]);

    unsigned int loopCount = 0;
    bool carry = false;
    while (otherCenter != center) {
      if (loopCount > 100) {
        otherCenter = center;
        carry = false;
        break;
      }
      loopCount++;
      if (isInSameMergedRing(mergedRings, center->GetIndex() + 1, otherCenter->GetIndex() + 1)) {
        otherCenter = center;
        break;
      }

      std::vector<unsigned long> tlist;
      FOR_NBORS_OF_ATOM (nbr, otherCenter) {
        if (symmetry_classes[nbr->GetIndex()] == duplicatedSymClass)
          tlist.push_back(nbr->GetIndex() + 1);
      }

      std::vector<unsigned long> map;
      for (unsigned int i = 0; i < p.map.size(); ++i) {
        if (std::find(tlist.begin(), tlist.end(), p.map[i]) != tlist.end())
          map.push_back(p.map[i]);
      }

      if (OBStereo::NumInversions(map) % 2)
        carry = !carry;

      otherCenter = mol->GetAtom(p.map[otherCenter->GetIndex()]);
    }

    std::vector<unsigned long> tlist;
    FOR_NBORS_OF_ATOM (nbr, center) {
      if (symmetry_classes[nbr->GetIndex()] == duplicatedSymClass)
        tlist.push_back(nbr->GetIndex() + 1);
    }

    std::vector<unsigned long> map;
    for (unsigned int i = 0; i < p.map.size(); ++i) {
      if (std::find(tlist.begin(), tlist.end(), p.map[i]) != tlist.end())
        map.push_back(p.map[i]);
    }

    if (carry)
      return (OBStereo::NumInversions(map) % 2) == 0;
    return OBStereo::NumInversions(map) % 2;
  }
  
  bool permutationInvertsCisTransCenter(const OBPermutation &p, OBBond *bond, 
      const std::vector<unsigned int> &symmetry_classes)
  {
    OBMol *mol = bond->GetParent();
      
    std::vector<unsigned int> list;
    std::vector<unsigned long> map;
    // begin
    FOR_NBORS_OF_ATOM (nbr, bond->GetBeginAtom()) {
      if (nbr->GetId() == bond->GetEndAtom()->GetId())
        continue;
      list.push_back(nbr->GetIndex() + 1);
    }
    for (unsigned int i = 0; i < p.map.size(); ++i) {
      if (std::find(list.begin(), list.end(), p.map[i]) != list.end())
        map.push_back(p.map[i]);
    }
    bool beginInverted = OBStereo::NumInversions(map) % 2;
    list.clear();
    map.clear();
    // end
    FOR_NBORS_OF_ATOM (nbr, bond->GetEndAtom()) {
      if (nbr->GetId() == bond->GetBeginAtom()->GetId())
        continue;
      list.push_back(nbr->GetIndex() + 1);
    }
    for (unsigned int i = 0; i < p.map.size(); ++i) {
      if (std::find(list.begin(), list.end(), p.map[i]) != list.end())
        map.push_back(p.map[i]);
    }
    bool endInverted = OBStereo::NumInversions(map) % 2;
 
    if (beginInverted ^ endInverted)
      return true;
    return false;
  }

  bool isTetrahedral(OBAtom *atom, const std::vector<StereogenicUnit> &units);


  struct StereoInverted {
    struct Entry {
      OBPermutation p;
      std::vector<OBAtom*> invertedAtoms;
      std::vector<OBBond*> invertedBonds;
    };

    std::vector<Entry> list;

    static StereoInverted compute(OBMol *mol, const std::vector<unsigned int> &symClasses, 
        const OBPermutationGroup &automorphisms)
    {
      StereoInverted result;

      // make a list of tetrahedral centers inverted by the automorphism permutations
      for (unsigned int i = 0; i < automorphisms.Size(); ++i) {
        Entry entry;
        entry.p = automorphisms.At(i);

        std::vector<OBAtom*>::iterator ia;
        for (OBAtom *atom = mol->BeginAtom(ia); atom; atom = mol->NextAtom(ia)) {
          if (!isPotentialTetrahedral(atom))
            continue;
          if (permutationInvertsTetrahedralCenter(entry.p, atom, symClasses))
            entry.invertedAtoms.push_back(atom);
        }
 
        std::vector<OBBond*>::iterator ib;
        for (OBBond *bond = mol->BeginBond(ib); bond; bond = mol->NextBond(ib)) {
          if (!isPotentialCisTrans(bond))
            continue;
          if (permutationInvertsCisTransCenter(entry.p, bond, symClasses))
            entry.invertedBonds.push_back(bond);
        }
     
        cout << "automorphism " << i+1 << "     "; entry.p.Print();
        cout << "  invertedAtoms: ";
        for (unsigned int l = 0; l < entry.invertedAtoms.size(); ++l)
          cout << entry.invertedAtoms[l]->GetId() << " ";
        cout << endl;
        cout << "  invertedBonds: ";
        for (unsigned int l = 0; l < entry.invertedBonds.size(); ++l)
          cout << entry.invertedBonds[l]->GetId() << " ";
        cout << endl;
 
        result.list.push_back(entry); 
      }
      return result;

    }
 

  };



  OBAtom* findAtomWithSymmetryClass(OBAtom *atom, unsigned int symClass, const std::vector<unsigned int> &symClasses)
  {
    OBAtom *ligandAtom = 0;
    FOR_NBORS_OF_ATOM (nbr, atom)
      if (symClasses.at(nbr->GetIndex()) == symClass)
        ligandAtom = &*nbr;
    return ligandAtom;
  }
 

  bool containsAtLeast_1true_1para(OBAtom *ligandAtom, OBAtom *atom, const std::vector<StereogenicUnit> &units)
  {
    OBMol *mol = atom->GetParent();
    OBBitVec ligand = getFragment(ligandAtom, atom);
    for (std::vector<StereogenicUnit>::const_iterator u2 = units.begin(); u2 != units.end(); ++u2) {
      if ((*u2).type == OBStereo::Tetrahedral) {
        if (ligand.BitIsOn((*u2).id))
          return true;
      } else if((*u2).type == OBStereo::CisTrans) {
        OBBond *bond2 = mol->GetBondById((*u2).id);
        OBAtom *begin = bond2->GetBeginAtom();
        OBAtom *end = bond2->GetEndAtom();
        if (ligand.BitIsOn(begin->GetId()) || ligand.BitIsOn(end->GetId()))
          return true;
      }
    }
    return false;
  }

  bool containsAtLeast_1true_2para(OBAtom *ligandAtom, OBAtom *atom, const std::vector<StereogenicUnit> &units)
  {
    OBMol *mol = atom->GetParent();
    // check if ligand contains at least:
    // - 1 true-stereocenter
    // - 2 para-stereocenters
    OBBitVec ligand = getFragment(ligandAtom, atom);
    bool foundTrueStereoCenter = false;
    int paraStereoCenterCount = 0;
    for (std::vector<StereogenicUnit>::const_iterator u2 = units.begin(); u2 != units.end(); ++u2) {
      if ((*u2).type == OBStereo::Tetrahedral) {
        if (ligand.BitIsOn((*u2).id)) {
          if ((*u2).para)
            paraStereoCenterCount++;
          else
            foundTrueStereoCenter = true;
        }
      } else if((*u2).type == OBStereo::CisTrans) {
        OBBond *bond = mol->GetBondById((*u2).id);
        OBAtom *begin = bond->GetBeginAtom();
        OBAtom *end = bond->GetEndAtom();
        if (ligand.BitIsOn(begin->GetId()) || ligand.BitIsOn(end->GetId()))
          if ((*u2).para)
            paraStereoCenterCount++;
          else
            foundTrueStereoCenter = true;
      }
    }

    if (foundTrueStereoCenter || paraStereoCenterCount >= 2) 
      return true;
    if (ligandAtom->IsInRing() && atom->IsInRing() && paraStereoCenterCount)
      return true;
    return false;
  }
  
  bool containsAtLeast_2true_2paraAssemblies(OBAtom *ligandAtom, OBAtom *atom, const std::vector<StereogenicUnit> &units)
  {
    OBMol *mol = atom->GetParent();
    std::vector<OBBitVec> mergedRings = mergeRings(mol);
    // check if ligand contains at least:
    // - 2 true-stereocenter
    // - 2 separate para-stereocenters assemblies
    OBBitVec ligand = getFragment(ligandAtom, atom);
    int trueStereoCenterCount = 0;
    std::vector<unsigned int> ringIndices;
    for (std::vector<StereogenicUnit>::const_iterator u2 = units.begin(); u2 != units.end(); ++u2) {
      if ((*u2).type == OBStereo::Tetrahedral) {
        if (ligand.BitIsOn((*u2).id)) {
          if ((*u2).para) {
            OBAtom *paraAtom = mol->GetAtomById((*u2).id);
            for (unsigned int ringIdx = 0; ringIdx < mergedRings.size(); ++ringIdx) {
              if (mergedRings.at(ringIdx).BitIsOn(paraAtom->GetIdx()))
                if (std::find(ringIndices.begin(), ringIndices.end(), ringIdx) == ringIndices.end())
                  ringIndices.push_back(ringIdx);
            }
          } else
            trueStereoCenterCount++;
        }
      } else if((*u2).type == OBStereo::CisTrans) {
        OBBond *bond = mol->GetBondById((*u2).id);
        OBAtom *begin = bond->GetBeginAtom();
        OBAtom *end = bond->GetEndAtom();
        if (ligand.BitIsOn(begin->GetId()) || ligand.BitIsOn(end->GetId()))
          if ((*u2).para) {
            for (unsigned int ringIdx = 0; ringIdx < mergedRings.size(); ++ringIdx) {
              if (mergedRings.at(ringIdx).BitIsOn(begin->GetIdx()) || mergedRings.at(ringIdx).BitIsOn(end->GetIdx()))
                if (std::find(ringIndices.begin(), ringIndices.end(), ringIdx) == ringIndices.end())
                  ringIndices.push_back(ringIdx);
            }
          } else
            trueStereoCenterCount++;
      }
    }

    if (trueStereoCenterCount >= 2 || ringIndices.size() >= 2)
      return true;
    return false;
  }


  std::vector<StereogenicUnit> FindStereogenicUnits(OBMol *mol, 
      const std::vector<unsigned int> &symClasses, const OBPermutationGroup &automorphisms)
  {
    std::vector<StereogenicUnit> units;
 
    // do quick test to see if there are any possible stereogenic units
    if (!mayHaveCisTransBond(mol) && !mayHaveChiralCenter(mol))
      return units;

    // make sure we have symmetry classes for all atoms
    if (symClasses.size() != mol->NumAtoms())
      return units;

    StereoInverted inverted = StereoInverted::compute(mol, symClasses, automorphisms);
    
    std::vector<unsigned long> doneAtoms, doneBonds;
    unsigned int lastSize = units.size();
    while (true) {
      std::vector<OBAtom*>::iterator ia;
      for (OBAtom *atom = mol->BeginAtom(ia); atom; atom = mol->NextAtom(ia)) {
        if (std::find(doneAtoms.begin(), doneAtoms.end(), atom->GetId()) != doneAtoms.end())
          continue;
        // consider only potential steroecenters
        if (!isPotentialTetrahedral(atom))
          continue;

        // A potential stereocenter is really a stereocenter if there exists no automorphic
        // permutation causing an inversion of the configuration of only the potential
        // stereogenic unit under consideration.
        bool foundPermutation = false; // invert __only__ configuration of atom
        for (unsigned int i = 0; i < inverted.list.size(); ++i) {
          const std::vector<OBAtom*> &atoms = inverted.list[i].invertedAtoms;
          if (atoms.size() != 1)
            continue;
          const std::vector<OBBond*> &bonds = inverted.list[i].invertedBonds;
          if (bonds.size())
            continue;
          if (atoms[0] == atom) {
            foundPermutation = true;
            break;
          }
        }

        int classification = classifyTetrahedralNbrSymClasses(symClasses, atom);

        if (!foundPermutation) {
          // true-stereocenter found
          bool isParaCenter = (classification == T1234) ? false : true;
          cout << "found(2) " << atom->GetId() << endl;
          units.push_back(StereogenicUnit(OBStereo::Tetrahedral, atom->GetId(), isParaCenter));
          doneAtoms.push_back(atom->GetId());
        } else {
          cout << "need to count... " << atom->GetId() << endl;
          // count ligand configurations:
          // If there exists at least oen automorphic permutation casuing the inversion of the 
          // configuration of only the stereogenic unit under consideration, then the potential
          // stereocenter can be a stereocenter if the number of topologically equivalent neighbors
          // (ligands) of potential stereogenic is less than or equal to the number of configurations
          // of these ligands.
          // 
          // In practise: 
          //    T1123 -> 1 true stereocenter OR 2 para stereocenters
          //    T1122 -> 1 true stereocenter OR 2 para stereocenters (for both)
          //    T1112 -> 2 true stereocenters OR 2 para stereocenter assemblies
          //    T1111 -> 2 true stereocenters OR 2 para stereocenter assemblies
          switch (classification) {
            case T1123:
              {
                unsigned int duplicatedSymClass = findDuplicatedSymmetryClass(atom, symClasses);
                OBAtom *ligandAtom = findAtomWithSymmetryClass(atom, duplicatedSymClass, symClasses);
                if (containsAtLeast_1true_2para(ligandAtom, atom, units)) {
                  units.push_back(StereogenicUnit(OBStereo::Tetrahedral, atom->GetId(), true));
                  doneAtoms.push_back(atom->GetId());
                }
              }
              break;
            case T1122:
              {
                unsigned int duplicatedSymClass1, duplicatedSymClass2;
                findDuplicatedSymmetryClasses(atom, symClasses, duplicatedSymClass1, duplicatedSymClass2);
                OBAtom *ligandAtom1 = findAtomWithSymmetryClass(atom, duplicatedSymClass1, symClasses);
                OBAtom *ligandAtom2 = findAtomWithSymmetryClass(atom, duplicatedSymClass2, symClasses);
                if (containsAtLeast_1true_2para(ligandAtom1, atom, units) &&
                    containsAtLeast_1true_2para(ligandAtom2, atom, units)) {
                  units.push_back(StereogenicUnit(OBStereo::Tetrahedral, atom->GetId(), true));
                  doneAtoms.push_back(atom->GetId());
                  cout << "found " << atom->GetId() << endl;
                }
              }
              break;
            case T1112:
            case T1111:
              {
                unsigned int duplicatedSymClass = findDuplicatedSymmetryClass(atom, symClasses);
                OBAtom *ligandAtom = findAtomWithSymmetryClass(atom, duplicatedSymClass, symClasses);
                if (containsAtLeast_2true_2paraAssemblies(ligandAtom, atom, units)) {
                  units.push_back(StereogenicUnit(OBStereo::Tetrahedral, atom->GetId(), true));
                  doneAtoms.push_back(atom->GetId());
                }
              }
              break;
          }
        }
      }

      std::vector<OBBond*>::iterator ib;
      for (OBBond *bond = mol->BeginBond(ib); bond; bond = mol->NextBond(ib)) {
        if (std::find(doneBonds.begin(), doneBonds.end(), bond->GetId()) != doneBonds.end())
          continue;
        if (!isPotentialCisTrans(bond))
          continue;
 
        // A double bond is a stereogenic bond if there exists no automorphic
        // permutation causing an inversion of the configuration of only the potential
        // stereogenic unit under consideration.
        bool foundPermutation = false; // invert __only__ configuration of atom
        for (unsigned int i = 0; i < inverted.list.size(); ++i) {
          const std::vector<OBAtom*> &atoms = inverted.list[i].invertedAtoms;
          if (atoms.size())
            continue;
          const std::vector<OBBond*> &bonds = inverted.list[i].invertedBonds;
          if (bonds.size() != 1)
            continue;
          if (bonds[0] == bond) {
            foundPermutation = true;
            break;
          }
        }

        int beginClassification = classifyCisTransNbrSymClasses(symClasses, bond, bond->GetBeginAtom());
        int endClassification = classifyCisTransNbrSymClasses(symClasses, bond, bond->GetEndAtom());
      
        if (!foundPermutation) {
          // true-stereocenter found
          bool isParaCenter = (beginClassification == C12) && (endClassification == C12) ? false : true;
          units.push_back(StereogenicUnit(OBStereo::CisTrans, bond->GetId(), isParaCenter));
          doneBonds.push_back(bond->GetId());
        } else {
          // count ligand configurations:
          bool beginValid = false;
          switch (beginClassification) {
            case C12:
              beginValid = true;
              break;
            case C11:
              {
                cout << "need to count bond... " << bond->GetId() << endl;
                // find the ligand
                OBAtom *ligandAtom = 0;
                FOR_NBORS_OF_ATOM (nbr, bond->GetBeginAtom()) {
                  if ((nbr->GetIdx() != bond->GetBeginAtomIdx()) && (nbr->GetIdx() != bond->GetEndAtomIdx())) {
                    ligandAtom = &*nbr;
                    break;
                  }
                }
                
                beginValid = containsAtLeast_1true_1para(ligandAtom, bond->GetBeginAtom(), units);
              }
              break;
          }

          if (!beginValid)      
            continue;

          bool endValid = false;
          switch (endClassification) {
            case C12:
              endValid = true;
              break;
            case C11:
              {
                cout << "need to count bond... " << bond->GetId() << endl;
                // find the ligand
                OBAtom *ligandAtom = 0;
                FOR_NBORS_OF_ATOM (nbr, bond->GetEndAtom()) {
                  if ((nbr->GetIdx() != bond->GetBeginAtomIdx()) && (nbr->GetIdx() != bond->GetEndAtomIdx())) {
                    ligandAtom = &*nbr;
                    break;
                  }
                }

                endValid = containsAtLeast_1true_1para(ligandAtom, bond->GetEndAtom(), units);
              }
              break;
          }

          if (endValid) {
            units.push_back(StereogenicUnit(OBStereo::CisTrans, bond->GetId(), true));
            doneBonds.push_back(bond->GetId());
          }
        } 
      }


      if (units.size() == lastSize)
        break;
      lastSize = units.size(); 
    }

    for (std::vector<StereogenicUnit>::iterator unit = units.begin(); unit != units.end(); ++unit) {
      if (unit->type == OBStereo::Tetrahedral) 
        cout << "Tetrahedral(center = " << unit->id << ", para = " << unit->para << ")" << endl;
      if (unit->type == OBStereo::CisTrans) 
        cout << "CisTrans(bond = " << unit->id << ", para = " << unit->para << ")" << endl;
    }

    return units;
  }
 
  unsigned int getSymClass(const StereogenicUnit &unit, OBMol *mol, const std::vector<unsigned int> &symmetry_classes)
  {
    if (unit.type == OBStereo::Tetrahedral) {
      OBAtom *center = mol->GetAtomById(unit.id);
      return symmetry_classes[center->GetIndex()];
    } else
    if (unit.type == OBStereo::CisTrans) {
      OBBond *bond = mol->GetBondById(unit.id);
      OBAtom *begin = bond->GetBeginAtom();
      OBAtom *end = bond->GetEndAtom();
      return std::min(symmetry_classes[begin->GetIndex()], symmetry_classes[end->GetIndex()]);
    }
    return std::numeric_limits<unsigned int>::max();
  }

struct IndexClassPair
{
  IndexClassPair(unsigned int _setIndex, unsigned int _symClass) : setIndex(_setIndex), symClass(_symClass) {}
  unsigned int setIndex;
  unsigned int symClass;
};

bool CompareIndexClassPair(const IndexClassPair &pair1, const IndexClassPair &pair2)
{
  return pair1.symClass < pair2.symClass;
}

std::vector<StereogenicUnit> orderSetBySymmetryClasses(OBMol *mol, const std::vector<StereogenicUnit> &set,
    const std::vector<unsigned int> &symmetry_classes)
{
  std::vector<IndexClassPair> pairs;
  for (unsigned int i = 0; i < set.size(); ++i) {
    const StereogenicUnit &unit = set[i];
    pairs.push_back(IndexClassPair(i, getSymClass(unit, mol, symmetry_classes)));
  }

  std::sort(pairs.begin(), pairs.end(), CompareIndexClassPair);
  
  std::vector<StereogenicUnit> ordered;
  for (unsigned int i = 0; i < pairs.size(); ++i) {
    ordered.push_back(set[pairs[i].setIndex]);
  }
  return ordered;
}


  bool shareSymClass(std::vector<StereogenicUnit> &set1, const std::vector<StereogenicUnit> &set2, OBMol *mol,
      const std::vector<unsigned int> &symmetry_classes)
  {
    if (!set1.size() || !set2.size())
      return false;
    
    std::vector<unsigned int> symClasses;
    
    for (unsigned int i = 0; i < set1.size(); ++i) {
      if (set1[i].type == OBStereo::Tetrahedral) {
        OBAtom *center = mol->GetAtomById(set1[i].id);
        symClasses.push_back(symmetry_classes[center->GetIndex()]);
      }
    }
    for (unsigned int i = 0; i < set2.size(); ++i) {
      if (set2[i].type == OBStereo::Tetrahedral) {
        OBAtom *center = mol->GetAtomById(set2[i].id);
        symClasses.push_back(symmetry_classes[center->GetIndex()]);
      }
    }

    std::sort(symClasses.begin(), symClasses.end());


    unsigned int symClassCount = std::unique(symClasses.begin(), symClasses.end()) - symClasses.begin();
    cout << "symClassCount = " << symClassCount << endl;
    return (symClassCount == 1);
  }

  bool shareUnit(std::vector<StereogenicUnit> &set1, const std::vector<StereogenicUnit> &set2)
  {
    for (unsigned int i = 0; i < set2.size(); ++i)
      for (unsigned int j = 0; j  < set1.size(); ++j)
        if ((set1[j].type == set2[i].type) && (set1[j].id == set2[i].id))
          return true;
    return false;
  }

  void merge(std::vector<StereogenicUnit> &set1, const std::vector<StereogenicUnit> &set2)
  {
    for (unsigned int i = 0; i < set2.size(); ++i) {
      bool found = false;
      for (unsigned int j = 0; j  < set1.size(); ++j) {
        if ((set1[j].type == set2[i].type) && (set1[j].id == set2[i].id))
          found = true;
      }

      if (!found)
        set1.push_back(set2[i]);
    } 
  }

  void mergeSets(std::vector<std::vector<StereogenicUnit> > &mergedSets, const std::vector<StereogenicUnit> &set,
      OBMol *mol, const std::vector<unsigned int> &symmetry_classes)
  {
    for (unsigned int i = 0; i  < mergedSets.size(); ++i) {
      if (!shareUnit(mergedSets[i], set))
        continue;
      if (!shareSymClass(mergedSets[i], set, mol, symmetry_classes))
        continue;
      merge(mergedSets[i], set);
    }   
  }

  bool isSubSet(std::vector<StereogenicUnit> &set1, const std::vector<StereogenicUnit> &set2)
  {
    for (unsigned int i = 0; i < set2.size(); ++i) {
      bool found = false;
      for (unsigned int j = 0; j  < set1.size(); ++j) {
        if ((set1[j].type == set2[i].type) && (set1[j].id == set2[i].id))
          found = true;
      }

      if (!found)
        return false;
    }
    return true;  
  }
  
  void removeDuplicates(std::vector<std::vector<StereogenicUnit> > &mergedSets,
      OBMol *mol, const std::vector<unsigned int> &symmetry_classes)
  {
    std::vector<std::vector<StereogenicUnit> > result;

    // remove real duplicates
    std::vector<unsigned int> doneIndexes;
    for (unsigned int i = 0; i  < mergedSets.size(); ++i) {
      if (std::find(doneIndexes.begin(), doneIndexes.end(), i) != doneIndexes.end())
        continue;
      
      for (unsigned int j = 0; j  < mergedSets.size(); ++j) {
        if (mergedSets[i].size() != mergedSets[j].size())
          continue;
        if (isSubSet(mergedSets[i], mergedSets[j]))
          doneIndexes.push_back(j);
      }

      result.push_back(mergedSets[i]);
    }
    mergedSets = result;

    // remove subsets
    result.clear();
    doneIndexes.clear();
    for (unsigned int i = 0; i  < mergedSets.size(); ++i) {
      if (std::find(doneIndexes.begin(), doneIndexes.end(), i) != doneIndexes.end())
        continue;
 
      bool subSet = false;
      for (unsigned int j = 0; j  < mergedSets.size(); ++j) {
        if (i == j)
          continue;
        if (isSubSet(mergedSets[j], mergedSets[i]))
          subSet = true;
      }

      if (!subSet)
        result.push_back(mergedSets[i]);
    }
    mergedSets = result;

  }

  std::vector<std::vector<StereogenicUnit> > FindInterdependentStereogenicUnits(OBMol *mol,
      const std::vector<StereogenicUnit> &units, const std::vector<unsigned int> &symClasses, 
      const OBPermutationGroup &automorphisms)
  {
    typedef std::vector<StereogenicUnit>::const_iterator UnitIter;
      
    std::vector< std::vector<StereogenicUnit> > sets;

    StereoInverted inverted = StereoInverted::compute(mol, symClasses, automorphisms);
   

    //DepictAutomorphisms(mol, automorphisms, "img/automorphism");
    // Order the stereogenic units by symmetry class
    std::vector<StereogenicUnit> ordered = orderSetBySymmetryClasses(mol, units, symClasses);

    for (UnitIter unit = ordered.begin(); unit != ordered.end(); ++unit) {
      if (unit->type != OBStereo::Tetrahedral)
        continue;
      OBAtom *atom = mol->GetAtomById(unit->id);
      cout << "UNIT: id = " << atom->GetId() << "  nbr_idxs = ";
      FOR_NBORS_OF_ATOM (nbr, atom)
        cout << nbr->GetIdx() << " ";
      cout << "  nbr_symmetry_classes = ";
      FOR_NBORS_OF_ATOM (nbr, atom)
        cout << symClasses[nbr->GetIndex()] << " ";
      cout << endl;
    }

    cout << "ordered.size = " << ordered.size() << endl;
    unsigned int lastSymClass = std::numeric_limits<unsigned int>::max();
    std::vector<unsigned long> doneAtoms, doneBonds, currentAtoms, currentBonds;
    for (UnitIter unit = ordered.begin(); unit != ordered.end(); ++unit) {
      // Some logic to delay storing done stereogenic units. Once all
      // units of the same symmetry class are done, the currentXX lists are
      // moved to the doneXX lists.
      unsigned int currentSymClass = getSymClass(*unit, mol, symClasses);
      if (currentSymClass != lastSymClass) {
        for (unsigned int i = 0;  i < currentAtoms.size(); ++i)
          doneAtoms.push_back(currentAtoms[i]);
        currentAtoms.clear();
      }
      lastSymClass = currentSymClass;

      // Skip this stereounit if it is already assigned to a set
      if (unit->type == OBStereo::Tetrahedral)
        if (std::find(doneAtoms.begin(), doneAtoms.end(), unit->id) != doneAtoms.end())
          continue;
      if (unit->type == OBStereo::CisTrans)
        if (std::find(doneBonds.begin(), doneBonds.end(), unit->id) != doneBonds.end())
          continue;
      
      cout << "Checking unit.id = " << unit->id << endl;
      // Select all permutations which invert the current unit
      std::vector<StereoInverted::Entry> selection;
      unsigned int minNumInversions = units.size();

      for (unsigned int i = 0; i < automorphisms.Size(); ++i) {
        const std::vector<OBAtom*> &atoms = inverted.list[i].invertedAtoms;
        const std::vector<OBBond*> &bonds = inverted.list[i].invertedBonds;
        unsigned int numInversions = atoms.size() + bonds.size();
        if (unit->type == OBStereo::Tetrahedral) {
          if (std::find(atoms.begin(), atoms.end(), mol->GetAtomById(unit->id)) != atoms.end()) {
            selection.push_back(inverted.list[i]);
            if (numInversions < minNumInversions)
              minNumInversions = numInversions;
          }
        } else
        if (unit->type == OBStereo::CisTrans) {
          if (std::find(bonds.begin(), bonds.end(), mol->GetBondById(unit->id)) != bonds.end()) {
            selection.push_back(inverted.list[i]);
            if (numInversions < minNumInversions)
              minNumInversions = numInversions;
          }
        }
      }

      // If the first selection is empty, there are no automorphisms casuing
      // inversion of configuration for the stereogenic unit and it belongs
      // to a singleton set.
      if (selection.empty()) {
        std::vector<StereogenicUnit> set;
        set.push_back(*unit);
        sets.push_back(set);
        if (unit->type == OBStereo::Tetrahedral)
          currentAtoms.push_back(unit->id);
        if (unit->type == OBStereo::CisTrans)
          currentBonds.push_back(unit->id);
        continue;
      }

      /*
      cout << "  First selection" << endl;
      for (unsigned int i = 0; i < selection.size(); ++i) {
        const StereoInverted::Entry &entry = selection[i];
        cout << "    selection " << i+1 << endl;
        cout << "      invertedAtoms: ";
        for (unsigned int l = 0; l < entry.invertedAtoms.size(); ++l)
          cout << entry.invertedAtoms[l]->GetId() << " ";
        cout << endl;
        cout << "      invertedBonds: ";
        for (unsigned int l = 0; l < entry.invertedBonds.size(); ++l)
          cout << entry.invertedBonds[l]->GetId() << " ";
        cout << endl;
      }
      */ 

      // Select the permutations which cause the minimal number of inversions
      std::vector<StereoInverted::Entry> finalSelection;
      for (unsigned int i = 0; i < selection.size(); ++i) {
        const StereoInverted::Entry &entry = selection[i];
        const std::vector<OBAtom*> &atoms = entry.invertedAtoms;
        const std::vector<OBBond*> &bonds = entry.invertedBonds;
        unsigned int numInversions = atoms.size() + bonds.size();
        if (numInversions == minNumInversions) {
          finalSelection.push_back(selection[i]);
        }
      }

      cout << "minNumInversions = " << minNumInversions << endl;
      cout << "  Final selection" << endl;
      for (unsigned int i = 0; i < finalSelection.size(); ++i) {
        const StereoInverted::Entry &entry = finalSelection[i];
        cout << "    selection " << i+1 << endl;
        cout << "      invertedAtoms: ";
        for (unsigned int l = 0; l < entry.invertedAtoms.size(); ++l)
          cout << entry.invertedAtoms[l]->GetId() << " ";
        cout << endl;
        cout << "      invertedBonds: ";
        for (unsigned int l = 0; l < entry.invertedBonds.size(); ++l)
          cout << entry.invertedBonds[l]->GetId() << " ";
        cout << endl;
      }
 
      std::vector<unsigned long> setAtoms, setBonds;
      std::vector<StereogenicUnit> set;
      for (unsigned int i = 0; i < finalSelection.size(); ++i) {
        const std::vector<OBAtom*> &atoms = finalSelection[i].invertedAtoms;
        for (unsigned int j = 0; j < atoms.size(); ++j) {
          if (std::find(doneAtoms.begin(), doneAtoms.end(), atoms[j]->GetId()) != doneAtoms.end())
            continue;
          if (std::find(setAtoms.begin(), setAtoms.end(), atoms[j]->GetId()) != setAtoms.end())
            continue;
          //if (std::find(currentAtoms.begin(), currentAtoms.end(), atoms[j]->GetId()) != currentAtoms.end())
          //  continue;
          set.push_back(StereogenicUnit(OBStereo::Tetrahedral, atoms[j]->GetId()));
          currentAtoms.push_back(atoms[j]->GetId());
          setAtoms.push_back(atoms[j]->GetId());
        }
        const std::vector<OBBond*> &bonds = finalSelection[i].invertedBonds;
        for (unsigned int j = 0; j < bonds.size(); ++j) {
          if (std::find(doneBonds.begin(), doneBonds.end(), bonds[j]->GetId()) != doneBonds.end())
            continue;
          set.push_back(StereogenicUnit(OBStereo::CisTrans, bonds[j]->GetId()));
          doneBonds.push_back(bonds[j]->GetId());
        }
      }
      sets.push_back(set);
 
      cout << "SET [ ";
      for (unsigned int k = 0; k < set.size(); ++k) {
        cout << set[k].id << " ";
      }
      cout << "]" << endl;
    }
    
    cout << "BEFORE SETS: ";
    for (unsigned int i = 0; i  < sets.size(); ++i) {
      cout << "[ ";
      const std::vector<StereogenicUnit> &set = sets[i];
      for (unsigned int k = 0; k < set.size(); ++k) {
        cout << set[k].id << " ";
      }
      cout << "]  ";
    }
    cout << endl;
   
    FOR_RINGS_OF_MOL (ring, mol) {
      if (ring->_path.size() != 4)
        continue;
      bool allSameSymClass = true;
      unsigned int symClass = symClasses[ring->_path[0] - 1];
      for (unsigned int i = 0; i < ring->_path.size(); ++i) {
        if (symClass != symClasses[ring->_path[i] - 1])
          allSameSymClass = false;
        if (!isTetrahedral(mol->GetAtom(ring->_path[i]), units))
          allSameSymClass = false;
      }
      if (!allSameSymClass)
        continue;

      std::vector<StereogenicUnit> set;
      for (unsigned int i = 0; i < ring->_path.size(); ++i)
        set.push_back(StereogenicUnit(OBStereo::Tetrahedral, mol->GetAtom(ring->_path[i])->GetId(), true));
      sets.push_back(set);
    }

    // Merge overlapping sets (with same symmetry classes!)
    for (unsigned int i = 0; i  < sets.size(); ++i) {
      const std::vector<StereogenicUnit> &iSet = sets[i];
      mergeSets(sets, iSet, mol, symClasses);
    }
    // Remove duplicated sets and subsets
    removeDuplicates(sets, mol, symClasses);

    cout << "FINAL SETS: ";
    for (unsigned int i = 0; i  < sets.size(); ++i) {
      cout << "[ ";
      const std::vector<StereogenicUnit> &set = sets[i];
      for (unsigned int k = 0; k < set.size(); ++k) {
        cout << set[k].id << " ";
      }
      cout << "]  ";
    }
    cout << endl;

    return sets;  
  }
 
  /**
   * Perform symmetry analysis.
   *
   * @return vector containing symmetry classes index by OBAtom::GetIndex().
   */
  std::vector<unsigned int> FindSymmetry(OBMol *mol)
  {
    OBGraphSym symmetry(mol);
    std::vector<unsigned int> symClasses;
    symmetry.GetSymmetry(symClasses);
    return symClasses;
  }

  ////////////////////////////////////////////////////////////////////////////
  //
  //  From0D
  //
  ////////////////////////////////////////////////////////////////////////////

  void StereoFrom0D(OBMol *mol)
  {
    if (mol->HasChiralityPerceived())
      return;
     
    obErrorLog.ThrowError(__FUNCTION__, "Ran OpenBabel::StereoFrom0D", obAuditMsg);

    std::vector<unsigned int> symClasses = FindSymmetry(mol);
    std::vector<StereogenicUnit> stereogenicUnits = FindStereogenicUnits(mol, symClasses);
    TetrahedralFrom0D(mol, stereogenicUnits);
    CisTransFrom0D(mol, stereogenicUnits);
    mol->SetChiralityPerceived();
  }

  std::vector<OBTetrahedralStereo*> TetrahedralFrom0D(OBMol *mol, 
      const std::vector<StereogenicUnit> &stereoUnits, bool addToMol)
  {
    std::vector<OBTetrahedralStereo*> configs;
    obErrorLog.ThrowError(__FUNCTION__, "Ran OpenBabel::TetrahedralFrom0D", obAuditMsg);

    // Delete any existing stereo objects that are not a member of 'centers'
    // and make a map of the remaining ones
    std::map<unsigned long, OBTetrahedralStereo*> existingMap;
    std::vector<OBGenericData*>::iterator data;
    std::vector<OBGenericData*> stereoData = mol->GetAllData(OBGenericDataType::StereoData);
    for (data = stereoData.begin(); data != stereoData.end(); ++data) {
      if (static_cast<OBStereoBase*>(*data)->GetType() == OBStereo::Tetrahedral) {
        OBTetrahedralStereo *ts = dynamic_cast<OBTetrahedralStereo*>(*data);
        unsigned long center = ts->GetConfig().center;
        // check if the center is really stereogenic
        bool isStereogenic = false;
        std::vector<StereogenicUnit>::const_iterator u;
        for (u = stereoUnits.begin(); u != stereoUnits.end(); ++u) {
          if ((*u).type == OBStereo::Tetrahedral)
            if ((*u).id == center)
              isStereogenic = true;
        }

        if (isStereogenic) {
          existingMap[center] = ts;
          configs.push_back(ts);
        } else {
          // According to OpenBabel, this is not a tetrahedral stereo
          obErrorLog.ThrowError(__FUNCTION__, "Removed spurious TetrahedralStereo object", obAuditMsg);
          mol->DeleteData(ts);
        }
      }
    }

    std::vector<StereogenicUnit>::const_iterator u;
    for (u = stereoUnits.begin(); u != stereoUnits.end(); ++u) {
      // skip non-tetrahedral units
      if ((*u).type != OBStereo::Tetrahedral)
        continue;
      // if there already exists a OBTetrahedralStereo object for this 
      // center, continue
      if (existingMap.find((*u).id) != existingMap.end())
        continue;

      OBAtom *center = mol->GetAtomById((*u).id);
 
      OBTetrahedralStereo::Config config;
      config.specified = false;
      config.center = (*u).id;
      FOR_NBORS_OF_ATOM(nbr, center) {
        if (config.from == OBStereo::NoRef)
          config.from = nbr->GetId();
        else
          config.refs.push_back(nbr->GetId());
      }

      if ((config.refs.size() == 2))
        config.refs.push_back(OBStereo::ImplicitRef); // need to add largest number on end to work

      OBTetrahedralStereo *th = new OBTetrahedralStereo(mol);
      th->SetConfig(config);
      
      configs.push_back(th);
      // add the data to the molecule if needed
      if (addToMol)
        mol->SetData(th);
    }

    return configs;    
  }

  std::vector<OBCisTransStereo*> CisTransFrom0D(OBMol *mol, 
      const std::vector<StereogenicUnit> &stereoUnits, 
      bool addToMol)
  {
    std::vector<OBCisTransStereo*> configs;
    obErrorLog.ThrowError(__FUNCTION__, "Ran OpenBabel::CisTransFrom0D", obAuditMsg);
 
    std::vector<unsigned long> bonds;
    for (std::vector<StereogenicUnit>::const_iterator u = stereoUnits.begin(); u != stereoUnits.end(); ++u)
      if ((*u).type == OBStereo::CisTrans)
        bonds.push_back((*u).id);

    // Delete any existing stereo objects that are not a member of 'bonds'
    // and make a map of the remaining ones
    std::map<unsigned long, OBCisTransStereo*> existingMap;
    std::vector<OBGenericData*>::iterator data;
    std::vector<OBGenericData*> stereoData = mol->GetAllData(OBGenericDataType::StereoData);
    for (data = stereoData.begin(); data != stereoData.end(); ++data) {
      if (static_cast<OBStereoBase*>(*data)->GetType() == OBStereo::CisTrans) {
        OBCisTransStereo *ct = dynamic_cast<OBCisTransStereo*>(*data);
        OBCisTransStereo::Config config = ct->GetConfig();
        // find the bond id from begin & end atom ids
        unsigned long id = OBStereo::NoRef;
        OBAtom *a = mol->GetAtomById(config.begin);
        if (!a)
          continue;
        FOR_BONDS_OF_ATOM (bond, a) {
          unsigned long beginId = bond->GetBeginAtom()->GetId();
          unsigned long endId = bond->GetEndAtom()->GetId();
          if ((beginId == config.begin && endId == config.end) ||
              (beginId == config.end && endId == config.begin)) {
            id = bond->GetId();
            break;
          }
        }

        if (std::find(bonds.begin(), bonds.end(), id) == bonds.end()) {
          // According to OpenBabel, this is not a cis trans stereo
          obErrorLog.ThrowError(__FUNCTION__, "Removed spurious CisTransStereo object", obAuditMsg);
          mol->DeleteData(ct);
        }
        else {
          existingMap[id] = ct;
          configs.push_back(ct);
        }
      }
    }

    std::vector<unsigned long>::iterator i;
    for (i = bonds.begin(); i != bonds.end(); ++i) {
      // If there already exists a OBCisTransStereo object for this 
      // bond, leave it alone
      if (existingMap.find(*i) != existingMap.end())
        continue;

      OBBond *bond = mol->GetBondById(*i);
      OBAtom *begin = bond->GetBeginAtom();
      OBAtom *end = bond->GetEndAtom();

      OBCisTransStereo::Config config;
      config.specified = false;
      // begin
      config.begin = begin->GetId();
      FOR_NBORS_OF_ATOM (nbr, begin) {
        if (nbr->GetId() == end->GetId())
          continue;
        config.refs.push_back(nbr->GetId());
      }
      if (config.refs.size() == 1) {
        config.refs.push_back(OBStereo::ImplicitRef);
      }
      // end
      config.end = end->GetId();
      FOR_NBORS_OF_ATOM (nbr, end) {
        if (nbr->GetId() == begin->GetId())
          continue;
        config.refs.push_back(nbr->GetId());
      }
      if (config.refs.size() == 3) {
        config.refs.push_back(OBStereo::ImplicitRef);
      }

      OBCisTransStereo *ct = new OBCisTransStereo(mol);
      ct->SetConfig(config);
      
      configs.push_back(ct);
      // add the data to the molecule if needed
      if (addToMol)
        mol->SetData(ct);
    }

    return configs;
  }



  ////////////////////////////////////////////////////////////////////////////
  //
  //  From3D
  //
  ////////////////////////////////////////////////////////////////////////////

  void StereoFrom3D(OBMol *mol, bool force)
  {
    if (mol->HasChiralityPerceived() && !force)
      return;
     
    obErrorLog.ThrowError(__FUNCTION__, "Ran OpenBabel::StereoFrom3D", obAuditMsg);

    mol->DeleteData(OBGenericDataType::StereoData);
    std::vector<unsigned int> symClasses = FindSymmetry(mol);
    std::vector<StereogenicUnit> stereogenicUnits = FindStereogenicUnits(mol, symClasses);
    TetrahedralFrom3D(mol, stereogenicUnits);
    CisTransFrom3D(mol, stereogenicUnits);
    mol->SetChiralityPerceived();
  }

  //! Calculate the "sign of a volume" given by a set of 4 coordinates
  double VolumeSign(const vector3 &a, const vector3 &b, const vector3 &c, const vector3 &d)
  {
    vector3 A, B, C;
    A = b - a;
    B = c - a;
    C = d - a;
    matrix3x3 m(A, B, C);
    return m.determinant();
  }

  std::vector<OBTetrahedralStereo*> TetrahedralFrom3D(OBMol *mol, 
      const std::vector<StereogenicUnit> &stereoUnits, bool addToMol)
  {
    std::vector<OBTetrahedralStereo*> configs;
    obErrorLog.ThrowError(__FUNCTION__, "Ran OpenBabel::TetrahedralFrom3D", obAuditMsg);

    // find all tetrahedral centers
    std::vector<unsigned long> centers;
    for (std::vector<StereogenicUnit>::const_iterator u = stereoUnits.begin(); u != stereoUnits.end(); ++u)
      if ((*u).type == OBStereo::Tetrahedral)
        centers.push_back((*u).id);
     
    std::vector<unsigned long>::iterator i;
    for (i = centers.begin(); i != centers.end(); ++i) {
      OBAtom *center = mol->GetAtomById(*i);
 
      // make sure we have at least 3 heavy atom neighbors
      // timvdm 28 Jun 2009: This is already checked in FindTetrahedralAtoms
      if (center->GetHvyValence() < 3) {
        std::stringstream errorMsg;
        errorMsg << "Cannot calculate a signed volume for an atom with a heavy atom valence of " 
                 << center->GetHvyValence() << std::endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
        continue;
      }
        
      OBTetrahedralStereo::Config config;
      config.center = *i;
      FOR_NBORS_OF_ATOM(nbr, center) {
        if (config.from == OBStereo::NoRef)
          config.from = nbr->GetId();
        else
          config.refs.push_back(nbr->GetId());
      }

      bool use_central_atom = false;
           
      // Create a vector with the coordinates of the neighbor atoms
      // and check for a bond that indicates unspecified stereochemistry
      std::vector<vector3> nbrCoords;
      OBAtom *from = mol->GetAtomById(config.from);
      OBBond *bond = mol->GetBond(from, center);
      if (bond->IsWedgeOrHash() && bond->GetBeginAtom()==center)
        config.specified = false;

      nbrCoords.push_back(from->GetVector());
      for (OBStereo::RefIter id = config.refs.begin(); id != config.refs.end(); ++id) {
        OBAtom *nbr = mol->GetAtomById(*id);
        nbrCoords.push_back(nbr->GetVector());
        OBBond *bond = mol->GetBond(nbr, center);
        if (bond->IsWedgeOrHash() && bond->GetBeginAtom()==center)
          config.specified = false;
      }
    
        // Checks for a neighbour having 0 co-ords (added hydrogen etc)
        /* FIXME: needed? if the molecule has 3D coords, additional 
         * hydrogens will get coords using OBAtom::GetNewBondVector
        for (std::vector<vector3>::iterator coord = nbrCoords.begin(); coord != nbrCoords.end(); ++coord) { 
          // are the coordinates zero to 6 or more significant figures
          if (coord->IsApprox(VZero, 1.0e-6)) {
            if (!use_central_atom) {
              use_central_atom = true;
            } else {
              obErrorLog.ThrowError(__FUNCTION__, 
                  "More than 2 neighbours have 0 co-ords when attempting 3D chiral calculation", obInfo);
            }
          }
        }
        */

      // If we have three heavy atoms we can use the chiral center atom itself for the fourth
      // will always give same sign (for tetrahedron), magnitude will be smaller.
      if ((config.refs.size() == 2) || use_central_atom) {
        nbrCoords.push_back(center->GetVector());
        config.refs.push_back(OBStereo::ImplicitRef); // need to add largest number on end to work
      }

      double sign = VolumeSign(nbrCoords[0], nbrCoords[1], nbrCoords[2], nbrCoords[3]);
      if (sign < 0.0)
        config.winding = OBStereo::AntiClockwise;

      OBTetrahedralStereo *th = new OBTetrahedralStereo(mol);
      th->SetConfig(config);
      
      configs.push_back(th);
      // add the data to the molecule if needed
      if (addToMol)
        mol->SetData(th);
    }

    return configs;    
  }

  std::vector<OBCisTransStereo*> CisTransFrom3D(OBMol *mol, 
      const std::vector<StereogenicUnit> &stereoUnits, bool addToMol)
  {
    std::vector<OBCisTransStereo*> configs;
    obErrorLog.ThrowError(__FUNCTION__, "Ran OpenBabel::CisTransFrom3D", obAuditMsg);
 
    // find all cis/trans bonds
    std::vector<unsigned long> bonds;
    for (std::vector<StereogenicUnit>::const_iterator u = stereoUnits.begin(); u != stereoUnits.end(); ++u)
      if ((*u).type == OBStereo::CisTrans)
        bonds.push_back((*u).id);

    std::vector<unsigned long>::iterator i;
    for (i = bonds.begin(); i != bonds.end(); ++i) {
      OBBond *bond = mol->GetBondById(*i);
      OBAtom *begin = bond->GetBeginAtom();
      OBAtom *end = bond->GetEndAtom();

      // Create a vector with the coordinates of the neighbor atoms
      std::vector<vector3> bondVecs;
      OBCisTransStereo::Config config;
      // begin
      config.begin = begin->GetId();
      FOR_NBORS_OF_ATOM (nbr, begin) {
        if (nbr->GetId() == end->GetId())
          continue;
        config.refs.push_back(nbr->GetId());
        bondVecs.push_back(nbr->GetVector() - begin->GetVector());
      }
      if (config.refs.size() == 1) {
        config.refs.push_back(OBStereo::ImplicitRef);
        vector3 pos;
        mol->GetAtomById(config.refs.at(0))->GetNewBondVector(pos, 1.0);
        bondVecs.push_back(pos - begin->GetVector());
      }
      // end
      config.end = end->GetId();
      FOR_NBORS_OF_ATOM (nbr, end) {
        if (nbr->GetId() == begin->GetId())
          continue;
        config.refs.push_back(nbr->GetId());
        bondVecs.push_back(nbr->GetVector() - end->GetVector());
      }
      if (config.refs.size() == 3) {
        config.refs.push_back(OBStereo::ImplicitRef);
        vector3 pos;
        mol->GetAtomById(config.refs.at(2))->GetNewBondVector(pos, 1.0);
        bondVecs.push_back(pos - end->GetVector());
      }

      // 0      3       0   3
      //  \    /         \ /      angle 0-*-3 & 1-*-2: 60 degrees (cis)
      //   C==C    -->    *       angle 0-*-1 & 2-*-3: 120 degrees (same bond atom)
      //  /    \         / \      angle 0-*-2 & 1-*-3: 180 degrees (trans)
      // 1      2       1   2
      double angle = vectorAngle(bondVecs[0], bondVecs[2]);
      if (IsNear(angle, 60.0, 10.0))
        config.shape = OBStereo::ShapeZ;
      if (IsNear(angle, 180.0, 10.0))
        config.shape = OBStereo::ShapeU;

      OBCisTransStereo *ct = new OBCisTransStereo(mol);
      ct->SetConfig(config);
      
      configs.push_back(ct);
      // add the data to the molecule if needed
      if (addToMol)
        mol->SetData(ct);
    }

    return configs;
  }

  ////////////////////////////////////////////////////////////////////////////
  //
  //  From2D
  //
  //  Reference:
  //  [1] T. Cieplak, J.L. Wisniewski, A New Effective Algorithm for the 
  //  Unambiguous Identification of the Stereochemical Characteristics of 
  //  Compounds During Their Registration in Databases. Molecules 2000, 6,
  //  915-926, http://www.mdpi.org/molecules/papers/61100915/61100915.htm
  ////////////////////////////////////////////////////////////////////////////

  void StereoFrom2D(OBMol *mol, std::map<OBBond*, enum OBStereo::BondDirection> *updown, bool force)
  {
    if (mol->HasChiralityPerceived() && !force)
      return;
      
    obErrorLog.ThrowError(__FUNCTION__, "Ran OpenBabel::StereoFrom2D", obAuditMsg);

    std::vector<unsigned int> symClasses = FindSymmetry(mol);
    std::vector<StereogenicUnit> stereogenicUnits = FindStereogenicUnits(mol, symClasses);
    mol->DeleteData(OBGenericDataType::StereoData);
    TetrahedralFrom2D(mol, stereogenicUnits);
    CisTransFrom2D(mol, stereogenicUnits, updown);
    mol->SetChiralityPerceived();
  }
 
  //! Calculate the "sign of a triangle" given by a set of 3 2D coordinates
  double TriangleSign(const vector3 &a, const vector3 &b, const vector3 &c)
  {
    // equation 6 from [1]
    return (a.x() - c.x()) * (b.y() - c.y()) - (a.y() - c.y()) * (b.x() - c.x());
  }

  std::vector<OBTetrahedralStereo*> TetrahedralFrom2D(OBMol *mol, 
      const std::vector<StereogenicUnit> &stereoUnits, bool addToMol)
  {
    std::vector<OBTetrahedralStereo*> configs;
    obErrorLog.ThrowError(__FUNCTION__, "Ran OpenBabel::TetrahedralFrom2D", obAuditMsg);
 
    // find all tetrahedral centers
    std::vector<unsigned long> centers;
    for (std::vector<StereogenicUnit>::const_iterator u = stereoUnits.begin(); u != stereoUnits.end(); ++u)
      if ((*u).type == OBStereo::Tetrahedral)
        centers.push_back((*u).id);
 
      
    std::vector<unsigned long>::iterator i;
    for (i = centers.begin(); i != centers.end(); ++i) {
      OBAtom *center = mol->GetAtomById(*i);
 
      // make sure we have at least 3 heavy atom neighbors
      if (center->GetHvyValence() < 3) {
        std::stringstream errorMsg;
        errorMsg << "Cannot calculate a signed volume for an atom with a heavy atom valence of " 
                 << center->GetHvyValence() << std::endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);
        continue;
      }
        
      OBTetrahedralStereo::Config config;
      config.center = *i;
        
      // find the hash, wedge and 2 plane atoms
      std::vector<OBAtom*> planeAtoms;
      std::vector<OBAtom*> wedgeAtoms;
      std::vector<OBAtom*> hashAtoms;
      FOR_BONDS_OF_ATOM(bond, center) {
        OBAtom *nbr = bond->GetNbrAtom(center);
        // hash bonds
        if (bond->IsHash()) {
          if (bond->GetBeginAtom()->GetId() == center->GetId()) {
            // this is a 'real' hash bond going from center to nbr
            hashAtoms.push_back(nbr);
          } else {
            // this is an 'inverted' hash bond going from nbr to center
            wedgeAtoms.push_back(nbr);
          }
        } else if (bond->IsWedge()) {
          // wedge bonds
          if (bond->GetBeginAtom()->GetId() == center->GetId()) {
            // this is a 'real' wedge bond going from center to nbr
            wedgeAtoms.push_back(nbr);
          } else {
            // this is an 'inverted' hash bond going from nbr to center
            hashAtoms.push_back(nbr);
          }
        } else if (bond->IsWedgeOrHash()) {
          config.specified = false;
          break;
        } else { 
          // plane bonds
          planeAtoms.push_back(nbr);
        }
      }

      bool success = true;
      
      using namespace std;
      if (!config.specified) {
        // unspecified
        FOR_NBORS_OF_ATOM (nbr, center)
          if (config.from == OBStereo::NoRef)
            config.from = nbr->GetId();
          else
            config.refs.push_back(nbr->GetId());
        while (config.refs.size() < 3)
          config.refs.push_back(OBStereo::ImplicitRef);
      } else
      if (planeAtoms.size() == 2) {
        if (hashAtoms.size() == 1 && wedgeAtoms.size() == 1) {
          // plane1 + plane2, hash, wedge
          config.from = wedgeAtoms[0]->GetId();
          config.refs.resize(3);
          config.refs[0] = planeAtoms[0]->GetId();
          config.refs[1] = planeAtoms[1]->GetId();
          config.refs[2] = hashAtoms[0]->GetId();
          double sign = TriangleSign(planeAtoms[0]->GetVector(), 
              planeAtoms[1]->GetVector(), hashAtoms[0]->GetVector());
          if (sign > 0.0)
            config.winding = OBStereo::AntiClockwise;
        } else if ((hashAtoms.size() + wedgeAtoms.size()) == 1) {
          // Either: plane1 + plane2 + hash *or* plane1 + plane2 + wedge
          OBAtom* stereoAtom;
          if (hashAtoms.size() == 1) {
            config.from = OBStereo::ImplicitRef;
            config.view = OBStereo::ViewFrom;
            stereoAtom = hashAtoms[0];
          }
          else { // wedgeAtoms.size() == 1
            config.towards = OBStereo::ImplicitRef;
            config.view = OBStereo::ViewTowards;
            stereoAtom = wedgeAtoms[0];
          }
          config.refs.resize(3);
          config.refs[0] = planeAtoms[0]->GetId();
          config.refs[1] = planeAtoms[1]->GetId();
          config.refs[2] = stereoAtom->GetId();
          double sign = TriangleSign(planeAtoms[0]->GetVector(), 
              planeAtoms[1]->GetVector(), stereoAtom->GetVector());
          if (sign > 0.0)
            config.winding = OBStereo::AntiClockwise;
        } else {
          success = false;
        }
      } else if (planeAtoms.size() == 3) {
        if ( (hashAtoms.size() + wedgeAtoms.size()) == 1) {
          // Either: plane1 + plane2 + plane3, hash
          //     or: plane1 + plane2 + plane3, wedge
          if (hashAtoms.size() == 1) {
            config.towards = hashAtoms[0]->GetId();
            config.view = OBStereo::ViewTowards;
          }
          else { // wedgeAtoms.size() == 1
            config.from = wedgeAtoms[0]->GetId();
            config.view = OBStereo::ViewFrom;
          }
          config.refs.resize(3);
          config.refs[0] = planeAtoms[0]->GetId();
          config.refs[1] = planeAtoms[1]->GetId();
          config.refs[2] = planeAtoms[2]->GetId();
          double sign = TriangleSign(planeAtoms[0]->GetVector(), 
              planeAtoms[1]->GetVector(), planeAtoms[2]->GetVector());
          if (sign > 0.0)
            config.winding = OBStereo::AntiClockwise;
        } else {
          success = false;
        }
      
      } else {
        success = false;
      }

      if (!success) {
//         std::stringstream errorMsg;
//         errorMsg << "Symmetry analysis found atom with id " << center->GetId() 
//             << " to be a tetrahedral atom but the wedge/hash bonds can't be interpreted." << std::endl
//             << " # in-plane bonds = " << planeAtoms.size() << std::endl
//             << " # wedge bonds = " << wedgeAtoms.size() << std::endl
//             << " # hash bonds = " << hashAtoms.size() << std::endl
//             << std::endl;
//         obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
        continue;
      }
 

      OBTetrahedralStereo *th = new OBTetrahedralStereo(mol);
      th->SetConfig(config);

      configs.push_back(th);
      // add the data to the molecule if needed
      if (addToMol)
        mol->SetData(th);
    }
   
    return configs;
  }

  std::vector<OBCisTransStereo*> CisTransFrom2D(OBMol *mol, 
      const std::vector<StereogenicUnit> &stereoUnits, 
      const std::map<OBBond*, enum OBStereo::BondDirection> *updown, bool addToMol)
  {
    std::vector<OBCisTransStereo*> configs;
    obErrorLog.ThrowError(__FUNCTION__, "Ran OpenBabel::CisTransFrom2D", obAuditMsg);

    // find all cis/trans bonds
    std::vector<unsigned long> bonds;
    for (std::vector<StereogenicUnit>::const_iterator u = stereoUnits.begin(); u != stereoUnits.end(); ++u)
      if ((*u).type == OBStereo::CisTrans)
        bonds.push_back((*u).id);
    
    std::vector<unsigned long>::iterator i;
    for (i = bonds.begin(); i != bonds.end(); ++i) {
      OBBond *bond = mol->GetBondById(*i);
      OBAtom *begin = bond->GetBeginAtom();
      OBAtom *end = bond->GetEndAtom();

      // Create a vector with the coordinates of the neighbor atoms
      std::vector<vector3> bondVecs;
      OBCisTransStereo::Config config;
      // begin
      config.begin = begin->GetId();
      FOR_NBORS_OF_ATOM (nbr, begin) {
        if (nbr->GetId() == end->GetId())
          continue;
        config.refs.push_back(nbr->GetId());
        bondVecs.push_back(nbr->GetVector());
      }
      if (config.refs.size() == 1) {
        config.refs.push_back(OBStereo::ImplicitRef);
        vector3 pos;
        begin->GetNewBondVector(pos, 1.0);
        bondVecs.push_back(pos);
      }
      // end
      config.end = end->GetId();
      FOR_NBORS_OF_ATOM (nbr, end) {
        if (nbr->GetId() == begin->GetId())
          continue;
        config.refs.push_back(nbr->GetId());
        bondVecs.push_back(nbr->GetVector());
      }
      if (config.refs.size() == 3) {
        config.refs.push_back(OBStereo::ImplicitRef);
        vector3 pos;
        end->GetNewBondVector(pos, 1.0);
        bondVecs.push_back(pos);
      }

      config.specified = true;
      if (updown) {
        std::map<OBBond*, enum OBStereo::BondDirection>::const_iterator ud_cit;
        ud_cit = updown->find(bond);
        if (ud_cit!=updown->end() && ud_cit->second==OBStereo::UnknownDir)
            config.specified = false;
      }
      if (config.specified==true) { // Work out the stereochemistry
        // 0      3       
        //  \    /        2 triangles: 0-1-b & 2-3-a
        //   a==b    -->  same sign: U
        //  /    \        opposite sign: Z
        // 1      2       
        /*
        double sign1 = TriangleSign(begin->GetVector(), end->GetVector(), bondVecs[0]);
        double sign2 = TriangleSign(begin->GetVector(), end->GetVector(), bondVecs[2]);
        */
        double sign1 = TriangleSign(bondVecs[0], bondVecs[1], end->GetVector());
        double sign2 = TriangleSign(bondVecs[2], bondVecs[3], begin->GetVector());
        double sign = sign1 * sign2;

        if (sign < 0.0) // opposite sign
          config.shape = OBStereo::ShapeZ;
      }

      OBCisTransStereo *ct = new OBCisTransStereo(mol);
      ct->SetConfig(config);

      configs.push_back(ct);
      // add the data to the molecule if needed
      if (addToMol)
        mol->SetData(ct);
    }

    return configs;
  }

  void TetStereoTo0D(OBMol &mol, 
      std::map<OBBond*, enum OBStereo::BondDirection> &updown,
      std::map<OBBond*, OBStereo::Ref> &from)
  {
    // Store the tetcenters for the second loop (below)
    std::set <unsigned long> tetcenters;
    std::vector<OBGenericData*> vdata = mol.GetAllData(OBGenericDataType::StereoData);
    for (std::vector<OBGenericData*>::iterator data = vdata.begin(); data != vdata.end(); ++data)
      if (((OBStereoBase*)*data)->GetType() == OBStereo::Tetrahedral) {
        OBTetrahedralStereo *ts = dynamic_cast<OBTetrahedralStereo*>(*data);
        OBTetrahedralStereo::Config cfg = ts->GetConfig();
        tetcenters.insert(cfg.center);
      }

    // This loop sets one bond of each tet stereo to up or to down (2D only)
    std::set <OBBond *> alreadyset;
    for (std::vector<OBGenericData*>::iterator data = vdata.begin(); data != vdata.end(); ++data)
      if (((OBStereoBase*)*data)->GetType() == OBStereo::Tetrahedral) {
        OBTetrahedralStereo *ts = dynamic_cast<OBTetrahedralStereo*>(*data);
        OBTetrahedralStereo::Config cfg = ts->GetConfig();

        OBBond* chosen = (OBBond*) NULL;
        OBAtom* center = mol.GetAtomById(cfg.center);
        bool nottet_flag = false;
        bool acyclic_flag = false;
        // Find the best candidate bond to set to up/down
        // 1. **Should not already be set**
        // 2. Should not be connected to a 2nd tet center
        // (this is acceptable, as the wedge is only at one end, but will only confuse things)
        // 3. Preferably is not in a cycle
        // 4. Preferably is a terminal H
        FOR_BONDS_OF_ATOM(b, center) {
          if (alreadyset.find(&*b) == alreadyset.end()) {
            if (chosen==NULL) chosen = &*b;
            OBAtom* nbr = chosen->GetNbrAtom(center);
            if (tetcenters.find(nbr->GetId()) == tetcenters.end()) { // Not a tetcenter
              if (nottet_flag==false) {
                chosen = &*b;
                nottet_flag = true;
              }
              if (!b->IsInRing()) {
                if (acyclic_flag==false) {
                  chosen = &*b;
                  acyclic_flag = true;
                }
                if (nbr->IsHydrogen()) {
                  chosen = &*b;
                  break;
                }
              }
            }
          }
        }
        if (chosen==NULL) { // There is a remote possibility of this but let's worry about 99.9% of cases first
          obErrorLog.ThrowError(__FUNCTION__, 
            "Failed to set stereochemistry as unable to find an available bond", obError);
          return;
        }
        alreadyset.insert(chosen);
        
        OBStereo::BondDirection bonddir = OBStereo::UnknownDir;
        if (cfg.specified) {
          // Determine whether this bond should be set hash or wedge (or indeed unknown)
          // (Code inspired by perception.cpp, TetrahedralFrom2D: plane1 + plane2 + plane3, wedge)
          OBTetrahedralStereo::Config test_cfg = cfg;
           
          // If there is an implicit ref; let's make that the 'from' atom
          // otherwise use the atom on the chosen bond
          bool implicit = true;
          if (test_cfg.from != OBStereo::ImplicitRef) {
            OBStereo::RefIter ri = std::find(test_cfg.refs.begin(), test_cfg.refs.end(), (unsigned long) OBStereo::ImplicitRef);
            if (ri!=test_cfg.refs.end())
              test_cfg = OBTetrahedralStereo::ToConfig(test_cfg, OBStereo::ImplicitRef);
            else {
              test_cfg = OBTetrahedralStereo::ToConfig(test_cfg, chosen->GetNbrAtom(center)->GetId());
              implicit = false;
            }
          }
          // -ve sign implies clockwise
          double sign = TriangleSign(mol.GetAtomById(test_cfg.refs[0])->GetVector(), 
              mol.GetAtomById(test_cfg.refs[1])->GetVector(), mol.GetAtomById(test_cfg.refs[2])->GetVector());

          // Things are inverted from the point of view of the ImplicitH which we
          // assume to be of opposite stereochemistry to the wedge/hash
          bool useup = !implicit;
          if (sign > 0) useup = !useup;
          // Set to UpBond (filled wedge from cfg.center to chosen_nbr) or DownBond
          bonddir = useup ? OBStereo::UpBond : OBStereo::DownBond;
        }
        updown[chosen] = bonddir;
        from[chosen] = cfg.center;
      }
  }
  set<OBBond*> GetUnspecifiedCisTrans(OBMol& mol)
  {
    // Get double bonds with unspecified CisTransStereo
    set<OBBond*> unspec_ctstereo;
    std::vector<OBGenericData*> vdata = mol.GetAllData(OBGenericDataType::StereoData);
    for (std::vector<OBGenericData*>::iterator data = vdata.begin(); data != vdata.end(); ++data)
      if (((OBStereoBase*)*data)->GetType() == OBStereo::CisTrans) {
        OBCisTransStereo *ct = dynamic_cast<OBCisTransStereo*>(*data);
        OBCisTransStereo::Config cfg = ct->GetConfig();
        if (!cfg.specified) {
          OBBond* dbl_bond = mol.GetBond(mol.GetAtomById(cfg.begin), mol.GetAtomById(cfg.end));
          unspec_ctstereo.insert(dbl_bond);
        }
      }
    return unspec_ctstereo;
  }
}

