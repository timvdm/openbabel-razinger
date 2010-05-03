/* -*-C++-*-

**********************************************************************
Copyright (C) 2005-2006, eMolecules, Inc. (www.emolecules.com)
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************

+======================================================================
|
| AUTHOR: Craig A. James, eMolecules, Inc.
|
| DESCRIPTION: CANONICALIZATION OF SMILES
|
|       This is a specialized SMILES canonicalization algorithm.  Although
|       it can be applied in the standard fashion to a whole molecule, 
|       its real job is to generate canonical SMILES for fragments, or
|       "subsets", of the atoms of a molecule.
|
|       For example, consider the first three atoms of Oc1ccccc1.  With
|       a "normal" SMILES canonicalizer, you couldn't generate a SMILES
|       for Occ, because it's not a valid molecule.  However, this system
|       can do exactly that, by taking both the whole molecule (which 
|       retains the aromaticity), and a "subset" bitmap that specifies
|       which atoms are to be included in the SMILES.
|
|       Canonicalization is carried out per Weininger et al (J. Chem. 
|       Inf. Comput. Sci., Vol. 29, No. 2, 1989, pp 97-101), with some
|       modifications to handle bond symmetries not foreseen by Weininger
|       in that paper.
|
|       WARNING - KNOWN BUG: These functions make use of a bitmap vector
|       to represent a "fragment" -- a subset of the atoms in a molecule.
|       But this means the bonds of the fragment are implicit, not explicit,
|       which is incorrect.  For example, if you want to break one bond of
|       cyclehexane (C1CCCCC1), all six atoms will still be there, so the
|       "fragment" will be cyclic.  This is relevant when generating fragment
|       SMILES for ring systems where breaking a bond can reduce the number
|       of ring without removing any atoms.  We need to add a pair of bit
|       vectors, the atoms AND the bonds, to represent a fragment.  (Note
|       that this is also an ambiguity in OpenBabel itself, which represents
|       a ring as a set of atoms. This is only valid if the ring is a member
|       of a SSSR.)
+======================================================================
*/

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/graphsym.h>

#include <openbabel/stereo/cistrans.h>
#include <openbabel/stereo/tetrahedral.h>

#include <openbabel/depict/depict.h>

#include <iterator> // std::istream_iterator
#include <cassert>

using namespace std;

namespace OpenBabel {

/***************************************************************************
* FUNCTION: CompareXXX
*
* DESCRIPTION:
*       Three functions for use by the sort() method of a vector.
***************************************************************************/

  OBGraphSym::OBGraphSym(OBMol* pmol, OBBitVec* frag_atoms)
  {
    _pmol = pmol;
    if (frag_atoms) {
      _frag_atoms = frag_atoms;
    } else {
      _frag_atoms = new OBBitVec(_pmol->NumAtoms());
      FOR_ATOMS_OF_MOL(a, _pmol)
        _frag_atoms->SetBitOn(a->GetIdx());
    }
  }
 

// NOTE: Copied from OpenBabel/mol.cpp

  // Destructor
  OBGraphSym::~OBGraphSym()
  {
    // Nothing to free as we only hold pointers
  }

  const unsigned int OBGraphSym::NoSymmetryClass = 0x7FFFFFFF;

  bool OBGraphSym::CompareUnsigned(const unsigned int &a,const unsigned int &b)
  {
    return(a<b);
  }

  bool OBGraphSym::ComparePairFirst(const pair<OBAtom*,unsigned int> &a,const pair<OBAtom*,unsigned int> &b)
  {
    return(a.first->GetIdx() < b.first->GetIdx());
  }

  bool OBGraphSym::ComparePairSecond(const pair<OBAtom*,unsigned int> &a,const pair<OBAtom*,unsigned int> &b)
  {
    return(a.second < b.second);
  }

  bool OBGraphSym::CompareBondPairSecond(const pair<OBBond*,unsigned int> &a,const pair<OBBond*,unsigned int> &b)
  {
    return(a.second < b.second);
  }


/***************************************************************************
* FUNCTION: GetValence
*
* DESCRIPTION:
*       Like OBAtom::GetValence(): Counts the number of neighbors, but
*       doesn't count atoms not in the fragment.
***************************************************************************/

  unsigned int OBGraphSym::GetValence(OBAtom *atom)
  {
    unsigned int count = 0;
    OBBond *bond;
    OBAtom *nbr;

    vector<OBEdgeBase*>::iterator bi;
    for (bond = atom->BeginBond(bi); bond; bond = atom->NextBond(bi)) {
      nbr = bond->GetNbrAtom(atom);
      if (_frag_atoms->BitIsSet(nbr->GetIdx()))
        count++;
    }
    return(count);
  }

/***************************************************************************
* FUNCTION: GetHvyValence
*
* DESCRIPTION:
*       Like OBAtom::GetHvyValence(): Counts the number non-hydrogen
*       neighbors, but doesn't count atoms not in the fragment.
***************************************************************************/

  unsigned int OBGraphSym::GetHvyValence(OBAtom *atom)
  {
    unsigned int count = 0;
    OBBond *bond;
    OBAtom *nbr;

    vector<OBEdgeBase*>::iterator bi;
    for (bond = atom->BeginBond(bi); bond; bond = atom->NextBond(bi)) {
      nbr = bond->GetNbrAtom(atom);
      if (_frag_atoms->BitIsSet(nbr->GetIdx()) && !(nbr->IsHydrogen()))
        count++;
    }
    
    return(count);
  }

/***************************************************************************
* FUNCTION: GetHvyBondSum
*
* DESCRIPTION:
*       Sums the bond order over the bonds from this atom to other atoms
*       in the fragment.  Single = 1, double = 2, triple = 3, aromatic = 1.6,
*       but sum is rounded to nearest integer.
*
*       This is used for fragment symmetry perception instead of the "implicit
*       valence" used by the standard OpenBabel symmetry perception.  It
*       has the same effect, but we don't have to worry about hydrogen counts,
*       EXCEPT for aromatic N, where the difference between n and [nH] is
*       critical.
***************************************************************************/

  unsigned int OBGraphSym::GetHvyBondSum(OBAtom *atom)
  {
    float count = 0.0f;
    OBBond *bond;
    OBAtom *nbr;

    vector<OBEdgeBase*>::iterator bi;
    for (bond = atom->BeginBond(bi); bond; bond = atom->NextBond(bi)) {
      nbr = bond->GetNbrAtom(atom);
      if (_frag_atoms->BitIsSet(nbr->GetIdx()) && !(nbr->IsHydrogen())) {
        if (bond->IsSingle())        count += 1.0f;
        else if (bond->IsDouble())   count += 2.0f;
        else if (bond->IsTriple())   count += 3.0f;
        else if (bond->IsAromatic()) count += 1.6f;
      }
    }
    if (atom->GetAtomicNum() == 7 && atom->IsAromatic() && atom->GetImplicitValence() == 3) {
      count += 1;         // [nH] - add another bond
    }
    return(int(count + 0.5));     // round to nearest int
  }


/***************************************************************************
* FUNCTION: GetGTDVector
*
* DESCRIPTION:
*       
*       Calculates the graph theoretical distance of each atom.
*       Vector is indexed from zero.
*
*       NOTE: Copied from OpenBabel/mol.cpp
*
*       NOTE: "Indexed from zero" means it's one off from the atom->GetIdx()
*       that's used to index atoms inside the molecule!
*
*       NOTE: This function is hard to decipher, and seems to be misnamed.
*       A "distance" should be be between two atoms, but there's more here
*       than that.  It seems to be doing a breadth-first search to find the
*       most-distant atom from each atom, and reporting the number of steps
*       (which happens to be the graph-theoretical distance) to that atom.
*       The name "Graph Theoretical Distance" is thus misleading.
***************************************************************************/

  bool OBGraphSym::GetGTDVector(vector<int> &gtd)
{
  gtd.clear();
  gtd.resize(_pmol->NumAtoms());
  
  int gtdcount, natom;
  OBBitVec used, curr, next;
  OBAtom *atom, *atom1;
  OBBond *bond;
  vector<OBNodeBase*>::iterator ai;
  vector<OBEdgeBase*>::iterator j;

  next.Clear();

  for (atom = _pmol->BeginAtom(ai); atom; atom = _pmol->NextAtom(ai)) {

    int idx = atom->GetIdx();
    if (!_frag_atoms->BitIsOn(idx)) {     // Not in this fragment?
      gtd[idx-1] = 0;
      continue;
    }

    gtdcount = 0;
    used.Clear();curr.Clear();
    used.SetBitOn(idx);
    curr.SetBitOn(idx);

    while (!curr.IsEmpty()) {
      next.Clear();
      for (natom = curr.NextBit(-1);natom != curr.EndBit();natom = curr.NextBit(natom)) {
        atom1 = _pmol->GetAtom(natom);
        if (!_frag_atoms->BitIsOn(atom1->GetIdx()))
          continue;
        for (bond = atom1->BeginBond(j);bond;bond = atom1->NextBond(j)) {
          int nbr_idx = bond->GetNbrAtomIdx(atom1);
          if (   _frag_atoms->BitIsOn(nbr_idx)
              && !used.BitIsOn(nbr_idx)
              && !curr.BitIsOn(nbr_idx)
              && !(bond->GetNbrAtom(atom1))->IsHydrogen())
            next.SetBitOn(nbr_idx);
        }
      }
      used |= next;
      curr = next;
      gtdcount++;
    }
    gtd[idx-1] = gtdcount;
  }

  return(true);
}

/***************************************************************************
* FUNCTION: FindRingAtoms
*
* DESCRIPTION:
*       Finds all atoms that are part of a ring in the current fragment.
*       We start with the whole molecule's rings, and eliminate any that
*       have atoms not in the subset.  For the rings that are left, mark
*       each atom of the ring as a ring atom.
*
*       Returns a bit vector where TRUE means it's a ring atom.
***************************************************************************/

  void OBGraphSym::FindRingAtoms(OBBitVec &ring_atoms)
{
  vector<OBRing*> sssRings;
  vector<OBRing*>::iterator ri;

  ring_atoms.Resize(_pmol->NumAtoms());
  ring_atoms.Clear();

  sssRings = _pmol->GetSSSR();
  for (ri = sssRings.begin(); ri != sssRings.end(); ri++) {
    OBRing *ring = *ri;
    OBBitVec bvtmp = *_frag_atoms & ring->_pathset;       // intersection: fragment and ring
    if (bvtmp == ring->_pathset)                        // all ring atoms in fragment?
      ring_atoms |= ring->_pathset;                     //   yes - add this ring's atoms 
  }
}


/***************************************************************************
* FUNCTION: GetGIVector
*
* DESCRIPTION:
*       NOTE: Copied from OpenBabel/mol.cpp
*
*       Calculates a set of graph invariant indexes using the graph theoretical
*       distance, number of connected heavy atoms, aromatic boolean, ring
*       boolean, atomic number, and summation of bond orders connected to the
*       atom.
*
*       We have to recalculate which atoms are in rings by taking the fragment's
*       atoms into account when we generate the graph invarients.
*
*       Vector is indexed from zero (not one, like atom->GetIdx()).
*
*       NOTE: This may need to be extended to include the bond-invariant properties,
*       particularly the size of all rings the bond is in (from a SSSR).
***************************************************************************/

  void OBGraphSym::GetGIVector(vector<unsigned int> &vid)
{
  // Prepare the vector...
  vid.clear();
  vid.resize(_pmol->NumAtoms());

  // The "graph theoretical distance" for each atom (see comments in the function)
  vector<int> v;
  GetGTDVector(v);

  // Compute the ring atoms for this particular fragment (set of atoms)
  OBBitVec ring_atoms;
  FindRingAtoms(ring_atoms);

  int i;
  OBAtom *atom;
  vector<OBNodeBase*>::iterator ai;
  for (i=0, atom = _pmol->BeginAtom(ai); atom; atom = _pmol->NextAtom(ai)) {
    vid[i] = 0;
    if (_frag_atoms->BitIsOn(atom->GetIdx())) {
      vid[i] = 
        v[i]                                                    // 10 bits: graph-theoretical distance
        | (GetHvyValence(atom)                <<10)  //  4 bits: heavy valence
        | (((atom->IsAromatic()) ? 1 : 0)                <<14)  //  1 bit:  aromaticity
        | (((ring_atoms.BitIsOn(atom->GetIdx())) ? 1 : 0)<<15)  //  1 bit:  ring atom
        | (atom->GetAtomicNum()                          <<16)  //  7 bits: atomic number
        | (GetHvyBondSum(atom)               <<23)  //  4 bits: heavy bond sum
        | ((7 + atom->GetFormalCharge())                 <<27); //  4 bits: formal charge
    }
    i++;
  }
}

/***************************************************************************
* FUNCTION: BreakChiralTies
*
* DESCRIPTION:
*       After the achiral symmetry analysis ChiralSymmetry() is done, but
*       before the "tie breaker" step (see CanonicalLabels(), below), there
*       may be two (or more) chiral centers in the same symmetry class that
*       actually have different chirality.  This function finds such chiral
*       centers, and compares their chirality.  If it finds that two atoms
*       in the same symmetry class have the same chirality, it leaves them
*       alone, but if they have opposite chirality, it breaks the tie
*       between the two atoms.
*
*       Actually, it's more subtle than that.  Suppose there are a bunch of
*       chiral atoms in the same symmetry class.  The the class is divided
*       into two classes: All atoms with one chirality go into one class,
*       and all atoms with the opposite chirality go in the other class.
*
* INPUTS:
*       pmol                the molecule
*       frag_atoms          atoms of the molecules in this fragment
*       atom_sym_classes    vector of atom/symclass pairs
*
***************************************************************************/

  // Translate atom ids to symmetry classes for a CisTrans config struct
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
    // end
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

  // Translate atom ids to symmetry classes for a Tetrahedral config struct
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

  // Returns true if the sets of StereogenicUnits contain the same symmetry classes
  bool setsContainSameSymmetryClasses(OBMol *mol, const std::vector<StereogenicUnit> &set1,
      const std::vector<StereogenicUnit> &set2, const std::vector<unsigned int> &symmetry_classes)
  {
    for (unsigned int i = 0; i < set1.size(); ++i) {
      bool found = false;
      if (set1[i].type == OBStereo::Tetrahedral) {
        for (unsigned int j = 0; j < set2.size(); ++j) {
          if (set2[i].type != OBStereo::Tetrahedral)
            continue;

          OBAtom *iAtom = mol->GetAtomById(set1[i].id);
          OBAtom *jAtom = mol->GetAtomById(set2[j].id);
          if (symmetry_classes[iAtom->GetIndex()] == symmetry_classes[jAtom->GetIndex()]) {
            found = true;
            break;        
          }
        }
      }

      if (!found)
        return false;
    }

    return true;
  }

  // compare sets of sets by size (used for sorting below)
  bool CompareSetSizes(const std::vector< std::vector<StereogenicUnit> > &set1,
      const std::vector< std::vector<StereogenicUnit> > &set2)
  {
    return set1.size() > set2.size();
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

  // defined in src/stereo/perception.cpp
  int classifyTetrahedralNbrSymClasses(const std::vector<unsigned int> &symClasses, OBAtom *atom);
  int classifyCisTransNbrSymClasses(const std::vector<unsigned int> &symClasses, OBBond *doubleBond, OBAtom *atom);
  std::vector<StereogenicUnit> orderSetBySymmetryClasses(OBMol *mol, const std::vector<StereogenicUnit> &set,
      const std::vector<unsigned int> &symmetry_classes);
  unsigned int findDuplicatedSymmetryClass(OBAtom *atom, const std::vector<unsigned int> &symClasses);
  void findDuplicatedSymmetryClasses(OBAtom *atom, const std::vector<unsigned int> &symClasses,
      unsigned int &duplicated1, unsigned int &duplicated2);
  OBBitVec getFragment(OBAtom *atom, OBAtom *skip);


  // Part of orderSetByFragment...
  void orderSetByFragmentRecursive(std::vector<StereogenicUnit> &ordered, OBAtom *atom, OBAtom *skip, 
      const std::vector<StereogenicUnit> &set, OBBitVec &fragment)
  {
    FOR_NBORS_OF_ATOM (nbr, atom) {
      if (nbr->GetId() == skip->GetId())
        continue;
      if (!fragment.BitIsSet(nbr->GetId())) {
        // Add the stereocenter
        for (unsigned int i = 0; i < set.size(); ++i) {
          const StereogenicUnit &unit = set[i];
          if (unit.type == OBStereo::Tetrahedral) {
            if (nbr->GetId() == unit.id)
              ordered.push_back(unit);
          }
        }

        fragment.SetBitOn(nbr->GetId());
        orderSetByFragmentRecursive(ordered, &*nbr, skip, set, fragment);
      }
    }
  }

  // Order the set of StereogenicUnits by traversing the fragment
  std::vector<StereogenicUnit> orderSetByFragment(OBAtom *atom, OBAtom *skip, const std::vector<StereogenicUnit> &set)
  {
    std::vector<StereogenicUnit> ordered;
    OBBitVec fragment;

    for (unsigned int i = 0; i < set.size(); ++i) {
      const StereogenicUnit &unit = set[i];
      if (unit.type == OBStereo::Tetrahedral) {
        if (atom->GetId() == unit.id)
          ordered.push_back(unit);
      }
    }

    fragment.SetBitOn(atom->GetId());
    orderSetByFragmentRecursive(ordered, atom, skip, set, fragment);

    return ordered;
  }

void orderSetByRingRecursive(OBMol *mol, const OBBitVec &fragment, OBAtom *atom, 
    const std::vector<StereogenicUnit> &set, std::vector<StereogenicUnit> &result)
{
  FOR_NBORS_OF_ATOM (nbr, atom) {
    if (fragment.BitIsSet(nbr->GetId())) {
      bool alreadyAdded = false;
      for (unsigned int i = 0; i < result.size(); ++i)
        if (nbr->GetId() == result[i].id)
          alreadyAdded = true;

      if (!alreadyAdded) {
        for (unsigned int i = 0; i < set.size(); ++i)
          if (nbr->GetId() == set[i].id)
            result.push_back(set[i]);
        orderSetByRingRecursive(mol, fragment, &*nbr, set, result);
      }
    }
  }
}

std::vector<StereogenicUnit> orderSetByRing(OBMol *mol, const OBBitVec &fragment, const std::vector<StereogenicUnit> &set)
{
  std::vector<StereogenicUnit> result;
  result.push_back(set[0]);
  orderSetByRingRecursive(mol, fragment, mol->GetAtomById(set[0].id), set, result);
  return result;
}

  
  int findDescriptor(const OBTetrahedralStereo::Config &config)
  {
    std::vector<unsigned long> refs = config.refs;
    refs.insert(refs.begin(), config.from);
    if (OBStereo::NumInversions(refs) % 2)
      return 1;
    else
      return 0; 
  }
                
std::vector<int> findDescriptorVector(OBMol *mol, const OBBitVec &fragment,
    const std::vector<StereogenicUnit> &orderedUnits, const std::vector<unsigned int> &symmetry_classes)
{
  std::vector<int> v;
  for (unsigned int i = 0; i < orderedUnits.size(); ++i) {
    const StereogenicUnit &unit = orderedUnits[i];
    if (unit.type == OBStereo::Tetrahedral) {
      if (!fragment.BitIsOn(unit.id))
        continue;
      if (classifyTetrahedralNbrSymClasses(symmetry_classes, mol->GetAtomById(unit.id)) != T1234)
        continue;

      OBStereoFacade facade(mol, false);
      if (!facade.HasTetrahedralStereo(unit.id))
        continue;
      OBTetrahedralStereo::Config config = facade.GetTetrahedralStereo(unit.id)->GetConfig(); 
      IdsToSymClasses(mol, config, symmetry_classes);
      v.push_back(findDescriptor(config));
    }
  }

  return v;
}

int findDescriptorVectorValue(OBMol *mol, const OBBitVec &fragment, 
    const std::vector<StereogenicUnit> &orderedUnits, const std::vector<unsigned int> &symmetry_classes)
{
  std::vector<int> v = findDescriptorVector(mol, fragment, orderedUnits, symmetry_classes);

  int value = 0;
  for (unsigned int i = 0; i < v.size(); ++i) {
    if (!v[i])
      continue;
    int power = v.size() - i - 1;
    value += std::pow(2, power);
  }

  return value;
}

  // Check if the atom is a tetrahedral center (i.e. there is a Tetrahedral StereogenicUnit 
  // in units with the same id)
  bool isTetrahedral(OBAtom *atom, const std::vector<StereogenicUnit> &units)
  {
    for (unsigned int i = 0; i < units.size(); ++i) {
      const StereogenicUnit &unit = units[i];
      if (unit.type != OBStereo::Tetrahedral)
        continue;
      if (unit.id == atom->GetId())
        return true;
    }
    return false;
  }

  // Count the number of unique symmetry classes
  unsigned int countClasses(const std::vector<unsigned int> &symmetry_classes)
  {
    std::vector<unsigned int> copy(symmetry_classes);
    std::sort(copy.begin(), copy.end());
    return std::unique(copy.begin(), copy.end()) - copy.begin();
  }



void BreakUnit(OBMol *mol, const StereogenicUnit &unit, std::vector<unsigned int> &symmetry_classes,
    const std::vector<StereogenicUnit> &stereoUnits, int breakUnitsClassification, bool &brokenTies, bool force = false)
{
  if (unit.type == OBStereo::Tetrahedral) {
    OBAtom *center = mol->GetAtomById(unit.id);
    // the center might already be resolved by numbering a previous breakUnit...
    if (classifyTetrahedralNbrSymClasses(symmetry_classes, center) == T1234) {
      cout << "XXXXXXXXXX" << endl;
      return;
    }

    // Get the config struct for the center
    OBStereoFacade stereoFacade(mol, false);
    if (!stereoFacade.HasTetrahedralStereo(unit.id))
      return;
    cout << "..." << endl;
    OBTetrahedralStereo::Config idConfig = stereoFacade.GetTetrahedralStereo(unit.id)->GetConfig();

    switch (breakUnitsClassification) {
      case T1123:
        {
          cout << "T1123" << endl;
          // Find the duplicated symmetry class and equivalent neighbor atoms
          unsigned int duplicatedSymClass = findDuplicatedSymmetryClass(center, symmetry_classes);
          cout << "duplicatedSymClass = " << duplicatedSymClass << endl;
          std::vector<OBAtom*> equivalentNbrs;
          FOR_NBORS_OF_ATOM (nbr, center)
            if (duplicatedSymClass == symmetry_classes[nbr->GetIndex()])
              equivalentNbrs.push_back(&*nbr);

          cout << "RESOLVING: " << center->GetId() << endl;

          cout << "equivalentNbrs.size()  = " << equivalentNbrs.size() << endl;

          assert( equivalentNbrs.size() == 2 );

          // Get the fragments for the equivalent neighbor atoms
          OBBitVec fragment1 = getFragment(equivalentNbrs[0], center);
          OBBitVec fragment2 = getFragment(equivalentNbrs[1], center);

          // Check if all stereocenters in the ligands are resolved
          bool allStereoResolved = true;
          for (unsigned int k = 0; k < mol->NumAtoms(); ++k) {
            OBAtom *fragAtom = mol->GetAtom(k+1);
            if (fragment1.BitIsOn(fragAtom->GetId()) || fragment2.BitIsOn(fragAtom->GetId())) {
              if (!isTetrahedral(fragAtom, stereoUnits))
                continue;
              if (classifyTetrahedralNbrSymClasses(symmetry_classes, fragAtom) != T1234) {
                allStereoResolved = false;
                cout << "not resolved: " << fragAtom->GetId() << endl; 
              }
            }
          }

          if (allStereoResolved) {
            cout << "All stereo resolved for fragments" << endl;
            // If all stereocenters in both ligands are resolved, we need to take 
            // the descriptor values for the ligands into account.
            int dv1, dv2; // descriptor vector values

            if (fragment1 == fragment2) {
              // special case: e.g. 1,3-hydroxy-5-methyl-cyclohexane
              dv1 = findDescriptorVectorValue(mol, fragment1, 
                  orderSetByFragment(equivalentNbrs[0], center, stereoUnits), symmetry_classes);
              dv2 = findDescriptorVectorValue(mol, fragment2, 
                  orderSetByFragment(equivalentNbrs[1], center, stereoUnits), symmetry_classes);
            }  else {
              std::vector<StereogenicUnit> ordered(orderSetBySymmetryClasses(mol, stereoUnits, symmetry_classes));
              dv1 = findDescriptorVectorValue(mol, fragment1, ordered, symmetry_classes);
              dv2 = findDescriptorVectorValue(mol, fragment2, ordered, symmetry_classes);
            }
            cout << "dv1 = " << dv1 << endl;
            cout << "dv2 = " << dv2 << endl;

            // The highest descriptor vector value goes first
            brokenTies = true;
            if (dv1 > dv2) {
              symmetry_classes[equivalentNbrs[0]->GetIndex()] -= 1;
            } else
              if (dv2 > dv1) {
                symmetry_classes[equivalentNbrs[1]->GetIndex()] -= 1;                
              } 

            // If the ligands have the same descriptor vector value, the 
            // center is not a real stereocenter and we mark it unspecified.
            if (dv1 == dv2) {
              if (force) {
                symmetry_classes[equivalentNbrs[0]->GetIndex()] -= 1;
              } else {
                idConfig.specified = false;
                stereoFacade.GetTetrahedralStereo(unit.id)->SetConfig(idConfig);
              }
            }
          } else {
            cout << "All stereo NOT resolved for fragments" << endl;
            // If all stereocenters in both ligand are not resolved yet,
            // we can still choose the descriptor for the center. We always
            // choose the 1 descriptor. This choosing happens by numbering 
            // the equivalent atoms.

            brokenTies = true;

            OBTetrahedralStereo::Config symConfig;
            // try numbering 1
            symmetry_classes[equivalentNbrs[0]->GetIndex()] -= 1;
            symConfig = idConfig;
            IdsToSymClasses(mol, symConfig, symmetry_classes);

            if (!findDescriptor(symConfig)) {
              symmetry_classes[equivalentNbrs[0]->GetIndex()] += 1; // restore
              // try numbering 2
              symmetry_classes[equivalentNbrs[1]->GetIndex()] -= 1;
              symConfig = idConfig;
              IdsToSymClasses(mol, symConfig, symmetry_classes);
            }

            assert(findDescriptor(symConfig));
          }
        }
        break;
      case T1122:
        cout << "T1122" << endl;
        // old BreakChiralTies will take care of this...
        break;
      case T1112:
        {
          cout << "T1112" << endl;
          // Find the duplicated symmetry class and equivalent neighbor atoms
          unsigned int duplicatedSymClass = findDuplicatedSymmetryClass(center, symmetry_classes);
          cout << "duplicatedSymClass = " << duplicatedSymClass << endl;
          std::vector<OBAtom*> equivalentNbrs;
          FOR_NBORS_OF_ATOM (nbr, center)
            if (duplicatedSymClass == symmetry_classes[nbr->GetIndex()])
              equivalentNbrs.push_back(&*nbr);

          cout << "RESOLVING: " << center->GetId() << endl;

          cout << "equivalentNbrs.size()  = " << equivalentNbrs.size() << endl;

          assert( equivalentNbrs.size() == 3 );

          // Get the fragments for the equivalent neighbor atoms
          OBBitVec fragment1 = getFragment(equivalentNbrs[0], center);
          OBBitVec fragment2 = getFragment(equivalentNbrs[1], center);
          OBBitVec fragment3 = getFragment(equivalentNbrs[2], center);

          // Check if all stereocenters in the ligands are resolved
          bool allStereoResolved = true;
          for (unsigned int k = 0; k < mol->NumAtoms(); ++k) {
            OBAtom *fragAtom = mol->GetAtom(k+1);
            if (fragment1.BitIsOn(fragAtom->GetId()) || fragment2.BitIsOn(fragAtom->GetId()) ||
                fragment3.BitIsOn(fragAtom->GetId())) {
              if (!isTetrahedral(fragAtom, stereoUnits))
                continue;
              if (classifyTetrahedralNbrSymClasses(symmetry_classes, fragAtom) != T1234) {
                allStereoResolved = false;
                cout << "not resolved: " << fragAtom->GetId() << endl; 
              }
            }
          }

          if (allStereoResolved) {
            cout << "All stereo resolved for fragments" << endl;
            // If all stereocenters in both ligands are resolved, we need to take 
            // the descriptor values for the ligands into account.
            int dv1, dv2, dv3; // descriptor vector values

            if ((fragment1 == fragment2) && (fragment1 == fragment3)) {
              dv1 = findDescriptorVectorValue(mol, fragment1, 
                  orderSetByFragment(equivalentNbrs[0], center, stereoUnits), symmetry_classes);
              dv2 = findDescriptorVectorValue(mol, fragment2, 
                  orderSetByFragment(equivalentNbrs[1], center, stereoUnits), symmetry_classes);
              dv3 = findDescriptorVectorValue(mol, fragment3, 
                  orderSetByFragment(equivalentNbrs[2], center, stereoUnits), symmetry_classes);
            }  else {
              std::vector<StereogenicUnit> ordered(orderSetBySymmetryClasses(mol, stereoUnits, symmetry_classes));
              dv1 = findDescriptorVectorValue(mol, fragment1, ordered, symmetry_classes);
              dv2 = findDescriptorVectorValue(mol, fragment2, ordered, symmetry_classes);
              dv3 = findDescriptorVectorValue(mol, fragment3, ordered, symmetry_classes);
            }
            cout << "dv1 = " << dv1 << endl;
            cout << "dv2 = " << dv2 << endl;
            cout << "dv3 = " << dv3 << endl;

            // The highest descriptor vector value goes first
            brokenTies = true;
            if (dv1 == dv2 || dv1 == dv3 || dv2 == dv3) {
              // If the ligands have the same descriptor vector value, the 
              // center is not a real stereocenter and we mark it unspecified.
              idConfig.specified = false;
              stereoFacade.GetTetrahedralStereo(unit.id)->SetConfig(idConfig);
            } else {
              // decrement lowest ranked
              if (dv1 > dv2 && dv1 > dv3) {
                symmetry_classes[equivalentNbrs[0]->GetIndex()] -= 1;
              } else if (dv2 > dv1 && dv2 > dv3) {
                symmetry_classes[equivalentNbrs[1]->GetIndex()] -= 1;                
              } else if (dv3 > dv1 && dv3 > dv2) {
                symmetry_classes[equivalentNbrs[2]->GetIndex()] -= 1;                
              }
              // increment highest ranked
              if (dv1 < dv2 && dv1 < dv3) {
                symmetry_classes[equivalentNbrs[0]->GetIndex()] += 1;
              } else if (dv2 < dv1 && dv2 < dv3) {
                symmetry_classes[equivalentNbrs[1]->GetIndex()] += 1;                
              } else if (dv3 < dv1 && dv3 < dv2) {
                symmetry_classes[equivalentNbrs[2]->GetIndex()] += 1;                
              }
            }

          } else {
            assert(0); // should not happen AFAIK
          }
        }
        break;
      default:
        // old BreakChiralTies will take care of this...
        break;
    }

  }


}



void OBGraphSym::NewBreakChiralTies(std::vector<std::pair<OBAtom*, unsigned int> > &atom_sym_classes)
{
  cout << "ENTER NewBreakChiralTies..." << endl;
  // Convert the atom/class pairs to an array indexed by atom index (from 0).
  // This is just for convenience in the next step.  Note that there
  // will be "holes" in this vector since it's a molecule fragment.
  vector<unsigned int> symmetry_classes(_pmol->NumAtoms(), std::numeric_limits<unsigned long>::max());
  vector<pair<OBAtom*,unsigned int> >::iterator api;
  for (api = atom_sym_classes.begin(); api < atom_sym_classes.end(); api++)
    symmetry_classes[api->first->GetIndex()] = api->second;

  unsigned int beforeCount = countClasses(symmetry_classes);
  //DepictSymmetryClasses(_pmol, symmetry_classes, "img/NewBreakChiralTies_ENTER");

  // FIXME 
  // Skip breaking chiral ties if there are no stereogenic units with defined configuration.
  // This is only for performance and not part of the algorithm.
  //if (_pmol->HasChiralityPerceived()) {
    bool hasAtLeastOneDefined = false;
    OBStereoFacade sf(_pmol, false);
    FOR_ATOMS_OF_MOL (atom, _pmol) {
      if (sf.HasTetrahedralStereo(atom->GetId())) {
        if (sf.GetTetrahedralStereo(atom->GetId())->GetConfig().specified) {
          hasAtLeastOneDefined = true;
          break;
        }
      }
    }
    if (!hasAtLeastOneDefined)
      return;
  //}

  // Compute the automorphisms and derived stereogenic units only once.
  // Computing this every time would result in different interdependent 
  // units since the symmetry classes would be different after breaking
  // a tie.
  if (!m_G.Size()) {
    // Find all automorphisms
    m_G = FindAutomorphisms(_pmol, symmetry_classes);
    // Find all types of stereogenic units (i.e. tetrahedral, cis/trans, ...)
    m_stereoUnits = FindStereogenicUnits(_pmol, symmetry_classes, m_G);
    // Find the interdependent sets of stereogenic units
    m_interdependentSets = FindInterdependentStereogenicUnits(_pmol,
      m_stereoUnits, symmetry_classes, m_G);
  }

  // Determine stereochemistry from coordinates if needed
  if (!_pmol->HasChiralityPerceived()) {
    switch (_pmol->GetDimension()) {
      case 2:
        _pmol->DeleteData(OBGenericDataType::StereoData);
        TetrahedralFrom2D(_pmol, m_stereoUnits);
        break;
      case 3:
        _pmol->DeleteData(OBGenericDataType::StereoData);
        TetrahedralFrom3D(_pmol, m_stereoUnits);
        break;
      default:
        TetrahedralFrom0D(_pmol, m_stereoUnits);
        break;
    }
  }
 
  // If there are no unresolved stereocenters, we're done
  // Unresolved: classification != T1234
  // This is only for performance and not part of the algorithm.
  bool foundTie = false;
  for (unsigned int i = 0; i < m_stereoUnits.size(); ++i) {
    const StereogenicUnit &unit = m_stereoUnits[i];
    if (unit.type == OBStereo::Tetrahedral) {
      OBAtom *atom = _pmol->GetAtomById(unit.id);
      if (classifyTetrahedralNbrSymClasses(symmetry_classes, atom) != T1234) {
        foundTie = true;
      }
    }
  }
  if (!foundTie)
    return;

  //////////////////////////////////////////////////////////////////////////////
  //
  // Start tie breaking for stereocenters that coincide with the symmetry axis
  //
  //////////////////////////////////////////////////////////////////////////////

  // Find duplicated sets of interdependent stereogenic units.
  // duplicated: identical sets (containing same elements) after ids are 
  // replaced by symmetry classes
  std::vector<unsigned int> doneIndexes;
  std::vector< std::vector< std::vector<StereogenicUnit> > > duplicatedSets;
  for (unsigned int i = 0; i < m_interdependentSets.size(); ++i) {
    // skip already handled units
    if (std::find(doneIndexes.begin(), doneIndexes.end(), i) != doneIndexes.end())
      continue;
    // add the i-th set
    const std::vector<StereogenicUnit> &iSet = m_interdependentSets[i];
    std::vector< std::vector<StereogenicUnit> > sets;
    sets.push_back(iSet);
    doneIndexes.push_back(i);
    for (unsigned int j = 0; j < m_interdependentSets.size(); ++j) {
      // skip already handled units
      if (std::find(doneIndexes.begin(), doneIndexes.end(), j) != doneIndexes.end())
        continue;
      const std::vector<StereogenicUnit> &jSet = m_interdependentSets[j];
      if (iSet.size() != jSet.size())
        continue;
      // add the j-th set if it contains the same elements 
      if (setsContainSameSymmetryClasses(_pmol, iSet, jSet, symmetry_classes)) {
        sets.push_back(jSet);
        doneIndexes.push_back(j);
      }
    }

    duplicatedSets.push_back(sets);
  }

  // Sort duplicated sets by their size
  // This is a key step in the algorithm since it ensures the most duplicated
  // units are resolved first. This way, any stereocenter depending on stereocenters
  // in it's duplicated ligands will be resolved after these ligand stereocenters.
  std::sort(duplicatedSets.begin(), duplicatedSets.end(), CompareSetSizes); 

  // DEBUG print the sorted duplicated sets
  for (unsigned int i = 0; i  < duplicatedSets.size(); ++i) {
    const std::vector< std::vector<StereogenicUnit> > &sets = duplicatedSets[i];
    cout << "Duplicated set:" << endl;
    for (unsigned int j = 0; j < sets.size(); ++j) {
      cout << "    Set " << j+1 << endl << "        ";
      for (unsigned int k = 0; k < sets[j].size(); ++k) {
        cout << sets[j][k].id << " ";
      }
      cout << endl;
    }
  }
  // DEBUG

  // Multiply all symmetry classes by 2. This enables us to break ties by
  // substracting 1 from a symmetry class and iterating to propagate the 
  // changes.
  for (unsigned int j = 0; j < symmetry_classes.size(); ++j)
    symmetry_classes[j] *= 3;

  // Sort the duplicated sets
  for (unsigned int i = 0; i  < duplicatedSets.size(); ++i) { // foreach duplicated set of sets
    std::vector< std::vector<StereogenicUnit> > &sets = duplicatedSets[i];
    for (unsigned int j = 0; j < sets.size(); ++j) { // foreach interdependent set
      std::vector<StereogenicUnit> &set = sets[j];
      set = orderSetBySymmetryClasses(_pmol, set, symmetry_classes);
    }
  }
 
  // Start the actual tie breaking
  bool brokenTies = false;
  for (unsigned int i = 0; i  < duplicatedSets.size(); ++i) { // foreach duplicated set of sets
    std::vector< std::vector<StereogenicUnit> > &sets = duplicatedSets[i];
  
    // Do tie breaking for cases where there is no duplicate and the single set consists
    // of more than 2 stereogenic units with the same symmetry class. (e.g. inositol)
    if (sets.size() == 1) {
      std::vector<StereogenicUnit> &set = sets[0];
      if (set.size() > 2) {
        bool allSameSymClass = true;
        

        if (set[0].type == OBStereo::Tetrahedral) {
          unsigned int symClass = symmetry_classes[_pmol->GetAtomById(set[0].id)->GetIndex()];
          for (unsigned int k = 1; k < set.size(); ++k) { // foreach stereogenic unit
            const StereogenicUnit &unit = set[k];
            if (unit.type == OBStereo::Tetrahedral) {
              OBAtom *center = _pmol->GetAtomById(unit.id);
              if (symmetry_classes[center->GetIndex()] != symClass)
                allSameSymClass = false;
            }
          }

          if (allSameSymClass) {
            int canonValue = -1;
            std::vector<unsigned int> symClassesCanon;

            OBBitVec fragment;
            FOR_RINGS_OF_MOL (ring, _pmol) {
              bool allUnitsInRing = true;
              for (unsigned int l = 0; l < set.size(); ++l) {
                const StereogenicUnit &unit = set[l];
                if (unit.type == OBStereo::Tetrahedral) {
                  OBAtom *center = _pmol->GetAtomById(unit.id);
                  if (std::find(ring->_path.begin(), ring->_path.end(), center->GetIndex()+1) == ring->_path.end())
                    allUnitsInRing = false;
                }
              }
              if (allUnitsInRing) {
                for (unsigned int l = 0; l < ring->_path.size(); ++l)
                  fragment.SetBitOn(_pmol->GetAtom(ring->_path[l])->GetId());
                break;
              }
            }
            /*
            for (unsigned int l = 0; l < set.size(); ++l) {
              const StereogenicUnit &unit = set[l];
              if (unit.type == OBStereo::Tetrahedral) {
                fragment.SetBitOn(unit.id);
              }
            }
            */

            set = orderSetByRing(_pmol, fragment, set);

            // Check starting from each unit in forward direction
            for (unsigned int j = 0; j < set.size(); ++j) {
              std::vector<unsigned int> symClassesCopy = symmetry_classes;
              for (unsigned int l = 0; l < symClassesCopy.size(); ++l)
                symClassesCopy[l] += set.size() + 1;

              for (unsigned int l = 0; l < set.size(); ++l) {
                const StereogenicUnit &unit = set[l];
                if (unit.type == OBStereo::Tetrahedral) {
                  OBAtom *center = _pmol->GetAtomById(unit.id);
                  symClassesCopy[center->GetIndex()] = l;
                }
              }
              Iterate(symClassesCopy);

              int value = findDescriptorVectorValue(_pmol, fragment, set, symClassesCopy);
              if (value > canonValue) {
                canonValue = value;
                symClassesCanon = symClassesCopy;
              }
              cout << "DEBUG " << value << endl;

              StereogenicUnit u0 = set[0];
              for (unsigned int l = 1; l < set.size(); ++l) {
                set[l-1] = set[l];
              }
              set[set.size()-1] = u0;
            }

            // Reverse the order of the set
            std::vector<StereogenicUnit> set2;
            for (unsigned int j = 0; j < set.size(); ++j)
              set2.push_back(set[set.size() - j - 1]);
            set = set2;

            // Check starting from each unit in backward direction
            // (Code is identical copy of code above)
            for (unsigned int j = 0; j < set.size(); ++j) {
              std::vector<unsigned int> symClassesCopy = symmetry_classes;
              for (unsigned int l = 0; l < symClassesCopy.size(); ++l)
                symClassesCopy[l] += set.size() + 1;

              for (unsigned int l = 0; l < set.size(); ++l) {
                const StereogenicUnit &unit = set[l];
                if (unit.type == OBStereo::Tetrahedral) {
                  OBAtom *center = _pmol->GetAtomById(unit.id);
                  symClassesCopy[center->GetIndex()] = l;
                }
              }
              Iterate(symClassesCopy);

              int value = findDescriptorVectorValue(_pmol, fragment, set, symClassesCopy);
              if (value > canonValue) {
                canonValue = value;
                symClassesCanon = symClassesCopy;
              }
              cout << "DEBUG " << value << endl;

              StereogenicUnit u0 = set[0];
              for (unsigned int l = 1; l < set.size(); ++l) {
                set[l-1] = set[l];
              }
              set[set.size()-1] = u0;
            }



            DepictSymmetryClasses(_pmol, symClassesCanon, "img/NewBreakChiralTies");
            symmetry_classes = symClassesCanon;

            for (unsigned int i = 0; i < atom_sym_classes.size(); ++i)
              atom_sym_classes[i].second =  symmetry_classes[atom_sym_classes[i].first->GetIndex()];

            return;
          }
        }
      }
    }

    // Select the unresolved units with lowest symmetry class
    // unresolved: classification != T1234
    int breakUnitsClassification;
    std::vector<StereogenicUnit> breakUnits;
    cout << "BREAK UNITS: ";
    for (unsigned int j = 0; j < sets.size(); ++j) { // foreach interdependent set
      const std::vector<StereogenicUnit> &set = sets[j];
      
      for (unsigned int k = 0; k < set.size(); ++k) { // foreach stereogenic unit
        const StereogenicUnit &unit = set[k];
        if (unit.type == OBStereo::Tetrahedral) {
          OBAtom *center = _pmol->GetAtomById(unit.id);
          int classification = classifyTetrahedralNbrSymClasses(symmetry_classes, center);
          if (classification != T1234) {
            breakUnitsClassification = classification; // will be the same for all units...
            breakUnits.push_back(unit);
            cout << unit.id << " ";
            break;          
          }
        }
      }
    }
    cout << endl;

    if (breakUnits.empty())
      continue;

    // Special case: handle 4 ring 
    bool breakOnlyOne = false;
    if (breakUnits.size() == 2) {
      if (breakUnits[0].type == OBStereo::Tetrahedral && breakUnits[1].type == OBStereo::Tetrahedral) {
        OBAtom *breakAtom1 = _pmol->GetAtomById(breakUnits[0].id);
        OBAtom *breakAtom2 = _pmol->GetAtomById(breakUnits[1].id);

        if (breakAtom1->IsConnected(breakAtom2))
          breakOnlyOne = true;
      }
    }

    for (unsigned int j = 0; j < breakUnits.size(); ++j) {
      const StereogenicUnit &unit = breakUnits[j];
      cout << "RESOLV BREAK UNIT  " << unit.id << endl;
      BreakUnit(_pmol, unit, symmetry_classes, m_stereoUnits, breakUnitsClassification, brokenTies);
      if (breakOnlyOne)
        break;
    }

    // Break only 1 tie
    if (brokenTies)
      break;
  }

  Iterate(symmetry_classes);
  for (unsigned int i = 0; i < atom_sym_classes.size(); ++i)
    atom_sym_classes[i].second =  symmetry_classes[atom_sym_classes[i].first->GetIndex()];
  
    
  //DepictSymmetryClasses(_pmol, symmetry_classes, "img/NewBreakChiralTies_EXIT");
  cout << "EXIT NewBreakChiralTies..." << endl;
  
  unsigned int afterCount = countClasses(symmetry_classes);

  if (afterCount > beforeCount)
    NewBreakChiralTies(atom_sym_classes);
}


void OBGraphSym::BreakChiralTies(vector<pair<OBAtom*, unsigned int> > &atom_sym_classes)
{
  vector<pair<OBAtom*,unsigned int> > vp1, vp2;

  // for keeping track of atoms we've already considered
  OBBitVec used_atoms;
  used_atoms.Clear();
  used_atoms.Resize(_pmol->NumAtoms());

  // Convert the atom/class pairs to an array indexed by atom idx.
  // This is just for convenience in the next step.  Note that there
  // will be "holes" in this vector since it's a molecule fragment.
  vector<unsigned int> index2sym_class(_pmol->NumAtoms());
  vector<pair<OBAtom*,unsigned int> >::iterator api;
  for (api = atom_sym_classes.begin(); api < atom_sym_classes.end(); api++)
    index2sym_class[api->first->GetIndex()] = api->second;
  
  // find all types of stereogenic units (i.e. tetrahedral, cis/trans, ...)
  vector<StereogenicUnit> stereoUnits = FindStereogenicUnits(_pmol, index2sym_class);

  // 
  // Tetrahedral atoms
  //
  for (vector<StereogenicUnit>::iterator unit1 = stereoUnits.begin(); unit1 != stereoUnits.end(); ++unit1) {
    if ((*unit1).type != OBStereo::Tetrahedral)
      continue;
    unsigned long id1 = (*unit1).id;
    OBAtom *atom1 = _pmol->GetAtomById(id1);
    int index1 = atom1->GetIndex();
    // We only want: unused and part of this fragment.
    if (!(*_frag_atoms)[index1+1])
      continue;
    if (used_atoms[index1])
      continue;
    used_atoms.SetBitOn(index1);
    //if (GetValence(atom) < 4) // Valence relative to this fragment
    //  continue;

    // Get its symmetry class
    int symclass = index2sym_class[index1];
 
    // Start the vector of "same class" atoms by adding this atom
    vector<OBAtom *> same_class;
    same_class.push_back(atom1);

    // Inner Loop: Find all other atoms with the same symmetry class
    for (vector<StereogenicUnit>::iterator unit2 = unit1+1; unit2 != stereoUnits.end(); ++unit2) {
      if ((*unit2).type != OBStereo::Tetrahedral)
        continue;
      unsigned long id2 = (*unit2).id;
      OBAtom *atom2 = _pmol->GetAtomById(id2);
      int index2 = atom2->GetIndex();
      if (used_atoms[index2])
        continue;
      if (index2sym_class[index2] == symclass) {
        same_class.push_back(atom2);
        used_atoms.SetBitOn(index2);
      }
    }
    
    // Unless at least two atoms in the class, there are no ties to break
    if (same_class.size() < 2)
      continue;

    /*cout << "BreakChiralTies: same_class = ";
    vector<OBAtom*>::iterator ia;
    for (ia = same_class.begin(); ia != same_class.end(); ia++)
      cout << (*ia)->GetIndex() << " ";
    cout << "\n";*/

    // find tetrahedral atoms using current symmetry classes
    if (!_pmol->HasChiralityPerceived()) {
      //cout << "!_pmol->HasChiralityPerceived()" << endl;
      switch (_pmol->GetDimension()) {
        case 2:
          _pmol->DeleteData(OBGenericDataType::StereoData);
          TetrahedralFrom2D(_pmol, stereoUnits);
          break;
        case 3:
          _pmol->DeleteData(OBGenericDataType::StereoData);
          TetrahedralFrom3D(_pmol, stereoUnits);
          break;
        default:
          TetrahedralFrom0D(_pmol, stereoUnits);
          break;
      }
    } 

    // false = don't perceive stereochemistry, we already did it explicitly above
    OBStereoFacade stereoFacade(_pmol, false); 
    // get the reference config
    OBTetrahedralStereo::Config refConfig;
    if (stereoFacade.HasTetrahedralStereo(id1)) {
      refConfig = stereoFacade.GetTetrahedralStereo(id1)->GetConfig();
    } else {
      refConfig.specified = false;
    }
    // convert id --> symmetry classes to compare
    IdsToSymClasses(_pmol, refConfig, index2sym_class);

    // Devide the centers in 3 groups based on comparing the Configs
    // symclass1: all centers who's Config struct matches refConfig 
    // symclass2: all centers who's Config struct does not match refConfig 
    // unspecified: all centers with unspecified stereochemistry
    vector<OBAtom*> symclass1, symclass2, unspecified;

    vector<OBAtom*>::iterator iatom;
    for (iatom = same_class.begin(); iatom != same_class.end(); iatom++) {
      if (!stereoFacade.HasTetrahedralStereo((*iatom)->GetId())) {
        // this happens for 0D when we don't want to delete data (see above) 
        // but have found new previously unidentified chiral centers
        unspecified.push_back(*iatom);
      } else {
        OBTetrahedralStereo::Config otherConfig = stereoFacade.GetTetrahedralStereo((*iatom)->GetId())->GetConfig();
        // unspecified is group 3
        if (!otherConfig.specified) {
          unspecified.push_back(*iatom);
          continue;
        }

        // if refConfig is still unspecified, assign otherConfig to it and 
        // otherConfig will go in symclass1
        if (!refConfig.specified)
          refConfig = otherConfig;

        IdsToSymClasses(_pmol, otherConfig, index2sym_class);
        // compare
        if (refConfig == otherConfig)
          symclass1.push_back(*iatom); // same, add it to symclass1
        else
          symclass2.push_back(*iatom); // different, add to symclass2
      }
    }

    // If there's nothing in symclass2, then we don't have to split 
    // the symmetry class.
    if (symclass1.empty())
      if (symclass2.empty() || unspecified.empty())
        continue;
    if (symclass2.empty() && unspecified.empty())
      continue;

    //
    // Make a copy of refConfig and sort it
    //
    OBTetrahedralStereo::Config orderedConfig;    // constructor initializes view & winding --> canonical
    // copy specified flag and center Ref
    orderedConfig.center = refConfig.center;
    orderedConfig.specified = refConfig.specified;

    std::vector<unsigned long> refs = refConfig.refs;
    refs.insert(refs.begin(), refConfig.from); // doesn't matter if view is from or towards, will be sorted anyway
    std::sort(refs.begin(), refs.end());
    for (unsigned int i = 0; i < refs.size(); ++i)
      if (i)
        orderedConfig.refs.push_back(refs.at(i));
      else
        orderedConfig.from = refs.at(0);
    // store the match/mismatch result
    bool refMatchesOrdered = (refConfig == orderedConfig) ? true : false;

    //cout << "refConfig = " << refConfig << endl;
    //cout << "orderedConfig = " << orderedConfig << endl;

    // Time to break the class in 3. 
    //                         
    //                         1-2-3  symclass1
    //                        /             \ 
    //   unspecified  3-2-1--a               symclass1 != symclass2
    //                        \             /
    //                         1-2-3  symclass2
    //
    //                           |
    //                           | triple all symmetry classes
    //                           V
    //
    //                         3-6-9  symclass1  -1 <----- mismatch -----+
    //                        /             \                            |
    //   unspecified  9-6-3--a               check which one matches orderedConfig
    //      (=)               \             /                            |
    //                         3-6-9  symclass2  +1 <------ match -------+
    //
    //                           |
    //                           | perform illustrated operations
    //                           V
    //
    //                         2-5-8  symclass1  -1 <----- mismatch -----+
    //                        /             \                            |
    //   unspecified  9-6-3--a               check which one matches orderedConfig
    //      (=)               \             /                            |
    //                         4-7-10 symclass2  +1 <------ match -------+
    //
    //
    // Triple all symmetry classes.
    for (int i = 0; i < atom_sym_classes.size(); i++) {
      atom_sym_classes[i].second *= 3;
      
      for (int j = 0; j < symclass1.size(); j++) {
        if (symclass1[j] == atom_sym_classes[i].first) {
          if (refMatchesOrdered)
            atom_sym_classes[i].second += 1; // symclass1 == orderedConfig
          else
            atom_sym_classes[i].second -= 1; // symclass1 != orderedConfig
        }
      }
      for (int j = 0; j < symclass2.size(); j++) {
        if (symclass2[j] == atom_sym_classes[i].first) {
          if (refMatchesOrdered)
            atom_sym_classes[i].second -= 1; // symclass1 == orderedConfig --> symclass2 != orderedConfig
          else
            atom_sym_classes[i].second += 1; // symclass1 != orderedConfig --> symclass2 == orderedConfig
        }
      }
    }
    
    // Now propagate the change across the whole molecule with the
    // extended sum-of-invariants.
    //ExtendInvariants(atom_sym_classes);

    ///cout << "AFTER ExtendInvariants" << endl;
    //for (int i = 0; i < atom_sym_classes.size(); i++) {
    //  cout << atom_sym_classes[i].first->GetIndex() << ": " << atom_sym_classes[i].second << endl;
    //}


  }
  // 
  // Cis/Trans bonds
  //
  for (vector<StereogenicUnit>::iterator unit1 = stereoUnits.begin(); unit1 != stereoUnits.end(); ++unit1) {
    if ((*unit1).type != OBStereo::CisTrans)
      continue;

    unsigned long id1 = (*unit1).id;
    OBBond *bond1 = _pmol->GetBondById(id1);
    unsigned int begin1 = bond1->GetBeginAtom()->GetIndex();
    unsigned int end1 = bond1->GetEndAtom()->GetIndex();
    
    // We only want: unused and part of this fragment.
    if (!(*_frag_atoms)[begin1+1] || !(*_frag_atoms)[end1+1])
      continue;
    if (used_atoms[begin1] || used_atoms[end1])
      continue;
    used_atoms.SetBitOn(begin1);
    used_atoms.SetBitOn(end1);

    // Get its symmetry class
    int beginSymClass = index2sym_class[begin1];
    int endSymClass = index2sym_class[end1];
 
    // Start the vector of "same class" bonds by adding this bond
    vector<OBBond*> same_class;
    same_class.push_back(bond1);

    // Inner Loop: Find all other bonds with the same symmetry classes
    for (vector<StereogenicUnit>::iterator unit2 = unit1 + 1; unit2 != stereoUnits.end(); ++unit2) {
      if ((*unit2).type != OBStereo::CisTrans)
        continue;

      unsigned long id2 = (*unit2).id;
      OBBond *bond2 = _pmol->GetBondById(id2);
      unsigned int begin2 = bond2->GetBeginAtom()->GetIndex();
      unsigned int end2 = bond2->GetEndAtom()->GetIndex();
      if (used_atoms[begin2] || used_atoms[end2])
        continue;
      if (((index2sym_class[begin2] == beginSymClass) && (index2sym_class[end2] == endSymClass)) ||
          ((index2sym_class[begin2] == endSymClass) && (index2sym_class[end2] == beginSymClass))) {
        same_class.push_back(bond2);
        used_atoms.SetBitOn(begin2);
        used_atoms.SetBitOn(end2);
      }
    }
    
    // Unless at least two atoms in the class, there are no ties to break
    if (same_class.size() < 2)
      continue;

    /*cout << "BreakChiralTies: same_class = ";
    vector<OBBond*>::iterator ib;
    for (ib = same_class.begin(); ib != same_class.end(); ib++)
      cout << (*ib)->GetIdx() << " ";
    cout << "\n";*/

    // find cis/trans bonds using current symmetry classes
    if (!_pmol->HasChiralityPerceived()) {
      //cout << "!_pmol->HasChiralityPerceived()" << endl;
      switch (_pmol->GetDimension()) {
        case 2:
          _pmol->DeleteData(OBGenericDataType::StereoData);
          CisTransFrom2D(_pmol, stereoUnits);
          break;
        case 3:
          _pmol->DeleteData(OBGenericDataType::StereoData);
          CisTransFrom3D(_pmol, stereoUnits);
          break;
        default:
          CisTransFrom0D(_pmol, stereoUnits);
          break;
      }
    }

    // false = don't perceive stereochemistry, we already did it explicitly above
    OBStereoFacade stereoFacade(_pmol, false); 
    // get the reference config
    OBCisTransStereo::Config refConfig;
    if (stereoFacade.HasCisTransStereo(id1)) {
      refConfig = stereoFacade.GetCisTransStereo(id1)->GetConfig();
    } else {
      refConfig.specified = false;
    }
    // convert id --> symmetry classes to compare
    IdsToSymClasses(_pmol, refConfig, index2sym_class);

    // Devide the cis/trans bonds in 3 groups based on comparing the Configs
    // symclass1: all bonds who's Config struct matches refConfig 
    // symclass2: all bonds who's Config struct does not match refConfig 
    // unspecified: all bonds with unspecified stereochemistry
    vector<OBBond*> symclass1, symclass2, unspecified;

    vector<OBBond*>::iterator ibond;
    for (ibond = same_class.begin(); ibond != same_class.end(); ibond++) {
      if (!stereoFacade.HasCisTransStereo((*ibond)->GetId())) {
        // this happens for 0D when we don't want to delete data (see above) 
        // but have found new previously unidentified chiral centers
        unspecified.push_back(*ibond);
      } else {
        OBCisTransStereo::Config otherConfig = stereoFacade.GetCisTransStereo((*ibond)->GetId())->GetConfig();
        // unspecified is group 3
        if (!otherConfig.specified) {
          unspecified.push_back(*ibond);
          continue;
        }

        // if refConfig is still unspecified, assign otherConfig to it and 
        // otherConfig will go in symclass1
        if (!refConfig.specified)
          refConfig = otherConfig;

        IdsToSymClasses(_pmol, otherConfig, index2sym_class);
        // compare
        if (refConfig == otherConfig)
          symclass1.push_back(*ibond); // same, add it to symclass1
        else
          symclass2.push_back(*ibond); // different, add to symclass2
      }
    }

    // if there is only 1 group found, leave the symmetry classes alone
    if (symclass1.empty())
      if (symclass2.empty() || unspecified.empty())
        continue;
    if (symclass2.empty() && unspecified.empty())
      continue;

    // Make a copy of refConfig and sort it
    OBCisTransStereo::Config orderedConfig = refConfig;
    std::sort(orderedConfig.refs.begin(), orderedConfig.refs.end());
    // store the match/mismatch result
    bool refMatchesOrdered = (refConfig == orderedConfig) ? true : false;

    //cout << "refConfig = " << refConfig << endl;
    //cout << "orderedConfig = " << orderedConfig << endl;

    // Time to break the class in 3. (see above for details)
    for (int i = 0; i < atom_sym_classes.size(); i++) {
      atom_sym_classes[i].second *= 3;
      
      for (int j = 0; j < symclass1.size(); j++) {
        if (symclass1[j]->GetBeginAtom() == atom_sym_classes[i].first) {
          if (refMatchesOrdered)
            atom_sym_classes[i].second += 1; // symclass1 == orderedConfig
          else
            atom_sym_classes[i].second -= 1; // symclass1 != orderedConfig
        }
      }
      for (int j = 0; j < symclass2.size(); j++) {
        if (symclass2[j]->GetBeginAtom() == atom_sym_classes[i].first) {
          if (refMatchesOrdered)
            atom_sym_classes[i].second -= 1; // symclass1 == orderedConfig --> symclass2 != orderedConfig
          else
            atom_sym_classes[i].second += 1; // symclass1 != orderedConfig --> symclass2 == orderedConfig
        }
      }
    }

    // Now propagate the change across the whole molecule with the
    // extended sum-of-invariants.
    //ExtendInvariants(atom_sym_classes);
  }
 

}

/***************************************************************************
* FUNCTION: CreateNewClassVector
*
* DESCRIPTION:
*       NOTE: Derived from OpenBabel/mol.cpp
*
*       Creates a new vector of symmetry classes based on an existing
*       vector.  (Helper routine to GetGIDVector.)  On return, vp2 will
*       have newly-extended connectivity sums, but the numbers (the class
*       IDs) are very large.
*
*       (Comments by CJ) This appears to compute the "extended connectivity
*       sums" similar to those described by Weininger, Morgan, etc. It uses
*       vp1 as its starting point (the current connectivity sums), and puts
*       the new sums in vp2.  Note that vp1 is modified along the way.
*
*       Note that, per Weininger's warning, this assumes the initial class
*       ID's are less than 100, which is a BAD assumption, e.g. OCC...CCN
*       would have more than 100 symmetry classes if the chain is more than
*       98 carbons long.  Should change this to use Weininger's product of
*       corresponding primes.
***************************************************************************/

  void OBGraphSym::CreateNewClassVector(vector<pair<OBAtom*,unsigned int> > &vp1,
                                        vector<pair<OBAtom*,unsigned int> > &vp2)
{
  int m,id;
  OBAtom *atom, *nbr;
  vector<OBEdgeBase*>::iterator nbr_iter;
  vector<unsigned int>::iterator k;
  vector<pair<OBAtom*,unsigned int> >::iterator vp_iter;

#if DEBUG
  cout << "CreateNewClassVector: START\n";
  print_vector_pairs("    ", vp1);
#endif

  // There may be fewer atoms than in the whole molecule, so we can't
  // index the vp1 array by atom->GetIdx().  Instead, create a quick
  // mapping vector of idx-to-index for vp1.
  vector<int> idx2index(_pmol->NumAtoms() + 1, -1);  // natoms + 1
  int index = 0;
  for (vp_iter = vp1.begin(); vp_iter != vp1.end(); vp_iter++) {
    int idx = vp_iter->first->GetIdx();
    idx2index[idx] = index++;
  }

  // vp2 will hold the newly-extended symmetry classes
  vp2.resize(vp1.size());
  vp2.clear();

  // Loop over original atoms.
  // Create a new extended varient for each atom.  Get its neighbors' class ID's,
  // sort them into ascending order, and create a sum of (c0 + c1*10^2 + c2*10^4 + ...)
  // which becomes the new class ID (where c0 is the current classID).
  
  for (vp_iter = vp1.begin(); vp_iter != vp1.end(); vp_iter++) {
    atom = vp_iter->first;
    id   = vp_iter->second;
    vector<unsigned int> vtmp;
    for (nbr = atom->BeginNbrAtom(nbr_iter); nbr; nbr = atom->NextNbrAtom(nbr_iter)) {
      int idx = nbr->GetIdx();
      if (_frag_atoms->BitIsOn(idx))
        vtmp.push_back(vp1[idx2index[idx]].second);
    }

    sort(vtmp.begin(),vtmp.end(),CompareUnsigned);
    for (m = 100, k = vtmp.begin(); k != vtmp.end(); k++, m*=100) 
      id += *k * m;
    vp2.push_back(pair<OBAtom*,unsigned int> (atom, id));
  }
#if DEBUG
  cout << "CreateNewClassVector: FINISH\n";
  print_vector_pairs("    ", vp2);
#endif

}

void OBGraphSym::CreateNewClassVector(OBMol *mol, vector<pair<OBAtom*,unsigned int> > &vp1,
                                        vector<pair<OBAtom*,unsigned int> > &vp2)
{
  int m,id;
  OBAtom *atom, *nbr;
  vector<OBEdgeBase*>::iterator nbr_iter;
  vector<unsigned int>::iterator k;
  vector<pair<OBAtom*,unsigned int> >::iterator vp_iter;

#if DEBUG
  cout << "CreateNewClassVector: START\n";
  print_vector_pairs("    ", vp1);
#endif

  // There may be fewer atoms than in the whole molecule, so we can't
  // index the vp1 array by atom->GetIdx().  Instead, create a quick
  // mapping vector of idx-to-index for vp1.
  vector<int> idx2index(mol->NumAtoms() + 1, -1);  // natoms + 1
  int index = 0;
  for (vp_iter = vp1.begin(); vp_iter != vp1.end(); vp_iter++) {
    int idx = vp_iter->first->GetIdx();
    idx2index[idx] = index++;
  }

  // vp2 will hold the newly-extended symmetry classes
  vp2.resize(vp1.size());
  vp2.clear();

  // Loop over original atoms.
  // Create a new extended varient for each atom.  Get its neighbors' class ID's,
  // sort them into ascending order, and create a sum of (c0 + c1*10^2 + c2*10^4 + ...)
  // which becomes the new class ID (where c0 is the current classID).
  
  for (vp_iter = vp1.begin(); vp_iter != vp1.end(); vp_iter++) {
    atom = vp_iter->first;
    id   = vp_iter->second;
    vector<unsigned int> vtmp;
    for (nbr = atom->BeginNbrAtom(nbr_iter); nbr; nbr = atom->NextNbrAtom(nbr_iter)) {
      int idx = nbr->GetIdx();
      vtmp.push_back(vp1[idx2index[idx]].second);
    }

    sort(vtmp.begin(),vtmp.end(),CompareUnsigned);
    for (m = 100, k = vtmp.begin(); k != vtmp.end(); k++, m*=100) 
      id += *k * m;
    vp2.push_back(pair<OBAtom*,unsigned int> (atom, id));
  }
#if DEBUG
  cout << "CreateNewClassVector: FINISH\n";
  print_vector_pairs("    ", vp2);
#endif

}

/***************************************************************************
* FUNCTION: CountAndRenumberClasses
*
* DESCRIPTION:
*       NOTE: Copied from OpenBabel/mol.cpp
*
*       Counts the number of unique symmetry classes in a list.
*
*       (NOTE: CJ -- It also appears to MODIFY the list.  It sorts it in order
*       of class ID, then renumbers the ID's zero through N-1.  See the comments
*       in CreateNewClassVector() about how it returns very large numbers for the
*       class IDs it creates.  These are replaced by lower, sequential numbers here.)
***************************************************************************/

  void OBGraphSym::CountAndRenumberClasses(vector<pair<OBAtom*,unsigned int> > &vp,
                                           unsigned int &count)
  {
    count = 1;
    vector<pair<OBAtom*,unsigned int> >::iterator k;

    sort(vp.begin(), vp.end(), ComparePairSecond);
    k = vp.begin();
    if (k != vp.end()) {
      unsigned int id = k->second;
      k->second = 1;
      ++k;
      for (;k != vp.end(); ++k) {
        if (k->second != id) {
          id = k->second;
          k->second = ++count;
        } else {
          k->second = count;
        }
      }
    }
  }


/***************************************************************************
* FUNCTION: ExtendInvariants
*
* DESCRIPTION:
*       This is the core of symmetry analysis.  Starting with a set of
*       classes on each atom, it "spreads" them using a sum-of-invariants
*       of each atom's class and its neighbors' classes.  This iterates
*       until a stable solution is found (further spreading doesn't
*       change the answer).
*       
* RETURNS: The number of distinct symmetry classes found.
***************************************************************************/

  int OBGraphSym::ExtendInvariants(vector<pair<OBAtom*, unsigned int> > &symmetry_classes, bool breakChiralTies)
  {
    unsigned int nclasses1, nclasses2;
    vector<pair<OBAtom*,unsigned int> > tmp_classes;

    // How many classes are we starting with?  (The "renumber" part isn't relevant.)
    CountAndRenumberClasses(symmetry_classes, nclasses1);

    int natoms = _pmol->NumAtoms();
    int nfragatoms = _frag_atoms->CountBits();

    // LOOP: Do extended sum-of-invarients until no further changes are
    // noted.  (Note: This is inefficient, as it re-computes extended sums
    // and re-sorts the entire list each time.  You can save a lot of time by
    // only recomputing and resorting within regions where there is a tie
    // initially.  But it's a lot more code.)

    if (nclasses1 < nfragatoms) {
      for (int i = 0; i < 100;i++) {  //sanity check - shouldn't ever hit this number
        CreateNewClassVector(symmetry_classes, tmp_classes);
        CountAndRenumberClasses(tmp_classes, nclasses2);
        symmetry_classes = tmp_classes;
        if (nclasses1 == nclasses2) break;
        nclasses1 = nclasses2;
      }
    }

    if (breakChiralTies) {
    /*cout << "BEFORE BreakChiralTies, nclasses1 = " << nclasses1 << endl;
      for (int i = 0; i < symmetry_classes.size(); i++) {
      cout << symmetry_classes[i].first->GetIndex() << ": " << symmetry_classes[i].second << endl;
      }*/

      NewBreakChiralTies(symmetry_classes);
      /*
      std::vector<unsigned int> symClasses(symmetry_classes.size());
      for (unsigned int i = 0; i < symmetry_classes.size(); ++i)
        symClasses[symmetry_classes[i].first->GetIndex()] = symmetry_classes[i].second;
        */

      //DepictSymmetryClasses(_pmol, symClasses, "img/NewBreakChiralTies");
      // handle cases where the symmetry axis doesn't coincide with stereocenter
      // NEEDS TO BE MERGED
      BreakChiralTies(symmetry_classes); 
      //BreakChiralTies(symmetry_classes);
      CreateNewClassVector(symmetry_classes, tmp_classes);
      CountAndRenumberClasses(tmp_classes, nclasses2);

    /*cout << "AFTER BreakChiralTies, nclasses2 = " << nclasses2 << endl;
      for (int i = 0; i < symmetry_classes.size(); i++) {
      cout << symmetry_classes[i].first->GetIndex() << " (" << symmetry_classes[i].first->GetType() 
      << "): " << symmetry_classes[i].second << endl;
      }*/
    } else {
      CreateNewClassVector(symmetry_classes, tmp_classes);
      CountAndRenumberClasses(tmp_classes, nclasses2);
    }



    if (nclasses1 != nclasses2) {
      symmetry_classes = tmp_classes;
      return ExtendInvariants(symmetry_classes, breakChiralTies);
    }

    return nclasses1;
  }

/***************************************************************************
* FUNCTION: CalculateSymmetry
*
* DESCRIPTION:
*       Calculates a set of canonical symmetry identifiers for a molecule.
*       Atoms with the same symmetry ID are symmetrically equivalent.  By
*       "canonical", we mean it generates a repeatable labelling of the
*       atoms, i.e. the same fragment will get the same symmetry labels in
*       any molecule in which it occurs.
*
*       Vector is indexed from zero, corresponding to (atom->GetIdx() - 1).
*
*       The bit vector "_frag_atoms" specifies a fragment of the molecule,
*       where each bit represents the presence or absence of the atom in
*       the fragment.  Symmetry is computed as though the fragment is the
*       only part that exists.
***************************************************************************/
  
  int OBGraphSym::CalculateSymmetry(vector<unsigned int> &atom_sym_classes, bool breakChiralTies)
  {
    vector<unsigned int> vgi;
    vector<OBNodeBase*>::iterator j;
    OBAtom *atom;

    // Get vector of graph invariants.  These are the starting "symmetry classes".
    GetGIVector(vgi);

    // Create a vector-of-pairs, associating each atom with its Class ID.
    std::vector<std::pair<OBAtom*, unsigned int> > symmetry_classes;
    for (atom = _pmol->BeginAtom(j); atom; atom = _pmol->NextAtom(j)) {
      int idx = atom->GetIdx();
      if (_frag_atoms->BitIsOn(idx))
        symmetry_classes.push_back(pair<OBAtom*, unsigned int> (atom, vgi[idx-1]));
    }

    // The heart of the matter: Do extended sum-of-invariants until no further
    // changes are noted. 
    int nclasses = ExtendInvariants(symmetry_classes, breakChiralTies);

    // Convert to a vector indexed by Index
    // Atoms not in the fragment will have a value of OBGraphSym::NoSymmetryClass
    atom_sym_classes.clear();
    atom_sym_classes.resize(_pmol->NumAtoms(), OBGraphSym::NoSymmetryClass);
    for (unsigned int i = 0; i < symmetry_classes.size(); ++i) {
      atom_sym_classes[symmetry_classes.at(i).first->GetIndex()] = symmetry_classes.at(i).second;
    }

    // Store the symmetry classes in an OBPairData
    stringstream temp;
    vector<unsigned int>::iterator sym_iter = atom_sym_classes.begin();
    if (sym_iter != atom_sym_classes.end())
      temp << (*sym_iter++);
    for (; sym_iter != atom_sym_classes.end(); ++sym_iter)
      temp << " " << (*sym_iter);
  
    OBPairData *symData = new OBPairData;
    symData->SetAttribute("OpenBabel Symmetry Classes");
    symData->SetOrigin(local); //will not show as sdf or cml property
    symData->SetValue(temp.str());
    _pmol->SetData(symData);
 
    return nclasses;
  }

  int OBGraphSym::GetSymmetry(vector<unsigned int> &symmetry_classes, bool breakChiralTies)
  {
    cout << "ENTER GetSymmetry" << endl;
    ClearSymmetry(); // For the moment just recalculate the symmetry classes

    // Check to see whether we have already calculated the symmetry classes
    OBPairData *pd = dynamic_cast<OBPairData*>(_pmol->GetData("OpenBabel Symmetry Classes"));

    int nclasses = 0;
    if (!pd) {
      nclasses = CalculateSymmetry(symmetry_classes, breakChiralTies);
    } else {
      istringstream iss(pd->GetValue());
      symmetry_classes.clear();
      copy(istream_iterator<unsigned int>(iss),
           istream_iterator<unsigned int>(),
           back_inserter<vector<unsigned int> >(symmetry_classes));
      // Now find the number of unique elements
      vector<unsigned int> copy_sym = symmetry_classes;
      sort(copy_sym.begin(), copy_sym.end());
      vector<unsigned int>::iterator end_pos = unique(copy_sym.begin(), copy_sym.end()); // Requires sorted elements
      nclasses = end_pos - copy_sym.begin();
    }

    return nclasses;      
  }

  // Clears perceived symmetry
  void OBGraphSym::ClearSymmetry()
  {
    _pmol->DeleteData("OpenBabel Symmetry Classes");
  }

/***************************************************************************
* FUNCTION: CanonicalLabels
*
* DESCRIPTION:
*       Generates a canonical labeling of the atoms of a molecule, and as
*       a side benefit, returns the symmetry classes of the molecule.
*
*       To create a canonical labeling, we need every node to have a unique
*       label.  The canonical symmetry classes (see CalculateSymmetry(),
*       above) are a good start, but the atoms in each symmetry class are
*       still indistinguishable.  For writing a canonical string, we need
*       to create an arbitrary, but canonical (repeatable) distinction
*       between the atoms in each symmetry class -- "break the ties" in the
*       symmetry values.
*
*       To break ties, we sort into symetry-class order, double all class
*       IDs, then arbitrarily subtract one from the first repeated symmetry
*       class, thus breaking the tie (see Weininger et al).  With this new
*       set of symmetry classes, we repeat the extended-connectivity sums
*       to "spread" the broken symmetry class, and check again.  This is
*       repeated until all symmetry is gone and every atom has a unique
*       label.
*
* RETURNS:
*       canonical_labels - a vector indexed by [ OBAtom::GetIdx() - 1].
***************************************************************************/

  void OBGraphSym::CanonicalLabels(vector<unsigned int> &canonical_labels)
  {
    vector<pair<OBAtom*,unsigned int> > vp1, vp2;
    vector<OBNodeBase*>::iterator j;
    unsigned int nclass1, nclass2; //number of classes
    int i;

    int nfragatoms = _frag_atoms->CountBits();
    int natoms = _pmol->NumAtoms();

    std::vector<unsigned int> symmetry_classes;
    nclass1 = GetSymmetry(symmetry_classes, true);
    for (int i = 0; i < symmetry_classes.size(); ++i) {
      if (symmetry_classes.at(i) != OBGraphSym::NoSymmetryClass)
        vp1.push_back(
            pair<OBAtom*, unsigned int>(_pmol->GetAtom(i+1), symmetry_classes[i]) );
    }
    CountAndRenumberClasses(vp1, nclass1);

    /*cout << "BEFORE TieBreaker: nclass1 = " << nclass1 << ", nfragatoms = " << nfragatoms << "\n";
      for (int i = 0; i < vp1.size(); i++) {
      cout << vp1[i].first->GetIndex() << ": " << vp1[i].second << endl;
      }*/
    //DepictSymmetryClasses(_pmol, symmetry_classes, "img/NewBreakChiralTies");

    // The symmetry classes are the starting point for the canonical labels
    if (nclass1 < nfragatoms) {
      int tie_broken = 1;
      while (tie_broken) {
        tie_broken = 0;
        int last_rank = -1;
        for (i = 0; i < vp1.size(); i++) {
          vp1[i].second *= 2;             // Double symmetry classes
          if (vp1[i].second == last_rank && !tie_broken) {
            vp1[i-1].second -= 1;         // Break a tie
            tie_broken = 1;
          }
          last_rank = vp1[i].second;
        }
        if (tie_broken) {
          for (i = 0; i < 100;i++) {  //sanity check - shouldn't ever hit this number
            CreateNewClassVector(vp1, vp2);
            CountAndRenumberClasses(vp2, nclass2);
            vp1 = vp2;
            if (nclass1 == nclass2) break;
            nclass1 = nclass2;
          }
        } else {
          CountAndRenumberClasses(vp1, nclass1);  // no more ties - undo the doublings
        }
      }
    }

    /*cout << "AFTER TieBreaker: nclass1 = " << nclass1 << ", nfragatoms = " << nfragatoms << "\n";
      for (int i = 0; i < vp1.size(); i++) {
      cout << vp1[i].first->GetIndex() << ": " << vp1[i].second << endl;
      }*/

    canonical_labels.resize(_pmol->NumAtoms(), OBGraphSym::NoSymmetryClass);
    for (int i = 0; i < vp1.size(); ++i) {
      canonical_labels[vp1.at(i).first->GetIndex()] = vp1.at(i).second;
    }

    //DepictSymmetryClasses(_pmol, canonical_labels, "img/NewBreakChiralTies");

  }

  void OBGraphSym::CanonicalLabels(OBMol *mol, const std::vector<unsigned int> &symmetry_classes, 
      std::vector<unsigned int> &canonical_labels)
  {
    vector<pair<OBAtom*,unsigned int> > vp1, vp2;
    vector<OBNodeBase*>::iterator j;
    unsigned int nclass1, nclass2; //number of classes
    int i;

    int natoms = mol->NumAtoms();

    for (int i = 0; i < symmetry_classes.size(); ++i) {
      if (symmetry_classes.at(i) != OBGraphSym::NoSymmetryClass)
        vp1.push_back(
            pair<OBAtom*, unsigned int>(mol->GetAtom(i+1), symmetry_classes[i]) );
    }
    OBGraphSym::CountAndRenumberClasses(vp1, nclass1);

    /*cout << "BEFORE TieBreaker: nclass1 = " << nclass1 << ", nfragatoms = " << nfragatoms << "\n";
      for (int i = 0; i < vp1.size(); i++) {
      cout << vp1[i].first->GetIndex() << ": " << vp1[i].second << endl;
      }*/

    // The symmetry classes are the starting point for the canonical labels
    if (nclass1 < natoms) {
      int tie_broken = 1;
      while (tie_broken) {
        tie_broken = 0;
        int last_rank = -1;
        for (i = 0; i < vp1.size(); i++) {
          vp1[i].second *= 2;             // Double symmetry classes
          if (vp1[i].second == last_rank && !tie_broken) {
            vp1[i-1].second -= 1;         // Break a tie
            tie_broken = 1;
          }
          last_rank = vp1[i].second;
        }
        if (tie_broken) {
          for (i = 0; i < 100;i++) {  //sanity check - shouldn't ever hit this number
            OBGraphSym::CreateNewClassVector(mol, vp1, vp2);
            OBGraphSym::CountAndRenumberClasses(vp2, nclass2);
            vp1 = vp2;
            if (nclass1 == nclass2) break;
            nclass1 = nclass2;
          }
        } else {
          OBGraphSym::CountAndRenumberClasses(vp1, nclass1);  // no more ties - undo the doublings
        }
      }
    }

    /*cout << "AFTER TieBreaker: nclass1 = " << nclass1 << ", nfragatoms = " << nfragatoms << "\n";
      for (int i = 0; i < vp1.size(); i++) {
      cout << vp1[i].first->GetIndex() << ": " << vp1[i].second << endl;
      }*/

    canonical_labels.resize(mol->NumAtoms(), OBGraphSym::NoSymmetryClass);
    for (int i = 0; i < vp1.size(); ++i) {
      canonical_labels[vp1.at(i).first->GetIndex()] = vp1.at(i).second;
    }
 
  }
 
      

  int OBGraphSym::Iterate(vector<unsigned int> &symClasses)
  {
    // Create a vector-of-pairs, associating each atom with its Class ID.
    vector<OBAtom*>::iterator j;
    std::vector<std::pair<OBAtom*, unsigned int> > symmetry_classes;
    for (OBAtom *atom = _pmol->BeginAtom(j); atom; atom = _pmol->NextAtom(j)) {
      int idx = atom->GetIdx();
      if (_frag_atoms->BitIsOn(idx))
        symmetry_classes.push_back(pair<OBAtom*, unsigned int> (atom, symClasses[idx-1]));
    }

    // The heart of the matter: Do extended sum-of-invariants until no further
    // changes are noted. 
    int nclasses = ExtendInvariants(symmetry_classes, false);

    // Convert to a vector indexed by Index
    // Atoms not in the fragment will have a value of OBGraphSym::NoSymmetryClass
    symClasses.clear();
    symClasses.resize(_pmol->NumAtoms(), OBGraphSym::NoSymmetryClass);
    for (unsigned int i = 0; i < symmetry_classes.size(); ++i) {
      symClasses[symmetry_classes.at(i).first->GetIndex()] = symmetry_classes.at(i).second;
    }

    return nclasses;
  }


} // namespace OpenBabel

//! \file graphsym.cpp
//! \brief XXXX
