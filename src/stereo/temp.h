namespace OpenBabel {

  void contractPermutations(const PermutationGroup &G, PermutationGroup &Gc, const std::vector<unsigned int> &stereoAtoms)
  {
    // contract the permutations (i.e. remove non-stereogenic atoms)
    //std::cout << "Contracted Permutations:" << std::endl;
    for (unsigned int g = 0; g < G.permutations.size(); ++g) {
      Permutation p;
      for (unsigned int j = 0; j < G.at(g).map.size(); ++j) {
        if (std::find(stereoAtoms.begin(), stereoAtoms.end(), G.at(g).map.at(j)) != stereoAtoms.end())
          p.map.push_back(G.at(g).map.at(j));
      }
      //p.print();
      if (p.map.size())
        Gc.add(p);
    } 
  }

  void signedPermutationMatrices(const PermutationGroup &Gc, const std::vector<Eigen::VectorXi> &stereoIndexVectors, 
      std::vector<Eigen::MatrixXi> &signedMatrices, const std::vector<unsigned int> &tetrahedralAtoms, int numCisTrans, int n)
  {
    int numTetrahedral = tetrahedralAtoms.size();
    //cout << "signedPermutationMatrices" << endl;
    for (unsigned int g = 0; g < Gc.permutations.size(); ++g) {
      const Permutation &p = Gc.permutations.at(g);

      Eigen::MatrixXi Ps = Eigen::MatrixXi::Zero(n, n);
 
      std::vector<unsigned int> tetrahedralMap, cistransMap;
      for (unsigned int i = 0; i < p.map.size(); ++i)
        if (std::find(tetrahedralAtoms.begin(), tetrahedralAtoms.end(), p.map.at(i)) != tetrahedralAtoms.end())
          tetrahedralMap.push_back(p.map.at(i));
        else
          cistransMap.push_back(p.map.at(i));
     
      if (numTetrahedral) {
        Permutation p2(tetrahedralMap);
        Eigen::MatrixXi P = p2.matrix();
        for (unsigned int j = 0; j < numTetrahedral; ++j)
          for (unsigned int k = 0; k < numTetrahedral; ++k)
            Ps(k, j) = P(k, j) * stereoIndexVectors.at(g)[k];
      }
      if (numCisTrans) {
        Permutation p2(cistransMap);
        Eigen::MatrixXi P = p2.matrix();
        for (unsigned int j = numTetrahedral; j < numTetrahedral + numCisTrans; ++j)
          for (unsigned int k = numTetrahedral; k < numTetrahedral + numCisTrans; ++k)
            Ps(k, j) = P(k - numTetrahedral, j - numTetrahedral) * stereoIndexVectors.at(g)[k];
      }

      signedMatrices.push_back( Ps );

      //std::cout << "signed permutation matrix:" << std::endl;
      //std::cout << signedMatrices.back() << std::endl;
    }
  }
        
  Eigen::VectorXi multiplyMatrixByRow(const Eigen::MatrixXi &matrix, const Eigen::VectorXi &row)
  {
    //Eigen::VectorXi Pirj = matrix * row/*.transpose()*/;
    Eigen::VectorXi Pirj = matrix * row;
    return Pirj;
  }

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

  void createStereoIndexVectors(std::vector<Eigen::VectorXi> &stereoIndexVectors, const PermutationGroup &Gc,
      const PermutationGroup &G, const std::vector<StereogenicUnit> &stereoUnits, OBMol *mol, 
      const std::vector<unsigned int> &symmetry_classes, int n)
  {      
    //cout << "createStereoIndexVectors" << endl;
    //
    // Construct the stereoindex vectors, one for each automorphism:
    //
    // - true-stereo atoms: +1
    // - para-stereo atoms:
    //   - identical ligands not interchanged: +1 (or even number of permutation)
    //   - identical liganss interchanged: -1 (or odd number of permutations)
    //
    for (unsigned int g = 0; g < Gc.permutations.size(); ++g) {
      const Permutation &p = Gc.permutations.at(g);
      const Permutation &pFull = G.permutations.at(g);
      //std::cout << "permutation: "; p.print();
      //std::cout << "matrix: " << endl;
      //std::cout << p.matrix() << endl;

      Eigen::VectorXi stereoIndex(n);
      for (unsigned int c = 0; c < n; ++c)
        stereoIndex[c] = 0;
      //
      // make sure all tetrahedral centers go before cis/trans atoms
      //
      
      // make sure we respect the order of the atoms when constructing the stereoindex vectors
      std::vector<unsigned int> indexes;
      for (vector<StereogenicUnit>::const_iterator unit = stereoUnits.begin(); unit != stereoUnits.end(); ++unit) {
        if (unit->type == OBStereo::Tetrahedral) {
          indexes.push_back(mol->GetAtomById(unit->id)->GetIndex());
        }
      }
      std::sort(indexes.begin(), indexes.end());
      
      for (vector<StereogenicUnit>::const_iterator unit = stereoUnits.begin(); unit != stereoUnits.end(); ++unit) {
        if (unit->type == OBStereo::Tetrahedral) {
          OBAtom *atom = mol->GetAtomById(unit->id);
          //std::cout << "atom id: " << atom->GetId() << std::endl;

          //compute the insert position
          unsigned int insertpos;
          for (unsigned int i = 0; i < indexes.size(); ++i)
            if (indexes.at(i) == atom->GetIndex()) {
              insertpos = i;
              break;
            }

          if (!unit->para) {
            //cout << "-> true stereocenter: id = " << unit->id << endl;
            stereoIndex[insertpos] = 1;
          } else {
            //cout << "-> para stereocenter: id = " << unit->id << endl;
            // need to count the number of permutations for
            stereoIndex[insertpos] = getStereoIndex(atom, pFull, symmetry_classes);
          }
        }
      }

      int numTetrahedral = indexes.size();
      indexes.clear();

      for (vector<StereogenicUnit>::const_iterator unit = stereoUnits.begin(); unit != stereoUnits.end(); ++unit) {
        if (unit->type == OBStereo::CisTrans) {
          OBBond *bond = mol->GetBondById(unit->id);
          indexes.push_back(bond->GetBeginAtom()->GetIndex());
          indexes.push_back(bond->GetEndAtom()->GetIndex());
        }
      } 
      std::sort(indexes.begin(), indexes.end());
      for (vector<StereogenicUnit>::const_iterator unit = stereoUnits.begin(); unit != stereoUnits.end(); ++unit) {
        if (unit->type == OBStereo::CisTrans) {
          OBBond *bond = mol->GetBondById(unit->id);
          //std::cout << "bond id: " << bond->GetId() << std::endl;

          //compute the insert position
          unsigned int insertpos_begin, insertpos_end;
          for (unsigned int i = 0; i < indexes.size(); ++i)
            if (indexes.at(i) == bond->GetBeginAtom()->GetIndex()) {
              insertpos_begin = i + numTetrahedral;
            } else if (indexes.at(i) == bond->GetEndAtom()->GetIndex()) {
              insertpos_end = i + numTetrahedral;
            }


          if (!unit->para) {
            stereoIndex[insertpos_begin] = 1;
            stereoIndex[insertpos_end] = 1;
          } else {
            // need to count the number of permutations for
            stereoIndex[insertpos_begin] = getStereoIndex(bond->GetBeginAtom(), pFull, symmetry_classes);
            stereoIndex[insertpos_end] = getStereoIndex(bond->GetEndAtom(), pFull, symmetry_classes);
          }
        }
      }

      //std::cout << "stereoIndex =" << std::endl;
      //std::cout << stereoIndex << std::endl;
      stereoIndexVectors.push_back(stereoIndex);
      //cout << "stereoindex:" << endl << stereoIndex << endl;
 
    }

    // debug: print out stereoindex vectors
    /*
    cout << "StereoIndexVectors:" << endl;
    for (unsigned int idx = 0; idx < stereoIndexVectors.size(); ++idx)
      cout << "stereoindex:" << endl << stereoIndexVectors.at(idx) << endl;
      */
  }


}
