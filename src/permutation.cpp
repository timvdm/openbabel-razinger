#include <openbabel/permutation.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/graphsym.h>

#include <graph.hh>

#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>

namespace OpenBabel {

  void callback(void *param, unsigned int n, const unsigned int *aut)
  {
    // construct vector with elements
    std::vector<unsigned int> elements;
    //cout << "Generator: ";
    for (unsigned int i = 0; i < n; ++i) {
      //cout << aut[i] + 1 << " ";
      elements.push_back(aut[i]+1);
    }
    //cout << endl;

    // add the generator
    OBPermutationGroup *generators = static_cast<OBPermutationGroup*>(param);
    generators->Add(OBPermutation(elements));
  }

  /**
   * Add the permutations' inverses to group @p G. These inverses automatically
   * belong to the automorphism group without the need of further testing.
   */
  void addInverses(OBPermutationGroup &G)
  {
    //cout << "addInverses..." << endl;
    unsigned int size = G.Size();
    for (unsigned int i = 0; i < size; ++i) {
      int n = G.At(i).map.size();
      if (!n)
        continue;
      Eigen::MatrixXi P = Eigen::MatrixXi::Zero(n, n);
      for (int j = 0; j < n; ++j) {
        P(j, G.At(i).map.at(j)-1) = 1;
      }
 
      OBPermutation inv_p(P.transpose());
      //cout << "matrix:" << endl; cout << G.permutations.at(i).matrix() << endl;
      //inv_p.print();
      if (!G.Contains(inv_p)) {
        //cout << "--------------> found inverse" << endl;
        G.Add(inv_p); 
      }
    }
  }

  /**
   * Add the all permutation products to group @p G. These inverses automatically
   * belong to the automorphism group without the need of further testing.
   */
  void addProducts(OBPermutationGroup &G)
  {
    //cout << "addProducts..." << endl;
    for (unsigned int i = 0; i < G.Size(); ++i) {
      for (unsigned int j = 0; j < G.Size(); ++j) {
        if (i >= j)
          continue;

        OBPermutation p = G.At(i) * G.At(j);
        //p.print();
        if (!G.Contains(p)) {
          //cout << "-------------> found product" << endl;
          G.Add(p);
        }
      }
    }
  }

  OBPermutationGroup FindAutomorphisms(OpenBabel::OBMol *obmol, const std::vector<unsigned int> &symClasses)
  {
    // construct the bliss graph
    bliss::Graph g;
    std::map<OpenBabel::OBAtom*, unsigned int> atom2index;
    FOR_ATOMS_OF_MOL (atom, obmol) {
      atom2index[&*atom] = g.add_vertex(symClasses.at(atom->GetIndex()));
    }
    FOR_BONDS_OF_MOL (bond, obmol) {
      g.add_edge(atom2index[bond->GetBeginAtom()], atom2index[bond->GetEndAtom()]);
    }

    // use bliss to get the automorphism group generators
    OBPermutationGroup generators;
    bliss::Stats stats;
    g.find_automorphisms(stats, &callback, &generators);

    //cout << "# Automorphisms:" << stats.group_size_approx << endl;
    unsigned long nAut = stats.group_size_approx;

    // construct the automorphism group
    OBPermutationGroup G;

    // add identity permutation
    std::vector<unsigned int> eElements;
    for (unsigned int i = 0; i < obmol->NumAtoms(); ++i)
      eElements.push_back(i+1);
    OBPermutation e(eElements);
    if (!G.Contains(e))
      G.Add(e);

    // add the generators
    for (unsigned int i = 0; i < generators.Size(); ++i) {
      if (!G.Contains(generators.At(i)))
        G.Add(generators.At(i));
    }

    // loop and add inverses and products until no more automorphisms are found
    unsigned int counter = 0;
    unsigned int lastSize = 0;
    while (G.Size() != lastSize) {
      lastSize = G.Size();

      addInverses(G);
      addProducts(G);
    
      counter++;
      if (counter > 1000)
        break;
    }

    /*
    if (nAut != G.size())
      cout << "ERROR: Not all " << nAut << " automorphisms are found!" << endl;
    else 
      cout << "SUCCESS: all " << nAut << " automorphisms are found!" << endl;
      */
    
    return G;
  }

} // namespace
