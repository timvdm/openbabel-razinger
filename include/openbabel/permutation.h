#ifndef OB_PERMUTATION_H
#define OB_PERMUTATION_H

#include <openbabel/babelconfig.h>

#include <map>
#include <vector>
#include <algorithm>

#include <Eigen/Core>

namespace OpenBabel {

  class OBMol;

  ///@addtogroup permutation Permutations 
  ///@{

  /**
   * @brief Class representing single permutation
   *
   * The OBPermutation class represents a single permutation of a series of 
   * numbers. The class is meant to be used with the findAutomorphisms
   * function (other uses might exist though).
   */
  struct OBAPI OBPermutation
  {
    /**
     * Default constructor.
     */
    OBPermutation()
    {}

    /**
     * Constructor taking a @p map as argument.
     */
    OBPermutation(const std::vector<unsigned int> &_map) : map(_map)
    {}

    /**
     * Constructor taking a @p matrix as argument.
     */
    OBPermutation(const Eigen::MatrixXi &matrix)
    {
      SetMatrix(matrix);
    }

    /**
     * Copy constructor.
     */
    OBPermutation(const OBPermutation &other)
    {
      this->map = other.map;
    }

    /**
     * Print the permutation to std::cout in shortened notation.
     */
    void Print() const
    {
      std::vector<unsigned int>::const_iterator i;
      for (i = map.begin(); i != map.end(); ++i)
        std::cout << *i << " ";
      std::cout << std::endl;    
    }

    OBPermutation Apply(const OBPermutation &input) const
    {
      OBPermutation p;
      if (input.map.size() != map.size())
        return p;
      p.map.resize(map.size());

      unsigned int index = 0;
      std::vector<unsigned int>::const_iterator i;
      for (i = map.begin(); i != map.end(); ++i, ++index) {
        p.map[index] = input.map.at(*i-1);
      }

      return p;   
    }

    /**
     * Compute the permutation matrix for this OBPermutation.
     */
    Eigen::MatrixXi GetMatrix() const
    {
      int n = map.size();
      Eigen::MatrixXi P = Eigen::MatrixXi::Zero(n, n);

      // make sure we can create matrices for contracted permutations 
      // (i.e. permutations where the values don't go from 1 to N, there
      // may be gaps)
      std::map<unsigned int, unsigned int> renum;
      unsigned int j = 0;
      unsigned int k = 1;
      while (renum.size() < n) {
        if (std::find(map.begin(), map.end(), k) != map.end()) {
          renum[k] = j;
          j++;
        }
      
        k++;
      }

      for (int i = 0; i < n; ++i) {
        P(renum[map.at(i)], i) = 1;
      }

      return P;
    }

    /**
     * Set the permutation matrix for this permutation. The matrix is
     * converted to a shortened notation permutation.
     */
    void SetMatrix(const Eigen::MatrixXi &m)
    {
      if (m.rows() != m.cols())
        return;

      int size = m.rows();
      map.resize(size);

      for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j) {
          if (m(i,j))
            map[i] = j + 1;
        }
    }

    unsigned long NumInversions() const
    {
      std::vector<unsigned int> invVec; // the inversion vector
      std::vector<unsigned int>::const_iterator i, j;
      for (i = map.begin(); i != map.end(); ++i) {
        int e = 0; // ith element
        // loop over elements to the right
        for (j = i; j != map.end(); ++j)
          // increment e if element to the right is lower
          if (*j < *i)
            e++;

        invVec.push_back(e);
      }

      unsigned long sum = 0;
      for (std::vector<unsigned int>::iterator k = invVec.begin(); k != invVec.end(); ++k)
        sum += *k;

      return sum;
    }


    /**
     * Multiply two permutation matrices and return the resulting permutation.
     * This is the same as applying the two permutations consecutively.
     */
    OBPermutation operator*(const OBPermutation &rhs) const
    {
      return Apply(rhs);
    }

    /**
     * Equality operator. Two permutations are equal if their shortened notations 
     * are the same (compared element-by-element).
     */
    bool operator==(const OBPermutation &rhs) const
    {
      if (map.size() != rhs.map.size())
        return false;

      std::vector<unsigned int>::const_iterator i1 = map.begin();
      std::vector<unsigned int>::const_iterator i2 = rhs.map.begin();
      for (; i1 != map.end(); ++i1, ++i2)
        if (*i1 != *i2)
          return false;

      return true;  
    }

    std::vector<unsigned int> map; //!< The actual mapping in shortened notation
  };

  /**
   * @brief Class representing a group of permutations
   *
   * The OBPermutationGroup class represents a group of OBPermutation objects.
   * The biggest benefit of this separate class over a std::vector is the 
   * contains function. This class is the return type of the findAutomorphisms
   * function.
   */
  struct OBAPI OBPermutationGroup
  {
    /**
     * Default constructor.
     */
    OBPermutationGroup()
    {}

    /**
     * Constructor taking vector of permutations as argument.
     */
    OBPermutationGroup(const std::vector<OBPermutation> &_permutations) : permutations(_permutations)
    {}

    /**
     * @return The number of permutations in this group.
     */
    unsigned int Size() const
    {
      return permutations.size();
    }

    /**
     * Add permutation @p p to this permutation group.
     */
    void Add(const OBPermutation &p)
    {
      permutations.push_back(p);
    }

    /**
     * @return A constant reference to the @p index-th permutation.
     */
    const OBPermutation& At(unsigned int index) const
    {
      return permutations.at(index);
    }

    /**
     * @return True if this permutation group contains permutation @p p.
     */
    bool Contains(const OBPermutation &p) const
    {
      std::vector<OBPermutation>::const_iterator i;
      for (i = permutations.begin(); i != permutations.end(); ++i)
        if (*i == p)
          return true;
      return false; 
    }

    std::vector<OBPermutation> permutations; //!< The actual permutations in the group
  };


  /**
   * Find all automorphisms of a molecule. The specified symmetry classes partition 
   * the atoms in groups of topological equivalent atoms. The permutations in the 
   * automorphisms group are all possible permutations which preserve connectivity.
   * Mathematically this can be expressed using the adjacency matrix (\f$A\f$) and
   * permutation matrix (\f$P\f$) in the equation:
     \f[
     P^T A P = A     
     \f]
   */
  OBAPI OBPermutationGroup FindAutomorphisms(OpenBabel::OBMol *obmol, const std::vector<unsigned int> &symClasses);


  ///@}

} // namespace

#endif
