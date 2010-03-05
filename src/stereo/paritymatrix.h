#include <Eigen/Core>

namespace OpenBabel {

  // Helper struct for stereoParityMatrix
  struct ParityGenerator
  {
    ParityGenerator(unsigned int _power) : power(_power), counter(0), retval(false)
    {
    }

    bool inc()
    {
      if (counter >= power) {
        retval = !retval;
        counter = 0;
      }

      counter++;
      return retval;
    }

    unsigned int power, counter;
    bool retval;
  };


  /**
   * Get the stereo parity matrix. This matrix has @p n columns and 2^n rows. 
   * The rows represent the row index in binary where 0=1 and 1=-1.
   *
   @verbatim
    1  1  1  1
    1  1  1 -1
    1  1 -1  1
    1  1 -1 -1
    1 -1  1  1
    1 -1  1 -1
   ...
   -1 -1 -1 -1
   @endverbatim
   *
   * @param n The number of stereocenters
   */
  Eigen::MatrixXi stereoParityMatrix(unsigned int n)
  {
    unsigned int N = pow(2, n);
    std::vector<ParityGenerator> generators;
    for (unsigned int j = 0; j < n; ++j) {
      generators.push_back(ParityGenerator(pow(2, n - j - 1)));
    }

    Eigen::MatrixXi parityMatrix(Eigen::MatrixXi::Zero(N, n));
    for (unsigned int i = 0; i < N; ++i) {
      for (unsigned int j = 0; j < n; ++j) {
        if (generators[j].inc())
          parityMatrix(i,j) = -1;
        else
          parityMatrix(i,j) = 1;
      }
    }

    return parityMatrix;
  }

}
