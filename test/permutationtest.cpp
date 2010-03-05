#include "obtest.h"
#include <openbabel/permutation.h>

using namespace std;
using namespace OpenBabel;

void testConstructors()
{
  // default constructor
  Permutation p1;

  std::vector<unsigned int> map(3);
  map[0] = 1;
  map[1] = 2;
  map[2] = 3;
  Permutation p2(map);
  OB_REQUIRE( p2.map.size() == 3 );
  OB_REQUIRE( p2.map[0] == 1 );
  OB_REQUIRE( p2.map[1] == 2 );
  OB_REQUIRE( p2.map[2] == 3 );
}

void testApply()
{
  std::vector<unsigned int> orig(4);
  orig[0] = 1;
  orig[1] = 2;
  orig[2] = 3;
  orig[3] = 4;
  std::vector<unsigned int> map(4);
  map[0] = 4;
  map[1] = 1;
  map[2] = 3;
  map[3] = 2;

  Permutation p(map);
  Permutation result = p.apply(orig);

  OB_REQUIRE( result.map.size() == 4 );
  OB_REQUIRE( result.map[0] == 4 );
  OB_REQUIRE( result.map[1] == 1 );
  OB_REQUIRE( result.map[2] == 3 );
  OB_REQUIRE( result.map[3] == 2 );
}

void testMatrix()
{
  std::vector<unsigned int> map(6);
  map[0] = 2;
  map[1] = 1;
  map[2] = 4;
  map[3] = 3;
  map[4] = 6;
  map[5] = 5;

  Permutation p(map);

  Eigen::MatrixXi m = p.matrix();

  // 0 1 0 0 0 0
  // 1 0 0 0 0 0
  // 0 0 0 1 0 0
  // 0 0 1 0 0 0
  // 0 0 0 0 0 1
  // 0 0 0 0 1 0
  OB_REQUIRE( m(0,0) == 0 );
  OB_REQUIRE( m(0,1) == 1 );
  OB_REQUIRE( m(0,2) == 0 );
  OB_REQUIRE( m(0,3) == 0 );
  OB_REQUIRE( m(0,4) == 0 );
  OB_REQUIRE( m(0,5) == 0 );

  OB_REQUIRE( m(1,0) == 1 );
  OB_REQUIRE( m(1,1) == 0 );
  OB_REQUIRE( m(1,2) == 0 );
  OB_REQUIRE( m(1,3) == 0 );
  OB_REQUIRE( m(1,4) == 0 );
  OB_REQUIRE( m(1,5) == 0 );

  OB_REQUIRE( m(2,0) == 0 );
  OB_REQUIRE( m(2,1) == 0 );
  OB_REQUIRE( m(2,2) == 0 );
  OB_REQUIRE( m(2,3) == 1 );
  OB_REQUIRE( m(2,4) == 0 );
  OB_REQUIRE( m(2,5) == 0 );

  OB_REQUIRE( m(3,0) == 0 );
  OB_REQUIRE( m(3,1) == 0 );
  OB_REQUIRE( m(3,2) == 1 );
  OB_REQUIRE( m(3,3) == 0 );
  OB_REQUIRE( m(3,4) == 0 );
  OB_REQUIRE( m(3,5) == 0 );

  OB_REQUIRE( m(4,0) == 0 );
  OB_REQUIRE( m(4,1) == 0 );
  OB_REQUIRE( m(4,2) == 0 );
  OB_REQUIRE( m(4,3) == 0 );
  OB_REQUIRE( m(4,4) == 0 );
  OB_REQUIRE( m(4,5) == 1 );

  OB_REQUIRE( m(5,0) == 0 );
  OB_REQUIRE( m(5,1) == 0 );
  OB_REQUIRE( m(5,2) == 0 );
  OB_REQUIRE( m(5,3) == 0 );
  OB_REQUIRE( m(5,4) == 1 );
  OB_REQUIRE( m(5,5) == 0 );

  Permutation p2(m);
  OB_REQUIRE( p2.map.size() == 6 );
  OB_REQUIRE( p2.map[0] == 2 );
  OB_REQUIRE( p2.map[1] == 1 );
  OB_REQUIRE( p2.map[2] == 4 );
  OB_REQUIRE( p2.map[3] == 3 );
  OB_REQUIRE( p2.map[4] == 6 );
  OB_REQUIRE( p2.map[5] == 5 );
}

void testNumInversions()
{
}


  
int main()
{
  testConstructors();
  testApply();
  testMatrix();
  testNumInversions();
  
  return 0;
}
