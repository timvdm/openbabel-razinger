#include "obtest.h"
#include <openbabel/stereo/stereoisomer.h>
#include <openbabel/stereo/stereo.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/permutation.h>
#include <openbabel/graphsym.h>

#include "../src/stereo/paritymatrix.h"
#include "../src/stereo/temp.h"

using namespace std;
using namespace OpenBabel;

std::string GetFilename(const std::string &filename)
{
  string path = TESTDATADIR + filename;
  return path;
}

static int testCount = 0;

bool doStereoisomerTest(OBMol &mol, int numEnantiomerPairs, int numDiastereomers)
{
  OBStereoisomer isomers(&mol);
  
  OB_COMPARE(isomers.numEnantiomerPairs(), numEnantiomerPairs);
  OB_COMPARE(isomers.numDiastereomers(), numDiastereomers);

  testCount++;

  return (isomers.numEnantiomerPairs() == numEnantiomerPairs) && (isomers.numDiastereomers() == numDiastereomers);
}


/**
 * Test detection of stereoisomers
 */
void test_Stereoisomer(int n = 0, const std::string &file = std::string())
{
  OBMol mol;
  OBConversion conv;
  OB_ASSERT( conv.SetInFormat("mol") );

  if (!file.empty()) {
    cout << "file = " << file << endl;
    OBFormat *format = conv.FormatFromExt(file.c_str());
    OB_REQUIRE( format );
    OB_ASSERT( conv.SetInFormat(format) );
    OB_ASSERT( conv.ReadFile(&mol, file) );
    OBStereoisomer isomers(&mol);
    return;
  }

  if (n > 0 && n < 70) {
    cout << "Razinger paper, fig. 7: structure " << n << endl;
    stringstream ss;
    ss << "stereo/razinger_fig7_" << n << ".mol";
    OB_ASSERT( conv.ReadFile(&mol, GetFilename(ss.str())) );
    OBStereoisomer isomers(&mol);
    return;
  }

  cout << "Razinger paper, fig. 3" << endl;
  OB_REQUIRE( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig3.mol")) );
  OB_ASSERT( doStereoisomerTest(mol, 3, 4) );
  
  cout << "Razinger paper, fig. 6" << endl;
  OB_REQUIRE( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig6.mol")) );
  OB_ASSERT( doStereoisomerTest(mol, 2, 3) );
 

  std::vector<int> failed;
  /*
   * J. Chem. Inf. Comput. Sci., Vol. 33, No. 6, 1993
   *
   * Figure 7. Test compounds
   */
  cout << "Razinger paper, fig. 7: structure 1" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_1.mol")) );
  { bool result = doStereoisomerTest(mol, 1, 0); OB_ASSERT( result ); if (!result) failed.push_back(1); }

  cout << "Razinger paper, fig. 7: structure 2" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_2.mol")) );
  { bool result = doStereoisomerTest(mol, 8, 4); OB_ASSERT( result ); if (!result) failed.push_back(2); }

  cout << "Razinger paper, fig. 7: structure 3" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_3.mol")) );
  { bool result = doStereoisomerTest(mol, 6, 4); OB_ASSERT( result ); if (!result) failed.push_back(3); }

  cout << "Razinger paper, fig. 7: structure 4" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_4.mol")) );
  { bool result = doStereoisomerTest(mol, 4, 0); OB_ASSERT( result ); if (!result) failed.push_back(4); }

  cout << "Razinger paper, fig. 7: structure 5" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_5.mol")) );
  { bool result = doStereoisomerTest(mol, 2, 0); OB_ASSERT( result ); if (!result) failed.push_back(5); }

  cout << "Razinger paper, fig. 7: structure 6" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_6.mol")) );
  { bool result = doStereoisomerTest(mol, 2, 1); OB_ASSERT( result ); if (!result) failed.push_back(6); }

  cout << "Razinger paper, fig. 7: structure 7" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_7.mol")) );
  { bool result = doStereoisomerTest(mol, 4, 0); OB_ASSERT( result ); if (!result) failed.push_back(7); }

  cout << "Razinger paper, fig. 7: structure 8" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_8.mol")) );
  { bool result = doStereoisomerTest(mol, 5, 0); OB_ASSERT( result ); if (!result) failed.push_back(8); }

  cout << "Razinger paper, fig. 7: structure 9" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_9.mol")) );
  { bool result = doStereoisomerTest(mol, 8, 0); OB_ASSERT( result ); if (!result) failed.push_back(9); }

  cout << "Razinger paper, fig. 7: structure 10" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_10.mol")) );
  { bool result = doStereoisomerTest(mol, 16, 0); OB_ASSERT( result ); if (!result) failed.push_back(10); }

  cout << "Razinger paper, fig. 7: structure 11" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_11.mol")) );
  { bool result = doStereoisomerTest(mol, 5, 0); OB_ASSERT( result ); if (!result) failed.push_back(11); }

  cout << "Razinger paper, fig. 7: structure 12" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_12.mol")) );
  { bool result = doStereoisomerTest(mol, 16, 0); OB_ASSERT( result ); if (!result) failed.push_back(12); }

  cout << "Razinger paper, fig. 7: structure 13" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_13.mol")) );
  { bool result = doStereoisomerTest(mol, 3, 4); OB_ASSERT( result ); if (!result) failed.push_back(13); }

  cout << "Razinger paper, fig. 7: structure 14" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_14.mol")) );
  { bool result = doStereoisomerTest(mol, 1, 1); OB_ASSERT( result ); if (!result) failed.push_back(14); }

  cout << "Razinger paper, fig. 7: structure 15" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_15.mol")) );
  { bool result = doStereoisomerTest(mol, 1, 2); OB_ASSERT( result ); if (!result) failed.push_back(15); }

  cout << "Razinger paper, fig. 7: structure 16" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_16.mol")) );
  { bool result = doStereoisomerTest(mol, 1, 1); OB_ASSERT( result ); if (!result) failed.push_back(16); }

  cout << "Razinger paper, fig. 7: structure 17" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_17.mol")) );
  { bool result = doStereoisomerTest(mol, 1, 1); OB_ASSERT( result ); if (!result) failed.push_back(17); }

  cout << "Razinger paper, fig. 7: structure 18" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_18.mol")) );
  { bool result = doStereoisomerTest(mol, 0, 2); OB_ASSERT( result ); if (!result) failed.push_back(18); }

  cout << "Razinger paper, fig. 7: structure 19" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_19.mol")) );
  { bool result = doStereoisomerTest(mol, 1, 2); OB_ASSERT( result ); if (!result) failed.push_back(19); }

  cout << "Razinger paper, fig. 7: structure 20" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_20.mol")) );
  { bool result = doStereoisomerTest(mol, 4, 2); OB_ASSERT( result ); if (!result) failed.push_back(20); }

  cout << "Razinger paper, fig. 7: structure 21" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_21.mol")) );
  { bool result = doStereoisomerTest(mol, 6, 4); OB_ASSERT( result ); if (!result) failed.push_back(21); }

  cout << "Razinger paper, fig. 7: structure 22" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_22.mol")) );
  { bool result = doStereoisomerTest(mol, 0, 4); OB_ASSERT( result ); if (!result) failed.push_back(22); }

  cout << "Razinger paper, fig. 7: structure 23" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_23.mol")) );
  { bool result = doStereoisomerTest(mol, 0, 5); OB_ASSERT( result ); if (!result) failed.push_back(23); }

  cout << "Razinger paper, fig. 7: structure 24" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_24.mol")) );
  { bool result = doStereoisomerTest(mol, 4, 2); OB_ASSERT( result ); if (!result) failed.push_back(24); }

  cout << "Razinger paper, fig. 7: structure 25" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_25.mol")) );
  { bool result = doStereoisomerTest(mol, 16, 7); OB_ASSERT( result ); if (!result) failed.push_back(25); }

  cout << "Razinger paper, fig. 7: structure 26" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_26.mol")) );
  { bool result = doStereoisomerTest(mol, 6, 8); OB_ASSERT( result ); if (!result) failed.push_back(26); }

  cout << "Razinger paper, fig. 7: structure 27" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_27.mol")) );
  { bool result = doStereoisomerTest(mol, 2, 0); OB_ASSERT( result ); if (!result) failed.push_back(27); }

  cout << "Razinger paper, fig. 7: structure 28" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_28.mol")) );
  { bool result = doStereoisomerTest(mol, 0, 2); OB_ASSERT( result ); if (!result) failed.push_back(28); }

  cout << "Razinger paper, fig. 7: structure 29" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_29.mol")) );
  { bool result = doStereoisomerTest(mol, 1, 2); OB_ASSERT( result ); if (!result) failed.push_back(29); }

  cout << "Razinger paper, fig. 7: structure 30" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_30.mol")) );
  { bool result = doStereoisomerTest(mol, 0, 2); OB_ASSERT( result ); if (!result) failed.push_back(30); }

  /*
  cout << "Razinger paper, fig. 7: structure 31" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_31.mol")) );
  { bool result = doStereoisomerTest(mol, 1, 0) );

  cout << "Razinger paper, fig. 7: structure 32" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_32.mol")) );
  { bool result = doStereoisomerTest(mol, 0, 2) );

  cout << "Razinger paper, fig. 7: structure 33" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_33.mol")) );
  { bool result = doStereoisomerTest(mol, 1, 0) );
  */
  
  cout << "Razinger paper, fig. 7: structure 34" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_34.mol")) );
  { bool result = doStereoisomerTest(mol, 1, 2); OB_ASSERT( result ); if (!result) failed.push_back(34); }
  
  cout << "Razinger paper, fig. 7: structure 35" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_35.mol")) );
  { bool result = doStereoisomerTest(mol, 0, 3); OB_ASSERT( result ); if (!result) failed.push_back(35); }

  cout << "Razinger paper, fig. 7: structure 36" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_36.mol")) );
  { bool result = doStereoisomerTest(mol, 2, 2); OB_ASSERT( result ); if (!result) failed.push_back(36); }

  cout << "Razinger paper, fig. 7: structure 37" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_37.mol")) );
  { bool result = doStereoisomerTest(mol, 2, 2); OB_ASSERT( result ); if (!result) failed.push_back(37); }

  cout << "Razinger paper, fig. 7: structure 38" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_38.mol")) );
  { bool result = doStereoisomerTest(mol, 2, 3); OB_ASSERT( result ); if (!result) failed.push_back(38); }

/*
  cout << "Razinger paper, fig. 7: structure 39" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_39.mol")) );
  { bool result = doStereoisomerTest(mol, 3, 1); OB_ASSERT( result ); if (!result) failed.push_back(39); }
  */  
  /*
  cout << "Razinger paper, fig. 7: structure 40" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_40.mol")) );
  { bool result = doStereoisomerTest(mol, 2, 3) );
*/
  cout << "Razinger paper, fig. 7: structure 41" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_41.mol")) );
  { bool result = doStereoisomerTest(mol, 1, 1); OB_ASSERT( result ); if (!result) failed.push_back(41); }

  cout << "Razinger paper, fig. 7: structure 42" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_42.mol")) );
  { bool result = doStereoisomerTest(mol, 2, 3); OB_ASSERT( result ); if (!result) failed.push_back(42); }

  cout << "Razinger paper, fig. 7: structure 43" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_43.mol")) );
  { bool result = doStereoisomerTest(mol, 2, 0); OB_ASSERT( result ); if (!result) failed.push_back(43); }

  cout << "Razinger paper, fig. 7: structure 44" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_44.mol")) );
  { bool result = doStereoisomerTest(mol, 1, 1); OB_ASSERT( result ); if (!result) failed.push_back(44); }

  cout << "Razinger paper, fig. 7: structure 45" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_45.mol")) );
  { bool result = doStereoisomerTest(mol, 1, 0); OB_ASSERT( result ); if (!result) failed.push_back(45); }

  cout << "Razinger paper, fig. 7: structure 46" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_46.mol")) );
  { bool result = doStereoisomerTest(mol, 0, 0); OB_ASSERT( result ); if (!result) failed.push_back(46); }
/*
  cout << "Razinger paper, fig. 7: structure 47" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_47.mol")) );
  { bool result = doStereoisomerTest(mol, 3, 0) );
  */
  
  cout << "Razinger paper, fig. 7: structure 48" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_48.mol")) );
  { bool result = doStereoisomerTest(mol, 0, 2); OB_ASSERT( result ); if (!result) failed.push_back(48); }

  cout << "Razinger paper, fig. 7: structure 49" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_49.mol")) );
  { bool result = doStereoisomerTest(mol, 4, 0); OB_ASSERT( result ); if (!result) failed.push_back(49); }

  cout << "Razinger paper, fig. 7: structure 50" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_50.mol")) );
  { bool result = doStereoisomerTest(mol, 1, 1); OB_ASSERT( result ); if (!result) failed.push_back(50); }

  cout << "Razinger paper, fig. 7: structure 51" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_51.mol")) );
  { bool result = doStereoisomerTest(mol, 4, 2); OB_ASSERT( result ); if (!result) failed.push_back(51); }

  cout << "Razinger paper, fig. 7: structure 52" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_52.mol")) );
  { bool result = doStereoisomerTest(mol, 2, 3); OB_ASSERT( result ); if (!result) failed.push_back(52); }
  
  /*
  cout << "Razinger paper, fig. 7: structure 53" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_53.mol")) );
  { bool result = doStereoisomerTest(mol, 0, 2) );

  cout << "Razinger paper, fig. 7: structure 54" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_54.mol")) );
  { bool result = doStereoisomerTest(mol, 0, 2) );

  cout << "Razinger paper, fig. 7: structure 55" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_55.mol")) );
  { bool result = doStereoisomerTest(mol, 4, 0) );

  cout << "Razinger paper, fig. 7: structure 56" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_56.mol")) );
  { bool result = doStereoisomerTest(mol, 8, 0) );
  */
  
  cout << "Razinger paper, fig. 7: structure 57" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_57.mol")) );
  { bool result = doStereoisomerTest(mol, 4, 2); OB_ASSERT( result ); if (!result) failed.push_back(57); }
  
  cout << "Razinger paper, fig. 7: structure 58" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_58.mol")) );
  { bool result = doStereoisomerTest(mol, 1, 0); OB_ASSERT( result ); if (!result) failed.push_back(58); }

  cout << "Razinger paper, fig. 7: structure 59" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_59.mol")) );
  { bool result = doStereoisomerTest(mol, 3, 1); OB_ASSERT( result ); if (!result) failed.push_back(59); }

  cout << "Razinger paper, fig. 7: structure 60" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_60.mol")) );
  { bool result = doStereoisomerTest(mol, 10, 0); OB_ASSERT( result ); if (!result) failed.push_back(60); }

  cout << "Razinger paper, fig. 7: structure 61" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_61.mol")) );
  { bool result = doStereoisomerTest(mol, 128, 0); OB_ASSERT( result ); if (!result) failed.push_back(61); }

  cout << "Razinger paper, fig. 7: structure 62" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_62.mol")) );
  { bool result = doStereoisomerTest(mol, 1, 0); OB_ASSERT( result ); if (!result) failed.push_back(62); }

  cout << "Razinger paper, fig. 7: structure 63" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_63.mol")) );
  { bool result = doStereoisomerTest(mol, 1, 0); OB_ASSERT( result ); if (!result) failed.push_back(63); }

  cout << "Razinger paper, fig. 7: structure 64" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_64.mol")) );
  { bool result = doStereoisomerTest(mol, 1, 1); OB_ASSERT( result ); if (!result) failed.push_back(64); }
  
  cout << "Razinger paper, fig. 7: structure 65" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_65.mol")) );
  { bool result = doStereoisomerTest(mol, 0, 3); OB_ASSERT( result ); if (!result) failed.push_back(65); }

  cout << "Razinger paper, fig. 7: structure 66" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_66.mol")) );
  { bool result = doStereoisomerTest(mol, 0, 14); OB_ASSERT( result ); if (!result) failed.push_back(66); }

  cout << "Razinger paper, fig. 7: structure 67" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_67.mol")) );
  { bool result = doStereoisomerTest(mol, 2, 3); OB_ASSERT( result ); if (!result) failed.push_back(67); }

  cout << "Razinger paper, fig. 7: structure 68" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_68.mol")) );
  { bool result = doStereoisomerTest(mol, 0, 3); OB_ASSERT( result ); if (!result) failed.push_back(68); }
  
  cout << "Razinger paper, fig. 7: structure 69" << endl;
  OB_ASSERT( conv.ReadFile(&mol, GetFilename("stereo/razinger_fig7_69.mol")) );
  { bool result = doStereoisomerTest(mol, 1, 7); OB_ASSERT( result ); if (!result) failed.push_back(69); }
  
  if (failed.size()) {
    cout << "FAILED TESTS: ";
    for (unsigned int i = 0; i < failed.size(); ++i)
      cout << failed.at(i) << " ";
    cout << " [" << failed.size() << " total]" << endl;
  }
  
  cout << "PASSED TESTS: " << testCount - failed.size() << "/" << testCount << endl;

//  conv.ReadString(&mol, "CC1CCC(C)CC1");
}

void testParityMatrix()
{
  Eigen::MatrixXi parityMatrix = stereoParityMatrix(4);
  OB_REQUIRE( parityMatrix.rows() == 16 );
  OB_REQUIRE( parityMatrix.cols() == 4 );

  // row 1
  OB_REQUIRE( parityMatrix(0,0) == 1 );
  OB_REQUIRE( parityMatrix(0,1) == 1 );
  OB_REQUIRE( parityMatrix(0,2) == 1 );
  OB_REQUIRE( parityMatrix(0,3) == 1 );
  // 2 
  OB_REQUIRE( parityMatrix(1,0) == 1 );
  OB_REQUIRE( parityMatrix(1,1) == 1 );
  OB_REQUIRE( parityMatrix(1,2) == 1 );
  OB_REQUIRE( parityMatrix(1,3) == -1 );
  // 3
  OB_REQUIRE( parityMatrix(2,0) == 1 );
  OB_REQUIRE( parityMatrix(2,1) == 1 );
  OB_REQUIRE( parityMatrix(2,2) == -1 );
  OB_REQUIRE( parityMatrix(2,3) == 1 );
  // 4
  OB_REQUIRE( parityMatrix(3,0) == 1 );
  OB_REQUIRE( parityMatrix(3,1) == 1 );
  OB_REQUIRE( parityMatrix(3,2) == -1 );
  OB_REQUIRE( parityMatrix(3,3) == -1 );
  // 5
  OB_REQUIRE( parityMatrix(4,0) == 1 );
  OB_REQUIRE( parityMatrix(4,1) == -1 );
  OB_REQUIRE( parityMatrix(4,2) == 1 );
  OB_REQUIRE( parityMatrix(4,3) == 1 );
  // 6
  OB_REQUIRE( parityMatrix(5,0) == 1 );
  OB_REQUIRE( parityMatrix(5,1) == -1 );
  OB_REQUIRE( parityMatrix(5,2) == 1 );
  OB_REQUIRE( parityMatrix(5,3) == -1 );
  // 7
  OB_REQUIRE( parityMatrix(6,0) == 1 );
  OB_REQUIRE( parityMatrix(6,1) == -1 );
  OB_REQUIRE( parityMatrix(6,2) == -1 );
  OB_REQUIRE( parityMatrix(6,3) == 1 );
  // 8
  OB_REQUIRE( parityMatrix(7,0) == 1 );
  OB_REQUIRE( parityMatrix(7,1) == -1 );
  OB_REQUIRE( parityMatrix(7,2) == -1 );
  OB_REQUIRE( parityMatrix(7,3) == -1 );
  // 9
  OB_REQUIRE( parityMatrix(8,0) == -1 );
  OB_REQUIRE( parityMatrix(8,1) == 1 );
  OB_REQUIRE( parityMatrix(8,2) == 1 );
  OB_REQUIRE( parityMatrix(8,3) == 1 );
  // 10
  OB_REQUIRE( parityMatrix(9,0) == -1 );
  OB_REQUIRE( parityMatrix(9,1) == 1 );
  OB_REQUIRE( parityMatrix(9,2) == 1 );
  OB_REQUIRE( parityMatrix(9,3) == -1 );
  // 11
  OB_REQUIRE( parityMatrix(10,0) == -1 );
  OB_REQUIRE( parityMatrix(10,1) == 1 );
  OB_REQUIRE( parityMatrix(10,2) == -1 );
  OB_REQUIRE( parityMatrix(10,3) == 1 );
  // 12
  OB_REQUIRE( parityMatrix(11,0) == -1 );
  OB_REQUIRE( parityMatrix(11,1) == 1 );
  OB_REQUIRE( parityMatrix(11,2) == -1 );
  OB_REQUIRE( parityMatrix(11,3) == -1 );
  // 13
  OB_REQUIRE( parityMatrix(12,0) == -1 );
  OB_REQUIRE( parityMatrix(12,1) == -1 );
  OB_REQUIRE( parityMatrix(12,2) == 1 );
  OB_REQUIRE( parityMatrix(12,3) == 1 );
  // 14
  OB_REQUIRE( parityMatrix(13,0) == -1 );
  OB_REQUIRE( parityMatrix(13,1) == -1 );
  OB_REQUIRE( parityMatrix(13,2) == 1 );
  OB_REQUIRE( parityMatrix(13,3) == -1 );
  // 15
  OB_REQUIRE( parityMatrix(14,0) == -1 );
  OB_REQUIRE( parityMatrix(14,1) == -1 );
  OB_REQUIRE( parityMatrix(14,2) == -1 );
  OB_REQUIRE( parityMatrix(14,3) == 1 );
  // 16
  OB_REQUIRE( parityMatrix(15,0) == -1 );
  OB_REQUIRE( parityMatrix(15,1) == -1 );
  OB_REQUIRE( parityMatrix(15,2) == -1 );
  OB_REQUIRE( parityMatrix(15,3) == -1 );
}

void testContractPermutations()
{
  // stereoAtoms = 1 2 5
  std::vector<unsigned int> stereoAtoms(3, 0);
  stereoAtoms[0] = 1;
  stereoAtoms[1] = 2;
  stereoAtoms[2] = 5;
  // create G:  1 2 3 4 5 6
  //            2 1 3 4 5 6
  //            2 3 1 4 5 6
  PermutationGroup G;
  std::vector<unsigned int> map(6, 0);
  map[0] = 1; map[1] = 2; map[2] = 3; map[3] = 4; map[4] = 5; map[5] = 6;
  G.add(Permutation(map));
  map[0] = 2; map[1] = 1; map[2] = 3; map[3] = 4; map[4] = 5; map[5] = 6;
  G.add(Permutation(map));
  map[0] = 2; map[1] = 3; map[2] = 1; map[3] = 4; map[4] = 5; map[5] = 6;
  G.add(Permutation(map));

  PermutationGroup Gc; // contracted automorphism permutation
  contractPermutations(G, Gc, stereoAtoms);
 
  OB_REQUIRE( Gc.size() == 3 );
  // 1 2 3 4 5 6 -> 1 2 5
  OB_REQUIRE( Gc.at(0).map.size() == 3 );
  OB_REQUIRE( Gc.at(0).map[0] == 1 );
  OB_REQUIRE( Gc.at(0).map[1] == 2 );
  OB_REQUIRE( Gc.at(0).map[2] == 5 );
  // 2 1 3 4 5 6 -> 2 1 5
  OB_REQUIRE( Gc.at(1).map.size() == 3 );
  OB_REQUIRE( Gc.at(1).map[0] == 2 );
  OB_REQUIRE( Gc.at(1).map[1] == 1 );
  OB_REQUIRE( Gc.at(1).map[2] == 5 );
  // 2 3 1 4 5 6 -> 2 1 5
  OB_REQUIRE( Gc.at(2).map.size() == 3 );
  OB_REQUIRE( Gc.at(2).map[0] == 2 );
  OB_REQUIRE( Gc.at(2).map[1] == 1 );
  OB_REQUIRE( Gc.at(2).map[2] == 5 );
}


void testSignedPermutationMatrices()
{
  // create Gc
  PermutationGroup Gc;
  std::vector<unsigned int> map(4);
  // 3 5 8 9
  map[0] = 3; map[1] = 5; map[2] = 8; map[3] = 9;
  Gc.add(Permutation(map));
  std::vector<unsigned int> tetrahedralAtoms = map;
  // 5 3 8 9
  map[0] = 5; map[1] = 3; map[2] = 8; map[3] = 9;
  Gc.add(Permutation(map));
  // 3 5 9 8
  map[0] = 3; map[1] = 5; map[2] = 9; map[3] = 8;
  Gc.add(Permutation(map));
  // 3 8 9 5
  map[0] = 3; map[1] = 8; map[2] = 9; map[3] = 5;
  Gc.add(Permutation(map));
  // 5 8 9 3
  map[0] = 5; map[1] = 8; map[2] = 9; map[3] = 3;
  Gc.add(Permutation(map));
  // 5 3 9 8
  map[0] = 5; map[1] = 3; map[2] = 9; map[3] = 8;
  Gc.add(Permutation(map));
  // 9 5 8 3
  map[0] = 9; map[1] = 5; map[2] = 8; map[3] = 3;
  Gc.add(Permutation(map));

  // create stereoIndexVectors
  std::vector<Eigen::VectorXi> stereoIndexVectors;

  std::vector<Eigen::MatrixXi> signedMatrices;
  signedPermutationMatrices(Gc, stereoIndexVectors, signedMatrices, tetrahedralAtoms, 0, 4);
}

void testMultiplyMatrixByRow()
{
  // 0 1 0 0  0 0
  // 1 0 0 0  0 0
  // 0 0 1 0  0 0
  // 0 0 0 1  0 0
  // 0 0 0 0 -1 0
  // 0 0 0 0  0 1
  Eigen::MatrixXi m = Eigen::MatrixXi::Zero(6, 6);
  m(0,1) = 1;
  m(1,0) = 1;
  m(2,2) = 1;
  m(3,3) = 1;
  m(4,4) = -1;
  m(5,5) = 1;

  Eigen::VectorXi row1 = Eigen::VectorXi::Zero(6);
  row1[0] = 1;
  row1[1] = 1;
  row1[2] = 1;
  row1[3] = 1;
  row1[4] = 1;
  row1[5] = 1;

  Eigen::VectorXi result = multiplyMatrixByRow(m, row1);

  OB_REQUIRE( result[0] == 1 );
  OB_REQUIRE( result[1] == 1 );
  OB_REQUIRE( result[2] == 1 );
  OB_REQUIRE( result[3] == 1 );
  OB_REQUIRE( result[4] == -1 );
  OB_REQUIRE( result[5] == 1 );

}

OBMol* readFig3()
{
  OBMol *mol = new OBMol;
  OBConversion conv;

  std::string file = GetFilename("stereo/razinger_fig3.mol");
  cout << file << endl;
  OBFormat *format = conv.FormatFromExt(file.c_str());
  OB_REQUIRE( format );
  OB_ASSERT( conv.SetInFormat(format) );
  OB_ASSERT( conv.ReadFile(mol, file) );
  return mol;
}

void testContractPermutationsForFig3()
{
  OBMol *mol = readFig3();

  OBGraphSym graphSym(mol);
  std::vector<unsigned int> symmetry_classes;
  graphSym.GetSymmetry(symmetry_classes, false);
  for (unsigned int i = 0; i < symmetry_classes.size(); ++i)
    cout << i+1 << ": " << symmetry_classes[i] << endl;

  PermutationGroup G = findAutomorphisms(mol, symmetry_classes);
  cout << "G.size = " << G.size() << endl;
  vector<StereogenicUnit> stereoUnits = FindStereogenicUnits(mol, symmetry_classes); // FIXME need to cache this

  OB_REQUIRE( stereoUnits.size() == 6 );

  // make a list of stereo atom indexes
  unsigned int numTetrahedral = 0;
  std::vector<unsigned int> stereoAtoms;
  for (vector<StereogenicUnit>::iterator unit = stereoUnits.begin(); unit != stereoUnits.end(); ++unit) {
    if (unit->type == OBStereo::Tetrahedral) {
      OBAtom *atom = mol->GetAtomById(unit->id);
      stereoAtoms.push_back(atom->GetIndex()+1);
      numTetrahedral++;
    }
  }
  unsigned int n = numTetrahedral;
  OB_REQUIRE( stereoAtoms.size() == 6 );
      
  PermutationGroup Gc; // contracted automorphism permutations
  contractPermutations(G, Gc, stereoAtoms);

  OB_REQUIRE( Gc.size() == 8 );

  // create Gc_ref
  PermutationGroup Gc_ref;
  std::vector<unsigned int> map(6);
  map[0] = 9; map[1] = 10; map[2] = 11; map[3] = 12; map[4] = 13; map[5] = 14;
  Gc_ref.add(Permutation(map)); // 1
  std::vector<unsigned int> tetrahedralAtoms = map;
  map[0] = 10; map[1] = 9; map[2] = 11; map[3] = 12; map[4] = 13; map[5] = 14;
  Gc_ref.add(Permutation(map)); // 2
  map[0] = 10; map[1] = 9; map[2] = 12; map[3] = 11; map[4] = 13; map[5] = 14;
  Gc_ref.add(Permutation(map)); // 3
  map[0] = 9; map[1] = 10; map[2] = 12; map[3] = 11; map[4] = 13; map[5] = 14;
  Gc_ref.add(Permutation(map)); // 4
  map[0] = 11; map[1] = 12; map[2] = 10; map[3] = 9; map[4] = 14; map[5] = 13;
  Gc_ref.add(Permutation(map)); // 5
  map[0] = 12; map[1] = 11; map[2] = 9; map[3] = 10; map[4] = 14; map[5] = 13;
  Gc_ref.add(Permutation(map)); // 6
  map[0] = 12; map[1] = 11; map[2] = 10; map[3] = 9; map[4] = 14; map[5] = 13;
  Gc_ref.add(Permutation(map)); // 7
  map[0] = 11; map[1] = 12; map[2] = 9; map[3] = 10; map[4] = 14; map[5] = 13;
  Gc_ref.add(Permutation(map)); // 8

  //9 10 11 12 13 14    1
  //9 10 12 11 13 14    4
  //10 9 11 12 13 14    2
  //11 12 9 10 14 13    8
  //10 9 12 11 13 14    3
  //11 12 10 9 14 13    5
  //12 11 9 10 14 13    6
  //12 11 10 9 14 13    7


  for (unsigned int i = 0; i < Gc_ref.size(); ++i) {
      cout << i << endl;
    OB_REQUIRE( Gc.contains(Gc_ref.at(i)) );
    if (!Gc.contains(Gc_ref.at(i))) {
    }
  }

  std::vector<Eigen::VectorXi> stereoIndexVectors;
  createStereoIndexVectors(stereoIndexVectors, Gc, G, stereoUnits, mol, symmetry_classes, n);

  OB_REQUIRE( stereoIndexVectors.size() == 8 );

  OB_REQUIRE( stereoIndexVectors.at(0)[0] == 1 );
  OB_REQUIRE( stereoIndexVectors.at(0)[1] == 1 );
  OB_REQUIRE( stereoIndexVectors.at(0)[2] == 1 );
  OB_REQUIRE( stereoIndexVectors.at(0)[3] == 1 );
  OB_REQUIRE( stereoIndexVectors.at(0)[4] == 1 );
  OB_REQUIRE( stereoIndexVectors.at(0)[5] == 1 );

  OB_REQUIRE( stereoIndexVectors.at(1)[0] == 1 );
  OB_REQUIRE( stereoIndexVectors.at(1)[1] == 1 );
  OB_REQUIRE( stereoIndexVectors.at(1)[2] == 1 );
  OB_REQUIRE( stereoIndexVectors.at(1)[3] == 1 );
  OB_REQUIRE( stereoIndexVectors.at(1)[4] == 1 );
  OB_REQUIRE( stereoIndexVectors.at(1)[5] == -1 );

  OB_REQUIRE( stereoIndexVectors.at(2)[0] == 1 );
  OB_REQUIRE( stereoIndexVectors.at(2)[1] == 1 );
  OB_REQUIRE( stereoIndexVectors.at(2)[2] == 1 );
  OB_REQUIRE( stereoIndexVectors.at(2)[3] == 1 );
  OB_REQUIRE( stereoIndexVectors.at(2)[4] == -1 );
  OB_REQUIRE( stereoIndexVectors.at(2)[5] == 1 );

  OB_REQUIRE( stereoIndexVectors.at(3)[0] == 1 );
  OB_REQUIRE( stereoIndexVectors.at(3)[1] == 1 );
  OB_REQUIRE( stereoIndexVectors.at(3)[2] == 1 );
  OB_REQUIRE( stereoIndexVectors.at(3)[3] == 1 );
  OB_REQUIRE( stereoIndexVectors.at(3)[4] == 1 );
  OB_REQUIRE( stereoIndexVectors.at(3)[5] == 1 );

  OB_REQUIRE( stereoIndexVectors.at(4)[0] == 1 );
  OB_REQUIRE( stereoIndexVectors.at(4)[1] == 1 );
  OB_REQUIRE( stereoIndexVectors.at(4)[2] == 1 );
  OB_REQUIRE( stereoIndexVectors.at(4)[3] == 1 );
  OB_REQUIRE( stereoIndexVectors.at(4)[4] == -1 );
  OB_REQUIRE( stereoIndexVectors.at(4)[5] == -1 );

  OB_REQUIRE( stereoIndexVectors.at(5)[0] == 1 );
  OB_REQUIRE( stereoIndexVectors.at(5)[1] == 1 );
  OB_REQUIRE( stereoIndexVectors.at(5)[2] == 1 );
  OB_REQUIRE( stereoIndexVectors.at(5)[3] == 1 );
  OB_REQUIRE( stereoIndexVectors.at(5)[4] == -1 );
  OB_REQUIRE( stereoIndexVectors.at(5)[5] == 1 );

  OB_REQUIRE( stereoIndexVectors.at(6)[0] == 1 );
  OB_REQUIRE( stereoIndexVectors.at(6)[1] == 1 );
  OB_REQUIRE( stereoIndexVectors.at(6)[2] == 1 );
  OB_REQUIRE( stereoIndexVectors.at(6)[3] == 1 );
  OB_REQUIRE( stereoIndexVectors.at(6)[4] == 1 );
  OB_REQUIRE( stereoIndexVectors.at(6)[5] == -1 );

  OB_REQUIRE( stereoIndexVectors.at(7)[0] == 1 );
  OB_REQUIRE( stereoIndexVectors.at(7)[1] == 1 );
  OB_REQUIRE( stereoIndexVectors.at(7)[2] == 1 );
  OB_REQUIRE( stereoIndexVectors.at(7)[3] == 1 );
  OB_REQUIRE( stereoIndexVectors.at(7)[4] == -1 );
  OB_REQUIRE( stereoIndexVectors.at(7)[5] == -1 );


  std::vector<Eigen::MatrixXi> signedMatrices;
  signedPermutationMatrices(Gc, stereoIndexVectors, signedMatrices, tetrahedralAtoms, 0, n);

  /*
  OB_REQUIRE( signedMatrices.at(0)(0,0) == 1 );
  OB_REQUIRE( signedMatrices.at(0)(1,1) == 1 );
  OB_REQUIRE( signedMatrices.at(0)(2,2) == 1 );
  OB_REQUIRE( signedMatrices.at(0)(3,3) == 1 );
  OB_REQUIRE( signedMatrices.at(0)(4,4) == 1 );
  OB_REQUIRE( signedMatrices.at(0)(5,5) == 1 );

  OB_REQUIRE( signedMatrices.at(1)(0,0) == 1 );
  OB_REQUIRE( signedMatrices.at(1)(1,1) == 1 );
  OB_REQUIRE( signedMatrices.at(1)(2,3) == 1 );
  OB_REQUIRE( signedMatrices.at(1)(3,2) == 1 );
  OB_REQUIRE( signedMatrices.at(1)(4,4) == 1 );
  OB_REQUIRE( signedMatrices.at(1)(5,5) == -1 );

  OB_REQUIRE( signedMatrices.at(2)(0,1) == 1 );
  OB_REQUIRE( signedMatrices.at(2)(1,0) == 1 );
  OB_REQUIRE( signedMatrices.at(2)(2,2) == 1 );
  OB_REQUIRE( signedMatrices.at(2)(3,3) == 1 );
  OB_REQUIRE( signedMatrices.at(2)(4,4) == -1 );
  OB_REQUIRE( signedMatrices.at(2)(5,5) == 1 );

  OB_REQUIRE( signedMatrices.at(3)(0,2) == 1 );
  OB_REQUIRE( signedMatrices.at(3)(1,3) == 1 );
  OB_REQUIRE( signedMatrices.at(3)(2,0) == 1 );
  OB_REQUIRE( signedMatrices.at(3)(3,1) == 1 );
  OB_REQUIRE( signedMatrices.at(3)(4,5) == 1 );
  OB_REQUIRE( signedMatrices.at(3)(5,4) == 1 );

  OB_REQUIRE( signedMatrices.at(4)(0,1) == 1 );
  OB_REQUIRE( signedMatrices.at(4)(1,0) == 1 );
  OB_REQUIRE( signedMatrices.at(4)(2,3) == 1 );
  OB_REQUIRE( signedMatrices.at(4)(3,2) == 1 );
  OB_REQUIRE( signedMatrices.at(4)(4,4) == -1 );
  OB_REQUIRE( signedMatrices.at(4)(5,5) == -1 );

  OB_REQUIRE( signedMatrices.at(5)(0,3) == 1 );
  OB_REQUIRE( signedMatrices.at(5)(1,2) == 1 );
  OB_REQUIRE( signedMatrices.at(5)(2,0) == 1 );
  OB_REQUIRE( signedMatrices.at(5)(3,1) == 1 );
  OB_REQUIRE( signedMatrices.at(5)(4,5) == -1 );
  OB_REQUIRE( signedMatrices.at(5)(5,4) == 1 );

  OB_REQUIRE( signedMatrices.at(6)(0,2) == 1 );
  OB_REQUIRE( signedMatrices.at(6)(1,3) == 1 );
  OB_REQUIRE( signedMatrices.at(6)(2,1) == 1 );
  OB_REQUIRE( signedMatrices.at(6)(3,0) == 1 );
  OB_REQUIRE( signedMatrices.at(6)(4,5) == 1 );
  OB_REQUIRE( signedMatrices.at(6)(5,4) == -1 );

  OB_REQUIRE( signedMatrices.at(7)(0,3) == 1 );
  OB_REQUIRE( signedMatrices.at(7)(1,2) == 1 );
  OB_REQUIRE( signedMatrices.at(7)(2,1) == 1 );
  OB_REQUIRE( signedMatrices.at(7)(3,0) == 1 );
  OB_REQUIRE( signedMatrices.at(7)(4,5) == -1 );
  OB_REQUIRE( signedMatrices.at(7)(5,4) == -1 );
  */








  delete mol;
}


int main(int argc, char **argv)
{
  // some static tests first
  testParityMatrix();
  testContractPermutations();
  //testSignedPermutationMatrices();
  testMultiplyMatrixByRow();

  // fig3
  testContractPermutationsForFig3();

  std::string file;
  if (argc > 2)
    file = argv[2];

  if (argc > 1)
    test_Stereoisomer(atoi(argv[1]), file);
  else
    test_Stereoisomer(0, file);
  
  return 0;
}
