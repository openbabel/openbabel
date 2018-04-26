"""Test mol2 fileformat convarsion

On Windows or Linux, you can run these tests at the commandline
in the build folder with:
"C:\Program Files\CMake 2.6\bin\ctest.exe" -C CTestTestfile.cmake
                                           -R pybindtest -VV

You could also "chdir" into build/test and run the test file directly:
python ../../test/testexample.py

In both cases, the test file is run directly from the source folder,
and so you can quickly develop the tests and try them out.
"""

from testbindings import *

class Mol2Test(PybelWrapper):
    def testResNumAppend(self):
        mol_string = """@<TRIPOS>MOLECULE
*****
 54 53 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
      1  N         -2.5810    2.2660    0.2260 N.3     1  LEU1       -0.1231
      2  CA        -2.8680    1.7080   -1.1400 C.3     1  LEU1        0.1693
      3  C         -3.9170    0.6460   -1.2420 C.2     1  LEU1        0.2592
      4  O         -3.7340   -0.2850   -2.0340 O.2     1  LEU1       -0.2718
      5  CB        -3.4130    2.5810   -2.2680 C.3     1  LEU1        0.0284
      6  CG        -2.5740    3.8140   -2.2850 C.3     1  LEU1        0.0019
      7  CD1       -3.3880    5.0870   -2.2100 C.3     1  LEU1        0.0001
      8  CD2       -1.3500    3.5120   -3.1300 C.3     1  LEU1        0.0001
      9  N         -5.2070    0.8630   -0.6750 N.am    2  MET2       -0.1965
     10  CA        -5.3060    0.6220    0.7720 C.3     2  MET2        0.1736
     11  C         -4.5610   -0.6800    0.7200 C.2     2  MET2        0.2616
     12  O         -4.9780   -1.5400   -0.0660 O.2     2  MET2       -0.2715
     13  CB        -6.6150    0.0530    1.5330 C.3     2  MET2        0.0427
     14  CG        -6.1450   -0.3270    2.9860 C.3     2  MET2        0.0599
     15  SD        -6.2790   -1.7760    4.0800 S.3     2  MET2       -0.1393
     16  CE        -4.5960   -1.5820    4.8950 C.3     2  MET2        0.0692
     17  N         -3.4960   -0.8220    1.5520 N.am    3  ALA3       -0.1963
     18  CA        -2.6430   -1.9960    1.4510 C.3     3  ALA3        0.1728
     19  C         -1.4070   -1.5100    0.5560 C.2     3  ALA3        0.2616
     20  O         -1.4660   -0.4290   -0.0430 O.2     3  ALA3       -0.2715
     21  CB        -3.4150   -3.2280    0.9710 C.3     3  ALA3        0.0334
     22  N         -0.1340   -2.1100    0.6280 N.am    4  ILE4       -0.1963
     23  CA         1.0440   -1.1230    0.7360 C.3     4  ILE4        0.1727
     24  C          2.8800   -1.1710    1.0320 C.2     4  ILE4        0.2616
     25  O          3.3760   -2.1310    1.6420 O.2     4  ILE4       -0.2715
     26  CB         0.6270   -0.4850    2.0700 C.3     4  ILE4        0.0290
     27  CG1        0.8620   -1.4560    3.2320 C.3     4  ILE4        0.0022
     28  CG2       -0.4700    0.5900    2.3490 C.3     4  ILE4        0.0023
     29  CD1       -0.2030   -1.4650    4.3310 C.3     4  ILE4        0.0001
     30  N          3.6870    0.0840    0.7950 N.am    5  ASN5       -0.1958
     31  CA         4.9100    0.9270    1.6770 C.3     5  ASN5        0.1819
     32  C          4.5770    2.6080    1.8700 C.2     5  ASN5        0.2621
     33  O          4.6660    3.3300    0.8450 O.2     5  ASN5       -0.2715
     34  CB         6.3520    0.9800    1.0220 C.3     5  ASN5        0.1195
     35  CG         7.4410    1.6010    1.9300 C.2     5  ASN5        0.2630
     36  OD2        8.3940    0.9540    2.3670 O.2     5  ASN5       -0.2716
     37  ND1        7.3890    2.9480    2.1640 N.am    5  ASN5       -0.0877
     38  N          4.4020    3.3260    3.0800 N.am    6  ILE6       -0.1963
     39  CA         3.0040    3.7620    3.6070 C.3     6  ILE6        0.1727
     40  C          2.5800    4.8180    4.6970 C.2     6  ILE6        0.2617
     41  O          2.9010    5.9940    4.5740 O.2     6  ILE6       -0.2715
     42  CB         2.2890    2.5250    3.9750 C.3     6  ILE6        0.0290
     43  CG1        2.1350    1.8150    2.7110 C.3     6  ILE6        0.0022
     44  CG2        3.0680    1.6830    5.0110 C.3     6  ILE6        0.0023
     45  CD1        1.7180    2.5110    1.3810 C.3     6  ILE6        0.0001
     46  N          1.8110    4.4030    5.7930 N.am    7  ASP7       -0.1936
     47  CA         1.6900    5.1670    7.0780 C.3     7  ASP7        0.2090
     48  C          3.1120    5.4260    7.5360 C.2     7  ASP7        0.3848
     49  O          3.6920    4.8270    8.4270 O.co2   7  ASP7       -0.2437
     50  CB         0.9230    4.2470    8.1040 C.3     7  ASP7        0.1472
     51  CG        -0.2800    4.8440    8.8050 C.2     7  ASP7        0.3672
     52  OD1       -0.2680    5.2030    9.9720 O.co2   7  ASP7       -0.2456
     53  OD2       -1.4170    4.8600    8.0790 O.co2   7  ASP7       -0.2456
     54  OXT        3.7770    6.3550    6.7920 O.co2   7  ASP7       -0.2437
@<TRIPOS>BOND
     1     1     2    1
     2     2     3    1
     3     2     5    1
     4     3     4    2
     5     8     6    1
     6     6     5    1
     7     6     7    1
     8     3     9   am
     9     9    10    1
    10    10    11    1
    11    10    13    1
    12    11    12    2
    13    15    14    1
    14    15    16    1
    15    14    13    1
    16    11    17   am
    17    17    18    1
    18    18    19    1
    19    18    21    1
    20    20    19    2
    21    19    22   am
    22    22    23    1
    23    27    29    1
    24    27    26    1
    25    23    26    1
    26    23    24    1
    27    26    28    1
    28    24    25    2
    29    24    30   am
    30    30    31    1
    31    31    32    1
    32    31    34    1
    33    33    32    2
    34    34    35    1
    35    37    35   am
    36    35    36    2
    37    32    38   am
    38    38    39    1
    39    43    45    1
    40    43    42    1
    41    39    42    1
    42    39    40    1
    43    42    44    1
    44    40    41    2
    45    40    46   am
    46    46    47    1
    47    47    50    1
    48    47    48    1
    49    50    51    1
    50    52    51   ar
    51    51    53   ar
    52    48    49   ar
    53    48    54   ar
"""

        mol = pybel.readstring("mol2", mol_string)
        for i in range(5):
            for res in mol.residues:
                self.assertEqual(res.name, "%s%i" % (res.name[:3], res.number))

def gettests():
    testsuite = []
    for myclass in [Mol2Test]:
        suite = unittest.TestLoader().loadTestsFromTestCase(myclass)
        testsuite.append(suite)
    return testsuite

if __name__ == "__main__":
    unittest.TextTestRunner().run(unittest.TestSuite(gettests()))
