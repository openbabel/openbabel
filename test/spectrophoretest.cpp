#include <vector>

#include "obtest.h"
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/spectrophore.h>



// Normal Spectrophore(TM): test all values
void test01(OpenBabel::OBMol* mol)
{
   OpenBabel::OBSpectrophore s;
   s.SetNormalization(OpenBabel::OBSpectrophore::NoNormalization);
   s.SetResolution(3.0);
   s.SetAccuracy(OpenBabel::OBSpectrophore::AngStepSize20);
   s.SetStereo(OpenBabel::OBSpectrophore::NoStereoSpecificProbes);
   std::vector<double> r = s.GetSpectrophore(mol);
   OB_REQUIRE( (r[0]  >  1.598) && (r[0]  <  1.600) );
   OB_REQUIRE( (r[1]  >  1.576) && (r[1]  <  1.578) );
   OB_REQUIRE( (r[2]  >  1.169) && (r[2]  <  1.171) );
   OB_REQUIRE( (r[3]  >  3.760) && (r[3]  <  3.762) );
   OB_REQUIRE( (r[4]  >  5.174) && (r[4]  <  5.176) );
   OB_REQUIRE( (r[5]  >  5.780) && (r[5]  <  5.782) );
   OB_REQUIRE( (r[6]  >  3.796) && (r[6]  <  3.798) );
   OB_REQUIRE( (r[7]  >  3.712) && (r[7]  <  3.714) );
   OB_REQUIRE( (r[8]  >  4.650) && (r[8]  <  4.652) );
   OB_REQUIRE( (r[9]  >  7.736) && (r[9]  <  7.738) );
   OB_REQUIRE( (r[10] >  7.949) && (r[10] <  7.951) );
   OB_REQUIRE( (r[11] >  4.868) && (r[11] <  4.870) );
   OB_REQUIRE( (r[12] >  2.707) && (r[12] <  2.709) );
   OB_REQUIRE( (r[13] >  3.470) && (r[13] <  3.472) );
   OB_REQUIRE( (r[14] >  6.698) && (r[14] <  6.700) );
   OB_REQUIRE( (r[15] >  9.485) && (r[15] <  9.487) );
   OB_REQUIRE( (r[16] >  7.667) && (r[16] <  7.669) );
   OB_REQUIRE( (r[17] >  8.881) && (r[17] <  8.883) );
   OB_REQUIRE( (r[18] >  4.899) && (r[18] <  4.901) );
   OB_REQUIRE( (r[19] >  7.478) && (r[19] <  7.480) );
   OB_REQUIRE( (r[20] >  9.323) && (r[20] <  9.325) );
   OB_REQUIRE( (r[21] > 10.292) && (r[21] < 10.294) );
   OB_REQUIRE( (r[22] > 12.955) && (r[22] < 12.957) );
   OB_REQUIRE( (r[23] > 10.334) && (r[23] < 10.336) );
   OB_REQUIRE( (r[24] >  4.020) && (r[24] <  4.022) );
   OB_REQUIRE( (r[25] >  3.813) && (r[25] <  3.815) );
   OB_REQUIRE( (r[26] >  2.946) && (r[26] <  2.948) );
   OB_REQUIRE( (r[27] >  6.380) && (r[27] <  6.382) );
   OB_REQUIRE( (r[28] > 11.003) && (r[28] < 11.005) );
   OB_REQUIRE( (r[29] >  8.278) && (r[29] <  8.280) );
   OB_REQUIRE( (r[30] >  6.548) && (r[30] <  6.550) );
   OB_REQUIRE( (r[31] >  7.135) && (r[31] <  7.137) );
   OB_REQUIRE( (r[32] >  8.612) && (r[32] <  8.614) );
   OB_REQUIRE( (r[33] > 13.181) && (r[33] < 13.183) );
   OB_REQUIRE( (r[34] > 13.743) && (r[34] < 13.745) );
   OB_REQUIRE( (r[35] >  9.083) && (r[35] <  9.085) );
   OB_REQUIRE( (r[36] >  0.458) && (r[36] <  0.460) );
   OB_REQUIRE( (r[37] >  0.641) && (r[37] <  0.643) );
   OB_REQUIRE( (r[38] >  2.171) && (r[38] <  2.173) );
   OB_REQUIRE( (r[39] >  2.752) && (r[39] <  2.754) );
   OB_REQUIRE( (r[40] >  2.347) && (r[40] <  2.349) );
   OB_REQUIRE( (r[41] >  2.604) && (r[41] <  2.606) );
   OB_REQUIRE( (r[42] >  1.613) && (r[42] <  1.615) );
   OB_REQUIRE( (r[43] >  3.165) && (r[43] <  3.167) );
   OB_REQUIRE( (r[44] >  3.390) && (r[44] <  3.392) );
   OB_REQUIRE( (r[45] >  3.131) && (r[45] <  3.133) );
   OB_REQUIRE( (r[46] >  4.104) && (r[46] <  4.106) );
   OB_REQUIRE( (r[47] >  2.874) && (r[47] <  2.876) );
}



// With increased accuracy
void test02(OpenBabel::OBMol* mol)
{
   OpenBabel::OBSpectrophore s;
   s.SetNormalization(OpenBabel::OBSpectrophore::NoNormalization);
   s.SetResolution(3.0);
   s.SetAccuracy(OpenBabel::OBSpectrophore::AngStepSize5);
   s.SetStereo(OpenBabel::OBSpectrophore::NoStereoSpecificProbes);
   std::vector<double> r = s.GetSpectrophore(mol);
   OB_REQUIRE( (r[0]  > 1.644) && (r[0]  < 1.645) );
   OB_REQUIRE( (r[12] > 2.724) && (r[12] < 2.725) );
   OB_REQUIRE( (r[24] > 4.043) && (r[24] < 4.044) );
   OB_REQUIRE( (r[36] > 0.458) && (r[36] < 0.459) );
}



// With modified resolution
void test03(OpenBabel::OBMol* mol)
{
   OpenBabel::OBSpectrophore s;
   s.SetNormalization(OpenBabel::OBSpectrophore::NoNormalization);
   s.SetResolution(10.0);
   s.SetAccuracy(OpenBabel::OBSpectrophore::AngStepSize20);
   s.SetStereo(OpenBabel::OBSpectrophore::NoStereoSpecificProbes);
   std::vector<double> r = s.GetSpectrophore(mol);
   OB_REQUIRE( (r[1]  > 0.146) && (r[1]  < 0.147) );
   OB_REQUIRE( (r[13] > 0.327) && (r[13] < 0.328) );
   OB_REQUIRE( (r[25] > 0.359) && (r[25] < 0.360) );
   OB_REQUIRE( (r[37] > 0.061) && (r[37] < 0.062) );
}



// With normalization towards zero mean
void test04(OpenBabel::OBMol* mol)
{
   OpenBabel::OBSpectrophore s;
   s.SetNormalization(OpenBabel::OBSpectrophore::NormalizationTowardsZeroMean);
   s.SetResolution(3.0);
   s.SetAccuracy(OpenBabel::OBSpectrophore::AngStepSize20);
   s.SetStereo(OpenBabel::OBSpectrophore::NoStereoSpecificProbes);
   std::vector<double> r = s.GetSpectrophore(mol);
   OB_REQUIRE( (r[2]  > -3.145) && (r[2]  < -3.144) );
   OB_REQUIRE( (r[14] > -1.152) && (r[14] < -1.151) );
   OB_REQUIRE( (r[26] > -4.949) && (r[26] < -4.948) );
   OB_REQUIRE( (r[38] > -0.267) && (r[38] < -0.266) );
}



// With normalization towards unit std
void test05(OpenBabel::OBMol* mol)
{
   OpenBabel::OBSpectrophore s;
   s.SetNormalization(OpenBabel::OBSpectrophore::NormalizationTowardsUnitStd);
   s.SetResolution(3.0);
   s.SetAccuracy(OpenBabel::OBSpectrophore::AngStepSize20);
   s.SetStereo(OpenBabel::OBSpectrophore::NoStereoSpecificProbes);
   std::vector<double> r = s.GetSpectrophore(mol);
   OB_REQUIRE( (r[3]  > 1.698) && (r[3]  < 1.699) );
   OB_REQUIRE( (r[15] > 3.147) && (r[15] < 3.148) );
   OB_REQUIRE( (r[27] > 1.823) && (r[27] < 1.824) );
   OB_REQUIRE( (r[39] > 2.539) && (r[39] < 2.540) );
}



// With normalization towards unit std and zero mean
void test06(OpenBabel::OBMol* mol)
{
   OpenBabel::OBSpectrophore s;
   s.SetNormalization(OpenBabel::OBSpectrophore::NormalizationTowardsZeroMeanAndUnitStd);
   s.SetResolution(3.0);
   s.SetAccuracy(OpenBabel::OBSpectrophore::AngStepSize20);
   s.SetStereo(OpenBabel::OBSpectrophore::NoStereoSpecificProbes);
   std::vector<double> r = s.GetSpectrophore(mol);
   OB_REQUIRE( (r[4]  >  0.388) && (r[4]  <  0.389) );
   OB_REQUIRE( (r[16] > -0.061) && (r[16] < -0.060) );
   OB_REQUIRE( (r[28] >  0.887) && (r[28] <  0.888) );
   OB_REQUIRE( (r[40] > -0.084) && (r[40] < -0.083) );
}



// Unique stereo probes
void test07(OpenBabel::OBMol* mol)
{
   OpenBabel::OBSpectrophore s;
   s.SetNormalization(OpenBabel::OBSpectrophore::NoNormalization);
   s.SetResolution(3.0);
   s.SetAccuracy(OpenBabel::OBSpectrophore::AngStepSize20);
   s.SetStereo(OpenBabel::OBSpectrophore::UniqueStereoSpecificProbes);
   std::vector<double> r = s.GetSpectrophore(mol);
   OB_REQUIRE( (r[0]  >  1.150) && (r[0]  <  1.151) );
   OB_REQUIRE( (r[1]  >  3.149) && (r[1]  <  3.150) );
   OB_REQUIRE( (r[2]  >  3.209) && (r[2]  <  3.210) );
   OB_REQUIRE( (r[3]  >  3.174) && (r[3]  <  3.175) );
   OB_REQUIRE( (r[4]  >  3.569) && (r[4]  <  3.570) );
   OB_REQUIRE( (r[5]  >  3.304) && (r[5]  <  3.305) );
   OB_REQUIRE( (r[6]  >  3.051) && (r[6]  <  3.052) );
   OB_REQUIRE( (r[7]  >  4.562) && (r[7]  <  4.563) );
   OB_REQUIRE( (r[8]  >  6.095) && (r[8]  <  6.096) );
   OB_REQUIRE( (r[9]  >  5.226) && (r[9]  <  5.227) );
   OB_REQUIRE( (r[10] >  4.518) && (r[10] <  4.519) );
   OB_REQUIRE( (r[11] >  5.049) && (r[11] <  5.050) );
   OB_REQUIRE( (r[12] >  6.967) && (r[12] <  6.968) );
   OB_REQUIRE( (r[13] >  5.714) && (r[13] <  5.715) );
   OB_REQUIRE( (r[14] >  4.741) && (r[14] <  4.742) );
   OB_REQUIRE( (r[15] >  5.948) && (r[15] <  5.949) );
   OB_REQUIRE( (r[16] >  5.890) && (r[16] <  5.891) );
   OB_REQUIRE( (r[17] >  7.335) && (r[17] <  7.336) );
   OB_REQUIRE( (r[18] >  2.015) && (r[18] <  2.016) );
   OB_REQUIRE( (r[19] >  8.571) && (r[19] <  8.572) );
   OB_REQUIRE( (r[20] >  5.878) && (r[20] <  5.879) );
   OB_REQUIRE( (r[21] >  7.161) && (r[21] <  7.162) );
   OB_REQUIRE( (r[22] >  5.824) && (r[22] <  5.825) );
   OB_REQUIRE( (r[23] >  5.744) && (r[23] <  5.745) );
   OB_REQUIRE( (r[24] >  5.676) && (r[24] <  5.677) );
   OB_REQUIRE( (r[25] >  8.814) && (r[25] <  8.815) );
   OB_REQUIRE( (r[26] >  9.347) && (r[26] <  9.348) );
   OB_REQUIRE( (r[27] >  8.610) && (r[27] <  8.611) );
   OB_REQUIRE( (r[28] >  6.540) && (r[28] <  6.541) );
   OB_REQUIRE( (r[29] >  9.610) && (r[29] <  9.611) );
   OB_REQUIRE( (r[30] > 13.544) && (r[30] < 13.545) );
   OB_REQUIRE( (r[31] >  8.755) && (r[31] <  8.756) );
   OB_REQUIRE( (r[32] >  9.704) && (r[32] <  9.705) );
   OB_REQUIRE( (r[33] >  9.714) && (r[33] <  9.715) );
   OB_REQUIRE( (r[34] >  9.948) && (r[34] <  9.949) );
   OB_REQUIRE( (r[35] > 10.934) && (r[35] < 10.935) );
   OB_REQUIRE( (r[36] >  2.432) && (r[36] <  2.433) );
   OB_REQUIRE( (r[37] >  5.466) && (r[37] <  5.467) );
   OB_REQUIRE( (r[38] >  7.455) && (r[38] <  7.456) );
   OB_REQUIRE( (r[39] >  5.584) && (r[39] <  5.585) );
   OB_REQUIRE( (r[40] >  5.687) && (r[40] <  5.688) );
   OB_REQUIRE( (r[41] >  5.350) && (r[41] <  5.351) );
   OB_REQUIRE( (r[42] >  5.879) && (r[42] <  5.880) );
   OB_REQUIRE( (r[43] >  8.811) && (r[43] <  8.812) );
   OB_REQUIRE( (r[44] > 10.500) && (r[44] < 10.501) );
   OB_REQUIRE( (r[45] >  9.302) && (r[45] <  9.303) );
   OB_REQUIRE( (r[46] >  8.525) && (r[46] <  8.526) );
   OB_REQUIRE( (r[47] >  9.723) && (r[47] <  9.724) );
   OB_REQUIRE( (r[48] > 12.615) && (r[48] < 12.616) );
   OB_REQUIRE( (r[49] >  9.912) && (r[49] <  9.913) );
   OB_REQUIRE( (r[50] >  8.245) && (r[50] <  8.246) );
   OB_REQUIRE( (r[51] > 10.045) && (r[51] < 10.046) );
   OB_REQUIRE( (r[52] > 10.673) && (r[52] < 10.674) );
   OB_REQUIRE( (r[53] > 12.036) && (r[53] < 12.037) );
   OB_REQUIRE( (r[54] >  0.405) && (r[54] <  0.406) );
   OB_REQUIRE( (r[55] >  2.217) && (r[55] <  2.218) );
   OB_REQUIRE( (r[56] >  1.789) && (r[56] <  1.790) );
   OB_REQUIRE( (r[57] >  2.424) && (r[57] <  2.425) );
   OB_REQUIRE( (r[58] >  2.259) && (r[58] <  2.260) );
   OB_REQUIRE( (r[59] >  1.872) && (r[59] <  1.873) );
   OB_REQUIRE( (r[60] >  1.794) && (r[60] <  1.795) );
   OB_REQUIRE( (r[61] >  2.502) && (r[61] <  2.503) );
   OB_REQUIRE( (r[62] >  3.007) && (r[62] <  3.008) );
   OB_REQUIRE( (r[63] >  2.563) && (r[63] <  2.564) );
   OB_REQUIRE( (r[64] >  1.916) && (r[64] <  1.917) );
   OB_REQUIRE( (r[65] >  2.708) && (r[65] <  2.709) );
   OB_REQUIRE( (r[66] >  3.421) && (r[66] <  3.422) );
   OB_REQUIRE( (r[67] >  2.706) && (r[67] <  2.707) );
   OB_REQUIRE( (r[68] >  2.581) && (r[68] <  2.582) );
   OB_REQUIRE( (r[69] >  2.838) && (r[69] <  2.839) );
   OB_REQUIRE( (r[70] >  3.403) && (r[70] <  3.404) );
   OB_REQUIRE( (r[71] >  3.637) && (r[71] <  3.638) );
}



// Mirror stereo probes
void test08(OpenBabel::OBMol* mol)
{
   OpenBabel::OBSpectrophore s;
   s.SetNormalization(OpenBabel::OBSpectrophore::NoNormalization);
   s.SetResolution(3.0);
   s.SetAccuracy(OpenBabel::OBSpectrophore::AngStepSize20);
   s.SetStereo(OpenBabel::OBSpectrophore::MirrorStereoSpecificProbes);
   std::vector<double> r = s.GetSpectrophore(mol);
   OB_REQUIRE( (r[0]  >  1.147) && (r[0]  <  1.148) );
   OB_REQUIRE( (r[1]  >  3.170) && (r[1]  <  3.171) );
   OB_REQUIRE( (r[2]  >  3.174) && (r[2]  <  3.175) );
   OB_REQUIRE( (r[3]  >  3.339) && (r[3]  <  3.340) );
   OB_REQUIRE( (r[4]  >  3.669) && (r[4]  <  3.670) );
   OB_REQUIRE( (r[5]  >  3.108) && (r[5]  <  3.109) );
   OB_REQUIRE( (r[6]  >  3.160) && (r[6]  <  3.161) );
   OB_REQUIRE( (r[7]  >  4.886) && (r[7]  <  4.887) );
   OB_REQUIRE( (r[8]  >  6.239) && (r[8]  <  6.240) );
   OB_REQUIRE( (r[9]  >  5.184) && (r[9]  <  5.185) );
   OB_REQUIRE( (r[10] >  4.438) && (r[10] <  4.439) );
   OB_REQUIRE( (r[11] >  5.232) && (r[11] <  5.233) );
   OB_REQUIRE( (r[12] >  6.973) && (r[12] <  6.974) );
   OB_REQUIRE( (r[13] >  5.545) && (r[13] <  5.546) );
   OB_REQUIRE( (r[14] >  4.827) && (r[14] <  4.828) );
   OB_REQUIRE( (r[15] >  5.804) && (r[15] <  5.805) );
   OB_REQUIRE( (r[16] >  5.809) && (r[16] <  5.810) );
   OB_REQUIRE( (r[17] >  7.107) && (r[17] <  7.108) );
   OB_REQUIRE( (r[18] >  1.986) && (r[18] <  1.987) );
   OB_REQUIRE( (r[19] >  8.473) && (r[19] <  8.474) );
   OB_REQUIRE( (r[20] >  5.681) && (r[20] <  5.682) );
   OB_REQUIRE( (r[21] >  7.358) && (r[21] <  7.359) );
   OB_REQUIRE( (r[22] >  5.874) && (r[22] <  5.875) );
   OB_REQUIRE( (r[23] >  5.997) && (r[23] <  5.998) );
   OB_REQUIRE( (r[24] >  5.596) && (r[24] <  5.597) );
   OB_REQUIRE( (r[25] >  9.362) && (r[25] <  9.363) );
   OB_REQUIRE( (r[26] >  9.389) && (r[26] <  9.390) );
   OB_REQUIRE( (r[27] >  9.566) && (r[27] <  9.567) );
   OB_REQUIRE( (r[28] >  6.162) && (r[28] <  6.163) );
   OB_REQUIRE( (r[29] >  8.745) && (r[29] <  8.746) );
   OB_REQUIRE( (r[30] > 12.694) && (r[30] < 12.695) );
   OB_REQUIRE( (r[31] >  8.925) && (r[31] <  8.926) );
   OB_REQUIRE( (r[32] >  8.854) && (r[32] <  8.855) );
   OB_REQUIRE( (r[33] >  9.599) && (r[33] <  9.600) );
   OB_REQUIRE( (r[34] > 10.087) && (r[34] < 10.088) );
   OB_REQUIRE( (r[35] > 10.454) && (r[35] < 10.455) );
   OB_REQUIRE( (r[36] >  2.371) && (r[36] <  2.372) );
   OB_REQUIRE( (r[37] >  5.658) && (r[37] <  5.659) );
   OB_REQUIRE( (r[38] >  7.224) && (r[38] <  7.225) );
   OB_REQUIRE( (r[39] >  5.556) && (r[39] <  5.557) );
   OB_REQUIRE( (r[40] >  5.617) && (r[40] <  5.618) );
   OB_REQUIRE( (r[41] >  5.240) && (r[41] <  5.241) );
   OB_REQUIRE( (r[42] >  5.576) && (r[42] <  5.577) );
   OB_REQUIRE( (r[43] >  9.171) && (r[43] <  9.172) );
   OB_REQUIRE( (r[44] > 10.523) && (r[44] < 10.524) );
   OB_REQUIRE( (r[45] >  9.109) && (r[45] <  9.110) );
   OB_REQUIRE( (r[46] >  8.459) && (r[46] <  8.460) );
   OB_REQUIRE( (r[47] >  9.413) && (r[47] <  9.414) );
   OB_REQUIRE( (r[48] > 12.418) && (r[48] < 12.419) );
   OB_REQUIRE( (r[49] > 10.191) && (r[49] < 10.192) );
   OB_REQUIRE( (r[50] >  8.444) && (r[50] <  8.445) );
   OB_REQUIRE( (r[51] > 10.295) && (r[51] < 10.296) );
   OB_REQUIRE( (r[52] > 10.572) && (r[52] < 10.573) );
   OB_REQUIRE( (r[53] > 12.187) && (r[53] < 12.188) );
   OB_REQUIRE( (r[54] >  0.371) && (r[54] <  0.372) );
   OB_REQUIRE( (r[55] >  2.238) && (r[55] <  2.239) );
   OB_REQUIRE( (r[56] >  2.022) && (r[56] <  2.023) );
   OB_REQUIRE( (r[57] >  2.330) && (r[57] <  2.331) );
   OB_REQUIRE( (r[58] >  2.009) && (r[58] <  2.010) );
   OB_REQUIRE( (r[59] >  1.968) && (r[59] <  1.969) );
   OB_REQUIRE( (r[60] >  1.978) && (r[60] <  1.979) );
   OB_REQUIRE( (r[61] >  2.682) && (r[61] <  2.683) );
   OB_REQUIRE( (r[62] >  2.877) && (r[62] <  2.878) );
   OB_REQUIRE( (r[63] >  2.494) && (r[63] <  2.495) );
   OB_REQUIRE( (r[64] >  1.869) && (r[64] <  1.870) );
   OB_REQUIRE( (r[65] >  2.450) && (r[65] <  2.451) );
   OB_REQUIRE( (r[66] >  3.385) && (r[66] <  3.386) );
   OB_REQUIRE( (r[67] >  2.565) && (r[67] <  2.566) );
   OB_REQUIRE( (r[68] >  2.589) && (r[68] <  2.590) );
   OB_REQUIRE( (r[69] >  3.011) && (r[69] <  3.012) );
   OB_REQUIRE( (r[70] >  3.380) && (r[70] <  3.381) );
   OB_REQUIRE( (r[71] >  3.848) && (r[71] <  3.849) );
}



// All stereo probes
void test09(OpenBabel::OBMol* mol)
{
   OpenBabel::OBSpectrophore s;
   s.SetNormalization(OpenBabel::OBSpectrophore::NoNormalization);
   s.SetResolution(3.0);
   s.SetAccuracy(OpenBabel::OBSpectrophore::AngStepSize20);
   s.SetStereo(OpenBabel::OBSpectrophore::AllStereoSpecificProbes);
   std::vector<double> r = s.GetSpectrophore(mol);
   OB_REQUIRE( (r[0]   >  1.150) && (r[0]   <  1.151) );
   OB_REQUIRE( (r[1]   >  3.149) && (r[1]   <  3.150) );
   OB_REQUIRE( (r[2]   >  3.209) && (r[2]   <  3.210) );
   OB_REQUIRE( (r[3]   >  3.174) && (r[3]   <  3.175) );
   OB_REQUIRE( (r[4]   >  3.569) && (r[4]   <  3.570) );
   OB_REQUIRE( (r[5]   >  3.304) && (r[5]   <  3.305) );
   OB_REQUIRE( (r[6]   >  3.051) && (r[6]   <  3.052) );
   OB_REQUIRE( (r[7]   >  4.562) && (r[7]   <  4.563) );
   OB_REQUIRE( (r[8]   >  6.095) && (r[8]   <  6.096) );
   OB_REQUIRE( (r[9]   >  5.226) && (r[9]   <  5.227) );
   OB_REQUIRE( (r[10]  >  4.518) && (r[10]  <  4.519) );
   OB_REQUIRE( (r[11]  >  5.049) && (r[11]  <  5.050) );
   OB_REQUIRE( (r[12]  >  6.967) && (r[12]  <  6.968) );
   OB_REQUIRE( (r[13]  >  5.714) && (r[13]  <  5.715) );
   OB_REQUIRE( (r[14]  >  4.741) && (r[14]  <  4.742) );
   OB_REQUIRE( (r[15]  >  5.948) && (r[15]  <  5.949) );
   OB_REQUIRE( (r[16]  >  5.890) && (r[16]  <  5.891) );
   OB_REQUIRE( (r[17]  >  7.335) && (r[17]  <  7.336) );
   OB_REQUIRE( (r[18]  >  1.147) && (r[18]  <  1.148) );
   OB_REQUIRE( (r[19]  >  3.170) && (r[19]  <  3.171) );
   OB_REQUIRE( (r[20]  >  3.174) && (r[20]  <  3.175) );
   OB_REQUIRE( (r[21]  >  3.339) && (r[21]  <  3.340) );
   OB_REQUIRE( (r[22]  >  3.669) && (r[22]  <  3.670) );
   OB_REQUIRE( (r[23]  >  3.108) && (r[23]  <  3.109) );
   OB_REQUIRE( (r[24]  >  3.160) && (r[24]  <  3.161) );
   OB_REQUIRE( (r[25]  >  4.886) && (r[25]  <  4.887) );
   OB_REQUIRE( (r[26]  >  6.239) && (r[26]  <  6.240) );
   OB_REQUIRE( (r[27]  >  5.184) && (r[27]  <  5.185) );
   OB_REQUIRE( (r[28]  >  4.438) && (r[28]  <  4.439) );
   OB_REQUIRE( (r[29]  >  5.232) && (r[29]  <  5.233) );
   OB_REQUIRE( (r[30]  >  6.973) && (r[30]  <  6.974) );
   OB_REQUIRE( (r[31]  >  5.545) && (r[31]  <  5.546) );
   OB_REQUIRE( (r[32]  >  4.827) && (r[32]  <  4.828) );
   OB_REQUIRE( (r[33]  >  5.804) && (r[33]  <  5.805) );
   OB_REQUIRE( (r[34]  >  5.809) && (r[34]  <  5.810) );
   OB_REQUIRE( (r[35]  >  7.107) && (r[35]  <  7.108) );
   OB_REQUIRE( (r[36]  >  2.015) && (r[36]  <  2.016) );
   OB_REQUIRE( (r[37]  >  8.571) && (r[37]  <  8.572) );
   OB_REQUIRE( (r[38]  >  5.878) && (r[38]  <  5.879) );
   OB_REQUIRE( (r[39]  >  7.161) && (r[39]  <  7.162) );
   OB_REQUIRE( (r[40]  >  5.824) && (r[40]  <  5.825) );
   OB_REQUIRE( (r[41]  >  5.744) && (r[41]  <  5.745) );
   OB_REQUIRE( (r[42]  >  5.676) && (r[42]  <  5.677) );
   OB_REQUIRE( (r[43]  >  8.814) && (r[43]  <  8.815) );
   OB_REQUIRE( (r[44]  >  9.347) && (r[44]  <  9.348) );
   OB_REQUIRE( (r[45]  >  8.610) && (r[45]  <  8.611) );
   OB_REQUIRE( (r[46]  >  6.540) && (r[46]  <  6.541) );
   OB_REQUIRE( (r[47]  >  9.610) && (r[47]  <  9.611) );
   OB_REQUIRE( (r[48]  > 13.544) && (r[48]  < 13.545) );
   OB_REQUIRE( (r[49]  >  8.755) && (r[49]  <  8.756) );
   OB_REQUIRE( (r[50]  >  9.704) && (r[50]  <  9.705) );
   OB_REQUIRE( (r[51]  >  9.714) && (r[51]  <  9.715) );
   OB_REQUIRE( (r[52]  >  9.948) && (r[52]  <  9.949) );
   OB_REQUIRE( (r[53]  > 10.934) && (r[53]  < 10.935) );
   OB_REQUIRE( (r[54]  >  1.986) && (r[54]  <  1.987) );
   OB_REQUIRE( (r[55]  >  8.473) && (r[55]  <  8.474) );
   OB_REQUIRE( (r[56]  >  5.681) && (r[56]  <  5.682) );
   OB_REQUIRE( (r[57]  >  7.358) && (r[57]  <  7.359) );
   OB_REQUIRE( (r[58]  >  5.874) && (r[58]  <  5.875) );
   OB_REQUIRE( (r[59]  >  5.997) && (r[59]  <  5.998) );
   OB_REQUIRE( (r[60]  >  5.596) && (r[60]  <  5.597) );
   OB_REQUIRE( (r[61]  >  9.362) && (r[61]  <  9.363) );
   OB_REQUIRE( (r[62]  >  9.389) && (r[62]  <  9.390) );
   OB_REQUIRE( (r[63]  >  9.566) && (r[63]  <  9.567) );
   OB_REQUIRE( (r[64]  >  6.162) && (r[64]  <  6.163) );
   OB_REQUIRE( (r[65]  >  8.745) && (r[65]  <  8.746) );
   OB_REQUIRE( (r[66]  > 12.694) && (r[66]  < 12.695) );
   OB_REQUIRE( (r[67]  >  8.925) && (r[67]  <  8.926) );
   OB_REQUIRE( (r[68]  >  8.854) && (r[68]  <  8.855) );
   OB_REQUIRE( (r[69]  >  9.599) && (r[69]  <  9.600) );
   OB_REQUIRE( (r[70]  > 10.087) && (r[70]  < 10.088) );
   OB_REQUIRE( (r[71]  > 10.454) && (r[71]  < 10.455) );
   OB_REQUIRE( (r[72]  >  2.432) && (r[72]  <  2.433) );
   OB_REQUIRE( (r[73]  >  5.466) && (r[73]  <  5.467) );
   OB_REQUIRE( (r[74]  >  7.455) && (r[74]  <  7.456) );
   OB_REQUIRE( (r[75]  >  5.584) && (r[75]  <  5.585) );
   OB_REQUIRE( (r[76]  >  5.687) && (r[76]  <  5.688) );
   OB_REQUIRE( (r[77]  >  5.350) && (r[77]  <  5.351) );
   OB_REQUIRE( (r[78]  >  5.879) && (r[78]  <  5.880) );
   OB_REQUIRE( (r[79]  >  8.811) && (r[79]  <  8.812) );
   OB_REQUIRE( (r[80]  > 10.500) && (r[80]  < 10.501) );
   OB_REQUIRE( (r[81]  >  9.302) && (r[81]  <  9.304) );
   OB_REQUIRE( (r[82]  >  8.525) && (r[82]  <  8.526) );
   OB_REQUIRE( (r[83]  >  9.723) && (r[83]  <  9.724) );
   OB_REQUIRE( (r[84]  > 12.615) && (r[84]  < 12.616) );
   OB_REQUIRE( (r[85]  >  9.912) && (r[85]  <  9.913) );
   OB_REQUIRE( (r[86]  >  8.245) && (r[86]  <  8.246) );
   OB_REQUIRE( (r[87]  > 10.045) && (r[87]  < 10.046) );
   OB_REQUIRE( (r[88]  > 10.673) && (r[88]  < 10.674) );
   OB_REQUIRE( (r[89]  > 12.036) && (r[89]  < 12.037) );
   OB_REQUIRE( (r[90]  >  2.371) && (r[90]  <  2.372) );
   OB_REQUIRE( (r[91]  >  5.658) && (r[91]  <  5.659) );
   OB_REQUIRE( (r[92]  >  7.224) && (r[92]  <  7.225) );
   OB_REQUIRE( (r[93]  >  5.556) && (r[93]  <  5.557) );
   OB_REQUIRE( (r[94]  >  5.617) && (r[94]  <  5.618) );
   OB_REQUIRE( (r[95]  >  5.240) && (r[95]  <  5.241) );
   OB_REQUIRE( (r[96]  >  5.576) && (r[96]  <  5.577) );
   OB_REQUIRE( (r[97]  >  9.171) && (r[97]  <  9.172) );
   OB_REQUIRE( (r[98]  > 10.523) && (r[98]  < 10.524) );
   OB_REQUIRE( (r[99]  >  9.109) && (r[99]  <  9.110) );
   OB_REQUIRE( (r[100] >  8.459) && (r[100] <  8.460) );
   OB_REQUIRE( (r[101] >  9.413) && (r[101] <  9.414) );
   OB_REQUIRE( (r[102] > 12.418) && (r[102] < 12.419) );
   OB_REQUIRE( (r[103] > 10.191) && (r[103] < 10.192) );
   OB_REQUIRE( (r[104] >  8.444) && (r[104] <  8.445) );
   OB_REQUIRE( (r[105] > 10.295) && (r[105] < 10.296) );
   OB_REQUIRE( (r[106] > 10.572) && (r[106] < 10.573) );
   OB_REQUIRE( (r[107] > 12.187) && (r[107] < 12.188) );
   OB_REQUIRE( (r[108] >  0.405) && (r[108] <  0.406) );
   OB_REQUIRE( (r[109] >  2.217) && (r[109] <  2.218) );
   OB_REQUIRE( (r[110] >  1.789) && (r[110] <  1.790) );
   OB_REQUIRE( (r[111] >  2.424) && (r[111] <  2.425) );
   OB_REQUIRE( (r[112] >  2.259) && (r[112] <  2.260) );
   OB_REQUIRE( (r[113] >  1.872) && (r[113] <  1.873) );
   OB_REQUIRE( (r[114] >  1.794) && (r[114] <  1.795) );
   OB_REQUIRE( (r[115] >  2.502) && (r[115] <  2.503) );
   OB_REQUIRE( (r[116] >  3.007) && (r[116] <  3.008) );
   OB_REQUIRE( (r[117] >  2.563) && (r[117] <  2.564) );
   OB_REQUIRE( (r[118] >  1.916) && (r[118] <  1.917) );
   OB_REQUIRE( (r[119] >  2.708) && (r[119] <  2.709) );
   OB_REQUIRE( (r[120] >  3.421) && (r[120] <  3.422) );
   OB_REQUIRE( (r[121] >  2.706) && (r[121] <  2.707) );
   OB_REQUIRE( (r[122] >  2.581) && (r[122] <  2.582) );
   OB_REQUIRE( (r[123] >  2.838) && (r[123] <  2.839) );
   OB_REQUIRE( (r[124] >  3.403) && (r[124] <  3.404) );
   OB_REQUIRE( (r[125] >  3.637) && (r[125] <  3.638) );
   OB_REQUIRE( (r[126] >  0.371) && (r[126] <  0.372) );
   OB_REQUIRE( (r[127] >  2.238) && (r[127] <  2.239) );
   OB_REQUIRE( (r[128] >  2.022) && (r[128] <  2.023) );
   OB_REQUIRE( (r[129] >  2.330) && (r[129] <  2.331) );
   OB_REQUIRE( (r[130] >  2.009) && (r[130] <  2.010) );
   OB_REQUIRE( (r[131] >  1.968) && (r[131] <  1.969) );
   OB_REQUIRE( (r[132] >  1.978) && (r[132] <  1.979) );
   OB_REQUIRE( (r[133] >  2.682) && (r[133] <  2.683) );
   OB_REQUIRE( (r[134] >  2.877) && (r[134] <  2.878) );
   OB_REQUIRE( (r[135] >  2.494) && (r[135] <  2.495) );
   OB_REQUIRE( (r[136] >  1.869) && (r[136] <  1.870) );
   OB_REQUIRE( (r[137] >  2.450) && (r[137] <  2.451) );
   OB_REQUIRE( (r[138] >  3.385) && (r[138] <  3.386) );
   OB_REQUIRE( (r[139] >  2.565) && (r[139] <  2.566) );
   OB_REQUIRE( (r[140] >  2.589) && (r[140] <  2.590) );
   OB_REQUIRE( (r[141] >  3.011) && (r[141] <  3.012) );
   OB_REQUIRE( (r[142] >  3.380) && (r[142] <  3.381) );
   OB_REQUIRE( (r[143] >  3.848) && (r[143] <  3.849) );
}



int spectrophoretest(int argc, char* argv[])
{
  int defaultchoice = 1;
  
  int choice = defaultchoice;

  if (argc > 1) {
    if(sscanf(argv[1], "%d", &choice) != 1) {
      printf("Couldn't parse that input as a number\n");
      return -1;
    }
  }

   // Create a test molecule
   OpenBabel::OBMol mol;
   OpenBabel::OBAtom* a[5];
   a[0] = mol.NewAtom(); a[0]->SetAtomicNum(6);  a[0]->SetVector(-0.013,  1.086,  0.008);
   a[1] = mol.NewAtom(); a[1]->SetAtomicNum(1);  a[1]->SetVector( 0.002, -0.004,  0.002);
   a[2] = mol.NewAtom(); a[2]->SetAtomicNum(9);  a[2]->SetVector( 1.300,  1.570, -0.002);
   a[3] = mol.NewAtom(); a[3]->SetAtomicNum(35); a[3]->SetVector(-0.964,  1.737, -1.585);
   a[4] = mol.NewAtom(); a[4]->SetAtomicNum(17); a[4]->SetVector(-0.857,  1.667,  1.491);
   OpenBabel::OBBond* b;
   for (int i(1); i < 5; ++i)
   {
      b = mol.NewBond();
      b->SetBegin(a[0]); b->SetEnd(a[i]); b->SetBondOrder(1);
   }
   
  switch(choice) {
  case 1:
    test01(&mol);
    test02(&mol);
    break;
  case 2:
    test03(&mol);
    test04(&mol);
    break;
  case 3:
    test05(&mol);
    test06(&mol);
    break;
  case 4:
    test07(&mol);
    test08(&mol);
    break;
  case 5:
    test09(&mol);
    break;
  default:
    std::cout << "Test number " << choice << " does not exist!\n";
    return -1;
  }

  return 0;
}
