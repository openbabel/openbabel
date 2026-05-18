#ifndef OB_STEREO_TEST_HELPERS_H
#define OB_STEREO_TEST_HELPERS_H

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

#include <sstream>
#include <string>

// Convert a 3D OBMol to canonical SMILES by doing a full SDF roundtrip so
// that StereoFrom3D is invoked and the SMILES reflects the 3D stereo.
inline std::string canSmiFrom3D(OpenBabel::OBMol& mol3D)
{
  OpenBabel::OBConversion conv;
  conv.SetInAndOutFormats("sdf", "can");

  std::ostringstream sdfBuf;
  conv.SetOutFormat("sdf");
  conv.Write(&mol3D, &sdfBuf);

  OpenBabel::OBMol mol2D;
  conv.SetInFormat("sdf");
  std::istringstream iss(sdfBuf.str());
  conv.Read(&mol2D, &iss);

  conv.SetOutFormat("can");
  std::string result = conv.WriteString(&mol2D, true);
  while (!result.empty() &&
         (result.back() == '\n' || result.back() == '\r' || result.back() == '\t'))
    result.pop_back();
  return result;
}

// Flip every tetrahedral chirality marker in a SMILES (@ <-> @@); returns
// the enantiomer.  E/Z bond directions (/ and \) are left alone because
// mirror reflection does not change cis/trans.
inline std::string makeEnantiomer(const std::string& smiles)
{
  std::string result;
  result.reserve(smiles.size() * 2);  // worst case: every '@' becomes "@@"
  for (size_t i = 0; i < smiles.size(); ++i) {
    if (smiles[i] != '@') {
      result += smiles[i];
      continue;
    }
    if (i + 1 < smiles.size() && smiles[i + 1] == '@') {
      result += '@';        // @@ -> @
      ++i;
    } else {
      result += "@@";       // @  -> @@
    }
  }
  return result;
}

// Flip only the Nth (0-indexed) tetrahedral chirality marker, producing one
// specific diastereomer.  For >= 2 stereocenters this yields a structure
// distinct from both the input and its enantiomer.  index must be less than
// the number of @ markers in the input.
inline std::string makeDiastereomer(const std::string& smiles, size_t index)
{
  std::string result = smiles;
  size_t count = 0;
  for (size_t i = 0; i < result.size(); ++i) {
    if (result[i] != '@')
      continue;
    bool isDouble = (i + 1 < result.size() && result[i + 1] == '@');
    if (count == index) {
      if (isDouble)
        result.erase(i, 1);            // @@ -> @
      else
        result.insert(i, 1, '@');      // @  -> @@
      break;
    }
    ++count;
    if (isDouble)
      ++i;
  }
  return result;
}

// Count tetrahedral chirality markers (each @ or @@ counts once).
inline size_t countChirality(const std::string& smiles)
{
  size_t count = 0;
  for (size_t i = 0; i < smiles.size(); ++i) {
    if (smiles[i] == '@') {
      ++count;
      if (i + 1 < smiles.size() && smiles[i + 1] == '@')
        ++i;
    }
  }
  return count;
}

#endif
