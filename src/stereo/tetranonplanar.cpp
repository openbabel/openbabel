#include <openbabel/stereo/tetranonplanar.h>
#include <openbabel/oberror.h>
#include <cassert>

namespace OpenBabel {

  OBTetraNonPlanarStereo::OBTetraNonPlanarStereo(OBMol *mol) : OBStereoBase(mol)
  {
  }

  OBTetraNonPlanarStereo::~OBTetraNonPlanarStereo()
  {
  }

  std::vector<unsigned long> OBTetraNonPlanarStereo::ToInternal(const std::vector<unsigned long> &refs, 
      unsigned long id, OBStereo::Winding winding, OBStereo::View view)
  {
    assert( refs.size() == 3 );
    assert( (winding == OBStereo::Clockwise) || (winding == OBStereo::AntiClockwise) );
    assert( (view == OBStereo::ViewFrom) || (view == OBStereo::ViewTowards) );

    // copy the refs
    std::vector<unsigned long> result(refs);
    // add the from/to atom to the front
    result.insert(result.begin(), id);

    bool odd = false;
    // clockwise <-> anti-clockwise : odd = true
    if (winding == OBStereo::AntiClockwise)
      odd = !odd;
    // ViewFrom <-> ViewTowards : odd = true
    if (view == OBStereo::ViewTowards)
      odd = !odd;

    if (odd)
      return Permutate(result, 2, 3);
    else
      return result;
  }
  
  std::vector<unsigned long> OBTetraNonPlanarStereo::ToView(const std::vector<unsigned long> &refs,
      unsigned long id, OBStereo::Winding winding, OBStereo::View view)
  {
    assert( refs.size() == 4 );
    assert( id != OBStereo::NoId );
    assert( (winding == OBStereo::Clockwise) || (winding == OBStereo::AntiClockwise) );
    assert( (view == OBStereo::ViewFrom) || (view == OBStereo::ViewTowards) );

    // copy the internal refs
    std::vector<unsigned long> result(refs);
 
    // keep track of the permuations by using the oddness
    bool odd = false;

    // find id
    if (refs.at(0) == id) {
      // view from/towards atom is same as internal, no permutation needed
      result.erase(result.begin());
      // no permuations performed --> odd = false
      odd = false;
    } else {
      // move id to front and remove it = 1 permutation
      for (int i = 1; i < 4; ++i) {
        if (refs.at(i) == id) {
          result.erase(result.begin() + i);
          break;
        }
      }
      // 1 permutation perfromed --> odd = true
      odd = true;
    }

    // clockwise <-> anti-clockwise : odd = true
    if (winding == OBStereo::AntiClockwise)
      odd = !odd;
    // ViewFrom <-> ViewTowards : odd = true
    if (view == OBStereo::ViewTowards)
      odd = !odd;

    // make sure we actually found id
    if (result.size() == 3) {
      if (odd)
        return Permutate(result, 1, 2);
      else
        return result;
    }

    obErrorLog.ThrowError(__FUNCTION__, 
        "OBTetraNonPlanarStereo::ToInternal : Paramter id not found in internal refs.", obError);
    return std::vector<unsigned long>();
  }
 
  int OBTetraNonPlanarStereo::NumInversions(const std::vector<unsigned long> &refs)
  {
    std::vector<int> invVec; // the inversion vector
    std::vector<unsigned long>::const_iterator i, j;
    for (i = refs.begin(); i != refs.end(); ++i) {
      int e = 0; // ith element
      // loop over elements to the right
      for (j = i; j != refs.end(); ++j)
        // increment e if element to the right is lower
        if (*j < *i)
          e++;
        
      invVec.push_back(e);
    }
    int sum = 0;
    for (std::vector<int>::iterator k = invVec.begin(); k != invVec.end(); ++k)
      sum += *k;
    return sum;
  }

  std::vector<unsigned long> OBTetraNonPlanarStereo::Permutate(const std::vector<unsigned long> &refs, int i, int j)
  {
    assert( i >= 0 && i < refs.size() );
    assert( j >= 0 && j < refs.size() );
    std::vector<unsigned long> result(refs);
    result[i] = refs.at(j);
    result[j] = refs.at(i);
    return result;   
  }

}

