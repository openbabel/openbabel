#include <openbabel/stereo/tetraplanar.h>
#include <cassert>

namespace OpenBabel {

  OBTetraPlanarStereo::OBTetraPlanarStereo(OBMol *mol) : OBStereoBase(mol)
  {
  }

  OBTetraPlanarStereo::~OBTetraPlanarStereo()
  {
  }

  std::vector<unsigned long> OBTetraPlanarStereo::ToInternal(const std::vector<unsigned long> &refs, 
      OBStereo::Shape shape)
  {
    assert( refs.size() == 4 );
    std::vector<unsigned long> result(refs);

    switch (shape) {
      case OBStereo::ShapeU:
        // same as internal, just copy
        return result;
      case OBStereo::ShapeZ:
        // normalize to U shape
        result[1] = refs.at(2);
        result[2] = refs.at(3);
        result[3] = refs.at(1);
        return result;
      case OBStereo::Shape4:
        // normalize to U shape
        result[1] = refs.at(2);
        result[2] = refs.at(1);
        return result;
    }
  
  }
  
  std::vector<unsigned long> OBTetraPlanarStereo::ToShape(const std::vector<unsigned long> &refs, 
      OBStereo::Shape shape)
  {
    assert( refs.size() == 4 );
    std::vector<unsigned long> result(refs);

    switch (shape) {
      case OBStereo::ShapeU:
        // same as internal, just copy
        return result;
      case OBStereo::ShapeZ:
        // convert to U shape
        result[1] = refs.at(3);
        result[2] = refs.at(1);
        result[3] = refs.at(2);
        return result;
      case OBStereo::Shape4:
        // normalize to U shape
        result[1] = refs.at(2);
        result[2] = refs.at(1);
        return result;
    }

  }
 
}

