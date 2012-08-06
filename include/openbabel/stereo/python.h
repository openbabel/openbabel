#include <openbabel/stereo/stereo.h>

namespace OpenBabel {

  struct OBTetrahedralConfig
  {
#ifndef SWIG
    static OBTetrahedralStereo::Config Convert(const OBTetrahedralConfig &other)
    {
      OBTetrahedralStereo::Config config;
      config.center = other.center;
      config.from = other.from_or_towards;
      config.refs = other.refs;
      config.winding = other.winding;
      config.view = other.view;
      config.specified = other.specified;
      return config;
    }
#endif

    bool operator==(const OBTetrahedralConfig &other) const
    {
      if (center != other.center)
        return false;
      if ((refs.size() != 3) || (other.refs.size() != 3))
        return false;
      // return true if either is unspecified (i.e. accidental)
      if (!specified || !other.specified)
        return true;

      // Convert both Config's refs to same from, winding and view while
      // avoiding having an ImplicitRef in the 'from' position of either
      OBTetrahedralStereo::Config thisConfig = Convert(*this), otherConfig = Convert(other);

      if (from_or_towards == OBStereo::ImplicitRef) {
        thisConfig = OBTetraNonPlanarStereo::ToConfig(thisConfig, refs[0], winding, view);
        otherConfig = OBTetraNonPlanarStereo::ToConfig(otherConfig, thisConfig.from, winding, view);
      }
      else if (other.from_or_towards == OBStereo::ImplicitRef) {
        otherConfig = OBTetraNonPlanarStereo::ToConfig(otherConfig, other.refs[0], winding, view);
        thisConfig = OBTetraNonPlanarStereo::ToConfig(thisConfig, otherConfig.from, winding, view);
      }
      else {
        otherConfig = OBTetraNonPlanarStereo::ToConfig(otherConfig, thisConfig.from, winding, view);
      }

      if (!OBStereo::ContainsSameRefs(thisConfig.refs, otherConfig.refs)) {
        if (OBStereo::ContainsRef(thisConfig.refs, OBStereo::ImplicitRef)) {
          // if both refs already contain ImplicitRef, return false
          if (OBStereo::ContainsRef(otherConfig.refs, OBStereo::ImplicitRef))
            return false;

          // example: *this       = 23H
          //          otherConfig = 234 --> 23H

          // for each ref in otherConfig
          for (unsigned int i = 0; i < otherConfig.refs.size(); ++i) {
            bool found = false;
            for (OBStereo::RefIter j = thisConfig.refs.begin(); j != thisConfig.refs.end(); ++j)
              if (otherConfig.refs.at(i) == *j)
                found = true;

            if (!found) {
              // the ref from otherConfig is not found in this config
              otherConfig.refs[i] = OBStereo::ImplicitRef;
              break;
            }
          }
        } else
          if (OBStereo::ContainsRef(otherConfig.refs, OBStereo::ImplicitRef)) {
            // if both refs already contain ImplicitRef, return false
            if (OBStereo::ContainsRef(thisConfig.refs, OBStereo::ImplicitRef))
              return false;

            // example: *this       = 234
            //          otherConfig = 23H --> 234

            // for each ref in *this
            for (unsigned int i = 0; i < thisConfig.refs.size(); ++i) {
              bool found = false;
              // for each refs in otherConfig
              for (OBStereo::RefIter j = otherConfig.refs.begin(); j != otherConfig.refs.end(); ++j)
                if (thisConfig.refs.at(i) == *j)
                  found = true;

              if (!found) {
                for (OBStereo::RefIter j = otherConfig.refs.begin(); j != otherConfig.refs.end(); ++j)
                  if (*j == OBStereo::ImplicitRef)
                    *j = thisConfig.refs.at(i);
                break;
              }
            }
          }
      }

      int Ni1 = OBStereo::NumInversions(thisConfig.refs);
      int Ni2 = OBStereo::NumInversions(otherConfig.refs);
      return ((Ni1 + Ni2) % 2 == 0);
    }

    bool operator!=(const OBTetrahedralConfig &other) const
    {
      return !(*this == other);
    }

    unsigned long center;
    unsigned long from_or_towards;
    OBStereo::Refs refs;
    OBStereo::Winding winding;
    OBStereo::View view;
    bool specified;
  };

  struct OBCisTransConfig
  {
#ifndef SWIG
    static OBCisTransStereo::Config Convert(const OBCisTransConfig &other)
    {
      OBCisTransStereo::Config config;
      config.begin = other.begin;
      config.end = other.end;
      config.refs = other.refs;
      config.shape = other.shape;
      config.specified = other.specified;
      return config;
    }
#endif

    bool operator==(const OBCisTransConfig &other) const
    {
      if ((begin != other.begin) && (begin != other.end))
        return false;
      if ((end != other.begin) && (end != other.end))
        return false;
      if ((refs.size() != 4) || (other.refs.size() != 4))
        return false;

      OBCisTransStereo::Config u1, u2;
      if (!OBStereo::ContainsSameRefs(refs, other.refs)) {
        // find a ref that occurs in both
        for (OBStereo::ConstRefIter i = refs.begin(); i != refs.end(); ++i)
          if (OBStereo::ContainsRef(other.refs, *i)) {
            u1 = OBTetraPlanarStereo::ToConfig(Convert(*this), *i, OBStereo::ShapeU); // refs[0] = u1.refs[0]
            u2 = OBTetraPlanarStereo::ToConfig(Convert(other), *i, OBStereo::ShapeU); // refs[0] = u2.refs[0]
          }

        // check if they actualy share an id...
        if (u1.refs.empty())
          return false;
      } else {
        // normalize the other Config struct
        u1 = OBTetraPlanarStereo::ToConfig(Convert(*this), refs.at(0), OBStereo::ShapeU); // refs[0] = u1.refs[0]
        u2 = OBTetraPlanarStereo::ToConfig(Convert(other), refs.at(0), OBStereo::ShapeU); // refs[0] = u2.refs[0]
        // both now start with the same ref
        return (u1.refs[2] == u2.refs[2]);
      }

      if ((u1.refs[2] == OBStereo::ImplicitRef) || (u2.refs[2] == OBStereo::ImplicitRef)) {
        // 1 2 H 4
        if ((u1.refs[3] == OBStereo::ImplicitRef) || (u2.refs[3] == OBStereo::ImplicitRef)) {
          return (u1.refs[1] == u2.refs[1]); // 1 2 H H
        } else {
          return (u1.refs[3] == u2.refs[3]); // 1 H H 4
        }
      } else
        return (u1.refs[2] == u2.refs[2]); // 1 2 3 4  &  1 H 3 4  &  1 2 3 H

      return false;
    }

    bool operator!=(const OBCisTransConfig &other) const
    {
      return !(*this == other);
    }

    unsigned long begin, end;
    OBStereo::Refs refs;
    OBStereo::Shape shape;
    bool specified;
  };

  struct OBSquarePlanarConfig
  {
#ifndef SWIG
    static OBSquarePlanarStereo::Config Convert(const OBSquarePlanarConfig &other)
    {
      OBSquarePlanarStereo::Config config;
      config.center = other.center;
      config.refs = other.refs;
      config.shape = other.shape;
      config.specified = other.specified;
      return config;
    }
#endif

    bool operator==(const OBSquarePlanarConfig &other) const
    {
      if (center != other.center)
        return false;
      if ((refs.size() != 4) || (other.refs.size() != 4))
        return false;

      OBSquarePlanarStereo::Config u1, u2;
      if (!OBStereo::ContainsSameRefs(refs, other.refs)) {
        // find a ref that occurs in both
        for (OBStereo::ConstRefIter i = refs.begin(); i != refs.end(); ++i)
          if (OBStereo::ContainsRef(other.refs, *i)) {
            u1 = OBTetraPlanarStereo::ToConfig(Convert(*this), *i, OBStereo::ShapeU); // refs[0] = u1.refs[0]
            u2 = OBTetraPlanarStereo::ToConfig(Convert(other), *i, OBStereo::ShapeU); // refs[0] = u2.refs[0]
          }

        // check if they actualy share an id...
        if (u1.refs.empty())
          return false;
      } else {
        // normalize the other Config struct
        u1 = OBTetraPlanarStereo::ToConfig(Convert(*this), refs.at(0), OBStereo::ShapeU); // refs[0] = u1.refs[0]
        u2 = OBTetraPlanarStereo::ToConfig(Convert(other), refs.at(0), OBStereo::ShapeU); // refs[0] = u2.refs[0]
        // both now start with the same ref
        return (u1.refs[2] == u2.refs[2]);
      }

      if ((u1.refs[2] == OBStereo::ImplicitRef) || (u2.refs[2] == OBStereo::ImplicitRef)) {
        // 1 2 H 4
        if ((u1.refs[3] == OBStereo::ImplicitRef) || (u2.refs[3] == OBStereo::ImplicitRef)) {
          return (u1.refs[1] == u2.refs[1]); // 1 2 H H
        } else {
          return (u1.refs[3] == u2.refs[3]); // 1 H H 4
        }
      } else
        return (u1.refs[2] == u2.refs[2]); // 1 2 3 4  &  1 H 3 4  &  1 2 3 H

      return false;
    }

    bool operator!=(const OBSquarePlanarConfig &other) const
    {
      return !(*this == other);
    }

    unsigned long center;
    OBStereo::Refs refs;
    OBStereo::Shape shape;
    bool specified;
  };


}
