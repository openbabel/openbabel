%include <openbabel/stereo/tetranonplanar.h>
%include <openbabel/stereo/tetraplanar.h>
%include <openbabel/stereo/tetrahedral.h>
%include <openbabel/stereo/cistrans.h>
%include <openbabel/stereo/squareplanar.h>
%include <openbabel/stereo/bindings.h>

%extend OpenBabel::OBTetrahedralStereo {

  void SetConfig(const OpenBabel::OBTetrahedralConfig &config)
  {
    self->SetConfig(OpenBabel::OBTetrahedralConfig::Convert(config));
  }

  OpenBabel::OBTetrahedralConfig GetConfig(OBStereo::Winding winding = OBStereo::Clockwise, OBStereo::View view = OBStereo::ViewFrom)
  {
    OpenBabel::OBTetrahedralStereo::Config cConfig = self->GetConfig(winding, view);
    
    OpenBabel::OBTetrahedralConfig pyConfig;
    pyConfig.center = cConfig.center;
    pyConfig.from_or_towards = cConfig.from;
    pyConfig.refs = cConfig.refs;
    pyConfig.winding = cConfig.winding;
    pyConfig.view = cConfig.view;
    pyConfig.specified = cConfig.specified;

    return pyConfig;
  }
  
  OpenBabel::OBTetrahedralConfig GetConfig(unsigned long from_or_towards, OBStereo::Winding winding = OBStereo::Clockwise, OBStereo::View view = OBStereo::ViewFrom)
  {
    OpenBabel::OBTetrahedralStereo::Config cConfig = self->GetConfig(from_or_towards, winding, view);
    
    OpenBabel::OBTetrahedralConfig pyConfig;
    pyConfig.center = cConfig.center;
    pyConfig.from_or_towards = cConfig.from;
    pyConfig.refs = cConfig.refs;
    pyConfig.winding = cConfig.winding;
    pyConfig.view = cConfig.view;
    pyConfig.specified = cConfig.specified;

    return pyConfig;
  }

}

%extend OpenBabel::OBCisTransStereo {

  void SetConfig(const OpenBabel::OBCisTransConfig &config)
  {
    self->SetConfig(OpenBabel::OBCisTransConfig::Convert(config));
  }

  OpenBabel::OBCisTransConfig GetConfig(OBStereo::Shape shape = OBStereo::ShapeU)
  {
    OpenBabel::OBCisTransStereo::Config cConfig = self->GetConfig(shape);
  
    OpenBabel::OBCisTransConfig pyConfig;
    pyConfig.begin = cConfig.begin;
    pyConfig.end = cConfig.end;
    pyConfig.refs = cConfig.refs;
    pyConfig.shape = cConfig.shape;
    pyConfig.specified = cConfig.specified;

    return pyConfig;
  }

  OpenBabel::OBCisTransConfig GetConfig(unsigned long start, OBStereo::Shape shape = OBStereo::ShapeU)
  {
    OpenBabel::OBCisTransStereo::Config cConfig = self->GetConfig(start, shape);
  
    OpenBabel::OBCisTransConfig pyConfig;
    pyConfig.begin = cConfig.begin;
    pyConfig.end = cConfig.end;
    pyConfig.refs = cConfig.refs;
    pyConfig.shape = cConfig.shape;
    pyConfig.specified = cConfig.specified;

    return pyConfig;
  }

}

%extend OpenBabel::OBSquarePlanarStereo {

  void SetConfig(const OpenBabel::OBSquarePlanarConfig &config)
  {
    self->SetConfig(OpenBabel::OBSquarePlanarConfig::Convert(config));
  }

  OpenBabel::OBSquarePlanarConfig GetConfig(OBStereo::Shape shape = OBStereo::ShapeU)
  {
    OpenBabel::OBSquarePlanarStereo::Config cConfig = self->GetConfig(shape);
  
    OpenBabel::OBSquarePlanarConfig pyConfig;
    pyConfig.center = cConfig.center;
    pyConfig.refs = cConfig.refs;
    pyConfig.shape = cConfig.shape;
    pyConfig.specified = cConfig.specified;

    return pyConfig;
  }

  OpenBabel::OBSquarePlanarConfig GetConfig(unsigned long start, OBStereo::Shape shape = OBStereo::ShapeU)
  {
    OpenBabel::OBSquarePlanarStereo::Config cConfig = self->GetConfig(start, shape);
  
    OpenBabel::OBSquarePlanarConfig pyConfig;
    pyConfig.center = cConfig.center;
    pyConfig.refs = cConfig.refs;
    pyConfig.shape = cConfig.shape;
    pyConfig.specified = cConfig.specified;

    return pyConfig;
  }

}
