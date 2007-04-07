@REM test OB with CML

@set BABEL=..\..\src\babel

@REM test input
REM CML2 with array
%BABEL% -icml cs2a.cml -omdl cs2a.mol 
REM 3D molecules in SDF 
REM CML2 with XML version
%BABEL% -isdf 3d.head.sdf -ocml 3d.head.2.cml -x2v 
REM CML1 with DOCTYPE
%BABEL% -isdf cs2a.mol -ocml cs2a.mol.cml -x1d
REM CML2 arrays with namespaces (large)
%BABEL% -isdf 3d.head.sdf -ocml 3d.head.2an.cml -x2an 

@REM roundtripping; arguments are fileroot; input format; input suffix
REM 2d MDL to CML and back again through all main variants
@call roundtrip.bat nsc2dmol mol mdl
REM 3d MDL to CML and back again through all main variants
@call roundtrip.bat nsc3dmol mol mdl



