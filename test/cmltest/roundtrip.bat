@echo off
REM test OB with CML
REM %1 is inputfile root
REM %2 is inputfile suffix
REM %3 is inputfile type
REM e.g. roundtrip foo mol mdl

set BABEL=..\..\src\babel

@echo on
REM CML1 output
%BABEL% -i%3 %1.%2 -ocml %1.1.cml -x1v
REM CML2 output
%BABEL% -i%3 %1.%2 -ocml %1.2.cml -x2v
REM CML1+array output
%BABEL% -i%3 %1.%2 -ocml %1.a1.cml -xa1v
REM CML2+array output
%BABEL% -i%3 %1.%2 -ocml %1.a2.cml -xa2v

REM roundtrip to MOL; should be identical
%BABEL% -icml %1.1.cml  -o%3 %1.1.%2 -x2v
%BABEL% -icml %1.2.cml  -o%3 %1.2.%2 -x2v
%BABEL% -icml %1.a1.cml -o%3 %1.a1.%2 -x2v
%BABEL% -icml %1.a2.cml -o%3 %1.a2.%2 -x2v



