set ver=0.1
%CSHARP%/Csc.exe /target:library /optimize /out:OBDotNet.dll *.cs
del /Q *.cs
set dist=OBDotNet-%ver%
rmdir /s /q %dist%
mkdir %dist%
copy OBDotNet.dll %dist%
copy openbabel.dll %dist%
copy IronPython_Instructions.txt %dist%
copy ..\*.obf %dist%
copy ..\*.dll %dist%
