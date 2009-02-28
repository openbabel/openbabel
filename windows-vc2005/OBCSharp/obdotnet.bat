set ver=0.1
copy OBDotNetAssemblyInfo.txt OBDotNetAssemblyInfo.cs
%CSHARP%/Csc.exe /target:library /keyfile:obdotnet.snk /optimize /out:OBDotNet.dll *.cs
del /Q *.cs
set dist=OBDotNet-%ver%
rmdir /s /q %dist%
mkdir %dist%
mkdir %dist%\data
copy ..\..\data %dist%\data
copy OBDotNet.dll %dist%
copy openbabelcsharp.dll %dist%
copy IronPython_Instructions.txt %dist%
copy ..\*.obf %dist%
copy ..\*.dll %dist%
