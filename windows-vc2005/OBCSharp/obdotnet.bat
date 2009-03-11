set ver=0.1
set dist=..\..\scripts\csharp\OBDotNet\bin\release
mkdir %dist%\data
copy ..\..\data %dist%\data
copy openbabelcsharp.dll %dist%
copy ..\*.obf %dist%
copy ..\*.dll %dist%
del ..\..\scripts\csharp\OBDotNet\*.cs
copy *.cs ..\..\scripts\csharp\OBDotNet
del *.cs