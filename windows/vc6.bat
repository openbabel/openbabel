@echo off
REM Batch file for installation of OpenBabel with Visual Studio
REM Opens Visual Studio. Click Compile to produce obabel.exe in \src folder.
copy babelconfig.h.vc6 ..\src\babelconfig.h
MSDEV.EXE obabel-vc6.dsw
