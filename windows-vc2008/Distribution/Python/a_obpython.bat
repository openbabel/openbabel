@echo off
echo a_obpython.bat PYTHON_DIR PYTHON_VER

if %1a==a GOTO:EOF
set PYTHON_DIR=%1
if %2a==a GOTO:EOF
set PYTHON_VER=%2


cd ..\..
call default_build.bat "-DPYTHON_EXECUTABLE=%PYTHON_DIR%\python.exe" "-DPYTHON_INCLUDE_DIR=%PYTHON_DIR%\include" "-DPYTHON_INCLUDE_PATH=%PYTHON_DIR%\include" "-DPYTHON_LIBRARY=%PYTHON_DIR%\libs\python%PYTHON_VER%.lib"
cd Distribution\Python

echo. 
echo ****************************************************************
echo Now open openbabel.sln, touch the openbabel-python.i and rebuild 
echo the Python subproject
echo ****************************************************************
