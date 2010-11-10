@echo off
echo b_obpython.bat PYTHON_DIR PYTHON_VER

set VER=1.6

if %1a==a GOTO:EOF
set PYTHON_DIR=%1
if %2a==a GOTO:EOF
set PYTHON_VER=%2
echo %PYTHON_DIR%\python setup.py bdist_wininst --bitmap=logo.bmp
%PYTHON_DIR%\python setup.py bdist_wininst --bitmap=logo.bmp
echo move dist\openbabel-python-%VER%.win32.exe dist\openbabel-python-%VER%.py%PYTHON_VER%.exe
move dist\openbabel-python-%VER%.win32.exe dist\openbabel-python-%VER%.py%PYTHON_VER%.exe
