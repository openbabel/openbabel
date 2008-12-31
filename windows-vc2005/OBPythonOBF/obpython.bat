set VER=1.4

if %1a==a EXIT 1 
set PYTHON_DIR=%1
if %2a==a EXIT 1 
set PYTHON_VER=%2
echo %PYTHON_DIR%\python setup.py bdist_wininst --bitmap=logo.bmp
%PYTHON_DIR%\python setup.py bdist_wininst --bitmap=logo.bmp
echo move dist\openbabel-python-%VER%.win32.exe dist\openbabel-python-%VER%.py%PYTHON_VER%.exe
move dist\openbabel-python-%VER%.win32.exe dist\openbabel-python-%VER%.py%PYTHON_VER%.exe
