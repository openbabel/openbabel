@echo off
rem Run coverage.py under Windows
rem Creates a coverage report and the annotated file "pybel.py,cover"

set pythondir=C:\Python25
set coverage=D:\Tools\coverage\coverage-2.78\coverage.py

%pythondir%\python %coverage% -e
%pythondir%\python %coverage% -x testpybel.py
%pythondir%\python %coverage% -r %pythondir%\Lib\site-packages\pybel.py
%pythondir%\python %coverage% -a -d . %pythondir%\Lib\site-packages\pybel.py

echo.
echo Check pybel.py,cover for lines starting with !
