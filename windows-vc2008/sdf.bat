@rem Drop chemical files on to this file in Windows Explorer to convert
@rem to the format in the name of this file. For instance, if this file
@rem is called sdf.bat, an sdf file will be produced. 
@rem You can copy and rename the file to any OpenBabel format.
@babel.exe %1 "%~ndp1.%~n0"
