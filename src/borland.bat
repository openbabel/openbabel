@REM compile with Borland free commandline compiler
@REM edit the directory
@SET BORLAND=C:\borland
 
REM %BORLAND%\bcc55\bin\bcc32 -DDATADIR=\".\" -w-8004 -w-8008 -w-8012 -w-8022 -w-8026 -w-8027 -w-8060 -w-8066 -w-8057 -I%BORLAND%\bcc55\include -L%BORLAND%\bcc55\Lib -ebabel.exe *.cpp math\*.cpp

%BORLAND%\bcc55\bin\bcc32 -DDATADIR=\".\" -w-8004 -w-8022 -w-8026 -w-8027  -w-8057 -I%BORLAND%\bcc55\include -L%BORLAND%\bcc55\Lib -ebabel.exe *.cpp math\*.cpp



