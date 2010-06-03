To create the zlib1 dll and lib, I did the following:

(1) Got zlib124.zip from http://www.winimage.com/zLibDll/index.html
(2) Opened zlib.dsw into MSVC++ 2008
(3) In Configuration manager for zlib, make a x64 target
(4) Make the DLL Release =>  zlib-1.2.4\projects\visualc6\Win32_DLL_Release\zlib1.dll and zlib1.lib

To create libinchi

(1) Got stdinchi-1-api.zip (InChI v 1.0.2) from http://www.iupac.org/inchi/download/index.html
(2) Open STDINCHI-1-API\STDINCHI_API\vc6_INCHI_DLL\vc6_INCHI_DLL.dsw
(3) Change afxres.h to windows.h in INCHI_DLL.rc
(4) In configuration manager, make a x64 target
(5) Creates Release\libstdinchi.lib and .dll

The libxml files were downloaded as binaries from the PHP project:

(1) Got libxml2-2.7.3-vc9-x64.zip from http://pecl2.php.net/downloads/php-windows-builds/php-libs/VC9/x64/
(2) Extracted libxml2.lib and libxml2.dll

libiconv.dll (required for XML) was also downloaded as a binary from the PHP project:

(1) Got libiconv-1.12-vc9-x64.zip from http://pecl2.php.net/downloads/php-windows-builds/php-libs/VC9/x64/
(2) Extracted libiconv.dll

