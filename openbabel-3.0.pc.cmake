prefix=@CMAKE_INSTALL_PREFIX@
exec_prefix=${prefix}
libdir=${exec_prefix}/@LIB@
includedir=${prefix}/include
pkgincludedir=${includedir}/openbabel-3.0

Name: Open Babel library
Description: libopenbabel3
Version: @BABEL_VERSION@
Libs: -L${libdir} -lopenbabel3
Cflags: -I${pkgincludedir}
