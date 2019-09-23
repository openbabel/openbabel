prefix=@CMAKE_INSTALL_PREFIX@
exec_prefix=${prefix}
libdir=${exec_prefix}/@LIB_INSTALL_DIR@
includedir=${prefix}/include
pkgincludedir=${includedir}/openbabel-2.0

Name: Open Babel library
Description: libopenbabel
Version: @BABEL_VERSION@
Libs: -L${libdir} -lopenbabel
Cflags: -I${pkgincludedir}
