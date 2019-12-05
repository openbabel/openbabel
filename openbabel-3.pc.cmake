prefix=@CMAKE_INSTALL_PREFIX@
exec_prefix=${prefix}
libdir=${exec_prefix}/@LIB_INSTALL_DIR@
includedir=${prefix}/include
pkgincludedir=${includedir}/openbabel@BABEL_MAJ_VER@

Name: Open Babel library
Description: libopenbabel
Version: @BABEL_VERSION@
Libs: -L${libdir} -lopenbabel
Cflags: -I${pkgincludedir}
