dnl @synopsis AC_DEFINE_DIR(VARNAME, DIR [, DESCRIPTION])
dnl
dnl This macro _AC_DEFINEs VARNAME to the expansion of the DIR
dnl variable, taking care of fixing up ${prefix} and such.
dnl
dnl Note that the 3 argument form is only supported with autoconf 2.13 and
dnl later (i.e. only where _AC_DEFINE supports 3 arguments).
dnl
dnl Examples:
dnl
dnl    AC_DEFINE_DIR(DATADIR, datadir)
dnl    AC_DEFINE_DIR(PROG_PATH, bindir, [Location of installed binaries])
dnl
dnl @version $Id$
dnl @author Guido Draheim <guidod@gmx.de>, original by Alexandre Oliva

AC_DEFUN([AC_DEFINE_DIR], [
  test "x$prefix" = xNONE && prefix="$ac_default_prefix"
  test "x$exec_prefix" = xNONE && exec_prefix='${prefix}'
  ac_define_dir=`eval echo [$]$2`
  ac_define_dir=`eval echo [$]ac_define_dir`
  ifelse($3, ,
    AC_DEFINE_UNQUOTED($1, "$ac_define_dir"),
    AC_DEFINE_UNQUOTED($1, "$ac_define_dir", $3))
])
