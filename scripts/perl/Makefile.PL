use ExtUtils::MakeMaker;

use lib 'inc';
use Devel::CheckLib;

$CC = 'g++';

$cmakesrcdir = @ENV{"SRC_DIR"};
$cmakebindir = @ENV{"OBJ_DIR"};

# Look for libopenbabel.so
if (-e "$cmakebindir/lib/libopenbabel.so" || -e "$cmakebindir/lib/libopenbabel.dylib") {
  check_lib_or_exit( lib => 'openbabel' , libpath => "$cmakebindir/lib");
}
elsif (-e "$cmakebindir/lib64/libopenbabel.so") {
  check_lib_or_exit( lib => 'openbabel' , libpath => "$cmakebindir/lib64");
}
else {
  check_lib_or_exit( lib => 'openbabel' , libpath => @ENV{"LD_LIBRARY_PATH"});
}

$ldfrom = "\$(OBJECT) -lopenbabel";
$ldfrom = "\$(OBJECT) -L$cmakebindir/lib -lopenbabel"
    if (-r "$cmakebindir/lib/libopenbabel.dylib") and (-s _) and (-B _);
$ldfrom = "\$(OBJECT) -L$cmakebindir/lib -lopenbabel"
    if (-r "$cmakebindir/lib/libopenbabel.so") and (-s _) and (-B _);
$ldfrom = "\$(OBJECT) -L$cmakebindir/lib64 -lopenbabel"
    if (-r "$cmakebindir/lib64/libopenbabel.so") and (-s _) and (-B _);

# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    'NAME'	=>  'Chemistry::OpenBabel',
    'AUTHOR'    =>  'Geoffrey Hutchison <openbabel-scripting@lists.sourceforge.net>',
    'VERSION'   => '1.5.0',
    'LDFROM'    => $ldfrom,
    'CC'        => $CC,
    'LD'        => '$(CC)',
    'INC'       => '-I../../include -I'.$cmakesrcdir.'/include',
    'OBJECT'    => 'openbabel-perl.o'
);
