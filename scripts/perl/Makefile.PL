use ExtUtils::MakeMaker;

use lib 'inc';
use Devel::CheckLib;

check_lib_or_exit( lib => 'openbabel' , libpath =>
@ENV{"LD_LIBRARY_PATH"});

$CC = 'g++';
$srcdir = '../../src';

# check if we're compiling in the source directory
# check first for the static, then change to dynamic if available
$ldfrom = "\$(OBJECT) $srcdir/.libs/libopenbabel.a -lz"
    if (-r "$srcdir/.libs/libopenbabel.a") and (-s _) and (-B _);
$ldfrom = "\$(OBJECT) -L$srcdir/.libs -lopenbabel -lz"
    if (-r "$srcdir/.libs/libopenbabel.dylib") and (-s _) and (-B _);
$ldfrom = "\$(OBJECT) -L$srcdir/.libs -lopenbabel -lz"
    if (-r "$srcdir/.libs/libopenbabel.so") and (-s _) and (-B _);

$cmake_src_dir = @ENV{"SRC_DIR"};

# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    'NAME'	=>  'Chemistry::OpenBabel',
    'AUTHOR'    =>  'Geoffrey Hutchison <openbabel-scripting@lists.sourceforge.net>',
    'VERSION'   => '1.3.0',
    'LDFROM'    => $ldfrom,
    'CC'        => $CC,
    'LD'        => '$(CC) -shared',
    'INC'       => '-I../../include -I'.$cmake_src_dir.'/include',
    'OBJECT'    => 'openbabel-perl.o'
);
