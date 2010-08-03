use ExtUtils::MakeMaker;

use lib 'inc';
use Devel::CheckLib;

check_lib_or_exit( lib => 'openbabel' , libpath =>
@ENV{"LD_LIBRARY_PATH"});

$CC = 'g++';

$cmake_src_dir = @ENV{"SRC_DIR"};
$cmake_bin_dir = @ENV{"OBJ_DIR"};

$ldfrom = "\$(OBJECT) -lopenbabel -lz";
$ldfrom = "\$(OBJECT) -L$(cmake_bin_dir)/lib -lopenbabel -lz"
    if (-r "$(cmake_bin_dir)/lib/libopenbabel.dylib") and (-s _) and (-B _);
$ldfrom = "\$(OBJECT) -L$(cmake_bin_dir)/lib -lopenbabel -lz"
    if (-r "$(cmake_bin_dir)/lib/libopenbabel.so") and (-s _) and (-B _);

# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    'NAME'	=>  'Chemistry::OpenBabel',
    'AUTHOR'    =>  'Geoffrey Hutchison <openbabel-scripting@lists.sourceforge.net>',
    'VERSION'   => '1.4.0',
    'LDFROM'    => $ldfrom,
    'CC'        => $CC,
    'LD'        => '$(CC)',
    'INC'       => '-I../../include -I'.$cmake_src_dir.'/include',
    'OBJECT'    => 'openbabel-perl.o'
);
