use ExtUtils::MakeMaker;

$CC = 'g++';
$srcdir = '../../src';

# check if we're compiling in the source directory
$ldfrom = '$(OBJECT) -lopenbabel';
$ldfrom = "\$(OBJECT) $srcdir/.libs/libopenbabel.a"
    if (-r "$srcdir/.libs/libopenbabel.a") and (-s _) and (-B _);

# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    'NAME'	=>  'Chemistry::OpenBabel',
    'AUTHOR'    =>  'Geoffrey Hutchison <openbabel-scripting@lists.sourceforge.net>',
    'VERSION'   => '1.0.0b5',
    'LDFROM'    => $ldfrom,
    'CC'        => $CC,
    'LD'        => '$(CC)',
    'INC'       => '-I../../src',
    'OBJECT'    => 'openbabel_perl.o'
);
