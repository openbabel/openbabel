use ExtUtils::MakeMaker;

$CC = 'g++';

# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    'NAME'	=>  'Chemistry::OpenBabel',
    'AUTHOR'    =>  'Geoffrey Hutchison <openbabel-scripting@lists.sourceforge.net>',
    'VERSION'   => '2.0.0b4',
    'LIBS'      => '-lopenbabel',
    'CC'        => $CC,
    'LD'        => '$(CC)',
    'INC'       => '-I../../src',
    'OBJECT'    =>  'openbabel_perl.o'
);
