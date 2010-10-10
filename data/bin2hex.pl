#!/usr/bin/env perl
# bin2hex.pl

########################################################################
# Copyright (C) 2001 by OpenEye Scientific Software, Inc.
# Some portions Copyright (c) 2002 by Geoffrey R. Hutchison
#
# This file is part of the Open Babel project.
# For more information, see <http://openbabel.org/>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation version 2 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
########################################################################


# Autoflush STDOUT
#STDOUT->autoflush(1);
$| = 1;

$argc = @ARGV;
if( $argc == 2 ) {
    $filename = @ARGV[0];
    $arrayname = @ARGV[1];

    $pos = index($filename,".");
    if( $pos > 0 ) {
        $guard = uc(substr($filename,0,$pos));
    } else {
        $guard = uc($filename);
    }
} elsif( $argc == 3 ) {
    $filename = @ARGV[0];
    $arrayname = @ARGV[1];
    $guard = @ARGV[2];
} else {
    print "usage:  bin2hex.pl <binaryfile> <arrayname>\n\n";
    exit;
}

$debug = 0;

open(F,$filename) || die "Error: Unable to open binary file!\n";

if( !$debug ) {
    print "/***************************************************************\n";
    print "This file is part of the Open Babel project.\n";
    print "This is copyright under the GNU General Public License (GPL)\n";
    print "For more information, see <http://openbabel.org/>\n\n";
    print "This file contains a binary representation of data tables\n";
    print " used by Open Babel. It is used as a fallback if the textual\n";
    print " data table is not found at runtime. It is not a normal header.\n";
    print "***************************************************************/\n";
    print "\n\n";

    print "#ifndef OB_" . $guard . "_H\n";
    print "#define OB_" . $guard . "_H\n\n";
    print "namespace OpenBabel\n{\n";
    print "static const char " . $arrayname . "[] = {\n ";
}

binmode(F);

$col = 0;
$init = 0;
$ignore = 0;
$newline = 1;
$spacerun = 0; # collapse runs of spaces

while( !eof(F) ) {
    $ch = ord(getc(F));
    if( $ch == 13 ) { # ignore \r characters
        $ch = 0;
    }

    if ( $spacerun ) {
        if ( $ch == 32 ) {
            $ch = 0;
        } else {
            $spacerun = 0;
        }
    } elsif ( $ch == 32) {
        $spacerun = 1;
    }

    if( $ignore ) {
        if( $ch == 10 ) { # found the \n after an '#' ignore comment
            $ignore = 0;
        }
        $ch = 0;
    } elsif( $newline ) { # just saw a \n -- do we see real text
        if( $ch == 10 ) {
            $ch = 0;
        } elsif ( $ch == 35 ) { # ignore anything after a '#' until a \n
            $ignore = 1;
            $ch = 0;
        } elsif( $ch == 32 ) { # space
            $ch = 0;
        } elsif( $ch == 9 ) { # tab
            $ch = 0;
        } elsif( $ch ) { # something else, so clear the blank-line boolean
            $newline = 0;
        }
    } elsif( $ch == 10 ) { # set the blank-line detector
        $newline = 1;
    }

    if( $ch ) {
        if( $debug ) {
            print chr($ch);
        } else {
            if( $init ) {
                print ",";
            } else {
                $init = 1;
            }
            if( $col >= 15 ) {
                print "\n ";
                $col = 0;
            }
            print sprintf("0x%02X",$ch);
            $col++;
        }
    }
}

if( !$debug ) {
    if( $col >= 15 ) {
        print ",\n0x00};\n\n";
    } else {
        print ",0x00};\n\n";
    }
    print "} // namespace OpenBabel\n";
    print "#endif // OB_" . $guard . "_H\n\n";
}

close(F);
exit;

