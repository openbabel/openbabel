#!/usr/bin/perl
# bin2hex.pl
# OpenEye Scientific Software
# September 2001

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
    print "#ifndef OB_" . $guard . "_H\n";
    print "#define OB_" . $guard . "_H\n\n";
    print "namespace OpenBabel\n{\n";
    print "static const char " . $arrayname . "[] = {\n";
}

binmode(F);

$col = 0;
$init = 0;
$ignore = 0;
$newline = 1;

while( !eof(F) ) {
    $ch = ord(getc(F));
    if( $ch == 13 ) {
        $ch = 0;
    }

    if( $ignore ) {
        if( $ch == 10 ) {
            $ignore = 0;
        }
        $ch = 0;
    } elsif( $newline ) {
        if( $ch == 10 ) {
            $ch = 0;
        } elsif ( $ch == 35 ) {
            $ignore = 1;
            $ch = 0;
        } elsif( $ch == 32 ) {
            $ch = 0;
        } elsif( $ch ) {
            $newline = 0;
        }
    } elsif( $ch == 10 ) {
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
                print "\n";
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

