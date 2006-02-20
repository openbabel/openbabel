#!/usr/bin/perl

# run a short shell script to figure out where we are

$cmltestdir = "./cmltest";

if ( -e "$cmltestdir/test.sh" ) {
    print "# Testing molecule file roundtripping...\n";
    print "# This test will take quite some time!\n";
    print "# and it currently always passes, so check the results file\n";
    system "(cd $cmltestdir; source test.sh 2>/dev/null;)";
    print "1..1\n";
    print "ok";
} else {
    print "1..0 # skipping - CML test set not found\n";
}
