#!/usr/bin/perl
# test Open Babel with roundtripping if the test suite exists
#   (just un-tar into this folder as "test-set" and it'll run)

if ( -d "test-set" ) {
    print "# Testing molecule file roundtripping...\n";
    print "# This test will take quite some time!\n";
    print "# and it currently always passes, so check the results file\n";
    system "(cd test-set; source test.sh;)";
} else {
    print "1..0 # skipping - roundtrip test set not found\n";
}
