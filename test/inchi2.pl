#!/usr/bin/perl
use Env qw(TESTDATADIR);

if (defined $TESTDATADIR) {
  $INCHIDIR="$TESTDATADIR/../inchi";
} else {
  $INCHIDIR="./inchi";
}

system("./inchiwrite $INCHIDIR/SamplesTechMan.sdf $INCHIDIR/SamplesTechMan.txt 2>/dev/null");

# If the program failed to execute, ignore the test
if ($? == -1) {
    print "1..0 skip because program would not run";
}

# Otherwise exit
exit($? >>8);
