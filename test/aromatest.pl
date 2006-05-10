#!/usr/bin/perl
use Env qw(TESTDATADIR);

if (defined $TESTDATADIR) {
    system("./aromatic $TESTDATADIR/aromatics.smi");
} else {
    system("./aromatic files/aromatics.smi");
}

# If the program failed to execute, ignore the test
if ($? == -1) {
    print "1..0 skip because program would not run";
}

# Otherwise exit
exit($? >>8);
