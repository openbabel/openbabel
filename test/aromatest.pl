#!/usr/bin/perl

system("./aromatic aromatics.smi");

# If the program failed to execute, ignore the test
if ($? == -1) {
    print "1..0 skip because program would not run";
}

# Otherwise exit
exit($? >>8);
