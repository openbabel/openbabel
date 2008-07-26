# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.pl'

######################### We start with some black magic to print on failure.

# Change 1..1 below to 1..last_test_to_print .
# (It may become useful if the test is moved to ./t subdirectory.)

BEGIN { $| = 1; print "1..7\n"; }
END {print "not ok 1\n" unless $loaded;}
use Chemistry::OpenBabel;
$loaded = 1;
print "ok 1\n";

######################### End of black magic.

# Insert your test code below (better if it prints "ok 13"
# (correspondingly "not ok 13") depending on the success of chunk 13
# of the test code):


###
### OBMol tests
###

my $emptyMol = new Chemistry::OpenBabel::OBMol;
my $testMol1 = new Chemistry::OpenBabel::OBMol;
print "ok 2\n";

$testMol1->ReserveAtoms(-1);
$testMol1->ReserveAtoms(0);
$testMol1->ReserveAtoms(2);
print "ok 3\n";

## atom component tests

if ($testMol1->NumAtoms() == 0) {
	print "ok 4\n";
} else {
	print "not ok 4\n";
}

$testMol1->NewAtom();
if ($testMol1->NumAtoms() == 1) {
	print "ok 5\n";
} else {
	print "not ok 5\n";
}

$testMol1->NewAtom();
$testMol1->AddBond(1, 2, 1);
if ($testMol1->NumBonds() == 1) {
	print "ok 6\n";
} else {
	print "not ok 6\n";
}

$testMol1->Clear();
if ($testMol1->NumAtoms() == 0) {
	print "ok 7\n";
} else {
	print "not ok 7\n";
}
