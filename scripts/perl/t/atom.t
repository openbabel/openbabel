# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.pl'

######################### We start with some black magic to print on failure.

# Change 1..1 below to 1..last_test_to_print .
# (It may become useful if the test is moved to ./t subdirectory.)

BEGIN { $| = 1; print "1..6\n"; }
END {print "not ok 1\n" unless $loaded;}
use Chemistry::OpenBabel;
$loaded = 1;
print "ok 1\n";

######################### End of black magic.

# Insert your test code below (better if it prints "ok 13"
# (correspondingly "not ok 13") depending on the success of chunk 13
# of the test code):


###
### OBAtom isolation tests (no connection to residue, bond, molecule...)
###

my $emptyAtom = new Chemistry::OpenBabel::OBAtom;
my $testAtom1 = new Chemistry::OpenBabel::OBAtom;
my $testAtom2 = new Chemistry::OpenBabel::OBAtom;
print "ok 2\n";

$testAtom1->SetIdx(0);
print $testAtom1->GetIdx(), "\n";
$testAtom1->SetIdx(-1);
print $testAtom1->GetIdx(), "\n";
$testAtom1->SetIdx(1);
print "ok 3\n";

$testAtom1->SetAtomicNum(0);
print $testAtom1->GetAtomicNum(), "\n";
$testAtom1->SetAtomicNum(-1);
print $testAtom1->GetAtomicNum(), "\n";
$testAtom1->SetAtomicNum(200);
print $testAtom1->GetAtomicNum(), "\n";
$testAtom1->SetAtomicNum(300);
print $testAtom1->GetAtomicNum(), "\n";
$testAtom1->SetAtomicNum(1);
print "ok 4\n";

$coordPtr = $testAtom1->GetCoordinate();
print "ok 5\n";

$testAtom1->SetCoordPtr($coordPtr);
print "ok 6\n"
