# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.pl'

######################### We start with some black magic to print on failure.

# Change 1..1 below to 1..last_test_to_print .
# (It may become useful if the test is moved to ./t subdirectory.)

BEGIN { $| = 1; print "1..11\n"; }
END {print "not ok 1\n" unless $loaded;}
use Chemistry::OpenBabel;
$loaded = 1;
print "ok 1\n";

######################### End of black magic.

# Insert your test code below (better if it prints "ok 13"
# (correspondingly "not ok 13") depending on the success of chunk 13
# of the test code):

my $obMol = new Chemistry::OpenBabel::OBMol;
print "ok 2\n";

$numAtoms = $obMol->NumAtoms();
if ($numAtoms == 0) {
	print "ok 3\n";
} else {
	print "not ok 3\n";
}

$obMol->NewAtom();
$numAtoms = $obMol->NumAtoms();
if ($numAtoms == 1) {
	print "ok 4\n";
} else {
	print "not ok 4\n";
}

$obMol->NewAtom();
$obMol->AddBond(1, 2, 1);
$numBonds = $obMol->NumBonds();
if ($numBonds == 1) {
	print "ok 5\n";
} else {
	print "not ok 5\n";
}

my $obConversion = new Chemistry::OpenBabel::OBConversion;
$obConversion->SetInAndOutFormats("smi", "mdl");
print "ok 6\n";

$obMol->Clear();
$obConversion->ReadString($obMol, "C1=CC=CS1");
print "ok 7\n";

$numAtoms = $obMol->NumAtoms();
if ($numAtoms == 5) {
    print "ok 8\n";
} else {
    print "not ok 8\n";
}

$obMol->AddHydrogens();
$numAtoms = $obMol->NumAtoms();
if ($numAtoms == 9) {
    print "ok 9\n";
} else {
    print "not ok 9\n";
}

my $outMDL = $obConversion->WriteString($obMol);
print "$outMDL\n";
print "ok 10\n";

my $mass = $Chemistry::OpenBabel::etab->GetMass(2);
print "mass: $mass\n";
print "ok 11\n";