# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.pl'

######################### We start with some black magic to print on failure.

# Change 1..1 below to 1..last_test_to_print .
# (It may become useful if the test is moved to ./t subdirectory.)

BEGIN { $| = 1; print "1..9\n"; }
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

my $obConversion = new Chemistry::OpenBabel::OBConversion;
$obConversion->SetInAndOutFormats("smi", "mdl");
print "ok 3\n";

$obConversion->ReadString($obMol, "C1=CC=CS1");
print "ok 4\n";

if ($obMol->NumAtoms() == 5) {
    print "ok 5\n";
} else {
    print "not ok 5\n";
}

$obMol->AddHydrogens();
if ($obMol->NumAtoms() == 9) {
    print "ok 6\n";
} else {
    print "not ok 6\n";
}

my $outMDL = $obConversion->WriteString($obMol);
print "ok 7\n";

$obConversion->WriteFile($obMol, "test.mdl");
if (-e "test.mdl") {
    print "ok 8\n"
} else {
    print "not ok 8\n";
}

$obConversion->SetInAndOutFormats("mdl", "mdl");
$obConversion->ReadFile($obMol, "test.mdl");
unlink "test.mdl";
print "ok 9\n";

# RegisterFormat
# FindFormat
# FormatFromExt
# FormatFromMIME
# GetNextFormat
