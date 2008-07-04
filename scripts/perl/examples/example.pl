#!/usr/bin/perl

use Chemistry::OpenBabel;

my $obMol = new Chemistry::OpenBabel::OBMol;

$obMol->NewAtom();
$numAtoms = $obMol->NumAtoms(); # now 1 atom

$obMol->NewAtom();
$obMol->AddBond(1, 2, 1);
$numBonds = $obMol->NumBonds(); # now 1 bond

$obMol->Clear();

my $obConversion = new Chemistry::OpenBabel::OBConversion;
$obConversion->SetInAndOutFormats("smi", "mdl");
$obConversion->ReadString($obMol, "C1=CC=CS1");

$numAtoms = $obMol->NumAtoms(); # now 5 atoms

$obMol->AddHydrogens();
$numAtoms = $obMol->NumAtoms(); # now 9 atoms

my $outMDL = $obConversion->WriteString($obMol);
