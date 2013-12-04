<?php
include_once '/usr/local/lib/openbabel.php';

$obMol = new OBMol;

$obMol->NewAtom();
$numAtoms = $obMol->NumAtoms(); # now 1 atom

$obMol->NewAtom();
$obMol->AddBond(1, 2, 1);
$numBonds = $obMol->NumBonds(); # now 1 bond

$obMol->Clear();

$obConversion = new OBConversion;
$obConversion->SetInAndOutFormats("smi", "svg");
$obConversion->ReadString($obMol, "C1=CC=CS1");

$numAtoms = $obMol->NumAtoms(); # now 5 atoms

$obMol->AddHydrogens();
$numAtoms = $obMol->NumAtoms(); # now 9 atoms

$outMDL = $obConversion->WriteString($obMol);

echo $outMDL
?>
