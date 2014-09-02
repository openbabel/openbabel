<?php
/********************************************************************
baphpel - PHP interface to OpenBabel
Copyright (C) 2013  Maciej Wojcikowski <maciek@wojcikowski.pl>

This file is part of  of the Open Babel project.

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

********************************************************************/
namespace baphpel;

include_once 'openbabel.php';

# define global elements
$_builder = new \OBBuilder;
$_obconv = new \OBConversion;

function vector2array($v) {
	$out = array();
	$num = $v -> size();
	for($i=0;$i<$num;$i++) {
		$out[] = $v -> get($i);
	}
	return $out;
}

function _formatstodict($list) {
	global $_obconv;
	$list = vector2array($_obconv -> $list());
	foreach($list as $f) {
		$fmt = explode(' -- ', preg_replace(array('/\[Read-only\]/'), '/\[Write-only\]/', $f));
		$out[trim($fmt[0])] = trim($fmt[1]);
	}
	return $out;
}

$informats = _formatstodict('GetSupportedInputFormat');
$outformats = _formatstodict('GetSupportedOutputFormat');

function _getplugins($findplugins, $names) {
	$find = explode('::', $findplugins);
	foreach($names as $x) {
		if(in_array($find[0], array('OBDescriptor', 'OBForceField'))) {
			$f = $find[0]::$find[1]($x);
		}
		else {
			$f1 = new $find[0](null);
			$f = $f1 -> $find[1]($x);
		}
		if($f) {
			$plugins[$x] = $f;
		}
	}
	return $plugins;
}

function _getpluginnames($ptype) {
	$plugins = new \vectorString;
	\OBPlugin::ListAsVector($ptype, null, $plugins);
	$plugins = vector2array($plugins);
	foreach($plugins as $p) {
		$out[] = strtolower(explode(' ', trim($p))[0]);
	}
	return $out;
}


/* A list of supported descriptors */
$descs = _getpluginnames('descriptors');
$_descdict = _getplugins('\OBDescriptor::FindType', $descs);

/* A list of supported forcefields */
$forcefields = _getpluginnames('forcefields');
$_forcefields = _getplugins('\OBForceField::FindForceField', $forcefields);

/* A list of supported fingerprint types */
$fps = _getpluginnames('fingerprints');
$_fingerprinters = _getplugins('\OBFingerprint::FindType', $fps);

/* A list of supported operations */
$operations = _getpluginnames('ops');
$_operations = _getplugins('\OBOp::FindType', $operations);

function readfile($format, $filename, $opt=array()) {
	# setup converter
	$OBConversion = new \OBConversion;
	$OBConversion -> SetInFormat($format);
	# set options
	if(!empty($opt)) {
		foreach($opt as $k => $v) {
			if(!empty($v)) {
				$OBConversion -> AddOption($k, $OBConversion::INOPTIONS, $v);
			}
			else {
				$OBConversion -> AddOption($k, $OBConversion::INOPTIONS);
			}
		}
	}
	
	class filereader implements \Iterator {
		private $position = 0;  
		private $notatend = null;
		private $OBMol = null;
		
		public function __construct($OBConversion, $filename) {
			$this -> OBConversion = $OBConversion;
			$this -> OBMol = new \OBMol;
			$this -> notatend = $this -> OBConversion -> ReadFile($this -> OBMol, $filename);
		}

		function rewind() {
			$this->position = 0;
		}

		function current() {
			return new Molecule($this -> OBMol);
		}

		function key() {
			return $this->position;
		}

		function next() {
			$this -> position++;
			$this -> OBMol = new \OBMol;
			$this -> notatend = $this -> OBConversion -> Read($this -> OBMol);
		}		

		function valid() {
			return $this -> notatend;
		}
	}
	
	return new filereader($OBConversion, $filename);
}

function readstring($format, $string, $opt=array()) {
	# setup converter
	$OBConversion = new \OBConversion;
	$OBConversion -> SetInFormat($format);
	# set options
	if(!empty($opt)) {
		foreach($opt as $k => $v) {
			if(!empty($v)) {
				$OBConversion -> AddOption($k, $OBConversion::INOPTIONS, $v);
			}
			else {
				$OBConversion -> AddOption($k, $OBConversion::INOPTIONS);
			}
		}
	}
	
	$OBMol = new \OBMol;
	$OBConversion->ReadString($OBMol, $string);
	return new Molecule($OBMol);
}

class Outputfile {
	public function __construct($format, $filename, $overwrite=False, $opt=array()) {
		$this -> format = $format;
		$this -> filename = $filename;
		if(!$overwrite && file_exists(self.filename)) {
			throw new \Exception($this -> filename." already exists. Use 'overwrite=True' to overwrite it.");
		}
		$this -> OBConversion = new \OBConversion;
		$formatok = $this -> OBConversion -> SetOutFormat($this -> format);
		if(!$formatok) {
			throw new \Exception($this -> format." is not a recognised Open Babel format");
		}
		foreach($opt as $key => $value) {
			if($v === null) {
				$this -> OBConversion -> AddOption($k, $this -> OBConversion -> OUTOPTIONS);
			}
			else {
				$this -> OBConversion -> AddOption($k, $this -> OBConversion -> OUTOPTIONS, $v);
			}
		}
		$this -> total = 0; # The total number of molecules written to the file
	}
	
	public function write($molecule) {
		if(!$this -> filename) {
			throw new \Exception("Outputfile instance is closed.");
		}

		if($this -> total == 0) {
			$this -> OBConversion -> WriteFile($molecule -> OBMol, $this -> filename);
		}
		else {
			$this -> OBConversion -> Write($molecule -> OBMol);
		}
		$this -> total += 1;
	}

	public function close() {
		$this -> OBConversion -> CloseOutFile();
		$this -> filename = null;
	}
}

class Molecule implements \Iterator {
	public $OBMol = null;
	private $position = 0;
	private $atoms = array();
	private $residues = array();
	
	public function __construct($OBMol) {
		$this -> OBMol = $OBMol;
	}
	
	public function __get($name) {
		switch($name) {
			case 'atoms':
				if(empty($this -> atoms)) { # generate atoms only once
					$out = array();
					$num = $this -> OBMol -> NumAtoms();
					for($i=1;$i<=$num;$i++) {
						$atom = $this -> OBMol -> GetAtom($i);
						if(empty($atom)) {
							echo $i.'!!!';
						}
						$out[] = new Atom($this -> OBMol -> GetAtom($i));
					}
					$this -> atoms = $out;
				}
				else {
					$out = $this -> atoms;
				}
			break;
			case 'residues':
				if(empty($this -> residues)) { # generate residues only once
					$out = array();
					$num = $this -> OBMol -> NumResidues();
					for($i=0;$i<$num;$i++) {
						$out[] = new Residue($this -> OBMol -> GetResidue($i));
					}
					$this -> residues = $out;
				}
				else {
					$out = $this -> residues;
				}
			break;
			case 'charge':
				$out = $this -> OBMol -> GetTotalCharge();
			break;
			case 'conformers':
				$out = $this -> OBMol -> GetConformers();
			break;
			case 'data':
				$out = new MoleculeData($this -> OBMol);
			break;
			case 'dim':
				$out = $this -> OBMol -> GetDimension();
			break;
			case 'energy':
				$out = $this -> OBMol -> GetEnergy();
			break;
			case 'exactmass':
				$out = $this -> OBMol -> GetExactMass();
			break;
			case 'formula':
				$out = $this -> OBMol -> GetFormula();
			break;
			case 'molwt':
				$out = $this -> OBMol -> GetMolWt();
			break;
			case 'spin':
				$out = $this -> OBMol -> GetTotalSpinMultiplicity();
			break;
			case 'sssr':
				$out = vector2array($this -> OBMol -> GetSSSR());
			break;
			case 'title':
				$out = $this -> OBMol -> GetTitle();
			break;
#			case 'unitcell':
#			
#			break;
			case '_exchange':
				if($this -> OBMol -> HasNonZeroCoords()) {
					return $this -> write('mol');	
				}
				else {
					return preg_split('/[\t\s]+/', trim($this -> write('can')))[0];
				}
			break;
		}
		return $out;
	}
	
	public function __set($name, $value) {
		switch($name) {
			case 'title':
				$out = $this -> OBMol -> SetTitle($value);
			break;
		}
		return $out;
	}
	
	# iterator part
	public function rewind() {
		$this->position = 0;
		# generate atoms
		$this-> __get('atoms');
	}

	public function current() {
		return $this->atoms[$this->position];
	}

	public function key() {
		return $this->position;
	}

	public function next() {
		$this->position++;
	}

	public function valid() {
		return isset($this->atoms[$this->position]);
	}
	# end iterator
	
	public function calcdesc($descnames = array()) {
		global $descs, $_descdict;
		if(empty($descnames)) {
			$descnames = $descs;
		}
		foreach($descnames as $descname) {
			if(isset($_descdict[$descname])) {
				$out[$descname] = $_descdict[$descname] -> Predict($this -> OBMol);
			}
			else {
				throw new \Exception($descname.' is not recognised Open Babel descriptor type');
			}
		}
		return $out;
	}

	public function calcfp($fptype="FP2") {
		global $_fingerprinters;
		$fp = new \vectorUnsignedInt;
		$fptype = strtolower($fptype);
		if(isset($_fingerprinters[$fptype])) {
			$fingerprinter = $_fingerprinters[$fptype];
			$fingerprinter -> GetFingerprint($this -> OBMol, $fp);
			return new Fingerprint($fp);
		}
		else {
			throw new \Exception($fptype.' is not a recognised Open Babel Fingerprint type');
		}
	}
	
	public function __toString() {
		return $this -> write();
	}
	
	public function write($format = 'smi', $filename = null, $overwrite = false, $opt = array()) {
		$OBConversion = new \OBConversion;
		$OBConversion -> SetOutFormat($format);
		# set options
		if(!empty($opt)) {
			foreach($opt as $k => $v) {
				if($v === null) {
					$OBConversion -> AddOption($k, $OBConversion::OUTOPTIONS, $v);
				}
				else {
					$OBConversion -> AddOption($k, $OBConversion::OUTOPTIONS);
				}
			}
		}
		
		if(!empty($filename)) {
			if(!file_exists($filename) || file_exists($filename) && $overwrite) {
				return $OBConversion->WriteFile($this -> OBMol, $filename);
			}
		}
		else {
			return $OBConversion->WriteString($this -> OBMol);
		}
	}

	public function localopt($forcefield = 'mmff94', $steps = 50) {
		global $_forcefields;
		$forcefield = strtolower($forcefield);
		if($this -> dim != 3) {
			$this -> make3D($forcefield);
		}
		$ff = $_forcefields[$forcefield];
		if(!$ff -> Setup($this -> OBMol)) {
			return false;
		}
		$ff -> SteepestDescent($steps);
		$ff -> GetCoordinates($this -> OBMol);
	}

	public function make3D($forcefield = 'mmff94', $steps = 50) {
		global $_builder;
		$forcefield = strtolower($forcefield);
		$_builder -> Build($this -> OBMol);
		$this -> addh();
		$this -> localopt($forcefield, $steps);
	}
	
	public function addh() {
		$this -> OBMol -> AddHydrogens();
	}
	
	public function removeh() {
		$this -> OBMol -> DeleteHydrogens();
	}
	
	public function convertdbonds() {
		$this -> OBMol -> ConvertDativeBonds();
	}
	
	public function draw() {
		return $this -> write('svg');
	}
}

class Atom {
	public function __construct($OBAtom) {
		$this -> OBAtom = $OBAtom;
	}
	
	public function __get($name) {
		switch($name) {
			case 'coords':
				$out = array($this -> OBAtom -> GetX(), $this -> OBAtom -> GetY(), $this -> OBAtom -> GetZ());
			break;
			case 'atomicmass':
				$out = $this -> OBAtom -> GetAtomicMass();
			break;
			case 'atomicnum':
				$out = $this -> OBAtom -> GetAtomicNum();
			break;
			case 'cidx':
				$out = $this -> OBAtom -> GetCIdx();
			break;
			case 'coordidx':
				$out = $this -> OBAtom -> GetCoordinateIdx();
			break;
			case 'exactmass':
				$out = $this -> OBAtom -> GetExactMass();
			break;
			case 'formalcharge':
				$out = $this -> OBAtom -> GetFormalCharge();
			break;
			case 'heavyvalence':
				$out = $this -> OBAtom -> GetHvyValence();
			break;
			case 'heterovalence':
				$out = $this -> OBAtom -> GetHeteroValence();
			break;
			case 'hyb':
				$out = $this -> OBAtom -> GetHyb();
			break;
			case 'idx':
				$out = $this -> OBAtom -> GetIdx();
			break;
			case 'implicitvalence':
				$out = $this -> OBAtom -> GetImplicitValence();
			break;
			case 'isotope':
				$out = $this -> OBAtom -> GetIsotope();
			break;
			case 'partialcharge':
				$out = $this -> OBAtom -> GetPartialCharge();
			break;
			case 'residue':
				$out = new Residue($this -> OBAtom -> GetResidue());
			break;
			case 'spin':
				$out = $this -> OBAtom -> GetSpinMultiplicity();
			break;
			case 'type':
				$out = $this -> OBAtom -> GetType();
			break;
			case 'valence':
				$out = $this -> OBAtom -> GetValence();
			break;
			case 'vector':
				$out = $this -> OBAtom -> GetVector();
			break;
		}
		return $out;
	}
	
	public function __toString() {
		
		return 'Atom: '.$this -> atomicnum.' ('.$this -> coords[0] .', '.$this -> coords[1] .', '.$this -> coords[2] .')';
	}
}

class Residue {
	private $atoms = array();
	
	public function __construct($OBResidue) {
		$this -> OBResidue = $OBResidue;
	}
	
	public function __get($name) {
		switch($name) {
			case 'idx':
				$out = $this -> OBResidue -> GetIdx();
			break;
			case 'num':
				$out = $this -> OBResidue -> GetNum();
			break;
			case 'name':
				$out = $this -> OBResidue -> GetName();
			break;
			case 'atoms':
				if(empty($this -> atoms)) { # generate atoms only once
					$num = $this -> OBResidue -> GetNumAtoms();
					$iter = $this -> OBResidue -> BeginAtoms();
					$out[] = new Atom($this -> OBResidue -> BeginAtom($iter));
					for($i=0;$i<$num-1;$i++) {
						$out[] = new Atom($this -> OBResidue -> NextAtom($iter));
					}
					$this -> atoms = $out;
				}
				else {
					$out = $this -> atoms;
				}
			break;
		}
		return $out;
	}
}

function tanimoto($fp1, $fp2) {
	return \OBFingerprint::Tanimoto($fp1 -> fp, $fp2 -> fp);
}

class Fingerprint {
	public $fp = null;
	
	public function __construct($fp) {
		$this -> fp = $fp;
	}
	
	public function __toString() {
		return implode(', ', vector2array($this -> fp));
	}
}

class Smarts {
	private $obsmarts = null;
	
	public function __construct($smartspattern) {
		$this -> obsmarts = new \OBSmartsPattern;
		if(!$this -> obsmarts -> Init($smartspattern)) {
			throw new Exception('Invalid SMARTS pattern');
		}
	}
	
	public function findall($mol) {
		$this -> obsmarts -> Match($mol -> OBMol);
		foreach(vector2array($this -> obsmarts -> GetUMapList()) as $match) {
			$out[] = vector2array($match);
		}
		return $out;
	}
}

class MoleculeData {
	private $_mol = null;
	
	public function __construct($obmol) {
		$this -> _mol = $obmol;
	}
	
	public function _data() {
		$data = vector2array($this -> _mol -> GetData());
		foreach($data as $d) {
			if(in_array($d -> GetDataType(), array(\openbabel::PairData, \openbabel::CommentData))) {
				$out[] = \openbabel::toPairData($d); 
				
			}
		}
		return $out;
	}
	
	public function keys() {
		foreach($this -> _data() as $d) {
			$out[] = $d -> GetAttribute();
		}
		return $out;
	}
	
	public function values() {
		foreach($this -> _data() as $d) {
			$out[] = $d -> GetValue();
		}
		return $out;
	}
	
	public function items() {
		foreach($this -> _data() as $d) {
			$out[$d -> GetAttribute()] = $d -> GetValue();
		}
		return $out;
	}
}
?>
