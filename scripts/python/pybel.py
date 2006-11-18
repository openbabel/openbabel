import openbabel as ob
import os.path

def readfile(format, filename):
    """Iterate over the molecules in a file.

    Required parameters:
       format
       filename

    You can access the first molecule in a file using:
        mol = readfile("smi", "myfile.smi").next()
        
    You can make a list of the molecules in a file using:
        mols = [mol for mol in readfile("smi", "myfile.smi")]
        
    You can iterate over the molecules in a file as shown in the
    following code snippet...

    >>> atomtotal = 0
    >>> for mol in readfile("sdf","head.sdf"):
    ...     atomtotal += len(mol.atoms)
    ...
    >>> print atomtotal
    43
    """
    obconversion = ob.OBConversion()
    formatok = obconversion.SetInFormat(format)
    if not formatok:
        raise ValueError,"%s is not a recognised OpenBabel format" % format

    obmol = ob.OBMol()
    notatend = obconversion.ReadFile(obmol,filename)
    while notatend:
        yield Molecule(obmol)
        obmol = ob.OBMol()
        notatend = obconversion.Read(obmol)

def readstring(format, string):
    """Read in a molecule from a string.

    Required parameters:
       format
       string

    >>> input = "C1=CC=CS1"
    >>> mymol = readstring("smi",input)
    >>> len(mymol.atoms)
    5
    """
    obmol = ob.OBMol()
    obconversion = ob.OBConversion()

    formatok = obconversion.SetInFormat(format)
    if not formatok:
        raise ValueError,"%s is not a recognised OpenBabel format" % format

    obconversion.ReadString(obmol, string)
    return Molecule(obmol)

class Outputfile(object):
    """Represent a file to which *output* is to be sent.
    
    Although it's possible to write a single molecule to a file by
    calling the write() method of a molecule, if multiple molecules
    are to be written to the same file you should use the Outputfile
    class.
    
    Required parameters:
       format
       filename
    Optional parameters:
       overwrite (default is False) -- if the output file already exists,
                                       should it be overwritten?
    Methods:
       write(molecule)
    """
    def __init__(self, format, filename, overwrite=False):
        self.format = format
        self.filename = filename
        if not overwrite and os.path.isfile(self.filename):
            raise IOError, "%s already exists. Use 'overwrite=False' to overwrite it." % self.filename
        self.obConversion = ob.OBConversion()
        formatok = self.obConversion.SetOutFormat(self.format)
        if not formatok:
            raise ValueError,"%s is not a recognised OpenBabel format" % format
        self.total = 0 # The total number of molecules written to the file
    
    def write(self, molecule):
        """Write a molecule to the output file.
        
        Required parameters:
           molecule
        """
        if self.total==0:
            self.obConversion.WriteFile(molecule.OBMol, self.filename)
        else:
            self.obConversion.Write(molecule.OBMol)
        self.total += 1


class Molecule(object):
    """Represent a Pybel molecule.

    Optional parameters:
       OBMol -- an Open Babel molecule (default is None)
    
    An empty Molecule is created if an Open Babel molecule is not provided.
    
    Attributes:
       atoms, charge, dim, energy, exactmass, flags, formula, 
       mod, molwt, spin, sssr, title.
    (refer to the Open Babel library documentation for more info).
    
    Methods:
       write()
      
    The original Open Babel molecule can be accessed using the attribute:
       OBMol
    """
    _getmethods = {
        'conformers':'GetConformers',
        # 'coords':'GetCoordinates', you can access the coordinates the atoms elsewhere
        # 'data':'GetData', has been removed
        'dim':'GetDimension',
        'energy':'GetEnergy',
        'exactmass':'GetExactMass',
        'flags':'GetFlags',
        'formula':'GetFormula',
        # 'internalcoord':'GetInternalCoord', # Causes SWIG warning
        'mod':'GetMod',
        'molwt':'GetMolWt',
        'sssr':'GetSSSR',
        'title':'GetTitle',
        'charge':'GetTotalCharge',
        'spin':'GetTotalSpinMultiplicity'
    }
    
    def __init__(self, OBMol=None):

        self.OBMol = OBMol
        if not self.OBMol:
            self.OBMol = ob.OBMol()

    def __getattr__(self, attr):
        """Return the value of an attribute

        Note: The values are calculated on-the-fly. You may want to store the value in
        a variable if you repeatedly access the same attribute.
        """
        # This function is not accessed in the case of OBMol
        if attr == "atoms":
            # Create an atoms attribute on-the-fly
            return [ Atom(self.OBMol.GetAtom(i+1),i+1) for i in range(self.OBMol.NumAtoms()) ]
        elif attr in self._getmethods:
            # Call the OB Method to find the attribute value
            return getattr(self.OBMol, self._getmethods[attr])()
        else:
            raise AttributeError, "Molecule has no attribute %s" % attr

    def __iter__(self):
        """Iterate over the Atoms of the Molecule.
        
        This allows constructions such as the following:
           for atom in mymol:
               print atom
        """
        for atom in self.atoms:
            yield atom

    def write(self, format="SMI", filename=None, overwrite=False):
        """Write the molecule to a file or return a string.
        
        Optional parameters:
           format -- default is "SMI"
           filename -- default is None
           overwite -- default is False

        If a filename is specified, the result is written to a file.
        Otherwise, a string is returned containing the result.
        The overwrite flag is ignored if a filename is not specified.
        It controls whether to overwrite an existing file.
        """

        obconversion = ob.OBConversion()
        formatok = obconversion.SetOutFormat(format)
        if not formatok:
            raise ValueError,"%s is not a recognised OpenBabel format" % format

        if filename:
            if not overwrite and os.path.isfile(filename):
                raise IOError, "%s already exists. Use 'overwrite=False' to overwrite it." % filename
            obconversion.WriteFile(self.OBMol,filename)
        else:
            return obconversion.WriteString(self.OBMol)

    def __str__(self):
        return self.write()


class Atom(object):
    """Represent a Pybel atom.

    Optional parameters:
       OBAtom -- an Open Babel Atom (default is None)
       index -- the index of the atom in the molecule (default is None)
     
    An empty Atom is created if an Open Babel atom is not provided.
    
    Attributes:
       atomicmass, atomicnum, cidx, coords, coordidx, exactmass,
       formalcharge, heavyvalence, heterovalence, hyb, idx,
       implicitvalence, index, isotope, partialcharge, spin, type,
       valence, vector.

    (refer to the Open Babel library documentation for more info).
    
    The original Open Babel atom can be accessed using the attribute:
       OBAtom
    """
    
    _getmethods = {
        'atomicmass':'GetAtomicMass',
        'atomicnum':'GetAtomicNum',
        'cidx':'GetCIdx',
        'coordidx':'GetCoordinateIdx',
        # 'data':'GetData', has been removed
        'exactmass':'GetExactMass',
        'formalcharge':'GetFormalCharge',
        'heavyvalence':'GetHvyValence',
        'heterovalence':'GetHeteroValence',
        'hyb':'GetHyb',
        'idx':'GetIdx',
        'implicitvalence':'GetImplicitValence',
        'isotope':'GetIsotope',
        'partialcharge':'GetPartialCharge',
        'spin':'GetSpinMultiplicity',
        'type':'GetType',
        'valence':'GetValence',
        'vector':'GetVector',
        }

    def __init__(self, OBAtom=None, index=None):
        if not OBAtom:
            OBAtom = ob.OBAtom()
        self.OBAtom = OBAtom
        # For the moment, I will remember the index of the atom in the molecule...
        # I'm not sure if this is useful, though.
        self.index = index
        
    def __getattr__(self, attr):
        if attr == "coords":
            return (self.OBAtom.GetX(), self.OBAtom.GetY(), self.OBAtom.GetZ())
        elif attr in self._getmethods:
            return getattr(self.OBAtom, self._getmethods[attr])()
        else:
            raise AttributeError, "Molecule has no attribute %s" % attr

    def __str__(self):
        """Create a string representation of the atom.

        >>> a = Atom()
        >>> print a
        Atom: 0 (0.0, 0.0, 0.0)
        """
        return "Atom: %d %s" % (self.atomicnum, self.coords.__str__())

class Smarts(object):
    """A Smarts Pattern Matcher

    Required parameters:
       smartspattern
    
    Methods:
       findall()
    
    Example:
    >>> mol = readstring("smi","CCN(CC)CC") # triethylamine
    >>> smarts = Smarts("[#6][#6]") # Matches an ethyl group
    >>> print smarts.findall(mol) 
    [(1, 2), (4, 5), (6, 7)]
    """
    def __init__(self,smartspattern):
        """Initialise with a SMARTS pattern."""
        self.obsmarts = ob.OBSmartsPattern()
        self.obsmarts.Init(smartspattern)
    def findall(self,molecule):
        """Find all matches of the SMARTS pattern to a particular molecule.
        
        Required parameters:
           molecule
        """
        self.obsmarts.Match(molecule.OBMol)
        return [x for x in self.obsmarts.GetUMapList()]
        
if __name__=="__main__":
    import doctest
    doctest.testmod()
    
