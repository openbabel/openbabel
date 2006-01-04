import openbabel as ob

def readstring(format,string):
    """Read in a molecule from a string.

    >>> input = "C1=CC=CS1"
    >>> mymol = readstring("smi",input)
    >>> len(mymol.atoms)
    5
    """
    obmol = ob.OBMol()
    obconversion = ob.OBConversion()

# TO DO: Validate the format before passing to SetInFormat
#        and print a helpful list of alternatives if not valid
    formatok = obconversion.SetInFormat(format)
    if not formatok:
        raise ValueError,"%s is not a recognised OpenBabel format" % format

    obconversion.ReadString(obmol,string)
    return Molecule(obmol)


class Molecule(object):
    """Represent a molecule."""

    _getmethods = {
        'conformers':'GetConformers',
        'coords':'GetCoordinates',
        'data':'GetData',
        'dim':'GetDimension',
        'energy':'GetEnergy',
        'exactmass':'GetExactMass',
        'flags':'GetFlags',
        'formula':'GetFormula',
        'internalcoord':'GetInternalCoord',
        'mod':'GetMod',
        'molwt':'GetMolWt',
        'sssr':'GetSSSR',
        'title':'GetTitle',
        'charge':'GetTotalCharge',
        'spin':'GetTotalSpinMultiplicity'
    }
    
    def __init__(self,obmol=None):

        self.OBMol = obmol
        if not self.OBMol:
            self.OBMol = ob.OBMol()

        for x,v in self._getmethods.iteritems():
            setattr( self,x,getattr(self,x) )
        self.atoms = self.atoms
    
    def __getattr__(self,attr):
        if attr == "atoms":
            listofatoms = [ Atom(self.OBMol.GetAtom(i+1),i+1) for i in range(self.OBMol.NumAtoms()) ]
            return listofatoms
        elif attr in self._getmethods:
            return getattr(self.OBMol,self._getmethods[attr])()
        else:
            raise AttributeError,"Cannot find %s" % attr

    def __iter__(self):
        """Iterate over the Atoms of the Molecule.
        
        This allows constructions such as the following:
           for atom in mymol:
               print atom
        """
        for atom in self.atoms:
            yield atom

    def write(self,format="SMI",filename=None):
        """Write the Molecule to a file or return a string."""

        obconversion = ob.OBConversion()
        formatok = obconversion.SetOutFormat(format)
        if not formatok:
            raise ValueError,"%s is not a recognised OpenBabel format" % format

        if filename:
            obconversion.WriteFile(self.OBMol,filename)
        else:
            return obconversion.WriteString(self.OBMol)

    def __str__(self):
        return self.write()


class Atom(object):
    """Represent an atom."""
    
    _getmethods = {
        'atomicmass':'GetAtomicMass',
        'atomicnum':'GetAtomicNum',
        'cidx':'GetCIdx',
        'coordidx':'GetCoordinateIdx',
        'data':'GetData',
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

    def __init__(self,OBAtom=None,index=None):
        # For the moment, I will remember the index of the atom in the molecule...
        # I'm not sure if this is useful, though.
        if not OBAtom:
            OBAtom = ob.OBAtom()
        self.OBAtom = OBAtom
        self.index = index
        
        for x,v in self._getmethods.iteritems():
            setattr( self,x,getattr(self,x) )
        self.coords = self.coords
    
    def __getattr__(self,attr):
        # Reminder to add corresponding __setattr__ methods
        if attr == "coords":
            return (self.OBAtom.GetX(),self.OBAtom.GetY(),self.OBAtom.GetZ())
        elif attr in self._getmethods:
            return getattr(self.OBAtom,self._getmethods[attr])()
        else:
            raise AttributeError,"Cannot find %s" % attr

    def __str__(self):
        """Create a string representation of the atom.

        >>> a = Atom()
        >>> print a
        Atom: 0 (0.0, 0.0, 0.0)
        """
        return "Atom: %d %s" % (self.atomicnum, self.coords.__str__())
            
        
if __name__=="__main__":
    import doctest,pyopenbabel
    doctest.testmod(pyopenbabel)
    
