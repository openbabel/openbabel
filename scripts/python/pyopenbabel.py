import openbabel as ob

class Molecule(object):
    """Represent a molecule."""

    _getmethods = {
        'conformers':'GetConformers',
        'coordinates':'GetCoordinates',
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
    
    def __init__(self):
        self.OBMol = ob.OBMol()
        for x,v in self._getmethods.iteritems():
            setattr( self,x,getattr(self,x) )
        self.atoms = self.atoms
    
    def __getattr__(self,attr):
        if attr == "atoms":
            listofatoms = [ Atom(self.OBMol.GetAtom(i+1),i+1) for i in range(0,self.OBMol.NumAtoms()) ]
            return listofatoms
        elif attr in self._getmethods:
            return getattr(self.OBMol,self._getmethods[attr])()
        else:
            raise AttributeError,"Cannot find %s" % attr

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
        line = []
        line.append("Atom: ")
        line.append("(%f,%f,%f)" % (self.x,self.y,self.z))        
        return "".join(line)
            
        
if __name__=="__main__":
    mol = Molecule()
    a = Atom()
    b = Atom()
     
