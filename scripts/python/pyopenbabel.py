import openbabel as ob

class Molecule(object):
    """Represent a molecule."""

    _properties = {
        'insertprops':'here'
    }
    
    def __init__(self):
        self.OBMol = ob.OBMol()
    def __getattr__(self,attr):
        if attr == "atoms":
            listofatoms = [ Atom(self.mol.GetAtom(i+1),i+1) for i in range(0,self.mol.NumAtoms()) ]
            return listofatoms
    def newatom(self):
        self.mol.NewAtom()

class Atom(object):
    """Represent an atom."""
    
    _properties = {
        'isotope':'GetIsotope',
        'hyb':'GetHyb',
        'cidx':'GetCIdx',
        'valence':'GetValence',
        'idx':'GetIdx',
        'atomicnum':'GetAtomicNum',
        'spin':'GetSpinMultiplicity',
        'formalcharge':'GetFormalCharge',
        'coordidx':'GetCoordinateIdx',
        'atomicmass':'GetAtomicMass',
        'heterovalence':'GetHeteroValence',
        'data':'GetData',
        'heavyvalence':'GetHvyValence',
        'vector':'GetVector',
        'type':'GetType',
        'exactmass':'GetExactMass',
        'implicitvalence':'GetImplicitValence',
        'partialcharge':'GetPartialCharge',
        }

    def __init__(self,OBAtom=None,index=None):
        # For the moment, I will remember the index of the atom in the molecule
        # I'm not sure if this is useful.
        if not OBAtom:
            OBAtom = ob.OBAtom()
        self.OBAtom = OBAtom
        self.index = index
        
        for x,v in self._properties.iteritems():
            setattr( self,x,getattr(self,x) )
        self.coords = self.coords
    
    def __getattr__(self,attr):
        # Reminder to add corresponding __setattr__ methods
        if attr == "coords":
            return (self.OBAtom.GetX(),self.OBAtom.GetY(),self.OBAtom.GetZ())
        elif attr in self._properties:
            return getattr(self.OBAtom,self._properties[attr])()

    def __str__(self):
        line = []
        line.append("Atom: ")
        line.append("(%f,%f,%f)" % (self.x,self.y,self.z))        
        return "".join(line)
            
        
if __name__=="__main__":
    mol = Molecule()
    a = Atom()
    b = Atom()
     
