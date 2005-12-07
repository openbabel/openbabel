import openbabel as ob

class Molecule(object):
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
    
    _notwanted = {
        'coord':'GetCoordinate',
        'nextatom':'GetNextAtom',
        'x':'GetX',
        'y':'GetY',
        'z':'GetZ',
        }
    _notimplemeneted = {
        'angle':'GetAngle',
        'distance':'GetDistance',
        }
    _needstwoargs = {
        'bond':'GetBond',
        }
    _needsthreeargs = {
        'newbondvector':'GetNewBondVector'
        }
    _segfaults = {
        'residue':'GetResidue',
        }
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
        if not OBAtom:
            OBAtom = ob.OBAtom()
        self.OBAtom = OBAtom
        self.index = index
        
        for x,v in self._properties.iteritems():
            # setattr(self,x,getattr(self.OBAtom,v)())
            setattr( self,x,getattr(self,x) )
        self.coords = self.coords
    
    def __getattr__(self,attr):
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
     
