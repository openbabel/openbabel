import math
import os.path
import tempfile
import openbabel as ob

try:
    import oasa
    import oasa.cairo_out
except ImportError: #pragma: no cover
    oasa = None

try:
    import Tkinter as tk
    import Image as PIL
    import ImageTk as piltk
except ImportError: #pragma: no cover
    tk = None

def _formatstodict(list):
    broken = [x.replace("[Read-only]", "").replace("[Write-only]","").split(" -- ") for x in list]
    broken = [(x,y.strip()) for x,y in broken]
    return dict(broken)
_obconv = ob.OBConversion()
informats = _formatstodict(_obconv.GetSupportedInputFormat())
outformats = _formatstodict(_obconv.GetSupportedOutputFormat())

def _getplugins(findplugin, names):
    plugins = dict([(x, findplugin(x)) for x in names if findplugin(x)])
    return plugins

descs = ['LogP', 'MR', 'TPSA']
_descdict = _getplugins(ob.OBDescriptor.FindType, descs)
fps = ['FP2', 'FP3', 'FP4']
_fingerprinters = _getplugins(ob.OBFingerprint.FindFingerprint, fps)
forcefields = ['UFF', 'MMFF94', 'Ghemical']
_forcefields = _getplugins(ob.OBForceField.FindType, forcefields)
operations = ['Gen3D']
_operations = _getplugins(ob.OBOp.FindType, operations)

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
    if not os.path.isfile(filename):
        raise IOError, "No such file: '%s'" % filename
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

    success = obconversion.ReadString(obmol, string)
    if not success:
        raise IOError, "Failed to convert '%s' to format '%s'" % (
            string, format)
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
       write(molecule), close()
    """
    def __init__(self, format, filename, overwrite=False):
        self.format = format
        self.filename = filename
        if not overwrite and os.path.isfile(self.filename):
            raise IOError, "%s already exists. Use 'overwrite=True' to overwrite it." % self.filename
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
        if not self.filename:
            raise IOError, "Outputfile instance is closed."

        if self.total==0:
            self.obConversion.WriteFile(molecule.OBMol, self.filename)
        else:
            self.obConversion.Write(molecule.OBMol)
        self.total += 1

    def close(self):
        """Close the Outputfile to further writing."""
        self.obConversion.CloseOutFile()
        self.filename = None

class Molecule(object):
    """Represent a Pybel molecule.

    Required parameter:
       OBMol -- an Open Babel OBMol
       or
       Molecule -- any type of cinfony Molecule (e.g. one from cinfony.rdkit)

    If a cinfony Molecule is provided it will be converted into a pybel Molecule.       
    
    Attributes:
       atoms, charge, data, dim, energy, exactmass, flags, formula, 
       mod, molwt, spin, sssr, title, unitcell.
    (refer to the Open Babel library documentation for more info).
    
    Methods:
       addh(), calcfp(), calcdesc(), draw(), localopt(), make3D(), removeh(), write() 
      
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
    _cinfony = True

    def __init__(self, OBMol):
        
        if hasattr(OBMol, "_cinfony"):
            a, b = OBMol._exchange
            if a == 0:
                mol = readstring("smi", b)
            else:
                mol = readstring("mol", b)
            OBMol = mol.OBMol

        self.OBMol = OBMol
 
    def __getattr__(self, attr):
        """Return the value of an attribute

        Note: The values are calculated on-the-fly. You may want to store the value in
        a variable if you repeatedly access the same attribute.
        """
        # This function is not accessed in the case of OBMol
        if attr == "atoms":
            # Create an atoms attribute on-the-fly
            return [ Atom(self.OBMol.GetAtom(i+1)) for i in range(self.OBMol.NumAtoms()) ]
        elif attr == "data":
            # Create a data attribute on-the-fly
            return MoleculeData(self.OBMol)
        elif attr == "unitcell":
            # Create a unitcell attribute on-the-fly
            unitcell = self.OBMol.GetData(ob.UnitCell)
            if unitcell:
                return ob.toUnitCell(unitcell)
            else:
                raise AttributeError, "Molecule has no attribute 'unitcell'"
        elif attr in self._getmethods:
            # Call the OB Method to find the attribute value
            return getattr(self.OBMol, self._getmethods[attr])()
        elif attr == "_exchange":
            if self.OBMol.HasNonZeroCoords():
                return (1, self.write("mol"))
            else:
                return (0, self.write("can").split()[0])
        else:
            raise AttributeError, "Molecule has no attribute '%s'" % attr

    def __iter__(self):
        """Iterate over the Atoms of the Molecule.
        
        This allows constructions such as the following:
           for atom in mymol:
               print atom
        """
        for atom in self.atoms:
            yield atom

    def calcdesc(self, descnames=[]):
        """Calculate descriptor values.

        Optional parameter:
           descnames -- a list of names of descriptors

        If descnames is not specified, the full list of Open Babel
        descriptors is calculated: LogP, PSA and MR.
        """
        if not descnames:
            descnames = descs
        ans = {}
        for descname in descnames:
            try:
                desc = _descdict[descname]
            except KeyError:
                raise ValueError, "%s is not a recognised Open Babel descriptor type" % descname
            ans[descname] = desc.Predict(self.OBMol)
        return ans
    
    def calcfp(self, fptype="FP2"):
        """Calculate a molecular fingerprint.
        
        Optional parameters:
           fptype -- the name of the Open Babel fingerprint type.

        If fptype is not specified, the FP2 fingerprint type
        is used. See the Open Babel library documentation for more
        details.
        """
        fp = ob.vectorUnsignedInt()
        try:
            fingerprinter = _fingerprinters[fptype]
        except KeyError:
            raise ValueError, "%s is not a recognised Open Babel Fingerprint type" % fptype
        fingerprinter.GetFingerprint(self.OBMol, fp)
        return Fingerprint(fp)

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
                raise IOError, "%s already exists. Use 'overwrite=True' to overwrite it." % filename
            obconversion.WriteFile(self.OBMol,filename)
        else:
            return obconversion.WriteString(self.OBMol)

    def localopt(self, forcefield="MMFF94", steps=500):
        """Locally optimize the coordinates.
        
        Optional parameters:
           forcefield -- default is "MMFF94"
           steps -- default is 500

        If the molecule does not have 2D or 3D coordinates, make3D() is
        called before the optimization. Note that the molecule needs
        to have explicit hydrogens. If not, call addh().
        """
        
        if not (self.OBMol.Has2D() or self.OBMol.Has3D()):
            self.make3D()
        ff = _forcefields[forcefield]
        ff.Setup(self.OBMol)
        ff.SteepestDescent(steps)
        ff.GetCoordinates(self.OBMol)
    
##    def globalopt(self, forcefield="MMFF94", steps=1000):
##        if not (self.OBMol.Has2D() or self.OBMol.Has3D()):
##            self.make3D()
##        self.localopt(forcefield, 250)
##        ff = _forcefields[forcefield]
##        numrots = self.OBMol.NumRotors()
##        if numrots > 0:
##            ff.WeightedRotorSearch(numrots, int(math.log(numrots + 1) * steps))
##        ff.GetCoordinates(self.OBMol)
    
    def make3D(self, forcefield="MMFF94", steps=50):
        """Generate 3D coordinates.
        
        Optional parameters:
           forcefield -- default is "MMFF94"
           steps -- default is 50

        Once coordinates are generated, hydrogens are added and a quick
        local optimization is carried out with 50 steps and the
        MMFF94 forcefield. Call localopt() if you want
        to improve the coordinates further.
        """
        _operations['Gen3D'].Do(self.OBMol)
        self.addh()
        self.localopt(forcefield, steps)

    def addh(self):
        """Add hydrogens."""
        self.OBMol.AddHydrogens()

    def removeh(self):
        """Remove hydrogens."""
        self.OBMol.DeleteHydrogens()
        
    def __str__(self):
        return self.write()

    def draw(self, show=True, filename=None, update=False, usecoords=False):
        """Create a 2D depiction of the molecule.

        Optional parameters:
          show -- display on screen (default is True)
          filename -- write to file (default is None)
          update -- update the coordinates of the atoms to those
                    determined by the structure diagram generator
                    (default is False)
          usecoords -- don't calculate 2D coordinates, just use
                       the current coordinates (default is False)

        OASA is used for 2D coordinate generation and depiction. Tkinter and
        Python Imaging Library are required for image display.
        """
        etab = ob.OBElementTable()

        if not oasa:
            errormessage = ("OASA not found, but is required for 2D structure "
                            "generation and depiction. OASA is part of BKChem. "
                            "See installation instructions for more "
                            "information.")
            raise ImportError, errormessage
        mol = oasa.molecule()
        for atom in self.atoms:
            v = mol.create_vertex()
            v.symbol = etab.GetSymbol(atom.atomicnum)
            v.charge = atom.formalcharge
            if usecoords:
                v.x, v.y, v.z = atom.coords[0] * 30., atom.coords[1] * 30., 0.0
            mol.add_vertex(v)

        for bond in ob.OBMolBondIter(self.OBMol):
            e = mol.create_edge()
            e.order = bond.GetBO()
            if bond.IsHash():
                e.type = "h"
            elif bond.IsWedge():
                e.type = "w"
            mol.add_edge(bond.GetBeginAtomIdx() - 1,
                         bond.GetEndAtomIdx() - 1,
                         e)
        # I'm sure there's a more elegant way to do the following, but here goes...
        # let's set the stereochemistry around double bonds
        self.write("smi") # Perceive UP/DOWNness
        for bond in ob.OBMolBondIter(self.OBMol):
            ends = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if bond.GetBO() == 2:
                stereobonds = [[b for b in ob.OBAtomBondIter(self.OBMol.GetAtom(x)) if b.GetIdx() != bond.GetIdx() and (b.IsUp() or b.IsDown())]
                               for x in ends]
                if stereobonds[0] and stereobonds[1]: # Needs to be defined at either end
                    if stereobonds[0][0].IsUp() == stereobonds[1][0].IsUp():
                        # Either both up or both down
                        stereo = oasa.stereochemistry.cis_trans_stereochemistry.SAME_SIDE
                    else:
                        stereo = oasa.stereochemistry.cis_trans_stereochemistry.OPPOSITE_SIDE
                    atomids = [(b[0].GetBeginAtomIdx(), b[0].GetEndAtomIdx()) for b in stereobonds]
                    extremes = []
                    for id, end in zip(ends, atomids):
                        if end[0] == id:
                            extremes.append(end[1])
                        else:
                            extremes.append(end[0])
                    center = mol.get_edge_between(mol.atoms[ends[0] - 1], mol.atoms[ends[1] - 1])
                    st = oasa.stereochemistry.cis_trans_stereochemistry(
                              center = center, value = stereo,
                              references = (mol.atoms[extremes[0] - 1], mol.atoms[ends[0] - 1],
                                            mol.atoms[ends[1] - 1], mol.atoms[extremes[1] - 1]))
                    mol.add_stereochemistry(st)
        
        mol.remove_unimportant_hydrogens()
        if not usecoords:
            oasa.coords_generator.calculate_coords(mol, bond_length=30)
            if update:
                newcoords = [(v.x / 30., v.y / 30., 0.0) for v in mol.vertices]
                for atom, newcoord in zip(ob.OBMolAtomIter(self.OBMol), newcoords):
                    atom.SetVector(*newcoord)
        if filename or show:
            maxx = max([v.x for v in mol.vertices])
            minx = min([v.x for v in mol.vertices])
            maxy = max([v.y for v in mol.vertices])
            miny = min([v.y for v in mol.vertices])
            maxcoord = max(maxx - minx, maxy - miny)
            fontsize = 16
            bondwidth = 6
            linewidth = 2
            if maxcoord > 270: # 300  - margin * 2
                for v in mol.vertices:
                    v.x *= 270. / maxcoord
                    v.y *= 270. / maxcoord
                fontsize *= math.sqrt(270. / maxcoord)
                bondwidth *= math.sqrt(270. / maxcoord)
                linewidth *= math.sqrt(270. / maxcoord)
            if filename:
                filedes = None
            else:
                filedes, filename = tempfile.mkstemp()
            
            canvas = oasa.cairo_out.cairo_out()
            canvas.show_hydrogens_on_hetero = True
            canvas.font_size = fontsize
            canvas.bond_width = bondwidth
            canvas.line_width = linewidth
            canvas.mol_to_cairo(mol, filename)
            if show:
                if not tk:
                    errormessage = ("Tkinter or Python Imaging "
                                    "Library not found, but is required for image "
                                    "display. See installation instructions for "
                                    "more information.")
                    raise ImportError, errormessage
                root = tk.Tk()
                root.title((hasattr(self, "title") and self.title)
                           or self.__str__().rstrip())
                frame = tk.Frame(root, colormap="new", visual='truecolor').pack()
                image = PIL.open(filename)
                imagedata = piltk.PhotoImage(image)
                label = tk.Label(frame, image=imagedata).pack()
                quitbutton = tk.Button(root, text="Close", command=root.destroy).pack(fill=tk.X)
                root.mainloop()
            if filedes:
                os.close(filedes)
                os.remove(filename)

class Atom(object):
    """Represent a Pybel atom.

    Required parameter:
       OBAtom -- an Open Babel OBAtom
        
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

    def __init__(self, OBAtom):
        self.OBAtom = OBAtom
        
    def __getattr__(self, attr):
        if attr == "coords":
            return (self.OBAtom.GetX(), self.OBAtom.GetY(), self.OBAtom.GetZ())
        elif attr in self._getmethods:
            return getattr(self.OBAtom, self._getmethods[attr])()
        else:
            raise AttributeError, "Atom has no attribute %s" % attr

    def __str__(self):
        """Create a string representation of the atom.

        >>> a = Atom()
        >>> print a
        Atom: 0 (0.0, 0.0, 0.0)
        """
        c = self.coords
        return "Atom: %d (%.2f %.2f %.2f)" % (self.atomicnum, c[0], c[1], c[2])

def _findbits(fp, bitsperint):
    """Find which bits are set in a list/vector.

    This function is used by the Fingerprint class.

    >>> _findbits([13, 71], 8)
    [1, 3, 4, 9, 10, 11, 15]
    """
    ans = []
    start = 1
    for x in fp:
        i = start
        while x > 0:
            if x % 2:
                ans.append(i)
            x >>= 1
            i += 1
        start += bitsperint
    return ans
        
class Fingerprint(object):
    """A Molecular Fingerprint.
    
    Required parameters:
       obFingerprint -- a vector calculated by OBFingerprint.FindFingerprint()

    Attributes:
       fp -- the original obFingerprint
       bits -- a list of bits set in the Fingerprint

    Methods:
       The "|" operator can be used to calculate the Tanimoto coeff. For example,
       given two Fingerprints 'a', and 'b', the Tanimoto coefficient is given by:
          tanimoto = a | b
    """
    def __init__(self, obFingerprint):
        self.fp = obFingerprint
    def __or__(self, other):
        return ob.OBFingerprint.Tanimoto(self.fp, other.fp)
    def __getattr__(self, attr):
        if attr == "bits":
            # Create a bits attribute on-the-fly
            return _findbits(self.fp, ob.OBFingerprint.Getbitsperint())
        else:
            raise AttributeError, "Fingerprint has no attribute %s" % attr
    def __str__(self):
        return ", ".join([str(x) for x in self.fp])

class Smarts(object):
    """A Smarts Pattern Matcher

    Required parameters:
       smartspattern
    
    Methods:
       findall(molecule)
    
    Example:
    >>> mol = readstring("smi","CCN(CC)CC") # triethylamine
    >>> smarts = Smarts("[#6][#6]") # Matches an ethyl group
    >>> print smarts.findall(mol) 
    [(1, 2), (4, 5), (6, 7)]

    The numbers returned are the indices (starting from 1) of the atoms
    that match the SMARTS pattern. In this case, there are three matches
    for each of the three ethyl groups in the molecule.
    """
    def __init__(self,smartspattern):
        """Initialise with a SMARTS pattern."""
        self.obsmarts = ob.OBSmartsPattern()
        success = self.obsmarts.Init(smartspattern)
        if not success:
            raise IOError, "Invalid SMARTS pattern"
    def findall(self,molecule):
        """Find all matches of the SMARTS pattern to a particular molecule.
        
        Required parameters:
           molecule
        """
        self.obsmarts.Match(molecule.OBMol)
        return [x for x in self.obsmarts.GetUMapList()]
        
class MoleculeData(object):
    """Store molecule data in a dictionary-type object
    
    Required parameters:
      obmol -- an Open Babel OBMol 

    Methods and accessor methods are like those of a dictionary except
    that the data is retrieved on-the-fly from the underlying OBMol.

    Example:
    >>> mol = readfile("sdf", 'head.sdf').next()
    >>> data = mol.data
    >>> print data
    {'Comment': 'CORINA 2.61 0041  25.10.2001', 'NSC': '1'}
    >>> print len(data), data.keys(), data.has_key("NSC")
    2 ['Comment', 'NSC'] True
    >>> print data['Comment']
    CORINA 2.61 0041  25.10.2001
    >>> data['Comment'] = 'This is a new comment'
    >>> for k,v in data.iteritems():
    ...    print k, "-->", v
    Comment --> This is a new comment
    NSC --> 1
    >>> del data['NSC']
    >>> print len(data), data.keys(), data.has_key("NSC")
    1 ['Comment'] False
    """
    def __init__(self, obmol):
        self._mol = obmol
    def _data(self):
        return [ob.toPairData(x) for x in self._mol.GetData() if x.GetDataType()==ob.PairData or x.GetDataType()==ob.CommentData]
    def _testforkey(self, key):
        if not key in self:
            raise KeyError, "'%s'" % key
    def keys(self):
        return [x.GetAttribute() for x in self._data()]
    def values(self):
        return [x.GetValue() for x in self._data()]
    def items(self):
        return zip(self.keys(), self.values())
    def __iter__(self):
        return iter(self.keys())
    def iteritems(self):
        return iter(self.items())
    def __len__(self):
        return len(self._data())
    def __contains__(self, key):
        return self._mol.HasData(key)
    def __delitem__(self, key):
        self._testforkey(key)
        self._mol.DeleteData(self._mol.GetData(key))
    def clear(self):
        for key in self:
            del self[key]
    def has_key(self, key):
        return key in self
    def update(self, dictionary):
        for k, v in dictionary.iteritems():
            self[k] = v
    def __getitem__(self, key):
        self._testforkey(key)
        return ob.toPairData(self._mol.GetData(key)).GetValue()
    def __setitem__(self, key, value):
        if key in self:
            pairdata = ob.toPairData(self._mol.GetData(key))
            pairdata.SetValue(str(value))
        else:
            pairdata = ob.OBPairData()
            pairdata.SetAttribute(key)
            pairdata.SetValue(str(value))
            pairdata.thisown = 0 # So that SWIG Proxy will not delete pairdata
            self._mol.SetData(pairdata)
    def __repr__(self):
        return dict(self.iteritems()).__repr__()
 
if __name__=="__main__": #pragma: no cover
    import doctest
    doctest.testmod(verbose=True)
