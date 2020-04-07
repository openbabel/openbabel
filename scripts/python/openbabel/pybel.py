# -*. coding: utf-8 -*-
# Copyright (c) 2008-2012, Noel O'Boyle; 2012, Adrià Cereto-Massagué
# All rights reserved.
#
# This file is part of Cinfony.
# The contents are covered by the terms of the GPL v2 license
# which is included in the file LICENSE_GPLv2.txt.

"""
pybel - A Cinfony module for accessing Open Babel

Global variables:
  ob - the underlying SWIG bindings for Open Babel
  informats - a dictionary of supported input formats
  outformats - a dictionary of supported output formats
  descs - a list of supported descriptors
  fps - a list of supported fingerprint types
  forcefields - a list of supported forcefields
"""

import sys
import os.path
import tempfile
import xml.etree.ElementTree as ET

if sys.platform[:4] == "java":
    import org.openbabel as ob
    import java.lang.System
    java.lang.System.loadLibrary("openbabel_java")
    _obfuncs = ob.openbabel_java
    _obconsts = ob.openbabel_javaConstants
    import javax
elif sys.platform[:3] == "cli":
    import System
    import clr
    clr.AddReference('System.Windows.Forms')
    clr.AddReference('System.Drawing')

    from System.Windows.Forms import Application, DockStyle, Form, PictureBox
    from System.Windows.Forms import PictureBoxSizeMode
    from System.Drawing import Image, Size

    _obdotnet = os.environ["OBDOTNET"]
    if _obdotnet[0] == '"':  # Remove trailing quotes
        _obdotnet = _obdotnet[1:-1]
    clr.AddReferenceToFileAndPath(os.path.join(_obdotnet, "OBDotNet.dll"))
    import OpenBabel as ob
    _obfuncs = ob.openbabel_csharp
    _obconsts = ob.openbabel_csharp
else:
    from . import openbabel as ob
    _obfuncs = _obconsts = ob
    try:
        if sys.version_info[0] >= 3:
            import tkinter as tk
        else:
            import Tkinter as tk
        from PIL import Image as PIL
        from PIL import ImageTk as piltk
    except ImportError:  # pragma: no cover
        tk = None


def _formatstodict(list):
    if sys.platform[:4] == "java":
        list = [list.get(i) for i in range(list.size())]
    broken = [x.replace("[Read-only]", "").replace("[Write-only]", "").split(
              " -- ") for x in list]
    broken = [(x, y.strip()) for x, y in broken]
    return dict(broken)


def _getplugins(findplugin, names):
    return dict([(x, findplugin(x)) for x in names if findplugin(x)])


def _getpluginnames(ptype):
    if sys.platform[:4] == "cli":
        plugins = ob.VectorString()
    else:
        plugins = ob.vectorString()
    ob.OBPlugin.ListAsVector(ptype, None, plugins)
    if sys.platform[:4] == "java":
        plugins = [plugins.get(i) for i in range(plugins.size())]
    return [x.split()[0] for x in plugins if x.strip()]

_obconv = ob.OBConversion()
_builder = ob.OBBuilder()

informats = _formatstodict(_obconv.GetSupportedInputFormat())
"""A dictionary of supported input formats"""
outformats = _formatstodict(_obconv.GetSupportedOutputFormat())
"""A dictionary of supported output formats"""

descs = _getpluginnames("descriptors")
"""A list of supported descriptors"""
_descdict = _getplugins(ob.OBDescriptor.FindType, descs)

fps = [_x.lower() for _x in _getpluginnames("fingerprints")]
"""A list of supported fingerprint types"""
_fingerprinters = _getplugins(ob.OBFingerprint.FindFingerprint, fps)

forcefields = [_x.lower() for _x in _getpluginnames("forcefields")]
"""A list of supported forcefields"""
_forcefields = _getplugins(ob.OBForceField.FindType, forcefields)

charges = [_x.lower() for _x in _getpluginnames("charges")]
"""A list of supported charge models"""
_charges = _getplugins(ob.OBChargeModel.FindType, charges)

operations = _getpluginnames("ops")
"""A list of supported operations"""
_operations = _getplugins(ob.OBOp.FindType, operations)

ipython_3d = False
"""Toggles 2D vs 3D molecule representations in IPython notebook"""


def readfile(format, filename, opt=None):
    """Iterate over the molecules in a file.

    Required parameters:
       format - see the informats variable for a list of available
                input formats
       filename

    Optional parameters:
       opt    - a dictionary of format-specific options
                For format options with no parameters, specify the
                value as None.

    You can access the first molecule in a file using the next() method
    of the iterator (or the next() keyword in Python 3):
        mol = readfile("smi", "myfile.smi").next() # Python 2
        mol = next(readfile("smi", "myfile.smi"))  # Python 3

    You can make a list of the molecules in a file using:
        mols = list(readfile("smi", "myfile.smi"))

    You can iterate over the molecules in a file as shown in the
    following code snippet:
    >>> atomtotal = 0
    >>> for mol in readfile("sdf", "head.sdf"):
    ...     atomtotal += len(mol.atoms)
    ...
    >>> print atomtotal
    43
    """
    if opt is None:
        opt = {}
    obconversion = ob.OBConversion()
    formatok = obconversion.SetInFormat(format)
    for k, v in opt.items():
        if v is None:
            obconversion.AddOption(k, obconversion.INOPTIONS)
        else:
            obconversion.AddOption(k, obconversion.INOPTIONS, str(v))
    if not formatok:
        raise ValueError("%s is not a recognised Open Babel format" % format)
    if not os.path.isfile(filename):
        raise IOError("No such file: '%s'" % filename)

    def filereader():
        obmol = ob.OBMol()
        notatend = obconversion.ReadFile(obmol, filename)
        while notatend:
            yield Molecule(obmol)
            obmol = ob.OBMol()
            notatend = obconversion.Read(obmol)
    return filereader()


def readstring(format, string, opt=None):
    """Read in a molecule from a string.

    Required parameters:
       format - see the informats variable for a list of available
                input formats
       string

    Optional parameters:
       opt    - a dictionary of format-specific options
                For format options with no parameters, specify the
                value as None.

    Example:
    >>> input = "C1=CC=CS1"
    >>> mymol = readstring("smi", input)
    >>> len(mymol.atoms)
    5
    """
    if opt is None:
        opt = {}

    obmol = ob.OBMol()
    obconversion = ob.OBConversion()

    formatok = obconversion.SetInFormat(format)
    if not formatok:
        raise ValueError("%s is not a recognised Open Babel format" % format)
    for k, v in opt.items():
        if v is None:
            obconversion.AddOption(k, obconversion.INOPTIONS)
        else:
            obconversion.AddOption(k, obconversion.INOPTIONS, str(v))

    success = obconversion.ReadString(obmol, string)
    if not success:
        raise IOError("Failed to convert '%s' to format '%s'" % (
            string, format))
    return Molecule(obmol)


class Outputfile(object):
    """Represent a file to which *output* is to be sent.

    Although it's possible to write a single molecule to a file by
    calling the write() method of a molecule, if multiple molecules
    are to be written to the same file you should use the Outputfile
    class.

    Required parameters:
       format - see the outformats variable for a list of available
                output formats
       filename

    Optional parameters:
       overwrite -- if the output file already exists, should it
                   be overwritten? (default is False)
       opt -- a dictionary of format-specific options
              For format options with no parameters, specify the
              value as None.

    Methods:
       write(molecule)
       close()
    """

    def __init__(self, format, filename, overwrite=False, opt=None):
        if opt is None:
            opt = {}
        self.format = format
        self.filename = filename
        if not overwrite and os.path.isfile(self.filename):
            raise IOError(
                "%s already exists. Use 'overwrite=True' to overwrite it." %
                self.filename)

        self.obConversion = ob.OBConversion()
        formatok = self.obConversion.SetOutFormat(self.format)
        if not formatok:
            raise ValueError("%s is not a recognised Open Babel format" %
                             format)
        if filename:
            if isinstance(filename, bytes):
                gzextension = b'.gz'
            else:
                gzextension = '.gz'
            if os.path.splitext(filename)[1] == gzextension:
                self.obconversion.AddOption('z', self.obConversion.GENOPTIONS)
        for k, v in opt.items():
            if v is None:
                self.obConversion.AddOption(k, self.obConversion.OUTOPTIONS)
            else:
                self.obConversion.AddOption(k, self.obConversion.OUTOPTIONS, str(v))
        self.total = 0  # The total number of molecules written to the file

    def write(self, molecule):
        """Write a molecule to the output file.

        Required parameters:
           molecule
        """
        if not self.filename:
            raise IOError("Outputfile instance is closed.")

        if self.total == 0:
            self.obConversion.WriteFile(molecule.OBMol, self.filename)
        else:
            self.obConversion.Write(molecule.OBMol)
        self.total += 1

    def close(self):
        """Close the Outputfile to further writing."""
        self.obConversion.CloseOutFile()
        self.filename = None

    def __enter__(self):
        """Called by with statement, returns itself"""
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Called by with statement, closes itself"""
        self.close()


class Molecule(object):
    """Represent a Pybel Molecule.

    Required parameter:
       OBMol -- an Open Babel OBMol or any type of cinfony Molecule

    Attributes:
       atoms, charge, conformers, data, dim, energy, exactmass, formula,
       molwt, spin, sssr, title, unitcell.
    (refer to the Open Babel library documentation for more info).

    Methods:
       addh(), calcfp(), calcdesc(), draw(), localopt(), make2D(), make3D()
       calccharges(), removeh(), write()

    The underlying Open Babel molecule can be accessed using the attribute:
       OBMol
    """
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

    @property
    def atoms(self):
        return [Atom(self.OBMol.GetAtom(i + 1))
                for i in range(self.OBMol.NumAtoms())]

    @property
    def residues(self):
        return [Residue(res) for res in ob.OBResidueIter(self.OBMol)]

    @property
    def charge(self):
        return self.OBMol.GetTotalCharge()

    @property
    def conformers(self):
        return self.OBMol.GetConformers()

    @property
    def data(self):
        return MoleculeData(self.OBMol)

    @property
    def dim(self):
        return self.OBMol.GetDimension()

    @property
    def energy(self):
        return self.OBMol.GetEnergy()

    @property
    def exactmass(self):
        return self.OBMol.GetExactMass()

    @property
    def formula(self):
        return self.OBMol.GetFormula()

    @property
    def molwt(self):
        return self.OBMol.GetMolWt()

    @property
    def spin(self):
        return self.OBMol.GetTotalSpinMultiplicity()

    @property
    def sssr(self):
        return self.OBMol.GetSSSR()

    def _gettitle(self):
        return self.OBMol.GetTitle()

    def _settitle(self, val):
        self.OBMol.SetTitle(val)
    title = property(_gettitle, _settitle)

    @property
    def unitcell(self):
        unitcell_index = _obconsts.UnitCell
        if sys.platform[:3] == "cli":
            unitcell_index = System.UInt32(unitcell_index)
        unitcell = self.OBMol.GetData(unitcell_index)
        if unitcell:
            if sys.platform[:3] != "cli":
                return _obfuncs.toUnitCell(unitcell)
            else:
                return unitcell.Downcast[ob.OBUnitCell]()
        else:
            raise AttributeError("Molecule has no attribute 'unitcell'")

    @property
    def clone(self):
        return Molecule(ob.OBMol(self.OBMol))

    @property
    def _exchange(self):
        if self.OBMol.HasNonZeroCoords():
            return (1, self.write("mol"))
        else:
            return (0, self.write("can").split()[0])

    def __iter__(self):
        """Iterate over the Atoms of the Molecule.

        This allows constructions such as the following:
           for atom in mymol:
               print atom
        """
        return iter(self.atoms)

    def _repr_svg_(self):
        """For IPython notebook, renders 2D pybel.Molecule SVGs."""

        # Returning None defers to _repr_javascript_
        if ipython_3d:
            return None

        # Open Babel returns a nested svg, which IPython unpacks and treats as
        # two SVGs, messing with the display location. This parses out the
        # inner svg before handing over to IPython.
        namespace = "http://www.w3.org/2000/svg"
        ET.register_namespace("", namespace)
        obsvg = self.clone.write("svg")
        tree = ET.fromstring(obsvg)
        svg = tree.find("{{{ns}}}g/{{{ns}}}svg".format(ns=namespace))
        return ET.tostring(svg).decode("utf-8")

    def _repr_html_(self):
        """For IPython notebook, renders 3D pybel.Molecule webGL objects."""

        # Returning None defers to _repr_svg_
        if not ipython_3d:
            return None

        try:
            import imolecule
        except ImportError:
            raise ImportError("Cannot import 3D rendering. Please install "
                              "with `pip install imolecule`.")
        return imolecule.draw(self.clone, format="pybel", display_html=False)

    def calcdesc(self, descnames=[]):
        """Calculate descriptor values.

        Optional parameter:
           descnames -- a list of names of descriptors

        If descnames is not specified, all available descriptors are
        calculated. See the descs variable for a list of available
        descriptors.
        """
        if not descnames:
            descnames = descs
        ans = {}
        for descname in descnames:
            try:
                desc = _descdict[descname]
            except KeyError:
                raise ValueError(("%s is not a recognised Open Babel "
                                  "descriptor type") % descname)
            ans[descname] = desc.Predict(self.OBMol)
        return ans

    def calcfp(self, fptype="FP2"):
        """Calculate a molecular fingerprint.

        Optional parameters:
           fptype -- the fingerprint type (default is "FP2"). See the
                     fps variable for a list of of available fingerprint
                     types.
        """
        if sys.platform[:3] == "cli":
            fp = ob.VectorUInt()
        else:
            fp = ob.vectorUnsignedInt()
        fptype = fptype.lower()
        try:
            fingerprinter = _fingerprinters[fptype]
        except KeyError:
            raise ValueError(
                "%s is not a recognised Open Babel Fingerprint type" % fptype)
        fingerprinter.GetFingerprint(self.OBMol, fp)
        return Fingerprint(fp)

    def calccharges(self, model="mmff94"):
        """Estimates atomic partial charges in the molecule.

        Optional parameters:
           model -- default is "mmff94". See the charges variable for a list
                    of available charge models (in shell, `obabel -L charges`)

        This method populates the `partialcharge` attribute of each atom
        in the molecule in place.
        """
        model = model.lower()
        try:
            charge_model = _charges[model]
        except KeyError:
            raise ValueError(
                "%s is not a recognised Open Babel Charge Model type" % model)
        success = charge_model.ComputeCharges(self.OBMol)
        if not success:
            errors = ob.obErrorLog.GetMessagesOfLevel(ob.obError)
            error = errors[-1] if errors else "Molecule failed to charge."
            raise Exception(error)
        return [atom.partialcharge for atom in self.atoms]

    def write(self, format="smi", filename=None, overwrite=False, opt=None):
        """Write the molecule to a file or return a string.

        Optional parameters:
           format -- see the informats variable for a list of available
                     output formats (default is "smi")
           filename -- default is None
           overwite -- if the output file already exists, should it
                       be overwritten? (default is False)
           opt -- a dictionary of format specific options
                  For format options with no parameters, specify the
                  value as None.

        If a filename is specified, the result is written to a file.
        Otherwise, a string is returned containing the result.

        To write multiple molecules to the same file you should use
        the Outputfile class.
        """
        if opt is None:
            opt = {}
        obconversion = ob.OBConversion()
        formatok = obconversion.SetOutFormat(format)
        if not formatok:
            raise ValueError("%s is not a recognised Open Babel format" %
                             format)
        if filename:
            if isinstance(filename, bytes):
                gzextension = b'.gz'
            else:
                gzextension = '.gz'
            if os.path.splitext(filename)[1] == gzextension:
                obconversion.AddOption('z', self.obConversion.GENOPTIONS)
        for k, v in opt.items():
            if v is None:
                obconversion.AddOption(k, obconversion.OUTOPTIONS)
            else:
                obconversion.AddOption(k, obconversion.OUTOPTIONS, str(v))

        if filename:
            if not overwrite and os.path.isfile(filename):
                raise IOError(("%s already exists. Use 'overwrite=True' to "
                               "overwrite it.") % filename)
            obconversion.WriteFile(self.OBMol, filename)
            obconversion.CloseOutFile()
        else:
            return obconversion.WriteString(self.OBMol)

    def localopt(self, forcefield="mmff94", steps=500):
        """Locally optimize the coordinates.

        Optional parameters:
           forcefield -- default is "mmff94". See the forcefields variable
                         for a list of available forcefields.
           steps -- default is 500

        If the molecule does not have any coordinates, make3D() is
        called before the optimization. Note that the molecule needs
        to have explicit hydrogens. If not, call addh().
        """
        forcefield = forcefield.lower()
        if self.dim != 3:
            self.make3D(forcefield)
        ff = _forcefields[forcefield]
        success = ff.Setup(self.OBMol)
        if not success:
            return
        ff.SteepestDescent(steps)
        ff.GetCoordinates(self.OBMol)

    def make2D(self):
        """Generate 2D coordinates."""
        _operations['gen2D'].Do(self.OBMol)

    def make3D(self, forcefield="mmff94", steps=50):
        """Generate 3D coordinates.

        Optional parameters:
           forcefield -- default is "mmff94". See the forcefields variable
                         for a list of available forcefields.
           steps -- default is 50

        Once coordinates are generated, hydrogens are added and a quick
        local optimization is carried out with 50 steps and the
        MMFF94 forcefield. Call localopt() if you want
        to improve the coordinates further.
        """
        forcefield = forcefield.lower()
        _builder.Build(self.OBMol)
        self.addh()
        self.localopt(forcefield, steps)

    def addh(self):
        """Add hydrogens."""
        self.OBMol.AddHydrogens()

    def removeh(self):
        """Remove hydrogens."""
        self.OBMol.DeleteHydrogens()

    def convertdbonds(self):
        """Convert Dative Bonds."""
        self.OBMol.ConvertDativeBonds()

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

        Tkinter and Python Imaging Library are required for image display.
        """
        obconversion = ob.OBConversion()
        formatok = obconversion.SetOutFormat("_png2")
        if not formatok:
            raise ImportError("PNG depiction support not found. You should "
                              "compile Open Babel with support for Cairo. See "
                              "installation instructions for more "
                              "information.")

        # Need to copy to avoid removing hydrogens from self
        workingmol = Molecule(ob.OBMol(self.OBMol))
        workingmol.removeh()

        if not usecoords:
            _operations['gen2D'].Do(workingmol.OBMol)
        if update:
            if workingmol.OBMol.NumAtoms() != self.OBMol.NumAtoms():
                raise RuntimeError("It is not possible to update the original "
                                   "molecule with the calculated coordinates, "
                                   "as the original molecule contains "
                                   "explicit hydrogens for which no "
                                   "coordinates have been calculated.")
            else:
                for i in range(workingmol.OBMol.NumAtoms()):
                    self.OBMol.GetAtom(i + 1).SetVector(
                        workingmol.OBMol.GetAtom(i + 1).GetVector())
        if filename:
            filedes = None
        else:
            if sys.platform[:3] == "cli" and show:
                raise RuntimeError("It is only possible to show the molecule "
                                   "if you provide a filename. The reason for "
                                   "this is that I kept having problems "
                                   "when using temporary files.")

            filedes, filename = tempfile.mkstemp()

        workingmol.write("_png2", filename=filename, overwrite=True)

        if show:
            if sys.platform[:4] == "java":
                image = javax.imageio.ImageIO.read(java.io.File(filename))
                frame = javax.swing.JFrame(visible=1)
                frame.getContentPane().add(
                    javax.swing.JLabel(javax.swing.ImageIcon(image)))
                frame.setSize(300, 300)
                frame.setDefaultCloseOperation(
                    javax.swing.WindowConstants.DISPOSE_ON_CLOSE)
                frame.show()
            elif sys.platform[:3] == "cli":
                form = _MyForm()
                form.setup(filename, self.title)
                Application.Run(form)
            else:
                if not tk:
                    raise ImportError("Tkinter or Python Imaging Library not "
                                      "found, but is required for image "
                                      "display. See installation instructions "
                                      "for more information.")
                root = tk.Tk()
                root.title((hasattr(self, "title") and self.title)
                           or self.__str__().rstrip())
                frame = tk.Frame(root, colormap="new",
                                 visual='truecolor').pack()
                image = PIL.open(filename)
                imagedata = piltk.PhotoImage(image)
                tk.Label(frame, image=imagedata).pack()
                tk.Button(root, text="Close", command=root.destroy).pack(
                    fill=tk.X)
                root.mainloop()
        if filedes:
            os.close(filedes)
            os.remove(filename)


class Atom(object):
    """Represent a Pybel atom.

    Required parameter:
       OBAtom -- an Open Babel OBAtom

    Attributes:
       atomicmass, atomicnum, cidx, coords, coordidx, degree, exactmass,
       formalcharge, heavydegree, heterodegree, hyb, idx,
       implicitvalence, isotope, partialcharge, residue, spin, type,
       vector.

    (refer to the Open Babel library documentation for more info).

    The original Open Babel atom can be accessed using the attribute:
       OBAtom
    """

    def __init__(self, OBAtom):
        self.OBAtom = OBAtom

    @property
    def coords(self):
        return (self.OBAtom.GetX(), self.OBAtom.GetY(), self.OBAtom.GetZ())

    @property
    def atomicmass(self):
        return self.OBAtom.GetAtomicMass()

    @property
    def atomicnum(self):
        return self.OBAtom.GetAtomicNum()

    @property
    def cidx(self):
        return self.OBAtom.GetCIdx()

    @property
    def coordidx(self):
        return self.OBAtom.GetCoordinateIdx()

    @property
    def degree(self):
        return self.OBAtom.GetExplicitDegree()

    @property
    def exactmass(self):
        return self.OBAtom.GetExactMass()

    @property
    def formalcharge(self):
        return self.OBAtom.GetFormalCharge()

    @property
    def heavydegree(self):
        return self.OBAtom.GetHvyDegree()
    
    @property
    def heavyvalence(self):
        raise AttributeError("This property has been renamed. Use Atom.heavydegree instead.")

    @property
    def heterodegree(self):
        return self.OBAtom.GetHeteroDegree()
    
    @property
    def heterovalence(self):
        raise AttributeError("This property has been renamed. Use Atom.heterodegree instead.")

    @property
    def hyb(self):
        return self.OBAtom.GetHyb()

    @property
    def idx(self):
        return self.OBAtom.GetIdx()

    @property
    def implicitvalence(self):
        return self.OBAtom.GetImplicitValence()

    @property
    def isotope(self):
        return self.OBAtom.GetIsotope()

    @property
    def partialcharge(self):
        return self.OBAtom.GetPartialCharge()

    @property
    def residue(self):
        return Residue(self.OBAtom.GetResidue())

    @property
    def spin(self):
        return self.OBAtom.GetSpinMultiplicity()

    @property
    def type(self):
        return self.OBAtom.GetType()

    @property
    def valence(self):
        raise AttributeError("This property has been renamed. Use Atom.degree instead.")

    @property
    def vector(self):
        return self.OBAtom.GetVector()

    def __str__(self):
        c = self.coords
        return "Atom: %d (%.2f %.2f %.2f)" % (self.atomicnum, c[0], c[1], c[2])


class Residue(object):
    """Represent a Pybel residue.

    Required parameter:
       OBResidue -- an Open Babel OBResidue

    Attributes:
       atoms, idx, name.

    (refer to the Open Babel library documentation for more info).

    The original Open Babel atom can be accessed using the attribute:
       OBResidue
    """

    def __init__(self, OBResidue):
        self.OBResidue = OBResidue

    @property
    def atoms(self):
        return [Atom(atom) for atom in ob.OBResidueAtomIter(self.OBResidue)]

    @property
    def idx(self):
        return self.OBResidue.GetIdx()

    @property
    def name(self):
        return self.OBResidue.GetName()

    def __iter__(self):
        """Iterate over the Atoms of the Residue.

        This allows constructions such as the following:
           for atom in residue:
               print atom
        """
        return iter(self.atoms)


def _findbits(fp, bitsperint):
    """Find which bits are set in a list/vector.

    This function is used by the Fingerprint class.

    >>> _findbits([13, 71], 8)
    [1, 3, 4, 9, 10, 11, 15]
    """
    ans = []
    start = 1
    if sys.platform[:4] == "java":
        fp = [fp.get(i) for i in range(fp.size())]
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
       fingerprint -- a vector calculated by OBFingerprint.FindFingerprint()

    Attributes:
       fp -- the underlying fingerprint object
       bits -- a list of bits set in the Fingerprint

    Methods:
       The "|" operator can be used to calculate the Tanimoto coeff. For
       example, given two Fingerprints 'a', and 'b', the Tanimoto coefficient
       is given by:
          tanimoto = a | b
    """

    def __init__(self, fingerprint):
        self.fp = fingerprint

    def __or__(self, other):
        return ob.OBFingerprint.Tanimoto(self.fp, other.fp)

    @property
    def bits(self):
        return _findbits(self.fp, ob.OBFingerprint.Getbitsperint())

    def __str__(self):
        fp = self.fp
        if sys.platform[:4] == "java":
            fp = [self.fp.get(i) for i in range(self.fp.size())]
        return ", ".join([str(x) for x in fp])


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

    def __init__(self, smartspattern):
        """Initialise with a SMARTS pattern."""
        self.obsmarts = ob.OBSmartsPattern()
        success = self.obsmarts.Init(smartspattern)
        if not success:
            raise IOError("Invalid SMARTS pattern")

    def findall(self, molecule):
        """Find all matches of the SMARTS pattern to a particular molecule.

        Required parameters:
           molecule
        """
        self.obsmarts.Match(molecule.OBMol)
        vector = self.obsmarts.GetUMapList()
        if sys.platform[:4] == "java":
            vector = [vector.get(i) for i in range(vector.size())]
        return list(vector)


class MoleculeData(object):
    """Store molecule data in a dictionary-type object

    Required parameters:
      obmol -- an Open Babel OBMol

    Methods and accessor methods are like those of a dictionary except
    that the data is retrieved on-the-fly from the underlying OBMol.

    Example:
    >>> mol = readfile("sdf", 'head.sdf').next() # Python 2
    >>> # mol = next(readfile("sdf", 'head.sdf')) # Python 3
    >>> data = mol.data
    >>> print data
    {'Comment': 'CORINA 2.61 0041  25.10.2001', 'NSC': '1'}
    >>> print len(data), data.keys(), data.has_key("NSC")
    2 ['Comment', 'NSC'] True
    >>> print data['Comment']
    CORINA 2.61 0041  25.10.2001
    >>> data['Comment'] = 'This is a new comment'
    >>> for k,v in data.items():
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
        data = self._mol.GetData()
        if sys.platform[:4] == "java":
            data = [data.get(i) for i in range(data.size())]
        answer = [x for x in data if
                  x.GetDataType() == _obconsts.PairData or
                  x.GetDataType() == _obconsts.CommentData]
        if sys.platform[:3] != "cli":
            answer = [_obfuncs.toPairData(x) for x in answer]
        return answer

    def _testforkey(self, key):
        if key not in self:
            raise KeyError("'%s'" % key)

    def keys(self):
        return [x.GetAttribute() for x in self._data()]

    def values(self):
        return [x.GetValue() for x in self._data()]

    def items(self):
        return iter(zip(self.keys(), self.values()))

    def __iter__(self):
        return iter(self.keys())

    def iteritems(self):  # Can remove for Python 3
        return self.items()

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
        for k, v in dictionary.items():
            self[k] = v

    def __getitem__(self, key):
        self._testforkey(key)
        answer = self._mol.GetData(key)
        if sys.platform[:3] != "cli":
            answer = _obfuncs.toPairData(answer)
        return answer.GetValue()

    def __setitem__(self, key, value):
        if key in self:
            if sys.platform[:3] != "cli":
                pairdata = _obfuncs.toPairData(self._mol.GetData(key))
            else:
                pairdata = self._mol.GetData(key).Downcast[ob.OBPairData]()
            pairdata.SetValue(str(value))
        else:
            pairdata = ob.OBPairData()
            pairdata.SetAttribute(key)
            pairdata.SetValue(str(value))
            self._mol.CloneData(pairdata)

    def __repr__(self):
        return dict(self.items()).__repr__()

if sys.platform[:3] == "cli":
    class _MyForm(Form):

        def __init__(self):
            Form.__init__(self)

        def setup(self, filename, title):
            # adjust the form's client area size to the picture
            self.ClientSize = Size(300, 300)
            self.Text = title

            self.filename = filename
            self.image = Image.FromFile(self.filename)
            pictureBox = PictureBox()
            # this will fit the image to the form
            pictureBox.SizeMode = PictureBoxSizeMode.StretchImage
            pictureBox.Image = self.image
            # fit the picture box to the frame
            pictureBox.Dock = DockStyle.Fill

            self.Controls.Add(pictureBox)
            self.Show()

if __name__ == "__main__":  # pragma: no cover
    import doctest
    doctest.testmod(verbose=True)
