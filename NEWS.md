# Open Babel 2.4.0 (2016-9-21)

This release represents a major update and should be a stable upgrade,
strongly recommended for all users.

Note that this release deprecates the babel executable in favor of obabel. A future release will remove babel entirely. For information on the differences, please see [the documentation](http://openbabel.org/docs/current/Command-line_tools/babel.html).

## New file formats

* DALTON output files (read only) and DALTON input files (read/write) (Casper Steinmann)
* JSON format used by ChemDoodle (read/write) (Matt Swain)
* JSON format used by PubChem (read/write) (Matt Swain)
* LPMD's atomic configuration file (read/write) (Joaquin Peralta)
* The format used by the CONTFF and POSFF files in MDFF (read/write) (Kirill Okhotnikov)
* ORCA output files (read only) and ORCA input files (write only) (Dagmar Lenk)
* ORCA-AICCM's extended XYZ format (read/write) (Dagmar Lenk)
* Painter format for custom 2D depictions (write only) (Noel O'Boyle)
* Siesta output files (read only) (Patrick Avery)
* Smiley parser for parsing SMILES according to the OpenSMILES specification (read only) (Tim Vandermeersch)
* STL 3D-printing format (write only) (Matt Harvey)
* Turbomole AOFORCE output (read only) (Mathias Laurin)
* A representation of the VDW surface as a point cloud (write only) (Matt Harvey)

## New file format capabilities and options

* AutoDock PDBQT: Options to preserve hydrogens and/or atom names (Matt Harvey)
* CAR: Improved space group support in .car files (kartlee)
* CDXML: Read/write isotopes (Roger Sayle)
* CIF: Extract charges (Kirill Okhotnikov)
* CIF: Improved support for space-groups and symmetries (Alexandr Fonari)
* DL_Poly: Cell information is now read (Kirill Okhotnikov)
* Gaussian FCHK: Parse alpha and beta orbitals (Geoff Hutchison)
* Gaussian out: Extract true enthalpy of formation, quadrupole, polarizability tensor, electrostatic potential fitting points and potential values, and more (David van der Spoel)
* MDL Mol: Read in atom class information by default and optionally write it out (Roger Sayle)
* MDL Mol: Support added for ZBO, ZCH and HYD extensions (Matt Swain)
* MDL Mol: Implement the MDL valence model on reading (Roger Sayle)
* MDL SDF: Option to write out an ASCII depiction as a property (Noel O'Boyle)
* mmCIF: Improved mmCIF reading (Patrick Fuller)
* mmCIF: Support for atom occupancy and atom_type (Kirill Okhotnikov)
* Mol2: Option to read UCSF Dock scores (Maciej Wójcikowski)
* MOPAC: Read z-matrix data and parse (and prefer) ESP charges (Geoff Hutchison)
* NWChem: Support sequential calculations by optionally overwriting earlier ones (Dmitriy Fomichev)
* NWChem: Extract info on MEP(IRC), NEB and quadrupole moments (Dmitriy Fomichev)
* PDB: Read/write PDB insertion codes (Steffen Möller)
* PNG: Options to crop the margin, and control the background and bond colors (Fredrik Wallner)
* PQR: Use a stored atom radius (if present) in preference to the generic element radius (Zhixiong Zhao)
* PWSCF: Extend parsing of lattice vectors (David Lonie)
* PWSCF: Support newer versions, and the 'alat' term (Patrick Avery)
* SVG: Option to avoid addition of hydrogens to fill valence (Lee-Ping)
* SVG: Option to draw as ball-and-stick (Jean-Noël Avila)
* VASP: Vibration intensities are calculated (Christian Neiss, Mathias Laurin)
* VASP: Custom atom element sorting on writing (Kirill Okhotnikov)

## Other new features and improvements

* 2D layout: Improved the choice of which bonds to designate as hash/wedge bonds around a stereo center (Craig James)
* 3D builder: Use bond length corrections based on bond order from Pyykko and Atsumi (http://dx.doi.org/10.1002/chem.200901472) (Geoff Hutchison)
* 3D generation: "--gen3d", allow user to specify the desired speed/quality (Geoff Hutchison)
* Aromaticity: Improved detection (Geoff Hutchison)
* Canonicalisation: Changed behaviour for multi-molecule SMILES. Now each molecule is canonicalized individually and then sorted. (Geoff Hutchison/Tim Vandermeersch)
* Charge models: "--print" writes the partial charges to standard output after calculation (Geoff Hutchison)
* Conformations: Confab, the systematic conformation generator, has been incorporated into Open Babel (David Hall/Noel O'Boyle)
* Conformations: Initial support for ring rotamer sampling (Geoff Hutchison)
* Conformer searching: Performance improvement by avoiding gradient calculation and optimising the default parameters (Geoff Hutchison)
* EEM charge model: Extend to use additional params from http://dx.doi.org/10.1186/s13321-015-0107-1 (Tomáš Raček)
* FillUnitCell operation: Improved behavior (Patrick Fuller)
* Find duplicates: The "--duplicate" option can now return duplicates instead of just removing them (Chris Morley)
* GAFF forcefield: Atom types updated to match Wang et al. J. Comp. Chem. 2004, 25, 1157 (Mohammad Ghahremanpour)
* New charge model: EQeq crystal charge equilibration method (a speed-optimized crystal-focused charge estimator, http://pubs.acs.org/doi/abs/10.1021/jz3008485) (David Lonie)
* New charge model: "fromfile" reads partial charges from a named file (Matt Harvey)
* New conversion operation: "changecell", for changing cell dimensions (Kirill Okhotnikov)
* New command-line utility: "obthermo", for extracting thermochemistry data from QM calculations (David van der Spoel)
* New fingerprint: ECFP (Geoff Hutchison/Noel O'Boyle/Roger Sayle)
* OBConversion: Improvements and API changes to deal with a long-standing memory leak (David Koes)
* OBAtom::IsHBondAcceptor(): Definition updated to take into account the atom environment (Stefano Forli)
* Performance: Faster ring-finding algorithm (Roger Sayle)
* Performance: Faster fingerprint similarity calculations if compiled with -DOPTIMIZE_NATIVE=ON (Noel O'Boyle/Jeff Janes)
* SMARTS matching: The "-s" option now accepts an integer specifying the number of matches required (Chris Morley)
* UFF: Update to use traditional Rappe angle potential (Geoff Hutchison)

## Language bindings

* Bindings: Support compiling only the bindings against system libopenbabel (Reinis Danne)
* Java bindings: Add example Scala program using the Java bindings (Reinis Danne)
* New bindings: PHP (Maciej Wójcikowski)
* PHP bindings: BaPHPel, a simplified interface (Maciej Wójcikowski)
* Python bindings: Add 3D depiction support for Jupyter notebook  (Patrick Fuller)
* Python bindings, Pybel: calccharges() and convertdbonds() added (Patrick Fuller, Björn Grüning)
* Python bindings, Pybel: compress output if filename ends with .gz (Maciej Wójcikowski)
* Python bindings, Pybel: Residue support (Maciej Wójcikowski)

## Development/Build/Install Improvements

* Version control: move to git and GitHub from subversion and SourceForge
* Continuous integration: Travis for Linux builds and Appveyor for Windows builds (David Lonie and Noel O'Boyle)
* Python installer: Improvements to the Python setup.py installer and "pip install openbabel" (David Hall, Matt Swain, Joshua Swamidass)
* Compilation speedup: Speed up compilation by combining the tests (Noel O'Boyle)
* MacOSX: Support compiling with libc++ on MacOSX (Matt Swain)

## Cast of contributors

Alexandr Fonari, Anders Steen Christensen, Andreas Kempe, arkose, Benoit Leblanc, Björn Grüning, Casper Steinmann, Chris Morley, Christoph Willing, Craig James, Dagmar Lenk, David Hall, David Koes, David Lonie, David van der Spoel, Dmitriy Fomichev, Fulvio Ciriaco, Fredrik Wallner, Geoff Hutchison, Heiko Becker, Itay Zandbank, Jean-Noel Avila, Jeff Janes, Joaquin Peralta, Joshua Swamidass, Julien Nabet, Karol Langner, Karthik Rajagopalan, Katsuhiko Nishimra, Kevin Horan, Kirill Okhotnikov, Lee-Ping, Matt Harvey, Maciej Wójcikowski, Marcus Hanwell, Mathias Laurin, Matt Swain, Mohamad Mohebifar, Mohammad Ghahremanpour, Noel O'Boyle, Patrick Avery, Patrick Fuller, Paul van Maaren, Peng Bai, Philipp Thiel, Reinis Danne, Ronald Cohen, Scott McKechnie, Stefano Forli, Steve Roughley, Steffen Moeller, Tim Vandermeersch, Tomas Racek, Tomáš Trnka, Tor Colvin, Torsten Sachse, Yi-Shu Tu, Zhixiong Zhao
