/**********************************************************************
Copyright (C) 2010 by Reinis Danne
Some portions Copyright (C) 2004 by Chris Morley

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/elements.h>
#include <openbabel/generic.h>
#include <openbabel/obiter.h>


#include <sstream>
#include <iomanip>

using namespace std;
namespace OpenBabel
{

class GROFormat : public OBMoleculeFormat
{
public:
  //Register this format type ID in the constructor
  GROFormat()
  {
    /* GRO is the file extension and is case insensitive. A MIME type can be
       added as an optional third parameter.
       Multiple file extensions can be registered by adding extra statements.*/
    OBConversion::RegisterFormat("gro",this);

    /* If there are any format specific options they should be registered here
       so that the commandline interface works properly.
       The first parameter is the option name. If it is a single letter it can
       be concatinated with other single letter options. For output options it
       can be multicharcter and is then written as --optionname on the command
       line. The third parameter is the number of parameters the option takes.
       Currently this is either 1 or 0 and if it is 0 can be omitted for output
       options. The parameter is always text and needs to be parsed to extract
       a number.

       Options can apply when writing - 4th parameter is
       OBConversion::OUTOPTIONS or can be omitted as shown. A single letter
       output option is preceded by -x on the command line.
       Or options can apply to the input format - the 4th parameter is
       OBConversion::INOPTIONS. They are then preceded by -a on the command
       line.

       Each option letter may be reused in other formats, but within the same
       group, INOPTIONS or OUTOPTIONS, must take the same number of parameters
       (0 or 1). There will be an error message when OpenBabel  runs if there
       are conflicts between formats. A list of formats currently used (which
       may not be comprehensive) is in docs/options.html.
    */
    OBConversion::RegisterOptionParam("f", this, 1);
    OBConversion::RegisterOptionParam("n", this);
    OBConversion::RegisterOptionParam("s", this, 0, OBConversion::INOPTIONS);

  }

  /* The first line of the description should be a brief identifier, <40 chars,
     because it is used in dropdown lists, etc. in some user interfaces.
     The rest is optional.

     Describe any format specific options here. This text is parsed to provide
     checkboxes, etc for the GUI (for details click the control menu),
     so please try to keep to a similar form.

     Write options are the most common, and the "Write" is optional.
     The option f takes a text parameter, so that it is essential that the
     option is registered in the constructor of the class.
     Finish the options with a blank line as shown, if there are more than one
     group of options, or if there are further comments after them.
  */
  virtual const char* Description() //required
  {
    return
    "GRO format\n"
    "This is GRO file format as used in Gromacs.\n"
    "Right now there is only limited support for element perception. It works "
    "for \nelements with one letter symbols if the atomtype starts with the "
    "same letter.\n\n"

    "Read Options e.g. -as\n"
    " s  Consider single bonds only\n"
    " b  Disable bonding entierly\n"
    ;
  };

  //Optional URL where the file format is specified
  virtual const char* SpecificationURL()
  {
    return "http://manual.gromacs.org/documentation/current/reference-manual/file-formats.html#gro";
  }

  /* This optional function is for formats which can contain more than one
     molecule. It is used to quickly position the input stream after the nth
     molecule without requiring to convert and discard all the n molecules.
     See obconversion.cpp for details.*/
  virtual int SkipObjects(int n, OBConversion* pConv)
  {
    string line = "";
    int natoms = 0;
    int nlines = 0;

    if (n == 0)
      n++;
    istream& ifs = *pConv->GetInStream();
    getline(ifs, line);
    ifs >> natoms;
    nlines = (natoms+3)*n;
    while (ifs && --nlines) {
      getline(ifs, line);
    }

    return ifs.good() ? 1 : -1;
  };

  ////////////////////////////////////////////////////
  /// Declarations for the "API" interface functions. Definitions are below
  virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
  virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
};
  ////////////////////////////////////////////////////
//Make an instance of the format class
GROFormat theGROFormat;

/////////////////////////////////////////////////////////////////

bool GROFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{
  OBMol* pmol = pOb->CastAndClear<OBMol>();
  if(pmol==NULL)
      return false;

  istream& ifs = *pConv->GetInStream();

  pmol->BeginModify();

  /** Parse the input stream and use the OpenBabel API to populate the OBMol **/

  char buffer[BUFF_SIZE];
  string line = "";
  stringstream errorMsg;

  int natoms = 0;
  string title = "";
  long int resid = 0; // 5
  string resname = ""; //5
  string atomtype = ""; //5
  //long int atomid = 0; //5
  double x = 0.0; // 8.3
  double y = 0.0; // 8.3
  double z = 0.0; // 8.3
  double vx = 0.0; // 8.4
  double vy = 0.0; // 8.4
  double vz = 0.0; // 8.4
  string tempstr = "";
  long int residx = 0;
  OBAtom* atom;
  OBResidue* res;
  OBVectorData* velocity;

  if (!ifs || ifs.peek() == EOF) {
    return false; // Trying to read past end of the file
  }

  if (!ifs.getline(buffer, BUFF_SIZE)) {
    errorMsg << "Problems reading a GRO file: "
             << "Cannot read the first line!";
    obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
    return false;
  }

  // Get the title
  title.assign(buffer);
  if (title.size() < 1) {
    title = pConv->GetTitle();
    pmol->SetTitle(title);
  } else {
    pmol->SetTitle(title);
  }

  if (!ifs.getline(buffer, BUFF_SIZE)) {
    errorMsg << "Problems reading a GRO file: "
             << "Cannot read the second line!";
    obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
    return false;
  }

  // Get the number of atoms
  stringstream(buffer) >> natoms;
  if (natoms < 1) {
    errorMsg << "Problems reading a GRO file: "
             << "There are no atoms in the file or the second line is"
             << " incorrectly written.";
    obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
    return false;
  }
  pmol->ReserveAtoms(natoms);

  // Read all atom records
  for (int i=1; i<=natoms; i++) {
    if (!ifs.getline(buffer,BUFF_SIZE)) {
      errorMsg << "Problems reading a GRO file: "
               << "Could not read line #" << i+2 << ", file error." << endl
               << " According to the second line, there should be " << natoms
               << " atoms, and therefore " << natoms+3 << " lines in the file.";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      return false;
    }

    line = buffer;

    // Get atom
    atom  = pmol->NewAtom();

    tempstr.assign(line,0,5);
    stringstream(tempstr) >> resid;

    resname.assign(line,5,5);
    Trim(resname);

    atomtype.assign(line,10,5);
    Trim(atomtype);

    // Not used, OB assigns its own indizes
    //tempstr.assign(line,15,5);
    //stringstream(tempstr) >> atomid;

    tempstr.assign(line,20,8);
    stringstream(tempstr) >> x;

    tempstr.assign(line,28,8);
    stringstream(tempstr) >> y;

    tempstr.assign(line,36,8);
    stringstream(tempstr) >> z;

    if (line.size() > 44) {
      tempstr.assign(line,44,8);
      stringstream(tempstr) >> vx;

      tempstr.assign(line,52,8);
      stringstream(tempstr) >> vy;

      tempstr.assign(line,60,8);
      stringstream(tempstr) >> vz;

      velocity = new OBVectorData();
      velocity->SetData(vx, vy, vz);
      velocity->SetAttribute("Velocity");
      velocity->SetOrigin(fileformatInput);
      atom->SetData(velocity);
    }

    // OB translates this and, e.g., OW1 turns into O3
    // Type conversion should be done explicitly if that needs to be
    // controlled.
    atom->SetType(atomtype);

    // Set coordinates of the atom, multiply by 10 to convert from nm to
    // angstrom
    atom->SetVector(x*10, y*10, z*10);

    if (resid == residx) {
      // Add atom to an existing residue
      res->AddAtom(atom);
    } else {
      // Create new residue and use that
      res = pmol->NewResidue();
      res->SetName(resname);
      res->SetNum(resid);
      res->AddAtom(atom);
      residx = resid;
    }

    // Atom type has to be set in residues as AtomID
    res->SetAtomID(atom, atomtype);

    // For determining the element this doesn't work:
    // atom->SetAtomicNum(OBElements::GetAtomicNum(atom->GetType()));
    //
    // So the element simbol should be found while reading in the file. It
    // could be possible to provide an option for external atomtype<->element
    // simbol map. Now such functionality is not implemented.

    // TODO: Make element perception optional and improve it
    string element = "";
    if (atomtype[0] == 'C') {
      if (atomtype.find_first_of("aloru",1) == 1) {
        element.assign(atomtype,0,2);
      } else {
        element.assign(atomtype,0,1);
      }
    } else if (atomtype[0] == 'N') {
      if (atomtype.find_first_of("abei",1) == 1) {
        element.assign(atomtype,0,2);
      } else {
        element.assign(atomtype,0,1);
      }
    } else {
      element.assign(atomtype,0,1);
    }
    atom->SetAtomicNum(OBElements::GetAtomicNum(element.data()));
  }

  // Get periodic box
  if (!ifs.getline(buffer,BUFF_SIZE)) {
    errorMsg << "Problems reading a GRO file: "
             << "Could not read box vectors!" << endl;
    obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
    return false;
  }

  double v1x=0.0, v2y=0.0, v3z=0.0;
  double v1y=0.0, v1z=0.0, v2x=0.0;
  double v2z=0.0, v3x=0.0, v3y=0.0;

  stringstream ss(buffer);
  if (!ss) {
    errorMsg << "Problems reading a GRO file: "
             << "Could not read box vectors!" << endl;
    obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
    return false;
  }

  ss >> v1x >> v2y >> v3z;
  if (ss) {
    ss >> v1y >> v1z >> v2x;
    ss >> v2z >> v3x >> v3y;
  }

  if (!(v1x == 0.0 && v2y == 0.0 && v3z == 0.0 &&
        v1y == 0.0 && v1z == 0.0 && v2x == 0.0 &&
        v2z == 0.0 && v3x == 0.0 && v3y == 0.0)) {
  // Set box vectors and convert nm to angstroms
  const vector3 v1(v1x*10,v1y*10,v1z*10);
  const vector3 v2(v2x*10,v2y*10,v2z*10);
  const vector3 v3(v3x*10,v3y*10,v3z*10);

  OBUnitCell* cell = new OBUnitCell();
  cell->SetData(v1, v2, v3);
  pmol->SetData(cell);
  }

  pmol->EndModify();

  // Input options -as and -ab
  if (!pConv->IsOption("b", OBConversion::INOPTIONS)) {
    pmol->ConnectTheDots();
    if (!pConv->IsOption("s", OBConversion::INOPTIONS)) {
      pmol->PerceiveBondOrders();
    }
  }
  pmol->SetChainsPerceived();

  /* For multi-molecule formats, leave the input stream at the start of the
     next molecule, ready for this routine to be called again.

   * Return true if ok. Returning false means discard the OBMol and stop
     converting, unless the -e option is set. With a multi-molecule inputstream
     this will skip the current molecule and continue with the next, if
     SkipObjects() has been defined. If it has not, and continuation after
     errors is still required, it is necessary to leave the input stream at the
     beginning of next object when returning false;*/
  return true;
}

////////////////////////////////////////////////////////////////

bool GROFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
  OBMol* pmol = dynamic_cast<OBMol*>(pOb);
  if(pmol==NULL)
      return false; // Stop converting

  ostream& ofs = *pConv->GetOutStream();

  /** Write the representation of the OBMol molecule to the output stream **/

  OBVectorData* vector;
  vector3 v;
  OBResidue* res;
  string tempstr = "";
  long int atIdx = 0;
  long int resIdx = 0;

  ofs << pmol->GetTitle() << endl;
  ofs << pmol->NumAtoms() << endl;
  ofs.setf(ios::fixed);

  FOR_ATOMS_OF_MOL(atom, pmol) {
    res = atom->GetResidue();
    resIdx = res->GetNum();
    // Check if residue index excedes the field width and should be wrapped
    if (resIdx > 99999) {
      ofs << setw(5) << resIdx - 100000*int(resIdx/100000);
    } else {
      ofs << setw(5) << res->GetNum();
    }
    ofs << setw(5) << left  << res->GetName();
    // Remove whitespace from AtomID left by other formats
    tempstr = res->GetAtomID(&(*atom));
    ofs << setw(5) << right << Trim(tempstr);
    atIdx = atom->GetIdx();
    // Check if atom index excedes the field width and should be wrapped
    if (atIdx > 99999) {
      ofs << setw(5) << atIdx - 100000*int(atIdx/100000);
    } else {
      ofs << setw(5) << atIdx;
    }
    ofs.precision(3);
    // OpenBabel uses angstrom, so converting to nm by dividing by 10
    ofs << setw(8) << atom->x()/10
        << setw(8) << atom->y()/10
        << setw(8) << atom->z()/10;
    if (atom->GetData("Velocity")) {
      vector = (OBVectorData*) atom->GetData("Velocity");
      v = vector->GetData();
      ofs.precision(4);
      ofs << setw(8) << v.x()
          << setw(8) << v.y()
          << setw(8) << v.z();
    }
    ofs << endl;
  }

  // On the last line of the file goes periodic box specification
  if (pmol->HasData(OBGenericDataType::UnitCell)) {
    OBUnitCell* cell = (OBUnitCell*) pmol->GetData(OBGenericDataType::UnitCell);
    matrix3x3 m = cell->GetCellMatrix();
    vector3 v1 = m.GetRow(0);
    vector3 v2 = m.GetRow(1);
    vector3 v3 = m.GetRow(2);

    // Gromacs itself uses precision of 5, so it should be fine
    ofs.precision(5);
    ofs << "   " << v1.x()/10
        << "   " << v2.y()/10
        << "   " << v3.z()/10;

    // If there is any non-zero value among others, then write them all
    const double TRESHOLD = 1.0e-8;
    if (fabs(v1.y()) > TRESHOLD || fabs(v1.z()) > TRESHOLD ||
        fabs(v2.x()) > TRESHOLD || fabs(v2.z()) > TRESHOLD ||
        fabs(v3.x()) > TRESHOLD || fabs(v3.y()) > TRESHOLD) {
      ofs << "   " << v1.y()/10
          << "   " << v1.z()/10
          << "   " << v2.x()/10
          << "   " << v2.z()/10
          << "   " << v3.x()/10
          << "   " << v3.y()/10;
    }
  } else {
    // Set to zero if there is no box data in the molecule
    ofs << "   0.00000   0.00000   0.00000";
  }
  ofs << endl;

  return true;
}

} //namespace OpenBabel

