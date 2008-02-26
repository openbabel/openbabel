/**********************************************************************
gausscubeformat.cpp - Read in Gaussian cube format files.

// Molekel - Molecular Visualization Program
// Copyright (C) 2006, 2007 Swiss National Supercomputing Centre (CSCS)

 Some Portions Copyright (c) 2007 by Geoffrey R. Hutchison
 Some Portions Copyright (C) 2008 by Marcus D. Hanwell

This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation: either version 2 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include <sstream>

#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/griddata.h>

using namespace std;

/// @warning seems that every position in cube files has always to be
/// multiplied by this constant even if according to:
/// reference: http://www.gaussian.com/g_ur/u_cubegen.htm
/// <Num points along first axis> > 0 means that the values
/// are already in Angstroms.
static const double BOHR_TO_ANGSTROM = 0.529177249;

namespace OpenBabel
{

  class OBGaussianCubeFormat : public OpenBabel::OBMoleculeFormat
  {
  public:
    // Constructor: register 'cube' and "CUBE" format.
    OBGaussianCubeFormat()
    {
        OpenBabel::OBConversion::RegisterFormat( "cube", this );
        OpenBabel::OBConversion::RegisterFormat( "CUBE", this );
        OpenBabel::OBConversion::RegisterFormat( "cub", this );
    }

    // Return description.
    virtual const char* Description() //required
    {
        return
        "Gaussian cube format\n"
        "Read only.\n"
        "b no bonds\n"
        "s no multiple bonds\n\n";
    }

    // Return a specification url, not really a specification since
    // I couldn't find it but close enough.
    virtual const char* SpecificationURL()
    {
        return "http://www.gaussian.com/g_ur/u_cubegen.htm";
    }

    // Return MIME type, NULL in this case.
    virtual const char* GetMIMEType() { return 0; };

    // Return read/write flag: read only.
    virtual unsigned int Flags()
    {
        return READONEONLY;
    };

    // Skip to object: used for multi-object file formats.
    virtual int SkipObjects( int n, OpenBabel::OBConversion* pConv ) { return 0; }

    /// The "API" interface functions
    virtual bool ReadMolecule( OpenBabel::OBBase* pOb, OpenBabel::OBConversion* pConv );
    /// Write: always returns false right now - read only
    virtual bool WriteMolecule( OpenBabel::OBBase* , OpenBabel::OBConversion* )
    {
        return false;
    }
};

//------------------------------------------------------------------------------

    // Global variable used to register Gaussian cube format.
    OBGaussianCubeFormat theGaussianCubeFormat;

//------------------------------------------------------------------------------
bool OBGaussianCubeFormat::ReadMolecule( OBBase* pOb, OBConversion* pConv )
{
    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol == 0)
      return false;

    istream& ifs = *pConv->GetInStream();

    const char* title = pConv->GetTitle();
    char buffer[BUFF_SIZE];

    stringstream errorMsg;

    int nAtoms = 0;
    int negAtoms = false;

    int line = 0; // Line counter - help in debugging...

    if (!ifs)
      return false; // We are attempting to read past the end of the file

    // The first two lines are comments - and so not needed.
    ++line;
    if (!ifs.getline(buffer, BUFF_SIZE))
    {
      obErrorLog.ThrowError(__FUNCTION__,
                            "Problem reading the Gaussian cube file: cannot read the first line (title/comments).", obWarning);
      return false;
    }
    ++line;
    if (!ifs.getline(buffer,BUFF_SIZE))
    {
      obErrorLog.ThrowError(__FUNCTION__,
                            "Problem reading the Gaussian cube file: cannot read the second line (title/comments).", obWarning);
      return false;
    }

    // This line contains the number of atoms included in the file followed by
    // the position of the origin of the cube.
    ++line;
    if (!ifs.getline(buffer, BUFF_SIZE))
    {
      obErrorLog.ThrowError(__FUNCTION__,
                            "Problem reading the Gaussian cube file: cannot read the third line. This should contain the number of atoms and the position of the origin of the cube.", obWarning);
      return false;
    }
    vector<string> vs;
    tokenize(vs, buffer);
    if (vs.size() < 4)
    {
      obErrorLog.ThrowError(__FUNCTION__,
                            "Problem reading the Gassian cube file: The third line must contain the number of atoms and the origin of the cube.", obWarning);
      return false;
    }
    char *endptr;
    nAtoms = strtol(static_cast<const char*>(vs.at(0).c_str()), &endptr, 10);
    if (endptr == static_cast<const char*>(vs.at(0).c_str()))
    {
      errorMsg << "Problems reading the Gaussian cube file: "
               << "Could not read line #3.\n"
               << "According to the specification this line should contain "
               << "four entries separated by white space.\n"
               << "OpenBabel could not interpret item #0 as an integer.";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      return false;
    }
    if (nAtoms < 0)
    {
      // If the number of atoms is negative it means there is data between
      // the atom positions and grid data.
      negAtoms = true;
      nAtoms *= -1;
    }
//    cerr << "Number of atoms: " << nAtoms << " Orig: " << vs.at(0) << endl;
    double x = strtod(static_cast<const char*>(vs.at(1).c_str()), &endptr);
    if (endptr == static_cast<const char*>(vs.at(1).c_str()))
    {
      errorMsg << "Problems reading the Gaussian cube file: "
               << "Could not read line #3.\n"
               << "According to the specification this line should contain "
               << "four entries separated by white space.\n"
               << "OpenBabel could not interpret item #1 as a double.";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      return false;
    }
    double y = strtod(static_cast<const char*>(vs.at(2).c_str()), &endptr);
    if (endptr == static_cast<const char*>(vs.at(2).c_str()))
    {
      errorMsg << "Problems reading the Gaussian cube file: "
               << "Could not read line #3.\n"
               << "According to the specification this line should contain "
               << "four entries separated by white space.\n"
               << "OpenBabel could not interpret item #2 as a double.";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      return false;
    }
    double z = strtod(static_cast<const char*>(vs.at(3).c_str()), &endptr);
    if (endptr == static_cast<const char*>(vs.at(3).c_str()))
    {
      errorMsg << "Problems reading the Gaussian cube file: "
               << "Could not read line #3.\n"
               << "According to the specification this line should contain "
               << "four entries separated by white space.\n"
               << "OpenBabel could not interpret item #3 as a double.";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      return false;
    }
    vector3 origin(x, y, z);
    origin *= BOHR_TO_ANGSTROM;
//    cerr << "Origin: " << origin.x() << ", " << origin.y() << ", " << origin.z() << endl;

    // The next three lines contain the number of voxels in x, y and z along
    // with the three cube cell axes.
    vector<int> voxels(3);
    vector<vector3> axes(3);
    bool angstroms = true;
    for (int i = 0; i < 3; i++)
    {
      ++line;
      if (!ifs.getline(buffer, BUFF_SIZE))
      {
        errorMsg << "Problem reading the Gaussian cube file: cannot read the "
                 << line << "line.\nAccording to the specification this line "
                 << "should contain header information on the size of the cube.\n";
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
        return false;
      }
      vs.clear();
      tokenize(vs, buffer);
      if (vs.size() < 4)
      {
        errorMsg << "Problem reading the Gaussian cube file: cannot read the "
                 << i+4 << "line.\nAccording to the specification this line "
                 << "should contain four elements (cube header).\n";
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
        return false;
      }
      voxels[i] = strtol(static_cast<const char*>(vs.at(0).c_str()), &endptr, 10);
      if (endptr == static_cast<const char*>(vs.at(0).c_str()))
      {
        errorMsg << "Problems reading the Gaussian cube file: "
                 << "Could not read line " << i+4 << ".\n"
                 << "According to the specification this line should contain "
                 << "four entries separated by white space.\n"
                 << "OpenBabel could not interpret item #0 as an integer.";
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
        return false;
      }
      if (voxels[i] < 0)
      {
        // If the number of voxels is negative then the cube is in Bohr
        angstroms = false;
        voxels[i] *= -1;
      }
      x = strtod(static_cast<const char*>(vs.at(1).c_str()), &endptr);
      if (endptr == static_cast<const char*>(vs.at(1).c_str()))
      {
        errorMsg << "Problems reading the Gaussian cube file: "
                 << "Could not read line " << i+4 << ".\n"
                 << "According to the specification this line should contain "
                 << "four entries separated by white space.\n"
                 << "OpenBabel could not interpret item #1 as a double.";
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
        return false;
      }
      y = strtod(static_cast<const char*>(vs.at(2).c_str()), &endptr);
      if (endptr == static_cast<const char*>(vs.at(2).c_str()))
      {
        errorMsg << "Problems reading the Gaussian cube file: "
                 << "Could not read line " << i+4 << ".\n"
                 << "According to the specification this line should contain "
                 << "four entries separated by white space.\n"
                 << "OpenBabel could not interpret item #2 as a double.";
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
        return false;
      }
      z = strtod(static_cast<const char*>(vs.at(3).c_str()), &endptr);
      if (endptr == static_cast<const char*>(vs.at(3).c_str()))
      {
        errorMsg << "Problems reading the Gaussian cube file: "
                 << "Could not read line " << i+4 << ".\n"
                 << "According to the specification this line should contain "
                 << "four entries separated by white space.\n"
                 << "OpenBabel could not interpret item #3 as a double.";
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
        return false;
      }
      axes[i] = vector3(x, y, z);
      axes[i] *= BOHR_TO_ANGSTROM;
    }

    pmol->BeginModify();
    pmol->SetDimension(3);
    pmol->ReserveAtoms(nAtoms);

    // Now to read in the atoms contained in the Gaussian cube
    for (int i = 0; i < nAtoms; i++)
    {
      ++line;
      if (!ifs.getline(buffer, BUFF_SIZE))
      {
        errorMsg << "Problem reading the Gaussian cube file: cannot read the "
                 << line << "line.\nAccording to the specification this line "
                 << "should contain atomic number and coordinates.\n";
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
        return false;
      }
      vs.clear();
      tokenize(vs, buffer);
      if (vs.size() < 5)
      {
        errorMsg << "Problem reading the Gaussian cube file: cannot read the "
                 << line << "line.\nAccording to the specification this line "
                 << "should contain five elements (atom information).\n";
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
        return false;
      }

      // Only contains atoms - no bonds
      OBAtom *atom = pmol->NewAtom();

      int atomicNum = strtol(static_cast<const char*>(vs.at(0).c_str()), &endptr, 10);
      if (endptr == static_cast<const char*>(vs.at(0).c_str()))
      {
        errorMsg << "Problems reading the Gaussian cube file: "
                 << "Could not read line " << line << ".\n"
                 << "According to the specification this line should contain "
                 << "five entries separated by white space.\n"
                 << "OpenBabel could not interpret item #0 as an integer.";
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
        return false;
      }

      atom->SetAtomicNum(atomicNum);

      // Read the atom coordinates
      x = strtod(static_cast<const char*>(vs.at(2).c_str()), &endptr);
      if (endptr == static_cast<const char*>(vs.at(2).c_str()))
      {
        errorMsg << "Problems reading the Gaussian cube file: "
                 << "Could not read line " << line << ".\n"
                 << "According to the specification this line should contain "
                 << "five entries separated by white space.\n"
                 << "OpenBabel could not interpret item #2 as a double.";
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
        return false;
      }
      y = strtod(static_cast<const char*>(vs.at(3).c_str()), &endptr);
      if (endptr == static_cast<const char*>(vs.at(3).c_str()))
      {
        errorMsg << "Problems reading the Gaussian cube file: "
                 << "Could not read line " << line << ".\n"
                 << "According to the specification this line should contain "
                 << "five entries separated by white space.\n"
                 << "OpenBabel could not interpret item #3 as a double.";
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
        return false;
      }
      z = strtod(static_cast<const char*>(vs.at(4).c_str()), &endptr);
      if (endptr == static_cast<const char*>(vs.at(4).c_str()))
      {
        errorMsg << "Problems reading the Gaussian cube file: "
                 << "Could not read line " << line << ".\n"
                 << "According to the specification this line should contain "
                 << "five entries separated by white space.\n"
                 << "OpenBabel could not interpret item #4 as a double.";
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
        return false;
      }
      // Convert units...
      x *= BOHR_TO_ANGSTROM;
      y *= BOHR_TO_ANGSTROM;
      z *= BOHR_TO_ANGSTROM;
      atom->SetVector(x, y, z);
    }

    // Now to read in the actual cube data
    OBGridData *gd = new OBGridData;
    gd->SetNumberOfPoints(voxels[0], voxels[1], voxels[2]);
    gd->SetLimits(origin, axes[0], axes[1], axes[2]);
    gd->SetUnit(angstroms ? OBGridData::ANGSTROM : OBGridData::BOHR);

    // If the number of atoms was negative then there is some data between the
    // atom data and the cube data.
    if (negAtoms)
    {
      ++line;
      if (!ifs.getline(buffer, BUFF_SIZE))
      {
        obErrorLog.ThrowError(__FUNCTION__,
                              "Problem reading the Gaussian cube file: cannot read the line between the atom data and cube data.", obError);
        return false;
      }
    }
    vector<double> values;
    int n = voxels[0]*voxels[1]*voxels[2];
    values.reserve(n);
    while (values.size() < n)
    {
      // Read in values until we have a complete row of data
      ++line;
      if (!ifs.getline(buffer, BUFF_SIZE))
      {
        errorMsg << "Problem reading the Gaussian cube file: cannot"
                 << " read line " << line
                 << " of the file. More data was expected.\n"
                 << "Values read in = " << values.size()
                 << " and expected number of values = "
                 << voxels[0]*voxels[1]*voxels[2] << endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
        return false;
      }
      tokenize(vs, buffer);
      if (vs.size() == 0)
      {
        errorMsg << "Problem reading the Gaussian cube file: cannot"
                 << " read line " << line
                 << ", there does not appear to be any data in it.\n"
                 << buffer << "\n";
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
        return false;
      }

      for (int i = 0; i < vs.size(); ++i)
      {
        values.push_back(strtod(static_cast<const char*>(vs.at(i).c_str()), &endptr));
      }
    }

    // Now to translate from a vector stored in x, y, z to one stored in z, y, x
    // Gaussian cube format order to OpenBabel GridData order!
    vector<double> nValues;
    nValues.resize(values.size());
    for (int k = 0; k < voxels[2]; ++k)
      for (int j = 0; j < voxels[1]; ++j)
        for (int i = 0; i < voxels[0]; ++i)
        {
          nValues[k*voxels[0]*voxels[1] + j*voxels[0] + i] = values[i*voxels[1]*voxels[2] + j*voxels[2] + k];
        }

    gd->SetValues(nValues);
    gd->SetOrigin(fileformatInput); // i.e., is this data from a file or determined by Open Babel

    pmol->EndModify();
    
    // clean out any remaining blank lines
    while(ifs.peek() != EOF && ifs.good() && 
          (ifs.peek() == '\n' || ifs.peek() == '\r'))
      ifs.getline(buffer,BUFF_SIZE);

    // Connect the dots and guess bond orders
    if (!pConv->IsOption("b", OBConversion::INOPTIONS))
      pmol->ConnectTheDots();
    if (!pConv->IsOption("s", OBConversion::INOPTIONS) && !pConv->IsOption("b", OBConversion::INOPTIONS))
      pmol->PerceiveBondOrders();

    pmol->SetData(gd);
    pmol->SetDimension(3); // always a 3D structure

    return true;
  }
}
