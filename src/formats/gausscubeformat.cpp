/**********************************************************************
gausscubeformat.cpp - Read in Gaussian cube format files.

// Molekel - Molecular Visualization Program
// Copyright (C) 2006, 2007 Swiss National Supercomputing Centre (CSCS)

 Some Portions Copyright (c) 2007 by Geoffrey R. Hutchison
 Some Portions Copyright (C) 2008 by Marcus D. Hanwell
 Some Portions Copyright (C) 2008 by Tim Vandermeersch

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

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
#include <cstdlib>

#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/elements.h>

#include <openbabel/griddata.h>

using namespace std;

/// @warning seems that every position in cube files has always to be
/// multiplied by this constant even if according to:
/// reference: https://www.gaussian.com/cubegen/
/// <Num points along first axis> > 0 means that the values
/// are already in Angstroms.
static const double BOHR_TO_ANGSTROM = 0.529177249;
static const double ANGSTROM_TO_BOHR = 1.0 / 0.529177249;

namespace OpenBabel
{

  class OBGaussianCubeFormat : public OpenBabel::OBMoleculeFormat
  {
  public:
    // Constructor: register 'cube' and "CUBE" format.
    OBGaussianCubeFormat()
    {
        OpenBabel::OBConversion::RegisterFormat( "cube", this );
        OpenBabel::OBConversion::RegisterFormat( "cub", this );
    }

    // Return description.
    virtual const char* Description() //required
    {
        return
        "Gaussian cube format\n"
        "A grid format for volume data used by Gaussian\n"
        "Open Babel supports reading and writing Gaussian cubes, including multiple\n"
        "grids in one file.\n\n"

        "Read Options e.g. -as\n"
        "  b no bonds\n"
        "  s no multiple bonds\n\n";
    }

    // Return a specification url, not really a specification since
    // I couldn't find it but close enough.
    virtual const char* SpecificationURL()
    {
        return "https://www.gaussian.com/cubegen/";
    }

    // Return MIME type, NULL in this case.
    virtual const char* GetMIMEType() { return 0; };

    // Skip to object: used for multi-object file formats.
    virtual int SkipObjects( int n, OpenBabel::OBConversion* pConv ) { return 0; }

    virtual unsigned int Flags()
    {
        return 0;
    };

    /// The "API" interface functions
    virtual bool ReadMolecule( OpenBabel::OBBase* pOb, OpenBabel::OBConversion* pConv );
    /// Write: always returns false right now - read only
    virtual bool WriteMolecule( OpenBabel::OBBase* pOb, OpenBabel::OBConversion* pConv );

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
    string cubeTitle = buffer; // save for later use

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

    int nCubes = 1;
    vector<OBGridData*> vgd;
    // If the number of atoms was negative then there is some data between the
    // atom data and the cube data, i.e. we have multiple cubes...
    if (negAtoms)
    {
      ++line;
      if (!ifs.getline(buffer, BUFF_SIZE))
      {
        obErrorLog.ThrowError(__FUNCTION__,
                              "Problem reading the Gaussian cube file: cannot read the line between the atom data and cube data.", obError);
        return false;
      }

      tokenize(vs, buffer);
      if (vs.size() < 1) return false; // timvdm 18/06/2008
      nCubes = strtol(static_cast<const char*>(vs.at(0).c_str()), &endptr, 10);
      if (endptr == static_cast<const char*>(vs.at(0).c_str()))
      {
        errorMsg << "Problems reading the Gaussian cube file: "
                 << "Could not read line " << line << ".\n"
                 << "According to the specification this line should contain "
                 << "the number of cubes and a number for each cube (orbital number).\n"
                 << "OpenBabel could not interpret item #0 as an integer.";
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
        return false;
      }

      vgd.reserve(nCubes);
      for (unsigned int i = 1; i < vs.size(); ++i)
      {
        int cubeNumber = strtol(static_cast<const char*>(vs.at(i).c_str()), &endptr, 10);
        if (endptr == static_cast<const char*>(vs.at(i).c_str()))
        {
          errorMsg << "Problems reading the Gaussian cube file: "
                   << "Could not read line " << line << ".\n"
                   << "According to the specification this line should contain "
                   << "the number of cubes and a number for each cube (orbital number).\n"
                   << "OpenBabel could not interpret item #" << vgd.size()
                   << " as an integer.";
          obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
          return false;
        }
        vgd.push_back(new OBGridData);
        stringstream cubeTitle;
        cubeTitle << "MO " << cubeNumber;
        vgd.at(vgd.size()-1)->SetAttribute(cubeTitle.str().c_str());
      }
      // Now if we have lots of cubes we need to read in more line(s)
      while (vgd.size() < nCubes)
      {
        ++line;
        if (!ifs.getline(buffer, BUFF_SIZE))
        {
          errorMsg << "Problem reading the Gaussian cube file: cannot"
                   << " read line " << line
                   << " of the file. More data was expected.\n"
                   << "Grid MO titles read in = " << vgd.size()
                   << " and expected number of values = "
                   << nCubes << endl;
          obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
          return false;
        }
        tokenize(vs, buffer);
        for (unsigned int i = 0; i < vs.size(); ++i)
        {
          int cubeNumber = strtol(static_cast<const char*>(vs.at(i).c_str()), &endptr, 10);
          if (endptr == static_cast<const char*>(vs.at(i).c_str()))
          {
            errorMsg << "Problems reading the Gaussian cube file: "
                     << "Could not read line " << line << ".\n"
                     << "According to the specification this line should contain "
                     << "the number of cubes and a number for each cube (orbital number).\n"
                     << "OpenBabel could not interpret item #" << vgd.size()
                     << " as an integer.";
            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
            return false;
          }
          vgd.push_back(new OBGridData);
          stringstream cubeTitle;
          cubeTitle << "MO " << cubeNumber;
          vgd.at(vgd.size()-1)->SetAttribute(cubeTitle.str().c_str());
        }
      }
    }
    else
    {
      // only one cube
      vgd.push_back(new OBGridData);
      vgd.at(0)->SetAttribute(cubeTitle.c_str());
    }

    // get all values as one vector<double>
    vector<double> values;
    int n = voxels[0]*voxels[1]*voxels[2];
    values.reserve(n*nCubes);
    while (values.size() < n*nCubes)
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
                 << n*nCubes << endl;
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

      for (unsigned int l = 0; l < vs.size(); ++l)
      {
        values.push_back(strtod(static_cast<const char*>(vs.at(l).c_str()), &endptr));
      }
    }

    // translate the vector<double> to vector<vector<double> >
    vector<vector<double> > vvd;
    vvd.resize(nCubes);
    for (unsigned int i = 0; i < values.size(); )
    {
      for (int j = 0; j < nCubes; ++j, ++i) // foreach cube
      {
        vvd[j].push_back(values[i]);
      }
    }

    for (int i = 0; i < nCubes; ++i) // foreach cube
    {
      vgd[i]->SetNumberOfPoints(voxels[0], voxels[1], voxels[2]);
      vgd[i]->SetLimits(origin, axes[0], axes[1], axes[2]);
      vgd[i]->SetUnit(angstroms ? OBGridData::ANGSTROM : OBGridData::BOHR);
      vgd[i]->SetOrigin(fileformatInput); // i.e., is this data from a file or determined by Open Babel
      vgd[i]->SetValues(vvd[i]); // set the values
      pmol->SetData(vgd[i]); // store the grids in the OBMol
    }

    pmol->EndModify();

    // clean out any remaining blank lines
    std::streampos ipos;
    do
    {
      ipos = ifs.tellg();
      ifs.getline(buffer,BUFF_SIZE);
    }
    while(strlen(buffer) == 0 && !ifs.eof() );
    ifs.seekg(ipos);

    // Connect the dots and guess bond orders
    if (!pConv->IsOption("b", OBConversion::INOPTIONS))
      pmol->ConnectTheDots();
    if (!pConv->IsOption("s", OBConversion::INOPTIONS) && !pConv->IsOption("b", OBConversion::INOPTIONS))
      pmol->PerceiveBondOrders();

    pmol->SetDimension(3); // always a 3D structure

    return true;
  }

//------------------------------------------------------------------------------
  bool OBGaussianCubeFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    char buffer[BUFF_SIZE];
    string str;
    stringstream errorMsg;

    // first two lines are comments
    str = mol.GetTitle();
    if (str.empty())
      ofs << "*****" << endl;
    else
      ofs << str << endl;

    ofs << endl; // line 2

    OBGridData *gd = (OBGridData*)mol.GetData(OBGenericDataType::GridData);
    if (gd == NULL) {
      errorMsg << "The molecule has no grid.";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      return false;
    }

    int nx, ny, nz;
    double origin[3], xAxis[3], yAxis[3], zAxis[3];
    gd->GetAxes(xAxis, yAxis, zAxis);
    gd->GetNumberOfPoints(nx, ny, nz);
    gd->GetOriginVector(origin);

    // line 3: number of atoms, origin x y z
    snprintf(buffer, BUFF_SIZE,"%5d%12.6f%12.6f%12.6f", - static_cast<signed int> (mol.NumAtoms()),
        origin[0]*ANGSTROM_TO_BOHR, origin[1]*ANGSTROM_TO_BOHR, origin[2]*ANGSTROM_TO_BOHR);
    ofs << buffer << endl;

    // line 4: number of points x direction, axis x direction x y z
    snprintf(buffer, BUFF_SIZE,"%5d%12.6f%12.6f%12.6f", nx,
        xAxis[0]*ANGSTROM_TO_BOHR, xAxis[1]*ANGSTROM_TO_BOHR, xAxis[2]*ANGSTROM_TO_BOHR);
    ofs << buffer << endl;

    // line 5: number of points y direction, axis y direction x y z
    snprintf(buffer, BUFF_SIZE,"%5d%12.6f%12.6f%12.6f", ny,
        yAxis[0]*ANGSTROM_TO_BOHR, yAxis[1]*ANGSTROM_TO_BOHR, yAxis[2]*ANGSTROM_TO_BOHR);
    ofs << buffer << endl;

    // line 6: number of points z direction, axis z direction x y z
    snprintf(buffer, BUFF_SIZE,"%5d%12.6f%12.6f%12.6f", nz,
        zAxis[0]*ANGSTROM_TO_BOHR, zAxis[1]*ANGSTROM_TO_BOHR, zAxis[2]*ANGSTROM_TO_BOHR);
    ofs << buffer << endl;

    // Atom lines: atomic number, ?, X, Y, Z
    FOR_ATOMS_OF_MOL (atom, mol) {
      double *coordPtr = atom->GetCoordinate();
      snprintf(buffer, BUFF_SIZE,"%5d%12.6f%12.6f%12.6f%12.6f", atom->GetAtomicNum(),
          static_cast<double>(atom->GetAtomicNum()),
          coordPtr[0]*ANGSTROM_TO_BOHR, coordPtr[1]*ANGSTROM_TO_BOHR, coordPtr[2]*ANGSTROM_TO_BOHR);
      ofs << buffer << endl;
    }

    vector<OBGenericData*> grids = pmol->GetAllData(OBGenericDataType::GridData);
    snprintf(buffer, BUFF_SIZE," %5lu", (unsigned long)grids.size());
    ofs << buffer << flush;
    for (unsigned int l = 1; l <= grids.size(); ++l)
    {
      snprintf(buffer, BUFF_SIZE," %3d", l);
      ofs << buffer << flush;
    }
    ofs << endl;

    for (unsigned int l = 0; l < grids.size(); ++l)
    {
      gd = static_cast<OBGridData*>(grids[l]);
      int mx, my, mz;
      gd->GetNumberOfPoints(mx, my, mz);

      if ((nx != mx) || (ny != my) || (nz != mz))
      {
          errorMsg << "Problem writing the Gaussian cube file: cube " << l
                   << " does not have the same dimentions as cube 0.\n"
                   << "This cube will be skipped.\n";
          obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      }
    }

    // The cube(s)
    double value;
    unsigned int count = 1;
    for (int i = 0; i < nx; ++i)
    {
      for (int j = 0; j < ny; ++j)
      {
        for (int k = 0; k < nz; ++k)
        {
          // The cube files are stored in staggered z
          for (unsigned int l = 0; l < grids.size(); ++l)
          {
            value = static_cast<OBGridData*>(grids[l])->GetValue(i, j, k);
            snprintf(buffer, BUFF_SIZE," %12.5E", value);
            if (count % 6 == 0)
              ofs << buffer << endl;
            else
              ofs << buffer;
            count++;
          }
        } // z-axis
      } // y-axis
    } // x-axis

    return true;
  }
}
