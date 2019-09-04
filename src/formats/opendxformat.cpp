/**********************************************************************
opendxformat.cpp - Read in OpenDX cube format files from APBS

 Copyright (c) 2007-2008 by Geoffrey R. Hutchison
 Some Portions Copyright (C) 2008 by Marcus D. Hanwell

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

#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/elements.h>

#include <openbabel/griddata.h>
#include <cstdlib>

using namespace std;

namespace OpenBabel
{

  // Ultimately, this should just read in an OBFloatGrid, rather than a "molecule"
  // since there's no actual molecular data. For now, we'll read in an OBMol with no atoms
  // and attach an OBGridData
  class OBOpenDXCubeFormat : public OpenBabel::OBMoleculeFormat
  {
  public:
    // Constructor: register "dx" extension code
    OBOpenDXCubeFormat()
    {
        OpenBabel::OBConversion::RegisterFormat( "dx", this );
    }

    // Return description.
    virtual const char* Description() //required
    {
        return
        "OpenDX cube format for APBS\n"
        "A volume data format for IBM's Open Source visualization software\n"
          "The OpenDX support is currently designed to read the OpenDX cube\n"
          "files from APBS. \n\n";
    }

    // Return a specification url, not really a specification since
    // I couldn't find it but close enough.
    virtual const char* SpecificationURL()
    {
        return "http://apbs.sourceforge.net/doc/user-guide/index.html#id504516";
    }

    // Return MIME type, NULL in this case.
    virtual const char* GetMIMEType() { return 0; };

    // Skip to object: used for multi-object file formats.
    virtual int SkipObjects( int n, OpenBabel::OBConversion* pConv ) { return 0; }

    virtual unsigned int Flags()
    {
      return READONEONLY | WRITEONEONLY | ZEROATOMSOK;
    };

    /// The "API" interface functions
    virtual bool ReadMolecule( OpenBabel::OBBase* pOb, OpenBabel::OBConversion* pConv );
    virtual bool WriteMolecule( OpenBabel::OBBase* pOb, OpenBabel::OBConversion* pConv );

};

//------------------------------------------------------------------------------

    // Global variable used to register OpenDX cube format.
    OBOpenDXCubeFormat theOpenDXCubeFormat;

//------------------------------------------------------------------------------
bool OBOpenDXCubeFormat::ReadMolecule( OBBase* pOb, OBConversion* pConv )
{
    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol == 0)
      return false;

    istream& ifs = *pConv->GetInStream();

    const char* title = pConv->GetTitle();
    char buffer[BUFF_SIZE];

    stringstream errorMsg;

    if (!ifs)
      return false; // We are attempting to read past the end of the file

    pmol->SetTitle(title);

    while (ifs.good() && ifs.getline(buffer,BUFF_SIZE)) {
      if (buffer[0] == '#')
        continue; // comment line
      if (EQn(buffer, "object", 6))
        break;
    }
    if (!ifs)
      return false; // ran out of lines

    vector<string> vs;
    tokenize(vs, buffer);

    // Number of grid points (voxels)
    vector<int> voxels(3);
    if (!EQn(buffer, "object", 6) || vs.size() != 8)
      return false;
    else {
      voxels[0] = atoi(vs[5].c_str());
      voxels[1] = atoi(vs[6].c_str());
      voxels[2] = atoi(vs[7].c_str());
    }

    double x, y, z;
    if (!ifs.getline(buffer, BUFF_SIZE) || !EQn(buffer, "origin", 6))
      return false;
    else {
      tokenize(vs, buffer);
      if (vs.size() != 4)
        return false;
      x = atof(vs[1].c_str());
      y = atof(vs[2].c_str());
      z = atof(vs[3].c_str());
    }
    vector3 origin(x, y, z);

    // now three lines with the x, y, and z axes
    vector<vector3> axes;
    for (unsigned int i = 0; i < 3; ++i) {
      if (!ifs.getline(buffer, BUFF_SIZE) || !EQn(buffer, "delta", 5))
        return false;
      else {
        tokenize(vs, buffer);
        if (vs.size() != 4)
          return false;
        x = atof(vs[1].c_str());
        y = atof(vs[2].c_str());
        z = atof(vs[3].c_str());
        axes.push_back(vector3(x, y, z));
      }
    }

    // Two remaining header lines before the data:
    /*
      object 2 class gridconnections counts nx ny nz
      object 3 class array type double rank 0 times n data follows
    */
    if (!ifs.getline(buffer, BUFF_SIZE) || !EQn(buffer, "object", 6))
      return false;
    if (!ifs.getline(buffer, BUFF_SIZE) || !EQn(buffer, "object", 6))
      return false;

    pmol->BeginModify();
    pmol->SetDimension(3);

    OBGridData *gd = new OBGridData;
    gd->SetAttribute("OpenDX");

    // get all values as one vector<double>
    char *endptr;
    vector<double> values;
    int n = voxels[0]*voxels[1]*voxels[2];
    int line = 0;
    values.reserve(n);
    while (ifs.getline(buffer, BUFF_SIZE))
    {
      ++line;
      if (EQn(buffer, "attribute", 9))
        break; // we're finished with reading data -- although we should probably have a voxel check in here too

      tokenize(vs, buffer);
      if (vs.size() == 0)
      {
        errorMsg << "Problem reading the OpenDX grid file: cannot"
                 << " read line " << line
                 << ", there does not appear to be any data in it.\n"
                 << buffer << "\n";
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
        return false;
      }

      for (unsigned int l = 0; l < vs.size(); ++l)
      {
        values.push_back(strtod(static_cast<const char*>(vs[l].c_str()), &endptr));
      }
    }

    gd->SetNumberOfPoints(voxels[0], voxels[1], voxels[2]);
    gd->SetLimits(origin, axes[0], axes[1], axes[2]);
    gd->SetUnit(OBGridData::ANGSTROM);
    gd->SetOrigin(fileformatInput); // i.e., is this data from a file or determined by Open Babel
    gd->SetValues(values); // set the values
    pmol->SetData(gd); // store the grids in the OBMol
    pmol->EndModify();

    // Trailing control lines
    /*
     attribute "dep" string "positions"
     object "regular positions regular connections" class field
     component "positions" value 1
     component "connections" value 2
     component "data" value 3
    */
    if (!ifs.getline(buffer, BUFF_SIZE) || !EQn(buffer, "object", 6))
      return false;
    if (!ifs.getline(buffer, BUFF_SIZE) || !EQn(buffer, "component", 9))
      return false;
    if (!ifs.getline(buffer, BUFF_SIZE) || !EQn(buffer, "component", 9))
      return false;
    if (!ifs.getline(buffer, BUFF_SIZE) || !EQn(buffer, "component", 9))
      return false;

    // clean out any remaining blank lines
    std::streampos ipos;
    do
    {
      ipos = ifs.tellg();
      ifs.getline(buffer,BUFF_SIZE);
    }
    while(strlen(buffer) == 0 && !ifs.eof() );
    ifs.seekg(ipos);

    return true;
  }

//------------------------------------------------------------------------------
  bool OBOpenDXCubeFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    char buffer[BUFF_SIZE];
    string str;
    stringstream errorMsg;

    OBGridData *gd = (OBGridData*)mol.GetData(OBGenericDataType::GridData);
    if (gd == NULL) {
      errorMsg << "The molecule has no grid.";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      return false;
    }

    // APBS-style OpenDX Multigrid
    // First some comments
    ofs << "# Data from Open Babel " << BABEL_VERSION << "\n";
    str = mol.GetTitle();
    if (str.empty())
      ofs << "# Molecule Title: *****" << "\n";
    else
      ofs << "# Molecule Title: " << str << "\n";

    int nx, ny, nz;
    double origin[3], xAxis[3], yAxis[3], zAxis[3];
    gd->GetAxes(xAxis, yAxis, zAxis);
    gd->GetNumberOfPoints(nx, ny, nz);
    gd->GetOriginVector(origin);

    // data line 1: # of points in x, y, z (nx, ny, nz)
    snprintf(buffer, BUFF_SIZE, "object 1 class gridpositions counts %5d %5d %5d", nx, ny, nz);
    ofs << buffer << "\n";

    // data line 2: origin (x, y, z)
    snprintf(buffer, BUFF_SIZE,"origin %12.6f %12.6f %12.6f",
        origin[0], origin[1], origin[2]);
    ofs << buffer << "\n";

    // data line 3: x-displacement
    snprintf(buffer, BUFF_SIZE,"delta %12.6f %12.6f %12.6f",
        xAxis[0], xAxis[1], xAxis[2]);
    ofs << buffer << "\n";

    // data line 4: y-displacement
    snprintf(buffer, BUFF_SIZE,"delta %12.6f %12.6f %12.6f",
        yAxis[0], yAxis[1], yAxis[2]);
    ofs << buffer << "\n";

    // data line 5: z-displacement
    snprintf(buffer, BUFF_SIZE,"delta %12.6f %12.6f %12.6f",
        zAxis[0], zAxis[1], zAxis[2]);
    ofs << buffer << "\n";

    // data line 6: # of points in x, y, z (nx, ny, nz)
    snprintf(buffer, BUFF_SIZE, "object 2 class gridconnections counts %5d %5d %5d", nx, ny, nz);
    ofs << buffer << "\n";

    // data line 7: total # of points
    snprintf(buffer, BUFF_SIZE, "object 3 class array type double rank 0 items %5d data follows", nx*ny*nz);
    ofs << buffer << "\n";

    // The cube(s)
    double value;
    unsigned int count = 1;
    for (int i = 0; i < nx; ++i)
    {
      for (int j = 0; j < ny; ++j)
      {
        for (int k = 0; k < nz; ++k)
        {
          value = gd->GetValue(i, j, k);
          snprintf(buffer, BUFF_SIZE," %12.5E", value);
          if (count % 3 == 0)
            ofs << buffer << "\n";
          else
            ofs << buffer;
          count++;
        } // z-axis
      } // y-axis
    } // x-axis

    if (count % 3 != 0)
      ofs << "\n";
    ofs << "attribute \"dep\" string \"positions\"\n";
    ofs << "object \"regular positions regular connections\" class field\n";
    ofs << "component \"positions\" value 1\n";
    ofs << "component \"connections\" value 2\n";
    ofs << "component \"data\" value 3\n";

    return true;
  }
}
