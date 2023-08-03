/**********************************************************************
obxyz - Open Babel XYZ transformations

Copyright (C) 2001-2006 Geoffrey R. Hutchison

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/


// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif

#include <openbabel/babelconfig.h>

#include <openbabel/mol.h>
#include <openbabel/math/matrix3x3.h>
#include <openbabel/obconversion.h>

using namespace std;
using namespace OpenBabel;

int main(int argc,char **argv)
{
  char c;
  int transX, transY, transZ, centerIdx;
  transX = transY = transZ = centerIdx = -1;

  char *FileIn = nullptr;
  char *program_name = argv[0];
  //   char *iext;

  OBConversion conv(&cin,&cout);
  OBFormat *pFormat = conv.FindFormat("smi"); // default format is SMILES

  // Still need to add:
  // rotate X, Y, or Z by set angle
  // translate X, Y, or Z by set amount

  // Parse options
  while ((c = getopt(argc, argv, "x:y:z:c:")) != -1)
    {
      switch (c)
        {
        case 'c': /// atom to be centered

          c = sscanf(optarg, "%d", &centerIdx);
          if (c != 1 )
            {
              cerr << program_name << ": unable to parse -c option" << endl;
              exit (-1);
            }
          break;

        case 'x': /// atom to be centered

          c = sscanf(optarg, "%d", &transX);
          if (c != 1 )
            {
              cerr << program_name << ": unable to parse -x option" << endl;
              exit (-1);
            }
          break;

        case 'y': /// atom to y-axis

          c = sscanf(optarg, "%d", &transY);
          if (c != 1 )
            {
              cerr << program_name << ": unable to parse -y option" << endl;
              exit (-1);
            }
          break;

        case 'z': /// atom to z-axis

          c = sscanf(optarg, "%d", &transZ);
          if (c != 1 )
            {
              cerr << program_name << ": unable to parse -z option" << endl;
              exit (-1);
            }
          break;

        }
    }

  ifstream ifs;
  FileIn  = argv[optind];
  if (FileIn != nullptr)
    {
      // Read the file
      ifs.open(FileIn);
      if (!ifs)
        {
          cerr << program_name << ": cannot read input file!" << endl;
          exit (-1);
        }
      conv.SetInStream(&ifs);

      // Find Input filetype
      pFormat = conv.FormatFromExt(FileIn);
      if (pFormat == nullptr)
        {
          cerr << program_name << ": cannot read input format!" << endl;
          return (-1);
        }
    }

  if (! conv.SetInAndOutFormats(pFormat, pFormat))
    {
      cerr << program_name << ": cannot read or write to this file format" << endl;
      return (-1);
    }

  OBMol mol;
  mol.Clear();
  conv.Read(&mol);
  if (mol.Empty())
    return(1);

  vector3 v;
  if (centerIdx != -1)
    {
      v = (mol.GetAtom(centerIdx))->GetVector();
      mol.Translate(-1.0f*v);
    }

  double alpha, beta, gamma;
  alpha = beta = gamma = 0.0f;
  if (transX != -1)
    {
      v = (mol.GetAtom(transX))->GetVector();

      // angle to rotate vector into XZ plane
      // (i.e., angle between x and y components)
      gamma = M_PI/2.0 - atan(v.x() / v.y());

      // angle to rotate vector from XZ plane to X axis
      // (by rotation along Y-axis)
      beta = atan(v.z() / sqrt(v.y()*v.y() + v.x()*v.x()) );
    }
  else if (transY != -1)
    {
      v = (mol.GetAtom(transY))->GetVector();

      // angle to rotate vector into YZ plane
      // (i.e., angle between x and y components)
      gamma = -1.0 * atan(v.x() / v.y());

      // angle to rotate vector from YZ plane to Y axis
      // (by rotation along X-axis)
      alpha = -1.0 * atan(v.z() / sqrt(v.y()*v.y() + v.x()*v.x()) );
    }
  else if (transZ != -1)
    {
      v = (mol.GetAtom(transZ))->GetVector();

      // angle to rotate vector into YZ plane
      // (i.e., angle between x and y components)
      gamma = -1.0 * atan(v.x() / v.y());

      // angle to rotate vector from YZ plane to Z axis
      // (by rotation along X axis)
      alpha = M_PI/2.0f - atan(v.z() / sqrt(v.y()*v.y() + v.x()*v.x()) );
    }

  if ( !IsNearZero(alpha) ||
       !IsNearZero(beta)  ||
       !IsNearZero(gamma) )
    {
      matrix3x3 mat;
      mat.SetupRotMat(alpha*RAD_TO_DEG, beta*RAD_TO_DEG, gamma*RAD_TO_DEG);
      double array[9];
      mat.GetArray(array);
      mol.Rotate(array);
    }

  conv.Write(&mol, &cout);

  return(0);
}
