/**********************************************************************
Copyright (C) 2003-2005 by Geoffrey R. Hutchison
 
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
#include <openbabel/obutil.h>
#include <openbabel/generic.h>
#include <cstdlib>

using namespace std;
using namespace OpenBabel;

#ifdef TESTDATADIR
    string ptestdatadir = TESTDATADIR;
    string punitcell_file = ptestdatadir + "unitcell.txt";
    string presults_file = ptestdatadir + "unitcell_results.txt";
#else
    string punitcell_file = "files/unitcell.txt";
    string presults_file = "files/unitcell_results.txt";
#endif

int unitcell(int argc, char* argv[])
{
  int defaultchoice = 1;
  
  int choice = defaultchoice;

  if (argc > 1) {
    if(sscanf(argv[1], "%d", &choice) != 1) {
      printf("Couldn't parse that input as a number\n");
      return -1;
    }
  }

    cout << "# Testing unit cell transformations ..." << endl;
    cout << "1..12" << endl;

    double a, b, c, alpha, beta, gamma;
    vector3 v1, v2, v3, v4, v5, v6, v7, v8, v9;
    vector3 coords1, coords2, tmpcoords;
    double x = 0.0, y = 0.0, z = 0.0;
    char buffer[BUFF_SIZE];
    std::ifstream ifs, results;
    OBUnitCell cell, cell2;
    vector<vector3> v3Return;
    vector<string> vs;
    unsigned int currTest = 1;

    if (!SafeOpen(ifs, punitcell_file.c_str()))
      {
	cout << "Bail out! Couldn't open test file 'unitcell.txt'" << endl;
        return(-1);
      }
    if (!SafeOpen(results, presults_file.c_str()))
      {
	cout << "Bail out! Couldn't open test file 'unitcell_results.txt'" << endl;
        return(-1);
      }

    // Get cell vectors
    ifs.getline(buffer,BUFF_SIZE);
    sscanf(buffer,"%lf %lf %lf",&x, &y, &z);
    v1.Set(x, y, z);

    ifs.getline(buffer,BUFF_SIZE);
    sscanf(buffer,"%lf %lf %lf",&x, &y, &z);
    v2.Set(x, y, z);

    ifs.getline(buffer,BUFF_SIZE);
    sscanf(buffer,"%lf %lf %lf",&x, &y, &z);
    v3.Set(x, y, z);

    cell.SetData(v1, v2, v3);

    // Get test coordinates
    ifs.getline(buffer,BUFF_SIZE); // blank

    // Cartesian, outside of uc
    ifs.getline(buffer,BUFF_SIZE);
    sscanf(buffer,"%lf %lf %lf",&x, &y, &z);
    coords1.Set(x, y, z);

    // Fractional, outside of uc
    ifs.getline(buffer,BUFF_SIZE);
    sscanf(buffer,"%lf %lf %lf",&x, &y, &z);
    coords2.Set(x, y, z);

    // Test
    a = cell.GetA();
    b = cell.GetB();
    c = cell.GetC();
    alpha = cell.GetAlpha();
    beta = cell.GetBeta();
    gamma = cell.GetGamma();

    results.getline(buffer,BUFF_SIZE);
    tokenize(vs,buffer);
    if (vs.size() != 6)
      {
	cout << "Bail out! Cannot parse results file 'unitcell_results.txt'" 
	     << endl;
	return(-1);
      }

    if ( IsNear(a, atof(vs[0].c_str()), 1.0e-3) )
      cout << "ok " << currTest++ << " # a distance" << endl;
    else
      cout << "not ok " << currTest++ << " # a distance" << endl;
    if ( IsNear(b, atof(vs[1].c_str()), 1.0e-3) )
      cout << "ok " << currTest++ << " # b distance" << endl;
    else
      cout << "not ok " << currTest++ << " # b distance" << endl;
    if ( IsNear(c, atof(vs[2].c_str()), 1.0e-3) )
      cout << "ok " << currTest++ << " # c distance" << endl;
    else
      cout << "not ok " << currTest++ << " # c distance" << endl;

    if ( IsNear(alpha, atof(vs[3].c_str()), 1.0e-3) )
      cout << "ok " << currTest++ << " # alpha angle" << endl;
    else
      cout << "not ok " << currTest++ << " # alpha angle" << alpha << endl;
    if ( IsNear(beta, atof(vs[4].c_str()), 1.0e-3) )
      cout << "ok " << currTest++ << " # beta angle" << endl;
    else
      cout << "not ok " << currTest++ << " # beta angle " << beta << endl;
    if ( IsNear(gamma, atof(vs[5].c_str()), 1.0e-3) )
      cout << "ok " << currTest++ << " # gamma angle" << endl;
    else
      cout << "not ok " << currTest++ << " # gamma angle: " << gamma << endl;


    // check to see if vector and a,b,c methods are equivalent
    cell2.SetData(a, b, c, alpha, beta, gamma);

    if ( IsNear(a, cell2.GetA(), 1.0e-3) )
      cout << "ok " << currTest++ << " # a distance" << endl;
    else
      cout << "not ok " << currTest++ << " # a distance" << endl;
    if ( IsNear(b, cell2.GetB(), 1.0e-3) )
      cout << "ok " << currTest++ << " # b distance" << endl;
    else
      cout << "not ok " << currTest++ << " # b distance" << endl;
    if ( IsNear(c, cell2.GetC(), 1.0e-3) )
      cout << "ok " << currTest++ << " # c distance" << endl;
    else
      cout << "not ok " << currTest++ << " # c distance" << endl;

    if ( IsNear(alpha, cell2.GetAlpha(), 1.0e-3) )
      cout << "ok " << currTest++ << " # a angle" << endl;
    else
      cout << "not ok " << currTest++ << " # a angle" << endl;
    if ( IsNear(beta, cell2.GetBeta(), 1.0e-3) )
      cout << "ok " << currTest++ << " # beta angle" << endl;
    else
      cout << "not ok " << currTest++ << " # beta angle" << endl;
    if ( IsNear(gamma, cell2.GetGamma(), 1.0e-3) )
      cout << "ok " << currTest++ << " # gamma angle" << endl;
    else
      cout << "not ok " << currTest++ << " # gamma angle" << endl;

    v3Return = cell2.GetCellVectors();
    v4 = v3Return[0];
    //    cout << v4 << endl;
    v5 = v3Return[1];
    //    cout << v5 << endl;
    v6 = v3Return[2];
    //    cout << v6 << endl;

    v9 = vector3(1.0f, 1.0f, 1.0f);
    v9 = cell2.FractionalToCartesian(v9);
    //    cout << v9 << endl;

    cout << "# Testing unit cell coordinate functions ..." << endl;
    cout << "13..24" << endl;

    // Cartesian wrapping
    results.getline(buffer,BUFF_SIZE);
    tokenize(vs,buffer);
    if (vs.size() != 3)
      {
	cout << "Bail out! Cannot parse results file 'unitcell_results.txt'" 
	     << endl;
	return(-1);
      }

    tmpcoords = cell2.WrapCartesianCoordinate(coords1);
    // cout << "Wrapped cartesian : " << tmpcoords << endl;
    if ( IsNear(tmpcoords[0], atof(vs[0].c_str()), 1.0e-3) )
      cout << "ok " << currTest++ << " # wrapped cartesian x" << endl;
    else
      cout << "not ok " << currTest++ << " # wrapped cartesian x" << endl;
    if ( IsNear(tmpcoords[1], atof(vs[1].c_str()), 1.0e-3) )
      cout << "ok " << currTest++ << " # wrapped cartesian y" << endl;
    else
      cout << "not ok " << currTest++ << " # wrapped cartesian y" << endl;
    if ( IsNear(tmpcoords[2], atof(vs[2].c_str()), 1.0e-3) )
      cout << "ok " << currTest++ << " # wrapped cartesian z" << endl;
    else
      cout << "not ok " << currTest++ << " # wrapped cartesian z" << endl;

    // Fractional wrapping
    results.getline(buffer,BUFF_SIZE);
    tokenize(vs,buffer);
    if (vs.size() != 3)
      {
	cout << "Bail out! Cannot parse results file 'unitcell_results.txt'" 
	     << endl;
	return(-1);
      }
    tmpcoords = cell2.WrapFractionalCoordinate(coords2);
    // cout << "Wrapped fractional: " << tmpcoords << endl;
    if ( IsNear(tmpcoords[0], atof(vs[0].c_str()), 1.0e-3) )
      cout << "ok " << currTest++ << " # wrapped fractional x" << endl;
    else
      cout << "not ok " << currTest++ << " # wrapped fractional x" << endl;
    if ( IsNear(tmpcoords[1], atof(vs[1].c_str()), 1.0e-3) )
      cout << "ok " << currTest++ << " # wrapped fractional y" << endl;
    else
      cout << "not ok " << currTest++ << " # wrapped fractional y" << endl;
    if ( IsNear(tmpcoords[2], atof(vs[2].c_str()), 1.0e-3) )
      cout << "ok " << currTest++ << " # wrapped fractional z" << endl;
    else
      cout << "not ok " << currTest++ << " # wrapped fractional z" << endl;

    // cart2frac
    results.getline(buffer,BUFF_SIZE);
    tokenize(vs,buffer);
    if (vs.size() != 3)
      {
	cout << "Bail out! Cannot parse results file 'unitcell_results.txt'" 
	     << endl;
	return(-1);
      }
    tmpcoords = cell2.CartesianToFractional(coords1);
    // cout << "cart2frac: " << tmpcoords << endl;
    if ( IsNear(tmpcoords[0], atof(vs[0].c_str()), 1.0e-3) )
      cout << "ok " << currTest++ << " # cart2frac x" << endl;
    else
      cout << "not ok " << currTest++ << " # cart2frac x" << endl;
    if ( IsNear(tmpcoords[1], atof(vs[1].c_str()), 1.0e-3) )
      cout << "ok " << currTest++ << " # cart2frac y" << endl;
    else
      cout << "not ok " << currTest++ << " # cart2frac y" << endl;
    if ( IsNear(tmpcoords[2], atof(vs[2].c_str()), 1.0e-3) )
      cout << "ok " << currTest++ << " # cart2frac z" << endl;
    else
      cout << "not ok " << currTest++ << " # cart2frac z" << endl;

    // frac2cart
    results.getline(buffer,BUFF_SIZE);
    tokenize(vs,buffer);
    if (vs.size() != 3)
      {
	cout << "Bail out! Cannot parse results file 'unitcell_results.txt'" 
	     << endl;
	return(-1);
      }
    tmpcoords = cell2.FractionalToCartesian(coords1);
    // cout << "frac2cart: " << tmpcoords << endl;
    if ( IsNear(tmpcoords[0], atof(vs[0].c_str()), 1.0e-3) )
      cout << "ok " << currTest++ << " # frac2cart x" << endl;
    else
      cout << "not ok " << currTest++ << " # frac2cart x" << endl;
    if ( IsNear(tmpcoords[1], atof(vs[1].c_str()), 1.0e-3) )
      cout << "ok " << currTest++ << " # frac2cart y" << endl;
    else
      cout << "not ok " << currTest++ << " # frac2cart y" << endl;
    if ( IsNear(tmpcoords[2], atof(vs[2].c_str()), 1.0e-3) )
      cout << "ok " << currTest++ << " # frac2cart z" << endl;
    else
      cout << "not ok " << currTest++ << " # frac2cart z" << endl;

    cout << "# cell volume " << cell.GetCellVolume() << endl;
    cout << "# lattice type " << cell.GetLatticeType() << endl;
    cout << "# cell volume " << cell2.GetCellVolume() << endl;
    cout << "# lattice type " << cell2.GetLatticeType() << endl;

    return(0); // success
}
