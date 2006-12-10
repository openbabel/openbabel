/**********************************************************************
Copyright (C) 2003-2005 by Geoffrey R. Hutchison
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>

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

namespace OpenBabel
{
bool SafeOpen(std::ifstream &fs, const char *filename);
bool SafeOpen(std::ofstream &fs, const char *filename);
}

using namespace std;
using namespace OpenBabel;

#ifdef TESTDATADIR
    string testdatadir = TESTDATADIR;
    string unitcell_file = testdatadir + "unitcell.txt";
    string results_file = testdatadir + "unitcell_results.txt";
#else
    string unitcell_file = "files/unitcell.txt";
    string results_file = "files/unitcell_results.txt";
#endif

int main(int argc,char *argv[])
{
    if (argc != 1)
    {
        cout << "Usage: unitcell" << endl;
        cout << "   Tests Open Babel unit cell conversions." << endl;
        return 0;
    }

    cout << "# Testing unit cell transformations ..." << endl;
    cout << "1..12" << endl;

    double a, b, c, alpha, beta, gamma;
    vector3 v1, v2, v3, v4, v5, v6, v7, v8, v9;
    double x = 0.0, y = 0.0, z = 0.0;
    char buffer[BUFF_SIZE];
    std::ifstream ifs, results;
    OBUnitCell cell, cell2;
    vector<vector3> v3Return;
    vector<string> vs;
    unsigned int currTest = 1;

    if (!SafeOpen(ifs, unitcell_file.c_str()))
      {
	cout << "Bail out! Couldn't open test file 'unitcell.txt'" << endl;
        return(-1);
      }
    if (!SafeOpen(results, results_file.c_str()))
      {
	cout << "Bail out! Couldn't open test file 'unitcell_results.txt'" << endl;
        return(-1);
      }

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
    matrix3x3 m = cell2.GetOrthoMatrix();
    v9 *= m;
    //    cout << v9 << endl;

    cout << "# cell volume " << cell.GetCellVolume() << endl;
    cout << "# lattice type " << cell.GetLatticeType() << endl;
    cout << "# cell volume " << cell2.GetCellVolume() << endl;
    cout << "# lattice type " << cell2.GetLatticeType() << endl;

    return(0); // success
}
