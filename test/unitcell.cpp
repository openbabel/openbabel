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

#include "mol.h"
#include "math/matrix3x3.h"

namespace OpenBabel
{
bool SafeOpen(std::ifstream &fs, char *filename);
bool SafeOpen(std::ofstream &fs, char *filename);
}

using namespace std;
using namespace OpenBabel;

bool TestUnitCell()
{
    double a, b, c, alpha, beta, gamma;
    vector3 v1, v2, v3, v4, v5, v6, v7, v8, v9;
    double x = 0.0, y = 0.0, z = 0.0;
    char buffer[BUFF_SIZE];
    std::ifstream ifs;
    OBUnitCell cell, cell2;
    vector<vector3> v3Return;

#ifdef TESTDATADIR

    string testdatadir = TESTDATADIR;
    string unitcell_file = testdatadir + "unitcell.txt";
#else

    string unitcell_file = "unitcell.txt";
#endif

    if (!SafeOpen(ifs, (char*)unitcell_file.c_str()))
        return(false);
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

    //    cerr << " a: " << a << " b: " << b << " c: " << c << " alpha: " << alpha
    //	 << " beta: " << beta << " gamma: " << gamma << endl;

    cell2.SetData(a, b, c, alpha, beta, gamma);
    v3Return = cell2.GetCellVectors();
    v4 = v3Return[0];
    v5 = v3Return[1];
    v6 = v3Return[2];

    if (cell.GetA() != cell2.GetA() || cell.GetB() != cell2.GetB() ||
	cell.GetC() != cell2.GetC() || cell.GetAlpha() != cell2.GetAlpha() ||
	cell.GetBeta() != cell2.GetBeta() || 
	cell.GetGamma() != cell2.GetGamma())
      return(false);

    //    cerr << "x: " << v1.x() << " y: " << v1.y() << " z: " << v1.z() << " " <<
    //      v1.length() << endl;
    //    cerr << "x: " << v2.x() << " y: " << v2.y() << " z: " << v2.z() << " " <<
    //  v2.length() << endl;
    //    cerr << "x: " << v3.x() << " y: " << v3.y() << " z: " << v3.z() << " " <<
    //  v3.length() << endl;
    //    cerr << endl;

    //    cerr << "x: " << v4.x() << " y: " << v4.y() << " z: " << v4.z() << " " <<
    //  v4.length() << endl;
    //    cerr << "x: " << v5.x() << " y: " << v5.y() << " z: " << v5.z() << " " <<
    //  v5.length() << endl;
    //    cerr << "x: " << v6.x() << " y: " << v6.y() << " z: " << v6.z() << " " <<
    //  v6.length() << endl;
    //    cerr << endl;

    v9 = vector3(1.0f, 1.0f, 1.0f);
    matrix3x3 m = cell2.GetOrthoMatrix();
    v9 *= m;

    //    cerr << "x: " << v9.x() << " y: " << v9.y() << " z: " << v9.z() << " " <<
    //  v9.length() << endl;    
    //    cerr << endl;

    const vector3 testV    ( 0.5, 1.0, 1.5 ) ;

    v7 = testV + v1 + v2 + v3;
    v8 = testV + v4 + v5 + v6;

    //    cerr << "x: " << v7.x() << " y: " << v7.y() << " z: " << v7.z() << " " <<
    //      v7.length() << endl;
    //    cerr << "x: " << v8.x() << " y: " << v8.y() << " z: " << v8.z() << " " <<
    //     v8.length() << endl;

    return(true);
}
