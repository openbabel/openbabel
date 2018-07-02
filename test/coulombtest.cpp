/**********************************************************************************
coulombtest.cpp - Test the Slater and Gaussian integrals

Copyright (C) 2018 by Mohammad M Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
                      David van der Spoel Group
                      Uppsala University

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************************/

#include "obtest.h"

#include <openbabel/slater_integrals.h>
#include <openbabel/gaussian_integrals.h>

using namespace std;
using namespace OpenBabel;

enum CoulombTestType{
    eEnergy,
    eForce 
};

static int test_Slater(double r, int  i, int  j, 
                       double  xi, double  xj, 
                       double ref_value, int test_number, 
                       CoulombTestType ctt)
{
    double test_value = 0;
    if (ctt == eEnergy)
        test_value = Coulomb_SS(r, i, j, xi, xj);   
    else
        test_value = -DCoulomb_SS(r, i, j, xi, xj);
        
    if (IsNear(test_value, ref_value, 1e-10))
    {
        cout << "ok Slater test " << test_number << endl;
        return(0);
    }
    else
    {
        cout << "not ok Slater test " << test_number  << endl;
        cout << "Test value: " << test_value << "Ref Value: " << ref_value << endl;
        return(-1);
    }
}


static int test_Gaussian(double r, double  xi, double  xj, 
                         double ref_value, int test_number, 
                         CoulombTestType ctt)
{
    double test_value = 0;
    if (ctt == eEnergy)
        test_value = Coulomb_GG(r, xi, xj);   
    else
        test_value = -DCoulomb_GG(r, xi, xj);
           
    if (IsNear(test_value, ref_value, 1e-10))
    {
        cout << "ok Slater test " << test_number << endl;
        return(0);
    }
    else
    {
      cout << "not ok Slater test " << test_number  << endl;
      cout << "Test value: " << test_value << "Ref Value: " << ref_value << endl;
      return(-1);
    }
}

int coulombtest(int argc, char* argv[])
{
  int defaultchoice = 1;  
  int choice = defaultchoice;

  if (argc > 1) {
    if(sscanf(argv[1], "%d", &choice) != 1) {
      printf("Couldn't parse that input as a number\n");
      return -1;
    }
  }

  switch(choice) {
  case 1:
      return test_Slater(0.1, 1, 2, 12.0, 16.0, 5.7637290944563961, choice, eEnergy);
      break;
  case 2:
      return test_Slater(0.1, 1, 3, 12.0, 17.0, 4.86819764284163, choice, eEnergy);
      break;    
  case 3:
      return test_Slater(0.1, 1, 3, 12.0, 17.0, -5.7309185468702326, choice, eForce);
      break;
  case 4:
      return test_Slater(0, 3, 3, 11.0, 8.0, 2.3476310210961322, choice, eEnergy); // Slater 3s-3s energy at r = 0
      break;    
  case 5:
      return test_Slater(0, 3, 3, 11.0, 8.0, -0, choice, eForce); //Slater 3s-3s force at r = 0
      break;
  case 6:
      return test_Slater(2, 3, 3, 11.0, 8.0, 0.49999998663716444, choice, eEnergy);
      break;    
  case 7:
      return test_Slater(2, 3, 3, 11.0, 8.0, -0.24999982269334378, choice, eForce);
      break;
  case 8:
      return test_Gaussian(0.1, 7.0, 7.0, 5.1607269555385402, choice, eEnergy);
      break;
  case 9:
      return test_Gaussian(0.1, 7.0, 7.0, -7.8917188840388235, choice, eForce);
      break;
  case 10:
      return test_Gaussian(0, 5.0, 8.0, 4.7843180998583428, choice, eEnergy); //Gaussian energy at r = 0
      break;
  case 11:
      return test_Gaussian(0, 5.0, 8.0, -0, choice, eForce); //Gaussian force at r = 0
      break;
  default:
      cout << "Test number " << choice << " does not exist!\n";
      return -1;
  }
  return 0;
}

