/**********************************************************************************
DSlater_1S_1S.cpp - Compute electrostatic force between two 1S Slater orbitals

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


#include <openbabel/slater_integrals.h>

namespace OpenBabel
{
double DSlater_1S_1S(double r, double xi, double xj)
{
    double S, rxi, rxj;

    rxi = rxj = S = 0;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0)
        {
            S = 0

            ;
        }
        else
        {
            S = -(-33*xi + 48*exp(2*r*xi)*xi - 36*r*pow(xi, 2) -

                  12*pow(r, 2)*pow(xi, 3))/(24*exp(2*r*xi)*r) +

                (-24 + 24*exp(2*r*xi) - 33*r*xi - 18*pow(r, 2)*pow(xi, 2) -

                 4*pow(r, 3)*pow(xi, 3))/(24*exp(2*r*xi)*pow(r, 2)) +

                (xi*(-24 + 24*exp(2*r*xi) - 33*r*xi - 18*pow(r, 2)*pow(xi, 2) -

                     4*pow(r, 3)*pow(xi, 3)))/(12*exp(2*r*xi)*r)

            ;
        }

    }
    else
    {
        if (r == 0)
        {
            S = 0

            ;
        }
        else
        {
            S = (exp(2*r*(xi + xj))*pow(pow(xi, 2) - pow(xj, 2), 3) +

                 exp(2*r*xj)*pow(xj, 4)*

                 (-3*pow(xi, 2) - r*pow(xi, 3) + pow(xj, 2) + r*xi*pow(xj, 2)) -

                 exp(2*r*xi)*pow(xi, 4)*

                 (pow(xi, 2)*(1 + r*xj) - pow(xj, 2)*(3 + r*xj)))/

                (exp(2*r*(xi + xj))*pow(r, 2)*pow(xi - xj, 3)*pow(xi + xj, 3)) +

                (2*(exp(2*r*(xi + xj))*pow(pow(xi, 2) - pow(xj, 2), 3) +

                      exp(2*r*xj)*pow(xj, 4)*

                      (-3*pow(xi, 2) - r*pow(xi, 3) + pow(xj, 2) + r*xi*pow(xj, 2)) -

                      exp(2*r*xi)*pow(xi, 4)*

                      (pow(xi, 2)*(1 + r*xj) - pow(xj, 2)*(3 + r*xj))))/

                (exp(2*r*(xi + xj))*r*pow(xi - xj, 3)*pow(xi + xj, 2)) -

                (2*exp(2*r*(xi + xj))*(xi + xj)*pow(pow(xi, 2) - pow(xj, 2), 3) +

                 exp(2*r*xj)*pow(xj, 4)*(-pow(xi, 3) + xi*pow(xj, 2)) +

                 2*exp(2*r*xj)*pow(xj, 5)*

                 (-3*pow(xi, 2) - r*pow(xi, 3) + pow(xj, 2) + r*xi*pow(xj, 2)) -

                 exp(2*r*xi)*pow(xi, 4)*(pow(xi, 2)*xj - pow(xj, 3)) -

                 2*exp(2*r*xi)*pow(xi, 5)*

                 (pow(xi, 2)*(1 + r*xj) - pow(xj, 2)*(3 + r*xj)))/

                (exp(2*r*(xi + xj))*r*pow(xi - xj, 3)*pow(xi + xj, 3))

            ;
        }

    }
    return S;
}

} // end namespace openbabel
