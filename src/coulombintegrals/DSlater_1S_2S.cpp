/**********************************************************************************
DSlater_1S_2S.cpp - Compute electrostatic force between 1S and 2S Slater orbitals

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

double DSlater_1S_2S(double r, double xi, double xj)
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
            S = -(-375*xi + 480*exp(2*r*xi)*xi - 540*r*pow(xi, 2) -

                  345*pow(r, 2)*pow(xi, 3) - 120*pow(r, 3)*pow(xi, 4) -

                  20*pow(r, 4)*pow(xi, 5))/(240*exp(2*r*xi)*r) +

                (-240 + 240*exp(2*r*xi) - 375*r*xi - 270*pow(r, 2)*pow(xi, 2) -

                 115*pow(r, 3)*pow(xi, 3) - 30*pow(r, 4)*pow(xi, 4) -

                 4*pow(r, 5)*pow(xi, 5))/(240*exp(2*r*xi)*pow(r, 2)) +

                (xi*(-240 + 240*exp(2*r*xi) - 375*r*xi - 270*pow(r, 2)*pow(xi, 2) -

                     115*pow(r, 3)*pow(xi, 3) - 30*pow(r, 4)*pow(xi, 4) -

                     4*pow(r, 5)*pow(xi, 5)))/(120*exp(2*r*xi)*r)

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
            S = (6*exp(2*r*(xi + xj))*pow(pow(xi, 2) - pow(xj, 2), 5) +

                 6*exp(2*r*xj)*pow(xj, 6)*

                 (-4*pow(xi, 4) - r*pow(xi, 5) - 5*pow(xi, 2)*pow(xj, 2) +

                  pow(xj, 4) + r*xi*pow(xj, 4)) -

                 exp(2*r*xi)*pow(xi, 4)*

                 (pow(xi, 6)*(6 + 9*r*xj + 6*pow(r, 2)*pow(xj, 2) +

                                  2*pow(r, 3)*pow(xj, 3)) -

                  3*pow(xi, 4)*pow(xj, 2)*

                  (10 + 15*r*xj + 10*pow(r, 2)*pow(xj, 2) +

                   2*pow(r, 3)*pow(xj, 3)) +

                  3*pow(xi, 2)*pow(xj, 4)*

                  (20 + 33*r*xj + 14*pow(r, 2)*pow(xj, 2) +

                   2*pow(r, 3)*pow(xj, 3)) -

                  pow(xj, 6)*(84 + 63*r*xj + 18*pow(r, 2)*pow(xj, 2) +

                                  2*pow(r, 3)*pow(xj, 3))))/

                (6*exp(2*r*(xi + xj))*pow(r, 2)*pow(xi - xj, 5)*pow(xi + xj, 5)) +

                (6*exp(2*r*(xi + xj))*pow(pow(xi, 2) - pow(xj, 2), 5) +

                 6*exp(2*r*xj)*pow(xj, 6)*

                 (-4*pow(xi, 4) - r*pow(xi, 5) - 5*pow(xi, 2)*pow(xj, 2) +

                  pow(xj, 4) + r*xi*pow(xj, 4)) -

                 exp(2*r*xi)*pow(xi, 4)*

                 (pow(xi, 6)*(6 + 9*r*xj + 6*pow(r, 2)*pow(xj, 2) +

                                  2*pow(r, 3)*pow(xj, 3)) -

                  3*pow(xi, 4)*pow(xj, 2)*

                  (10 + 15*r*xj + 10*pow(r, 2)*pow(xj, 2) +

                   2*pow(r, 3)*pow(xj, 3)) +

                  3*pow(xi, 2)*pow(xj, 4)*

                  (20 + 33*r*xj + 14*pow(r, 2)*pow(xj, 2) +

                   2*pow(r, 3)*pow(xj, 3)) -

                  pow(xj, 6)*(84 + 63*r*xj + 18*pow(r, 2)*pow(xj, 2) +

                                  2*pow(r, 3)*pow(xj, 3))))/

                (3*exp(2*r*(xi + xj))*r*pow(xi - xj, 5)*pow(xi + xj, 4)) -

                (12*exp(2*r*(xi + xj))*(xi + xj)*pow(pow(xi, 2) - pow(xj, 2), 5) +

                 6*exp(2*r*xj)*pow(xj, 6)*(-pow(xi, 5) + xi*pow(xj, 4)) +

                 12*exp(2*r*xj)*pow(xj, 7)*

                 (-4*pow(xi, 4) - r*pow(xi, 5) - 5*pow(xi, 2)*pow(xj, 2) +

                  pow(xj, 4) + r*xi*pow(xj, 4)) -

                 exp(2*r*xi)*pow(xi, 4)*

                 (pow(xi, 6)*(9*xj + 12*r*pow(xj, 2) + 6*pow(r, 2)*pow(xj, 3)) -

                  3*pow(xi, 4)*pow(xj, 2)*

                  (15*xj + 20*r*pow(xj, 2) + 6*pow(r, 2)*pow(xj, 3)) +

                  3*pow(xi, 2)*pow(xj, 4)*

                  (33*xj + 28*r*pow(xj, 2) + 6*pow(r, 2)*pow(xj, 3)) -

                  pow(xj, 6)*(63*xj + 36*r*pow(xj, 2) + 6*pow(r, 2)*pow(xj, 3))) -

                 2*exp(2*r*xi)*pow(xi, 5)*

                 (pow(xi, 6)*(6 + 9*r*xj + 6*pow(r, 2)*pow(xj, 2) +

                                  2*pow(r, 3)*pow(xj, 3)) -

                  3*pow(xi, 4)*pow(xj, 2)*

                  (10 + 15*r*xj + 10*pow(r, 2)*pow(xj, 2) +

                   2*pow(r, 3)*pow(xj, 3)) +

                  3*pow(xi, 2)*pow(xj, 4)*

                  (20 + 33*r*xj + 14*pow(r, 2)*pow(xj, 2) +

                   2*pow(r, 3)*pow(xj, 3)) -

                  pow(xj, 6)*(84 + 63*r*xj + 18*pow(r, 2)*pow(xj, 2) +

                                  2*pow(r, 3)*pow(xj, 3))))/

                (6*exp(2*r*(xi + xj))*r*pow(xi - xj, 5)*pow(xi + xj, 5))

            ;
        }

    }
    return S;
}


double DSlater_2S_1S(double r, double xi, double xj)
{
    return DSlater_1S_2S(r, xj, xi);
}

