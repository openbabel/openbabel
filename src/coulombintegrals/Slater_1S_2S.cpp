/**********************************************************************************
Slater_1S_2S.cpp - Compute electrostatic energy between 1S and 2S Slater orbitals

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
double Slater_1S_2S(double r, double xi, double xj)
{
    double S, rxi, rxj;

    rxi = rxj = S = 0;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0)
        {
            S = (7*xi)/16

            ;
        }
        else
        {
            S = (1/r)*((-240 + 240*exp(2*rxi) - 375*rxi - 270*pow(rxi, 2) - 115*pow(rxi, 3) -

                          30*pow(rxi, 4) - 4*pow(rxi, 5))/(240*exp(2*rxi))

                         );
        }

    }
    else
    {
        if (r == 0)
        {
            S = (xi*xj*(pow(xi, 4) + 5*pow(xi, 3)*xj + 10*pow(xi, 2)*pow(xj, 2) +

                        10*xi*pow(xj, 3) + 2*pow(xj, 4)))/(2*pow(xi + xj, 5))

            ;
        }
        else
        {
            S = (1/r)*((6*exp(2*(rxi + rxj))*pow(pow(rxi, 2) - pow(rxj, 2), 5) +

                          6*exp(2*rxj)*pow(rxj, 6)*

                          (-4*pow(rxi, 4) - pow(rxi, 5) - 5*pow(rxi, 2)*pow(rxj, 2) +

                           pow(rxj, 4) + rxi*pow(rxj, 4)) -

                          exp(2*rxi)*pow(rxi, 4)*

                          (pow(rxi, 6)*(6 + 9*rxj + 6*pow(rxj, 2) + 2*pow(rxj, 3)) -

                           3*pow(rxi, 4)*pow(rxj, 2)*

                           (10 + 15*rxj + 10*pow(rxj, 2) + 2*pow(rxj, 3)) +

                           3*pow(rxi, 2)*pow(rxj, 4)*

                           (20 + 33*rxj + 14*pow(rxj, 2) + 2*pow(rxj, 3)) -

                           pow(rxj, 6)*(84 + 63*rxj + 18*pow(rxj, 2) + 2*pow(rxj, 3))))/

                         (6*exp(2*(rxi + rxj))*pow(rxi - rxj, 5)*pow(rxi + rxj, 5))

                         );
        }

    }
    return S;
}


double Slater_2S_1S(double r, double xi, double xj)
{
    return Slater_1S_2S(r, xj, xi);
}

}// end namespace openbabel
