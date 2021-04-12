/**********************************************************************************
Slater_1S_1S.cpp - Compute electrostatic energy between two 1S Slater orbitals

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
double Slater_1S_1S(double r, double xi, double xj)
{
    double S, rxi, rxj;

    rxi = rxj = S = 0;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0)
        {
            S = (5*xi)/8

            ;
        }
        else
        {
            S = (1/r)*((-24 + 24*exp(2*rxi) - 33*rxi - 18*pow(rxi, 2) - 4*pow(rxi, 3))/

                         (24*exp(2*rxi))

                         );
        }

    }
    else
    {
        if (r == 0)
        {
            S = (xi*xj*(pow(xi, 2) + 3*xi*xj + pow(xj, 2)))/pow(xi + xj, 3)

            ;
        }
        else
        {
            S = (1/r)*((exp(2*(rxi + rxj))*pow(pow(rxi, 2) - pow(rxj, 2), 3) +

                          exp(2*rxj)*pow(rxj, 4)*

                          (-3*pow(rxi, 2) - pow(rxi, 3) + pow(rxj, 2) + rxi*pow(rxj, 2)) -

                          exp(2*rxi)*pow(rxi, 4)*

                          (pow(rxi, 2)*(1 + rxj) - pow(rxj, 2)*(3 + rxj)))/

                         (exp(2*(rxi + rxj))*pow(rxi - rxj, 3)*pow(rxi + rxj, 3))

                         );
        }

    }
    return S;
}

} // end namespace openbabel
