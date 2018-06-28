/**********************************************************************************
Slater_2S_2S.cpp - Compute electrostatic energy between two 2S Slater orbitals

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

double Slater_2S_2S(double r, double xi, double xj)
{
    double S, rxi, rxj;

    rxi = rxj = S = 0;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0)
        {
            S = (93*xi)/256

            ;
        }
        else
        {
            S = (1/r)*((-80640 + 80640*exp(2*rxi) - 131985*rxi - 102690*pow(rxi, 2) -

                          49980*pow(rxi, 3) - 16800*pow(rxi, 4) - 4032*pow(rxi, 5) -

                          672*pow(rxi, 6) - 64*pow(rxi, 7))/(80640*exp(2*rxi))

                         );
        }

    }
    else
    {
        if (r == 0)
        {
            S = (xi*xj*(pow(xi, 6) + 7*pow(xi, 5)*xj + 21*pow(xi, 4)*pow(xj, 2) +

                        35*pow(xi, 3)*pow(xj, 3) + 21*pow(xi, 2)*pow(xj, 4) +

                        7*xi*pow(xj, 5) + pow(xj, 6)))/(2*pow(xi + xj, 7))

            ;
        }
        else
        {
            S = (1/r)*((6*exp(2*(rxi + rxj))*pow(pow(rxi, 2) - pow(rxj, 2), 7) -

                          exp(2*rxi)*pow(rxi, 6)*

                          (21*pow(rxi, 4)*pow(rxj, 4)*(6 + 11*rxj + 2*pow(rxj, 2)) -

                           2*pow(rxj, 8)*(90 + 54*rxj + 12*pow(rxj, 2) + pow(rxj, 3)) +

                           pow(rxi, 8)*(6 + 9*rxj + 6*pow(rxj, 2) + 2*pow(rxj, 3)) +

                           pow(rxi, 2)*pow(rxj, 6)*

                           (-390 - 69*rxj + 18*pow(rxj, 2) + 4*pow(rxj, 3)) -

                           pow(rxi, 6)*pow(rxj, 2)*

                           (42 + 63*rxj + 42*pow(rxj, 2) + 4*pow(rxj, 3))) +

                          exp(2*rxj)*pow(rxj, 6)*

                          (-24*pow(rxi, 10) - 2*pow(rxi, 11) - 69*pow(rxi, 7)*pow(rxj, 2) +

                           6*pow(rxj, 8) + 9*rxi*pow(rxj, 8) +

                           4*pow(rxi, 9)*(-27 + pow(rxj, 2)) +

                           18*pow(rxi, 8)*(-10 + pow(rxj, 2)) +

                           6*pow(rxi, 2)*pow(rxj, 6)*(-7 + pow(rxj, 2)) -

                           42*pow(rxi, 4)*pow(rxj, 4)*(-3 + pow(rxj, 2)) +

                           pow(rxi, 3)*pow(rxj, 6)*(-63 + 2*pow(rxj, 2)) +

                           6*pow(rxi, 6)*pow(rxj, 2)*(-65 + 7*pow(rxj, 2)) +

                           pow(rxi, 5)*(231*pow(rxj, 4) - 4*pow(rxj, 6))))/

                         (6*exp(2*(rxi + rxj))*pow(rxi - rxj, 7)*pow(rxi + rxj, 7))

                         );
        }

    }
    return S;
}
