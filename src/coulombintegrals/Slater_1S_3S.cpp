/**********************************************************************************
Slater_1S_3S.cpp - Compute electrostatic energy between 1S and 3S Slater orbitals

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

double Slater_1S_3S(double r, double xi, double xj)
{
    double S, rxi, rxj;

    rxi = rxj = S = 0;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0)
        {
            S = (41*xi)/128

            ;
        }
        else
        {
            S = (1/r)*((-120960 + 120960*exp(2*rxi) - 203175*rxi - 164430*pow(rxi, 2) -

                          84420*pow(rxi, 3) - 30240*pow(rxi, 4) - 7728*pow(rxi, 5) -

                          1344*pow(rxi, 6) - 128*pow(rxi, 7))/(120960*exp(2*rxi))

                         );
        }

    }
    else
    {
        if (r == 0)
        {
            S = (xi*xj*(pow(xi, 6) + 7*pow(xi, 5)*xj + 21*pow(xi, 4)*pow(xj, 2) +

                        35*pow(xi, 3)*pow(xj, 3) + 35*pow(xi, 2)*pow(xj, 4) +

                        21*xi*pow(xj, 5) + 3*pow(xj, 6)))/(3*pow(xi + xj, 7))

            ;
        }
        else
        {
            S = (1/r)*((45*exp(2*(rxi + rxj))*pow(pow(rxi, 2) - pow(rxj, 2), 7) +

                          15*exp(2*rxj)*pow(rxj, 8)*

                          (-15*pow(rxi, 6) - 3*pow(rxi, 7) - 63*pow(rxi, 4)*pow(rxj, 2) -

                           7*pow(rxi, 5)*pow(rxj, 2) - 21*pow(rxi, 2)*pow(rxj, 4) +

                           7*pow(rxi, 3)*pow(rxj, 4) + 3*pow(rxj, 6) + 3*rxi*pow(rxj, 6)) +

                          exp(2*rxi)*pow(rxi, 4)*

                          (-10*pow(rxi, 2)*pow(rxj, 8)*

                           (135 + 333*rxj + 228*pow(rxj, 2) + 75*pow(rxj, 3) +

                            13*pow(rxj, 4) + pow(rxj, 5)) +

                           2*pow(rxj, 10)*(945 + 945*rxj + 420*pow(rxj, 2) +

                                                 105*pow(rxj, 3) + 15*pow(rxj, 4) + pow(rxj, 5)) -

                           pow(rxi, 10)*(45 + 75*rxj + 60*pow(rxj, 2) + 30*pow(rxj, 3) +

                                             10*pow(rxj, 4) + 2*pow(rxj, 5)) +

                           5*pow(rxi, 8)*pow(rxj, 2)*

                           (63 + 105*rxj + 84*pow(rxj, 2) + 42*pow(rxj, 3) +

                            14*pow(rxj, 4) + 2*pow(rxj, 5)) -

                           5*pow(rxi, 6)*pow(rxj, 4)*

                           (189 + 315*rxj + 252*pow(rxj, 2) + 132*pow(rxj, 3) +

                            36*pow(rxj, 4) + 4*pow(rxj, 5)) +

                           5*pow(rxi, 4)*pow(rxj, 6)*

                           (315 + 513*rxj + 468*pow(rxj, 2) + 204*pow(rxj, 3) +

                            44*pow(rxj, 4) + 4*pow(rxj, 5))))/

                         (45*exp(2*(rxi + rxj))*pow(rxi - rxj, 7)*pow(rxi + rxj, 7))

                         );
        }

    }
    return S;
}


double Slater_3S_1S(double r, double xi, double xj)
{
    return Slater_1S_3S(r, xj, xi);
}

