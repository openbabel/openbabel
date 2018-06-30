/**********************************************************************************
DSlater_1S_3S.cpp - Compute electrostatic force between 1S and 3S Slater orbitals

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
double DSlater_1S_3S(double r, double xi, double xj)
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
            S = -(-203175*xi + 241920*exp(2*r*xi)*xi - 328860*r*pow(xi, 2) -

                  253260*pow(r, 2)*pow(xi, 3) - 120960*pow(r, 3)*pow(xi, 4) -

                  38640*pow(r, 4)*pow(xi, 5) - 8064*pow(r, 5)*pow(xi, 6) -

                  896*pow(r, 6)*pow(xi, 7))/(120960*exp(2*r*xi)*r) +

                (-120960 + 120960*exp(2*r*xi) - 203175*r*xi -

                 164430*pow(r, 2)*pow(xi, 2) - 84420*pow(r, 3)*pow(xi, 3) -

                 30240*pow(r, 4)*pow(xi, 4) - 7728*pow(r, 5)*pow(xi, 5) -

                 1344*pow(r, 6)*pow(xi, 6) - 128*pow(r, 7)*pow(xi, 7))/

                (120960*exp(2*r*xi)*pow(r, 2)) +

                (xi*(-120960 + 120960*exp(2*r*xi) - 203175*r*xi -

                     164430*pow(r, 2)*pow(xi, 2) - 84420*pow(r, 3)*pow(xi, 3) -

                     30240*pow(r, 4)*pow(xi, 4) - 7728*pow(r, 5)*pow(xi, 5) -

                     1344*pow(r, 6)*pow(xi, 6) - 128*pow(r, 7)*pow(xi, 7)))/

                (60480*exp(2*r*xi)*r)

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
            S = (45*exp(2*r*(xi + xj))*pow(pow(xi, 2) - pow(xj, 2), 7) +

                 15*exp(2*r*xj)*pow(xj, 8)*

                 (-15*pow(xi, 6) - 3*r*pow(xi, 7) - 63*pow(xi, 4)*pow(xj, 2) -

                  7*r*pow(xi, 5)*pow(xj, 2) - 21*pow(xi, 2)*pow(xj, 4) +

                  7*r*pow(xi, 3)*pow(xj, 4) + 3*pow(xj, 6) + 3*r*xi*pow(xj, 6)) +

                 exp(2*r*xi)*pow(xi, 4)*

                 (-10*pow(xi, 2)*pow(xj, 8)*

                  (135 + 333*r*xj + 228*pow(r, 2)*pow(xj, 2) +

                   75*pow(r, 3)*pow(xj, 3) + 13*pow(r, 4)*pow(xj, 4) +

                   pow(r, 5)*pow(xj, 5)) +

                  2*pow(xj, 10)*(945 + 945*r*xj + 420*pow(r, 2)*pow(xj, 2) +

                                       105*pow(r, 3)*pow(xj, 3) + 15*pow(r, 4)*pow(xj, 4) +

                                       pow(r, 5)*pow(xj, 5)) -

                  pow(xi, 10)*(45 + 75*r*xj + 60*pow(r, 2)*pow(xj, 2) +

                                   30*pow(r, 3)*pow(xj, 3) + 10*pow(r, 4)*pow(xj, 4) +

                                   2*pow(r, 5)*pow(xj, 5)) +

                  5*pow(xi, 8)*pow(xj, 2)*

                  (63 + 105*r*xj + 84*pow(r, 2)*pow(xj, 2) +

                   42*pow(r, 3)*pow(xj, 3) + 14*pow(r, 4)*pow(xj, 4) +

                   2*pow(r, 5)*pow(xj, 5)) -

                  5*pow(xi, 6)*pow(xj, 4)*

                  (189 + 315*r*xj + 252*pow(r, 2)*pow(xj, 2) +

                   132*pow(r, 3)*pow(xj, 3) + 36*pow(r, 4)*pow(xj, 4) +

                   4*pow(r, 5)*pow(xj, 5)) +

                  5*pow(xi, 4)*pow(xj, 6)*

                  (315 + 513*r*xj + 468*pow(r, 2)*pow(xj, 2) +

                   204*pow(r, 3)*pow(xj, 3) + 44*pow(r, 4)*pow(xj, 4) +

                   4*pow(r, 5)*pow(xj, 5))))/

                (45*exp(2*r*(xi + xj))*pow(r, 2)*pow(xi - xj, 7)*pow(xi + xj, 7))

                + (2*(45*exp(2*r*(xi + xj))*pow(pow(xi, 2) - pow(xj, 2), 7) +

                        15*exp(2*r*xj)*pow(xj, 8)*

                        (-15*pow(xi, 6) - 3*r*pow(xi, 7) - 63*pow(xi, 4)*pow(xj, 2) -

                         7*r*pow(xi, 5)*pow(xj, 2) - 21*pow(xi, 2)*pow(xj, 4) +

                         7*r*pow(xi, 3)*pow(xj, 4) + 3*pow(xj, 6) + 3*r*xi*pow(xj, 6))

                        + exp(2*r*xi)*pow(xi, 4)*(-10*pow(xi, 2)*pow(xj, 8)*

                                                        (135 + 333*r*xj + 228*pow(r, 2)*pow(xj, 2) +

                                                         75*pow(r, 3)*pow(xj, 3) + 13*pow(r, 4)*pow(xj, 4) +

                                                         pow(r, 5)*pow(xj, 5)) +

                                                        2*pow(xj, 10)*(945 + 945*r*xj + 420*pow(r, 2)*pow(xj, 2) +

                                                                             105*pow(r, 3)*pow(xj, 3) + 15*pow(r, 4)*pow(xj, 4) +

                                                                             pow(r, 5)*pow(xj, 5)) -

                                                        pow(xi, 10)*(45 + 75*r*xj + 60*pow(r, 2)*pow(xj, 2) +

                                                                         30*pow(r, 3)*pow(xj, 3) + 10*pow(r, 4)*pow(xj, 4) +

                                                                         2*pow(r, 5)*pow(xj, 5)) +

                                                        5*pow(xi, 8)*pow(xj, 2)*

                                                        (63 + 105*r*xj + 84*pow(r, 2)*pow(xj, 2) +

                                                         42*pow(r, 3)*pow(xj, 3) + 14*pow(r, 4)*pow(xj, 4) +

                                                         2*pow(r, 5)*pow(xj, 5)) -

                                                        5*pow(xi, 6)*pow(xj, 4)*

                                                        (189 + 315*r*xj + 252*pow(r, 2)*pow(xj, 2) +

                                                         132*pow(r, 3)*pow(xj, 3) + 36*pow(r, 4)*pow(xj, 4) +

                                                         4*pow(r, 5)*pow(xj, 5)) +

                                                        5*pow(xi, 4)*pow(xj, 6)*

                                                        (315 + 513*r*xj + 468*pow(r, 2)*pow(xj, 2) +

                                                         204*pow(r, 3)*pow(xj, 3) + 44*pow(r, 4)*pow(xj, 4) +

                                                         4*pow(r, 5)*pow(xj, 5)))))/

                (45*exp(2*r*(xi + xj))*r*pow(xi - xj, 7)*pow(xi + xj, 6)) -

                (90*exp(2*r*(xi + xj))*(xi + xj)*pow(pow(xi, 2) - pow(xj, 2), 7) +

                 15*exp(2*r*xj)*pow(xj, 8)*

                 (-3*pow(xi, 7) - 7*pow(xi, 5)*pow(xj, 2) +

                  7*pow(xi, 3)*pow(xj, 4) + 3*xi*pow(xj, 6)) +

                 30*exp(2*r*xj)*pow(xj, 9)*

                 (-15*pow(xi, 6) - 3*r*pow(xi, 7) - 63*pow(xi, 4)*pow(xj, 2) -

                  7*r*pow(xi, 5)*pow(xj, 2) - 21*pow(xi, 2)*pow(xj, 4) +

                  7*r*pow(xi, 3)*pow(xj, 4) + 3*pow(xj, 6) + 3*r*xi*pow(xj, 6)) +

                 exp(2*r*xi)*pow(xi, 4)*

                 (-10*pow(xi, 2)*pow(xj, 8)*

                  (333*xj + 456*r*pow(xj, 2) + 225*pow(r, 2)*pow(xj, 3) +

                   52*pow(r, 3)*pow(xj, 4) + 5*pow(r, 4)*pow(xj, 5)) +

                  2*pow(xj, 10)*(945*xj + 840*r*pow(xj, 2) +

                                       315*pow(r, 2)*pow(xj, 3) + 60*pow(r, 3)*pow(xj, 4) +

                                       5*pow(r, 4)*pow(xj, 5)) -

                  pow(xi, 10)*(75*xj + 120*r*pow(xj, 2) +

                                   90*pow(r, 2)*pow(xj, 3) + 40*pow(r, 3)*pow(xj, 4) +

                                   10*pow(r, 4)*pow(xj, 5)) +

                  5*pow(xi, 8)*pow(xj, 2)*

                  (105*xj + 168*r*pow(xj, 2) + 126*pow(r, 2)*pow(xj, 3) +

                   56*pow(r, 3)*pow(xj, 4) + 10*pow(r, 4)*pow(xj, 5)) -

                  5*pow(xi, 6)*pow(xj, 4)*

                  (315*xj + 504*r*pow(xj, 2) + 396*pow(r, 2)*pow(xj, 3) +

                   144*pow(r, 3)*pow(xj, 4) + 20*pow(r, 4)*pow(xj, 5)) +

                  5*pow(xi, 4)*pow(xj, 6)*

                  (513*xj + 936*r*pow(xj, 2) + 612*pow(r, 2)*pow(xj, 3) +

                   176*pow(r, 3)*pow(xj, 4) + 20*pow(r, 4)*pow(xj, 5))) +

                 2*exp(2*r*xi)*pow(xi, 5)*

                 (-10*pow(xi, 2)*pow(xj, 8)*

                  (135 + 333*r*xj + 228*pow(r, 2)*pow(xj, 2) +

                   75*pow(r, 3)*pow(xj, 3) + 13*pow(r, 4)*pow(xj, 4) +

                   pow(r, 5)*pow(xj, 5)) +

                  2*pow(xj, 10)*(945 + 945*r*xj + 420*pow(r, 2)*pow(xj, 2) +

                                       105*pow(r, 3)*pow(xj, 3) + 15*pow(r, 4)*pow(xj, 4) +

                                       pow(r, 5)*pow(xj, 5)) -

                  pow(xi, 10)*(45 + 75*r*xj + 60*pow(r, 2)*pow(xj, 2) +

                                   30*pow(r, 3)*pow(xj, 3) + 10*pow(r, 4)*pow(xj, 4) +

                                   2*pow(r, 5)*pow(xj, 5)) +

                  5*pow(xi, 8)*pow(xj, 2)*

                  (63 + 105*r*xj + 84*pow(r, 2)*pow(xj, 2) +

                   42*pow(r, 3)*pow(xj, 3) + 14*pow(r, 4)*pow(xj, 4) +

                   2*pow(r, 5)*pow(xj, 5)) -

                  5*pow(xi, 6)*pow(xj, 4)*

                  (189 + 315*r*xj + 252*pow(r, 2)*pow(xj, 2) +

                   132*pow(r, 3)*pow(xj, 3) + 36*pow(r, 4)*pow(xj, 4) +

                   4*pow(r, 5)*pow(xj, 5)) +

                  5*pow(xi, 4)*pow(xj, 6)*

                  (315 + 513*r*xj + 468*pow(r, 2)*pow(xj, 2) +

                   204*pow(r, 3)*pow(xj, 3) + 44*pow(r, 4)*pow(xj, 4) +

                   4*pow(r, 5)*pow(xj, 5))))/

                (45*exp(2*r*(xi + xj))*r*pow(xi - xj, 7)*pow(xi + xj, 7))

            ;
        }

    }
    return S;
}


double DSlater_3S_1S(double r, double xi, double xj)
{
    return DSlater_1S_3S(r, xj, xi);
}

} //end namespace openbabel
