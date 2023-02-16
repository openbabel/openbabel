/**********************************************************************************
Slater_2S_3S.cpp - Compute electrostatic energy between 2S and 3S Slater orbitals

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
double Slater_2S_3S(double r, double xi, double xj)
{
    double S, rxi, rxj;

    rxi = rxj = S = 0;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0)
        {
            S = (451*xi)/1536

            ;
        }
        else
        {
            S = (1/r)*((-4354560 + 4354560*exp(2*rxi) - 7430535*rxi - 6151950*pow(rxi, 2) -

                          3275370*pow(rxi, 3) - 1251180*pow(rxi, 4) - 361368*pow(rxi, 5) -

                          80640*pow(rxi, 6) - 13824*pow(rxi, 7) - 1728*pow(rxi, 8) -

                          128*pow(rxi, 9))/(4.35456e6*exp(2*rxi))

                         );
        }

    }
    else
    {
        if (r == 0)
        {
            S = (xi*xj*(2*pow(xi, 8) + 18*pow(xi, 7)*xj + 72*pow(xi, 6)*pow(xj, 2) +

                        168*pow(xi, 5)*pow(xj, 3) + 252*pow(xi, 4)*pow(xj, 4) +

                        252*pow(xi, 3)*pow(xj, 5) + 108*pow(xi, 2)*pow(xj, 6) +

                        27*xi*pow(xj, 7) + 3*pow(xj, 8)))/(6*pow(xi + xj, 9))

            ;
        }
        else
        {
            S = (1/r)*((90*exp(2*(rxi + rxj))*pow(pow(rxi, 2) - pow(rxj, 2), 9) +

                          5*exp(2*rxj)*pow(rxj, 8)*

                          (-90*pow(rxi, 12) - 6*pow(rxi, 13) + 18*pow(rxj, 10) +

                           27*rxi*pow(rxj, 10) + 18*pow(rxi, 2)*pow(rxj, 8)*

                           (-9 + pow(rxj, 2)) - 162*pow(rxi, 4)*pow(rxj, 6)*

                           (-4 + pow(rxj, 2)) - 198*pow(rxi, 10)*(5 + pow(rxj, 2)) -

                           108*pow(rxi, 6)*pow(rxj, 4)*(36 + pow(rxj, 2)) +

                           2*pow(rxi, 5)*pow(rxj, 6)*(675 + pow(rxj, 2)) -

                           18*pow(rxi, 7)*pow(rxj, 4)*(-81 + 2*pow(rxj, 2)) +

                           3*pow(rxi, 3)*pow(rxj, 8)*(-81 + 2*pow(rxj, 2)) -

                           pow(rxi, 11)*(495 + 2*pow(rxj, 2)) +

                           9*pow(rxi, 9)*pow(rxj, 2)*(-233 + 4*pow(rxj, 2)) +

                           6*pow(rxi, 8)*pow(rxj, 2)*(-1063 + 90*pow(rxj, 2))) -

                          2*exp(2*rxi)*pow(rxi, 6)*

                          (-90*pow(rxi, 6)*pow(rxj, 6)*

                           (42 + 65*rxj + 76*pow(rxj, 2) + 22*pow(rxj, 3) + 2*pow(rxj, 4)) -

                           2*pow(rxj, 12)*(2970 + 2475*rxj + 900*pow(rxj, 2) +

                                                 180*pow(rxj, 3) + 20*pow(rxj, 4) + pow(rxj, 5)) +

                           10*pow(rxi, 8)*pow(rxj, 4)*

                           (162 + 270*rxj + 216*pow(rxj, 2) + 122*pow(rxj, 3) +

                            22*pow(rxj, 4) + pow(rxj, 5)) -

                           5*pow(rxi, 4)*pow(rxj, 8)*

                           (-639 - 3555*rxj - 1452*pow(rxj, 2) - 174*pow(rxj, 3) +

                            6*pow(rxj, 4) + 2*pow(rxj, 5)) +

                           pow(rxi, 12)*(45 + 75*rxj + 60*pow(rxj, 2) + 30*pow(rxj, 3) +

                                             10*pow(rxj, 4) + 2*pow(rxj, 5)) -

                           pow(rxi, 10)*pow(rxj, 2)*

                           (405 + 675*rxj + 540*pow(rxj, 2) + 270*pow(rxj, 3) +

                            90*pow(rxj, 4) + 8*pow(rxj, 5)) +

                           pow(rxi, 2)*pow(rxj, 10)*

                           (-21615 - 9075*rxj - 300*pow(rxj, 2) + 490*pow(rxj, 3) +

                            110*pow(rxj, 4) + 8*pow(rxj, 5))))/

                         (90*exp(2*(rxi + rxj))*pow(rxi - rxj, 9)*pow(rxi + rxj, 9))

                         );
        }

    }
    return S;
}


double Slater_3S_2S(double r, double xi, double xj)
{
    return Slater_2S_3S(r, xj, xi);
}

} // end namespace openbabel
