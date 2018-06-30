/**********************************************************************************
Slater_3S_3S.cpp - Compute electrostatic energy between two 3S Slater orbitals

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
double Slater_3S_3S(double r, double xi, double xj)
{
    double S, rxi, rxj;

    rxi = rxj = S = 0;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0)
        {
            S = (793*xi)/3072

            ;
        }
        else
        {
            S = (1/r)*((-1437004800 + 1437004800*exp(2*rxi) - 2503064025*rxi -

                          2132118450*pow(rxi, 2) - 1180664100*pow(rxi, 3) -

                          476506800*pow(rxi, 4) - 148856400*pow(rxi, 5) -

                          37255680*pow(rxi, 6) - 7603200*pow(rxi, 7) - 1267200*pow(rxi, 8) -

                          168960*pow(rxi, 9) - 16896*pow(rxi, 10) - 1024*pow(rxi, 11))/

                         (1.4370048e9*exp(2*rxi))

                         );
        }

    }
    else
    {
        if (r == 0)
        {
            S = (xi*xj*(pow(xi, 10) + 11*pow(xi, 9)*xj + 55*pow(xi, 8)*pow(xj, 2) +

                        165*pow(xi, 7)*pow(xj, 3) + 330*pow(xi, 6)*pow(xj, 4) +

                        462*pow(xi, 5)*pow(xj, 5) + 330*pow(xi, 4)*pow(xj, 6) +

                        165*pow(xi, 3)*pow(xj, 7) + 55*pow(xi, 2)*pow(xj, 8) +

                        11*xi*pow(xj, 9) + pow(xj, 10)))/(3*pow(xi + xj, 11))

            ;
        }
        else
        {
            S = (1/r)*((135*exp(2*(rxi + rxj))*pow(pow(rxi, 2) - pow(rxj, 2), 11) +

                          exp(2*rxj)*pow(rxj, 8)*

                          (-150*pow(rxi, 18) - 6*pow(rxi, 19) + 135*pow(rxj, 14) +

                           225*rxi*pow(rxj, 14) + 10*pow(rxi, 17)*(-165 + pow(rxj, 2)) -

                           30*pow(rxi, 16)*(330 + pow(rxj, 2)) +

                           45*pow(rxi, 3)*pow(rxj, 12)*(-55 + 2*pow(rxj, 2)) +

                           45*pow(rxi, 2)*pow(rxj, 12)*(-33 + 4*pow(rxj, 2)) +

                           pow(rxi, 9)*pow(rxj, 6)*

                           (234135 - 4950*pow(rxj, 2) - 34*pow(rxj, 4)) -

                           5*pow(rxi, 7)*pow(rxj, 8)*

                           (6237 - 1242*pow(rxj, 2) + 2*pow(rxj, 4)) +

                           3*pow(rxi, 5)*pow(rxj, 10)*

                           (4125 - 330*pow(rxj, 2) + 2*pow(rxj, 4)) +

                           15*pow(rxi, 4)*pow(rxj, 10)*

                           (495 - 132*pow(rxj, 2) + 2*pow(rxj, 4)) -

                           165*pow(rxi, 6)*pow(rxj, 8)*

                           (135 - 60*pow(rxj, 2) + 2*pow(rxj, 4)) -

                           5*pow(rxi, 13)*pow(rxj, 2)*

                           (43875 - 3438*pow(rxj, 2) + 22*pow(rxj, 4)) +

                           5*pow(rxi, 11)*pow(rxj, 4)*

                           (7695 - 2442*pow(rxj, 2) + 22*pow(rxj, 4)) +

                           15*pow(rxi, 8)*pow(rxj, 6)*

                           (-33 - 3564*pow(rxj, 2) + 26*pow(rxj, 4)) +

                           pow(rxi, 15)*(-32175 - 3690*pow(rxj, 2) + 34*pow(rxj, 4)) +

                           15*pow(rxi, 10)*pow(rxj, 4)*

                           (-32277 + 1364*pow(rxj, 2) + 66*pow(rxj, 4)) +

                           15*pow(rxi, 14)*(-3003 - 2932*pow(rxj, 2) + 94*pow(rxj, 4)) -

                           15*pow(rxi, 12)*pow(rxj, 2)*

                           (28119 - 5252*pow(rxj, 2) + 154*pow(rxj, 4))) +

                          exp(2*rxi)*pow(rxi, 8)*

                          (-5*pow(rxi, 2)*pow(rxj, 12)*

                           (-84357 - 43875*rxj - 8796*pow(rxj, 2) - 738*pow(rxj, 3) -

                            6*pow(rxj, 4) + 2*pow(rxj, 5)) -

                           3*pow(rxi, 14)*(45 + 75*rxj + 60*pow(rxj, 2) + 30*pow(rxj, 3) +

                                                 10*pow(rxj, 4) + 2*pow(rxj, 5)) -

                           55*pow(rxi, 8)*pow(rxj, 6)*

                           (-405 - 567*rxj - 972*pow(rxj, 2) - 90*pow(rxj, 3) +

                            18*pow(rxj, 4) + 2*pow(rxj, 5)) +

                           55*pow(rxi, 6)*pow(rxj, 8)*

                           (9 - 4257*rxj - 372*pow(rxj, 2) + 222*pow(rxj, 3) +

                            42*pow(rxj, 4) + 2*pow(rxj, 5)) +

                           3*pow(rxj, 14)*(15015 + 10725*rxj + 3300*pow(rxj, 2) +

                                                 550*pow(rxj, 3) + 50*pow(rxj, 4) + 2*pow(rxj, 5)) +

                           5*pow(rxi, 12)*pow(rxj, 2)*

                           (297 + 495*rxj + 396*pow(rxj, 2) + 198*pow(rxj, 3) +

                            66*pow(rxj, 4) + 2*pow(rxj, 5)) +

                           pow(rxi, 10)*pow(rxj, 4)*

                           (-7425 - 12375*rxj - 9900*pow(rxj, 2) - 6210*pow(rxj, 3) -

                            390*pow(rxj, 4) + 34*pow(rxj, 5)) -

                           pow(rxi, 4)*pow(rxj, 10)*

                           (-484155 + 38475*rxj + 78780*pow(rxj, 2) + 17190*pow(rxj, 3) +

                            1410*pow(rxj, 4) + 34*pow(rxj, 5))))/

                         (135*exp(2*(rxi + rxj))*pow(rxi - rxj, 11)*pow(rxi + rxj, 11))

                         );
        }

    }
    return S;
}
} // end namespace openbabel
