/**********************************************************************************
DSlater_2S_3S.cpp - Compute electrostatic force between 2S and 3S Slater orbitals

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

double DSlater_2S_3S(double r, double xi, double xj)
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
            S = -(-7430535*xi + 8709120*exp(2*r*xi)*xi - 12303900*r*pow(xi, 2) -

                  9826110*pow(r, 2)*pow(xi, 3) - 5004720*pow(r, 3)*pow(xi, 4) -

                  1806840*pow(r, 4)*pow(xi, 5) - 483840*pow(r, 5)*pow(xi, 6) -

                  96768*pow(r, 6)*pow(xi, 7) - 13824*pow(r, 7)*pow(xi, 8) -

                  1152*pow(r, 8)*pow(xi, 9))/(4.35456e6*exp(2*r*xi)*r) +

                (-4354560 + 4354560*exp(2*r*xi) - 7430535*r*xi -

                 6151950*pow(r, 2)*pow(xi, 2) - 3275370*pow(r, 3)*pow(xi, 3) -

                 1251180*pow(r, 4)*pow(xi, 4) - 361368*pow(r, 5)*pow(xi, 5) -

                 80640*pow(r, 6)*pow(xi, 6) - 13824*pow(r, 7)*pow(xi, 7) -

                 1728*pow(r, 8)*pow(xi, 8) - 128*pow(r, 9)*pow(xi, 9))/

                (4.35456e6*exp(2*r*xi)*pow(r, 2)) +

                (xi*(-4354560 + 4354560*exp(2*r*xi) - 7430535*r*xi -

                     6151950*pow(r, 2)*pow(xi, 2) - 3275370*pow(r, 3)*pow(xi, 3) -

                     1251180*pow(r, 4)*pow(xi, 4) - 361368*pow(r, 5)*pow(xi, 5) -

                     80640*pow(r, 6)*pow(xi, 6) - 13824*pow(r, 7)*pow(xi, 7) -

                     1728*pow(r, 8)*pow(xi, 8) - 128*pow(r, 9)*pow(xi, 9)))/

                (2.17728e6*exp(2*r*xi)*r)

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
            S = (90*exp(2*r*(xi + xj))*pow(pow(xi, 2) - pow(xj, 2), 9) +

                 5*exp(2*r*xj)*pow(xj, 8)*

                 (-90*pow(r, 2)*pow(xi, 12) - 6*pow(r, 3)*pow(xi, 13) +

                  18*pow(xj, 10) + 27*r*xi*pow(xj, 10) +

                  18*pow(xi, 2)*pow(xj, 8)*(-9 + pow(r, 2)*pow(xj, 2)) -

                  162*pow(xi, 4)*pow(xj, 6)*(-4 + pow(r, 2)*pow(xj, 2)) -

                  198*pow(xi, 10)*(5 + pow(r, 2)*pow(xj, 2)) -

                  108*pow(xi, 6)*pow(xj, 4)*(36 + pow(r, 2)*pow(xj, 2)) +

                  2*r*pow(xi, 5)*pow(xj, 6)*(675 + pow(r, 2)*pow(xj, 2)) -

                  18*r*pow(xi, 7)*pow(xj, 4)*(-81 + 2*pow(r, 2)*pow(xj, 2)) +

                  3*r*pow(xi, 3)*pow(xj, 8)*(-81 + 2*pow(r, 2)*pow(xj, 2)) -

                  r*pow(xi, 11)*(495 + 2*pow(r, 2)*pow(xj, 2)) +

                  9*r*pow(xi, 9)*pow(xj, 2)*(-233 + 4*pow(r, 2)*pow(xj, 2)) +

                  6*pow(xi, 8)*pow(xj, 2)*(-1063 + 90*pow(r, 2)*pow(xj, 2))) -

                 2*exp(2*r*xi)*pow(xi, 6)*

                 (-90*pow(xi, 6)*pow(xj, 6)*

                  (42 + 65*r*xj + 76*pow(r, 2)*pow(xj, 2) +

                   22*pow(r, 3)*pow(xj, 3) + 2*pow(r, 4)*pow(xj, 4)) -

                  2*pow(xj, 12)*(2970 + 2475*r*xj + 900*pow(r, 2)*pow(xj, 2) +

                                       180*pow(r, 3)*pow(xj, 3) + 20*pow(r, 4)*pow(xj, 4) +

                                       pow(r, 5)*pow(xj, 5)) +

                  10*pow(xi, 8)*pow(xj, 4)*

                  (162 + 270*r*xj + 216*pow(r, 2)*pow(xj, 2) +

                   122*pow(r, 3)*pow(xj, 3) + 22*pow(r, 4)*pow(xj, 4) +

                   pow(r, 5)*pow(xj, 5)) -

                  5*pow(xi, 4)*pow(xj, 8)*

                  (-639 - 3555*r*xj - 1452*pow(r, 2)*pow(xj, 2) -

                   174*pow(r, 3)*pow(xj, 3) + 6*pow(r, 4)*pow(xj, 4) +

                   2*pow(r, 5)*pow(xj, 5)) +

                  pow(xi, 12)*(45 + 75*r*xj + 60*pow(r, 2)*pow(xj, 2) +

                                   30*pow(r, 3)*pow(xj, 3) + 10*pow(r, 4)*pow(xj, 4) +

                                   2*pow(r, 5)*pow(xj, 5)) -

                  pow(xi, 10)*pow(xj, 2)*

                  (405 + 675*r*xj + 540*pow(r, 2)*pow(xj, 2) +

                   270*pow(r, 3)*pow(xj, 3) + 90*pow(r, 4)*pow(xj, 4) +

                   8*pow(r, 5)*pow(xj, 5)) +

                  pow(xi, 2)*pow(xj, 10)*

                  (-21615 - 9075*r*xj - 300*pow(r, 2)*pow(xj, 2) +

                   490*pow(r, 3)*pow(xj, 3) + 110*pow(r, 4)*pow(xj, 4) +

                   8*pow(r, 5)*pow(xj, 5))))/

                (90*exp(2*r*(xi + xj))*pow(r, 2)*pow(xi - xj, 9)*pow(xi + xj, 9))

                + (90*exp(2*r*(xi + xj))*pow(pow(xi, 2) - pow(xj, 2), 9) +

                   5*exp(2*r*xj)*pow(xj, 8)*

                   (-90*pow(r, 2)*pow(xi, 12) - 6*pow(r, 3)*pow(xi, 13) +

                    18*pow(xj, 10) + 27*r*xi*pow(xj, 10) +

                    18*pow(xi, 2)*pow(xj, 8)*(-9 + pow(r, 2)*pow(xj, 2)) -

                    162*pow(xi, 4)*pow(xj, 6)*(-4 + pow(r, 2)*pow(xj, 2)) -

                    198*pow(xi, 10)*(5 + pow(r, 2)*pow(xj, 2)) -

                    108*pow(xi, 6)*pow(xj, 4)*(36 + pow(r, 2)*pow(xj, 2)) +

                    2*r*pow(xi, 5)*pow(xj, 6)*(675 + pow(r, 2)*pow(xj, 2)) -

                    18*r*pow(xi, 7)*pow(xj, 4)*(-81 + 2*pow(r, 2)*pow(xj, 2)) +

                    3*r*pow(xi, 3)*pow(xj, 8)*(-81 + 2*pow(r, 2)*pow(xj, 2)) -

                    r*pow(xi, 11)*(495 + 2*pow(r, 2)*pow(xj, 2)) +

                    9*r*pow(xi, 9)*pow(xj, 2)*(-233 + 4*pow(r, 2)*pow(xj, 2)) +

                    6*pow(xi, 8)*pow(xj, 2)*(-1063 + 90*pow(r, 2)*pow(xj, 2))) -

                   2*exp(2*r*xi)*pow(xi, 6)*

                   (-90*pow(xi, 6)*pow(xj, 6)*

                    (42 + 65*r*xj + 76*pow(r, 2)*pow(xj, 2) +

                     22*pow(r, 3)*pow(xj, 3) + 2*pow(r, 4)*pow(xj, 4)) -

                    2*pow(xj, 12)*(2970 + 2475*r*xj + 900*pow(r, 2)*pow(xj, 2) +

                                         180*pow(r, 3)*pow(xj, 3) + 20*pow(r, 4)*pow(xj, 4) +

                                         pow(r, 5)*pow(xj, 5)) +

                    10*pow(xi, 8)*pow(xj, 4)*

                    (162 + 270*r*xj + 216*pow(r, 2)*pow(xj, 2) +

                     122*pow(r, 3)*pow(xj, 3) + 22*pow(r, 4)*pow(xj, 4) +

                     pow(r, 5)*pow(xj, 5)) -

                    5*pow(xi, 4)*pow(xj, 8)*

                    (-639 - 3555*r*xj - 1452*pow(r, 2)*pow(xj, 2) -

                     174*pow(r, 3)*pow(xj, 3) + 6*pow(r, 4)*pow(xj, 4) +

                     2*pow(r, 5)*pow(xj, 5)) +

                    pow(xi, 12)*(45 + 75*r*xj + 60*pow(r, 2)*pow(xj, 2) +

                                     30*pow(r, 3)*pow(xj, 3) + 10*pow(r, 4)*pow(xj, 4) +

                                     2*pow(r, 5)*pow(xj, 5)) -

                    pow(xi, 10)*pow(xj, 2)*

                    (405 + 675*r*xj + 540*pow(r, 2)*pow(xj, 2) +

                     270*pow(r, 3)*pow(xj, 3) + 90*pow(r, 4)*pow(xj, 4) +

                     8*pow(r, 5)*pow(xj, 5)) +

                    pow(xi, 2)*pow(xj, 10)*

                    (-21615 - 9075*r*xj - 300*pow(r, 2)*pow(xj, 2) +

                     490*pow(r, 3)*pow(xj, 3) + 110*pow(r, 4)*pow(xj, 4) +

                     8*pow(r, 5)*pow(xj, 5))))/

                (45*exp(2*r*(xi + xj))*r*pow(xi - xj, 9)*pow(xi + xj, 8)) -

                (180*exp(2*r*(xi + xj))*(xi + xj)*pow(pow(xi, 2) - pow(xj, 2), 9) +

                 5*exp(2*r*xj)*pow(xj, 8)*

                 (-180*r*pow(xi, 12) - 18*pow(r, 2)*pow(xi, 13) -

                  396*r*pow(xi, 10)*pow(xj, 2) -

                  4*pow(r, 2)*pow(xi, 11)*pow(xj, 2) +

                  1080*r*pow(xi, 8)*pow(xj, 4) +

                  72*pow(r, 2)*pow(xi, 9)*pow(xj, 4) -

                  216*r*pow(xi, 6)*pow(xj, 6) -

                  72*pow(r, 2)*pow(xi, 7)*pow(xj, 6) -

                  324*r*pow(xi, 4)*pow(xj, 8) +

                  4*pow(r, 2)*pow(xi, 5)*pow(xj, 8) + 27*xi*pow(xj, 10) +

                  36*r*pow(xi, 2)*pow(xj, 10) +

                  12*pow(r, 2)*pow(xi, 3)*pow(xj, 10) +

                  2*pow(xi, 5)*pow(xj, 6)*(675 + pow(r, 2)*pow(xj, 2)) -

                  18*pow(xi, 7)*pow(xj, 4)*(-81 + 2*pow(r, 2)*pow(xj, 2)) +

                  3*pow(xi, 3)*pow(xj, 8)*(-81 + 2*pow(r, 2)*pow(xj, 2)) -

                  pow(xi, 11)*(495 + 2*pow(r, 2)*pow(xj, 2)) +

                  9*pow(xi, 9)*pow(xj, 2)*(-233 + 4*pow(r, 2)*pow(xj, 2))) +

                 10*exp(2*r*xj)*pow(xj, 9)*

                 (-90*pow(r, 2)*pow(xi, 12) - 6*pow(r, 3)*pow(xi, 13) +

                  18*pow(xj, 10) + 27*r*xi*pow(xj, 10) +

                  18*pow(xi, 2)*pow(xj, 8)*(-9 + pow(r, 2)*pow(xj, 2)) -

                  162*pow(xi, 4)*pow(xj, 6)*(-4 + pow(r, 2)*pow(xj, 2)) -

                  198*pow(xi, 10)*(5 + pow(r, 2)*pow(xj, 2)) -

                  108*pow(xi, 6)*pow(xj, 4)*(36 + pow(r, 2)*pow(xj, 2)) +

                  2*r*pow(xi, 5)*pow(xj, 6)*(675 + pow(r, 2)*pow(xj, 2)) -

                  18*r*pow(xi, 7)*pow(xj, 4)*(-81 + 2*pow(r, 2)*pow(xj, 2)) +

                  3*r*pow(xi, 3)*pow(xj, 8)*(-81 + 2*pow(r, 2)*pow(xj, 2)) -

                  r*pow(xi, 11)*(495 + 2*pow(r, 2)*pow(xj, 2)) +

                  9*r*pow(xi, 9)*pow(xj, 2)*(-233 + 4*pow(r, 2)*pow(xj, 2)) +

                  6*pow(xi, 8)*pow(xj, 2)*(-1063 + 90*pow(r, 2)*pow(xj, 2))) -

                 2*exp(2*r*xi)*pow(xi, 6)*

                 (-90*pow(xi, 6)*pow(xj, 6)*

                  (65*xj + 152*r*pow(xj, 2) + 66*pow(r, 2)*pow(xj, 3) +

                   8*pow(r, 3)*pow(xj, 4)) -

                  2*pow(xj, 12)*(2475*xj + 1800*r*pow(xj, 2) +

                                       540*pow(r, 2)*pow(xj, 3) + 80*pow(r, 3)*pow(xj, 4) +

                                       5*pow(r, 4)*pow(xj, 5)) +

                  10*pow(xi, 8)*pow(xj, 4)*

                  (270*xj + 432*r*pow(xj, 2) + 366*pow(r, 2)*pow(xj, 3) +

                   88*pow(r, 3)*pow(xj, 4) + 5*pow(r, 4)*pow(xj, 5)) -

                  5*pow(xi, 4)*pow(xj, 8)*

                  (-3555*xj - 2904*r*pow(xj, 2) - 522*pow(r, 2)*pow(xj, 3) +

                   24*pow(r, 3)*pow(xj, 4) + 10*pow(r, 4)*pow(xj, 5)) +

                  pow(xi, 12)*(75*xj + 120*r*pow(xj, 2) +

                                   90*pow(r, 2)*pow(xj, 3) + 40*pow(r, 3)*pow(xj, 4) +

                                   10*pow(r, 4)*pow(xj, 5)) -

                  pow(xi, 10)*pow(xj, 2)*

                  (675*xj + 1080*r*pow(xj, 2) + 810*pow(r, 2)*pow(xj, 3) +

                   360*pow(r, 3)*pow(xj, 4) + 40*pow(r, 4)*pow(xj, 5)) +

                  pow(xi, 2)*pow(xj, 10)*

                  (-9075*xj - 600*r*pow(xj, 2) + 1470*pow(r, 2)*pow(xj, 3) +

                   440*pow(r, 3)*pow(xj, 4) + 40*pow(r, 4)*pow(xj, 5))) -

                 4*exp(2*r*xi)*pow(xi, 7)*

                 (-90*pow(xi, 6)*pow(xj, 6)*

                  (42 + 65*r*xj + 76*pow(r, 2)*pow(xj, 2) +

                   22*pow(r, 3)*pow(xj, 3) + 2*pow(r, 4)*pow(xj, 4)) -

                  2*pow(xj, 12)*(2970 + 2475*r*xj + 900*pow(r, 2)*pow(xj, 2) +

                                       180*pow(r, 3)*pow(xj, 3) + 20*pow(r, 4)*pow(xj, 4) +

                                       pow(r, 5)*pow(xj, 5)) +

                  10*pow(xi, 8)*pow(xj, 4)*

                  (162 + 270*r*xj + 216*pow(r, 2)*pow(xj, 2) +

                   122*pow(r, 3)*pow(xj, 3) + 22*pow(r, 4)*pow(xj, 4) +

                   pow(r, 5)*pow(xj, 5)) -

                  5*pow(xi, 4)*pow(xj, 8)*

                  (-639 - 3555*r*xj - 1452*pow(r, 2)*pow(xj, 2) -

                   174*pow(r, 3)*pow(xj, 3) + 6*pow(r, 4)*pow(xj, 4) +

                   2*pow(r, 5)*pow(xj, 5)) +

                  pow(xi, 12)*(45 + 75*r*xj + 60*pow(r, 2)*pow(xj, 2) +

                                   30*pow(r, 3)*pow(xj, 3) + 10*pow(r, 4)*pow(xj, 4) +

                                   2*pow(r, 5)*pow(xj, 5)) -

                  pow(xi, 10)*pow(xj, 2)*

                  (405 + 675*r*xj + 540*pow(r, 2)*pow(xj, 2) +

                   270*pow(r, 3)*pow(xj, 3) + 90*pow(r, 4)*pow(xj, 4) +

                   8*pow(r, 5)*pow(xj, 5)) +

                  pow(xi, 2)*pow(xj, 10)*

                  (-21615 - 9075*r*xj - 300*pow(r, 2)*pow(xj, 2) +

                   490*pow(r, 3)*pow(xj, 3) + 110*pow(r, 4)*pow(xj, 4) +

                   8*pow(r, 5)*pow(xj, 5))))/

                (90*exp(2*r*(xi + xj))*r*pow(xi - xj, 9)*pow(xi + xj, 9))

            ;
        }

    }
    return S;
}


double DSlater_3S_2S(double r, double xi, double xj)
{
    return DSlater_2S_3S(r, xj, xi);
}

