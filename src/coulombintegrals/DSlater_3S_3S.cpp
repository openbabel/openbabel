/**********************************************************************************
DSlater_3S_3S.cpp - Compute electrostatic force two 3S Slater orbitals

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

#pragma warning(disable : 4146) // warning C4146: unary minus operator applied to unsigned type, result still unsigned

namespace OpenBabel
{
double DSlater_3S_3S(double r, double xi, double xj)
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
            S = -(-2503064025*xi + 2874009600*exp(2*r*xi)*xi -

                  4264236900*r*pow(xi, 2) - 3541992300*pow(r, 2)*pow(xi, 3) -

                  1906027200*pow(r, 3)*pow(xi, 4) -

                  744282000*pow(r, 4)*pow(xi, 5) - 223534080*pow(r, 5)*pow(xi, 6) -

                  53222400*pow(r, 6)*pow(xi, 7) - 10137600*pow(r, 7)*pow(xi, 8) -

                  1520640*pow(r, 8)*pow(xi, 9) - 168960*pow(r, 9)*pow(xi, 10) -

                  11264*pow(r, 10)*pow(xi, 11))/(1.4370048e9*exp(2*r*xi)*r) +

                (-1437004800 + 1437004800*exp(2*r*xi) - 2503064025*r*xi -

                 2132118450*pow(r, 2)*pow(xi, 2) -

                 1180664100*pow(r, 3)*pow(xi, 3) - 476506800*pow(r, 4)*pow(xi, 4) -

                 148856400*pow(r, 5)*pow(xi, 5) - 37255680*pow(r, 6)*pow(xi, 6) -

                 7603200*pow(r, 7)*pow(xi, 7) - 1267200*pow(r, 8)*pow(xi, 8) -

                 168960*pow(r, 9)*pow(xi, 9) - 16896*pow(r, 10)*pow(xi, 10) -

                 1024*pow(r, 11)*pow(xi, 11))/(1.4370048e9*exp(2*r*xi)*pow(r, 2))

                + (xi*(-1437004800 + 1437004800*exp(2*r*xi) - 2503064025*r*xi -

                       2132118450*pow(r, 2)*pow(xi, 2) -

                       1180664100*pow(r, 3)*pow(xi, 3) -

                       476506800*pow(r, 4)*pow(xi, 4) - 148856400*pow(r, 5)*pow(xi, 5) -

                       37255680*pow(r, 6)*pow(xi, 6) - 7603200*pow(r, 7)*pow(xi, 7) -

                       1267200*pow(r, 8)*pow(xi, 8) - 168960*pow(r, 9)*pow(xi, 9) -

                       16896*pow(r, 10)*pow(xi, 10) - 1024*pow(r, 11)*pow(xi, 11)))/

                (7.185024e8*exp(2*r*xi)*r)

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
            S = (135*exp(2*r*(xi + xj))*pow(pow(xi, 2) - pow(xj, 2), 11) +

                 exp(2*r*xj)*pow(xj, 8)*

                 (-150*pow(r, 4)*pow(xi, 18) - 6*pow(r, 5)*pow(xi, 19) +

                  135*pow(xj, 14) + 225*r*xi*pow(xj, 14) +

                  10*pow(r, 3)*pow(xi, 17)*(-165 + pow(r, 2)*pow(xj, 2)) -

                  30*pow(r, 2)*pow(xi, 16)*(330 + pow(r, 2)*pow(xj, 2)) +

                  45*r*pow(xi, 3)*pow(xj, 12)*(-55 + 2*pow(r, 2)*pow(xj, 2)) +

                  45*pow(xi, 2)*pow(xj, 12)*(-33 + 4*pow(r, 2)*pow(xj, 2)) +

                  r*pow(xi, 9)*pow(xj, 6)*

                  (234135 - 4950*pow(r, 2)*pow(xj, 2) -

                   34*pow(r, 4)*pow(xj, 4)) -

                  5*r*pow(xi, 7)*pow(xj, 8)*

                  (6237 - 1242*pow(r, 2)*pow(xj, 2) + 2*pow(r, 4)*pow(xj, 4)) +

                  3*r*pow(xi, 5)*pow(xj, 10)*

                  (4125 - 330*pow(r, 2)*pow(xj, 2) + 2*pow(r, 4)*pow(xj, 4)) +

                  15*pow(xi, 4)*pow(xj, 10)*

                  (495 - 132*pow(r, 2)*pow(xj, 2) + 2*pow(r, 4)*pow(xj, 4)) -

                  165*pow(xi, 6)*pow(xj, 8)*

                  (135 - 60*pow(r, 2)*pow(xj, 2) + 2*pow(r, 4)*pow(xj, 4)) -

                  5*r*pow(xi, 13)*pow(xj, 2)*

                  (43875 - 3438*pow(r, 2)*pow(xj, 2) + 22*pow(r, 4)*pow(xj, 4))

                  + 5*r*pow(xi, 11)*pow(xj, 4)*

                  (7695 - 2442*pow(r, 2)*pow(xj, 2) + 22*pow(r, 4)*pow(xj, 4))

                  + 15*pow(xi, 8)*pow(xj, 6)*(-33 - 3564*pow(r, 2)*pow(xj, 2) +

                                                        26*pow(r, 4)*pow(xj, 4)) +

                  r*pow(xi, 15)*(-32175 - 3690*pow(r, 2)*pow(xj, 2) +

                                     34*pow(r, 4)*pow(xj, 4)) +

                  15*pow(xi, 10)*pow(xj, 4)*

                  (-32277 + 1364*pow(r, 2)*pow(xj, 2) +

                   66*pow(r, 4)*pow(xj, 4)) +

                  15*pow(xi, 14)*(-3003 - 2932*pow(r, 2)*pow(xj, 2) +

                                        94*pow(r, 4)*pow(xj, 4)) -

                  15*pow(xi, 12)*pow(xj, 2)*

                  (28119 - 5252*pow(r, 2)*pow(xj, 2) + 154*pow(r, 4)*pow(xj, 4))

                 ) + exp(2*r*xi)*pow(xi, 8)*

                 (-5*pow(xi, 2)*pow(xj, 12)*

                  (-84357 - 43875*r*xj - 8796*pow(r, 2)*pow(xj, 2) -

                   738*pow(r, 3)*pow(xj, 3) - 6*pow(r, 4)*pow(xj, 4) +

                   2*pow(r, 5)*pow(xj, 5)) -

                  3*pow(xi, 14)*(45 + 75*r*xj + 60*pow(r, 2)*pow(xj, 2) +

                                       30*pow(r, 3)*pow(xj, 3) + 10*pow(r, 4)*pow(xj, 4) +

                                       2*pow(r, 5)*pow(xj, 5)) -

                  55*pow(xi, 8)*pow(xj, 6)*

                  (-405 - 567*r*xj - 972*pow(r, 2)*pow(xj, 2) -

                   90*pow(r, 3)*pow(xj, 3) + 18*pow(r, 4)*pow(xj, 4) +

                   2*pow(r, 5)*pow(xj, 5)) +

                  55*pow(xi, 6)*pow(xj, 8)*

                  (9 - 4257*r*xj - 372*pow(r, 2)*pow(xj, 2) +

                   222*pow(r, 3)*pow(xj, 3) + 42*pow(r, 4)*pow(xj, 4) +

                   2*pow(r, 5)*pow(xj, 5)) +

                  3*pow(xj, 14)*(15015 + 10725*r*xj + 3300*pow(r, 2)*pow(xj, 2) +

                                       550*pow(r, 3)*pow(xj, 3) + 50*pow(r, 4)*pow(xj, 4) +

                                       2*pow(r, 5)*pow(xj, 5)) +

                  5*pow(xi, 12)*pow(xj, 2)*

                  (297 + 495*r*xj + 396*pow(r, 2)*pow(xj, 2) +

                   198*pow(r, 3)*pow(xj, 3) + 66*pow(r, 4)*pow(xj, 4) +

                   2*pow(r, 5)*pow(xj, 5)) +

                  pow(xi, 10)*pow(xj, 4)*

                  (-7425 - 12375*r*xj - 9900*pow(r, 2)*pow(xj, 2) -

                   6210*pow(r, 3)*pow(xj, 3) - 390*pow(r, 4)*pow(xj, 4) +

                   34*pow(r, 5)*pow(xj, 5)) -

                  pow(xi, 4)*pow(xj, 10)*

                  (-484155 + 38475*r*xj + 78780*pow(r, 2)*pow(xj, 2) +

                   17190*pow(r, 3)*pow(xj, 3) + 1410*pow(r, 4)*pow(xj, 4) +

                   34*pow(r, 5)*pow(xj, 5))))/

                (135*exp(2*r*(xi + xj))*pow(r, 2)*pow(xi - xj, 11)*

                 pow(xi + xj, 11)) + (2*(135*exp(2*r*(xi + xj))*

                                               pow(pow(xi, 2) - pow(xj, 2), 11) +

                                               exp(2*r*xj)*pow(xj, 8)*

                                               (-150*pow(r, 4)*pow(xi, 18) - 6*pow(r, 5)*pow(xi, 19) +

                                                135*pow(xj, 14) + 225*r*xi*pow(xj, 14) +

                                                10*pow(r, 3)*pow(xi, 17)*(-165 + pow(r, 2)*pow(xj, 2)) -

                                                30*pow(r, 2)*pow(xi, 16)*(330 + pow(r, 2)*pow(xj, 2)) +

                                                45*r*pow(xi, 3)*pow(xj, 12)*(-55 + 2*pow(r, 2)*pow(xj, 2)) +

                                                45*pow(xi, 2)*pow(xj, 12)*(-33 + 4*pow(r, 2)*pow(xj, 2)) +

                                                r*pow(xi, 9)*pow(xj, 6)*

                                                (234135 - 4950*pow(r, 2)*pow(xj, 2) -

                                                 34*pow(r, 4)*pow(xj, 4)) -

                                                5*r*pow(xi, 7)*pow(xj, 8)*

                                                (6237 - 1242*pow(r, 2)*pow(xj, 2) + 2*pow(r, 4)*pow(xj, 4))

                                                + 3*r*pow(xi, 5)*pow(xj, 10)*

                                                (4125 - 330*pow(r, 2)*pow(xj, 2) + 2*pow(r, 4)*pow(xj, 4))

                                                + 15*pow(xi, 4)*pow(xj, 10)*

                                                (495 - 132*pow(r, 2)*pow(xj, 2) + 2*pow(r, 4)*pow(xj, 4)) -

                                                165*pow(xi, 6)*pow(xj, 8)*

                                                (135 - 60*pow(r, 2)*pow(xj, 2) + 2*pow(r, 4)*pow(xj, 4)) -

                                                5*r*pow(xi, 13)*pow(xj, 2)*

                                                (43875 - 3438*pow(r, 2)*pow(xj, 2) +

                                                 22*pow(r, 4)*pow(xj, 4)) +

                                                5*r*pow(xi, 11)*pow(xj, 4)*

                                                (7695 - 2442*pow(r, 2)*pow(xj, 2) +

                                                 22*pow(r, 4)*pow(xj, 4)) +

                                                15*pow(xi, 8)*pow(xj, 6)*

                                                (-33 - 3564*pow(r, 2)*pow(xj, 2) + 26*pow(r, 4)*pow(xj, 4))

                                                + r*pow(xi, 15)*(-32175 - 3690*pow(r, 2)*pow(xj, 2) +

                                                                     34*pow(r, 4)*pow(xj, 4)) +

                                                15*pow(xi, 10)*pow(xj, 4)*

                                                (-32277 + 1364*pow(r, 2)*pow(xj, 2) +

                                                 66*pow(r, 4)*pow(xj, 4)) +

                                                15*pow(xi, 14)*(-3003 - 2932*pow(r, 2)*pow(xj, 2) +

                                                                      94*pow(r, 4)*pow(xj, 4)) -

                                                15*pow(xi, 12)*pow(xj, 2)*

                                                (28119 - 5252*pow(r, 2)*pow(xj, 2) +

                                                 154*pow(r, 4)*pow(xj, 4))) +

                                               exp(2*r*xi)*pow(xi, 8)*

                                               (-5*pow(xi, 2)*pow(xj, 12)*

                                                (-84357 - 43875*r*xj - 8796*pow(r, 2)*pow(xj, 2) -

                                                 738*pow(r, 3)*pow(xj, 3) - 6*pow(r, 4)*pow(xj, 4) +

                                                 2*pow(r, 5)*pow(xj, 5)) -

                                                3*pow(xi, 14)*(45 + 75*r*xj + 60*pow(r, 2)*pow(xj, 2) +

                                                                     30*pow(r, 3)*pow(xj, 3) + 10*pow(r, 4)*pow(xj, 4) +

                                                                     2*pow(r, 5)*pow(xj, 5)) -

                                                55*pow(xi, 8)*pow(xj, 6)*

                                                (-405 - 567*r*xj - 972*pow(r, 2)*pow(xj, 2) -

                                                 90*pow(r, 3)*pow(xj, 3) + 18*pow(r, 4)*pow(xj, 4) +

                                                 2*pow(r, 5)*pow(xj, 5)) +

                                                55*pow(xi, 6)*pow(xj, 8)*

                                                (9 - 4257*r*xj - 372*pow(r, 2)*pow(xj, 2) +

                                                 222*pow(r, 3)*pow(xj, 3) + 42*pow(r, 4)*pow(xj, 4) +

                                                 2*pow(r, 5)*pow(xj, 5)) +

                                                3*pow(xj, 14)*(15015 + 10725*r*xj +

                                                                     3300*pow(r, 2)*pow(xj, 2) + 550*pow(r, 3)*pow(xj, 3) +

                                                                     50*pow(r, 4)*pow(xj, 4) + 2*pow(r, 5)*pow(xj, 5)) +

                                                5*pow(xi, 12)*pow(xj, 2)*

                                                (297 + 495*r*xj + 396*pow(r, 2)*pow(xj, 2) +

                                                 198*pow(r, 3)*pow(xj, 3) + 66*pow(r, 4)*pow(xj, 4) +

                                                 2*pow(r, 5)*pow(xj, 5)) +

                                                pow(xi, 10)*pow(xj, 4)*

                                                (-7425 - 12375*r*xj - 9900*pow(r, 2)*pow(xj, 2) -

                                                 6210*pow(r, 3)*pow(xj, 3) - 390*pow(r, 4)*pow(xj, 4) +

                                                 34*pow(r, 5)*pow(xj, 5)) -

                                                pow(xi, 4)*pow(xj, 10)*

                                                (-484155 + 38475*r*xj + 78780*pow(r, 2)*pow(xj, 2) +

                                                 17190*pow(r, 3)*pow(xj, 3) + 1410*pow(r, 4)*pow(xj, 4) +

                                                 34*pow(r, 5)*pow(xj, 5)))))/

                (135*exp(2*r*(xi + xj))*r*pow(xi - xj, 11)*pow(xi + xj, 10)) -

                (270*exp(2*r*(xi + xj))*(xi + xj)*

                 pow(pow(xi, 2) - pow(xj, 2), 11) +

                 exp(2*r*xj)*pow(xj, 8)*

                 (-600*pow(r, 3)*pow(xi, 18) - 30*pow(r, 4)*pow(xi, 19) -

                  60*pow(r, 3)*pow(xi, 16)*pow(xj, 2) +

                  20*pow(r, 4)*pow(xi, 17)*pow(xj, 2) + 225*xi*pow(xj, 14) +

                  360*r*pow(xi, 2)*pow(xj, 14) +

                  180*pow(r, 2)*pow(xi, 3)*pow(xj, 14) +

                  30*pow(r, 2)*pow(xi, 17)*(-165 + pow(r, 2)*pow(xj, 2)) -

                  60*r*pow(xi, 16)*(330 + pow(r, 2)*pow(xj, 2)) +

                  45*pow(xi, 3)*pow(xj, 12)*(-55 + 2*pow(r, 2)*pow(xj, 2)) +

                  r*pow(xi, 9)*pow(xj, 6)*

                  (-9900*r*pow(xj, 2) - 136*pow(r, 3)*pow(xj, 4)) -

                  5*r*pow(xi, 7)*pow(xj, 8)*

                  (-2484*r*pow(xj, 2) + 8*pow(r, 3)*pow(xj, 4)) +

                  3*r*pow(xi, 5)*pow(xj, 10)*

                  (-660*r*pow(xj, 2) + 8*pow(r, 3)*pow(xj, 4)) +

                  15*pow(xi, 4)*pow(xj, 10)*

                  (-264*r*pow(xj, 2) + 8*pow(r, 3)*pow(xj, 4)) -

                  165*pow(xi, 6)*pow(xj, 8)*

                  (-120*r*pow(xj, 2) + 8*pow(r, 3)*pow(xj, 4)) -

                  5*r*pow(xi, 13)*pow(xj, 2)*

                  (-6876*r*pow(xj, 2) + 88*pow(r, 3)*pow(xj, 4)) +

                  5*r*pow(xi, 11)*pow(xj, 4)*

                  (-4884*r*pow(xj, 2) + 88*pow(r, 3)*pow(xj, 4)) +

                  15*pow(xi, 8)*pow(xj, 6)*

                  (-7128*r*pow(xj, 2) + 104*pow(r, 3)*pow(xj, 4)) +

                  r*pow(xi, 15)*(-7380*r*pow(xj, 2) + 136*pow(r, 3)*pow(xj, 4)) +

                  15*pow(xi, 10)*pow(xj, 4)*

                  (2728*r*pow(xj, 2) + 264*pow(r, 3)*pow(xj, 4)) +

                  15*pow(xi, 14)*(-5864*r*pow(xj, 2) +

                                        376*pow(r, 3)*pow(xj, 4)) -

                  15*pow(xi, 12)*pow(xj, 2)*

                  (-10504*r*pow(xj, 2) + 616*pow(r, 3)*pow(xj, 4)) +

                  pow(xi, 9)*pow(xj, 6)*

                  (234135 - 4950*pow(r, 2)*pow(xj, 2) - 34*pow(r, 4)*pow(xj, 4))

                  - 5*pow(xi, 7)*pow(xj, 8)*(6237 - 1242*pow(r, 2)*pow(xj, 2) +

                                                       2*pow(r, 4)*pow(xj, 4)) +

                  3*pow(xi, 5)*pow(xj, 10)*

                  (4125 - 330*pow(r, 2)*pow(xj, 2) + 2*pow(r, 4)*pow(xj, 4)) -

                  5*pow(xi, 13)*pow(xj, 2)*

                  (43875 - 3438*pow(r, 2)*pow(xj, 2) + 22*pow(r, 4)*pow(xj, 4))

                  + 5*pow(xi, 11)*pow(xj, 4)*(7695 - 2442*pow(r, 2)*pow(xj, 2) +

                                                        22*pow(r, 4)*pow(xj, 4)) +

                  pow(xi, 15)*(-32175 - 3690*pow(r, 2)*pow(xj, 2) +

                                   34*pow(r, 4)*pow(xj, 4))) +

                 2*exp(2*r*xj)*pow(xj, 9)*

                 (-150*pow(r, 4)*pow(xi, 18) - 6*pow(r, 5)*pow(xi, 19) +

                  135*pow(xj, 14) + 225*r*xi*pow(xj, 14) +

                  10*pow(r, 3)*pow(xi, 17)*(-165 + pow(r, 2)*pow(xj, 2)) -

                  30*pow(r, 2)*pow(xi, 16)*(330 + pow(r, 2)*pow(xj, 2)) +

                  45*r*pow(xi, 3)*pow(xj, 12)*(-55 + 2*pow(r, 2)*pow(xj, 2)) +

                  45*pow(xi, 2)*pow(xj, 12)*(-33 + 4*pow(r, 2)*pow(xj, 2)) +

                  r*pow(xi, 9)*pow(xj, 6)*

                  (234135 - 4950*pow(r, 2)*pow(xj, 2) - 34*pow(r, 4)*pow(xj, 4))

                  - 5*r*pow(xi, 7)*pow(xj, 8)*(6237 - 1242*pow(r, 2)*pow(xj, 2) +

                                                         2*pow(r, 4)*pow(xj, 4)) +

                  3*r*pow(xi, 5)*pow(xj, 10)*

                  (4125 - 330*pow(r, 2)*pow(xj, 2) + 2*pow(r, 4)*pow(xj, 4)) +

                  15*pow(xi, 4)*pow(xj, 10)*

                  (495 - 132*pow(r, 2)*pow(xj, 2) + 2*pow(r, 4)*pow(xj, 4)) -

                  165*pow(xi, 6)*pow(xj, 8)*

                  (135 - 60*pow(r, 2)*pow(xj, 2) + 2*pow(r, 4)*pow(xj, 4)) -

                  5*r*pow(xi, 13)*pow(xj, 2)*

                  (43875 - 3438*pow(r, 2)*pow(xj, 2) + 22*pow(r, 4)*pow(xj, 4))

                  + 5*r*pow(xi, 11)*pow(xj, 4)*

                  (7695 - 2442*pow(r, 2)*pow(xj, 2) + 22*pow(r, 4)*pow(xj, 4)) +

                  15*pow(xi, 8)*pow(xj, 6)*

                  (-33 - 3564*pow(r, 2)*pow(xj, 2) + 26*pow(r, 4)*pow(xj, 4)) +

                  r*pow(xi, 15)*(-32175 - 3690*pow(r, 2)*pow(xj, 2) +

                                     34*pow(r, 4)*pow(xj, 4)) +

                  15*pow(xi, 10)*pow(xj, 4)*

                  (-32277 + 1364*pow(r, 2)*pow(xj, 2) + 66*pow(r, 4)*pow(xj, 4))

                  + 15*pow(xi, 14)*(-3003 - 2932*pow(r, 2)*pow(xj, 2) +

                                          94*pow(r, 4)*pow(xj, 4)) -

                  15*pow(xi, 12)*pow(xj, 2)*

                  (28119 - 5252*pow(r, 2)*pow(xj, 2) + 154*pow(r, 4)*pow(xj, 4)))

                 + exp(2*r*xi)*pow(xi, 8)*(-5*pow(xi, 2)*pow(xj, 12)*

                                                 (-43875*xj - 17592*r*pow(xj, 2) - 2214*pow(r, 2)*pow(xj, 3) -

                                                  24*pow(r, 3)*pow(xj, 4) + 10*pow(r, 4)*pow(xj, 5)) -

                                                 3*pow(xi, 14)*(75*xj + 120*r*pow(xj, 2) +

                                                                      90*pow(r, 2)*pow(xj, 3) + 40*pow(r, 3)*pow(xj, 4) +

                                                                      10*pow(r, 4)*pow(xj, 5)) -

                                                 55*pow(xi, 8)*pow(xj, 6)*

                                                 (-567*xj - 1944*r*pow(xj, 2) - 270*pow(r, 2)*pow(xj, 3) +

                                                  72*pow(r, 3)*pow(xj, 4) + 10*pow(r, 4)*pow(xj, 5)) +

                                                 55*pow(xi, 6)*pow(xj, 8)*

                                                 (-4257*xj - 744*r*pow(xj, 2) + 666*pow(r, 2)*pow(xj, 3) +

                                                  168*pow(r, 3)*pow(xj, 4) + 10*pow(r, 4)*pow(xj, 5)) +

                                                 3*pow(xj, 14)*(10725*xj + 6600*r*pow(xj, 2) +

                                                                      1650*pow(r, 2)*pow(xj, 3) + 200*pow(r, 3)*pow(xj, 4) +

                                                                      10*pow(r, 4)*pow(xj, 5)) +

                                                 5*pow(xi, 12)*pow(xj, 2)*

                                                 (495*xj + 792*r*pow(xj, 2) + 594*pow(r, 2)*pow(xj, 3) +

                                                  264*pow(r, 3)*pow(xj, 4) + 10*pow(r, 4)*pow(xj, 5)) +

                                                 pow(xi, 10)*pow(xj, 4)*

                                                 (-12375*xj - 19800*r*pow(xj, 2) - 18630*pow(r, 2)*pow(xj, 3) -

                                                  1560*pow(r, 3)*pow(xj, 4) + 170*pow(r, 4)*pow(xj, 5)) -

                                                 pow(xi, 4)*pow(xj, 10)*

                                                 (38475*xj + 157560*r*pow(xj, 2) + 51570*pow(r, 2)*pow(xj, 3) +

                                                  5640*pow(r, 3)*pow(xj, 4) + 170*pow(r, 4)*pow(xj, 5))) +

                 2*exp(2*r*xi)*pow(xi, 9)*

                 (-5*pow(xi, 2)*pow(xj, 12)*

                  (-84357 - 43875*r*xj - 8796*pow(r, 2)*pow(xj, 2) -

                   738*pow(r, 3)*pow(xj, 3) - 6*pow(r, 4)*pow(xj, 4) +

                   2*pow(r, 5)*pow(xj, 5)) -

                  3*pow(xi, 14)*(45 + 75*r*xj + 60*pow(r, 2)*pow(xj, 2) +

                                       30*pow(r, 3)*pow(xj, 3) + 10*pow(r, 4)*pow(xj, 4) +

                                       2*pow(r, 5)*pow(xj, 5)) -

                  55*pow(xi, 8)*pow(xj, 6)*

                  (-405 - 567*r*xj - 972*pow(r, 2)*pow(xj, 2) -

                   90*pow(r, 3)*pow(xj, 3) + 18*pow(r, 4)*pow(xj, 4) +

                   2*pow(r, 5)*pow(xj, 5)) +

                  55*pow(xi, 6)*pow(xj, 8)*

                  (9 - 4257*r*xj - 372*pow(r, 2)*pow(xj, 2) +

                   222*pow(r, 3)*pow(xj, 3) + 42*pow(r, 4)*pow(xj, 4) +

                   2*pow(r, 5)*pow(xj, 5)) +

                  3*pow(xj, 14)*(15015 + 10725*r*xj + 3300*pow(r, 2)*pow(xj, 2) +

                                       550*pow(r, 3)*pow(xj, 3) + 50*pow(r, 4)*pow(xj, 4) +

                                       2*pow(r, 5)*pow(xj, 5)) +

                  5*pow(xi, 12)*pow(xj, 2)*

                  (297 + 495*r*xj + 396*pow(r, 2)*pow(xj, 2) +

                   198*pow(r, 3)*pow(xj, 3) + 66*pow(r, 4)*pow(xj, 4) +

                   2*pow(r, 5)*pow(xj, 5)) +

                  pow(xi, 10)*pow(xj, 4)*

                  (-7425 - 12375*r*xj - 9900*pow(r, 2)*pow(xj, 2) -

                   6210*pow(r, 3)*pow(xj, 3) - 390*pow(r, 4)*pow(xj, 4) +

                   34*pow(r, 5)*pow(xj, 5)) -

                  pow(xi, 4)*pow(xj, 10)*

                  (-484155 + 38475*r*xj + 78780*pow(r, 2)*pow(xj, 2) +

                   17190*pow(r, 3)*pow(xj, 3) + 1410*pow(r, 4)*pow(xj, 4) +

                   34*pow(r, 5)*pow(xj, 5))))/

                (135*exp(2*r*(xi + xj))*r*pow(xi - xj, 11)*pow(xi + xj, 11))

            ;
        }

    }
    return S;
}

} // end namescpace openbabel
