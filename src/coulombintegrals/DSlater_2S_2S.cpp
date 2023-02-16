/**********************************************************************************
DSlater_2S_2S.cpp - Compute electrostatic force two 2S Slater orbitals

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
double DSlater_2S_2S(double r, double xi, double xj)
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
            S = -(-131985*xi + 161280*exp(2*r*xi)*xi - 205380*r*pow(xi, 2) -

                  149940*pow(r, 2)*pow(xi, 3) - 67200*pow(r, 3)*pow(xi, 4) -

                  20160*pow(r, 4)*pow(xi, 5) - 4032*pow(r, 5)*pow(xi, 6) -

                  448*pow(r, 6)*pow(xi, 7))/(80640*exp(2*r*xi)*r) +

                (-80640 + 80640*exp(2*r*xi) - 131985*r*xi -

                 102690*pow(r, 2)*pow(xi, 2) - 49980*pow(r, 3)*pow(xi, 3) -

                 16800*pow(r, 4)*pow(xi, 4) - 4032*pow(r, 5)*pow(xi, 5) -

                 672*pow(r, 6)*pow(xi, 6) - 64*pow(r, 7)*pow(xi, 7))/

                (80640*exp(2*r*xi)*pow(r, 2)) +

                (xi*(-80640 + 80640*exp(2*r*xi) - 131985*r*xi -

                     102690*pow(r, 2)*pow(xi, 2) - 49980*pow(r, 3)*pow(xi, 3) -

                     16800*pow(r, 4)*pow(xi, 4) - 4032*pow(r, 5)*pow(xi, 5) -

                     672*pow(r, 6)*pow(xi, 6) - 64*pow(r, 7)*pow(xi, 7)))/

                (40320*exp(2*r*xi)*r)

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
            S = (6*exp(2*r*(xi + xj))*pow(pow(xi, 2) - pow(xj, 2), 7) -

                 exp(2*r*xi)*pow(xi, 6)*

                 (21*pow(xi, 4)*pow(xj, 4)*

                  (6 + 11*r*xj + 2*pow(r, 2)*pow(xj, 2)) -

                  2*pow(xj, 8)*(90 + 54*r*xj + 12*pow(r, 2)*pow(xj, 2) +

                                      pow(r, 3)*pow(xj, 3)) +

                  pow(xi, 8)*(6 + 9*r*xj + 6*pow(r, 2)*pow(xj, 2) +

                                  2*pow(r, 3)*pow(xj, 3)) +

                  pow(xi, 2)*pow(xj, 6)*

                  (-390 - 69*r*xj + 18*pow(r, 2)*pow(xj, 2) +

                   4*pow(r, 3)*pow(xj, 3)) -

                  pow(xi, 6)*pow(xj, 2)*

                  (42 + 63*r*xj + 42*pow(r, 2)*pow(xj, 2) +

                   4*pow(r, 3)*pow(xj, 3))) +

                 exp(2*r*xj)*pow(xj, 6)*

                 (-24*pow(r, 2)*pow(xi, 10) - 2*pow(r, 3)*pow(xi, 11) -

                  69*r*pow(xi, 7)*pow(xj, 2) + 6*pow(xj, 8) + 9*r*xi*pow(xj, 8) +

                  4*r*pow(xi, 9)*(-27 + pow(r, 2)*pow(xj, 2)) +

                  18*pow(xi, 8)*(-10 + pow(r, 2)*pow(xj, 2)) +

                  6*pow(xi, 2)*pow(xj, 6)*(-7 + pow(r, 2)*pow(xj, 2)) -

                  42*pow(xi, 4)*pow(xj, 4)*(-3 + pow(r, 2)*pow(xj, 2)) +

                  r*pow(xi, 3)*pow(xj, 6)*(-63 + 2*pow(r, 2)*pow(xj, 2)) +

                  6*pow(xi, 6)*pow(xj, 2)*(-65 + 7*pow(r, 2)*pow(xj, 2)) +

                  pow(xi, 5)*(231*r*pow(xj, 4) - 4*pow(r, 3)*pow(xj, 6))))/

                (6*exp(2*r*(xi + xj))*pow(r, 2)*pow(xi - xj, 7)*pow(xi + xj, 7)) +

                (6*exp(2*r*(xi + xj))*pow(pow(xi, 2) - pow(xj, 2), 7) -

                 exp(2*r*xi)*pow(xi, 6)*

                 (21*pow(xi, 4)*pow(xj, 4)*

                  (6 + 11*r*xj + 2*pow(r, 2)*pow(xj, 2)) -

                  2*pow(xj, 8)*(90 + 54*r*xj + 12*pow(r, 2)*pow(xj, 2) +

                                      pow(r, 3)*pow(xj, 3)) +

                  pow(xi, 8)*(6 + 9*r*xj + 6*pow(r, 2)*pow(xj, 2) +

                                  2*pow(r, 3)*pow(xj, 3)) +

                  pow(xi, 2)*pow(xj, 6)*

                  (-390 - 69*r*xj + 18*pow(r, 2)*pow(xj, 2) +

                   4*pow(r, 3)*pow(xj, 3)) -

                  pow(xi, 6)*pow(xj, 2)*

                  (42 + 63*r*xj + 42*pow(r, 2)*pow(xj, 2) +

                   4*pow(r, 3)*pow(xj, 3))) +

                 exp(2*r*xj)*pow(xj, 6)*

                 (-24*pow(r, 2)*pow(xi, 10) - 2*pow(r, 3)*pow(xi, 11) -

                  69*r*pow(xi, 7)*pow(xj, 2) + 6*pow(xj, 8) + 9*r*xi*pow(xj, 8) +

                  4*r*pow(xi, 9)*(-27 + pow(r, 2)*pow(xj, 2)) +

                  18*pow(xi, 8)*(-10 + pow(r, 2)*pow(xj, 2)) +

                  6*pow(xi, 2)*pow(xj, 6)*(-7 + pow(r, 2)*pow(xj, 2)) -

                  42*pow(xi, 4)*pow(xj, 4)*(-3 + pow(r, 2)*pow(xj, 2)) +

                  r*pow(xi, 3)*pow(xj, 6)*(-63 + 2*pow(r, 2)*pow(xj, 2)) +

                  6*pow(xi, 6)*pow(xj, 2)*(-65 + 7*pow(r, 2)*pow(xj, 2)) +

                  pow(xi, 5)*(231*r*pow(xj, 4) - 4*pow(r, 3)*pow(xj, 6))))/

                (3*exp(2*r*(xi + xj))*r*pow(xi - xj, 7)*pow(xi + xj, 6)) -

                (12*exp(2*r*(xi + xj))*(xi + xj)*pow(pow(xi, 2) - pow(xj, 2), 7) -

                 exp(2*r*xi)*pow(xi, 6)*

                 (21*pow(xi, 4)*pow(xj, 4)*(11*xj + 4*r*pow(xj, 2)) -

                  2*pow(xj, 8)*(54*xj + 24*r*pow(xj, 2) +

                                      3*pow(r, 2)*pow(xj, 3)) +

                  pow(xi, 8)*(9*xj + 12*r*pow(xj, 2) + 6*pow(r, 2)*pow(xj, 3)) +

                  pow(xi, 2)*pow(xj, 6)*

                  (-69*xj + 36*r*pow(xj, 2) + 12*pow(r, 2)*pow(xj, 3)) -

                  pow(xi, 6)*pow(xj, 2)*

                  (63*xj + 84*r*pow(xj, 2) + 12*pow(r, 2)*pow(xj, 3))) -

                 2*exp(2*r*xi)*pow(xi, 7)*

                 (21*pow(xi, 4)*pow(xj, 4)*

                  (6 + 11*r*xj + 2*pow(r, 2)*pow(xj, 2)) -

                  2*pow(xj, 8)*(90 + 54*r*xj + 12*pow(r, 2)*pow(xj, 2) +

                                      pow(r, 3)*pow(xj, 3)) +

                  pow(xi, 8)*(6 + 9*r*xj + 6*pow(r, 2)*pow(xj, 2) +

                                  2*pow(r, 3)*pow(xj, 3)) +

                  pow(xi, 2)*pow(xj, 6)*

                  (-390 - 69*r*xj + 18*pow(r, 2)*pow(xj, 2) +

                   4*pow(r, 3)*pow(xj, 3)) -

                  pow(xi, 6)*pow(xj, 2)*

                  (42 + 63*r*xj + 42*pow(r, 2)*pow(xj, 2) +

                   4*pow(r, 3)*pow(xj, 3))) +

                 exp(2*r*xj)*pow(xj, 6)*

                 (-48*r*pow(xi, 10) - 6*pow(r, 2)*pow(xi, 11) -

                  69*pow(xi, 7)*pow(xj, 2) + 36*r*pow(xi, 8)*pow(xj, 2) +

                  8*pow(r, 2)*pow(xi, 9)*pow(xj, 2) +

                  84*r*pow(xi, 6)*pow(xj, 4) - 84*r*pow(xi, 4)*pow(xj, 6) +

                  9*xi*pow(xj, 8) + 12*r*pow(xi, 2)*pow(xj, 8) +

                  4*pow(r, 2)*pow(xi, 3)*pow(xj, 8) +

                  4*pow(xi, 9)*(-27 + pow(r, 2)*pow(xj, 2)) +

                  pow(xi, 3)*pow(xj, 6)*(-63 + 2*pow(r, 2)*pow(xj, 2)) +

                  pow(xi, 5)*(231*pow(xj, 4) - 12*pow(r, 2)*pow(xj, 6))) +

                 2*exp(2*r*xj)*pow(xj, 7)*

                 (-24*pow(r, 2)*pow(xi, 10) - 2*pow(r, 3)*pow(xi, 11) -

                  69*r*pow(xi, 7)*pow(xj, 2) + 6*pow(xj, 8) + 9*r*xi*pow(xj, 8) +

                  4*r*pow(xi, 9)*(-27 + pow(r, 2)*pow(xj, 2)) +

                  18*pow(xi, 8)*(-10 + pow(r, 2)*pow(xj, 2)) +

                  6*pow(xi, 2)*pow(xj, 6)*(-7 + pow(r, 2)*pow(xj, 2)) -

                  42*pow(xi, 4)*pow(xj, 4)*(-3 + pow(r, 2)*pow(xj, 2)) +

                  r*pow(xi, 3)*pow(xj, 6)*(-63 + 2*pow(r, 2)*pow(xj, 2)) +

                  6*pow(xi, 6)*pow(xj, 2)*(-65 + 7*pow(r, 2)*pow(xj, 2)) +

                  pow(xi, 5)*(231*r*pow(xj, 4) - 4*pow(r, 3)*pow(xj, 6))))/

                (6*exp(2*r*(xi + xj))*r*pow(xi - xj, 7)*pow(xi + xj, 7))

            ;
        }

    }
    return S;
}

} // end namespace openbabel
