/**********************************************************************************
gaussian_integrals.cpp - coulombintegrals for Slater s-type orbitals

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


#include <stdio.h>
#include <cmath>
#include <iostream>

#include <openbabel/slater_integrals.h>

using namespace std;

#define SLATER_MAX 3

namespace OpenBabel
{
typedef double t_slater_SS_func (double r, double xi, double xj);
typedef double t_slater_NS_func (double r, double xi);

t_slater_SS_func (*Slater_SS[SLATER_MAX][SLATER_MAX]) = {
    {  Slater_1S_1S,  Slater_1S_2S,  Slater_1S_3S},
    {  Slater_2S_1S,  Slater_2S_2S,  Slater_2S_3S},
    {  Slater_3S_1S,  Slater_3S_2S,  Slater_3S_3S}
};

t_slater_SS_func (*DSlater_SS[SLATER_MAX][SLATER_MAX]) = {
    {  DSlater_1S_1S,  DSlater_1S_2S,  DSlater_1S_3S},
    {  DSlater_2S_1S,  DSlater_2S_2S,  DSlater_2S_3S},
    {  DSlater_3S_1S,  DSlater_3S_2S,  DSlater_3S_3S}
};

t_slater_NS_func (*Slater_NS[SLATER_MAX])  = {Nuclear_1S,  Nuclear_2S,  Nuclear_3S};

t_slater_NS_func (*DSlater_NS[SLATER_MAX]) = {DNuclear_1S,  DNuclear_2S,  DNuclear_3S};


double Nuclear_1S(double r, double xi)
{
    double S = 0;
    S = 1/r - (1 + r*xi)/(exp(2*r*xi)*r)

    ;

    return S;
}

double Nuclear_2S(double r, double xi)
{
    double S = 0;
    S = 1/r - (6 + 9*r*xi + 6*pow(r, 2)*pow(xi, 2) + 2*pow(r, 3)*pow(xi, 3))/

        (6*exp(2*r*xi)*r)

    ;

    return S;
}

double Nuclear_3S(double r, double xi)
{
    double S = 0;
    S = 1/r - (45 + 75*r*xi + 60*pow(r, 2)*pow(xi, 2) +

                 30*pow(r, 3)*pow(xi, 3) + 10*pow(r, 4)*pow(xi, 4) +

                 2*pow(r, 5)*pow(xi, 5))/(45*exp(2*r*xi)*r)

    ;

    return S;
}

double DNuclear_1S(double r, double xi)
{
    double S = 0;
    S = pow(r, -2) - (1 + 2*r*xi + 2*pow(r, 2)*pow(xi, 2))/

        (exp(2*r*xi)*pow(r, 2))

    ;

    return S;
}

double DNuclear_2S(double r, double xi)
{
    double S = 0;
    S = pow(r, -2) - (3 + 6*r*xi + 6*pow(r, 2)*pow(xi, 2) +

                          4*pow(r, 3)*pow(xi, 3) + 2*pow(r, 4)*pow(xi, 4))/

        (3*exp(2*r*xi)*pow(r, 2))

    ;

    return S;
}

double DNuclear_3S(double r, double xi)
{
    double S = 0;
    S = pow(r, -2) - (45 + 90*r*xi + 90*pow(r, 2)*pow(xi, 2) +

                          60*pow(r, 3)*pow(xi, 3) + 30*pow(r, 4)*pow(xi, 4) +

                          12*pow(r, 5)*pow(xi, 5) + 4*pow(r, 6)*pow(xi, 6))/

        (45*exp(2*r*xi)*pow(r, 2))

    ;

    return S;
}

double Coulomb_SS(double r, int i, int j, double xi, double xj)
{   
    double S = 0;
    if ((i > SLATER_MAX) || (j > SLATER_MAX))
    {
        cout << "Slater-Slater integral" << i << j << "not supported. Thus, they will reduce to the SLATER_MAX" << SLATER_MAX << endl;
    }    
    if (i > SLATER_MAX)
    {
        i = SLATER_MAX;
    }
    if (j > SLATER_MAX)
    {
        j = SLATER_MAX;
    }    
    if ((i > 0) && (j > 0))
    {
        S = Slater_SS[i-1][j-1](r, xi, xj);
    }
    else if (j > 0)
    {
        if (r == 0)
        {
            S = xj/j;
        }
        else
        {
            S = Slater_NS[j-1](r, xj);
        }
    }
    else if (i > 0)
    {
        if (r == 0)
        {
            S = xi/i;
        }
        else
        {
            S = Slater_NS[i-1](r, xi);
        }
    }
    else
    {
        if (r == 0)
        {
            return 0;
        }
        else
        {
            return 1/r;
        }
    }
    return S;
}

double Nuclear_SS(double r, int i, double xi)
{
    double S = 0;
    if (xi == 0)
    {
        if (r == 0)
        {
            return 0;
        }
        else
        {
            return 1/r;
        }
    }
    else if (r == 0)
    {
        return xi/i;
    }
    else if (i <= 0)
    {
        return 1/r;
    }
    else
    {
        S  = Slater_NS[i-1](r, xi);
        return S;
    }
}

double DCoulomb_SS(double r, int i, int j, double xi, double xj)
{
    double S = 0;    
    if ((i > SLATER_MAX) || (j > SLATER_MAX))
    {
        cout << "Slater-Slater integral" << i << j << "not supported. Thus, they will reduce to the SLATER_MAX" << SLATER_MAX << endl;
    }    
    if (i > SLATER_MAX)
    {
        i = SLATER_MAX;
    }
    if (j > SLATER_MAX)
    {
        j = SLATER_MAX;
    }
    if ((i > 0) && (j > 0))
    {
        S = DSlater_SS[i-1][j-1](r, xi, xj);
    }
    else if (j > 0)
    {
        if (r == 0)
        {
            S = 0;
        }
        else
        {
            S = DSlater_NS[j-1](r, xj);
        }
    }
    else if (i > 0)
    {
        if (r == 0)
        {
            S = 0;
        }
        else
        {
            S = DSlater_NS[i-1](r, xi);
        }
    }
    else
    {
        if (r == 0)
        {
            return 0;
        }
        else
        {
            return -1/(r*r);
        }
    }
    return S;
}

double DNuclear_SS(double r, int i, double xi)
{
    double S = 0;
    if (i > SLATER_MAX)
    {
        cout << "Slater-Slater integral" << i << "not supported. Thus, they will reduce to the SLATER_MAX" << SLATER_MAX << endl;
        i = SLATER_MAX;
    }    
    if (r == 0)
    {
        return 0;
    }
    else
    {
        if ((xi == 0) || (i <= 0))
        {
            return -1/(r*r);
        }
        else
        {
            S  = DSlater_NS[i-1](r, xi);
            return S;
        }
    }
}

} // end namespace openbabel
