/**********************************************************************************
gaussian_integrals.cpp - coulombintegrals for Gaussian s-type orbitals

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

#include <math.h>
#include <stdio.h>

#include <openbabel/obutil.h>
#include <openbabel/gaussian_integrals.h>

static double sqr(double x)
{
    return x*x;
}

double Nuclear_GG(double r, double zeta)
{
    /* This routine may be called with zeta 0.
     * In that case it is a simple 1/r interaction.
     */
    if (zeta == 0)
    {
        if (r == 0)
        {
            return 0;
        }
        else
        {
            return 1.0/r;
        }
    }
    else if (r == 0)
    {
        return 2.0*zeta/sqrt(M_PI);
    }
    else
    {
        return erf(zeta*r)/r;
    }
}

double Coulomb_GG(double r, double zi, double zj)
{
    double zeff;

    /* This routine may be called with one or both zeta 0.
     * In that case it is either a Nuclear_GG or simple 1/r interaction.
     */
    if ((zi == 0) && (zj == 0))
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
    else
    {
        if (zi == 0)
        {
            zeff = zj;
        }
        else if (zj == 0)
        {
            zeff = zi;
        }
        else
        {
            zeff = zi*zj/sqrt(sqr(zi)+sqr(zj));
        }

        return Nuclear_GG(r, zeff);
    }
}

double DNuclear_GG(double r, double z)
{
    double rz = r * z;
    double dngg;

    if (z == 0)
    {
        if (r == 0)
        {
            return 0;
        }
        else
        {
            return -1.0/sqr(r);
        }
    }
    else if (r == 0)
    {
        dngg = 0;
    }
    else
    {
        dngg = (-1)*((2.0/sqrt(M_PI))*exp(-sqr(rz))*(z/r) - erf(rz)/sqr(r));
    }

    return dngg;
}

double DCoulomb_GG(double r, double zi, double zj)
{
    double zeff;

    if ((zi == 0) && (zj == 0))
    {
        if (r == 0)
        {
            return 0;
        }
        else
        {
            return -1/sqr(r);
        }
    }
    else
    {
        if (zi == 0)
        {
            zeff = zj;
        }
        else if (zj == 0)
        {
            zeff = zi;
        }
        else
        {
            zeff = zi*zj/sqrt(sqr(zi)+sqr(zj));
        }

        return DNuclear_GG(r, zeff);
    }
}
