/**********************************************************************
cpcomplex.cpp - Handle CpComplex class.

Copyright (C) 2023 by Jesus N. M.

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include <openbabel/atom.h>
#include <openbabel/cpcomplex.h>

using namespace std;

namespace OpenBabel {

    //
    // CpComplex member functions
    //
    //OBAtom* CpComplex::BeginAtomCp(OBAtomIterator& i)
    //{
    //    i = _cpAtoms.begin();
    //    return i == _cpAtoms.end() ? nullptr : (OBAtom*)*i;
    //}

    //OBAtom* CpComplex::NextAtomCp(OBAtomIterator& i)
    //{
    //    ++i;
    //    return i == _cpAtoms.end() ? nullptr : (OBAtom*)*i;
    //}
    void CpComplex::FindCentroid()
    {
        double sumX = 0.0, sumY = 0.0, sumZ = 0.0;
        for (int i = 0; i < _cpAtoms.size(); i++) {
            sumX += _cpAtoms[i]->GetX();
            sumY += _cpAtoms[i]->GetY();
            sumZ += _cpAtoms[i]->GetZ();
        }
        sumX = sumX / _cpAtoms.size();
        sumY = sumY / _cpAtoms.size();
        sumZ = sumZ / _cpAtoms.size();

        center.Set(sumX, sumY, sumZ);
    }


    //Zero based access method to vector
    //unsigned int CpComplex::GetCarbonIdx(int i) const
    //{
    //    if ((unsigned)i < 0 || (unsigned)i >= idx_carbons.size())
    //    {
    //        obErrorLog.ThrowError(__FUNCTION__, "Requested CarbonIdx Out of Range", obDebug);
    //    }

    //    return(idx_carbons[i]);
    //}

    vector3 CpComplex::GetCircleCoord(unsigned int i) {
        if ((unsigned)i < 0 || (unsigned)i >= circlePath.size()) {
            obErrorLog.ThrowError(__FUNCTION__, "Requested CircleCoord Out of Range", obDebug);
        }

        return (circlePath[i]);
    }

    void CpComplex::SetCircleCoord(unsigned int i, vector3 _v) {
        if ((unsigned)i < 0 || (unsigned)i >= circlePath.size()) {
            obErrorLog.ThrowError(__FUNCTION__, "Requested CircleCoord Out of Range", obDebug);
        }

        circlePath[i] = _v;
    }

    void CpComplex::SetCircleCoord(unsigned int i, double _vx, double _vy, double _vz) {
        if ((unsigned)i < 0 || (unsigned)i >= circlePath.size()) {
            obErrorLog.ThrowError(__FUNCTION__, "Requested CircleCoord Out of Range", obDebug);
        }

        circlePath[i].Set(_vx, _vy, _vz);
    }

}
