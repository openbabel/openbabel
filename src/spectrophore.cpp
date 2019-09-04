/**********************************************************************
Spectrophore.cpp - Spectrophore(TM) calculator

Copyright (C) 2005-2010 by Silicos NV

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

The algorithm in this software has been covered by patent WO2009146735.
However, Silicos NV and the inventors of the above mentioned patent assure
that no patent infringment claims will be issued against individuals or
institutions that use this software, as is stipulated by the GNU General
Public License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include <openbabel/spectrophore.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/obiter.h>

namespace OpenBabel
{


OBSpectrophore::OBSpectrophore(void)
:  _resolution(3.0)
,  _property(NULL)
,  _radii(NULL)
,  _oricoor(NULL)
,  _coor(NULL)
,  _beginProbe(0)
,  _endProbe(0)
,  _numberOfProbes(0)
{
   SetAccuracy(OBSpectrophore::AngStepSize20);
   SetStereo(OBSpectrophore::NoStereoSpecificProbes);
   SetNormalization(OBSpectrophore::NoNormalization);
}



OBSpectrophore::OBSpectrophore(const OBSpectrophore& s)
:  _resolution(s._resolution)
,  _property(NULL)
,  _radii(NULL)
,  _oricoor(NULL)
,  _coor(NULL)
,  _beginProbe(s._beginProbe)
,  _endProbe(s._endProbe)
,  _numberOfProbes(s._numberOfProbes)
,  _spectro(s._spectro)
{
   SetAccuracy(s.GetAccuracy());
   SetStereo(s.GetStereo());
   SetNormalization(s.GetNormalization());
}



OBSpectrophore&
OBSpectrophore::operator=(const OBSpectrophore& s)
{
   if (this != &s)
   {
      _resolution = s._resolution;
      _accuracy = s._accuracy;
      _beginProbe = s._beginProbe;
      _endProbe = s._endProbe;
      _numberOfProbes = s._numberOfProbes;
      _spectro = s._spectro;
      SetAccuracy(s.GetAccuracy());
      SetStereo(s.GetStereo());
      SetNormalization(s.GetNormalization());
   }
   return *this;
}



OBSpectrophore::~OBSpectrophore(void)
{
   _property = NULL;
   _radii = NULL;
   _oricoor = NULL;
   _coor = NULL;
}



void
OBSpectrophore::SetNormalization(const OBSpectrophore::NormalizationOption n)
{
   _normalization = n;
}



OBSpectrophore::NormalizationOption
OBSpectrophore::GetNormalization(void) const
{
    return _normalization;
}



void
OBSpectrophore::SetResolution(const double r)
{
    if (r > 0.0)
    {
        _resolution = r;
    }
    else
    {
        _resolution = 3.0;
    }
}



double
OBSpectrophore::GetResolution(void) const
{
    return _resolution;
}



void
OBSpectrophore::SetAccuracy(const OBSpectrophore::AccuracyOption a)
{
   _accuracy = a;

    // Update the rotation list
    _rotationStepList.clear();
    switch (_accuracy)
    {
      case OBSpectrophore::AngStepSize1:
         _rotationStepList.push_back(1);
         break;
      case OBSpectrophore::AngStepSize2:
         _rotationStepList.push_back(2);
         _rotationStepList.push_back(5);
         _rotationStepList.push_back(36);
         break;
      case OBSpectrophore::AngStepSize5:
         _rotationStepList.push_back(5);
         _rotationStepList.push_back(36);
         break;
      case OBSpectrophore::AngStepSize10:
         _rotationStepList.push_back(10);
         _rotationStepList.push_back(15);
         _rotationStepList.push_back(36);
         break;
      case OBSpectrophore::AngStepSize15:
         _rotationStepList.push_back(15);
         _rotationStepList.push_back(20);
         _rotationStepList.push_back(36);
         break;
      case OBSpectrophore::AngStepSize20:
         _rotationStepList.push_back(20);
         _rotationStepList.push_back(30);
         _rotationStepList.push_back(36);
         _rotationStepList.push_back(45);
         break;
      case OBSpectrophore::AngStepSize30:
         _rotationStepList.push_back(30);
         _rotationStepList.push_back(36);
         _rotationStepList.push_back(45);
         break;
      case OBSpectrophore::AngStepSize36:
         _rotationStepList.push_back(36);
         _rotationStepList.push_back(45);
         _rotationStepList.push_back(60);
         break;
      case OBSpectrophore::AngStepSize45:
         _rotationStepList.push_back(45);
         _rotationStepList.push_back(60);
         break;
      case OBSpectrophore::AngStepSize60:
         _rotationStepList.push_back(60);
         break;
    }
}



OBSpectrophore::AccuracyOption
OBSpectrophore::GetAccuracy(void) const
{
    return _accuracy;
}



void
OBSpectrophore::SetStereo(const OBSpectrophore::StereoOption stereo)
{
   _stereoFlag = stereo;
   _setBox();
}



OBSpectrophore::StereoOption
OBSpectrophore::GetStereo(void) const
{
   return _stereoFlag;
}



std::vector<double>
OBSpectrophore::GetSpectrophore(OpenBabel::OBMol* mol)
{
   // Clear the return variable
   _spectro.clear();

   // Atoms
   _nAtoms = mol->NumAtoms();
   if (_nAtoms < 3)
   {
      std::cerr << "OBSpectrophore::GetSpectrophore() error: not enough atoms in molecule" << std::endl;;
      return _spectro;
   }

   // Coordinate and property arrays
   _oricoor = new double*[_nAtoms];
   _coor = new double*[_nAtoms];
   double** REF1 = new double*[_nAtoms];
   double** REF2 = new double*[_nAtoms];
   _property = new double*[_nAtoms];
   for (unsigned int i = 0; i < _nAtoms; ++i)
   {
      _oricoor[i] = new double[3];
      _coor[i] = new double[3];
      REF1[i] = new double[3];
      REF2[i] = new double[3];
      _property[i] = new double[N_PROPERTIES];
   }

    // Atom radii and coordinates
   _radii = new double[_nAtoms];
   _getMoleculeData(mol);

   // Atom properties
   _calculateProperties(mol);

   // Shift molecule to its center of gravity and orient it in a standard way
   _orient();

   // Calculate the first spectrum to initiate values
   unsigned int sphoreSize(N_PROPERTIES * _numberOfProbes);
   std::vector<double> ENERGY(sphoreSize);
   std::vector<double> SPHORE(sphoreSize);

   // Properties
   _getBox(_oricoor);
   _getEnergies(_oricoor, &(ENERGY[0]));
   _initiateSpectrophore(&(ENERGY[0]), &(SPHORE[0]));

   // Rotate
   double psi;
   double cos_psi;
   double sin_psi;
   double theta;
   double cos_theta;
   double sin_theta;
   double phi;
   double cos_phi;
   double sin_phi;

   for (unsigned int i = 0; i < _rotationStepList.size(); ++i)
   {
      int rotationStep(_rotationStepList[i]);

      for (int iTheta = 0; iTheta < 180; iTheta += rotationStep)
      {
         theta = 0.017453292519943 * iTheta;
         cos_theta = cos(theta);
         sin_theta = sin(theta);
         _rotateY(_oricoor, REF1, cos_theta, sin_theta);

         for (int iPsi = 0; iPsi < 360; iPsi += rotationStep)
         {
            psi = 0.017453292519943 * iPsi;
            cos_psi = cos(psi);
            sin_psi = sin(psi);
            _rotateZ(REF1, REF2, cos_psi, sin_psi);

            for (int iPhi = 0; iPhi < 360; iPhi += rotationStep)
            {
               phi = 0.017453292519943 * iPhi;
               cos_phi = cos(phi);
               sin_phi = sin(phi);
               _rotateX(REF2, _coor, cos_phi, sin_phi);

               // Calculate energies
               _getBox(_coor);
               _getEnergies(_coor, &(ENERGY[0]));
               _updateSpectrophore(&(ENERGY[0]), &(SPHORE[0]));
            }
         }
      }
   }


   // Cleanup
   for (unsigned int i = 0; i < _nAtoms; ++i)
   {
      delete[] _property[i];
      delete[] _oricoor[i];
      delete[] REF1[i];
      delete[] REF2[i];
      delete[] _coor[i];
      _property[i] = NULL;
      _oricoor[i] = NULL;
      REF1[i] = NULL;
      REF2[i] = NULL;
      _coor[i] = NULL;
   }
   delete[] REF1;
   delete[] REF2;
   delete[] _radii;
   _radii = NULL;


   // Modify the actual sphore data
   _spectro.resize(sphoreSize);
   for (unsigned int i = 0; i < sphoreSize; ++i)
   {
      _spectro[i] = -100 * SPHORE[i];
   }

   // Normalisation
   double mean[N_PROPERTIES];
   double std[N_PROPERTIES];
   unsigned int m;
   for (unsigned int i(0); i < N_PROPERTIES; ++i)
   {
      mean[i] = 0.0;
      for (unsigned int n(0); n < 12; ++n)
      {
         m = (i * 12) + n;
         mean[i] += _spectro[m];
      }
      mean[i] /= 12.0;
      std[i] = 0.0;
      for (unsigned int n(0); n < 12; ++n)
      {
         m = (i * 12) + n;
         std[i] += (_spectro[m] - mean[i]) * (_spectro[m] - mean[i]);
      }
      std[i] /= 11.0;
      std[i] = sqrt(std[i]);
   }
   if ((_normalization == OBSpectrophore::NormalizationTowardsZeroMean) ||
       (_normalization == OBSpectrophore::NormalizationTowardsZeroMeanAndUnitStd))
   {
      for (unsigned int i(0); i < N_PROPERTIES; ++i)
      {
         for (unsigned int n(0); n < 12; ++n)
         {
            m = (i * 12) + n;
            _spectro[m] -= mean[i];
         }
      }
   }
   if ((_normalization == OBSpectrophore::NormalizationTowardsUnitStd) ||
       (_normalization == OBSpectrophore::NormalizationTowardsZeroMeanAndUnitStd))
   {
      for (unsigned int i(0); i < N_PROPERTIES; ++i)
      {
         for (unsigned int n(0); n < 12; ++n)
         {
            m = (i * 12) + n;
            _spectro[m] /= std[i];
         }
      }
   }

   // Return
   return _spectro;
}



void
OBSpectrophore::_getEnergies(double** c, double* e)
{
   double d;
   double x2;
   double y2;
   double z2;

   // Distance between atom and each boxpoint
   for (unsigned int boxPoint = 0; boxPoint < 12; ++boxPoint)
   {
      // Reset value at box point
      for (unsigned int prop = 0; prop < N_PROPERTIES; ++prop)
      {
         _boxPoint[boxPoint].v[prop] = 0.0;
      }

      // Calculate squared distances between atoms and boxpoints
      for (unsigned int atom = 0; atom < _nAtoms; ++atom)
      {
         x2 = _boxPoint[boxPoint].x - c[atom][0];
         y2 = _boxPoint[boxPoint].y - c[atom][1];
         z2 = _boxPoint[boxPoint].z - c[atom][2];
         d = sqrt((x2 * x2) + (y2 * y2) + (z2 * z2));

         // Calculate the potential at each box point
         for (unsigned int prop = 0; prop < N_PROPERTIES; ++prop)
         {
            _boxPoint[boxPoint].v[prop] += _property[atom][prop] / d;
         }
      }
   }

   // Reset to zero
   for (unsigned int i = 0; i < N_PROPERTIES * _numberOfProbes; ++i)
   {
      e[i] = 0.0;
   }

   unsigned int index;
   for (unsigned int prop = 0; prop < N_PROPERTIES; ++prop)
   {
      for (unsigned int boxPoint = 0; boxPoint < 12; ++boxPoint)
      {
         // Calculate for each probe the total interaction energy
         for (unsigned int probe = _beginProbe; probe < _endProbe; ++probe)
         {
            index = prop * _numberOfProbes + (probe - _beginProbe);
            e[index] += _boxPoint[boxPoint].v[prop] *_probe[probe].value[boxPoint];
         }
      }
   }
}



void
OBSpectrophore::_initiateSpectrophore(double* e, double* s)
{
    for (unsigned int i = 0; i < N_PROPERTIES * _numberOfProbes; ++i)
    {
        s[i] = e[i];
    }
}



void
OBSpectrophore::_updateSpectrophore(double* ENERGY, double* SPHORE)
{
    for (unsigned int i = 0; i < N_PROPERTIES * _numberOfProbes; ++i)
    {
        SPHORE[i] = MIN(ENERGY[i], SPHORE[i]);
    }
}



void
OBSpectrophore::_getMoleculeData(OpenBabel::OBMol* mol)
{
   unsigned int a(0);
   unsigned int n;
   for (OpenBabel::OBMolAtomIter atom(mol); atom; ++atom)
   {
      // Coordinates
      _oricoor[a][0] = atom->GetX();
      _oricoor[a][1] = atom->GetY();
      _oricoor[a][2] = atom->GetZ();

      // Radii
      n = (unsigned int) atom->GetAtomicNum();
      switch (n)
      {
         case 1:  // H
            _radii[a] = 1.20;
            break;
         case 3:  // Li
            _radii[a] = 1.82;
            break;
         case 5:  // B
            _radii[a] = 2.00;
            break;
         case 6:  // C
            _radii[a] = 1.70;
            break;
         case 7:  // N
            _radii[a] = 1.55;
            break;
         case 8:  // O
            _radii[a] = 1.52;
            break;
         case 9:  // F
            _radii[a] = 1.47;
            break;
         case 11: // Na
            _radii[a] = 2.27;
            break;
         case 12: // Mg
            _radii[a] = 1.73;
            break;
         case 14: // Si
            _radii[a] = 2.10;
            break;
         case 15: // P
            _radii[a] = 1.80;
            break;
         case 16: // S
            _radii[a] = 1.80;
            break;
        case 17: // Cl
            _radii[a] = 1.75;
            break;
         case 19: // K
            _radii[a] = 2.75;
            break;
         case 20: // Ca
            _radii[a] = 2.00;
            break;
         case 26: // Fe
            _radii[a] = 1.10;
            break;
         case 29: // Cu
            _radii[a] = 1.40;
            break;
         case 30: // Zn
            _radii[a] = 1.39;
            break;
         case 35: // Br
            _radii[a] = 1.85;
            break;
         case 53: // I
            _radii[a] = 1.98;
            break;
         default:
            _radii[a] = 1.50;
      }

      // Counter
      ++a;
   }
}


void
OBSpectrophore::_getBox(double** c)
{
   // Calculate molecular extends
   double xm(c[0][0] - _radii[0]);
   double xp(c[0][0] + _radii[0]);
   double ym(c[0][1] - _radii[0]);
   double yp(c[0][1] + _radii[0]);
   double zm(c[0][2] - _radii[0]);
   double zp(c[0][2] + _radii[0]);
   for (unsigned int i = 1; i < _nAtoms; ++i)
   {
      xm = MIN(c[i][0] - _radii[i], xm);
      xp = MAX(c[i][0] + _radii[i], xp);
      ym = MIN(c[i][1] - _radii[i], ym);
      yp = MAX(c[i][1] + _radii[i], yp);
      zm = MIN(c[i][2] - _radii[i], zm);
      zp = MAX(c[i][2] + _radii[i], zp);
   }
   xm -= _resolution;
   xp += _resolution;
   ym -= _resolution;
   yp += _resolution;
   zm -= _resolution;
   zp += _resolution;
   double xh((xp + xm) / 2.0);
   double yh((yp + ym) / 2.0);
   double zh((zp + zm) / 2.0);

   // Box points
   _boxPoint[0].x =  xh;
   _boxPoint[1].x =  xp;
   _boxPoint[2].x =  xh;
   _boxPoint[3].x =  xm;
   _boxPoint[4].x =  xm;
   _boxPoint[5].x =  xp;
   _boxPoint[6].x =  xm;
   _boxPoint[7].x =  xp;
   _boxPoint[8].x =  xp;
   _boxPoint[9].x =  xh;
   _boxPoint[10].x = xm;
   _boxPoint[11].x = xh;

   _boxPoint[0].y =  ym;
   _boxPoint[1].y =  yh;
   _boxPoint[2].y =  yp;
   _boxPoint[3].y =  yh;
   _boxPoint[4].y =  ym;
   _boxPoint[5].y =  ym;
   _boxPoint[6].y =  yp;
   _boxPoint[7].y =  yp;
   _boxPoint[8].y =  yh;
   _boxPoint[9].y =  ym;
   _boxPoint[10].y = yh;
   _boxPoint[11].y = yp;

   _boxPoint[0].z =  zp;
   _boxPoint[1].z =  zp;
   _boxPoint[2].z =  zp;
   _boxPoint[3].z =  zp;
   _boxPoint[4].z =  zh;
   _boxPoint[5].z =  zh;
   _boxPoint[6].z =  zh;
   _boxPoint[7].z =  zh;
   _boxPoint[8].z =  zm;
   _boxPoint[9].z =  zm;
   _boxPoint[10].z = zm;
   _boxPoint[11].z = zm;
}



void
OBSpectrophore::_rotateX(double** oc, double** nc, const double c, const double s)
{
    for (unsigned int i = 0; i < _nAtoms; ++i)
    {
        nc[i][0] =  oc[i][0];
        nc[i][1] =  oc[i][1] * c + oc[i][2] * s;
        nc[i][2] = -oc[i][1] * s + oc[i][2] * c;
    }
}



void
OBSpectrophore::_rotateY(double** oc, double** nc, const double c, const double s)
{
    for (unsigned int i = 0; i < _nAtoms; ++i)
    {
        nc[i][0] =  oc[i][0] * c + oc[i][2] * s;
        nc[i][1] =  oc[i][1];
        nc[i][2] = -oc[i][0] * s + oc[i][2] * c;
    }
}



void
OBSpectrophore::_rotateZ(double** oc, double** nc, const double c, const double s)
{
    for (unsigned int i = 0; i < _nAtoms; ++i)
    {
        nc[i][0] =  oc[i][0] * c + oc[i][1] * s;
        nc[i][1] = -oc[i][0] * s + oc[i][1] * c;
        nc[i][2] =  oc[i][2];
    }
}



void
OBSpectrophore::_orient(void)
{
    // Center molecule around its COG
    double COG[3] = {0.0, 0.0, 0.0};
    for (unsigned int i = 0; i < _nAtoms; ++i)
    {
        for (unsigned int j = 0; j < 3; ++j)
        {
            COG[j] += _oricoor[i][j];
        }
    }
    for (unsigned int j = 0; j < 3; ++j)
    {
        COG[j] /= _nAtoms;
    }
    for (unsigned int i = 0; i < _nAtoms; ++i)
    {
        for (unsigned int j = 0; j < 3; ++j)
        {
            _oricoor[i][j] -= COG[j];
        }
    }

    // Determine atom that is furthest away from origin
    double maxDistance(0.0);
    int maxAtom(0);
    double d;
    for (unsigned int i = 0; i < _nAtoms; ++i)
    {
        d = _oricoor[i][0]*_oricoor[i][0] +  _oricoor[i][1]*_oricoor[i][1] + _oricoor[i][2]*_oricoor[i][2];
        if (d > maxDistance)
        {
            maxDistance = d;
            maxAtom = i;
        }
    }

    // Rotate all atoms along z-axis
    double angle(-atan2(_oricoor[maxAtom][1], _oricoor[maxAtom][0]));
    double x;
    double y;
    for (unsigned int i = 0; i < _nAtoms; ++i)
    {
        x = cos(angle) * _oricoor[i][0] - sin(angle) * _oricoor[i][1];
        y = sin(angle) * _oricoor[i][0] + cos(angle) * _oricoor[i][1];
        _oricoor[i][0] = x;
        _oricoor[i][1] = y;
    }

    // Rotate all atoms along y-axis to place the maxAtom on z
    angle = -atan2(_oricoor[maxAtom][0], _oricoor[maxAtom][2]);
    double z;
    for (unsigned int i = 0; i < _nAtoms; ++i)
    {
        x = cos(angle) * _oricoor[i][0] + sin(angle) * _oricoor[i][2];
        z = cos(angle) * _oricoor[i][2] - sin(angle) * _oricoor[i][0];
        _oricoor[i][0] = x;
        _oricoor[i][2] = z;
    }

    // Center molecule again around its COG
    COG[0] = 0.0;
    COG[1] = 0.0;
    COG[2] = 0.0;
    for (unsigned int i = 0; i < _nAtoms; ++i)
    {
        for (unsigned int j = 0; j < 3; ++j)
        {
            COG[j] += _oricoor[i][j];
        }
    }
    for (unsigned int j = 0; j < 3; ++j)
    {
        COG[j] /= _nAtoms;
    }
    for (unsigned int i = 0; i < _nAtoms; ++i)
    {
        for (unsigned int j = 0; j < 3; ++j)
        {
            _oricoor[i][j] -= COG[j];
        }
    }
}



void
OBSpectrophore::_setBox(void)
{
   // Define boundaries
   switch (_stereoFlag)
   {
      case OBSpectrophore::NoStereoSpecificProbes:  // No stereo
         _beginProbe = 0;
         _endProbe = 12;
         _numberOfProbes = 12;
         break;
      case OBSpectrophore::UniqueStereoSpecificProbes:  // Mirror stereo
         _beginProbe = 12;
         _endProbe = 30;
         _numberOfProbes = 18;
         break;
      case OBSpectrophore::MirrorStereoSpecificProbes:  // Unique stereo
         _beginProbe = 30;
         _endProbe = 48;
         _numberOfProbes = 18;
         break;
      case OBSpectrophore::AllStereoSpecificProbes:  // All stereo
         _beginProbe = 12;
         _endProbe = 48;
         _numberOfProbes = 36;
         break;
   }

   // Dodecapole - non-stereo - probe 1
   _probe[0].value[0] = +1;
   _probe[0].value[1] = +1;
   _probe[0].value[2] = -1;
   _probe[0].value[3] = -1;
   _probe[0].value[4] = -1;
   _probe[0].value[5] = +1;
   _probe[0].value[6] = +1;
   _probe[0].value[7] = -1;
   _probe[0].value[8] = -1;
   _probe[0].value[9] = -1;
   _probe[0].value[10] = +1;
   _probe[0].value[11] = +1;

   // Dodecapole - non-stereo - probe 2
   _probe[1].value[0] = +1;
   _probe[1].value[1] = +1;
   _probe[1].value[2] = -1;
   _probe[1].value[3] = -1;
   _probe[1].value[4] = +1;
   _probe[1].value[5] = -1;
   _probe[1].value[6] = -1;
   _probe[1].value[7] = +1;
   _probe[1].value[8] = -1;
   _probe[1].value[9] = -1;
   _probe[1].value[10] = +1;
   _probe[1].value[11] = +1;

   // Dodecapole - non-stereo - probe 3
   _probe[2].value[0] = +1;
   _probe[2].value[1] = +1;
   _probe[2].value[2] = +1;
   _probe[2].value[3] = -1;
   _probe[2].value[4] = -1;
   _probe[2].value[5] = -1;
   _probe[2].value[6] = -1;
   _probe[2].value[7] = -1;
   _probe[2].value[8] = -1;
   _probe[2].value[9] = +1;
   _probe[2].value[10] = +1;
   _probe[2].value[11] = +1;

   // Dodecapole - non-stereo - probe 4
   _probe[3].value[0] = +1;
   _probe[3].value[1] = +1;
   _probe[3].value[2] = +1;
   _probe[3].value[3] = -1;
   _probe[3].value[4] = -1;
   _probe[3].value[5] = -1;
   _probe[3].value[6] = -1;
   _probe[3].value[7] = -1;
   _probe[3].value[8] = +1;
   _probe[3].value[9] = +1;
   _probe[3].value[10] = -1;
   _probe[3].value[11] = +1;

   // Dodecapole - non-stereo - probe 5
   _probe[4].value[0] = +1;
   _probe[4].value[1] = +1;
   _probe[4].value[2] = +1;
   _probe[4].value[3] = -1;
   _probe[4].value[4] = -1;
   _probe[4].value[5] = +1;
   _probe[4].value[6] = -1;
   _probe[4].value[7] = +1;
   _probe[4].value[8] = -1;
   _probe[4].value[9] = -1;
   _probe[4].value[10] = +1;
   _probe[4].value[11] = -1;

   // Dodecapole - non-stereo - probe 6
   _probe[5].value[0] = +1;
   _probe[5].value[1] = +1;
   _probe[5].value[2] = +1;
   _probe[5].value[3] = -1;
   _probe[5].value[4] = +1;
   _probe[5].value[5] = -1;
   _probe[5].value[6] = +1;
   _probe[5].value[7] = -1;
   _probe[5].value[8] = -1;
   _probe[5].value[9] = -1;
   _probe[5].value[10] = +1;
   _probe[5].value[11] = -1;

   // Dodecapole - non-stereo - probe 7
   _probe[6].value[0] = +1;
   _probe[6].value[1] = +1;
   _probe[6].value[2] = +1;
   _probe[6].value[3] = -1;
   _probe[6].value[4] = +1;
   _probe[6].value[5] = -1;
   _probe[6].value[6] = +1;
   _probe[6].value[7] = -1;
   _probe[6].value[8] = +1;
   _probe[6].value[9] = -1;
   _probe[6].value[10] = -1;
   _probe[6].value[11] = -1;

   // Dodecapole - non-stereo - probe 8
   _probe[7].value[0] = +1;
   _probe[7].value[1] = +1;
   _probe[7].value[2] = +1;
   _probe[7].value[3] = +1;
   _probe[7].value[4] = -1;
   _probe[7].value[5] = -1;
   _probe[7].value[6] = -1;
   _probe[7].value[7] = -1;
   _probe[7].value[8] = +1;
   _probe[7].value[9] = -1;
   _probe[7].value[10] = +1;
   _probe[7].value[11] = -1;

   // Dodecapole - non-stereo - probe 9
   _probe[8].value[0] = +1;
   _probe[8].value[1] = +1;
   _probe[8].value[2] = +1;
   _probe[8].value[3] = +1;
   _probe[8].value[4] = -1;
   _probe[8].value[5] = -1;
   _probe[8].value[6] = -1;
   _probe[8].value[7] = -1;
   _probe[8].value[8] = +1;
   _probe[8].value[9] = +1;
   _probe[8].value[10] = -1;
   _probe[8].value[11] = -1;

   // Dodecapole - non-stereo - probe 10
   _probe[9].value[0] = +1;
   _probe[9].value[1] = +1;
   _probe[9].value[2] = +1;
   _probe[9].value[3] = +1;
   _probe[9].value[4] = +1;
   _probe[9].value[5] = -1;
   _probe[9].value[6] = -1;
   _probe[9].value[7] = +1;
   _probe[9].value[8] = -1;
   _probe[9].value[9] = -1;
   _probe[9].value[10] = -1;
   _probe[9].value[11] = -1;

   // Dodecapole - non-stereo - probe 11
   _probe[10].value[0] = +1;
   _probe[10].value[1] = +1;
   _probe[10].value[2] = +1;
   _probe[10].value[3] = +1;
   _probe[10].value[4] = +1;
   _probe[10].value[5] = +1;
   _probe[10].value[6] = -1;
   _probe[10].value[7] = -1;
   _probe[10].value[8] = -1;
   _probe[10].value[9] = -1;
   _probe[10].value[10] = -1;
   _probe[10].value[11] = -1;

   // Dodecapole - non-stereo - probe 12
   _probe[11].value[0] = +1;
   _probe[11].value[1] = +1;
   _probe[11].value[2] = +1;
   _probe[11].value[3] = -1;
   _probe[11].value[4] = -1;
   _probe[11].value[5] = +1;
   _probe[11].value[6] = -1;
   _probe[11].value[7] = -1;
   _probe[11].value[8] = -1;
   _probe[11].value[9] = +1;
   _probe[11].value[10] = -1;
   _probe[11].value[11] = +1;

   // Dodecapole - mirror-stereo - probe 1
   _probe[12].value[0] = +1;
   _probe[12].value[1] = +1;
   _probe[12].value[2] = -1;
   _probe[12].value[3] = -1;
   _probe[12].value[4] = -1;
   _probe[12].value[5] = -1;
   _probe[12].value[6] = +1;
   _probe[12].value[7] = +1;
   _probe[12].value[8] = -1;
   _probe[12].value[9] = +1;
   _probe[12].value[10] = +1;
   _probe[12].value[11] = -1;

   // Dodecapole - mirror-stereo - probe 2
   _probe[13].value[0] = +1;
   _probe[13].value[1] = +1;
   _probe[13].value[2] = +1;
   _probe[13].value[3] = -1;
   _probe[13].value[4] = -1;
   _probe[13].value[5] = -1;
   _probe[13].value[6] = -1;
   _probe[13].value[7] = -1;
   _probe[13].value[8] = +1;
   _probe[13].value[9] = -1;
   _probe[13].value[10] = +1;
   _probe[13].value[11] = +1;

   // Dodecapole - mirror-stereo - probe 3
   _probe[14].value[0] = +1;
   _probe[14].value[1] = +1;
   _probe[14].value[2] = +1;
   _probe[14].value[3] = -1;
   _probe[14].value[4] = -1;
   _probe[14].value[5] = -1;
   _probe[14].value[6] = -1;
   _probe[14].value[7] = +1;
   _probe[14].value[8] = -1;
   _probe[14].value[9] = +1;
   _probe[14].value[10] = +1;
   _probe[14].value[11] = -1;

   // Dodecapole - mirror-stereo - probe 4
   _probe[15].value[0] = +1;
   _probe[15].value[1] = +1;
   _probe[15].value[2] = +1;
   _probe[15].value[3] = -1;
   _probe[15].value[4] = -1;
   _probe[15].value[5] = -1;
   _probe[15].value[6] = +1;
   _probe[15].value[7] = -1;
   _probe[15].value[8] = -1;
   _probe[15].value[9] = +1;
   _probe[15].value[10] = -1;
   _probe[15].value[11] = +1;

   // Dodecapole - mirror-stereo - probe 5
   _probe[16].value[0] = +1;
   _probe[16].value[1] = +1;
   _probe[16].value[2] = +1;
   _probe[16].value[3] = -1;
   _probe[16].value[4] = -1;
   _probe[16].value[5] = -1;
   _probe[16].value[6] = +1;
   _probe[16].value[7] = -1;
   _probe[16].value[8] = -1;
   _probe[16].value[9] = +1;
   _probe[16].value[10] = +1;
   _probe[16].value[11] = -1;

   // Dodecapole - mirror-stereo - probe 6
   _probe[17].value[0] = +1;
   _probe[17].value[1] = +1;
   _probe[17].value[2] = +1;
   _probe[17].value[3] = -1;
   _probe[17].value[4] = -1;
   _probe[17].value[5] = -1;
   _probe[17].value[6] = +1;
   _probe[17].value[7] = -1;
   _probe[17].value[8] = +1;
   _probe[17].value[9] = -1;
   _probe[17].value[10] = +1;
   _probe[17].value[11] = -1;

   // Dodecapole - mirror-stereo - probe 7
   _probe[18].value[0] = +1;
   _probe[18].value[1] = +1;
   _probe[18].value[2] = +1;
   _probe[18].value[3] = -1;
   _probe[18].value[4] = -1;
   _probe[18].value[5] = -1;
   _probe[18].value[6] = +1;
   _probe[18].value[7] = -1;
   _probe[18].value[8] = +1;
   _probe[18].value[9] = +1;
   _probe[18].value[10] = -1;
   _probe[18].value[11] = -1;

   // Dodecapole - mirror-stereo - probe 8
   _probe[19].value[0] = +1;
   _probe[19].value[1] = +1;
   _probe[19].value[2] = +1;
   _probe[19].value[3] = -1;
   _probe[19].value[4] = -1;
   _probe[19].value[5] = -1;
   _probe[19].value[6] = +1;
   _probe[19].value[7] = +1;
   _probe[19].value[8] = -1;
   _probe[19].value[9] = +1;
   _probe[19].value[10] = -1;
   _probe[19].value[11] = -1;

   // Dodecapole - mirror-stereo - probe 9
   _probe[20].value[0] = +1;
   _probe[20].value[1] = +1;
   _probe[20].value[2] = +1;
   _probe[20].value[3] = -1;
   _probe[20].value[4] = -1;
   _probe[20].value[5] = -1;
   _probe[20].value[6] = +1;
   _probe[20].value[7] = +1;
   _probe[20].value[8] = +1;
   _probe[20].value[9] = -1;
   _probe[20].value[10] = -1;
   _probe[20].value[11] = -1;

   // Dodecapole - mirror-stereo - probe 10
   _probe[21].value[0] = +1;
   _probe[21].value[1] = +1;
   _probe[21].value[2] = +1;
   _probe[21].value[3] = -1;
   _probe[21].value[4] = -1;
   _probe[21].value[5] = +1;
   _probe[21].value[6] = +1;
   _probe[21].value[7] = -1;
   _probe[21].value[8] = -1;
   _probe[21].value[9] = -1;
   _probe[21].value[10] = -1;
   _probe[21].value[11] = +1;

   // Dodecapole - mirror-stereo - probe 11
   _probe[22].value[0] = +1;
   _probe[22].value[1] = +1;
   _probe[22].value[2] = +1;
   _probe[22].value[3] = -1;
   _probe[22].value[4] = -1;
   _probe[22].value[5] = +1;
   _probe[22].value[6] = +1;
   _probe[22].value[7] = -1;
   _probe[22].value[8] = -1;
   _probe[22].value[9] = -1;
   _probe[22].value[10] = +1;
   _probe[22].value[11] = -1;

   // Dodecapole - mirror-stereo - probe 12
   _probe[23].value[0] = +1;
   _probe[23].value[1] = +1;
   _probe[23].value[2] = +1;
   _probe[23].value[3] = -1;
   _probe[23].value[4] = -1;
   _probe[23].value[5] = +1;
   _probe[23].value[6] = +1;
   _probe[23].value[7] = -1;
   _probe[23].value[8] = -1;
   _probe[23].value[9] = +1;
   _probe[23].value[10] = -1;
   _probe[23].value[11] = -1;

   // Dodecapole - mirror-stereo - probe 13
   _probe[24].value[0] = +1;
   _probe[24].value[1] = +1;
   _probe[24].value[2] = +1;
   _probe[24].value[3] = -1;
   _probe[24].value[4] = -1;
   _probe[24].value[5] = +1;
   _probe[24].value[6] = +1;
   _probe[24].value[7] = +1;
   _probe[24].value[8] = -1;
   _probe[24].value[9] = -1;
   _probe[24].value[10] = -1;
   _probe[24].value[11] = -1;

   // Dodecapole - mirror-stereo - probe 14
   _probe[25].value[0] = +1;
   _probe[25].value[1] = +1;
   _probe[25].value[2] = +1;
   _probe[25].value[3] = -1;
   _probe[25].value[4] = +1;
   _probe[25].value[5] = -1;
   _probe[25].value[6] = -1;
   _probe[25].value[7] = +1;
   _probe[25].value[8] = +1;
   _probe[25].value[9] = -1;
   _probe[25].value[10] = -1;
   _probe[25].value[11] = -1;

   // Dodecapole - mirror-stereo - probe 15
   _probe[26].value[0] = +1;
   _probe[26].value[1] = +1;
   _probe[26].value[2] = +1;
   _probe[26].value[3] = -1;
   _probe[26].value[4] = +1;
   _probe[26].value[5] = -1;
   _probe[26].value[6] = +1;
   _probe[26].value[7] = -1;
   _probe[26].value[8] = -1;
   _probe[26].value[9] = -1;
   _probe[26].value[10] = -1;
   _probe[26].value[11] = +1;

   // Dodecapole - mirror-stereo - probe 16
   _probe[27].value[0] = +1;
   _probe[27].value[1] = +1;
   _probe[27].value[2] = +1;
   _probe[27].value[3] = -1;
   _probe[27].value[4] = +1;
   _probe[27].value[5] = -1;
   _probe[27].value[6] = +1;
   _probe[27].value[7] = +1;
   _probe[27].value[8] = -1;
   _probe[27].value[9] = -1;
   _probe[27].value[10] = -1;
   _probe[27].value[11] = -1;

   // Dodecapole - mirror-stereo - probe 17
   _probe[28].value[0] = +1;
   _probe[28].value[1] = +1;
   _probe[28].value[2] = +1;
   _probe[28].value[3] = +1;
   _probe[28].value[4] = +1;
   _probe[28].value[5] = -1;
   _probe[28].value[6] = -1;
   _probe[28].value[7] = -1;
   _probe[28].value[8] = -1;
   _probe[28].value[9] = -1;
   _probe[28].value[10] = -1;
   _probe[28].value[11] = +1;

   // Dodecapole - mirror-stereo - probe 18
   _probe[29].value[0] = +1;
   _probe[29].value[1] = +1;
   _probe[29].value[2] = +1;
   _probe[29].value[3] = +1;
   _probe[29].value[4] = +1;
   _probe[29].value[5] = -1;
   _probe[29].value[6] = -1;
   _probe[29].value[7] = -1;
   _probe[29].value[8] = -1;
   _probe[29].value[9] = -1;
   _probe[29].value[10] = +1;
   _probe[29].value[11] = -1;

   // Dodecapole - unique-stereo - probe 1
   _probe[30].value[0] = +1;
   _probe[30].value[1] = +1;
   _probe[30].value[2] = -1;
   _probe[30].value[3] = -1;
   _probe[30].value[4] = +1;
   _probe[30].value[5] = -1;
   _probe[30].value[6] = +1;
   _probe[30].value[7] = -1;
   _probe[30].value[8] = +1;
   _probe[30].value[9] = -1;
   _probe[30].value[10] = -1;
   _probe[30].value[11] = +1;

   // Dodecapole - unique-stereo - probe 2
   _probe[31].value[0] = +1;
   _probe[31].value[1] = +1;
   _probe[31].value[2] = +1;
   _probe[31].value[3] = -1;
   _probe[31].value[4] = -1;
   _probe[31].value[5] = -1;
   _probe[31].value[6] = -1;
   _probe[31].value[7] = -1;
   _probe[31].value[8] = +1;
   _probe[31].value[9] = +1;
   _probe[31].value[10] = +1;
   _probe[31].value[11] = -1;

   // Dodecapole - unique-stereo - probe 3
   _probe[32].value[0] = +1;
   _probe[32].value[1] = +1;
   _probe[32].value[2] = +1;
   _probe[32].value[3] = -1;
   _probe[32].value[4] = -1;
   _probe[32].value[5] = +1;
   _probe[32].value[6] = -1;
   _probe[32].value[7] = -1;
   _probe[32].value[8] = -1;
   _probe[32].value[9] = -1;
   _probe[32].value[10] = +1;
   _probe[32].value[11] = +1;

   // Dodecapole - unique-stereo - probe 4
   _probe[33].value[0] = +1;
   _probe[33].value[1] = +1;
   _probe[33].value[2] = +1;
   _probe[33].value[3] = -1;
   _probe[33].value[4] = +1;
   _probe[33].value[5] = -1;
   _probe[33].value[6] = -1;
   _probe[33].value[7] = -1;
   _probe[33].value[8] = -1;
   _probe[33].value[9] = +1;
   _probe[33].value[10] = -1;
   _probe[33].value[11] = +1;

   // Dodecapole - unique-stereo - probe 5
   _probe[34].value[0] = +1;
   _probe[34].value[1] = +1;
   _probe[34].value[2] = +1;
   _probe[34].value[3] = -1;
   _probe[34].value[4] = +1;
   _probe[34].value[5] = -1;
   _probe[34].value[6] = -1;
   _probe[34].value[7] = -1;
   _probe[34].value[8] = -1;
   _probe[34].value[9] = -1;
   _probe[34].value[10] = +1;
   _probe[34].value[11] = +1;

   // Dodecapole - unique-stereo - probe 6
   _probe[35].value[0] = +1;
   _probe[35].value[1] = +1;
   _probe[35].value[2] = +1;
   _probe[35].value[3] = -1;
   _probe[35].value[4] = +1;
   _probe[35].value[5] = -1;
   _probe[35].value[6] = -1;
   _probe[35].value[7] = -1;
   _probe[35].value[8] = +1;
   _probe[35].value[9] = -1;
   _probe[35].value[10] = +1;
   _probe[35].value[11] = -1;

   // Dodecapole - unique-stereo - probe 7
   _probe[36].value[0] = +1;
   _probe[36].value[1] = +1;
   _probe[36].value[2] = +1;
   _probe[36].value[3] = -1;
   _probe[36].value[4] = +1;
   _probe[36].value[5] = -1;
   _probe[36].value[6] = -1;
   _probe[36].value[7] = -1;
   _probe[36].value[8] = +1;
   _probe[36].value[9] = -1;
   _probe[36].value[10] = -1;
   _probe[36].value[11] = +1;

   // Dodecapole - unique-stereo - probe 8
   _probe[37].value[0] = +1;
   _probe[37].value[1] = +1;
   _probe[37].value[2] = +1;
   _probe[37].value[3] = -1;
   _probe[37].value[4] = +1;
   _probe[37].value[5] = +1;
   _probe[37].value[6] = -1;
   _probe[37].value[7] = -1;
   _probe[37].value[8] = -1;
   _probe[37].value[9] = -1;
   _probe[37].value[10] = -1;
   _probe[37].value[11] = +1;

   // Dodecapole - unique-stereo - probe 9
   _probe[38].value[0] = +1;
   _probe[38].value[1] = +1;
   _probe[38].value[2] = +1;
   _probe[38].value[3] = -1;
   _probe[38].value[4] = +1;
   _probe[38].value[5] = +1;
   _probe[38].value[6] = -1;
   _probe[38].value[7] = -1;
   _probe[38].value[8] = +1;
   _probe[38].value[9] = -1;
   _probe[38].value[10] = -1;
   _probe[38].value[11] = -1;

   // Dodecapole - unique-stereo - probe 10
   _probe[39].value[0] = +1;
   _probe[39].value[1] = +1;
   _probe[39].value[2] = +1;
   _probe[39].value[3] = -1;
   _probe[39].value[4] = +1;
   _probe[39].value[5] = -1;
   _probe[39].value[6] = -1;
   _probe[39].value[7] = +1;
   _probe[39].value[8] = -1;
   _probe[39].value[9] = +1;
   _probe[39].value[10] = -1;
   _probe[39].value[11] = -1;

   // Dodecapole - unique-stereo - probe 11
   _probe[40].value[0] = +1;
   _probe[40].value[1] = +1;
   _probe[40].value[2] = +1;
   _probe[40].value[3] = -1;
   _probe[40].value[4] = +1;
   _probe[40].value[5] = -1;
   _probe[40].value[6] = -1;
   _probe[40].value[7] = +1;
   _probe[40].value[8] = -1;
   _probe[40].value[9] = -1;
   _probe[40].value[10] = +1;
   _probe[40].value[11] = -1;

   // Dodecapole - unique-stereo - probe 12
   _probe[41].value[0] = +1;
   _probe[41].value[1] = +1;
   _probe[41].value[2] = +1;
   _probe[41].value[3] = -1;
   _probe[41].value[4] = +1;
   _probe[41].value[5] = -1;
   _probe[41].value[6] = -1;
   _probe[41].value[7] = +1;
   _probe[41].value[8] = -1;
   _probe[41].value[9] = -1;
   _probe[41].value[10] = -1;
   _probe[41].value[11] = +1;

   // Dodecapole - unique-stereo - probe 13
   _probe[42].value[0] = +1;
   _probe[42].value[1] = +1;
   _probe[42].value[2] = +1;
   _probe[42].value[3] = -1;
   _probe[42].value[4] = +1;
   _probe[42].value[5] = +1;
   _probe[42].value[6] = -1;
   _probe[42].value[7] = +1;
   _probe[42].value[8] = -1;
   _probe[42].value[9] = -1;
   _probe[42].value[10] = -1;
   _probe[42].value[11] = -1;

   // Dodecapole - unique-stereo - probe 14
   _probe[43].value[0] = +1;
   _probe[43].value[1] = +1;
   _probe[43].value[2] = +1;
   _probe[43].value[3] = -1;
   _probe[43].value[4] = +1;
   _probe[43].value[5] = +1;
   _probe[43].value[6] = -1;
   _probe[43].value[7] = -1;
   _probe[43].value[8] = -1;
   _probe[43].value[9] = -1;
   _probe[43].value[10] = +1;
   _probe[43].value[11] = -1;

   // Dodecapole - unique-stereo - probe 15
   _probe[44].value[0] = +1;
   _probe[44].value[1] = +1;
   _probe[44].value[2] = +1;
   _probe[44].value[3] = -1;
   _probe[44].value[4] = +1;
   _probe[44].value[5] = -1;
   _probe[44].value[6] = +1;
   _probe[44].value[7] = -1;
   _probe[44].value[8] = -1;
   _probe[44].value[9] = +1;
   _probe[44].value[10] = -1;
   _probe[44].value[11] = -1;

   // Dodecapole - unique-stereo - probe 16
   _probe[45].value[0] = +1;
   _probe[45].value[1] = +1;
   _probe[45].value[2] = +1;
   _probe[45].value[3] = -1;
   _probe[45].value[4] = +1;
   _probe[45].value[5] = +1;
   _probe[45].value[6] = +1;
   _probe[45].value[7] = -1;
   _probe[45].value[8] = -1;
   _probe[45].value[9] = -1;
   _probe[45].value[10] = -1;
   _probe[45].value[11] = -1;

   // Dodecapole - unique-stereo - probe 17
   _probe[46].value[0] = +1;
   _probe[46].value[1] = +1;
   _probe[46].value[2] = +1;
   _probe[46].value[3] = +1;
   _probe[46].value[4] = +1;
   _probe[46].value[5] = -1;
   _probe[46].value[6] = -1;
   _probe[46].value[7] = -1;
   _probe[46].value[8] = +1;
   _probe[46].value[9] = -1;
   _probe[46].value[10] = -1;
   _probe[46].value[11] = -1;

   // Dodecapole - unique-stereo - probe 18
   _probe[47].value[0] = +1;
   _probe[47].value[1] = +1;
   _probe[47].value[2] = +1;
   _probe[47].value[3] = +1;
   _probe[47].value[4] = +1;
   _probe[47].value[5] = -1;
   _probe[47].value[6] = -1;
   _probe[47].value[7] = -1;
   _probe[47].value[8] = -1;
   _probe[47].value[9] = +1;
   _probe[47].value[10] = -1;
   _probe[47].value[11] = -1;
}



void
OBSpectrophore::_calculateProperties(OpenBabel::OBMol* mol)
{
   //
   // PROPERTY 1: ATOMIC PARTIAL CHARGES [0]
   //

   // CHI and ETA
   unsigned int dim(_nAtoms + 1);
   std::vector<double> CHI(dim);
   double** ETA = new double*[dim];
   double** ETA2 = new double*[dim];      // A copy for the electrophilicity later on
   for (unsigned int i = 0; i < dim; ++i)
   {
      ETA[i] = new double[dim];
      ETA2[i] = new double[dim];
   }
   double totalCharge(0.0);
   unsigned int i(0);
   unsigned int n;
   double hardness;
   double electronegativity;
   for (OpenBabel::OBMolAtomIter atom(mol); atom; ++atom)
   {
      n = (unsigned int) atom->GetAtomicNum();
      switch (n)
      {
         case 1:  // H
            hardness = 0.65971;
            electronegativity = 0.20606;
            break;
         case 3:  // Li
            hardness = 0.32966;
            electronegativity = 0.36237;
            break;
         case 5:  // B
            hardness = 0.32966;
            electronegativity = 0.36237;
            break;
         case 6:  // C
            hardness = 0.32966;
            electronegativity = 0.36237;
            break;
         case 7:  // N
            hardness = 0.34519;
            electronegativity = 0.49279;
            break;
         case 8:  // O
            hardness = 0.54428;
            electronegativity = 0.73013;
            break;
         case 9:  // F
            hardness = 0.72664;
            electronegativity = 0.72052;
            break;
         case 11: // Na
            hardness = 0.32966;
            electronegativity = 0.36237;
            break;
         case 12: // Mg
            hardness = 0.32966;
            electronegativity = 0.36237;
            break;
         case 14: // Si
            hardness = 0.32966;
            electronegativity = 0.36237;
            break;
         case 15: // P
            hardness = 0.32966;
            electronegativity = 0.36237;
            break;
         case 16: // S
            hardness = 0.20640;
            electronegativity = 0.62020;
            break;
         case 17: // Cl
            hardness = 0.32966;
            electronegativity = 0.36237;
            break;
         case 19: // K
            hardness = 0.32966;
            electronegativity = 0.36237;
            break;
         case 20: // Ca
            hardness = 0.32966;
            electronegativity = 0.36237;
            break;
         case 26: // Fe
            hardness = 0.32966;
            electronegativity = 0.36237;
            break;
         case 29: // Cu
            hardness = 0.32966;
            electronegativity = 0.36237;
            break;
         case 30: // Zn
            hardness = 0.32966;
            electronegativity = 0.36237;
            break;
         case 35: // Br
            hardness = 0.54554;
            electronegativity = 0.70052;
            break;
         case 53: // I
            hardness = 0.30664;
            electronegativity = 0.68052;
            break;
         default:
            hardness = 0.65971;
            electronegativity = 0.20606;
            break;
      }

      CHI[i] = -electronegativity;
      ETA[i][i] = 2.0 * hardness;

      // Adjust the total molecular charge
      totalCharge += atom->GetFormalCharge();

      // Increment
      ++i;
   }

   // Complete CHI
   CHI[_nAtoms] = totalCharge;

   // Complete ETA
   double d;
   for (unsigned int r = 0; r < _nAtoms; ++r)
   {
      for (unsigned int c = r + 1; c < _nAtoms; ++c)
      {
         d =  (_oricoor[r][0] - _oricoor[c][0]) * (_oricoor[r][0] - _oricoor[c][0]);
         d += (_oricoor[r][1] - _oricoor[c][1]) * (_oricoor[r][1] - _oricoor[c][1]);
         d += (_oricoor[r][2] - _oricoor[c][2]) * (_oricoor[r][2] - _oricoor[c][2]);
         ETA[r][c] = 0.529176 / sqrt(d);     // 0.529176: Angstrom to au
         ETA[c][r] = ETA[r][c];
      }
   }
   for (unsigned int i = 0; i < dim; ++i)
   {
      ETA[i][_nAtoms] = -1.0;
      ETA[_nAtoms][i] = +1.0;
   }
   ETA[_nAtoms][_nAtoms] = 0.0;

   // Make ETA2 a copy of ETA
   for (unsigned int r(0); r < dim; ++r)
   {
      for (unsigned int c(0); c < dim; ++c)
      {
         ETA2[r][c] = ETA[r][c];
      }
   }

   // Solve the matrix equation
   _solveMatrix(ETA, &(CHI[0]), dim);    // CHI will contain the values

   // Add values to property matrix
   for (unsigned int i = 0; i < _nAtoms; ++i)
   {
      _property[i][0] = CHI[i];
   }
   double CHIeq2 = CHI[dim - 1] * CHI[dim - 1]; // For electrophilicity later on

   // Clean
   for (unsigned int i = 0; i < dim; ++i)
   {
      delete[] ETA[i];
      ETA[i] = NULL;
   }
   delete[] ETA;
   ETA = NULL;



   //
   // PROPERTY 2: ATOMIC LIPOPHILICITY VALUES [1]
   //
   unsigned int a(0);
   for (OpenBabel::OBMolAtomIter atom(mol); atom; ++atom)
   {
      n = (unsigned int) atom->GetAtomicNum();
      switch (n)
      {
         case 1:		// H
            if (atom->GetExplicitDegree())
            {
               if (atom->IsNonPolarHydrogen())
               // Non-polar H
               {
                  _property[a][1] = -0.018;
               }
               else
               // Polar H
               {
                  _property[a][1] = -0.374;
               }
            }
            else
            {
               _property[a][1] = -0.175;
            }
            break;
        case 6:		// C
            _property[a][1] = 0.271;
            break;
        case 7:		// N
            _property[a][1] = -0.137;
            break;
        case 8:		// O
            _property[a][1] = -0.321;
            break;
        case 9:		// F
            _property[a][1] = 0.217;
            break;
        case 16:	// S
            _property[a][1] = 0.385;
            break;
        case 17:	// Cl
            _property[a][1] = 0.632;
            break;
        case 35:	// Br
            _property[a][1] = 0.815;
            break;
        case 53:	// I
            _property[a][1] = 0.198;
            break;
        default:	// The rest
            _property[a][1] = -0.175;
            break;
        }

        // Increment
        ++a;
    }



   //
   // PROPERTY 3: ATOMIC SHAPE DEVIATIONS [2]
   //
   for (unsigned int a(0); a < _nAtoms; ++a)
   {
      for (unsigned int c(0); c < 3; ++c)
      {
         _coor[a][c] = _oricoor[a][c];
      }
   }
   double COG[3]; // Center of geometry
   for (unsigned int c(0); c < 3; ++c)
   {
      COG[c] = _coor[0][c];
      for (unsigned int a(1); a < _nAtoms; ++a)
      {
         COG[c] += _coor[a][c];
      }
      COG[c] /= _nAtoms;
   }

   // Shift molecules to COG and calculate individual distances
   std::vector<double> distance(_nAtoms);
   double averageDistance(0.0);
   for (unsigned int a(0); a < _nAtoms; ++a)
   {
      distance[a] = 0.0;
      for (unsigned int c(0); c < 3; ++c)
      {
         _coor[a][c] -= COG[c];
         distance[a] += _coor[a][c] * _coor[a][c];
      }
      distance[a] = sqrt(distance[a]);
      averageDistance += distance[a];
   }
   averageDistance /= _nAtoms;

   // Set relative to average distance
   for (unsigned int a(0); a < _nAtoms; ++a)
   {
      distance[a] -= averageDistance;
      distance[a] *= averageDistance;
   }

   // Add property
   for (unsigned int a(0); a < _nAtoms; ++a)
   {
      _property[a][2] = distance[a];
   }



   //
   // PROPERTY 4: ATOMIC ELECTROPHILICITY [3]
   //
   // CHI and ETA2
   for (unsigned int i(0); i < dim; ++i)
   {
      CHI[i] = 1.0;
   }
   CHI[_nAtoms] = 0.0;

   // Complete ETA2
   for (unsigned int i = 0; i < dim; ++i)
   {
      ETA2[i][_nAtoms] = 0.0;
      ETA2[_nAtoms][i] = 1.0;
   }
   ETA2[_nAtoms][_nAtoms] = -1.0;

   // Solve the matrix equation
   _solveMatrix(ETA2, &(CHI[0]), dim);    // CHI will contain the values

   // Add values to property matrix
   for (unsigned int i = 0; i < _nAtoms; ++i)
   {
      _property[i][3] = CHI[i] * CHIeq2;
   }

   // Clean
   for (unsigned int i = 0; i < dim; ++i)
   {
      delete[] ETA2[i];
      ETA2[i] = NULL;
   }
   delete[] ETA2;
   ETA2 = NULL;

   //
   // Return
   //
   return;
}



void
OBSpectrophore::_solveMatrix(double** A, double* B, unsigned int dim)
{
   std::vector<int> temp(dim);
   _luDecompose(A, temp, dim);
   _luSolve(A, temp, B, dim);
}



void
OBSpectrophore::_luDecompose(double** A, std::vector<int>& I, unsigned int dim)
{
   unsigned int i, j, k, kMax, iMax;
   std::vector<double> vScales(dim, 0);
   double maxVal = 0, dummy = 0;
   double * pRowi = NULL;

   // first find the highest pivot element in each row and store it for implicit scaling
   for (i = 0; i < dim; ++i)
   {
      maxVal = 0.0;
      for (j = 0; j < dim; ++j)
      {
         if ((dummy=fabs(A[i][j])) > maxVal)
         {
            maxVal = dummy;
         }
      }
      if (maxVal == 0)
      {
         std::cerr << "OBSpectrophore: Warning singular matrix..." << std::endl;
      }

      vScales[i] = 1.0 / maxVal;
   }

   std::vector<double> colJ(dim); // variable to store local copy of column

   // loop over columns
   for (j = 0; j < dim; ++j)
   {
      // make a local copy of column j
      for (i = 0; i < dim; ++i) colJ[i] = A[i][j];
      for (i = 0; i < dim; ++i)
      {
         pRowi = A[i];
         dummy = pRowi[j];
         kMax = i < j ? i : j;
         for (k = 0; k < kMax; ++k) dummy -= pRowi[k] * colJ[k];
         colJ[i] = dummy;
         pRowi[j] = colJ[i];
      }

      // search largest pivot element
      maxVal = 0.0;
      iMax = j;
      for (i = j + 1; i < dim; ++i)
      {
         if ((dummy = fabs(colJ[i]) * vScales[i]) >= maxVal)
         {
            maxVal = dummy;
            iMax = i;
         }
      }

      // check if we need to interchange rows
      if (j != iMax) // if current column index is not the maximal row index we need to interchange
      {
         // std::cerr << "Swap rows: " << iMax << " <-> " << j << std::endl;
         _swapRows(A, iMax, j, dim);
         vScales[iMax] = vScales[j];
      }
      // store row index in I
      I[j] = iMax;

      // finally divide by the pivot element
      if (j != dim - 1)
      {
         dummy = 1.0 / A[j][j]; // A.GetValueAt(j,j);
         for (i = j + 1; i < dim; ++i) A[i][j] *= dummy;
      }


   } // next column

   return;
}



void
OBSpectrophore::_luSolve(double** A, std::vector<int>& I, double* B, unsigned int dim)
{
   for (unsigned int i = 0; i < dim; ++i) _swapRows(B, i, I[i]);

   // forward substitution pass
   for (unsigned int k = 0; k < dim; ++k)
   {
      for (unsigned int i = k+1; i < dim; ++i)
      {
         B[i] -= A[i][k] * B[k];
      }
   }

   // do the backsubstitution
   for (int i = dim - 1; i >= 0; --i)
   {
      B[i] /= A[i][i];
      for (int k = 0; k < i; ++k)
      {
         B[k] -= A[k][i] * B[i];
      }
   }

   return;
}



void
OBSpectrophore::_swapRows(double** _pMatrix, unsigned int i, unsigned int j, unsigned int nCols)
{
   double dummy;
   for (unsigned int k = 0; k < nCols; ++k)         // loop over all columns
   {
      dummy = _pMatrix[i][k];
      _pMatrix[i][k] = _pMatrix[j][k];
      _pMatrix[j][k] = dummy;
   }
   return;
}



void
OBSpectrophore::_swapRows(double* _pMatrix, unsigned int i, unsigned int j)
{
   double dummy;
   dummy = _pMatrix[i];
   _pMatrix[i] = _pMatrix[j];
   _pMatrix[j] = dummy;
   return;
}

} // end namespace OpenBabel

//! \file spectrophore.cpp
//! \brief Spectrophore(TM) calculator. Implementation of OBSpectrophore.

