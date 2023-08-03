/**********************************************************************
Copyright (C) 2005 by Chris Morley

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
#ifndef OB_KINETICS_H
#define OB_KINETICS_H

#include <openbabel/generic.h>

namespace OpenBabel
{

const unsigned RateData = 55555;
const unsigned ThermoData = 55556;

/// \class OBRateData kinetics.h <openbabel/kinetics.h>
/// \brief Holds rate constant data for OBReaction
/**The data is in a form used by CHEMKIN, at least the original
   free versions developed by Sandia Labs. Cantera is a GPL chemical
   kinetics program with similar capabilities.
   Used by chemkin format and cml reaction format
**/


class OBRateData : public OBGenericData
{
protected:
  double Rates[3];
  double LoRates[3];
  double TroeParams[4];
  std::map<std::string,double> Efficiencies;
public:
  virtual OBGenericData* Clone(OBBase* parent) const{return new OBRateData(*this);}
  enum rate_type {A, n, E};
  enum reaction_type {ARRHENIUS=55555, LINDERMANN, TROE, SRI, THREEBODY};
  reaction_type ReactionType;
  OBRateData():OBGenericData("Rate data", RateData)
  {
    Rates[0]=Rates[1]=Rates[2]=0;
    LoRates[0]=LoRates[1]=LoRates[2]=0;
    TroeParams[0]=TroeParams[1]=TroeParams[2]=TroeParams[3]=0;
    ReactionType = ARRHENIUS;
  }

  double GetRate(rate_type n) const
  {
    return Rates[n];
  }

  void SetRate(rate_type n, const double val)
  {
    if(n<3)
      Rates[n] = val;
  }

  double GetLoRate(rate_type n) const
  {
    return LoRates[n];
  }

  void SetLoRate(rate_type n, const double val)
  {
    if(n<3)
      LoRates[n] = val;
  }

  double GetTroeParam(int n) const
  {
    return TroeParams[n];
  }

  void SetTroeParams(int n, const double val)
  {
    if(n<4)
      TroeParams[n] = val;
  }

  void SetEfficiency(std::string id, double Eff)
  {
    Efficiencies[id] = Eff;
  }

  double GetEfficiency(std::string id)
  {
    return Efficiencies[id]; //will be 0 if not found
  }

  bool GetNextEff(std::string& id, double& Eff)
  {
    //Supply id empty to begin, then id is the*last* id
    std::map<std::string, double>::iterator itr;
    if(id.empty())
      itr = Efficiencies.begin();
    else
    {
      itr = Efficiencies.find(id);
      if(itr!=Efficiencies.end())
        ++itr;
    }
    if(itr==Efficiencies.end())
      return false;
    id    = itr->first;
    Eff = itr->second;
    return true;
  }
};

//******************************************************************************
/// \class OBNasaThermoData kinetics.h <openbabel/kinetics.h>
/// \brief Thermodynamic data in old style NASA polynomial form for OBMol
/**This is a venerable data format used to describe specific heats, enthalpies
   and entropies, particularly in the gas phase and at high temperatures.
   There is a standard datafile with fixed format (for punched cards!) which
   can be read and written to this OBMol extension using the thermo format. It
   is also used in chemkin format and in cmlreact format
   For a brief description of the meaning of the coefficients see
   http://www.me.berkeley.edu/gri_mech/data/nasa_plnm.html
   The first 7 coefficients are for the high temperature range MidT to HiT;
   and the second 7 are for the low temperature range LoT to MidT
   Note that there is a more modern NASA polynomial with more terms, which
   is not supported here.
**/
class OBNasaThermoData : public OBGenericData
{
protected:
  double Coeffs[14];
  double LoT, MidT, HiT;
  char phase;
public:
  OBNasaThermoData(): LoT(300),MidT(1000),HiT(3000),phase('G')
  {	_type = ThermoData;	_attr = "Nasa thermo data";}

  virtual OBGenericData* Clone(OBBase* parent) const{return new OBNasaThermoData(*this);}

  double GetCoeff(unsigned n) const
  {
    return Coeffs[n];
  }

  void SetCoeff(unsigned n, const double val)
  {
    if(n<14)
      Coeffs[n] = val;
  }
  double GetLoT() const {return LoT;}
  double GetMidT() const {return MidT;}
  double GetHiT() const {return HiT;}
  void SetLoT(double val){LoT=val;}
  void SetMidT(double val){MidT=val;}
  void SetHiT(double val){HiT=val;}

  char GetPhase() const {return phase;}
  void SetPhase(char ph){phase=ph;}
};

} //namespace OpenBabel

#endif //OB_KINETICS_H

//! \file kinetics.h
//! \brief OBRateData and OBNasaThermoData classes
