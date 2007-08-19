  /**********************************************************************
Copyright (C) 2005 by Chris Morley
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#include "openbabel/babelconfig.h"
#include <string>
#include <iomanip>
#include <map>
#include <set>
#include "openbabel/mol.h"
#include "openbabel/oberror.h"
#include "openbabel/obconversion.h"
#include "openbabel/reaction.h"
#include "openbabel/kinetics.h"
#include "openbabel/obmolecformat.h"

using namespace std;
using std::tr1::shared_ptr;

namespace OpenBabel
{
  
class ChemKinFormat : public OBFormat
{
public:
  ChemKinFormat()
  {
    OBConversion::RegisterFormat("ck",this);
    OBConversion::RegisterOptionParam("s", this); //no params
    OBConversion::RegisterOptionParam("t", this);
  }

  virtual const char* Description()
  {
      return
"ChemKin format\n"
"Input Options e.g. -ai\n"
"f <file> File with standard thermo data: default therm.dat\n"
"z Use standard thermo only\n\n"
"Output options e.g. -xs\n"
"s Simple output: reactions only\n"
"t Do not include species thermo data\n\n";
  };

  virtual const char* TargetClassDescription()
  {
      return OBReaction::ClassDescription();
  };

  const type_info& GetType()
  {
    return typeid(OBReaction*);
  };
private:
  shared_ptr<OBMol> CheckSpecies(string& name, string& ln, bool inreactions);
  bool ReadThermo(OBConversion* pConv);
  bool ReadStdThermo(const string& datafilename);
  OBFormat* GetThermoFormat();
  bool CheckAllMolsHaveThermo();
  bool WriteReactionLine(OBReaction* pReact, OBConversion* pConv);
  bool WriteHeader(OBConversion* pConv);
private:
  typedef map<string,shared_ptr<OBMol> > MolMap;
  typedef set<shared_ptr<OBMol> > MolSet;
  MolMap IMols; //used on input
  MolSet OMols; //used on output
  stringstream ss;//used on output

  ////////////////////////////////////////////////////
  /// The "API" interface functions
  virtual bool WriteMolecule(OBBase* pReact, OBConversion* pConv);

  ////////////////////////////////////////////////////
  /// The "Convert" interface functions
  virtual bool ReadChemObject(OBConversion* pConv);
  
  virtual bool WriteChemObject(OBConversion* pConv)
  {
    OBBase* pOb=pConv->GetChemObject();
    OBReaction* pReact = dynamic_cast<OBReaction*>(pOb);
    bool ret=false;
    if(pReact!=NULL)
    {
      ret=WriteMolecule(pReact,pConv);

      std::string auditMsg = "OpenBabel::Write reaction ";
      std::string description(Description());
            auditMsg += description.substr( 0, description.find('\n') );
            obErrorLog.ThrowError(__FUNCTION__,
                                  auditMsg,
                                  obAuditMsg);

      vector<shared_ptr<OBMol> >::iterator itr;
    }
    delete pOb;
    return ret;
  }
};

//Make an instance of the format class
ChemKinFormat theChemKinFormat;

/////////////////////////////////////////////////////////////////
bool ChemKinFormat::ReadChemObject(OBConversion* pConv)
{
  //Called once. AddChemObject is called when each OBReaction is completed. 
  obErrorLog.ThrowError(__FUNCTION__, "OpenBabel:read reactions with Chemkin Format", obAuditMsg);
  istream &ifs = *pConv->GetInStream();

  OBReaction* pReact = NULL;
  OBRateData* pRD = NULL;
  IMols.clear();
  bool ThermoLoaded=false;
  double EUnitsFactor = 1.0;
  double AUnitsFactor = 1.0;

  vector<string>::iterator itr;
  unsigned i;
  string ln;
  bool doingspecies=false;
  bool inreactions=false;
  bool lookforend=false;

  //Special species name
  if(IMols.empty())
  {
    //Three body reactions and similar are stored in OBReaction with M as reactant and product
    shared_ptr<OBMol> sp(new OBMol);
    sp.get()->SetTitle("M");
    IMols["M"] = sp;
  }

//Loop through whole file. Save OBReaction objects when found. Save the last after this loop
  while(ifs.good())
  {
    if(!getline(ifs,ln))
      continue;
    if(Trim(ln).empty() || ln[0]=='!')
      continue;
    //Line before comment is made uppercase
    transform(ln.begin(), find(ln.begin(),ln.end(),'!'), ln.begin(),toupper);
    
    vector<string> toks;
    string::size_type eqpos, commentpos;
    eqpos = ln.find('=');
    commentpos = ln.find('!');
    if(eqpos==string::npos || (commentpos!=string::npos && commentpos < eqpos) )
    {
      //line is not a reaction
      tokenize(toks, ln, " \t\n\r/\\");
      if(lookforend)
      {
        if(find(toks.begin(),toks.end(),"END")!=toks.end())
          lookforend = false;
        continue;
      }

      if(pRD && toks[0]=="LOW")
      {
        if(pRD->ReactionType != OBRateData::TROE)
          pRD->ReactionType = OBRateData::LINDERMANN;
        unsigned n;
        for(n=0;n<3;++n)
        {
          double val = atof(toks[n+1].c_str());
          if(n==0)
            val /= pow(AUnitsFactor,pReact->NumReactants());
          else if(n==2)
            val /= EUnitsFactor;
          pRD->SetLoRate((OBRateData::rate_type)n, val );
        }
      }
      else if(pRD && toks[0]=="TROE")
      {
        pRD->ReactionType = OBRateData::TROE;
        for(i=0;i<4;++i)
          pRD->SetTroeParams(i, atof(toks[i+1].c_str()));
      }
      else if(toks[0]=="DUPLICATE")
      {}
      else if(pReact && toks[0]=="TS")
      {
        //Defines the molecule which is a transition state for a reaction
        //This is not a ChemKin keyword. Used for Mesmer.
        pReact->SetTransitionState(CheckSpecies(toks[1], ln, inreactions));
      }

      else if(toks[0]=="THERMO")
      {
        //Read following data using Thermo format
        if(!pConv->IsOption("z",OBConversion::INOPTIONS))
        {
          pConv->AddOption("e", OBConversion::INOPTIONS); //stops on END
          ReadThermo(pConv);
          pConv->RemoveOption("e", OBConversion::INOPTIONS);
        }
        else
        {
          lookforend=true;
        }
      }
      else if(doingspecies || toks[0]=="SPECIES" || toks[0]=="SPEC")
      {
        //Add all species to IMols
        vector<string>::iterator itr;
        itr=toks.begin();
        if(!doingspecies) ++itr; //ignore "SPECIES"
        doingspecies=true;
        for(;itr!=toks.end();++itr)
        {
          if(*itr=="END")
          {
            doingspecies=false;
            break;
          }
          shared_ptr<OBMol> sp(new OBMol);
          sp.get()->SetTitle(*itr);
          IMols[*itr] = sp;
        }
        continue;
      }
      else if(toks[0]=="REACTIONS" || toks[0]=="REAC")
      {
        inreactions=true;

        //Need to check here whether thermo data has been input and if not
        //load it from therm.dat
        if(!CheckAllMolsHaveThermo())
        {
          string stdthermo("therm.dat"); //default
          const char* pstd = pConv->IsOption("f",OBConversion::INOPTIONS);
          if(pstd)
            stdthermo=pstd;
          if(!ReadStdThermo(stdthermo))
            return false;
        }
        //TODO units on this line
        string EKeywords[6] ={"CAL/MOLE","KCAL/MOLE","JOULES/MOLE","KJOULES/MOLE","KELVINS","EVOLTS"};
        double EFactor[6]   ={   1.0    ,   0.001  ,    4.1816    ,   0.041816   ,   1.98  , 0.0};
        double AvFactor = 6.023E23;
        //MOLECULES
        for(i=1;i<toks.size();++i)
        {
          for(int j=0;j<6;++j)
            if(toks[i]==EKeywords[j])
              EUnitsFactor = EFactor[j];
          if(toks[i]=="MOLECULES")
            AUnitsFactor = AvFactor;
        }
        continue;
      }
      else if(inreactions && toks[0]=="END") //of reactions
      {
        IMols.clear();
        break;
      }		
      else if(inreactions && pRD)
      {
        //3-body efficiencies
        for(i=0;i<toks.size();++i)
          pRD->SetEfficiency(toks[i++], atof(toks[i].c_str()));
      }
      //other unrecognized lines before"REACTIONS" are ignored
    }
    else 
    {
      /* Line is a reaction
      A single reaction line is input here on each call
      and the input stream left ready for the next.
      Lines like the following are handled
      Label A + B => C + D 1E-12 0.2 2300 !comment
      H2 = 2H 1e-8 0 112000 comment: has A n E
      2H + M => H2 + M 1e-16 comment: has A only
      Label A+B = C+D comment: has no rates
      */
      //Save the previous reaction
      if(pReact)
        if(!pConv->AddChemObject(pReact->DoTransformations(pConv->GetOptions(OBConversion::GENOPTIONS))))
          return false;

      pReact = new OBReaction;
      pRD = new OBRateData; //to store rate constant data
      pReact->SetData(pRD);

      int n=0;
      shared_ptr<OBMol> sp;
    
      if(ln[eqpos-1]=='<' || ln[eqpos+1]!='>')
      {
        //Reaction is reversible: contains <=> or =
        pReact->SetReversible();
      }

      //Replace each (+M) by M
      string::size_type pos;
      while((pos = ln.find("(+M)")) != string::npos)
        ln.replace(pos, 4, " M  ");
      
      //Do reactants
      const char delim[] = " \t\n\r+<>=/";
      tokenize(toks, ln.substr(0,eqpos), delim);
      for(itr=toks.begin();itr!=toks.end();++itr)
      {		
        //See if first token is a label
        if(itr==toks.begin()&& isalpha((*itr)[0])  
          && (pConv->IsOption("l",OBConversion::INOPTIONS)
          || !IMols.count(*itr))) // not a species name?		
        {
          pReact->SetTitle(*itr);
          continue;
        }
        if(isalpha((*itr)[0]))
        {
          if(*itr == "M")
            pRD->ReactionType = OBRateData::THREEBODY;
          sp = CheckSpecies(*itr, ln, inreactions);
          if(!sp.get())
            return false;
          pReact->AddReactant(sp);
          continue;
        }
        else
        {
          if(isalpha((*itr)[1]))
          {
            //species multiplier (single digit)
            unsigned mult = atoi(itr->c_str());
            sp = CheckSpecies(itr->substr(1), ln, inreactions);
            if(!sp.get())
              return false;
            for(i=0;i<mult;++i)
              pReact->AddReactant(sp);
            continue;
          }
          else
          {
            obErrorLog.ThrowError(__FUNCTION__, 
              "In " + ln  +
              "\nThe species multiplier has currently to be a single digit integer",
              obError);
            return false; //incorrect multiplier
          }
        }
      }

      //Do products
      //Remove the +
      if(commentpos!=string::npos)
        tokenize(toks, ln.substr(eqpos+1,commentpos-eqpos-1), delim);
      else
        tokenize(toks, ln.substr(eqpos+1), delim);

      for(itr=toks.begin();itr!=toks.end();++itr)
      {		
        if(isalpha((*itr)[0]))
        {
          sp = CheckSpecies(*itr, ln, inreactions);
          if(!sp.get())
            return false;
          pReact->AddProduct(sp);
        }
        else
        {
          unsigned len = itr->size();
          if(len>1 && isalpha((*itr)[1]) && len>2 && !isalpha((*itr)[2])) //2OH 2E ~ 2 21 214 2E+10  
          {
            //species multiplier (single digit)
            unsigned mult = atoi(itr->c_str());
            sp = CheckSpecies(itr->substr(1), ln, inreactions);
            if(!sp.get())
              return false;
            for(i=0;i<mult;++i)
              pReact->AddProduct(sp);
          }
          else
          {
            //Rate parameters	
            //Correct for e.g. 1.23E+07 being split after E
            char lastch = (*itr)[itr->size()-1];
            string num(*itr);
            if(lastch=='E' || lastch=='e')
              num += *(++itr); 

            stringstream ss(num);
            double val;
            ss >> val;
            if(n==0)
              val /= pow(AUnitsFactor,pReact->NumReactants());
            else if(n==2)
              val /= EUnitsFactor;
            pRD->SetRate((OBRateData::rate_type)n++, val);
            if(!ss)
            {
              //not numeric: put into comment (better than doing nothing)
              pReact->SetComment(*itr);
              break;
            }
          }
        }
      }
      if(commentpos!=string::npos)
        pReact->SetComment(ln.substr(commentpos+1));
    }
  }

  bool ret;
  if(pReact && pReact->NumReactants() + pReact->NumProducts()>0)
    ret = pConv->AddChemObject(pReact->DoTransformations(pConv->GetOptions(OBConversion::GENOPTIONS)));
  else
    ret = pConv->AddChemObject(NULL);
  return ret;
}
///////////////////////////////////////////////////////////////

shared_ptr<OBMol> ChemKinFormat::CheckSpecies(string& name, string& ln, bool inreactions)
{
  MolMap::iterator mapitr = IMols.find(name);
  if(mapitr==IMols.end())
  {
    //unknown species
    if(inreactions)
    {
      obErrorLog.ThrowError(__FUNCTION__,
        name + " not recognized as a species in\n" + ln, obError);
      shared_ptr<OBMol> sp;
      return sp; //empty
    }
    else
    {
      // There was no REACTIONS section in input file and probably no SPECIES section.
      // Unknown species that appear in a reaction can be made here with just a name.
      shared_ptr<OBMol> sp(new OBMol);
      sp->SetTitle(name.c_str());
      return sp;
    }
  }
  else
    return mapitr->second;
}


//////////////////////////////////////////////////////////////////
bool ChemKinFormat::ReadThermo(OBConversion* pConv)
{
  /*	Reads molecule using thermoformat.
      Finds mol in IMols with same name
       and combines the one with OBNasaThermoData with it.
      Continue with all molecules.
      Construct index if pIndex!=NULL.
  */
  OBFormat* pThermFormat = OBConversion::FindFormat("therm");
  if(!pThermFormat)
  {
    obErrorLog.ThrowError(__FUNCTION__,
    "Thermo format needed but not available", obError);
    return false;
  }
  else
  {
    pConv->SetInFormat(pThermFormat);
    pConv->AddOption("e", OBConversion::INOPTIONS); //stops on END

    OBMol thmol;
    while(pConv->Read(&thmol))
    {
      MolMap::iterator mapitr = IMols.find(thmol.GetTitle());
      if(mapitr!=IMols.end())
      {
        shared_ptr<OBMol> psnewmol(OBMoleculeFormat::MakeCombinedMolecule(mapitr->second.get(),&thmol));
        IMols.erase(mapitr);
        IMols[thmol.GetTitle()] = psnewmol;
      }
      thmol.Clear();
    }
    pConv->SetInFormat(this);
  }
  pConv->RemoveOption("e", OBConversion::INOPTIONS);
  return true;
}

/////////////////////////////////////////////////////////////////
bool ChemKinFormat::ReadStdThermo(const string& datafilename)
{
  OBMoleculeFormat::NameIndexType index;
  OBFormat* pThermFormat = GetThermoFormat();

  //Get the index of std thermo file, which may involve it being prepared
  if(!pThermFormat || !OBMoleculeFormat::ReadNameIndex(index, datafilename, pThermFormat))
    return false;
  
  string missing; // list of molecules which do not have thermodata
  OBConversion StdThermConv;
  ifstream stdthermo;
  OpenDatafile(stdthermo, datafilename);
  if(!stdthermo) 
  {
    obErrorLog.ThrowError(__FUNCTION__,
    datafilename + " was not found", obError);
    return false;
  }
  StdThermConv.SetInFormat(pThermFormat);
  StdThermConv.SetInStream(&stdthermo);

  MolMap::iterator mapitr;
  for(mapitr=IMols.begin();mapitr!=IMols.end();++mapitr)
  {
    //Look up each molecules's name in index, move the the returned seek position,
    //read the molecule and combine it with the one in Imols
    OBMoleculeFormat::NameIndexType::iterator itr = index.find(mapitr->first);
    if(itr!=index.end())
    {
      OBMol thmol;
      stdthermo.seekg(itr->second);
      StdThermConv.Read(&thmol);
      shared_ptr<OBMol> psnewmol(OBMoleculeFormat::MakeCombinedMolecule(mapitr->second.get(),&thmol));
      IMols[thmol.GetTitle()] = psnewmol;
    }	
    else
      if(mapitr->first!="M")
        missing += mapitr->first + ',';
  }
  if(!missing.empty())
  {
    obErrorLog.ThrowError(__FUNCTION__,
    datafilename + " does not contain thermodata for " + missing, obError);
    return false;
  }
  return true;
}

//////////////////////////////////////////////////////////
bool ChemKinFormat::CheckAllMolsHaveThermo()
{
  MolMap::iterator mapitr;
  for(mapitr=IMols.begin();mapitr!=IMols.end();++mapitr)
  {
    if(!mapitr->second->GetData(ThermoData) && mapitr->first!="M")
      return false;
  }
  return true;
}

/////////////////////////////////////////////////////////////////
bool ChemKinFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
  //It's really a reaction, not a molecule. Called separately for each reaction.
  //Cast output object to the class type need, i.e. OBReaction
  OBReaction* pReact = dynamic_cast<OBReaction*>(pOb);
  if(pReact==NULL)
      return false;

  //Read in reaction, store mols in OMols, write reaction to stringstream ss.
  if(pConv->GetOutputIndex()==1)
  {
    OMols.clear();
    ss.str("");
  }
  if(!WriteReactionLine(pReact, pConv))
    return false;

  //At end, construct ELEMENTS and SPECIES and output to ofs followed by ss
  if(pConv->IsLast())
  {
    ostream& ofs = *pConv->GetOutStream();
    if(!pConv->IsOption("s"))
    {
      if(!WriteHeader(pConv))
        return false;
      ofs << "REACTIONS\n";
    }
    ofs  << ss.rdbuf() << endl;     
    if(!pConv->IsOption("s"))
      ofs << "END" << endl;
  }
  return true;
}

//////////////////////////////////////////////////////////////////
bool ChemKinFormat::WriteHeader(OBConversion* pConv)
{
  ostream& ofs = *pConv->GetOutStream();

  set<string> elements;
  vector<string> species;
  MolSet::iterator itr;
  for(itr= ++OMols.begin();itr!=OMols.end();++itr) //not first species "M"
  {
    const char* title = (*itr)->GetTitle();
    if(strcmp(title, "M"))
      species.push_back(title);
    FOR_ATOMS_OF_MOL(atom, itr->get())
      elements.insert(etab.GetSymbol(atom->GetAtomicNum()));
  }
  if(!elements.empty())
  {
    ofs << "ELEMENTS\n";
    copy(elements.begin(),elements.end(), ostream_iterator<string>(ofs," "));
    ofs << "\nEND\n";
  }
  else
    obErrorLog.ThrowError(__FUNCTION__, "No element data available", obWarning);

  ofs << "SPECIES\n";
  vector<string>::iterator sitr;
  int maxlen=0;
  for(sitr= species.begin();sitr!=species.end();++sitr)
    if(sitr->size()>maxlen) maxlen = sitr->size();

  int n=0;
  for(sitr=species.begin();sitr!=species.end();++sitr, ++n)
  {
    if(n > 80 / maxlen)
    {
      ofs << '\n';
      n=0;
    }
    ofs << setw(maxlen+1) << *sitr;
  }
  ofs << "\nEND\n";

  if(!pConv->IsOption("t"))
  {
    OBFormat* pFormat = OBConversion::FindFormat("therm");
    if(!pFormat)
    {
      obErrorLog.ThrowError(__FUNCTION__,
      "Thermo format needed but not available", obError);
      return false;
    }
    else
    {
      stringstream thermss;
      thermss << "THERMO ALL\n";
      thermss << "   300.000  1000.000  5000.000\n";
      OBConversion ConvThermo(*pConv);
      ConvThermo.SetOutFormat(pFormat);
      ConvThermo.SetOutStream(&thermss);
      int ntherm=0;
      for(itr= OMols.begin();itr!=OMols.end();++itr)
      {
        const char* title = (*itr)->GetTitle();
        if(strcmp(title, "M"))
          if(ConvThermo.Write(itr->get())) 
            ++ntherm;
      }

      thermss << "END\n";
      if(ntherm)
        ofs << thermss.str(); //but don't output unless there was some thermo data
    }
  }
  return true;
}

//////////////////////////////////////////////////////////////////
bool ChemKinFormat::WriteReactionLine(OBReaction* pReact, OBConversion* pConv)
{
  //Get rate data so that we know what kind of reaction it is
  OBRateData* pRD = static_cast<OBRateData*>(pReact->GetData(RateData));

  ss << pReact->GetTitle() << '\t';

  if(!pRD)obErrorLog.ThrowError(__FUNCTION__, "Reaction " + pReact->GetTitle()
     + " has no rate data", obWarning);

  int i;
  for(i=0;i<pReact->NumReactants();++i)
  {
    shared_ptr<OBMol> psMol = pReact->GetReactant(i);
//    if(strcasecmp(psMol->GetTitle(),"M"))
    OMols.insert(psMol);
    if (i)
      ss << " + ";
    ss << setw(3) << left << psMol->GetTitle();
  }

  /*
  3-body
  H + H + M <=> H2 + M  May have efficiencies
  Lindemann
  O + CO (+M) <=> CO2 (+M) Has LOW/ and may have efficiencies. Troe[0]=0
  Troe
  H + CH3 (+M) <=> CH4 (+M) Has LOW/ and TROE/ and may have efficiencies
  SRI
  */

  string mstring;
  if(pRD)
  {
    switch(pRD->ReactionType)
    {
    case OBRateData::TROE:
    case OBRateData::SRI:
    case OBRateData::LINDERMANN:
      mstring = " (+M) ";
    }
  }
  if(mstring.empty() && pReact->NumReactants()<3)
    ss << "     ";

  ss << mstring;

  if(pReact->IsReversible())
    ss << "\t <=> \t";
  else
    ss << "\t => \t";

  for(i=0;i<pReact->NumProducts();++i)
  {
    shared_ptr<OBMol> psMol = pReact->GetProduct(i);
    if(strcasecmp(psMol->GetTitle(),"M"))
      OMols.insert(psMol);
    if (i)
      ss << " + ";
    ss << setw(3) << left << psMol->GetTitle();
  }
  if(mstring.empty() && pReact->NumProducts()<3)
    ss << "     ";

  ss << mstring;

  if(pRD)
  {
    ss << " \t" << scientific << setprecision(3) << pRD->GetRate(OBRateData::A) << ' ' 
      << fixed << pRD->GetRate(OBRateData::n)	<< ' ' 
      << setprecision(1) << pRD->GetRate(OBRateData::E) 
      << " \t" << pReact->GetComment() << endl;

    switch(pRD->ReactionType)
    {
    case OBRateData::TROE:
      ss << "\tTROE / " << pRD->GetTroeParam(0) << ' ' << pRD->GetTroeParam(1) 
        << ' ' << pRD->GetTroeParam(2);
      if(pRD->GetTroeParam(3))
        ss << ' ' <<pRD->GetTroeParam(3);
      ss << '/' << endl;
      //fallthrough
    case OBRateData::LINDERMANN:
      ss << "\tLOW / " << scientific << setprecision(3) << pRD->GetLoRate(OBRateData::A) << ' '
        << fixed << pRD->GetLoRate(OBRateData::n) << ' '
        << setprecision(1) << pRD->GetLoRate(OBRateData::E) << '/' << endl;
      //fallthrough
    case OBRateData::THREEBODY:
      string id;
      double eff;
      int neffs=0;
      while(pRD->GetNextEff(id,eff))
      {
        if(!neffs) ss << '\t';
        ss << id << "/ " << eff << "/ ";
        ++neffs;
      }
      if(neffs)
        ss << endl;
    }
  }
  return true;
}

OBFormat* ChemKinFormat::GetThermoFormat()
{
  OBFormat* pThermFormat = OBConversion::FindFormat("therm");
  if(!pThermFormat)
  {
    obErrorLog.ThrowError(__FUNCTION__,
    "Thermo format needed but not available", obError);
    return NULL;
  }
  return pThermFormat;
}
}
/*
LINDEMANN FALLOFF FORM
This treatment applies if no specific falloff parameters are given. 
At pressures intermediate to the high and low pressure limits, 
the rate constant is given by the Lindemann formula:

                    k_inf
         k = ----------------
             1  +  k_inf/k_o[M]

In cases where no high pressure limit rate constant parameters are given
(i.e., the collider M as a reactant is not in parenthesis),
the reaction is in the low pressure limit.

TROE FALLOFF FORM
A more refined treatment of pressure effects than Lindemann is employed 
using the TROE parameters. The falloff parameter F_cent for a unimolecular 
reaction is calculated from the values of a, b, c, and d by the formula of Troe

      F_cent  =  (1-a) exp(-T/b)  +  a exp(-T/c)  +  exp(-d/T)

which gives the temperature dependence of F_cent, the factor by which 
the rate constant of a given unimolecular reaction at temperature T and 
reduced pressure P_r = k_o[M]/k_inf of 1.0 is less than the value k_inf/2 which
it would have if unimolecular reactions behaved according to the Lindemann formula.

The broadening factor F, which is 1 for the Lindemann case where no parameters 
for F_cent are provided, is computed from F_cent by

                                log F_cent
      log F = ---------------------------------------------
              1 + [(log P_r + C)/(N - 0.14{log P_r + C})]^2


      with N = 0.75 - 1.27log F_cent  and  C = -0.4 - 0.67log F_cent.

The rate coefficient, k, is then given by multiplying the Lindemann formula by F.
a,b,c,d = a, T***, T*, T**

See also http://gems.mines.edu/~reactionxml/Fall-off2.pdf
*/
/*
TODO
Make case independent
Isotopes like D in ELEMENTS list
ELEM and END is optional
*/