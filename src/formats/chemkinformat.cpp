  /**********************************************************************
Copyright (C) 2005-2007 by Chris Morley

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
#include "openbabel/babelconfig.h"

#include <string>
#include <iomanip>
#include <map>
#include <set>
#include <iterator>
#include <locale>

#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/elements.h>
#include <openbabel/obiter.h>

#include "openbabel/oberror.h"
#include "openbabel/obconversion.h"
#include "openbabel/reaction.h"
#include "openbabel/kinetics.h"
#include "openbabel/obmolecformat.h"

#include <cstdlib>

using namespace std;

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
    Init();
  }

  virtual const char* Description()
  {
      return
"ChemKin format\n"
"Read Options e.g. -aL\n"
" f <file> File with standard thermo data: default therm.dat\n"
" z Use standard thermo only\n"
" L Reactions have labels (Usually optional)\n"
"\n"
"Write Options e.g. -xs\n"
" s Simple output: reactions only\n"
" t Do not include species thermo data\n"
" 0 Omit reactions with zero rates\n"
"\n";
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
  void              Init();
  ///\return -1 eof or error; +1 reactionline found; 0 otherwise
  int               ReadLine(istream& ifs);
  bool              ReadHeader(istream& ifs, OBConversion* pConv);
  bool              ParseReactionLine(OBReaction* pReact, OBConversion* pConv);
  bool              ReadReactionQualifierLines(istream& ifs, OBReaction* pReact);
  obsharedptr<OBMol> CheckSpecies(string& name, string& ln, bool MustBeKnown);
  bool              ReadThermo(OBConversion* pConv);
  bool              ReadStdThermo(const string& datafilename);
  OBFormat*         GetThermoFormat();
  bool              CheckAllMolsHaveThermo();
  bool              WriteReactionLine(OBReaction* pReact, OBConversion* pConv);
  bool              WriteHeader(OBConversion* pConv);
private:
  typedef map<string,obsharedptr<OBMol> > MolMap;
  typedef set<obsharedptr<OBMol> > MolSet;
  //used on input
  MolMap IMols;
  string ln;
  bool SpeciesListed;
  double AUnitsFactor, EUnitsFactor;
  string comment;
  //used on output
  MolSet OMols;
  stringstream ss;

  ////////////////////////////////////////////////////
  /// The "API" interface functions
  virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
  virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

  ////////////////////////////////////////////////////
  /// The "Convert" interface functions
  virtual bool ReadChemObject(OBConversion* pConv)
  {
    std::string auditMsg = "OpenBabel::Read ChemKinFormat";
    std::string description(Description());
    auditMsg += description.substr(0,description.find('\n'));
    obErrorLog.ThrowError(__FUNCTION__,
              auditMsg,
              obAuditMsg);
    //Makes a new OBReaction
    OBReaction* pReact = new OBReaction;
    bool ret=ReadMolecule(pReact,pConv); //call the "API" read function

    if(ret) //Do transformation and return molecule
      return pConv->AddChemObject(pReact->DoTransformations(pConv->GetOptions(OBConversion::GENOPTIONS),pConv))!=0;
    else
        pConv->AddChemObject(NULL);
    return false;
  }

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
    }
    delete pOb;
    return ret;
  }
};

//Make an instance of the format class
ChemKinFormat theChemKinFormat;

/////////////////////////////////////////////////////////////////
bool ChemKinFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{
  /*Badly name function. It handles OBReaction objects.

  This format can be used with Chemkin files or with simple lists of
  reactions.

  In the simple case:
  Reversible reactions use = or <=> Irreversible reactions use =>
      A + B = C
      D     => E + F
  Rate parameters A, n and E can follow, separated by spaces.
      A+B=C 1.2E12 1.0 5000
  The species can have names containing spaces
      first reactant + second reactant => first product + second product
  The reactions can optionally have a label
      label first reactant => first product
  In ambiguous cases like
      label A + B = C
  the program assumes that "label A" is the first reactant. To force the
  correct interetation, use the -al option "Reactions have labels".
  Even in the simple case, each reaction can be followed by lines with
  additional information in the CheKin style, like low pressure rates,
  third body efficiencies, etc.

  On the first call, the function reads any header sections, ELEMENTS,
  SPECIES,  THERMO that are present and populates the IMols list of
  molecules. It then reads the first reaction and returns true.
  Subsequent calls also read a single reaction. A return of false
  indicates there are no more reactions.
  */
  OBReaction* pReact = dynamic_cast<OBReaction*>(pOb);
  if(!pReact)
    return false;

  istream& ifs = *pConv->GetInStream();

  if(pConv->IsFirstInput())
  {
    Init();
    if(!ReadHeader(ifs, pConv))
    {
      obErrorLog.ThrowError(__FUNCTION__, "Unexpected end of file or file reading error", obError);
      return false;
    }
  }

  if(!ifs                                      //possibly EOF
   || !ReadLine(ifs)                           //not a reaction line
   || !ParseReactionLine(pReact, pConv)        //faulty parse
   || !ReadReactionQualifierLines(ifs, pReact))//END or erroneous line found
   return false;

  //return true if reaction has either reactants or products
  return pReact->NumReactants() + pReact->NumProducts()>0;
}

/////////////////////////////////////////////////
void ChemKinFormat::Init()
{    //initialize the member variables used during input
    ln.clear();
    AUnitsFactor = 1.0;
    EUnitsFactor = 1.0;
    SpeciesListed=false;
    IMols.clear();
    //Special species name
    obsharedptr<OBMol> sp(new OBMol);
    sp.get()->SetTitle("M");
    IMols["M"] = sp;
}

/////////////////////////////////////////////////
// Uses the line from the member variable ln, probably from previous call to
// ReadMolecule, if it is not empty. Otherwise read a line and extract the comment.
// Returns -1 on eof or error; +1 if line is a reaction; 0 otherwise.
int ChemKinFormat::ReadLine(istream& ifs )
{
  while(ln.empty())
  {
    if(!getline(ifs,ln))
      return -1;
    //discard lines that are empty or contain just a comment
    if(Trim(ln).empty() || ln[0]=='!')
      ln.clear();
    comment.clear();
  }
  string::size_type eqpos, commentpos;
  commentpos = ln.find('!');
  //Extract and remove comment
  if(commentpos!=string::npos)
  {
    comment = ln.substr(commentpos+1);
    ln.erase(commentpos);
  }

  eqpos = ln.find('=');
  //eof may have been set, but we need ReadMolecule() to be called again to process this line
  ifs.clear();
  return eqpos==string::npos? 0 : +1;
}

//////////////////////////////////////////////////////
bool ChemKinFormat::ReadHeader(istream& ifs, OBConversion* pConv)
{
  bool doingspecies=false;
  //loop for each line until a reaction line is found
  while(ifs)
  {
    if(int ret=ReadLine(ifs)!=0)
      return ret>0; //reaction line found: there may have been no header

    vector<string> toks;
    tokenize(toks, ln, " \t\n\r/\\");
    ln.clear(); //have to clear line when it has been dealt with

    if(doingspecies || !strcasecmp(toks[0].c_str(),"SPECIES") || !strcasecmp(toks[0].c_str(),"SPEC"))
    {
      SpeciesListed = true; //Means that molecules in reactions must have been specified in SPECIES

      vector<string>::iterator itr;
      itr=toks.begin();
      if(!doingspecies) ++itr; //ignore "SPECIES"
      doingspecies=true;
      for(;itr!=toks.end();++itr)
      {
        if(*itr=="END" || *itr=="end")
        {
          doingspecies=false;
          break;
        }
        //Add all species to IMols
        obsharedptr<OBMol> sp(new OBMol);
        sp.get()->SetTitle(*itr);
        IMols[*itr] = sp;
      }
    }

    else if(!strcasecmp(toks[0].c_str(),"THERMO"))
    {
      //Read following data using Thermo format
      if(!pConv->IsOption("z",OBConversion::INOPTIONS))
      {
        pConv->AddOption("e", OBConversion::INOPTIONS); //stops on END
        ReadThermo(pConv);
        pConv->RemoveOption("e", OBConversion::INOPTIONS);
      }
    }

    else if(!strcasecmp(toks[0].c_str(),"REACTIONS") || !strcasecmp(toks[0].c_str(),"REAC"))
    {
      //Units may be specified on this line
      string EKeywords[6] ={"CAL/MOLE","KCAL/MOLE","JOULES/MOLE","KJOULES/MOLE","KELVINS","EVOLTS"};
      double EFactor[6]   ={   1.0    ,   0.001  ,    4.1816    ,   0.041816   ,   1.98  , 0.0};
      double AvFactor = 6.023E23;

      for (unsigned int i=1; i<toks.size(); ++i)
      {
        for(int j=0;j<6;++j)
          if(!strcasecmp(toks[i].c_str(), EKeywords[j].c_str()))
            EUnitsFactor = EFactor[j];
        if(!strcasecmp(toks[i].c_str(),"MOLECULES"))
          AUnitsFactor = AvFactor;
      }

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

    }

  // Anthing not in a SPECIES or THERMO section is ignored.
  // This includes the ELEMENTS section
  }
  return false; //failed file read
}

//////////////////////////////////////////////
bool ChemKinFormat::ParseReactionLine(OBReaction* pReact, OBConversion* pConv)
{
  /* Line is a reaction
  Lines like the following are handled
  Label A + B => C + D 1E-12 0.2 2300 !comment
  H2 = 2H 1e-8 0 112000 comment: has A n E
  2H + M => H2 + M 1e-16 comment: has A only
  Label A+B = C+D comment: has no rates
  */
  OBRateData* pRD = new OBRateData; //to store rate constant data. Attach only if rate data found

  int n=0;
  obsharedptr<OBMol> sp;

  string::size_type eqpos = ln.find('=');

  bool r1=false, r2=false;
  //Ensure divider between reactants and products is just '='
  if(eqpos>0 && ln[eqpos-1]=='<')
  {
    ln[eqpos-1] = ' ';
    r1=true;
  }
  if(eqpos < ln.size()-1 && ln[eqpos+1]=='>')
  {
    ln[eqpos+1] = ' ';
    r2=true;
  }
  if(r1 || !r2)
  {
    //Reaction is reversible: contains <=> or =
    pReact->SetReversible();
  }

  //Replace each (+M) by M
  string::size_type pos;
  while((pos = ln.find("(+M)")) != string::npos)
    ln.replace(pos, 4, " +M ");
  while((pos = ln.find("(+m)")) != string::npos)
    ln.replace(pos, 4, " +M ");

  //Do reactants
  vector<string> toks;
  vector<string>::iterator itr;
  string temp = ln.substr(0, eqpos);
  tokenize(toks, temp, "+");
  //(ln is cleared later)

  for(itr=toks.begin();itr!=toks.end();++itr)
  {
    Trim(*itr);
    if(itr==toks.begin())
    {
      /*First token can contain a label, and reactant can contain spaces
        label reactant1      +  1
        first reactant       +  case 2
        label first reactant +  case 3
        reactant1            +  case 4
        label                =  case 5
        reactant1            =  case 6
        (1 and 2) and (5 and 6) are ambiguous if -al option not set. Assume 2 or 6 and issue a warning
      */
      vector<string> firstr;
      tokenize(firstr, *itr, " \t");
      if(isalpha(firstr[0][0]))
      {
        //Starts with letter, so could be a label. Further tests...
        if(pConv->IsOption("L",OBConversion::INOPTIONS)//this option mandates a label
          || firstr.size()>2                           //case 3 above
          || (SpeciesListed && !IMols.count(*itr)))    // there is a species list and it is not a species name
        {
          pReact->SetTitle(firstr[0]);             //Add label to OBReaction
          Trim(toks[0].erase(0, firstr[0].size()));//Remove label leaving only first reactant
        }

        //Ambiguous cases
        else if(firstr.size()==2 || (firstr.size()==1 && toks.size()==1))
          obErrorLog.ThrowError(__FUNCTION__,
            "In " + ln +
            "\nThe string " + firstr[0] + " has been assumed NOT to be a label\n"
            "If it should be, use the -aL option which mandates labels on reactions.\n"
            "A species missing from the SPECIES section, if one is used, can also give this error",
            obWarning);
      }
    }

    if(isalpha((*itr)[0]))
    {
      if(*itr == "m")
        *itr="M";
      if(*itr == "M")
        pRD->ReactionType = OBRateData::THREEBODY;
      sp = CheckSpecies(*itr, ln, SpeciesListed);
      if(!sp.get())
      {
        ln.clear();
        return false;
      }
      pReact->AddReactant(sp);
      continue;
    }
    else
    {
      if(isalpha((*itr)[1]))
      {
        //species multiplier (single digit)
        unsigned mult = atoi(itr->c_str());
        string temp = itr->substr(1);
        sp = CheckSpecies(temp, ln, SpeciesListed);
        if(!sp.get())
        {
          ln.clear();
          return false;
        }
        for (unsigned int i=0; i<mult; ++i)
          pReact->AddReactant(sp);
        continue;
      }
      else
      {
        obErrorLog.ThrowError(__FUNCTION__,
          "In " + ln  +
          "\nThe species multiplier must be a single digit integer",
          obError);
        ln.clear();
        return false; //incorrect multiplier
      }
    }
  }

  //Do products
  temp = ln.substr(eqpos+1);
  tokenize(toks, temp, "+");
  if(toks.size()>0)
  {
    /*
      product1
      2product1
      first product
      last product   7.7 8.8
      product   7.7 8.8
      product   7.7 8.8E
      12 0
    */
    //Combine tokens erroneously split at + in 8.8E+12
    for(int i = toks.size()-1;i>0 && isdigit(toks[i][0]);--i)//break when starts with letter
    {
      char lastch = toks[i-1][toks[i-1].size()-1];
      if(lastch=='E' || lastch=='e')
      {
        toks[i-1] += toks[i];
        toks.pop_back();
      }
    }

    //Split the last token and separate the rate parameters
    vector<string> lastp;
    tokenize(lastp, toks[toks.size()-1], " \t");
    toks[toks.size()-1].clear();
    bool HasRateData=false;
    n=0;
    for(itr=lastp.begin();itr!=lastp.end();++itr)
    {
      Trim(*itr);
      //copy species names and species names with multiplier (single digit) back into orig token
      unsigned len = itr->size();

      //Separate products like 2C2H6 2O 2ETR 2H2 from rate params like 2E25 2E+12 -1 .00
      if(isalpha((*itr)[0]) ||            //to be a product: 1st char is a letter or...
        (len>=2 && isalpha((*itr)[1])     //both second char (it there is one) is a letter and
        && (len<=2 || isalpha((*itr)[2])  //  third char(if there is one) is a letter
        || toupper((*itr)[1])!='E')))     //  or the second char is not 'E'
      {
        toks[toks.size()-1] += ' ' + *itr;
        continue;
      }

      //Read in rate parameters
      stringstream ss(*itr);
      locale cLocale("C");
      ss.imbue(cLocale);

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
      HasRateData=true;
    }
    //Rate parameters were specified, so OBReaction needs to have the OBRateData attached
    if(HasRateData)
      pReact->SetData(pRD);
    else if(SpeciesListed) //a true ChemKin file
      obErrorLog.ThrowError(__FUNCTION__,
            "In " + ln + "\nNo rate data found.", obWarning);

    //Read in product species
    for(itr=toks.begin();itr!=toks.end();++itr)
    {
      Trim(*itr);
      if(isalpha((*itr)[0]))
      {
        if(*itr == "m")
          *itr="M";

        sp = CheckSpecies(*itr, ln, SpeciesListed);
        if(!sp.get())
        {
          ln.clear();
          return false;
        }
        pReact->AddProduct(sp);
      }
      else
      {
        if(itr->size()>1 && isalpha((*itr)[1]))
        {
          //species multiplier (single digit)
          unsigned mult = atoi(itr->c_str());
          string temp = itr->substr(1);
          sp = CheckSpecies(temp, ln, SpeciesListed);
          if(!sp.get())
          {
            ln.clear();
            return false;
          }
          for (unsigned int j=0; j<mult; ++j)
            pReact->AddProduct(sp);
        }
        else
          obErrorLog.ThrowError(__FUNCTION__,
            "In " + ln + "\nError in products or rate parameters.", obError);
      }
    }
  }
  pReact->SetComment(comment);
  ln.clear();
  return true;
}

/////////////////////////////////////////////
bool ChemKinFormat::ReadReactionQualifierLines(istream& ifs, OBReaction* pReact)
{
  OBRateData* pRD = (OBRateData*)pReact->GetData("Rate data");

  while(ifs)
  {
    if(int ret=ReadLine(ifs)!=0)
      return ret>0; //The next reaction has been found

    vector<string> toks;
    tokenize(toks, ln, " \t\n\r/\\");
    ln.clear(); //have to clear line when it has been dealt with

    if(pRD && !strcasecmp(toks[0].c_str(),"LOW"))
    {
      if(pRD->ReactionType != OBRateData::TROE)
        pRD->ReactionType = OBRateData::LINDERMANN;
      unsigned n;
      for(n=0;n<3;++n)
      {
        double val = atof(toks[n+1].c_str());
        if(n==0)
          val /= pow(AUnitsFactor, pReact->NumReactants());
        else if(n==2)
          val /= EUnitsFactor;
        pRD->SetLoRate((OBRateData::rate_type)n, val );
      }
    }
    else if(pRD && !strcasecmp(toks[0].c_str(),"TROE"))
    {
      pRD->ReactionType = OBRateData::TROE;
      for(int i=0;i<4;++i)
        pRD->SetTroeParams(i, atof(toks[i+1].c_str()));
    }

    else if(!strcasecmp(toks[0].c_str(),"DUPLICATE"))
    {}

    else if(pReact && !strcasecmp(toks[0].c_str(),"TS"))
    {
      //Defines the molecule which is a transition state for a reaction
      //This is not a ChemKin keyword. Used for Mesmer.
      pReact->SetTransitionState(CheckSpecies(toks[1], ln, SpeciesListed));
    }

    else if(pRD && strcasecmp(toks[0].c_str(),"END") && toks.size()%2==0)
    {
      //not "END". Has an even number of tokens.
      //3-body efficiencies
      for(int i=0;i<toks.size()-1;++i)//also incremented in body to retrieve id,val pairs
      {
        string sp(toks[i++]);
        pRD->SetEfficiency(sp, atof(toks[i].c_str()));
      }
    }
  }
  return (bool)ifs;
}

///////////////////////////////////////////////////////////////
obsharedptr<OBMol> ChemKinFormat::CheckSpecies(string& name, string& ln, bool MustBeKnown)
{
  MolMap::iterator mapitr = IMols.find(name);
  if(mapitr==IMols.end())
  {
    //unknown species
    if(MustBeKnown)
    {
      obErrorLog.ThrowError(__FUNCTION__,
        name + " not recognized as a species in\n" + ln, obError);
      obsharedptr<OBMol> sp;
      return sp; //empty
    }
    else
    {
      // There was no REACTIONS section in input file and probably no SPECIES section.
      // Unknown species that appear in a reaction can be made here with just a name.
      obsharedptr<OBMol> sp(new OBMol);
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
        obsharedptr<OBMol> psnewmol(OBMoleculeFormat::MakeCombinedMolecule(mapitr->second.get(),&thmol));
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
      obsharedptr<OBMol> psnewmol(OBMoleculeFormat::MakeCombinedMolecule(mapitr->second.get(),&thmol));
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

  WriteReactionLine(pReact, pConv);

  //At end, construct ELEMENTS and SPECIES and output to ofs followed by ss
  if(pConv->IsLast())
  {
    ostream& ofs = *pConv->GetOutStream();
    if(!pConv->IsOption("s")) //Simple output option - reactions only
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
  for(itr= OMols.begin();itr!=OMols.end();++itr)
  {
    const char* title = (*itr)->GetTitle();
    if(strcmp(title, "M"))
      species.push_back(title);
    FOR_ATOMS_OF_MOL(atom, itr->get())
      elements.insert(OBElements::GetSymbol(atom->GetAtomicNum()));
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
  unsigned int maxlen=0;
  for(sitr= species.begin();sitr!=species.end();++sitr)
    if(sitr->size()>maxlen) maxlen = sitr->size();

  unsigned int n=0;
  for(sitr=species.begin();sitr!=species.end();++sitr, ++n)
  {
    if(maxlen>0 && n > 80 / maxlen)
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

  //If -0 option set, omit reactions with zero rates. However, number of reactions converted remains the same.
  if(pConv->IsOption("0"))
    if(!pRD || pRD->GetRate(OBRateData::A)==0.0)
      return false;

  ss << pReact->GetTitle() << '\t';

  if(!pRD && !pConv->IsOption("s"))
    obErrorLog.ThrowError(__FUNCTION__, "Reaction " + pReact->GetTitle()
     + " has no rate data", obWarning);

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

  int i;
  for(i=0;i<pReact->NumReactants();++i)
  {
    obsharedptr<OBMol> psMol = pReact->GetReactant(i);
//    if(strcasecmp(psMol->GetTitle(),"M"))
    OMols.insert(psMol);

    //If reactant has no title use its formula
    if(*psMol->GetTitle()=='\0')
      psMol->SetTitle(psMol->GetSpacedFormula(1,"").c_str());

    //write species name but, if M, only if (+M) is not going to be output
    if(mstring.empty() || strcasecmp(psMol->GetTitle(),"M"))
    {
      if (i)
        ss << " + ";
      ss << setw(3) << left << psMol->GetTitle();
    }
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

  if(mstring.empty() && pReact->NumReactants()<3)
    ss << "     ";

  ss << mstring;

  if(pReact->IsReversible())
    ss << "\t <=> \t";
  else
    ss << "\t => \t";

  for(i=0;i<pReact->NumProducts();++i)
  {
    obsharedptr<OBMol> psMol = pReact->GetProduct(i);
    if(strcasecmp(psMol->GetTitle(),"M"))
      OMols.insert(psMol);

    //If product has no title use its formula
    if(*psMol->GetTitle()=='\0')
      psMol->SetTitle(psMol->GetSpacedFormula(1,"").c_str());

    //write species name but, if M, only if (+M) is not going to be output
    if(mstring.empty() || strcasecmp(psMol->GetTitle(),"M"))
    {
      if (i)
        ss << " + ";
      ss << setw(3) << left << psMol->GetTitle();
    }
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
      ss << "\tTROE / " << setprecision(3) << pRD->GetTroeParam(0) << ' '
        << pRD->GetTroeParam(1) << ' ' << pRD->GetTroeParam(2);
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
        ss << id << "/ " << setprecision(2) << eff << "/ ";
        ++neffs;
      }
      if(neffs)
        ss << endl;
    }
  }
  else //simple option
    ss << pReact->GetComment() << endl;

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
} //namespace
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
@todo
Make case independent
Isotopes like D in ELEMENTS list
ELEM and END is optional
*/
