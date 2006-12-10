/**********************************************************************
Copyright (C) 2005,2006 Chris Morley

Based on the IUPAC InChI reference software, which is distributed
under the GNU LGPL:
Copyright (C) 2005 The International Union of Pure and Applied Chemistry
IUPAC International Chemical Identifier (InChI) (contact:secretariat@iupac.org)

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/obmolecformat.h>

#include "inchi_api.h"
#include <sstream>
#include <set>
#include <vector>

using namespace std;
namespace OpenBabel
{
extern string GetInChI(istream& is);

class InChIFormat : public OBMoleculeFormat
{
public:
  InChIFormat()
  {
    OBConversion::RegisterFormat("inchi",this);
    OBConversion::RegisterOptionParam("n", this, 0, OBConversion::INOPTIONS);
    OBConversion::RegisterOptionParam("t", this);
  }

  virtual const char* Description()
  {
      return 
"InChI format\n \
IUPAC/NIST molecular identifier\n \
Write options, e.g. -xat \n \
 t add molecule name\n \
 a output auxilliary information\n \
 u output only unique molecules\n \
 U output only unique molecules and sort them\n \
 e compare first molecule to others\n \
 w don't warn on undef stereo or charge rearrangement\n\n \
Note that when reading an InChI the stereochemical info is currently ignored.\n \
Input options, e.g. -at\n \
 n molecule name follows InChI on same line \n \
 a add InChI string to molecule name \n\n"
  ;
  };

  virtual const char* SpecificationURL()
  { return "http://www.iupac.org/inchi/";};

  virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
  virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
  virtual int  SkipObjects(int n, OBConversion* pConv);

  static char   CompareInchi(const char* Inchi1, const char* Inchi2);
  static string InChIErrorMessage(const char ch);

private:
  OBAtom* GetCommonAtom(OBBond* pb1, OBBond* pb2);

  ///Compare std::strings with embedded numbers so that 
  // "a6b" (or "a06b") is less than "a15b"
  // and "CH4" is less than "C2H6"
  // and "CH4" is less than "ClH" (hydrogen chloride)
  struct InchiLess
    : public binary_function<const string&, const string&, bool>
  {
    bool operator()(const string& s1, const string& s2) const
    {
      string::const_iterator p1, p2;
      p1=s1.begin(); p2=s2.begin();
      while( p1!=s1.end() && p2!=s2.end() )
      {
        if(iscntrl(*p1) || iscntrl(*p2) || isspace(*p1) || isspace(*p2))
          return false; //stop comparison at whitespace. Identical up to here 
        int n1=-1,n2=-1;
        if(isdigit(*p1))
          {
            n1 = atoi(&*p1);
            //skip over number
            while(p1!=s1.end() && isdigit(*p1++)); --p1;
          }
        if(isdigit(*p2))
          {
            n2 = atoi(&*p2);
            while(p2!=s2.end() && isdigit(*p2++)); --p2;
          }
        if(n1<0 && n2 < 0)
          {
            //neither numbers
            if(*p1 != *p2)
        return *p1 < *p2;
          }
        else if(n1>=0 && n2>0)
          {
            //both numbers
            if(n1!=n2)
        return n1 < n2;
          }
        else if(n1>0)
          return islower(*p2)!=0;
        else if(n2>0)
          return !islower(*p1);

        ++p1; ++p2; // iterate
      } // while loop
      return false; //identical
    }
  };

  typedef	set<string, InchiLess> nSet;
  nSet allInchi;
  string firstInchi;
  string firstID;
};

//Make an instance of the format class
InChIFormat theInChIFormat;

/////////////////////////////////////////////////////////////////
bool InChIFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{
  OBMol* pmol = dynamic_cast<OBMol*>(pOb);
  if(pmol==NULL) return false;
  istream &ifs = *pConv->GetInStream();

  //Extract InChI from input stream even if it is split
  string inchi;
  do
  {
    inchi = GetInChI(ifs);
    if(inchi.empty())
      return false; //eof
  }while(inchi.size()<9); //ignore empty "InChI=" or "InChI=1/"

  if(pConv->IsFirstInput())
    clog << "Note that any stereochemical info in the InChI is not currently read." << endl;

  //Set up input struct
  inchi_InputINCHI inInchi;
  inInchi.szOptions="";
  inInchi.szInChI = const_cast<char*>( inchi.c_str()); //***KLUDGE***
  inchi_OutputStruct out;

  //Call the conversion routine in InChI code
  int ret = GetStructFromINCHI( &inInchi, &out );

  if (ret!=inchi_Ret_OKAY)
  {
    obErrorLog.ThrowError("InChI code", out.szMessage, obWarning);
  }

  //Read name if requested e.g InChI=1/CH4/h1H4 methane 
  //OR InChI=1/CH4/h1H4 "First alkane"  Quote can be any punct char and
  //uses upto the end of the line if second quote is not found
  if(pConv->IsOption("n",OBConversion::INOPTIONS))
  {
    string name;
    if(getline(ifs, name))
      pmol->SetTitle(Trim(name));
  }

  //Option to add InChI text to title
  if(pConv->IsOption("a",OBConversion::INOPTIONS))
  {
    string title(pmol->GetTitle());
    title += ' ' + inchi;
    pmol->SetTitle(title);
  }

  //Translate the returned structure into OBMol
  pmol->SetDimension(0);
  pmol->BeginModify();
  int i;
  //Make all atoms first because we need pointers to later ones
  for(i=0;i<out.num_atoms;++i)
  {
    pmol->NewAtom();
  }
  for(i=0;i<out.num_atoms;++i)
  {
    OBAtom* patom = pmol->GetAtom(i+1); //index starts at 1
    inchi_Atom* piat = &out.atom[i];
    int iso=0;
    patom->SetAtomicNum(etab.GetAtomicNum(piat->elname,iso));
    patom->SetIsotope(iso);
    if(piat->isotopic_mass)
      patom->SetIsotope(piat->isotopic_mass - ISOTOPIC_SHIFT_FLAG +
          (int)(isotab.GetExactMass(patom->GetAtomicNum())+0.5));

    patom->SetSpinMultiplicity(piat->radical);
    patom->SetFormalCharge(piat->charge);
    patom->SetVector(piat->x,piat->y,piat->z);

    int j;
    for(j=0;j<piat->num_bonds;++j)
    {
      pmol->AddBond(i+1, piat->neighbor[j]+1, piat->bond_type[j]);	
    }
  }
  
  //***TODO 0D stereo, implicit H isotopes
  /*for(i=0;i<out.num_stereo0D;++i)
  {
    inchi_Stereo0D stereo = out.stereo0D;
    switch(stereo.type)
    {
    case INCHI_StereoType_DoubleBond:
      break;
    case INCHI_StereoType_Tetrahedral:
      OBAtom* patom = pmol->GetAtom(stereo.central_atom+1);
      if(!patom)
        return false;

      break;
    case INCHI_StereoType_Allene:
      default:
      obErrorLog.ThrowError("InChI code", "Unsupported stereo type has been ignored.", obWarning);
    }
  }
  
  
  vector<inchi_Stereo0D> stereoVec;
  //Tetrahedral stereo
  OBAtom* patom;
  vector<OBNodeBase*>::iterator itr;
  if(mol.IsChiral())
  {
    for(patom = mol.BeginAtom(itr);patom;patom = mol.NextAtom(itr))
    {
      if(patom->IsChiral())
      {
        inchi_Stereo0D stereo;
        stereo.central_atom = patom->GetIdx()-1;
        stereo.type = INCHI_StereoType_Tetrahedral;
        int i=0;
        vector<OBEdgeBase*>::iterator itr;
        OBBond *pbond;
        for (pbond = patom->BeginBond(itr);pbond;pbond = patom->NextBond(itr),++i)
          stereo.neighbor[i] = pbond->GetNbrAtomIdx(patom)-1;
        //if(i!=4)...error

        stereo.parity = INCHI_PARITY_UNKNOWN;
        if(patom->IsPositiveStereo() || patom->IsClockwise())
          stereo.parity = INCHI_PARITY_ODD;
        if(patom->IsNegativeStereo() || patom->IsAntiClockwise())
          stereo.parity = INCHI_PARITY_EVEN;
        stereoVec.push_back(stereo);
      }
    }
  }
  */
  pmol->EndModify();
  FreeStructFromINCHI( &out );
  return true;
}

int InChIFormat::SkipObjects(int n, OBConversion* pConv)
{
  istream& ifs = *pConv->GetInStream();
  string inchi;
  while(ifs.good() && n)
  {
    inchi = GetInChI(ifs);
    if(inchi.size()>=8)//ignore empty "InChI=" or "InChI=1/"
      --n;
  }
  return ifs.good() ? 1 : -1;
}

/////////////////////////////////////////////////////////////////
bool InChIFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
  OBMol* pmol = dynamic_cast<OBMol*>(pOb);
  if(pmol==NULL) return false;
  
  //Use a copy of the molecule - we will be modifying it by adding H 
  OBMol mol = *pmol;

  stringstream molID;
  if(strlen(mol.GetTitle())==0)
    molID << '#' << pConv->GetOutputIndex() << ' ';
  else
    molID << mol.GetTitle() << ' ';
  if(pConv->GetOutputIndex()==1)
    firstID=molID.str();
  
  mol.AddHydrogens(false,false); //so stereo works

  inchi_Input inp;
  memset(&inp,0,sizeof(inchi_Input));
  bool Is0D=true;
  OBAtom* patom;
  vector<inchi_Atom> inchiAtoms(mol.NumAtoms());
  vector<OBNodeBase*>::iterator itr;
  for (patom = mol.BeginAtom(itr);patom;patom = mol.NextAtom(itr))
  {
    //OB atom index starts at 1; inchi atom index starts at 0
    inchi_Atom& iat = inchiAtoms[patom->GetIdx()-1];
    memset(&iat,0,sizeof(inchi_Atom));
    iat.x = patom->GetX();
    iat.y = patom->GetY();
    iat.z = patom->GetZ();
    if(iat.x!=0 || iat.y!=0 || iat.z!=0)
      Is0D=false;
    
    int nbonds = 0;
    vector<OBEdgeBase*>::iterator itr;
    OBBond *pbond;
    for (pbond = patom->BeginBond(itr);pbond;pbond = patom->NextBond(itr))
    {
      // do each bond only once. Seems necessary to avoid problems with stereo
      if(pbond->GetNbrAtomIdx(patom)<patom->GetIdx()) continue;
      iat.neighbor[nbonds]      = pbond->GetNbrAtomIdx(patom)-1;
      int bo = pbond->GetBO();
      if(bo==5)
        bo=4;
      iat.bond_type[nbonds]     = bo;

      S_CHAR bs = INCHI_BOND_STEREO_NONE;

      if(pbond->IsWedge())
        bs = INCHI_BOND_STEREO_SINGLE_1UP;
      if(pbond->IsHash())
        bs = INCHI_BOND_STEREO_SINGLE_1DOWN;

      iat.bond_stereo[nbonds++] = bs;
      if(nbonds>MAXVAL)
      {
#ifdef HAVE_SSTREAM
        stringstream errorMsg;
#else
        strstream errorMsg;
#endif
        errorMsg << "Too many bonds to " << iat.elname << " atom" << endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
        return false;
      }
    }
  
    strcpy(iat.elname,etab.GetSymbol(patom->GetAtomicNum()));
    iat.num_bonds = nbonds;
    iat.num_iso_H[0] = -1; //Let inchi add implicit Hs
    if(patom->GetIsotope())
    {
      iat.isotopic_mass = ISOTOPIC_SHIFT_FLAG +
        patom->GetIsotope() - (int)(isotab.GetExactMass(patom->GetAtomicNum())+0.5);
    }
    else
      iat.isotopic_mass = 0 ;
    iat.radical = patom->GetSpinMultiplicity();
    //InChI doesn't recognize spin miltiplicity of 4 or 5 (as used in OB for CH and C atom)
    if(iat.radical>=4)
    {
      iat.radical=0;
      iat.num_iso_H[0] = 0; //no implicit hydrogens
    }
    iat.charge  = patom->GetFormalCharge();
  }
  
  inp.atom = &inchiAtoms[0];

  vector<inchi_Stereo0D> stereoVec;
  if(Is0D)
  {
    //Tetrahedral stereo
    mol.FindChiralCenters();
    OBAtom* patom;
    vector<OBNodeBase*>::iterator itr;
    if(mol.IsChiral())
    {
      for(patom = mol.BeginAtom(itr);patom;patom = mol.NextAtom(itr))
      {
        if(patom->IsChiral())
        {
          inchi_Stereo0D stereo;
          stereo.central_atom = patom->GetIdx()-1;
          stereo.type = INCHI_StereoType_Tetrahedral;
          int i=0;
          vector<OBEdgeBase*>::iterator itr;
          OBBond *pbond;
          for (pbond = patom->BeginBond(itr);pbond;pbond = patom->NextBond(itr),++i)
            stereo.neighbor[i] = pbond->GetNbrAtomIdx(patom)-1;
          //if(i!=4)...error

          stereo.parity = INCHI_PARITY_UNKNOWN;
          if(patom->IsPositiveStereo() || patom->IsClockwise())
            stereo.parity = INCHI_PARITY_ODD;
          if(patom->IsNegativeStereo() || patom->IsAntiClockwise())
            stereo.parity = INCHI_PARITY_EVEN;
          stereoVec.push_back(stereo);
        }
      }
    }
    
    //Double bond stereo
    //Currently does not handle cumulenes
    set<OBBond*> UpDown;
    set<OBBond*>::iterator uditr;
    OBBond *pbond;
    vector<OBEdgeBase*>::iterator bitr;
    for (pbond = mol.BeginBond(bitr);pbond;pbond = mol.NextBond(bitr))
    {
     if(pbond->IsUp() || pbond->IsDown())
       UpDown.insert(pbond);
    }
    if(!UpDown.empty())
    {
      //For all double bonds look for up/down bonds attached
      for (pbond = mol.BeginBond(bitr);pbond;pbond = mol.NextBond(bitr))
      {
        if(pbond->IsDouble())
        {
          OBAtom* patom, *pA=NULL, *pB=NULL;
          OBBond* pX, *pY;;
          for(uditr=UpDown.begin();uditr!=UpDown.end();++uditr)
          {
            patom = GetCommonAtom(pbond, *uditr);
            if(patom && !pA)
            {
              //first one in pA
              pA = patom;
              pX = *uditr;
            }
            else if(patom && pA)
            {
              //second one in pB
              pB = patom;
              pY = *uditr;
              break;
            }
          }
          if(pA && pB)
          {
            inchi_Stereo0D stereo;
            stereo.central_atom = NO_ATOM;
            stereo.type = INCHI_StereoType_DoubleBond;
            stereo.neighbor[0]= pX->GetNbrAtomIdx(pA)-1;
            stereo.neighbor[1] = pA->GetIdx()-1;
            stereo.neighbor[2] = pB->GetIdx()-1;
            stereo.neighbor[3]= pY->GetNbrAtomIdx(pB)-1;
            if((pX->IsUp() && pY->IsUp())||(pX->IsDown() && pY->IsDown()))
              stereo.parity = INCHI_PARITY_ODD;
            else
              stereo.parity = INCHI_PARITY_EVEN;
            stereoVec.push_back(stereo);
          }
        }
      }
    }
  }

//	*inp.szOptions = '\0';
  inp.num_atoms = mol.NumAtoms();
  inp.num_stereo0D = stereoVec.size();
  if(inp.num_stereo0D>0)
    
    inp.stereo0D = &stereoVec[0];

  inchi_Output inout;
  memset(&inout,0,sizeof(inchi_Output));

  int ret = GetINCHI(&inp, &inout);

  if(ret!=inchi_Ret_OKAY)
  {
    string mes(inout.szMessage);
    if(!(pConv->IsOption("w") 
         && (mes=="Omitted undefined stereo" || mes=="Charges were rearranged")))
    {
      obErrorLog.ThrowError("InChI code", molID.str() + ':' + mes, obWarning);
    }
    if(ret!=inchi_Ret_WARNING)
    {
      FreeINCHI(&inout);
      return false;
    }
  }
  
  string ostring = inout.szInChI;
  if(pConv->IsOption("t"))
  {
    ostring += ' ';
    ostring +=  mol.GetTitle();
  }
    
  ostream &ofs = *pConv->GetOutStream();

  if(pConv->IsOption("U"))
  {
    if(pConv->GetOutputIndex()==1)
      allInchi.clear();
    //Just add to set and don't output, except at the end
    allInchi.insert(ostring);
    
    if(pConv->IsLast())
    {
      nSet::iterator itr;
      for(itr=allInchi.begin();itr!=allInchi.end();++itr)
        ofs << *itr << endl;
    }
    return true;
  }
  else if(pConv->IsOption("u"))
  {
    if(pConv->GetOutputIndex()==1)
      allInchi.clear();
    if(!allInchi.insert(ostring).second)
      return true; //no output if already in set
  }

  ofs << ostring << endl;

  if (pConv->IsOption("a"))
    ofs << inout.szAuxInfo << endl;
  
  if(pConv->IsOption("e"))
  {
    if(pConv->GetOutputIndex()==1)
      firstInchi = inout.szInChI;
    else
    {
      ofs << "Molecules " << firstID << "and " << molID.str();
      ofs << InChIErrorMessage(CompareInchi(firstInchi.c_str(), inout.szInChI)) << endl; 
    }
  }

  FreeINCHI(&inout);
  return true;
}

//////////////////////////////////////////////////////////
char InChIFormat::CompareInchi(const char* Inchi1, const char* Inchi2)
{
  //Returns 0 if identical or an char identifying the layer where they first differed
  string s1(Inchi1), s2(Inchi2);

  //Remove anything after the end of the Inchi
  string::size_type pos;
  pos = s1.find_first_of(" \t\n");
  if(pos!=string::npos)
    s1.erase(pos);
  pos = s2.find_first_of(" \t\n");
  if(pos!=string::npos)
    s2.erase(pos);

  vector<string> layers1, layers2;
  tokenize(layers1,s1,"/\n");
  tokenize(layers2,s2,"/\n");
  int i;
  if(layers1.size()<layers2.size())
    layers1.swap(layers2); //layers1 is the longest
  
  for(i=1;i<layers2.size();++i)
  {
    if(layers1[i]!=layers2[i])
    {
      char ch = '+';
      if(i>1) //not formula layer
        ch=layers1[i][0];
      return ch;
    }
  }
  if(layers1.size()==layers2.size())	
    return 0;
  else
    return layers1[i][0];
}

string InChIFormat::InChIErrorMessage(const char ch)
{
  string s;
  switch (ch) 
  {
  case 0:
    s = " are identical";
    break;
  case '+':
    s = " have different formulae";
    break;
  case 'c':
    s = " have different connection tables";
    break;
  case 'h':
    s = " have different bond orders, or radical character";
    break;
  case 'q':
    s = " have different charges";
    break;
  case 'p':
    s = " have different numbers of attached protons";
    break;
  case 'b':
    s = " have different double bond stereochemistry";
    break;
  case 'm':
  case 't':
    s = " have different sp3 stereochemistry";
    break;
  case 'i':
    s = " have different isotopic composition";
    break;
  default:
    s = " are different";
  }
  return s;
}

OBAtom* InChIFormat::GetCommonAtom(OBBond* pb1, OBBond* pb2)
{
  OBAtom* pa1 = pb1->GetBeginAtom();
  if(pa1==pb2->GetBeginAtom() || pa1==pb2->GetEndAtom())
    return pa1;
  pa1 = pb1->GetEndAtom();
  if(pa1==pb2->GetBeginAtom() || pa1==pb2->GetEndAtom())
    return pa1;
  return NULL; //not adjacent bonds
}

//******************************************************
class InChICompareFormat : public OBMoleculeFormat
{
public:
  InChICompareFormat()
  {
      OBConversion::RegisterFormat("k",this);
  }

  virtual const char* Description() //required
  {
      return 
"Compares first molecule with others using InChI\n \
e.g. babel first.smi second.mol third.cml -ok\n \
Same as -oinchi -xet\n \
";
  };
  
  virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    pConv->AddOption("e",OBConversion::OUTOPTIONS);
    pConv->AddOption("t",OBConversion::OUTOPTIONS);
    return theInChIFormat.WriteMolecule(pOb,pConv);
  };
  
  virtual unsigned int Flags() { return NOTREADABLE;};
};

InChICompareFormat theInChICompareFormat;

//*********************************************************
class TestFormat : public OBMoleculeFormat
{
public:
  //Register this format type ID
  TestFormat() 
  {
    OBConversion::RegisterFormat("test",this);
    OBConversion::RegisterOptionParam("O", this, 1);
    OBConversion::RegisterOptionParam("m", this);
  }

  virtual const char* Description()
  { return
"Test format\n \
Does a round trip conversion and compares before and after using InChI.\n \
Uses the input format unless -xO set\n \
Option e.g. -xOsmi\n \
 O<ext> Test the format that has the specified ID\n \
 m      Output message for each successful molecule\n\n \
Also use opts appropriate for the target format.\n \
";
  };

  virtual unsigned int Flags() { return NOTREADABLE;}
  virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

};
/////////////////////////////////////////////////////////////////
TestFormat theTestFormat;
/////////////////////////////////////////////////////////////////
bool TestFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
  OBMol* pmol = dynamic_cast<OBMol*>(pOb);

  ostream &ofs = *pConv->GetOutStream();
  OBMol &mol = *pmol;
  istream* origInStream = pConv->GetInStream(); 
  ostream* origOutStream = pConv->GetOutStream(); 
  stringstream ssinchi1, ssinchi2;
  
  static int nMols;
  static int nFailures;
  if(pConv->GetOutputIndex()==1)
    nMols=nFailures=0;
  nMols++;

  //mol has already been input using the target format

  //Send its InChI to a stringstream
  OBFormat* pInchi = OBConversion::FindFormat("inchi");
  if(!pInchi)
  {	
    obErrorLog.ThrowError(__FUNCTION__, "InChIFormat needs to be installed to use TestFormat", obWarning);
    return false;
  }
  pConv->AddOption("w",OBConversion::OUTOPTIONS);//no trivial warnings
  pConv->SetOutFormat(pInchi);
  if(!pConv->Write(pmol,&ssinchi1))
    return false;

  //Use a new OBConversion to write and then read using the target format
  //(Using the same OBConversion gave problems with CMLFormat)
  OBConversion NewConv(*pConv);
  NewConv.SetAuxConv(NULL); //until proper copy constructor for OBConversion written
  
  const char* pTargetExt = pConv->IsOption("O");
#ifdef HAVE_SSTREAM
  stringstream errorMsg;
#else
  strstream errorMsg;
#endif
  if(pTargetExt)
  {
    OBFormat* pTargetFormat = OBConversion::FindFormat(pTargetExt);
    if(!pTargetFormat)
    {
      errorMsg << pTargetExt <<  " format is not available" << endl;
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      return false;
    }
    if(!NewConv.SetInFormat(pTargetFormat))
    {
      errorMsg << pTargetExt << " format being tested needs to be readable" << endl;
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      return false;
    }
  }

  if(!NewConv.SetOutFormat(NewConv.GetInFormat()))
  {
    errorMsg << "The input format being tested needs also to be writeable" << endl;
    obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
    return false;
  }
  //Output with the target format to a stringstream ss
  stringstream ss;
  NewConv.SetOneObjectOnly();
  NewConv.SetInStream(&ss); //Quirk in CMLFormat: needs to haveinput stream set before writing
  if(!NewConv.Write(pmol,&ss))
    return false;

  //Read it back in again using the target format
  OBMol Remol;
  if(!NewConv.Read(&Remol))
    return false;

  //Take the InChI of the reconverted molecule
  pConv->SetOutFormat(pInchi);
  if(!pConv->Write(&Remol,&ssinchi2))
    return false;
  
  pConv->SetInStream(origInStream);
  pConv->SetOutStream(origOutStream);
  pConv->SetOutFormat(this);

  char ch = InChIFormat::CompareInchi(ssinchi1.str().c_str(), ssinchi2.str().c_str());

  if(!ch && pConv->IsOption("m"))
  {
    stringstream molID;
    if(strlen(mol.GetTitle())==0)
      molID << "Mol #" << nMols;
    else
      molID << mol.GetTitle();
    ofs << molID.str() << " in " << pConv->GetInFilename();
    if(ch)
    {
      ofs << " and its conversion" << InChIFormat::InChIErrorMessage(ch) << endl;
      ++nFailures;
    }
    else
      ofs << " ok" << endl;
  }

  char s = nFailures==1 ? ' ' : 's';
  if(pConv->IsLast())
    ofs << nFailures << " failure" << s << endl;
  return true;
}

}//namespace OpenBabel

