/**********************************************************************
Copyright (C) 2011 by Chris Morley

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <happyhttp.h>
#include <limits>

using namespace std;
namespace OpenBabel
{
#undef max

class ResolverFormat : public OBMoleculeFormat
{
public:
  //Register this format type ID in the constructor
  ResolverFormat()
  {
    OBConversion::RegisterFormat("web",this);
  }

  virtual const char* Description() //required
  {
    return
    "NIH Web Resolver(Names,InChIKeys)\n"
    "As an input format,\n"
    "it converts IUPAC names via OPSIN, trivial chemical names\n"
    "via its own and ChemSpider's database, and full standard InChIKeys.\n"
    "An input file can have multiple lines, each with a name\n"
    "(which can contain spaces) or full InChiKey (with or without InChIKey=)\n"
    "Lines starting with # are comments and are ignored.\n"
    "The title of the molecule optionally follows after a tab,\n"
    "otherwise uses input string if less than 72 characters.\n"
    "For a single value use the GUI or:\n"
    "  echo \"value text\" | obabel -ires -oxxx\n"
    "(No quotes under Windows)\n"
    " \n "
    "As an output format it produces the IUPAC name.\n"
    "The resolver seems to fail on rather a high proportion\n"
    "of molecules, which is slightly mitigated by removing\n"
    "any uncharged non-bonded fragments from the molecule before\n"
    "submission."
    ;
  };

  virtual const char* SpecificationURL(){return
     "http://cactus.nci.nih.gov/blog/?p=841\n"
     "http://opsin.ch.cam.ac.uk/";
  };

  ////////////////////////////////////////////////////
  virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
  virtual int SkipObjects(int n, OBConversion* pConv);
  virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
private:
  void RemoveSmallUnchargedFragments(OBMol* pmol);
};
  ////////////////////////////////////////////////////

//Make an instance of the format class
ResolverFormat theResolverFormat;

/////////////////////////////////////////////////////////////////

bool ResolverFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{
  OBMol* pmol = pOb->CastAndClear<OBMol>();
  if(pmol==NULL)
    return false;
  istream& ifs = *pConv->GetInStream();
  string ln, title;
  while(getline(ifs, ln))
  {
    if(ln[0]=='#') //comment
      continue;

    //Read title (after tab)
    string::size_type pos;
    pos = ln.find('\t');
    if(pos!=string::npos)
    {
      title = ln.substr(pos).c_str();
      ln.erase(pos);
    }
    else //no title specified
    {
      //use the value as the molecule title, if shorter than 72 chars
      if(ln.size()<72)
        title = ln;
    }

    //Replace spaces
    while((pos=ln.find_first_of(' '))!=string::npos)
      ln.replace(pos, 1, "%20");

    //Use the input text as part of a URL to call web service. Put result in a stream.
    stringstream ss;
    try {
      happyhttp::HttpToStream(&ss, "cactus.nci.nih.gov", "/chemical/structure/" + ln
        + "/sdf?resolver=name_by_opsin,name_by_cir,stdinchikey,name_by_chemspider"
        + "&requester=OpenBabel232");
    }
    catch(happyhttp::Wobbly w) {
      string msg("While resolving "+title+' ');
      obErrorLog.ThrowError(__FUNCTION__, msg + w.what(), obError);
    }
    ss.seekp(0);

    //Convert from return sdf to OBMol
    OBConversion chemConv;
    if(!chemConv.SetInFormat("sdf") || !chemConv.Read(pOb, &ss))
      return false;
    if(!title.empty()) //override nih provided title(the formula)
      pOb->SetTitle(title.c_str());
    return true;
  }
  return false;
}

////////////////////////////////////////////////////////////////////////
int ResolverFormat::SkipObjects(int n, OBConversion* pConv)
{
  if(n==0) return 1; //already points after current line
  istream& ifs = *pConv->GetInStream();
  if (ifs.eof())
    return -1;

  int i=0;
  while(i<n && ifs.good())
    {
      if(ifs.peek()!='#')
        i++;
      ifs.ignore(numeric_limits<streamsize>::max(),'\n');
    }
  return ifs ? 1 : -1;
}

////////////////////////////////////////////////////////////////////////
bool ResolverFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
  OBMol* pmol = dynamic_cast<OBMol*>(pOb);
  if(!pmol)
    return false;
  ostream &ofs = *pConv->GetOutStream();

  RemoveSmallUnchargedFragments(pmol);
  
  OBConversion sconv;
  sconv.SetOutFormat("smi");
  sconv.AddOption("n"); // no molecule name
  string smiles = sconv.WriteString(pmol, true);

  //Escape some characters because used in a URL
  string::size_type pos;
  while( (pos = smiles.find_first_of("#[]"))!=string::npos)
  {
    if(smiles[pos]=='#')
      smiles.replace(pos,1,"%23");
    if(smiles[pos]=='[')
      smiles.replace(pos,1,"%5B");
    if(smiles[pos]==']')
      smiles.replace(pos,1,"%5D");
  }

  //Use the SMILES as part of a URL to call web service. Put result in a stream.
  stringstream ss;
  try {
    happyhttp::HttpToStream(&ss, "cactus.nci.nih.gov", "/chemical/structure/" + smiles
      + "/iupac_name");
  }
  catch(happyhttp::Wobbly w) {
    ss << pOb->GetTitle()<< ' ' << w.what();
  }
  ss.seekp(0);
  ofs << ss.rdbuf() << endl;
  return true;
}

class charged {
private:
  OBMol* _pmol;
public:
  charged(OBMol* pmol) : _pmol(pmol){}
  bool operator()(int idx) {
    return (_pmol->GetAtom(idx))->GetFormalCharge();
  }
};

void ResolverFormat::RemoveSmallUnchargedFragments(OBMol* pmol)
{
  //e.g. from  O.[Na+].CC(=O)[O-]  remove water but not Na+
  vector<vector<int> > cfl;
  vector<vector<int> >::iterator it,biggestfrag;
  pmol->ContigFragList(cfl);
  if (!cfl.empty() && cfl.size() != 1)
  {
    vector<int> unwanted;
    for(it=cfl.begin()+1;it!=cfl.end();++it)
    {
      if(count_if(it->begin(), it->end(), charged(pmol))>0)
        continue;
      //note the atoms in uncharged fragments
      unwanted.insert(unwanted.end(), it->begin(), it->end());
    }
    sort(unwanted.begin(), unwanted.end()); //smallest first
    for(;!unwanted.empty();unwanted.pop_back())//delete from top down to avoid index invalidation
      pmol->DeleteAtom(pmol->GetAtom(unwanted.back()));
  }
}

} //namespace OpenBabel

