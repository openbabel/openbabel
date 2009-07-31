/**********************************************************************
Copyright (C) 2005 by Chris Morley
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#include <openbabel/babelconfig.h>

#include <vector>
#include <string>
#include <iomanip>

#include <openbabel/obmolecformat.h>
#include <openbabel/fingerprint.h>

using namespace std;
namespace OpenBabel
{

  /// \brief Constructs and displays fingerprints. For details see OBFingerprint class
  class FingerprintFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    FingerprintFormat() {OBConversion::RegisterFormat("fpt",this);}

    virtual const char* Description() //required
    { return
      "Fingerprint format\n"
      "Constructs and displays fingerprints and (for multiple input objects)\n"
      "the Tanimoto coefficient and whether a superstructure of the first object\n"
      "Output options e.g. -xfFP3 -xN128\n"
      " f<id> fingerprint type\n"
      " N# fold to specified number of bits, 32, 64, 128, etc.\n"
      " h  hex output when multiple molecules\n"
      " s  describe each set bit\n"
      " u  describe each unset bit\n"
;
    };

    virtual unsigned int Flags(){return NOTREADABLE;};
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

  private:
    vector<unsigned int> firstfp;
    string firstname;
    bool IsPossibleSubstructure(vector<unsigned int>Mol, vector<unsigned int>Frag);
  };

  ////////////////////////////////////////////////////
  //Make an instance of the format class
  FingerprintFormat theFingerprintFormat;

  //*******************************************************************
  bool FingerprintFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    ostream &ofs = *pConv->GetOutStream();

    bool hexoutput=false;
    if(pConv->IsOption("h") || (pConv->GetOutputIndex()==1 && pConv->IsLast()))
      hexoutput=true;

    string fpid;
    int nbits=0;
    const char* p=pConv->IsOption("f");
    if(p)
      {
        fpid=p;
        fpid = fpid.substr(0,fpid.find('"'));
      }

    OBFingerprint* pFP = OBFingerprint::FindFingerprint(fpid.c_str());
    if(!pFP)
      {
        stringstream errorMsg;
        errorMsg << "Fingerprint type '" << fpid << "' not available" << endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
        return false;
      }

    p=pConv->IsOption("N");
    if(p)
      nbits = atoi(p);
    if(nbits<0)
      obErrorLog.ThrowError(__FUNCTION__,
      "The number of bits to fold to, in the-xN option, should be >=0", obWarning);

    vector<unsigned int> fptvec;
    if(!pFP->GetFingerprint(pOb, fptvec, nbits))
      return false;
	
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol)
      ofs << ">" << pmol->GetTitle();
    
    // checkmol-type output
    if(pConv->IsOption("s") || pConv->IsOption("u"))
    {
      if(nbits!=0)
      {
        obErrorLog.ThrowError(__FUNCTION__,
        "The fingerprint must be unfolded when describing bits.", obError);
        return false;
      }
      string descr = pFP->DescribeBits(fptvec, pConv->IsOption("s")!=NULL);
      if(descr=="")
        obErrorLog.ThrowError(__FUNCTION__,
        "Bit descriptions are not available for this fingerprint type", obError, onceOnly);

      ofs << descr;
      return true;
    }

    if(hexoutput && pConv->GetOutputIndex()<=1)
      {
        unsigned int i, bitsset=0;
        for (i=0;i<fptvec.size();++i)
          {
            int wd = fptvec[i];
            for(;wd;wd=wd<<1)//count bits set by shifting into sign bit until word==0
              if(wd<0) ++bitsset;
          }
        ofs  << "   " << bitsset << " bits set "; 
      }

    if(pConv->GetOutputIndex()<=1)
      {
        //store the fingerprint and name of first molecule
        firstfp=fptvec;
        if(pmol)
          firstname=pmol->GetTitle();
        if(firstname.empty())
          firstname = "first mol";


    }
    else
      {
        ofs << "   Tanimoto from " << firstname << " = " << OBFingerprint::Tanimoto(firstfp, fptvec);
        if(IsPossibleSubstructure(fptvec,firstfp))
          ofs << "\nPossible superstructure of " << firstname;
      }
    ofs << endl;
	
    int i;

    if(hexoutput)
      {
        for(i=fptvec.size()-1;i>=0;i--)
          {
            ofs << hex << setfill('0') << setw(8) << fptvec[i] << " " ;
            if((fptvec.size()-i)%6==0)
              ofs <<endl;
          }
        ofs << dec << endl;
      }
    return true;
  }

  bool FingerprintFormat::IsPossibleSubstructure(vector<unsigned int>Mol, vector<unsigned int>Frag)
  {
    //Returns false if Frag is definitely NOT a substructure of Mol
    unsigned int i;
    for (i=0;i<Mol.size();++i)
      if((Mol[i] & Frag[i]) ^ Frag[i]) return false;
    return true;
  }

} //namespace OpenBabel
