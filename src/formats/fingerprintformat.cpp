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
#include <cstdlib>

#include <openbabel/obmolecformat.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/elements.h>

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
      "Generate or display molecular fingerprints.\n"
"This format constructs and displays fingerprints and (for multiple input\n"
"objects) the Tanimoto coefficient and whether a superstructure of the first\n"
"object.\n\n"

"A list of available fingerprint types can be obtained by::\n\n"

"  babel -L fingerprints\n\n"

"The current default type FP2 is is of the Daylight type, indexing a molecule\n"
"based on the occurrence of linear fragment up to 7 atoms in length. To use a\n"
"fingerprint type other than the default, use the ``-xf`` option, for example::\n\n"

"  babel infile.xxx -ofpt -xfFP3\n\n"

"For a single molecule the fingerprint is output in hexadecimal form\n"
"(intended mainly for debugging).\n\n"

"With multiple molecules the hexadecimal form is output only if the ``-xh``\n"
"option is specified. But in addition the Tanimoto coefficient between the\n"
"first molecule and each of the subsequent ones is displayed. If the first\n"
"molecule is a substructure of the target molecule a note saying this is\n"
"also displayed.\n\n"

"The Tanimoto coefficient is defined as::\n\n"

" Number of bits set in (patternFP & targetFP) / Number of bits in (patternFP | targetFP)\n\n"

"where the boolean operations between the fingerprints are bitwise.\n\n"

"The Tanimoto coefficient has no absolute meaning and depends on the design of the fingerprint.\n\n"


"Use the ``-xs`` option to describe the bits that are set in the fingerprint.\n"
"The output depends on the fingerprint type. For Fingerprint FP4, each bit\n"
"corresponds to a particular chemical feature, which are specified as SMARTS\n"
"patterns in :file:`SMARTS_InteLigand.txt`, and the output is a tab-separated\n"
"list of the features of a molecule. For instance, a well-known molecule\n"
"gives::\n\n"

" Primary_carbon: Carboxylic_acid: Carboxylic_ester: Carboxylic_acid_derivative:\n"
" Vinylogous_carbonyl_or_carboxyl_derivative: Vinylogous_ester: Aromatic:\n"
" Conjugated_double_bond: C_ONS_bond: 1,3-Tautomerizable: Rotatable_bond: CH-acidic:\n\n"

"For the path-based fingerprint FP2, the output from the ``-xs`` option is\n"
"instead a list of the chemical fragments used to set bits, e.g.::\n\n"

" $ obabel -:\"CCC(=O)Cl\" -ofpt -xs -xf FP2\n"
" >\n"
" 0 6 1 6 <670>\n"
" 0 6 1 6 1 6 <260>\n"
" 0 8 2 6 <623>\n"
" ...etc\n\n"

"where the first digit is 0 for linear fragments but is a bond order\n"
"for cyclic fragments. The remaining digits indicate the atomic number\n"
"and bond order alternatively. Note that a bond order of 5 is used for\n"
"aromatic bonds. For example, bit 623 above is the linear fragment O=C\n"
"(8 for oxygen, 2 for double bond and 6 for carbon).\n\n"

      "Write Options e.g. -xfFP3 -xN128\n"
      " f<id> fingerprint type\n"
      " N# fold to specified number of bits, 32, 64, 128, etc.\n"
      " h  hex output when multiple molecules\n"
      " o  hex output only\n"
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
    bool WriteHex(ostream &ofs, vector<unsigned int> fptvec);
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

    if(pConv->IsOption("o"))
      return WriteHex(ofs, fptvec);

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

      ofs << '\n' << descr;
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

    if(hexoutput)
      {
        WriteHex(ofs, fptvec);
        ofs << endl;
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

  bool FingerprintFormat::WriteHex(ostream &ofs, vector<unsigned int> fptvec)
  {
    for(int i=fptvec.size()-1;i>=0;i--)
    {
      ofs << hex << setfill('0') << setw(8) << fptvec[i] << " " ;
      if((fptvec.size()-i)%6==0)
        ofs <<endl;
    }
    ofs << dec << flush;
    return true;
  }
} //namespace OpenBabel
