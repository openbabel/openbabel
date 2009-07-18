/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey Hutchison
Portions Copyright (C) 2004-2006 by Chris Morley

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
#include <openbabel/babelconfig.h>

#ifdef _WIN32
#pragma warning (disable : 4786)
#endif

#include <ctime>
#include <vector>
#include <iomanip>
#include <map>
#include <algorithm>
#include <openbabel/obmolecformat.h>
#include <openbabel/stereo/stereo.h>
#include <openbabel/stereo/cistrans.h>
#include <openbabel/stereo/tetrahedral.h>
#include <openbabel/alias.h>
#include <openbabel/tokenst.h>

using namespace std;
namespace OpenBabel
{

  //MDLFormat is a base class which is never instantiated.
  //MOLFormat and SDFormat are derived from it and have their own constructors.
  //SDFormat has its own WriteMolecule which sets the "sd" option.

  class MDLFormat : public OBMoleculeFormat
  {
    public:
      virtual const char* Description()
      { 
        return "MDL MOL format\n"
               "Reads and writes V2000 and V3000 versions\n"
               "Write Options, e.g. -x3\n"
            /* " 2  output V2000 (default) or\n" */
               " 3  output V3000 not V2000 (used for >999 atoms/bonds) \n"
               " m  write no properties\n";
      }

      virtual const char* SpecificationURL()
      {
        return "http://www.mdl.com/downloads/public/ctfile/ctfile.jsp";
      }

      virtual const char* GetMIMEType() 
      { 
        return "chemical/x-mdl-molfile"; 
      }

      virtual unsigned int Flags() { return DEFAULTFORMAT; }
      virtual const char* TargetClassDescription() { return OBMol::ClassDescription(); }

      virtual int SkipObjects(int n, OBConversion* pConv)
      {
        if (n == 0) 
          n++;
        istream& ifs = *pConv->GetInStream();
        do {
         ignore(ifs, "$$$$\n");   
        } while(ifs && --n);
        return ifs.good() ? 1 : -1;
      }

      ////////////////////////////////////////////////////
      /// The "API" interface functions
      virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
      virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

      ////////////////////////////////////////////////////
      //V3000 routines
    private:
      bool ReadV3000Block(istream& ifs, OBMol& mol, OBConversion* pConv,bool DoMany);
      bool ReadV3000Line(istream& ifs, vector<string>& vs);
      bool ReadAtomBlock(istream& ifs,OBMol& mol, OBConversion* pConv);
      bool ReadBondBlock(istream& ifs,OBMol& mol, OBConversion* pConv);
      bool ReadCollectionBlock(istream& ifs,OBMol& mol, OBConversion* pConv);
      bool WriteV3000(ostream& ofs,OBMol& mol, OBConversion* pConv);
    private:
      enum Parity {
        NotStereo, Clockwise, AntiClockwise, Unknown
      };
      bool  HasProperties;
      string GetTimeDate();
      void GetUpDown(OBMol& mol, map<OBBond*, OBStereo::BondDirection> &updown, set<OBBond*> &stereodbl);
      void GetParity(OBMol& mol, map<OBAtom*, Parity> &parity,
           map<OBBond*, OBStereo::BondDirection> &updown);
      map<int,int> indexmap; //relates index in file to index in OBMol
      vector<string> vs;
  };

  //**************************************
  class MOLFormat : public MDLFormat
  {
    public:
      //Register this format type ID
      MOLFormat() 
      {
        OBConversion::RegisterFormat("mol",this, "chemical/x-mdl-molfile");
        OBConversion::RegisterFormat("mdl",this, "chemical/x-mdl-molfile");
        OBConversion::RegisterOptionParam("2", this);
        OBConversion::RegisterOptionParam("3", this);
      }
  };
  
  //Make an instance of the format class
  MOLFormat theMOLFormat;

  //*************************************
  class SDFormat : public MDLFormat
  {
    public:
      SDFormat()
      {
        OBConversion::RegisterFormat("sd",this, "chemical/x-mdl-sdfile");
        OBConversion::RegisterFormat("sdf",this, "chemical/x-mdl-sdfile");
      }

      virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv)
      {
        //The sd option ensures that a $$$$ is written at the end of the file
        pConv->AddOption("sd", OBConversion::OUTOPTIONS);
        return MDLFormat::WriteMolecule(pOb, pConv);
      }
  };

  //Make an instance of the format class
  SDFormat theSDFormat;

  /////////////////////////////////////////////////////////////////
  bool MDLFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = pOb->CastAndClear<OBMol>();

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    bool setDimension = false; // did we extract the 'dimensional code' from line 2?
    stringstream errorMsg;
    string clearError; // empty string to clear the warning buffer

    int i, natoms, nbonds;
    char buffer[BUFF_SIZE];
    string comment;
    string r1, r2;
    map<OBBond*, OBStereo::BondDirection> updown;

    // Attempting to read past the end of the file -- don't bother
    if ( !ifs.good() || ifs.peek() == EOF ) 
      return false;

    std::string line;
    // 
    // The Header Block
    //

    // line1: molecule name
    if (!std::getline(ifs, line)) {
      errorMsg << "WARNING: Problems reading a MDL file\n";
      errorMsg << "Cannot read title line\n";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
      return(false);
    }
    mol.SetTitle(line);
  
    // line 2: IIPPPPPPPPMMDDYYHHmmddSSssssssssssEEEEEEEEEEEERRRRRR
    //
    //          0...1    I = user's initials
    //          2...9    P = program name
    //         10..19    M/D/Y,H:m = date/time
    //         20..21    d = dimentional code
    //         22..23    S = scaling factor
    //         24..33    s = scaling facter (double format 10.5)
    //         34..45    E = energy
    //         46..51    R = internal registry number
    if (!std::getline(ifs, line)) {
      errorMsg << "WARNING: Problems reading a MDL file\n";
      errorMsg << "Cannot read creator/dimension line line\n";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
      return false;
    }
    
    if (line.size() > 21) {
      string dim = line.substr(20, 2);
      if (dim == "3D") {
        mol.SetDimension(3);
        setDimension = true;
      } else 
      if (dim == "2D") {
        mol.SetDimension(2);
        setDimension = true;
      } else
      if (dim == "0D") {
        mol.SetDimension(0);
        setDimension = true;
      }
    }

    // line 3: comment line
    if (!std::getline(ifs, line)) {
      errorMsg << "WARNING: Problems reading a MDL file\n";
      errorMsg << "Cannot read comment line\n";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
      return false;
    }

    if (!line.empty())
      comment = line;

    //
    // Connection Table (Ctab)
    //

    // line 1: counts line
    if (!std::getline(ifs, line)) {
      errorMsg << "WARNING: Problems reading a MDL file\n";
      errorMsg << "Cannot read atom and bond count\n";
      errorMsg << "File ended prematurely\n";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
      return false;
    }
    if (line.size() < 6) { // error from Joe Bedell, Sigma-Aldrich
      errorMsg << "WARNING: Problems reading a MDL file\n";
      errorMsg << "Cannot read atom and bond count\n";
      errorMsg << "Expected standard 6 character atom and bond count\n";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
      return(false);
    }

    natoms = atoi((line.substr(0, 3)).c_str());
    nbonds = atoi((line.substr(3, 3)).c_str());

    mol.BeginModify();
    if(line.find("V3000") != string::npos) {
      // V3000
      indexmap.clear();
      if(!ReadV3000Block(ifs, mol, pConv, false)) 
        return false;
      //ifs.getline(buffer,BUFF_SIZE); //M END line
    } else {
      // V2000
      mol.ReserveAtoms(natoms);
      double x,y,z;
      string symbol;
      OBAtom atom;

      //
      // Atom Block
      //
      int massdiff, charge;
      vector<int> massDiffs, charges;
      for (i = 0; i < natoms; ++i) {
        if (!std::getline(ifs, line)) {
          errorMsg << "WARNING: Problems reading a MDL file\n";
	        errorMsg << "Not enough atoms to match atom count (" << natoms << ") in counts line\n";
	        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
          return false;
	      }

        // xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee
        //
        // 0...30   x y z = atom coordinates
        // 31..33   aaa = atom symbol
        // 34..35   dd = mass difference: -3, -2, -1, 0, 1, 2, 3, 4 ('M  ISO' lines take precedence)
        // 36..38   ccc = charge  ('M  CGH' and 'M  RAD' lines take precedence)
        // 39..41   sss = atom stereo parity (ignored when read)
        //          ... = query/reaction related
        massdiff = charge = 0;
        if (line.size() < 34) {
	        errorMsg << "WARNING: Problems reading a MDL file\n";
	        errorMsg << "Missing data following atom specification in atom block\n";
	        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
          return false;
        }
        // coordinates
        x = atof(line.substr(0, 10).c_str());
        y = atof(line.substr(10, 10).c_str());
        z = atof(line.substr(20, 10).c_str());
        atom.SetVector(x, y, z);
        // symbol & isotope
        symbol = line.substr(31, 3);
        Trim(symbol);
        atom.SetAtomicNum(etab.GetAtomicNum(symbol.c_str()));
        // mass difference
        if (line.size() >= 35)
          massdiff = atoi(line.substr(34, 2).c_str());
        massDiffs.push_back(massdiff);
        // charge      
        if (line.size() >= 38)
          charge = atoi(line.substr(36, 3).c_str());
        charges.push_back(charge);

        if (!mol.AddAtom(atom))
          return false;
        atom.Clear();
      }

      //
      // Bond Block
      //
      int stereo = 0;
      unsigned int begin, end, order, flag;
      for (i = 0;i < nbonds; ++i) {
        flag = 0;
        if (!std::getline(ifs, line)) {
          errorMsg << "WARNING: Problems reading a MDL file\n";
	        errorMsg << "Not enough bonds to match bond count (" << nbonds << ") in counts line\n";
	        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
                return false;
	      }
	      begin = end = order = 0;
	      // 111222tttsssxxxrrrccc
        //
        // 111 = first atom number
        // 222 = second atom number
        // ttt = bond type (1-3, 4 = aromatic, 4-8 = query)
        // sss = bond stereo (
        // ... = query/topology
        if (line.size() >= 9) {
	        begin = atoi((line.substr(0, 3)).c_str());
	        end = atoi((line.substr(3, 3)).c_str());
	        order = atoi((line.substr(6, 3)).c_str());
	      }
        if (begin == 0 || end == 0 || order == 0 || begin > mol.NumAtoms() || end > mol.NumAtoms()) {
	        errorMsg << "WARNING: Problems reading a MDL file\n";
	        errorMsg << "Invalid bond specification, atom numbers or bond order are wrong.\n";
	        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
          return false;
	      }

        order = (order == 4) ? 5 : order;
        if (line.size() >= 12) {  //handle wedge/hash data
          stereo = atoi((line.substr(9, 3)).c_str());
          if (stereo) {
            switch (stereo) {
              case 1:
                // single bond: wedge
                flag |= OBBond::Wedge;
                break;
              case 3:
                // double bond: either cis or trans
                flag |= OBBond::CisOrTrans;
              case 4:
                // single bond: either wedge or hash (unspecified)
                flag |= OBBond::WedgeOrHash;
                break;
              case 6:
                // single bond: hash
                flag |= OBBond::Hash;
                break;
              default:
                // single bonds: not stereo
                // double bonds: use x,y,z coordinates
                break;
            }
          }
        }

        if (!mol.AddBond(begin,end,order,flag)) {
	        errorMsg << "WARNING: Problems reading a MDL file\n";
	        errorMsg << "Invalid bond specification\n";
	        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
	        return false;
	      }
        if (stereo) {
          OBStereo::BondDirection bd;
          switch (stereo) {
            case 1: 
              bd = OBStereo::UpBond;
              break;
            case 6:
              bd = OBStereo::DownBond;
              break;
            case 4:
              bd = OBStereo::UnknownDir;
              break;
            default:
              bd = OBStereo::NotStereo;
              break;
          }
          if (bd != OBStereo::NotStereo)
            updown[mol.GetBond(begin, end)] = bd;
        }
      }

      // 
      // Properties Block (currently only M RAD and M CHG M ISO)
      //
      bool foundISO = false, foundCHG = false;
      while (std::getline(ifs, line)) {
        if (line.substr(0, 4) == "$$$$")
          return true;
        if (line.substr(0, 6) == "M  END")
          break;
        if (line.substr(0, 6) == "S  SKP") {
          int i = atoi(line.substr(6, line.size() - 6).c_str());
          for(; i > 0; --i)
            if (ifs.good()) // check for EOL, suggested by Dalke
              std::getline(ifs, line);
        }
        if (line.at(0) == 'A') { //alias
          int atomnum = atoi(line.substr(2, line.size() - 2).c_str());
          for(; i > 0; --i)
          std::getline(ifs, line);
          if(line.at(0) != '?' && line.at(0) != '*') {
            AliasData* ad = new AliasData();
            ad->SetAlias(line);
            ad->SetOrigin(fileformatInput);
            OBAtom* at = mol.GetAtom(atomnum);
            if (at) {
              at->SetData(ad);
              at->SetAtomicNum(0);
              //The alias has now been added as a dummy atom with a AliasData object.
              ad->Expand(mol, atomnum); //Make chemically meaningful, if possible.
            }
          }
          continue;
        }

        if ((line.substr(0, 6) != "M  CHG") && (line.substr(0, 6) != "M  RAD") &&
            (line.substr(0, 6) != "M  ISO"))
          continue;
        int n = -1;
        if (line.size() >= 9)
          n = atoi((line.substr(6, 3)).c_str()); //entries on this line
        if (n <= 0 || n > 8 || 6+n*8 > line.size()) { //catch ill-formed line
          obErrorLog.ThrowError(__FUNCTION__, "Error in line:\n" + line, obError);
          return false;
        }
        int pos = 10;
        for (; n > 0; n--, pos += 8) {
          int atomnumber = atoi((line.substr(pos,3)).c_str());
          OBAtom *at;
          if (atomnumber==0 || (at=mol.GetAtom(atomnumber))==NULL) {
            obErrorLog.ThrowError(__FUNCTION__, "Error in line:\n" + line, obError);
            return false;
          }

          at = mol.GetAtom(atomnumber); //atom numbers start at 1
          int value = atoi((line.substr(pos+4,3)).c_str());
          if (line.substr(3, 3) == "RAD") {
            at->SetSpinMultiplicity(value);
            foundCHG = true;
          } else if (line.substr(3, 3) == "CHG") {
            at->SetFormalCharge(value);
            foundCHG = true;
          } else if (line.substr(3, 3) == "ISO") {
            if (value)
              at->SetIsotope(value);
            foundISO = true;
          }
        }
        // Lines setting several other properties are not implemented
      }

      // if no 'M  ISO' properties are found, use the mass differences from the atom block
      if (!foundISO)
        FOR_ATOMS_OF_MOL (a, mol) {
          int massDifference = massDiffs.at(a->GetIndex());
          if (massDifference)
            a->SetIsotope(etab.GetMass(a->GetAtomicNum()) + massDifference);
        }
      // if no 'M  CHG' or 'M  RAD' properties are found, use the charges from the atom block
      if (!foundCHG)
        FOR_ATOMS_OF_MOL (a, mol) {
          charge = charges.at(a->GetIndex());
          switch (charge) {
            case 0: break;
            case 3: a->SetFormalCharge(1); break;
            case 2: a->SetFormalCharge(2); break;
            case 1: a->SetFormalCharge(3); break;
            case 5: a->SetFormalCharge(-1); break;
            case 6: a->SetFormalCharge(-2); break;
            case 7: a->SetFormalCharge(-3); break;
          }
        }
    }

    mol.AssignSpinMultiplicity();
    mol.EndModify();
        
    if (comment.length()) {
      OBCommentData *cd = new OBCommentData;
      cd->SetData(comment);
      cd->SetOrigin(fileformatInput);
      mol.SetData(cd);
    }
        
    //Get property lines
    while (std::getline(ifs, line)) {
      if (line.find("<") != string::npos) {
        size_t lt = line.find("<")+1;
        size_t rt = line.find_last_of(">");
        string attr = line.substr(lt, rt - lt);

        // sometimes we can hit more data than BUFF_SIZE, so we'll use a std::string
        string buff;
        while (std::getline(ifs, line)) {
          Trim(line);
          if (line.size()) {
            buff.append(line);
            buff += "\n";
          } else
            break;
        }
        Trim(buff);

        OBPairData *dp = new OBPairData;
        dp->SetAttribute(attr);
        dp->SetValue(buff);
        dp->SetOrigin(fileformatInput);
        mol.SetData(dp);
            
        if(!strcasecmp(attr.c_str(),"NAME") && *mol.GetTitle()=='\0')
          mol.SetTitle(buff);
      }
      if (line.substr(0, 4) ==  "$$$$") 
        break;
      if (line.substr(0, 4) == "$MOL") 
        break;
    }

    if (mol.Has3D()) {
      if (!setDimension)
        mol.SetDimension(3);
      // use 3D coordinates to determine stereochemistry
      StereoFrom3D(&mol);
    } else 
    if (mol.Has2D()) {
      if (!setDimension)
        mol.SetDimension(2);
      // use 2D coordinates + hash/wedge to determine stereochemistry
      StereoFrom2D(&mol);
    } else {
    if (!setDimension)
      mol.SetDimension(0);
    // use up/down to determine stereochemistry
    StereoFrom0D(&mol, &updown);
    }

    return true;
  }

  /////////////////////////////////////////////////////////////////
  bool MDLFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    if (pConv->GetOutputIndex()==1)
      HasProperties = false;

    // 
    // Header Block
    //

    string dimension("2D");
    if(mol.GetDimension()==3)
      dimension = "3D";

    // line 1: molecule name
    ofs << mol.GetTitle() <<  endl;

    // line 2: Program name, date/time, dimensions code
    ofs << " OpenBabel" << GetTimeDate() <<  dimension << endl; //line2

    // line 3: comment
    if (mol.HasData(OBGenericDataType::CommentData)) {
      OBCommentData *cd = (OBCommentData*)mol.GetData(OBGenericDataType::CommentData);
      ofs << cd->GetData(); 
    } 
    ofs << endl;
       
    //
    // Atom Block
    //

    if(pConv->IsOption("3") || mol.NumAtoms() > 999 || mol.NumBonds() > 999) {
      if (!WriteV3000(ofs, mol, pConv)) 
        return false;
    } else {
      //The rest of the function is the same as the original
      char buff[BUFF_SIZE];  

      if (mol.NumAtoms() > 999 || mol.NumBonds() > 999) { // Three digits!
        stringstream errorMsg;
        errorMsg << "MDL Molfile conversion failed: Molecule is too large to convert." << endl;
        errorMsg << "  File format (v2000) is limited to 999 atoms or bonds." << endl;
        errorMsg << "  Molecule size: " << mol.NumAtoms() << " atoms ";
        errorMsg << "and " << mol.NumBonds() << " bonds." << endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
        return false;
      }

      // Check to see if there are any untyped aromatic bonds (GetBO == 5)
      // These must be kekulized first
      FOR_BONDS_OF_MOL(b, mol) {
        if (b->GetBO() == 5) {
          mol.Kekulize();
          break;
        }
      }

      // Calculate up/downness of cis/trans bonds
      map<OBBond*, OBStereo::BondDirection> updown;
      set<OBBond*> stereodbl;
      GetUpDown(mol, updown, stereodbl);

      // Calculate parity of atoms (see Appendix A of ctfile.pdf)
      map<OBAtom*, Parity> parity;
      GetParity(mol, parity, updown);
                      
      // The counts line:
      // aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
      //
      // aaa = number of atoms
      // bbb = number of bonds
      // lll = number of atom lists (query)
      // ccc = chiral flag
      // ... = obsolete
      // mmm = no longer supported (default=999)
      //                         aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
      ofs << setw(3) << mol.NumAtoms() << setw(3) << mol.NumBonds();
      ofs << "  0  0  0  0  0  0  0  0999 V2000" << endl;

      OBAtom *atom;
      vector<OBAtom*>::iterator i;
      int charge = 0; 
      for (atom = mol.BeginAtom(i); atom; atom = mol.NextAtom(i)) {
        // convert charge
        switch (atom->GetFormalCharge()) {
          case 1: charge = 3; break;
          case 2: charge = 2; break;
          case 3: charge = 1; break;
          case -1: charge = 5; break;
          case -2: charge = 6; break;
          case -3: charge = 7; break;
          default: charge = 0; break;
        }
        Parity stereo = NotStereo;
        if (parity.find(atom) != parity.end())
          stereo = parity[atom];
        snprintf(buff, BUFF_SIZE, "%10.4f%10.4f%10.4f %-3s%2d%3d%3d%3d%3d",
                 atom->GetX(), atom->GetY(), atom->GetZ(),
                 atom->GetAtomicNum() ? etab.GetSymbol(atom->GetAtomicNum()) : "* ",
                 0, charge, stereo, 0, 0);    
          ofs << buff << endl;
        }

        //so the bonds come out sorted
        OBAtom *nbr;
        OBBond *bond;
        vector<OBBond*>::iterator j;
        for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i)) {
          for (nbr = atom->BeginNbrAtom(j);nbr;nbr = atom->NextNbrAtom(j))
            if (atom->GetIdx() < nbr->GetIdx()) {
              bond = (OBBond*) *j;
                                        
              int stereo = 0; //21Jan05 CM
              if(dimension == "2D") {
                if (bond->IsWedge()) 
                  stereo = 1;
                else if (bond->IsHash()) 
                  stereo = 6;
                else if (bond->IsWedgeOrHash())
                  stereo = 4;

                // For Cis/Trans bonds, set the stereo
                if (updown.find(bond) != updown.end())
                  stereo = updown[bond];
                // For Cis/Trans double bonds, set the stereo
                if (stereodbl.find(bond) != stereodbl.end())
                  stereo = 3;
              }

              ofs << setw(3) << bond->GetBeginAtomIdx(); // begin atom number
              ofs << setw(3) << bond->GetEndAtomIdx(); // end atom number
              ofs << setw(3) << bond->GetBO(); // bond type
              ofs << setw(3) << stereo; // bond stereo
              ofs << "  0  0" << endl;
            }
        }

        vector<OBAtom*> rads, isos, chgs;
        vector<OBAtom*>::iterator itr;
        for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
          {
            if(atom->GetSpinMultiplicity())
              rads.push_back(atom);
            if(atom->GetIsotope())
              isos.push_back(atom);
            if(atom->GetFormalCharge())
              chgs.push_back(atom);

            if(atom->HasData(AliasDataType))
            {
              AliasData* ad = static_cast<AliasData*>(atom->GetData(AliasDataType));
              if(!ad->IsExpanded()) //do nothing with an expanded alias
                ofs << "A  " << atom->GetIdx() << '\n' << ad->GetAlias() << endl;
            }

          }
        if(rads.size())
          {
            ofs << "M  RAD" << setw(3) << rads.size();
            for(itr=rads.begin();itr!=rads.end();++itr)
              ofs << setw(4) << (*itr)->GetIdx() << setw(4) << (*itr)->GetSpinMultiplicity();
            ofs << endl;
          }                     
        if(isos.size())
          {
            ofs << "M  ISO" << setw(3) << isos.size();
            for(itr=isos.begin();itr!=isos.end();++itr)
              ofs << setw(4) << (*itr)->GetIdx() << setw(4) << (*itr)->GetIsotope();
            ofs << endl;
          }                     
        if(chgs.size())
          {
            int counter = 0;
            for( itr=chgs.begin(); itr != chgs.end(); ++itr, counter++ ) {
              if( counter % 8 == 0 ) {
                if( counter > 0 ) ofs << endl;
                ofs << "M  CHG" << setw(3) << min(static_cast<unsigned long int>(chgs.size() - counter), static_cast<unsigned long int>(8));
              }
              ofs << setw(4) << (*itr)->GetIdx() << setw(4) << (*itr)->GetFormalCharge();
            }
            ofs << endl;
          }
      }

    ofs << "M  END" << endl;

    //For SD files only, write properties unless option m
    if(pConv->IsOption("sd") && !pConv->IsOption("m"))
      {
        vector<OBGenericData*>::iterator k;
        vector<OBGenericData*> vdata = mol.GetData();
        for (k = vdata.begin();k != vdata.end();k++)
          {
            if ((*k)->GetDataType() == OBGenericDataType::PairData)
              {
                HasProperties = true;
                //Since partial charges are not output
                //in this format, don't need the annotation
                if((*k)->GetAttribute()!="PartialCharges")
                {
                  ofs << ">  <" << (*k)->GetAttribute() << ">" << endl;
                  ofs << ((OBPairData*)(*k))->GetValue() << endl << endl;
                }
              }
          }
      }
        
    //Unless option no$$$$ is set, $$$$ is always written between molecules and
    //at the end any if properties have been output in any molecule,
    //or if the sd option is set.
    if(!pConv->IsOption("no$$$$"))
      if(!pConv->IsLast()  || HasProperties  || pConv->IsOption("sd"))  
        ofs << "$$$$" << endl;

    return(true);
  }


  //////////////////////////////////////////////////////
  bool MDLFormat::ReadV3000Block(istream& ifs, OBMol& mol, OBConversion* pConv,bool DoMany)
  {
    do
      {
        if(!ReadV3000Line(ifs,vs)) return false;
        if(vs[2]=="LINKNODE"){continue;} //not implemented
        if(vs[2]!="BEGIN") return false;

        if(vs[3]=="CTAB")
          {
            if(!ReadV3000Line(ifs,vs) || vs[2]!="COUNTS") return false;
            int natoms = atoi(vs[3].c_str());
            //int nbonds = atoi(vs[4].c_str());
            //int chiral = atoi(vs[7].c_str()); 
            //number of s groups, number of 3D contraints, chiral flag and regno not yet implemented
            mol.ReserveAtoms(natoms);

            ReadV3000Block(ifs,mol,pConv,true);//go for contained blocks        
            if(!ReadV3000Line(ifs,vs) || (vs[1]!="END" && vs[3]!="CTAB")) return false;
            return true;
          }
        else if(vs[3]=="ATOM")
          ReadAtomBlock(ifs,mol,pConv);
        else if(vs[3]=="BOND")
          ReadBondBlock(ifs,mol,pConv);
        else if(vs[3]=="COLLECTION")
          ReadCollectionBlock(ifs,mol,pConv);
          
        /*
          else if(vs[3]=="3D")
          //not currently implemented
          else if(vs[3]=="SGROUP")
          //not currently implemented
          else if(vs[3]=="RGROUP")
          //not currently implemented
        */
      }while(DoMany && ifs.good());
    //  if(is3D){mol.SetDimension(3);cout<<"SetDim to 3"<<endl;}
    //  else if(is2D){mol.SetDimension(2);cout<<"SetDim to 2"<<endl;}
    return true;
  }

  //////////////////////////////////////////////////////
  bool MDLFormat::ReadV3000Line(istream& ifs, vector<string>& vs)
  {
    char buffer[BUFF_SIZE];
    if(!ifs.getline(buffer,BUFF_SIZE)) return false;
    tokenize(vs,buffer," \t\n\r");
    if (vs.size() < 2) return false; // timvdm 18/06/2008
    if(vs[0]!="M" || (vs[1]!="V30" && vs[1]!="END")) return false;
        
    if(buffer[strlen(buffer)-1] == '-') //continuation char
      {
        //Read continuation line iteratively and add parsed tokens (without M V30) to vs
        vector<string> vsx;
        if(!ReadV3000Line(ifs,vsx)) return false;
        vs.insert(vs.end(),vsx.begin()+3,vsx.end());
      }
    return true;
  }

  //////////////////////////////////////////////////////
  bool MDLFormat::ReadAtomBlock(istream& ifs,OBMol& mol, OBConversion* pConv)
  {     
    OBAtom atom;
    bool chiralWatch=false;
    int obindex;
    for(obindex=1;;obindex++)
      {
        if(!ReadV3000Line(ifs,vs)) return false;
        if(vs[2]=="END") break;
                
        indexmap[atoi(vs[2].c_str())] = obindex;
        atom.SetVector(atof(vs[4].c_str()), atof(vs[5].c_str()), atof(vs[6].c_str()));
        //      if(abs(atof(vs[6].c_str()))>0)is3D=true;
        //      if(abs(atof(vs[4].c_str()))>0)is2D=true;
        //      if(abs(atof(vs[5].c_str()))>0)is2D=true;
        char type[5];
        strncpy(type,vs[3].c_str(),4);
        int iso=0;
        atom.SetAtomicNum(etab.GetAtomicNum(type,iso));
        if(iso)
          atom.SetIsotope(iso);
        atom.SetType(type); //takes a char not a const char!
        //mapping vs[7] not implemented
                
        //Atom properties
        vector<string>::iterator itr;
        for(itr=vs.begin()+8;itr!=vs.end();itr++)
          {
            string::size_type pos = (*itr).find('=');
            if (pos==string::npos) return false;
            int val = atoi((*itr).substr(pos+1).c_str());

            if((*itr).substr(0,pos)=="CHG")
              {
                atom.SetFormalCharge(val);
              }
            else if((*itr).substr(0,pos)=="RAD")
              {
                atom.SetSpinMultiplicity(val);
              }
            else if((*itr).substr(0,pos)=="CFG")
              {
                //Stereo configuration: 0 none; 1 odd parity; 2 even parity; (3 either parity)
                //Reversed 12Aug05 as advised by Nick England
               /* @todo
                if(val==2) atom.SetAntiClockwiseStereo();
                else if(val==1) atom.SetClockwiseStereo();
                else if(val==3) atom.SetChiral();
                chiralWatch=true;
                */
              }
            else if((*itr).substr(0,pos)=="MASS")
              {
                if(val) atom.SetIsotope(val);
              }
            else if((*itr).substr(0,pos)=="VAL")
              {
                //@todo Abnormal valence: 0 normal;-1 zero
              }
            //Several query properties unimplemented
            //Unknown properties ignored
          }
        if(!mol.AddAtom(atom)) return false;
        /*
        if(chiralWatch)
          _mapcd[mol.GetAtom(mol.NumAtoms())]= new OBChiralData; // fill the map with chrial data for each chiral atom
        */
        atom.Clear();
      }
    return true;
  }

  //////////////////////////////////////////////////////
  bool MDLFormat::ReadBondBlock(istream& ifs,OBMol& mol, OBConversion* pConv)
  {
    for(;;)
      {
        if(!ReadV3000Line(ifs,vs)) return false;
        if(vs[2]=="END") break;
                
        unsigned flag=0;

        int order = atoi(vs[3].c_str());
        if(order==4) order=5;

        int obstart = indexmap[atoi(vs[4].c_str())];
        int obend = indexmap[atoi(vs[5].c_str())];

        vector<string>::iterator itr;
        for(itr=vs.begin()+6;itr!=vs.end();itr++)
          {
            string::size_type pos = (*itr).find('=');
            if (pos==string::npos) return false;
            int val = atoi((*itr).substr(pos+1).c_str());

            if((*itr).substr(0,pos)=="CFG")
              {
                //@todo Bond Configuration 2 or 3D??
                if (val == 1) 
                  {
                    flag |= OB_WEDGE_BOND;
                  }
                else if (val == 3) 
                  {
                    flag |= OB_HASH_BOND;
                  }
              }
          }
        if (!mol.AddBond(obstart,obend,order,flag)) return false;
         
        /*
        // after adding a bond to atom "obstart"
        // search to see if atom is bonded to a chiral atom
        map<OBAtom*,OBChiralData*>::iterator ChiralSearch;
        ChiralSearch = _mapcd.find(mol.GetAtom(obstart));
        if (ChiralSearch!=_mapcd.end())
          {
            (ChiralSearch->second)->AddAtomRef(obend, input);
          }
        // after adding a bond to atom "obend"
        // search to see if atom is bonded to a chiral atom
        ChiralSearch = _mapcd.find(mol.GetAtom(obend));
        if (ChiralSearch!=_mapcd.end())
          {
            (ChiralSearch->second)->AddAtomRef(obstart, input);
          }
        */
      }
    return true;
  }

////////////////////////////////////////////////////////////
  bool MDLFormat::ReadCollectionBlock(istream& ifs,OBMol& mol, OBConversion* pConv)
  {
    //Not currently implemented
    obErrorLog.ThrowError(__FUNCTION__, 
      "COLLECTION blocks are not currently implemented and their contents ae ignored.", obWarning);
    for(;;)
    {
      if(!ReadV3000Line(ifs,vs))
        return false;
      if(vs[2]=="END")
        break;
    }
    return true;
  }

  //////////////////////////////////////////////////////////
  bool MDLFormat::WriteV3000(ostream& ofs,OBMol& mol, OBConversion* pConv)
  {
    // Check to see if there are any untyped aromatic bonds (GetBO == 5)
    // These must be kekulized first
    FOR_BONDS_OF_MOL(b, mol)
      {
        if (b->GetBO() == 5)
          {
            mol.Kekulize();
            break;
          }
      }
  
  
    ofs << "  0  0  0     0  0            999 V3000" << endl; //line 4
    ofs << "M  V30 BEGIN CTAB" <<endl;
    ofs << "M  V30 COUNTS " << mol.NumAtoms() << " " << mol.NumBonds() 
        << " 0 0 " << mol.IsChiral() << endl;
        
    ofs << "M  V30 BEGIN ATOM" <<endl;
    OBAtom *atom;
    int index=1;
    vector<OBAtom*>::iterator i;
    for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
      {  
        ofs     << "M  V30 "
                << index++ << " "
                << etab.GetSymbol(atom->GetAtomicNum()) << " "
                << atom->GetX() << " "
                << atom->GetY() << " "
                << atom->GetZ()
                << " 0";
        if(atom->GetFormalCharge()!=0)
          ofs << " CHG=" << atom->GetFormalCharge();
        if(atom->GetSpinMultiplicity()!=0)
          ofs << " RAD=" << atom->GetSpinMultiplicity();
        /*
        if(atom->IsChiral())
          {
            // MOLV3000 uses 1234 unless an H then 123H
         
            OBChiralData* cd=(OBChiralData*)atom->GetData(OBGenericDataType::ChiralData);
            if(!cd){ //if no Chiral Data Set, need to make one!
              cd=new OBChiralData;
              atom->SetData(cd);
            }
            if (atom->GetHvyValence()==3)
              {
                OBAtom *nbr;
                int Hid = (mol.NumAtoms()+1) ;// max Atom ID +1 
                vector<unsigned int> nbr_atms;
                vector<OBBond*>::iterator i;
                for (nbr = atom->BeginNbrAtom(i);nbr;nbr = atom->NextNbrAtom(i))
                  {
                    if (nbr->IsHydrogen()){Hid=nbr->GetIdx();continue;}
                    nbr_atms.push_back(nbr->GetIdx());
                  }
                sort(nbr_atms.begin(),nbr_atms.end());
                nbr_atms.push_back(Hid);
                cd->SetAtom4Refs(nbr_atms,output);   
              } 
            else if (atom->GetHvyValence()==4)
              {
                vector<unsigned int> nbr_atms;
                int n;
                for(n=1;n<5;n++)nbr_atms.push_back(n);
                cd->SetAtom4Refs(nbr_atms,output); 
              }
            double vol=0;         
            if (mol.HasNonZeroCoords())
              {
                vol=CalcSignedVolume(mol,atom);
                if (vol > 0.0)atom->SetClockwiseStereo();
                else if(vol < 0.0)atom->SetAntiClockwiseStereo();
                CorrectChirality(mol,atom,calcvolume,output);
              }
            else {            
              CorrectChirality(mol,atom); // will set the stereochem based on input/output atom4refs
            }

            int cfg=3; // if we don't know, then it's unspecified
            if(atom->IsClockwise())cfg=1;
            else if(atom->IsAntiClockwise())cfg=2;
                        
            ofs << " CFG=" << cfg;
          }
        */
        if(atom->GetIsotope()!=0)
          ofs << " MASS=" << atom->GetIsotope();
        ofs << endl;
      }
    ofs << "M  V30 END ATOM" <<endl;

    ofs << "M  V30 BEGIN BOND" <<endl;
    //so the bonds come out sorted
    index=1;
    OBAtom *nbr;
    OBBond *bond;
    vector<OBBond*>::iterator j;
    for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
      {
        for (nbr = atom->BeginNbrAtom(j);nbr;nbr = atom->NextNbrAtom(j))
          {
            if (atom->GetIdx() < nbr->GetIdx())
              {
                bond = (OBBond*) *j;
                ofs << "M  V30 "
                    << index++ << " "
                    << bond->GetBO() << " "
                    << bond->GetBeginAtomIdx() << " "
                    << bond->GetEndAtomIdx();
                //@todo do the following stereo chemistry properly
                int cfg=0;
                if(bond->IsWedge()) cfg=1;
                if(bond->IsHash()) cfg=6;
                if(bond->IsWedgeOrHash()) cfg=4;
                if(cfg) ofs << " CFG=" << cfg;
                ofs << endl;
              }
          }
      }
    ofs << "M  V30 END BOND" <<endl;
    ofs << "M  V30 END CTAB" <<endl;
    return true;
  }

  string MDLFormat::GetTimeDate()
  {
    char td[11];
    //returns MMDDYYHHmm
    struct tm* ts;
    time_t long_time;
    time( &long_time );
    ts = localtime( &long_time ); 
    snprintf(td, 11, "%02d%02d%02d%02d%02d", ts->tm_mon+1, ts->tm_mday, 
             ((ts->tm_year>=100)? ts->tm_year-100 : ts->tm_year),
             ts->tm_hour, ts->tm_min);
    return string(td);
  }

  void MDLFormat::GetUpDown(OBMol& mol, map<OBBond*, OBStereo::BondDirection> &updown,
                            set<OBBond*> &stereodbl)
  {
    // IMPROVEME: rewrite to loop over double bonds only

    // Get CisTransStereos
    std::vector<OBCisTransStereo*> cistrans, unvisited_cistrans;
    std::vector<OBGenericData*> vdata = mol.GetAllData(OBGenericDataType::StereoData);
    for (std::vector<OBGenericData*>::iterator data = vdata.begin(); data != vdata.end(); ++data)
      if (((OBStereoBase*)*data)->GetType() == OBStereo::CisTrans) {
        OBCisTransStereo *ct = dynamic_cast<OBCisTransStereo*>(*data);
        cistrans.push_back(ct);
      }
    unvisited_cistrans = cistrans;

    // Initialise two opposite configurations for up/downness
    bool use_alt_config;
    vector<OBStereo::BondDirection> config(4), alt_config(4);
    config[0] = OBStereo::UpBond;   config[3] = config[0];
    config[1] = OBStereo::DownBond; config[2] = config[1];
    alt_config[0] = config[1]; alt_config[3] = alt_config[0];
    alt_config[1] = config[0]; alt_config[2] = alt_config[1];

    // Find bonds in a BFS manner and set their up/downness
    vector<OBCisTransStereo*>::iterator ChiralSearch;
    vector<unsigned long>::iterator lookup; 
    FOR_BONDBFS_OF_MOL(b, mol) {
      if (updown.find(&(*b)) == updown.end()) 
      { // This bond has not yet been set
        vector<OBAtom*> bond_atoms(2);
        bond_atoms[0] = b->GetBeginAtom();
        bond_atoms[1] = b->GetEndAtom();
        // Is this a stereo bond?
        for (vector<OBAtom*>::iterator bond_atom = bond_atoms.begin(); bond_atom!=bond_atoms.end(); ++bond_atom)
        {
          for (ChiralSearch = unvisited_cistrans.begin(); ChiralSearch != unvisited_cistrans.end(); ChiralSearch++)
          {
            OBCisTransStereo::Config cfg = (*ChiralSearch)->GetConfig(OBStereo::ShapeU);
            lookup = std::find(cfg.refs.begin(), cfg.refs.end(), (*bond_atom)->GetId());
            if (lookup != cfg.refs.end() && (cfg.begin == b->GetNbrAtom(*bond_atom)->GetId() ||
                                               cfg.end == b->GetNbrAtom(*bond_atom)->GetId()) )
            { // We have a stereo bond. Now get the five bonds involved...
              OBBond* dbl_bond = mol.GetBond(mol.GetAtomById(cfg.begin), mol.GetAtomById(cfg.end));
              std::vector<OBBond *> refbonds(4, (OBBond*)NULL);
              if (cfg.refs[0] != OBStereo::ImplicitId) // Could be a hydrogen
                refbonds[0] = mol.GetBond(mol.GetAtomById(cfg.refs[0]), mol.GetAtomById(cfg.begin));
              if (cfg.refs[1] != OBStereo::ImplicitId) // Could be a hydrogen
                refbonds[1] = mol.GetBond(mol.GetAtomById(cfg.refs[1]), mol.GetAtomById(cfg.begin));
              
              if (cfg.refs[2] != OBStereo::ImplicitId) // Could be a hydrogen
                refbonds[2] = mol.GetBond(mol.GetAtomById(cfg.refs[2]), mol.GetAtomById(cfg.end));
              if (cfg.refs[3] != OBStereo::ImplicitId) // Could be a hydrogen
                refbonds[3] = mol.GetBond(mol.GetAtomById(cfg.refs[3]), mol.GetAtomById(cfg.end));              
              

              if (cfg.specified) {

                // If any of the bonds have been previously set, now set them all
                // in agreement
                use_alt_config = false;
                for (int i=0; i<4; ++i)
                  if (updown.find(refbonds[i]) != updown.end()) // We have already set this one (conjugated bond)
                    if (updown[refbonds[i]] != config[i])
                    {
                      use_alt_config = true;
                      break;
                    }
              }
              
              // Set the configuration
              stereodbl.insert(dbl_bond);
              if (cfg.specified) {
                for(int i=0;i<4;i++)
                  if (refbonds[i] != NULL)
                    updown[refbonds[i]] = use_alt_config ? alt_config[i] : config[i];
              }
              else { // Cis/Trans unknown
                for(int i=0;i<4;++i)
                  if (updown.find(refbonds[i]) == updown.end())
                    updown[refbonds[i]] = OBStereo::UnknownDir;
              }

              unvisited_cistrans.erase(ChiralSearch);
              break; // Break out of the ChiralSearch (could break out of the outer loop too...let's see)
            }
          }
        }
      }
    }
  }
  void MDLFormat::GetParity(OBMol& mol, map<OBAtom*, MDLFormat::Parity> &parity,
     map<OBBond*, OBStereo::BondDirection> &updown)
  {
    // Get TetrahedralStereos
    std::vector<OBTetrahedralStereo*> tet;
    std::vector<OBGenericData*> vdata = mol.GetAllData(OBGenericDataType::StereoData);
    for (std::vector<OBGenericData*>::iterator data = vdata.begin(); data != vdata.end(); ++data)
      if (((OBStereoBase*)*data)->GetType() == OBStereo::Tetrahedral) {
        OBTetrahedralStereo *ts = dynamic_cast<OBTetrahedralStereo*>(*data);
        
        OBTetrahedralStereo::Config cfg = ts->GetConfig();

        enum Parity atomparity = Unknown;
        if (cfg.specified) {
          // If, when looking towards the maxref, the remaining refs increase in number
          // clockwise, parity is 1 (Parity::Clockwise)
          OBStereo::Refs refs = cfg.refs;
          unsigned long maxref = std::max(*(std::max_element(refs.begin(), refs.end())), cfg.from);
          
          // Get a new cfg and refs looking towards the maxref
          cfg = ts->GetConfig(maxref, OBStereo::Clockwise, OBStereo::ViewTowards);
          int inversions = OBStereo::NumInversions(cfg.refs);
          
          // If they were in increasing order, inversions would be 0 or some even value      
          if (inversions % 2 == 0)
            atomparity = Clockwise;
          else
            atomparity = AntiClockwise;
        }
        parity[mol.GetAtomById(cfg.center)] = atomparity;
        
        // Set Bond Up?
      }

  }
  
}
