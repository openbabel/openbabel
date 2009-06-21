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
#include <openbabel/chiral.h>
#include <openbabel/alias.h>

using namespace std;
namespace OpenBabel
{

  //MDLFormat is a base class which is never instantiated.
  //MOLFormat and SDFormat are derived from it and have their own constructors.
  //SDFormat has its own WriteMolecule which sets the "sd" option.

  class MDLFormat : public OBMoleculeFormat
  {
    map<OBAtom*,OBChiralData*> _mapcd; // map of ChiralAtoms and their data
  public:
    virtual const char* Description()
    { return
        "MDL MOL/SDF format\n"
        "Reads and writes V2000 and V3000 versions\n"
        "Write Options, e.g. -x3\n"
        " 3  output V3000 not V2000 (used for >999 atoms/bonds)\n"
        " m  write no properties\n\n"
        "Input Options, e.g. -av\n"
        " v  AtomValue in V line is NOT a partial charge\n\n";
    };

    virtual const char* SpecificationURL()
    {return "http://www.mdl.com/downloads/public/ctfile/ctfile.jsp";};

    virtual const char* GetMIMEType() 
    { return "chemical/x-mdl-molfile"; };

    virtual unsigned int Flags() { return DEFAULTFORMAT;};
    virtual const char* TargetClassDescription(){return OBMol::ClassDescription();};

    virtual int SkipObjects(int n, OBConversion* pConv)
    {
      if(n==0) n++;
      string temp;
      istream& ifs = *pConv->GetInStream();
      do
        {
          getline(ifs,temp,'$');
          if(ifs.good())
            getline(ifs, temp);
        }while(ifs.good() && temp.substr(0,3)=="$$$" && --n);
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
    bool  HasProperties;
    char* GetTimeDate(char* td);

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
      pConv->AddOption("sd",OBConversion::OUTOPTIONS);
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
    _mapcd.clear();
    bool chiralWatch=false;
    bool setDimension = false;
    stringstream errorMsg;
    string clearError; // empty string to clear the warning buffer

    // Allows addition of further disconnected atoms to an existing molecule
    int offset = mol.NumAtoms(); 

    int i,natoms,nbonds;
    char buffer[BUFF_SIZE];
    string comment;
    string r1,r2;

 	  // Attempting to read past the end of the file -- don't bother
		if( !ifs.good() || ifs.peek() == EOF ) return( false );

    if (!ifs.getline(buffer,BUFF_SIZE)) {
      errorMsg << "WARNING: Problems reading a MDL file\n";
      errorMsg << "Cannot read title line\n";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
      return(false);
    }
    mol.SetTitle(buffer);
  
    if (!ifs.getline(buffer,BUFF_SIZE)) {
      errorMsg << "WARNING: Problems reading a MDL file\n";
      errorMsg << "Cannot read creator/dimension line line\n";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
      return(false);
    }
    
    if (strlen(buffer) > 20) {
      char* dimension = buffer+20;
      dimension[2]='\0'; //truncate after 2D
      if(strcmp(dimension,"3D") == 0)
        {
          mol.SetDimension(3);
          setDimension = true;
        }
    }

    if (!ifs.getline(buffer,BUFF_SIZE)) {
      errorMsg << "WARNING: Problems reading a MDL file\n";
      errorMsg << "Cannot read comment line\n";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
      return(false);
    }

    if (strlen(buffer) > 0) {
      comment = buffer;
    }

    if (!ifs.getline(buffer,BUFF_SIZE)) {
      errorMsg << "WARNING: Problems reading a MDL file\n";
      errorMsg << "Cannot read atom and bond count\n";
      errorMsg << "File ended prematurely\n";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
      return(false);
    }
    r1 = buffer;
    if (r1.size() < 6) { // error from Joe Bedell, Sigma-Aldrich
      errorMsg << "WARNING: Problems reading a MDL file\n";
      errorMsg << "Cannot read atom and bond count\n";
      errorMsg << "Expected standard 6 character atom and bond count\n";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
      return(false);
    }

    natoms = atoi((r1.substr(0,3)).c_str());
    nbonds = atoi((r1.substr(3,3)).c_str());

    mol.BeginModify();
    if(r1.find("V3000")!=string::npos)
      {
        indexmap.clear();
        if(!ReadV3000Block(ifs,mol,pConv,false)) return false;
        //              ifs.getline(buffer,BUFF_SIZE); //M END line
      }
    else
      {
        mol.ReserveAtoms(natoms);
        double x,y,z;
        char type[8];
        vector3 v;
        OBAtom atom;
        int charge, scanArgs, stereo;

        for (i = 0;i < natoms;i++) {
          if (!ifs.getline(buffer,BUFF_SIZE)) {
            errorMsg << "WARNING: Problems reading a MDL file\n";
            errorMsg << "Not enough atoms to match atom count (" << natoms << ") in counts line\n";
            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
            return(false);
          }

          scanArgs = sscanf(buffer,"%lf %lf %lf %5s %*d %d %d",&x,&y,&z,type,&charge, &stereo);
          if (scanArgs <4) {
            errorMsg << "WARNING: Problems reading a MDL file\n";
            errorMsg << "Missing data following atom specification in atom block\n";
            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
            return(false);
          }
          v.SetX(x);v.SetY(y);v.SetZ(z);
          atom.SetVector(x, y, z);
          int iso=0;
          atom.SetAtomicNum(etab.GetAtomicNum(type,iso));
          //                    atom.SetType(type);
          if(iso)
            atom.SetIsotope(iso);

          if (scanArgs >= 5)
            {
              switch (charge)
                {
                case 0: break;
                case 3: atom.SetFormalCharge(1); break;
                case 2: atom.SetFormalCharge(2); break;
                case 1: atom.SetFormalCharge(3); break;
                case 5: atom.SetFormalCharge(-1); break;
                case 6: atom.SetFormalCharge(-2); break;
                case 7: atom.SetFormalCharge(-3); break;
                }
            }

          if (scanArgs == 6) // set a stereo mark
            {
              //Stereo configuration: 0 none; 1 odd parity; 2 even parity; 3 unspecified)
              if (stereo == 2)
                {
                  chiralWatch=true;
                  atom.SetAntiClockwiseStereo();
                }
              else if (stereo == 1)
                {
                  chiralWatch=true;
                  atom.SetClockwiseStereo();
                }
              else if(stereo == 3)
                {
                  chiralWatch=true;
                  atom.SetChiral();
                }
            }

          if (!mol.AddAtom(atom))
            return(false);
          if(chiralWatch)  // fill the map with data for each chiral atom
            _mapcd[mol.GetAtom(mol.NumAtoms())] = new OBChiralData;
          atom.Clear();
        }

        unsigned int start,end,order,flag;
        for (i = 0;i < nbonds;i++) {
          flag = 0;
          if (!ifs.getline(buffer,BUFF_SIZE)) {
            errorMsg << "WARNING: Problems reading a MDL file\n";
            errorMsg << "Not enough bonds to match bond count (" << nbonds << ") in counts line\n";
            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
            return(false);
          }
          r1 = buffer;
          start = end = order = 0;
          if (r1.size() >= 9) {
            start = atoi((r1.substr(0,3)).c_str());
            end = atoi((r1.substr(3,3)).c_str());
            order = atoi((r1.substr(6,3)).c_str());
          }
          if (start == 0 || end == 0 || order == 0 ||
              start > mol.NumAtoms() || end > mol.NumAtoms()) {
            errorMsg << "WARNING: Problems reading a MDL file\n";
            errorMsg << "Invalid bond specification, atom numbers or bond order are wrong.\n";
            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
            return false;
          }

          order = (order == 4) ? 5 : order;
          if (r1.size() >= 12) {  //handle wedge/hash data
            stereo = atoi((r1.substr(9,3)).c_str());
            if (stereo) {
              if (stereo == 1) flag |= OB_WEDGE_BOND;
              if (stereo == 6) flag |= OB_HASH_BOND;
            }
          }

          if (!mol.AddBond(start+offset,end+offset,order,flag)) {
            errorMsg << "WARNING: Problems reading a MDL file\n";
            errorMsg << "Invalid bond specification\n";
            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
            return(false);
          }

          // after adding a bond to atom # "start+offset"
          // search to see if atom is bonded to a chiral atom
          map<OBAtom*,OBChiralData*>::iterator ChiralSearch;
          ChiralSearch = _mapcd.find(mol.GetAtom(start+offset));
          if (ChiralSearch!=_mapcd.end())
            {
              (ChiralSearch->second)->AddAtomRef(end+offset, input);
            }
          // after adding a bond to atom # "end + offset"
          // search to see if atom is bonded to a chiral atom
          ChiralSearch = _mapcd.find(mol.GetAtom(end+offset));
          if (ChiralSearch!=_mapcd.end())
            {
              (ChiralSearch->second)->AddAtomRef(start+offset, input);
            }
        }

        //Read Properties block, currently only A, V, S, M RAD, M CHG, M ISO

        while(ifs.getline(buffer,BUFF_SIZE))
          {
            if(!strncmp(buffer,"$$$$",4))
              return true;
            if(!strncmp(buffer,"M  END",6))
              break;
            if(!strncmp(buffer,"S  SKP",6))
              {
                int i = atoi(buffer+6);
                for(;i>0;--i)
                  if (ifs.good()) // check for EOL, suggested by Dalke
                    ifs.getline(buffer,BUFF_SIZE);
                continue;
              }
            if(buffer[0]=='A') //alias
              {
                int atomnum = atoi(buffer+2);
                ifs.getline(buffer,BUFF_SIZE);
                if(*buffer!='?' && *buffer!='*')
                  {
                    AliasData* ad = new AliasData();
                    ad->SetAlias(buffer);
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

            if(buffer[0]=='V') //Atom Value
              {
                //Put contents into a OBPairData with attribute "AtomValue" attached to the OBAtom,
                //and, if the -av option is NOT set, interpret as a partial charge.
                int atomnum = atoi(buffer+2);
                string val(buffer+7);
                Trim(val);
                OBAtom* at = mol.GetAtom(atomnum);
                if(at) {
                  OBPairData *dp = new OBPairData;
                  dp->SetAttribute("AtomValue");
                  dp->SetValue(val);
                  dp->SetOrigin(fileformatInput);
                  at->SetData(dp);
                  if(!pConv->IsOption("v",OBConversion::INOPTIONS)) {
                    double pchg = atof(val.c_str());
                    if(fabs(pchg)<4.0)
                      at->SetPartialCharge(pchg);
                  }
                }
              }

            if(strncmp(buffer,"M  CHG",6) && strncmp(buffer,"M  RAD",6) && strncmp(buffer,"M  ISO",6))
              continue;
            r1 = buffer;
            int n = -1;
            if (r1.size() >= 9)
              n = atoi((r1.substr(6,3)).c_str()); //entries on this line
            if(n<=0 || n>8 || 6+n*8>r1.size()) //catch ill-formed line
              {
                obErrorLog.ThrowError(__FUNCTION__, "Error in line:\n" + r1, obError);
                return false;
              }
            int pos = 10;
            for(;n>0;n--,pos+=8)
              {
                int atomnumber = atoi((r1.substr(pos,3)).c_str());
                OBAtom* at;
                if (atomnumber==0 || (at=mol.GetAtom(atomnumber+offset))==NULL)
                  {
                    obErrorLog.ThrowError(__FUNCTION__, "Error in line:\n" + r1, obError);
                    return false;
                  }

                at=mol.GetAtom(atomnumber+offset); //atom numbers start at 1
                int value = atoi((r1.substr(pos+4,3)).c_str());
                if(r1.substr(3,3)=="RAD")
                  at->SetSpinMultiplicity(value);
                else if(r1.substr(3,3)=="CHG")
                  at->SetFormalCharge(value);
                else if(r1.substr(3,3)=="ISO")
                  at->SetIsotope(value);
                //Although not done here,according to the specification, 
                //previously set formal charges should be reset to zero
              }
            // Lines setting several other properties are not implemented
          }
      }
    mol.AssignSpinMultiplicity();

    mol.EndModify();
        
    //NE add the OBChiralData stored inside the _mapcd to the atoms now after end
    // modify so they don't get lost.
    if(_mapcd.size()>0)
      {
        OBAtom* atom;
        OBChiralData* cd;
        map<OBAtom*,OBChiralData*>::iterator ChiralSearch;
        for(ChiralSearch=_mapcd.begin();ChiralSearch!=_mapcd.end();ChiralSearch++)
          {
            atom=ChiralSearch->first;
            cd=ChiralSearch->second;
            atom->SetData(cd);
          }    
      }

    if (comment.length())
      {
        OBCommentData *cd = new OBCommentData;
        cd->SetData(comment);
        cd->SetOrigin(fileformatInput);
        mol.SetData(cd);
      }
        
    //Get property lines
    while (ifs.getline(buffer,BUFF_SIZE))
      {
        if (strstr(buffer,"<"))
          {
            string buff(buffer);
            size_t lt=buff.find("<")+1;
            size_t rt = buff.find_last_of(">");
            string attr = buff.substr(lt,rt-lt);

            // sometimes we can hit more data than BUFF_SIZE, so we'll use a std::string
            string line;
            buff.clear();
            while (getline(ifs, line))
              {
                Trim(line);
                if (line.size())
                  {
                    buff.append(line);
                    buff += "\n";
                  }
                else
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
        if (!strncmp(buffer,"$$$$",4)) break;
        if (!strncmp(buffer,"$MOL",4)) break;
      }

    if (!setDimension && mol.Has3D())
      mol.SetDimension(3);
    else if (!setDimension && !mol.Has3D())
      mol.SetDimension(2);

    return(true);

  }

  /////////////////////////////////////////////////////////////////
  bool MDLFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    if(pConv->GetOutputIndex()==1)
      HasProperties=false;

    char dimension[3] = "2D";
    if(mol.GetDimension()==3)
      dimension[0]='3';

    mol.FindChiralCenters(); // Needed to mark centers as chiral from formats like xyz

    ofs << mol.GetTitle() <<  endl; //line 1

    char td[11];
    ofs << " OpenBabel" << GetTimeDate(td) <<  dimension << endl; //line2

    if (mol.HasData(OBGenericDataType::CommentData))
      {
        OBCommentData *cd = (OBCommentData*)mol.GetData(OBGenericDataType::CommentData);
        ofs << cd->GetData() << endl; //line 3
      }
    else
      ofs << endl;
        
    if(pConv->IsOption("3") || mol.NumAtoms() > 999 || mol.NumBonds() > 999)
      {
        if(!WriteV3000(ofs,mol,pConv)) return false;
      }

    else
      {
        //The rest of the function is the same as the original
        char buff[BUFF_SIZE];  

        if (mol.NumAtoms() > 999 || mol.NumBonds() > 999) // Three digits!
          {
            stringstream errorMsg;
            errorMsg << "MDL Molfile conversion failed: Molecule is too large to convert." << endl;
            errorMsg << "  File format (v2000) is limited to 999 atoms or bonds." << endl;
            errorMsg << "  Molecule size: " << mol.NumAtoms() << " atoms ";
            errorMsg << "and " << mol.NumBonds() << " bonds." << endl;
            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
            //      delete pOb;
            return(false);
          }

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

        snprintf(buff, BUFF_SIZE, "%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d V2000",
                 mol.NumAtoms(),mol.NumBonds(),0,0,0,0,0,0,0,0,999);
        ofs << buff << endl;

        OBAtom *atom;
        vector<OBAtom*>::iterator i;
        int charge=0; //Charge is now output in a M CHG statement
        for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i)) {
          /*          switch (atom->GetFormalCharge()) {
                      case 1: charge = 3; break;
                      case 2: charge = 2; break;
                      case 3: charge = 1; break;
                      case -1: charge = 5; break;
                      case -2: charge = 6; break;
                      case -3: charge = 7; break;
                      default:
                      charge=0; break;
                      }
          */
          snprintf(buff, BUFF_SIZE, "%10.4f%10.4f%10.4f %-3s%2d%3d%3d%3d%3d",
                   atom->GetX(),
                   atom->GetY(),
                   atom->GetZ(),
                   atom->GetAtomicNum() ? etab.GetSymbol(atom->GetAtomicNum()) : "* ",
                   0,charge,0,0,0);    
          ofs << buff << endl;
        }

        //so the bonds come out sorted
        OBAtom *nbr;
        OBBond *bond;
        vector<OBBond*>::iterator j;
        for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
          for (nbr = atom->BeginNbrAtom(j);nbr;nbr = atom->NextNbrAtom(j))
            if (atom->GetIdx() < nbr->GetIdx()) {
              bond = (OBBond*) *j;
                                        
              int stereo=0; //21Jan05 CM
              if(strcmp(dimension,"2D")==0)
                {
                  int flag = bond->GetFlags();
                  if (flag & OB_WEDGE_BOND) stereo=1;
                  if (flag & OB_HASH_BOND ) stereo=6;
                }
              snprintf(buff, BUFF_SIZE, "%3d%3d%3d%3d%3d%3d",
                       bond->GetBeginAtomIdx(),
                       bond->GetEndAtomIdx(),
                       bond->GetBO(),
                       stereo,0,0);
              ofs << buff << endl;
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
        //      if(fabs(atof(vs[6].c_str()))>0)is3D=true;
        //      if(fabs(atof(vs[4].c_str()))>0)is2D=true;
        //      if(fabs(atof(vs[5].c_str()))>0)is2D=true;
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
                if(val==2) atom.SetAntiClockwiseStereo();
                else if(val==1) atom.SetClockwiseStereo();
                else if(val==3) atom.SetChiral();
                chiralWatch=true;
              }
            else if((*itr).substr(0,pos)=="MASS")
              {
                if(val) atom.SetIsotope(val);
              }
            else if((*itr).substr(0,pos)=="VAL")
              {
                //TODO Abnormal valence: 0 normal;-1 zero
              }
            //Several query properties unimplemented
            //Unknown properties ignored
          }
        if(!mol.AddAtom(atom)) return false;
        if(chiralWatch)_mapcd[mol.GetAtom(mol.NumAtoms())]= new OBChiralData; // fill the map with chrial data for each chiral atom
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
                //TODO Bond Configuration 2 or 3D??
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
                //TODO do the following stereo chemistry properly
                int cfg=0;
                if(bond->IsWedge()) cfg=1;
                if(bond->IsHash()) cfg=3;
                if(cfg) ofs << " CFG=" << cfg;
                ofs << endl;
              }
          }
      }
    ofs << "M  V30 END BOND" <<endl;
    ofs << "M  V30 END CTAB" <<endl;
    return true;
  }

  char* MDLFormat::GetTimeDate(char* td)
  {
    //returns MMDDYYHHmm
    struct tm* ts;
    time_t long_time;
    time( &long_time );
    ts = localtime( &long_time ); 
    snprintf(td, 11, "%02d%02d%02d%02d%02d", ts->tm_mon+1, ts->tm_mday, 
             ((ts->tm_year>=100)? ts->tm_year-100 : ts->tm_year),
             ts->tm_hour, ts->tm_min);
    return td;
  }


}
