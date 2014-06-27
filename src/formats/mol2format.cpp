/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2007 by Geoffrey R. Hutchison
Some portions Copyright (C) 2004 by Chris Morley

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

using namespace std;
namespace OpenBabel
{

  //The routine WriteSmiOrderedMol2() in the original mol2.cpp is presumably
  //another output format, but was not made available in version 100.1.2, nor
  //is it here.

  class MOL2Format : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    MOL2Format()
    {
      OBConversion::RegisterFormat("mol2",this, "chemical/x-mol2");
      OBConversion::RegisterFormat("ml2",this);
      OBConversion::RegisterFormat("sy2",this);
      OBConversion::RegisterOptionParam("l", NULL, 0, OBConversion::OUTOPTIONS);
    }

    virtual const char* Description() //required
    {
      return
        "Sybyl Mol2 format\n"
        "Write Options e.g. -xl\n"
        "  l               Output ignores residue information (only ligands)\n\n";
    };

    virtual const char* SpecificationURL()
    {
      return "http://www.tripos.com/data/support/mol2.pdf";
    }; //optional

    virtual const char* GetMIMEType()
    { return "chemical/x-mol2"; };

    virtual int SkipObjects(int n, OBConversion* pConv);

    //*** This section identical for most OBMol conversions ***
    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
  };
  //***

  //Make an instance of the format class
  MOL2Format theMOL2Format;

  /////////////////////////////////////////////////////////////////
  bool MOL2Format::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;

    //Old code follows...
    bool foundAtomLine = false;
    char buffer[BUFF_SIZE];
    char *comment = NULL;
    string str,str1;
    vector<string> vstr;
    int len;

    mol.BeginModify();

    for (;;)
      {
        if (!ifs.getline(buffer,BUFF_SIZE))
          return(false);
        if (EQn(buffer,"@<TRIPOS>MOLECULE",17))
          break;
      }

    // OK, just read MOLECULE line
    int lcount;
    int natoms,nbonds;
    bool hasPartialCharges = true;
    for (lcount=0;;lcount++)
      {
        if (!ifs.getline(buffer,BUFF_SIZE))
          return(false);
        if (EQn(buffer,"@<TRIPOS>ATOM",13))
          {
            foundAtomLine = true;
            break;
          }

        if (lcount == 0)
          {
            tokenize(vstr,buffer);
            if (!vstr.empty())
              mol.SetTitle(buffer);
          }
        else if (lcount == 1)
          sscanf(buffer,"%d%d",&natoms,&nbonds);
        else if (lcount == 3) // charge descriptions
          {
            // Annotate origin of partial charges
            OBPairData *dp = new OBPairData;
            dp->SetAttribute("PartialCharges");
            dp->SetValue(buffer);
            dp->SetOrigin(fileformatInput);
            mol.SetData(dp);

            if (strncasecmp(buffer, "NO_CHARGES", 10) == 0)
              hasPartialCharges = false;
          }
        else if (lcount == 4) //energy (?)
          {
            tokenize(vstr,buffer);
            if (!vstr.empty() && vstr.size() == 3)
              if (vstr[0] == "Energy")
                mol.SetEnergy(atof(vstr[2].c_str()));
          }
        else if (lcount == 5) //comment
          {
            if ( buffer[0] )
              {
                len = (int) strlen(buffer)+1;
                //! @todo allow better multi-line comments
                // which don't allow ill-formed data to consume memory
                // Thanks to Andrew Dalke for the pointer
                if (comment != NULL)
                  delete [] comment;
                comment = new char [len];
                memcpy(comment,buffer,len);
              }
          }
      }

    if (!foundAtomLine)
      {
        mol.EndModify();
        mol.Clear();
        obErrorLog.ThrowError(__FUNCTION__, "Unable to read Mol2 format file. No atoms found.", obWarning);
        return(false);
      }

    mol.ReserveAtoms(natoms);

    int i;
    vector3 v;
    OBAtom atom;
    double x,y,z,pcharge;
    char temp_type[BUFF_SIZE], resname[BUFF_SIZE], atmid[BUFF_SIZE];
    int elemno, resnum = -1;

    ttab.SetFromType("SYB");
    for (i = 0;i < natoms;i++)
      {
        if (!ifs.getline(buffer,BUFF_SIZE))
          return(false);
        sscanf(buffer," %*s %1024s %lf %lf %lf %1024s %d %1024s %lf",
               atmid, &x,&y,&z, temp_type, &resnum, resname, &pcharge);

        atom.SetVector(x, y, z);
        atom.SetFormalCharge(0);

        // Handle "CL" and "BR" and other mis-typed atoms
        str = temp_type;
        if (strncmp(temp_type, "CL", 2) == 0) {
          str = "Cl";
        } else  if (strncmp(temp_type,"BR",2) == 0) {
          str = "Br";
        } else if (strncmp(temp_type,"S.o2", 4) == 0) {
          str = "S.O2";
        } else if (strncmp(temp_type,"S.o", 3) == 0) {
          str = "S.O";
        } else if (strncmp(temp_type,"SI", 2) == 0) {
          str = "Si";
        // The following cases are entries which are not in openbabel/data/types.txt
        // and should probably be added there
        } else if (strncmp(temp_type,"S.1", 3) == 0) {
          str = "S.2"; // no idea what the best type might be here
        } else if (strncmp(temp_type,"P.", 2) == 0) {
          str = "P.3";
        } else if (strncasecmp(temp_type,"Ti.", 3) == 0) { // e.g. Ti.th
          str = "Ti";
        } else if (strncasecmp(temp_type,"Ru.", 3) == 0) { // e.g. Ru.oh
          str = "Ru";
        // Fixes PR#3557898
        } else if (strncmp(temp_type, "N.4", 3) == 0) {
          atom.SetFormalCharge(1);
        }

        ttab.SetToType("ATN");
        ttab.Translate(str1,str);
        elemno = atoi(str1.c_str());
        ttab.SetToType("IDX");

        // We might have missed some SI or FE type things above, so here's
        // another check
        if( !elemno && isupper(temp_type[1]) )
          {
            temp_type[1] = (char)tolower(temp_type[1]);
            str = temp_type;
            ttab.SetToType("ATN");
            ttab.Translate(str1,str);
            elemno = atoi(str1.c_str());
            ttab.SetToType("IDX");
          }
        // One last check if there isn't a period in the type,
        // it's a malformed atom type, but it may be the element symbol
        // GaussView does this (PR#1739905)
        if ( !elemno ) {
          std::stringstream errorMsg;
          errorMsg << "This Mol2 file is non-standard. Problem with molecule: "
                   << mol.GetTitle()
                   << " Cannot interpret atom types correctly, instead attempting to interpret atom type: "
                   << str << " as elements instead.";
          obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);

          string::size_type dotPos = str.find('.');
          if (dotPos == string::npos) {
            elemno = etab.GetAtomicNum(str.c_str());
          }
        }

        atom.SetAtomicNum(elemno);
        ttab.SetToType("INT");
        ttab.Translate(str1,str);
        atom.SetType(str1);
        atom.SetPartialCharge(pcharge);
        // MMFF94 has different atom types for Cu(I) and Cu(II)
        // as well as for Fe(II) and Fe(III), so the correct formal
        // charge is needed for correct atom type assignment
        if (str1 == "Cu" || str1 == "Fe")
          atom.SetFormalCharge((int)pcharge);
        if (!mol.AddAtom(atom))
          return(false);
        if (!IsNearZero(pcharge))
          hasPartialCharges = true;

        // Add residue information if it exists
        if (resnum != -1 && resnum != 0 &&
            strlen(resname) != 0 && strncmp(resname,"<1>", 3) != 0)
          {
            OBResidue *res  = (mol.NumResidues() > 0) ?
              mol.GetResidue(mol.NumResidues()-1) : NULL;
            if (res == NULL || res->GetName() != resname ||
                static_cast<int>(res->GetNum()) != resnum)
              {
                vector<OBResidue*>::iterator ri;
                for (res = mol.BeginResidue(ri) ; res ; res = mol.NextResidue(ri))
                  if (res->GetName() == resname &&
                      static_cast<int>(res->GetNum()) == resnum)
                    break;

                if (res == NULL)
                  {
                    res = mol.NewResidue();
                    res->SetName(resname);
                    res->SetNum(resnum);
                  }
              }
            OBAtom *atomPtr = mol.GetAtom(mol.NumAtoms());
            res->AddAtom(atomPtr);
            res->SetAtomID(atomPtr, atmid);
          } // end adding residue info
      }

    for (;;)
      {
        if (!ifs.getline(buffer,BUFF_SIZE))
          return(false);
        str = buffer;
        if (!strncmp(buffer,"@<TRIPOS>BOND",13))
          break;
      }

    int start,end,order;
    for (i = 0; i < nbonds; i++)
      {
        if (!ifs.getline(buffer,BUFF_SIZE))
          return(false);

        sscanf(buffer,"%*d %d %d %1024s",&start,&end,temp_type);
        str = temp_type;
        order = 1;
        if (str == "ar" || str == "AR" || str == "Ar")
          order = 5;
        else if (str == "AM" || str == "am" || str == "Am")
          order = 1;
        else
          order = atoi(str.c_str());

        mol.AddBond(start,end,order);
      }

    // Make a pass to ensure that there are no double bonds
    // between atoms which are also involved in aromatic bonds
    // as that may ill-condition kekulization (fixes potential
    // issues with molecules like CEWYIM30 (MMFF94 validation suite)
    // Patch by Paolo Tosco 2012-06-07
    int idx1, idx2;
    bool idx1arom, idx2arom;
    FOR_BONDS_OF_MOL(bond, mol) {
      if (bond->GetBO() != 2)
          continue;
      idx1 = bond->GetBeginAtom()->GetIdx();
      idx2 = bond->GetEndAtom()->GetIdx();
      idx1arom = idx2arom = false;
      FOR_BONDS_OF_MOL(bond2, mol) {
        if (&*bond == &*bond2)
          continue;
        if ((bond2->GetBeginAtom()->GetIdx() == idx1 || bond2->GetEndAtom()->GetIdx() == idx1)
          && bond2->GetBO() == 5)
          idx1arom = true;
        else if ((bond2->GetBeginAtom()->GetIdx() == idx2 || bond2->GetEndAtom()->GetIdx() == idx2)
          && bond2->GetBO() == 5)
          idx2arom = true;
        if (idx1arom && idx2arom) {
          bond->SetBO(1);
          break;
        }
      }
    }

    // Now that bonds are added, make a pass to "de-aromatize" carboxylates
    // and (di)thiocarboxylates
    // Fixes PR#3092368
    OBAtom *carboxylCarbon, *oxysulf;
    FOR_BONDS_OF_MOL(bond, mol)
      {
        if (bond->GetBO() != 5)
          continue;

        if (bond->GetBeginAtom()->IsCarboxylOxygen() || bond->GetBeginAtom()->IsThiocarboxylSulfur()) {
          carboxylCarbon = bond->GetEndAtom();
          oxysulf = bond->GetBeginAtom();
        } else if (bond->GetEndAtom()->IsCarboxylOxygen() || bond->GetEndAtom()->IsThiocarboxylSulfur()) {
          carboxylCarbon = bond->GetBeginAtom();
          oxysulf = bond->GetEndAtom();
        } else // not a carboxylate
          continue;

        if (carboxylCarbon->HasDoubleBond()) { // we've already picked a double bond
          bond->SetBO(1); // this should be a single bond, not "aromatic"
          continue;
        }

        // We need to choose a double bond
        if (oxysulf->ExplicitHydrogenCount() == 1 || oxysulf->GetFormalCharge() == -1) { // single only
          bond->SetBO(1);
          continue;
        } else
          bond->SetBO(2); // we have to pick one, let's use this one
      }

    // Make a pass to fix aromatic bond orders and formal charges
    // involving nitrogen and oxygen atoms - before this patch
    // the aromaticity of a molecule as simple as pyridinium
    // cation could not be correctly perceived
    // Patch by Paolo Tosco 2012-06-07
    OBAtom *carbon, *partner, *boundToNitrogen;
    OBBitVec bv;

    bv.SetBitOn(nbonds);
    bv.Clear();
    FOR_BONDS_OF_MOL(bond, mol)
    {
      if (bv[bond->GetIdx()] || (bond->GetBO() != 5))
        continue;

      // only bother for 6 membered rings (e.g., pyridinium)
      // 5-membered rings like pyrrole, imidazole, or triazole are OK with nH
      if ( (bond->FindSmallestRing())->Size() != 6 )
        continue;

      if ((bond->GetBeginAtom()->IsCarbon() && bond->GetEndAtom()->IsNitrogen())
        || (bond->GetBeginAtom()->IsNitrogen() && bond->GetEndAtom()->IsCarbon())) {
        carbon = (bond->GetBeginAtom()->IsCarbon() ? bond->GetBeginAtom() : bond->GetEndAtom());
        int min_n_h_bonded = 100;
        int min_idx = mol.NumAtoms() + 1;
        FOR_BONDS_OF_ATOM(bond2, carbon) {
          if (bond2->GetBO() != 5)
            continue;
          partner = (bond2->GetBeginAtom() == carbon ? bond2->GetEndAtom() : bond2->GetBeginAtom());
          if (partner->IsNitrogen() && partner->GetValence() == 3 && partner->GetFormalCharge() == 0) {
            int n_h_bonded = 0;
            FOR_BONDS_OF_ATOM(bond3, partner) {
              boundToNitrogen = (bond3->GetBeginAtom() == partner ? bond3->GetEndAtom() : bond3->GetBeginAtom());
              if (boundToNitrogen->IsHydrogen())
                n_h_bonded++;
            }
            if (n_h_bonded < min_n_h_bonded || (n_h_bonded == min_n_h_bonded && partner->GetIdx() < min_idx)) {
              min_n_h_bonded = n_h_bonded;
              min_idx = partner->GetIdx();
            }
          }
        }
        FOR_BONDS_OF_ATOM(bond2, carbon) {
          if (bond2->GetBO() != 5)
            continue;
          partner = (bond2->GetBeginAtom() == carbon ? bond2->GetEndAtom() : bond2->GetBeginAtom());
          if (partner->IsNitrogen() && partner->GetValence() == 3 && partner->GetFormalCharge() == 0) {
            int n_ar_bond = 0;
            FOR_BONDS_OF_ATOM(bond3, partner) {
              boundToNitrogen = (bond3->GetBeginAtom() == partner ? bond3->GetEndAtom() : bond3->GetBeginAtom());
              if (boundToNitrogen->IsOxygen() && boundToNitrogen->GetValence() == 1) {
                n_ar_bond = -1;
                break;
              }
              if (bond3->GetBO() == 5)
                ++n_ar_bond;
            }
            if (n_ar_bond == -1)
              continue;
            if (partner->GetIdx() == min_idx) {
              partner->SetFormalCharge(1);
              if (n_ar_bond == 1) {
                bond2->SetBO(2);
              }
            }
            else if (n_ar_bond == 1) {
              bond2->SetBO(1);
            }
          }
          bv.SetBitOn(bond2->GetIdx());
        }
      } else if ((bond->GetBeginAtom()->IsCarbon() && bond->GetEndAtom()->IsOxygen())
        || (bond->GetBeginAtom()->IsOxygen() && bond->GetEndAtom()->IsCarbon())) {
        OBAtom *atom1, *atom2;
        atom1 = bond->GetBeginAtom();
        atom2 = bond->GetEndAtom();
        // set formal charges for pyrilium
        // (i.e., this bond is a 6-membered ring, aromatic, and C-O)
        if (atom1->IsOxygen() && atom1->IsInRingSize(6))
          atom1->SetFormalCharge(1);
        else if (atom2->IsOxygen() && atom2->IsInRingSize(6))
          atom2->SetFormalCharge(1);
      }
    }
    // Suggestion by Liu Zhiguo 2008-01-26
    // Mol2 files define atom types -- there is no need to re-perceive
    mol.SetAtomTypesPerceived();
    mol.EndModify();

    //must add generic data after end modify - otherwise it will be blown away
    if (comment)
      {
        OBCommentData *cd = new OBCommentData;
        cd->SetData(comment);
        cd->SetOrigin(fileformatInput);
        mol.SetData(cd);
        delete [] comment;
        comment = NULL;
      }
    if (hasPartialCharges)
      mol.SetPartialChargesPerceived();

    /* Disabled due to PR#3048758 -- seekg is very slow with gzipped mol2
    // continue untill EOF or untill next molecule record
    streampos pos;
    for(;;)
      {
        pos = ifs.tellg();
        if (!ifs.getline(buffer,BUFF_SIZE))
          break;
        if (EQn(buffer,"@<TRIPOS>MOLECULE",17))
          break;
      }

    ifs.seekg(pos); // go back to the end of the molecule
    */

    return(true);
  }

  ////////////////////////////////////////////////////////////////

  bool MOL2Format::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;
    bool ligandsOnly = pConv->IsOption("l", OBConversion::OUTOPTIONS)!=NULL;

    //The old code follows....
    string str,str1;
    char buffer[BUFF_SIZE],label[BUFF_SIZE];
    char rnum[BUFF_SIZE],rlabel[BUFF_SIZE];

    ofs << "@<TRIPOS>MOLECULE" << endl;
    str = mol.GetTitle();
    if (str.empty())
      ofs << "*****" << endl;
    else
      ofs << str << endl;

    snprintf(buffer, BUFF_SIZE," %d %d 0 0 0", mol.NumAtoms(),mol.NumBonds());
    ofs << buffer << endl;
    ofs << "SMALL" << endl; // TODO: detect if we have protein, biopolymer, etc.

    OBPairData *dp = (OBPairData*)mol.GetData("PartialCharges");
    if (dp != NULL) {
        // Tripos spec says:
        // NO_CHARGES, DEL_RE, GASTEIGER, GAST_HUCK, HUCKEL, PULLMAN,
        // GAUSS80_CHARGES, AMPAC_CHARGES, MULLIKEN_CHARGES, DICT_ CHARGES,
        // MMFF94_CHARGES, USER_CHARGES
      if (strcasecmp(dp->GetValue().c_str(),"Mulliken") == 0)
        ofs << "MULLIKEN_CHARGES" << endl;
      else if (strcasecmp(dp->GetValue().c_str(),"MMFF94") == 0)
        ofs << "MMFF94_CHARGES" << endl;
      else if (strcasecmp(dp->GetValue().c_str(),"ESP") == 0)
        ofs << "USER_CHARGES" << endl;
      else if (strcasecmp(dp->GetValue().c_str(),"Gasteiger") == 0)
        ofs << "GASTEIGER" << endl;
      else // ideally, code should pick from the Tripos types
        ofs << "USER_CHARGES" << endl;
    }
    else { // No idea what these charges are... all our code sets "PartialCharges"
        ofs << "GASTEIGER" << endl;
    }

    //    ofs << "Energy = " << mol.GetEnergy() << endl;

    if (mol.HasData(OBGenericDataType::CommentData))
      {
        ofs << "****\n"; // comment line printed, so we need to add "no status bits set"
        OBCommentData *cd = (OBCommentData*)mol.GetData(OBGenericDataType::CommentData);
        ofs << cd->GetData();
      }

    ofs << endl;
    ofs << "@<TRIPOS>ATOM" << endl;

    OBAtom *atom;
    OBResidue *res;

    vector<OBAtom*>::iterator i;
    vector<int> labelcount;
    labelcount.resize( etab.GetNumberOfElements() );

    ttab.SetFromType("INT");
    ttab.SetToType("SYB");

    for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
      {

        //
        //  Use sequentially numbered atom names if no residues
        //

        snprintf(label,BUFF_SIZE, "%s%d",
                 etab.GetSymbol(atom->GetAtomicNum()),
                 ++labelcount[atom->GetAtomicNum()]);
        strcpy(rlabel,"<1>");
        strcpy(rnum,"1");

        str = atom->GetType();
        ttab.Translate(str1,str);

        //
        //  Use original atom names if there are residues
        //

        if (!ligandsOnly && (res = atom->GetResidue()) )
          {
            // use original atom names defined by residue
            snprintf(label,BUFF_SIZE,"%s",(char*)res->GetAtomID(atom).c_str());
            // make sure that residue name includes its number
            snprintf(rlabel,BUFF_SIZE,"%s%d",res->GetName().c_str(), res->GetNum());
            snprintf(rnum,BUFF_SIZE,"%d",res->GetNum());
          }

        snprintf(buffer,BUFF_SIZE,"%7d %-6s   %9.4f %9.4f %9.4f %-5s %3s  %-8s %9.4f",
                 atom->GetIdx(),label,
                 atom->GetX(),atom->GetY(),atom->GetZ(),
                 str1.c_str(),
                 rnum,rlabel,
                 atom->GetPartialCharge());
        ofs << buffer << endl;
      }

    ofs << "@<TRIPOS>BOND" << endl;
    OBBond *bond;
    vector<OBBond*>::iterator j;
    OBSmartsPattern pat;
    string s1, s2;
    for (bond = mol.BeginBond(j);bond;bond = mol.NextBond(j))
      {
        s1 = bond->GetBeginAtom()->GetType();
        s2 = bond->GetEndAtom()->GetType();
        if (bond->IsAromatic() || s1 == "O.co2" || s2 == "O.co2")
          strcpy(label,"ar");
        else if (bond->IsAmide())
          strcpy(label,"am");
        else
          snprintf(label,BUFF_SIZE,"%d",bond->GetBO());

        snprintf(buffer, BUFF_SIZE,"%6d %5d %5d   %2s",
                 bond->GetIdx()+1,bond->GetBeginAtomIdx(),bond->GetEndAtomIdx(),
                 label);
        ofs << buffer << endl;
      }
    // NO trailing blank line (PR#1868929).
    //    ofs << endl;

    return(true);
  }

  int MOL2Format::SkipObjects(int n, OBConversion* pConv)
  {
    const char txt[] = "@<TRIPOS>MOLECULE";
    istream& ifs = *pConv->GetInStream();
    if(!ifs)
      return -1;
    if(n>0 && ifs.peek()==txt[0])
      ifs.ignore(); // move past '@' so that next mol will be found
    do {
      ignore(ifs, txt);
    } while(ifs && (--n)>0);

    if(!ifs.eof())
      ifs.seekg(1-sizeof(txt), ios::cur);//1 for '/0'
    char ch = ifs.peek();
   return 1;
  }

} //namespace OpenBabel
