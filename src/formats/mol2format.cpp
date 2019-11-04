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
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/elements.h>
#include <openbabel/generic.h>
#include <openbabel/kekulize.h>
#include <openbabel/obfunctions.h>
#include <openbabel/data.h>
#include <cstdlib>

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
      OBConversion::RegisterOptionParam("c", this, 0, OBConversion::INOPTIONS);
      OBConversion::RegisterOptionParam("c", this, 0, OBConversion::OUTOPTIONS);
      OBConversion::RegisterOptionParam("l", this, 0, OBConversion::OUTOPTIONS);
      OBConversion::RegisterOptionParam("u", this, 0, OBConversion::OUTOPTIONS);
    }

    virtual const char* Description() //required
    {
      return
        "Sybyl Mol2 format\n"
        "Read Options e.g. -ac\n"
        "  c               Read UCSF Dock scores saved in comments preceding molecules\n\n"
        "Write Options e.g. -xl\n"
        "  l               Output ignores residue information (only ligands)\n"
        "  c               Write UCSF Dock scores saved in comments preceding molecules\n"
        "  u               Do not write formal charge information in UNITY records\n\n";
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
  
  // Helper function for ReadMolecule
  // \return Is this atom a sulfur in a (di)thiocarboxyl (-CS2, -COS, CS2H or COSH) group?
  static bool IsThiocarboxylSulfur(OBAtom* queryatom)
  {
    if (queryatom->GetAtomicNum() != OBElements::Sulfur)
      return(false);
    if (queryatom->GetHvyDegree() != 1)
      return(false);

    OBAtom *atom = NULL;
    OBBond *bond;
    OBBondIterator i;

    for (bond = queryatom->BeginBond(i); bond; bond = queryatom->NextBond(i))
      if ((bond->GetNbrAtom(queryatom))->GetAtomicNum() == OBElements::Carbon)
      {
        atom = bond->GetNbrAtom(queryatom);
        break;
      }
    if (!atom)
      return(false);
    if (!(atom->CountFreeSulfurs() == 2)
      && !(atom->CountFreeOxygens() == 1 && atom->CountFreeSulfurs() == 1))
      return(false);

    //atom is connected to a carbon that has a total
    //of 2 attached free sulfurs or 1 free oxygen and 1 free sulfur
    return(true);
  }

  static bool IsOxygenOrSulfur(OBAtom *atom)
  {
    switch (atom->GetAtomicNum()) {
    case 8: case 16: return true;
    default: return false;
    }
  }

  static unsigned int GetAtomicNumAndIsotope(const char* symbol, int *isotope)
  {
    const char* p = symbol;
    switch (p[0]) {
    case 'D':
      if (p[1] == '\0') {
        *isotope = 2;
        return 1;
      }
      break;
    case 'T':
      if (p[1] == '\0') {
        *isotope = 3;
        return 1;
      }
      break;
    }
    return OBElements::GetAtomicNum(symbol);
  }

  //read from ifs until next rti is found and return it
  static string read_until_rti(istream & ifs) 
  {
      char buffer[BUFF_SIZE];
      for (;;)
      {
        if (!ifs.getline(buffer,BUFF_SIZE))
          return "";
        if (!strncmp(buffer,"@<TRIPOS>",9))
          return string(buffer);
      }      
  }
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

    // Prevent reperception
    mol.SetChainsPerceived();

    mol.BeginModify();

    for (;;)
      {
        if (!ifs.getline(buffer,BUFF_SIZE))
          return(false);
        if (pConv->IsOption("c", OBConversion::INOPTIONS)!=NULL && EQn(buffer,"###########",10))
          {
            char attr[32], val[32];
            sscanf(buffer, "########## %[^:]:%s", attr, val);
            OBPairData *dd = new OBPairData;
            dd->SetAttribute(attr);
            dd->SetValue(val);
            dd->SetOrigin(fileformatInput);
            mol.SetData(dd);
          }
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
    int isotope = 0;
    bool has_explicit_hydrogen = false;
    bool has_residue_information = false;

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
          // check if it's "Du" or "Xx" and the element is in the atom name
          if (str == "Du" || str == "Xx") {
            str = atmid;
            for (unsigned int i = 0; i < str.length(); ++i)
              if (!isalpha(str[i])) {
                str.erase(i);
                break; // we've erased the end of the string
              }
          }


          std::stringstream errorMsg;
          errorMsg << "This Mol2 file is non-standard. Problem with molecule: "
                   << mol.GetTitle()
                   << " Cannot interpret atom types correctly, instead attempting to interpret atom type: "
                   << str << " as elements instead.";
          obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);

          string::size_type dotPos = str.find('.');
          if (dotPos == string::npos) {
            elemno = GetAtomicNumAndIsotope(str.c_str(), &isotope);
          }
        }

        atom.SetAtomicNum(elemno);
        if (isotope)
          atom.SetIsotope(isotope);
        else if (elemno == 1)
          has_explicit_hydrogen = true;
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
            has_residue_information = true;
            OBResidue *res  = (mol.NumResidues() > 0) ?
              mol.GetResidue(mol.NumResidues()-1) : NULL;
            if (res == NULL || res->GetName() != resname ||
                res->GetNum() != resnum)
              {
                vector<OBResidue*>::iterator ri;
                for (res = mol.BeginResidue(ri) ; res ; res = mol.NextResidue(ri))
                  if (res->GetName() == resname &&
                      res->GetNum() == resnum)
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

    string nextrti;
    do { nextrti = read_until_rti(ifs); } 
    while(nextrti != "@<TRIPOS>UNITY_ATOM_ATTR" && nextrti != "@<TRIPOS>BOND" && nextrti.length() > 0);

    if(nextrti == "@<TRIPOS>UNITY_ATOM_ATTR")
    { //read in formal charge information, must be done before Kekulization
        int aid = 0, num = 0;
        while (ifs.peek() != '@' && ifs.getline(buffer,BUFF_SIZE))
        {
          sscanf(buffer,"%d %d",&aid, &num);
          for(int i = 0; i < num; i++) 
          {
            if (!ifs.getline(buffer,BUFF_SIZE))
              return(false);
            if(strncmp(buffer, "charge", 6) == 0)
            {
              int charge = 0;
              sscanf(buffer,"%*s %d",&charge);
              if(aid <= mol.NumAtoms()) 
              {
                OBAtom *atom = mol.GetAtom(aid);
                atom->SetFormalCharge(charge);
              }
            }
          }
        }
    }

    while(nextrti != "@<TRIPOS>BOND" && nextrti.length() > 0)
      nextrti = read_until_rti(ifs);

    if(nextrti != "@<TRIPOS>BOND")
      return false;

    int start, end;
    bool needs_kekulization = false;
    for (i = 0; i < nbonds; i++)
      {
        if (!ifs.getline(buffer,BUFF_SIZE))
          return(false);

        sscanf(buffer,"%*d %d %d %1024s",&start,&end,temp_type);
        str = temp_type;
        unsigned int flags = 0;
        int order;
        if (str == "ar" || str == "AR" || str == "Ar") {
          order = 1;
          flags = OB_AROMATIC_BOND;
          needs_kekulization = true;
        }
        else if (str == "AM" || str == "am" || str == "Am")
          order = 1;
        else
          order = atoi(str.c_str());

        mol.AddBond(start, end, order, flags);
      }

    // TODO: Add a test case for the statement below of Paolo Tosco
    //       - I am currently assuming that is not a problem for the
    //         the current kekulization code, but it needs to be
    //         checked
    //   "Make a pass to ensure that there are no double bonds
    //    between atoms which are also involved in aromatic bonds
    //    as that may ill-condition kekulization (fixes potential
    //    issues with molecules like CEWYIM30 (MMFF94 validation suite)"

    mol.SetAromaticPerceived(); // don't trigger reperception

    if (has_explicit_hydrogen) {
      FOR_ATOMS_OF_MOL(atom, mol) {
        unsigned int total_valence = atom->GetTotalDegree();
        switch (atom->GetAtomicNum()) {
        case 8:
          if (total_valence != 1) continue;
          if (strcmp(atom->GetType(), "O2") != 0) continue; // TODO: the O.co2 type is lost by this point
          {
            OBAtomBondIter bit(&*atom);
            if (!bit->IsAromatic() && bit->GetBondOrder() == 1)
              atom->SetFormalCharge(-1); // set -1 charge on dangling O.co2
          }
          break;
        case 17: // Cl
          if (total_valence == 0)
            atom->SetFormalCharge(-1);
          break;
        }
      }
    }

    // Kekulization is neccessary if an aromatic bond is present
    if (needs_kekulization) {
      // "de-aromatize" carboxylates and (di)thiocarboxylates
      // The typical case (in our test suite anyway) is a carboxylate binding to
      // a metal ion. The two O's have charges of -0.5 and the bond orders are aromatic
      FOR_ATOMS_OF_MOL(atom, mol) {
        OBAtom* oxygenOrSulfur = &*atom;
        // Look first for a terminal O/S
        if (!IsOxygenOrSulfur(oxygenOrSulfur) || oxygenOrSulfur->GetTotalDegree() != 1) continue;
        OBAtomBondIter bitA(oxygenOrSulfur);
        OBBond *bondA = &*bitA;
        if (!bondA->IsAromatic()) continue;
        // Look for the carbon
        OBAtom *carbon = bondA->GetNbrAtom(oxygenOrSulfur);
        if (carbon->GetAtomicNum() != 6) continue;
        // Look for the other oxygen or sulfur
        OBAtom* otherOxygenOrSulfur = (OBAtom*)0;
        OBBond* bondB = (OBBond*)0;
        FOR_BONDS_OF_ATOM(bitB, carbon) {
          if (&*bitB == bondA || !bitB->IsAromatic()) continue;
          OBAtom* nbr = bitB->GetNbrAtom(carbon);
          if (IsOxygenOrSulfur(nbr) && nbr->GetTotalDegree() == 1) {
            otherOxygenOrSulfur = nbr;
            bondB = &*bitB;
          }
        }
        if(otherOxygenOrSulfur->GetFormalCharge() != 0) continue; //formal charge already set on one
        if (!otherOxygenOrSulfur) continue;

        // Now set as C(=O)O
        bondA->SetAromatic(false);
        oxygenOrSulfur->SetFormalCharge(-1);

        bondB->SetAromatic(false);
        bondB->SetBondOrder(2);
      }

      // First of all, set the atoms at the ends of the aromatic bonds to also
      // be aromatic. This information is required for OBKekulize.
      FOR_BONDS_OF_MOL(bond, mol) {
        if (bond->IsAromatic()) {
          bond->GetBeginAtom()->SetAromatic();
          bond->GetEndAtom()->SetAromatic();
        }
      }
      bool ok = OBKekulize(&mol);
      if (!ok) {
        stringstream errorMsg;
        errorMsg << "Failed to kekulize aromatic bonds in MOL2 file";
        std::string title = mol.GetTitle();
        if (!title.empty())
          errorMsg << " (title is " << title << ")";
        errorMsg << endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
        // return false; Should we return false for a kekulization failure?
      }
    }

    mol.EndModify();

    // Suggestion by Liu Zhiguo 2008-01-26
    // Mol2 files define atom types -- there is no need to re-perceive
    mol.SetAtomTypesPerceived();

    if (has_residue_information)
      mol.SetChainsPerceived();

    if (!has_explicit_hydrogen) {
      // Guess how many hydrogens are present on each atom based on typical valencies
      // TODO: implement the MOL2 valence model (if it exists)
      FOR_ATOMS_OF_MOL(matom, mol) {
        if (matom->GetImplicitHCount() == 0)
          OBAtomAssignTypicalImplicitHydrogens(&*matom);
      }
    }

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
    bool skipFormalCharge = pConv->IsOption("u", OBConversion::OUTOPTIONS)!=NULL;

    //The old code follows....
    string str,str1;
    char buffer[BUFF_SIZE],label[BUFF_SIZE];
    char rnum[BUFF_SIZE],rlabel[BUFF_SIZE];

    //Check if UCSF Dock style coments are on
    if(pConv->IsOption("c", OBConversion::OUTOPTIONS)!=NULL) {
        vector<OBGenericData*>::iterator k;
        vector<OBGenericData*> vdata = mol.GetData();
        ofs << endl;
        for (k = vdata.begin();k != vdata.end();++k) {
            if ((*k)->GetDataType() == OBGenericDataType::PairData
            && (*k)->GetOrigin()!=local //internal OBPairData is not written
            && (*k)->GetAttribute()!="PartialCharges")
            {
                ofs << "##########\t" << (*k)->GetAttribute() << ":\t" << ((OBPairData*)(*k))->GetValue() << endl;
            }
        }
        ofs << endl;
    }

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
    std::map<int, int> labelcount;

    ttab.SetFromType("INT");
    ttab.SetToType("SYB");

    bool hasFormalCharges = false;
    for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
      {

        //
        //  Use sequentially numbered atom names if no residues
        //

        snprintf(label,BUFF_SIZE, "%s%d",
                 OBElements::GetSymbol(atom->GetAtomicNum()),
                 ++labelcount[atom->GetAtomicNum()]);
        strcpy(rlabel,"<1>");
        strcpy(rnum,"1");

        str = atom->GetType();
        ttab.Translate(str1,str);

        if (atom->GetFormalCharge() != 0) hasFormalCharges = true;
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

    //store formal charge info; put before bonds so we don't have
    //to read past the end of the molecule to realize it is there
    if(hasFormalCharges && !skipFormalCharge) {
      //dkoes - to enable roundtriping of charges
      ofs << "@<TRIPOS>UNITY_ATOM_ATTR\n";
      for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
      {
        int charge = atom->GetFormalCharge();
        if (charge != 0) 
        {
          ofs << atom->GetIdx() << " 1\n"; //one attribute
          ofs << "charge " << charge << "\n"; //namely charge
        }
      }
    }

    ofs << "@<TRIPOS>BOND" << endl;
    OBBond *bond;
    vector<OBBond*>::iterator j;
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
          snprintf(label,BUFF_SIZE,"%d",bond->GetBondOrder());

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
