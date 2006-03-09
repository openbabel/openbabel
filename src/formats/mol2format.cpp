/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison
Some portions Copyright (C) 2004 by Chris Morley
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#include "babelconfig.h"

#include "mol.h"
#include "obconversion.h"
#include "obmolecformat.h"

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
    }

    virtual const char* Description() //required
    {
        return
            "Sybyl Mol2 format\n \
            No comments yet\n";
    };

  virtual const char* SpecificationURL()
  {
    return "http://www.tripos.com/data/support/mol2.pdf";
  }; //optional

  virtual const char* GetMIMEType() 
  { return "chemical/x-mol2"; };

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

    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
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

    int lcount;
    int natoms,nbonds;
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
        else if (lcount == 4) //energy
        {
            tokenize(vstr,buffer);
            if (!vstr.empty() && vstr[0] == "Energy" && vstr.size() == 3)
                mol.SetEnergy(atof(vstr[2].c_str()));
        }
        else if (lcount == 5) //comment
        {
            if ( buffer[0] )
            {
                len = (int) strlen(buffer)+1;
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
    bool hasPartialCharges=false;
    double x,y,z,pcharge;
    char temp_type[BUFF_SIZE], resname[BUFF_SIZE], atmid[BUFF_SIZE];
    int elemno, resnum = -1;

    ttab.SetFromType("SYB");
    for (i = 0;i < natoms;i++)
    {
        if (!ifs.getline(buffer,BUFF_SIZE))
            return(false);
        sscanf(buffer," %*s %s %lf %lf %lf %s %d %s %lf",
               atmid, &x,&y,&z, temp_type, &resnum, resname, &pcharge);

        atom.SetVector(x, y, z);

        str = temp_type;
        ttab.SetToType("ATN");
        ttab.Translate(str1,str);
        elemno = atoi(str1.c_str());

        // Handle "CL" and "BR" atom types!
        if( !elemno && isupper(temp_type[1]) )
        {
            temp_type[1] = (char)tolower(temp_type[1]);
            str = temp_type;
            ttab.Translate(str1,str);
            elemno = atoi(str1.c_str());
        }
        atom.SetAtomicNum(elemno);
        ttab.SetToType("INT");
        ttab.Translate(str1,str);
        atom.SetType(str1);
        atom.SetPartialCharge(pcharge);
        if (!mol.AddAtom(atom))
            return(false);
        if (!IsNearZero(pcharge))
            hasPartialCharges = true;

        // Add residue information if it exists
        if (resnum != -1 && resname != "")
        {
            OBResidue *res  = (mol.NumResidues() > 0) ?
                              mol.GetResidue(mol.NumResidues()-1) : NULL;
            if (res == NULL || res->GetName() != resname ||
		static_cast<int>(res->GetNum()) != resnum)
            {
                vector<OBResidue*>::iterator ri;
                for (res = mol.BeginResidue(ri) ; res ; res = mol.NextResidue(ri))
                    if (res->GetName() == resname &&
                            static_cast<int>(res->GetNum())
                            == resnum)
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

        sscanf(buffer,"%*d %d %d %s",&start,&end,temp_type);
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

 
    // update neighbour bonds information for each atom.
    vector<OBNodeBase*>::iterator apos;
    vector<OBEdgeBase*>::iterator bpos;
    OBAtom* patom;
    OBBond* pbond;
    
    for (patom = mol.BeginAtom(apos); patom; patom = mol.NextAtom(apos))
      {
	patom->ClearBond();
	for (pbond = mol.BeginBond(bpos); pbond; pbond = mol.NextBond(bpos))
	  {
	    if (patom == pbond->GetBeginAtom() || patom == pbond->GetEndAtom())
	      {
		patom->AddBond(pbond);
	      }
	  }
      }

    mol.EndModify();

    //must add generic data after end modify - otherwise it will be blown away
    if (comment)
    {
        OBCommentData *cd = new OBCommentData;
        cd->SetData(comment);
        mol.SetData(cd);
        delete [] comment;
        comment = NULL;
    }
    if (hasPartialCharges)
        mol.SetPartialChargesPerceived();

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

    sprintf(buffer, " %d %d 0 0 0", mol.NumAtoms(),mol.NumBonds());
    ofs << buffer << endl;
    ofs << "SMALL" << endl;

    //if (mol.HasPartialChargesPerceived()) ofs << "GASTEIGER" << endl;
    //  else                                  ofs << "NO_CHARGES" << endl;

    ofs << "GASTEIGER" << endl;
    ofs << "Energy = " << mol.GetEnergy() << endl;

    if (mol.HasData(OBGenericDataType::CommentData))
    {
        OBCommentData *cd = (OBCommentData*)mol.GetData(OBGenericDataType::CommentData);
        ofs << cd->GetData();
    }

    ofs << endl;
    ofs << "@<TRIPOS>ATOM" << endl;

    OBAtom *atom;
    OBResidue *res;

    vector<OBNodeBase*>::iterator i;
    vector<int> labelcount;
    labelcount.resize( etab.GetNumberOfElements() );

    ttab.SetFromType("INT");
    ttab.SetToType("SYB");

    for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
    {

        //
        //  Use sequentially numbered atom names if no residues
        //

        sprintf(label,"%s%d",
                etab.GetSymbol(atom->GetAtomicNum()),
                ++labelcount[atom->GetAtomicNum()]);
	strcpy(rlabel,"<1>");
	strcpy(rnum,"1");

        str = atom->GetType();
        ttab.Translate(str1,str);

        //
        //  Use original atom names if there are residues
        //

        if ( (res = atom->GetResidue()) )
        {
	  // use original atom names defined by residue
	  sprintf(label,"%s",(char*)res->GetAtomID(atom).c_str());
	  // make sure that residue name includes its number
	  sprintf(rlabel,"%s%d",res->GetName().c_str(), res->GetNum());
	  sprintf(rnum,"%d",res->GetNum());
        }

        sprintf(buffer,"%7d%1s%-6s%12.4f%10.4f%10.4f%1s%-5s%4s%1s %-8s%10.4f",
                atom->GetIdx(),"",label,
                atom->GetX(),atom->GetY(),atom->GetZ(),
                "",str1.c_str(),
                rnum,"",rlabel,
                atom->GetPartialCharge());
        ofs << buffer << endl;
    }

    ofs << "@<TRIPOS>BOND" << endl;
    OBBond *bond;
    vector<OBEdgeBase*>::iterator j;
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
	sprintf(label,"%d",bond->GetBO());
      sprintf(buffer, "%6d%6d%6d%3s%2s",
	      bond->GetIdx()+1,bond->GetBeginAtomIdx(),bond->GetEndAtomIdx(),
	      "",label);
      ofs << buffer << endl;
    }
    ofs << endl;

    return(true);
}

} //namespace OpenBabel
