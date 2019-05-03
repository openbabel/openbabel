/**********************************************************************
Copyright (C) 2002-2006 by Dr. Alex M. Clark and Geoffrey Hutchison
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
#include <openbabel/elements.h>
#include <cstdlib>

#include <sstream>

using namespace std;
namespace OpenBabel
{

  class CRK2DFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    CRK2DFormat()
    {
      OBConversion::RegisterFormat("crk2d", this, "chemical/x-crk2d");
    }

    virtual const char* Description() //required
    {
      return
        "Chemical Resource Kit diagram(2D)\n"
        "No comments yet\n";
    };

    virtual const char* SpecificationURL()
    {return "http://crk.sourceforge.net/";}; //optional

    virtual const char* GetMIMEType()
    { return "chemical/x-crk2d"; };

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
      return READONEONLY;
    };

    //*** This section identical for most OBMol conversions ***
    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

    static bool ReadCRK(std::istream &ifs,OBMol &mol,const char *classTag);
    static void WriteCRK(std::ostream &ofs,OBMol &mol,bool GroupCharges);

  };


  //Make an instance of the format class
  CRK2DFormat theCRK2DFormat;

  /////////////////////////////////////////////////////////////////
  bool CRK2DFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    mol.SetTitle( pConv->GetTitle()); //default title is the filename

    char buffer[BUFF_SIZE];//CM extra buffer

    if (!ifs.getline(buffer,BUFF_SIZE))
      {
        obErrorLog.ThrowError(__FUNCTION__, "File is empty!", obError);
        return(false);
      }
    if (!strstr(buffer,"<Property"))
      {
        obErrorLog.ThrowError(__FUNCTION__, "Not valid CRK XML", obWarning);
        return false;
      }
    if (!strstr(buffer,"\"DiagramStructure\""))
      {
        obErrorLog.ThrowError(__FUNCTION__,"Not CRK DiagramStructure (2D)", obWarning);
        return false;
      }

    mol.SetDimension(2);
    return ReadCRK(ifs,mol,"Structure2D");
  }

  ////////////////////////////////////////////////////////////////

  bool CRK2DFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    ofs << "<Property Type=\"DiagramStructure\">" <<  endl;
    ofs << " <Structure2D>" << endl;

    WriteCRK(ofs,mol,true);

    ofs << " </Structure2D>" << endl;
    ofs << "</Property>" << endl;

    return true;
  }

  //******************************************************
  class CRK3DFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    CRK3DFormat()
    {
      OBConversion::RegisterFormat("crk3d", this, "chemical/x-crk3d");
    }

    virtual const char* Description() //required
    {
      return
        "Chemical Resource Kit 3D format\n"
        "No comments yet\n";
    };

    virtual const char* SpecificationURL()
    {return "http://crk.sourceforge.net/";}; //optional

    virtual const char* GetMIMEType()
    { return "chemical/x-crk3d"; };

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
      return READONEONLY;
    };

    //*** This section identical for most OBMol conversions ***
    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
  };
  //***

  //Make an instance of the format class
  CRK3DFormat theCRK3DFormat;

  /////////////////////////////////////////////////////////////////
  bool CRK3DFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    mol.SetTitle( pConv->GetTitle()); //default title is the filename

    char buffer[BUFF_SIZE];//CM extra buffer

    if (!ifs.getline(buffer,BUFF_SIZE))
      {
        obErrorLog.ThrowError(__FUNCTION__, "File is empty!", obError);
        return(false);
      }
    if (!strstr(buffer,"<Property"))
      {
        obErrorLog.ThrowError(__FUNCTION__, "Not valid CRK XML", obWarning);
        return false;
      }
    if (!strstr(buffer,"\"ModelStructure\"") && !strstr(buffer,"\"XRayStructure\""))
      {
        obErrorLog.ThrowError(__FUNCTION__,"Not CRK ModelStructure or XRayStructure (3D).", obWarning);
        return false;
      }

    return CRK2DFormat::ReadCRK(ifs,mol,"Structure3D");
  }

  ////////////////////////////////////////////////////////////////

  bool CRK3DFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    ofs << "<Property Type=\"ModelStructure\">" <<  endl;
    ofs << " <Structure3D>" << endl;

    CRK2DFormat::WriteCRK(ofs,mol,true);

    ofs << " </Structure3D>" << endl;
    ofs << "</Property>" << endl;

    return true;
  }

  //**************************************************************
  bool CRK2DFormat::ReadCRK(std::istream &ifs,OBMol &mol,const char *classTag)
  {
    bool foundClass=false;

#define MAX_ATOMS 1000

    int numAtoms=0;
    int statomID[MAX_ATOMS];

#define MAX_BONDS 1000

    int numBonds=0;
    int stbondFrom[MAX_BONDS],stbondTo[MAX_BONDS],stbondStyle[MAX_BONDS];
    double stbondOrder[MAX_BONDS];

    bool inAtom=false,inBond=false;
    int atomID,atomNumber;
    double atomX,atomY,atomZ,atomCharge;
    int bondFrom,bondTo,bondStyle;
    double bondOrder = 0.0f;
    char buffer[BUFF_SIZE];//was global

    mol.BeginModify();

    while (ifs.getline(buffer,BUFF_SIZE))
      {
        if (strstr(buffer,classTag) && foundClass == false)
          foundClass=true;
        else if (strstr(buffer,classTag) && foundClass == true)
          break;
        else if (strstr(buffer,"<Atom"))
          {
            atomID=0;
            char *tag=strstr(buffer,"ID=\"");
            if (tag)
              atomID=atoi(tag+4);
            if (atomID>0)
              {
                inAtom=true;
                atomNumber=0;
                atomX= atomY= atomZ= atomCharge =0.0;
              }
            else
              continue; // atomID <= 0
          }
        else if (strstr(buffer,"<Bond"))
          {
            inBond=true;
            bondFrom=bondTo=bondStyle=0;
            bondOrder=0;
          }
        else if (strstr(buffer,"</Atom>"))
          {
            if (inAtom && numAtoms<MAX_ATOMS)
              {
                OBAtom atm;
                atm.Clear();

                statomID[numAtoms++]=atomID;

                atm.SetAtomicNum(atomNumber);
                atm.SetVector(atomX,atomY,atomZ);
                atm.SetFormalCharge((int)atomCharge);

                if (!mol.AddAtom(atm))
                  {
                    obErrorLog.ThrowError(__FUNCTION__, "Unable to add atom.", obWarning);
                    return false;
                  }
              }
            inAtom=false;
          }
        else if (strstr(buffer,"</Bond>"))
          {
            if (inBond && numBonds<MAX_BONDS)
              {
                stbondFrom[numBonds]=bondFrom;
                stbondTo[numBonds]=bondTo;
                stbondOrder[numBonds]=bondOrder;
                stbondStyle[numBonds]=bondStyle;
                numBonds++;
              }
            inBond=false;
          }
        else
          {
            char *tag;
            if (inAtom)
              {
                tag=strstr(buffer,"<X>");
                if (tag)
                  atomX=atof(tag+3);
                tag=strstr(buffer,"<Y>");
                if (tag)
                  atomY=atof(tag+3);
                tag=strstr(buffer,"<Z>");
                if (tag)
                  atomZ=atof(tag+3);
                tag=strstr(buffer,"<Element>");
                if (tag)
                  {
                    char element[3]="\0\0";
                    element[0]=tag[9];
                    if (tag[10]>='a' && tag[10]<='z')
                      element[1]=tag[10];
                    atomNumber=OBElements::GetAtomicNum(element);
                  }
                tag=strstr(buffer,"<Charge>");
                if (tag)
                  atomCharge=atof(tag+8);
              }
            if (inBond)
              {
                tag=strstr(buffer,"<From>");
                if (tag)
                  bondFrom=atoi(tag+6);
                tag=strstr(buffer,"<To>");
                if (tag)
                  bondTo=atoi(tag+4);
                tag=strstr(buffer,"<Order>");
                if (tag)
                  bondOrder=atof(tag+7);
                tag=strstr(buffer,"<Style>");
                if (tag)
                  bondStyle=atoi(tag+7);
              }
          }
      }

    for(int n=0;n<numBonds;n++)
      {
        int fromIdx=0,toIdx=0;
        for(int i=0;i<numAtoms;i++)
          {
            if (stbondFrom[n]==statomID[i])
              fromIdx=i+1;
            if (stbondTo[n]==statomID[i])
              toIdx=i+1;
          }

        if (fromIdx>0 && toIdx>0)
          {
            OBAtom *from=mol.GetAtom(fromIdx),*to=mol.GetAtom(toIdx);

            int order=1;
            if (stbondOrder[n]==2)
              order=2;
            else if (stbondOrder[n]==3)
              order=3;
            else if (stbondOrder[n]==1.5)
              order=5;

            OBBond bnd;
            bnd.Set(n+1,from,to,order,0);

            if (stbondStyle[n]==1)
              bnd.SetWedge();
            if (stbondStyle[n]==2)
              bnd.SetHash();
            if (stbondOrder[n]==1.5)
              bnd.SetAromatic();

            if (!mol.AddBond(bnd))
              {
                obErrorLog.ThrowError(__FUNCTION__, "Unable to add bond.", obWarning);
                return false;
              }
          }
        else
          {
            stringstream errorMsg;
            errorMsg << "Unassigned bond ID (" << stbondFrom[n]
                     << " " << stbondTo[n] << "), source may be invalid.";
            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
            return false;
          }
      }

    mol.EndModify();

    // we likely have an </Property> line to gobble up
    if (ifs.peek() != EOF && ifs.good())
      {
        ifs.getline(buffer,BUFF_SIZE);
        if (strstr(buffer,"</Property>") == 0)
          return false; // something messed up
      }

    return foundClass;
  }

  void CRK2DFormat::WriteCRK(std::ostream &ofs,OBMol &mol,bool GroupCharges)
  {
    double groupCharge=0;
    if (GroupCharges)
      for(unsigned int n=1;n<=mol.NumAtoms();n++)
        groupCharge+=mol.GetAtom(n)->GetFormalCharge();

    ofs << "  <Group Charge=\"" << groupCharge << "\" Spin=\"0\">" << endl;

    for(unsigned int n=1;n<=mol.NumAtoms();n++)
      {
        OBAtom *atm=mol.GetAtom(n);

        int id=atm->GetIdx(),atomnum=atm->GetAtomicNum();
        double x=atm->GetX(),y=atm->GetY(),z=atm->GetZ();
        const char *element=OBElements::GetSymbol(atomnum);
        double charge=0;
        if (!GroupCharges)
          charge=atm->GetFormalCharge();

        //ofs << ((int)atm) << endl;

        ofs << "   <Atom ID=\"" << id << "\">" << endl;
        ofs << "    <X>" << x << "</X>" << endl;
        ofs << "    <Y>" << y << "</Y>" << endl;
        ofs << "    <Z>" << z << "</Z>" << endl;
        ofs << "    <Element>" << element << "</Element>" << endl;
        if (charge!=0)
          ofs << "    <Charge>" << charge << "</Charge>" << endl;
        ofs << "   </Atom>" << endl;
      }

    for(unsigned int m=0;m<mol.NumBonds();m++)
      {
        OBBond *bnd=mol.GetBond(m);

        int from=bnd->GetBeginAtom()->GetIdx(),to=bnd->GetEndAtom()->GetIdx();
        double order=bnd->GetBondOrder();
        if (bnd->IsAromatic())
          order=1.5;
        int style=0;
        if (bnd->IsHash())
          style=1;
        if (bnd->IsWedge())
          style=2;

        ofs << "   <Bond>" << endl;
        ofs << "    <From>" << from << "</From>" << endl;
        ofs << "    <To>" << to << "</To>" << endl;
        ofs << "    <Order>" << order << "</Order>" << endl;
        ofs << "    <Style>" << style << "</Style>" << endl;
        ofs << "   </Bond>" << endl;
      }

    ofs << "  </Group>" << endl;
  }

} //namespace OpenBabel
