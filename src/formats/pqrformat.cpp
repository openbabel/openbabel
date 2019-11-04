/**********************************************************************
Copyright (C) 2008 Geoffrey R. Hutchison
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
#include <openbabel/elements.h>
#include <openbabel/generic.h>
#include <openbabel/obiter.h>
#include <openbabel/data.h>



#include <vector>
#include <map>

#include <sstream>
#include <cstdlib>

using namespace std;
namespace OpenBabel
{

  class PQRFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    PQRFormat()
    {
      OBConversion::RegisterFormat("pqr",this, "chemical/x-pqr");
    }

    virtual const char* Description() //required
    {
      return
        "PQR format\n"
        "Read Options e.g. -as\n"
        "  s  Output single bonds only\n"
        "  b  Disable bonding entirely\n\n";
    };

    virtual const char* SpecificationURL()
    { return "";};

    virtual const char* GetMIMEType()
    { return "chemical/x-pqr"; };

    //*** This section identical for most OBMol conversions ***
    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual int SkipObjects(int n, OBConversion* pConv);
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

  };
  //***

  //Make an instance of the format class
  PQRFormat thePQRFormat;

  ////////////////////////////////////////////////////////////////
  /// Utility functions
  static bool parseAtomRecord(char *buffer, OBMol &mol,int /*chainNum*/);
  static double parseAtomRadius(char *buffer, OBMol &mol);
  static double parseAtomCharge(char *buffer, OBMol &mol);

  /////////////////////////////////////////////////////////////////
  int PQRFormat::SkipObjects(int n, OBConversion* pConv)
  {
    if (n == 0)
      ++ n;
    istream &ifs = *pConv->GetInStream();
    char buffer[BUFF_SIZE];
    while (n && ifs.getline(buffer,BUFF_SIZE))
      {
        if (EQn(buffer,"ENDMDL",6))
          -- n;
      }

    return ifs.good() ? 1 : -1;
  }
  /////////////////////////////////////////////////////////////////
  bool PQRFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title = pConv->GetTitle();

    int chainNum = 1;
    char buffer[BUFF_SIZE];
    vector<double> charges, radii;
    string line, key, value;

    mol.SetTitle(title);
    mol.SetChainsPerceived(); // It's a PDB-like file, we read all chain/res info.

    mol.BeginModify();
    while (ifs.good() && ifs.getline(buffer,BUFF_SIZE))
      {
        if (EQn(buffer,"ENDMDL",6))
          break;
        if (EQn(buffer,"END",3)) {
          // eat anything until the next ENDMDL
          while (ifs.getline(buffer,BUFF_SIZE) && !EQn(buffer,"ENDMDL",6));
          break;
        }
        if (EQn(buffer,"TER",3)) {
          chainNum++;
          continue;
        }
        if (EQn(buffer,"ATOM",4) || EQn(buffer,"HETATM",6))
          {
            if( ! parseAtomRecord(buffer,mol,chainNum))
              {
                stringstream errorMsg;
                errorMsg << "WARNING: Problems reading a PQR file\n"
                         << "  Problems reading a ATOM/HETATM record.\n";
                obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obError);
              }

            // Read in the partial charge and radius too
            charges.push_back( parseAtomCharge(buffer, mol) );
            radii.push_back( parseAtomRadius(buffer, mol) );
            continue;
          }
        }

    if (!mol.NumAtoms()) { // skip the rest of this processing
      mol.EndModify();
      return(false);
    }

    // Use residue definitions to assign bond orders
    resdat.AssignBonds(mol);

    mol.EndModify();

    /*Now assign hetatm bonds based on distance*/
    if (!pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.ConnectTheDots();

    if (!pConv->IsOption("s",OBConversion::INOPTIONS) && !pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();

    FOR_ATOMS_OF_MOL(a, mol) {
      // WARNING: Atom index issue here
      a->SetPartialCharge(charges[a->GetIdx() - 1]);

      cerr << " charge : " << charges[a->GetIdx() - 1] << endl;

      if (!a->HasData("Radius")) {
        std::ostringstream s;
        s << radii[ a->GetIdx()-1 ];
        OBPairData *p = new OBPairData;
        p->SetAttribute("Radius");
        p->SetValue( s.str() );
        a->SetData(p);
      }

      cerr << " radius : " << radii[a->GetIdx() - 1] << endl;

    }
    mol.SetPartialChargesPerceived();
    mol.SetChainsPerceived();

    // clean out remaining blank lines
    std::streampos ipos;
    do
    {
      ipos = ifs.tellg();
      ifs.getline(buffer,BUFF_SIZE);
    }
    while(strlen(buffer) == 0 && !ifs.eof() );
    ifs.seekg(ipos);

    return(true);
  }

  static double parseAtomCharge(char *buffer, OBMol &mol)
  // In PQR format, either:
  // Field name, atom number, atom name, residue name, residue number
  //    x y z charge radius element
  // OR
  // Field, atom number, atom name, chain id, residue number, X, Y, Z, chg, rad, ele
  {
    vector<string> vs;
    tokenize(vs,buffer);

    OBAtom *atom = mol.GetAtom(mol.NumAtoms());

    if (vs.size() == 11)//add element, Zhixiong Zhao
      return atof(vs[8].c_str());
    else if (vs.size() == 12)
      return atof(vs[9].c_str());

    return 0.0;
  }

  static double parseAtomRadius(char *buffer, OBMol &mol)
  {
    vector<string> vs;
    tokenize(vs,buffer);

    OBAtom *atom = mol.GetAtom(mol.NumAtoms());

    if (vs.size() == 11)
      return atof(vs[9].c_str());
    else if (vs.size() == 12)
      return atof(vs[10].c_str());

    return 0.0;
  }

  static bool parseAtomRecord(char *buffer, OBMol &mol,int /*chainNum*/)
  /* ATOMFORMAT "(i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3,2f6.2,1x,i3)" */
  {
    string sbuf = &buffer[6];
    if (sbuf.size() < 48)
      return(false);

    bool hetatm = (EQn(buffer,"HETATM",6)) ? true : false;

    /* serial number */
    string serno = sbuf.substr(0,5);
    //SerialNum(the_atom) = atoi(tmp_str);

    /* atom name */
    string atmid = sbuf.substr(6,4);

    /* chain */
    char chain = sbuf.substr(15,1)[0];

    /* element */
    string element;
    if (sbuf.size() > 71)
      element = sbuf.substr(70,2);
    else
      element = "  ";

    //trim spaces on the right and left sides
    while (!atmid.empty() && atmid[0] == ' ')
      atmid = atmid.substr(1,atmid.size()-1);

    while (!atmid.empty() && atmid[atmid.size()-1] == ' ')
      atmid = atmid.substr(0,atmid.size()-1);

    /* residue name */

    string resname = sbuf.substr(11,3);
    if (resname == "   ")
      resname = "UNK";
    else
      {
        while (!resname.empty() && resname[0] == ' ')
          resname = resname.substr(1,resname.size()-1);

        while (!resname.empty() && resname[resname.size()-1] == ' ')
          resname = resname.substr(0,resname.size()-1);
      }

    string type;
    if (EQn(buffer,"ATOM",4))
      {
        type = atmid.substr(0,2);
        if (isdigit(type[0])) {
          // sometimes non-standard files have, e.g 11HH
          if (!isdigit(type[1])) type = atmid.substr(1,1);
          else type = atmid.substr(2,1);
        } else if ((sbuf[6] == ' ' &&
                    strncasecmp(type.c_str(), "Zn", 2) != 0 &&
                    strncasecmp(type.c_str(), "Fe", 2) != 0) ||
                   isdigit(type[1]))    //type[1] is digit in Platon
          type = atmid.substr(0,1);     // one-character element


        if (resname.substr(0,2) == "AS" || resname[0] == 'N')
          {
            if (atmid == "AD1")
              type = "O";
            if (atmid == "AD2")
              type = "N";
          }
        if (resname.substr(0,3) == "HIS" || resname[0] == 'H')
          {
            if (atmid == "AD1" || atmid == "AE2")
              type = "N";
            if (atmid == "AE1" || atmid == "AD2")
              type = "C";
          }
        if (resname.substr(0,2) == "GL" || resname[0] == 'Q')
          {
            if (atmid == "AE1")
              type = "O";
            if (atmid == "AE2")
              type = "N";
          }
        // fix: #2002557
        if (atmid[0] == 'H' &&
            (atmid[1] == 'D' || atmid[1] == 'E' ||
             atmid[1] == 'G' || atmid[1] == 'H')) // HD, HE, HG, HH, ..
          type = "H";
      }
    else //must be hetatm record
      {
        if (isalpha(element[1]) && (isalpha(element[0]) || (element[0] == ' ')))
          {
            if (isalpha(element[0]))
              type = element.substr(0,2);
            else
              type = element.substr(1,1);
            if (type.size() == 2)
              type[1] = tolower(type[1]);
          }
        else
          {
            if (isalpha(atmid[0]))
              {
              if (atmid.size() > 2 && (atmid[2] == '\0' || atmid[2] == ' '))
                type = atmid.substr(0,2);
              else if (atmid[0] == 'A') // alpha prefix
                type = atmid.substr(1, atmid.size() - 1);
              else
                type = atmid.substr(0,1);
              }
            else if (atmid[0] == ' ')
              type = atmid.substr(1,1); // one char element
            else
              type = atmid.substr(1,2);

            // Some cleanup steps
            if (atmid == resname)
              {
                type = atmid;
                if (type.size() == 2)
                  type[1] = tolower(type[1]);
              }
            else
              if (resname == "ADR" || resname == "COA" || resname == "FAD" ||
                  resname == "GPG" || resname == "NAD" || resname == "NAL" ||
                  resname == "NDP" || resname == "ABA")
                {
                  if (type.size() > 1)
                    type = type.substr(0,1);
                  //type.erase(1,type.size()-1);
                }
              else
                if (isdigit(type[0]))
                  {
                    type = type.substr(1,1);
                  }
                else
                  if (type.size() > 1 && isdigit(type[1]))
                    type = type.substr(0,1);
                  else
                    if (type.size() > 1 && isalpha(type[1])) {
                      if (type[0] == 'O' && type[1] == 'H')
                        type = type.substr(0,1); // no "Oh" element (e.g. 1MBN)
                      else if(isupper(type[1]))
                        {
                          type[1] = tolower(type[1]);
                        }
                    }
          }

      }

    OBAtom atom;
    /* X, Y, Z */
    string xstr = sbuf.substr(24,8);
    string ystr = sbuf.substr(32,8);
    string zstr = sbuf.substr(40,8);
    vector3 v(atof(xstr.c_str()),atof(ystr.c_str()),atof(zstr.c_str()));
    atom.SetVector(v);

    // useful for debugging unknown atom types (e.g., PR#1577238)
    //    cout << mol.NumAtoms() + 1 << " " << atmid << " type: " << type << endl;
    atom.SetAtomicNum(OBElements::GetAtomicNum(type.c_str()));

    /* residue sequence number */
    string resnum = sbuf.substr(16,4);
    OBResidue *res  = (mol.NumResidues() > 0) ? mol.GetResidue(mol.NumResidues()-1) : NULL;
    if (res == NULL || res->GetName() != resname
        || res->GetNumString() != resnum)
      {
        vector<OBResidue*>::iterator ri;
        for (res = mol.BeginResidue(ri) ; res ; res = mol.NextResidue(ri))
          if (res->GetName() == resname
              && res->GetNumString() == resnum
              && static_cast<int>(res->GetChain()) == chain)
            break;

        if (res == NULL)
          {
            res = mol.NewResidue();
            res->SetChain(chain);
            res->SetName(resname);
            res->SetNum(resnum);
          }
      }

    if (!mol.AddAtom(atom))
      return(false);
    else
      {
        OBAtom *atom = mol.GetAtom(mol.NumAtoms());

        res->AddAtom(atom);
        res->SetSerialNum(atom, atoi(serno.c_str()));
        res->SetAtomID(atom, sbuf.substr(6,4));
        res->SetHetAtom(atom, hetatm);

        return(true);
      }
  }

  bool PQRFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    unsigned int i;
    char buffer[BUFF_SIZE];
    char type_name[10], padded_name[10];
    char the_res[10];
    char the_chain = ' ';
    const char *element_name;
    int res_num;
    bool het=true;
    int model_num = 0;
    if (!pConv->IsLast() || pConv->GetOutputIndex() > 1)
      { // More than one molecule record
      model_num = pConv->GetOutputIndex(); // MODEL 1-based index
      snprintf(buffer, BUFF_SIZE, "MODEL %8d", model_num);
      ofs << buffer << endl;
      }

    if (strlen(mol.GetTitle()) > 0)
      snprintf(buffer, BUFF_SIZE, "COMPND    %s ",mol.GetTitle());
    else
      snprintf(buffer, BUFF_SIZE, "COMPND    UNNAMED");
    ofs << buffer << endl;

    snprintf(buffer, BUFF_SIZE, "AUTHOR    GENERATED BY OPEN BABEL %s",BABEL_VERSION);
    ofs << buffer << endl;

    // before we write any records, we should check to see if any coord < -1000
    // which will cause errors in the formatting

    double minX, minY, minZ;
    minX = minY = minZ = -999.0f;
    FOR_ATOMS_OF_MOL(a, mol)
      {
        if (a->GetX() < minX)
          minX = a->GetX();
        if (a->GetY() < minY)
          minY = a->GetY();
        if (a->GetZ() < minZ)
          minZ = a->GetZ();
      }

    vector3 transV = VZero;
    if (minX < -999.0)
      transV.SetX( -1.0*minX - 900.0 );
    if (minY < -999.0)
      transV.SetY( -1.0*minY - 900.0 );
    if (minZ < -999.0)
      transV.SetZ( -1.0*minZ - 900.0 );

    // if minX, minY, or minZ was never changed, shift will be 0.0f
    // otherwise, move enough so that smallest coord is > -999.0f
    mol.Translate(transV);

    OBAtom *atom;
    OBResidue *res;
    for (i = 1; i <= mol.NumAtoms(); i++)
      {
        atom = mol.GetAtom(i);
        strncpy(type_name, OBElements::GetSymbol(atom->GetAtomicNum()), sizeof(type_name));
        type_name[sizeof(type_name) - 1] = '\0';

        //two char. elements are on position 13 and 14 one char. start at 14
        if (strlen(type_name) > 1)
          type_name[1] = toupper(type_name[1]);
        else
          {
            char tmp[10];
            strncpy(tmp, type_name, 9); // make sure to null-terminate
            snprintf(type_name, sizeof(type_name), " %-3s", tmp);
          }

        if ( (res = atom->GetResidue()) != 0 )
          {
            het = res->IsHetAtom(atom);
            snprintf(the_res,4,"%s",(char*)res->GetName().c_str());
            snprintf(type_name,5,"%s",(char*)res->GetAtomID(atom).c_str());
            the_chain = res->GetChain();

            //two char. elements are on position 13 and 14 one char. start at 14
            if (strlen(OBElements::GetSymbol(atom->GetAtomicNum())) == 1)
              {
                if (strlen(type_name) < 4)
                  {
                    char tmp[16];
                    strncpy(tmp, type_name, 15); // make sure to null-terminate
                    snprintf(padded_name, sizeof(padded_name), " %-3s", tmp);
                    strncpy(type_name,padded_name,4);
                    type_name[4] = '\0';
                  }
                else
                  {
                    type_name[4] = '\0';
                  }
              }
            res_num = res->GetNum();
          }
        else
          {
            strcpy(the_res,"UNK");
            snprintf(padded_name,sizeof(padded_name), "%s",type_name);
            strncpy(type_name,padded_name,4);
            type_name[4] = '\0';
            res_num = 1;
          }

        element_name = OBElements::GetSymbol(atom->GetAtomicNum());
        //snprintf(buffer, BUFF_SIZE, "%s%5d %-4s %-3s %c%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s  \n",
        snprintf(buffer, BUFF_SIZE, "%s%5d %-4s %-3s %c%4d    %8.3f%8.3f%8.3f %11.8f%8.3f %2s  \n",
                 het?"HETATM":"ATOM  ",
                 i,
                 type_name,
                 the_res,
		 the_chain,
                 res_num,
                 atom->GetX(),
                 atom->GetY(),
                 atom->GetZ(),
                 atom->GetPartialCharge(),
                 atom->HasData("Radius")//use atom radius data,Zhixiong Zhao
				 	?atof(atom->GetData("Radius")->GetValue().c_str())
					:OBElements::GetVdwRad(atom->GetAtomicNum()),
                 element_name);
        ofs << buffer;
      }

    OBAtom *nbr;
    vector<OBBond*>::iterator k;
    for (i = 1; i <= mol.NumAtoms(); i ++)
      {
        atom = mol.GetAtom(i);
        if (atom->GetExplicitDegree() == 0)
          continue; // no need to write a CONECT record -- no bonds

        snprintf(buffer, BUFF_SIZE, "CONECT%5d", i);
        ofs << buffer;
        // Write out up to 4 real bonds per line PR#1711154
        int currentValence = 0;
        for (nbr = atom->BeginNbrAtom(k);nbr;nbr = atom->NextNbrAtom(k))
          {
            snprintf(buffer, BUFF_SIZE, "%5d", nbr->GetIdx());
            ofs << buffer;
            if (++currentValence % 4 == 0) {
              // Add the trailing space to finish this record
              ofs << "                                       \n";
              // write the start of a new CONECT record
              snprintf(buffer, BUFF_SIZE, "CONECT%5d", i);
              ofs << buffer;
            }
          }

        // Add trailing spaces
        int remainingValence = atom->GetExplicitDegree() % 4;
        for (int count = 0; count < (4 - remainingValence); count++) {
          snprintf(buffer, BUFF_SIZE, "     ");
          ofs << buffer;
        }
        ofs << "                                       \n";
      }

    snprintf(buffer, BUFF_SIZE, "MASTER        0    0    0    0    0    0    0    0 ");
    ofs << buffer;
    snprintf(buffer, BUFF_SIZE, "%4d    0 %4d    0\n",mol.NumAtoms(),mol.NumAtoms());
    ofs << buffer;

    ofs << "END\n";
    if (model_num) {
      ofs << "ENDMDL" << endl;
    }

    return(true);
  }


} //namespace OpenBabel
