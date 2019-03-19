/**********************************************************************

Copyright (C) 2011-2019 by Culgi B.V., Leiden, The Netherlands
Authors: Shyamal K Nath, Ruben Serral Gracia and Paul Becherer for COF
Some portions copyright (C) 2004 by Chris Morley for template

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
#include <openbabel/obiter.h>
#include <openbabel/builder.h>
#include <openbabel/kekulize.h>
#include <openbabel/generic.h>
#include <iomanip>
#include <cstdlib>

using namespace std;
namespace OpenBabel
{
  class COFFormat : public OBMoleculeFormat
  {
  public:
    COFFormat()
    {
      OBConversion::RegisterFormat("cof", this);
    }

    virtual const char* Description()
    {
      return "Culgi object file format\n"
      "Culgi format\n"
      "No options currently \n";
    };

    /* Flags() can return be any of the following combined by | or be omitted if none apply
     NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY  DEFAULTFORMAT
     READBINARY  WRITEBINARY  READXML  ZEROATOMSOK*/
    virtual unsigned int Flags()
    {
      return READONEONLY | WRITEONEONLY;
    };

    /// Declarations for the "API" interface functions. Definitions are below
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

  };
  ////////////////////////////////////////////////////

  //Make an instance of the format class
  COFFormat theCOFFormat;

  /////////////////////////////////////////////////////////////////

  bool COFFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    stringstream errorMsg;
    char buffer[BUFF_SIZE];

    istream& ifs = *pConv->GetInStream();

    // Extract Culgi version, always on the first line
    vector<string> vs;
    string name;
    ifs.getline(buffer,BUFF_SIZE);
    const char *delimstr="\t\n\r"; // Only tabs can act as delimiters in a line, spaces cannot.
    const char *delimstrversion=" \t\n\r."; // Except when extracting the version number, as that line has a special format
    tokenize(vs, buffer, delimstrversion);
    if (vs.size() != 4)
    {
      obErrorLog.ThrowError(__FUNCTION__, "Unable to read Culgi Object File. First line incorrectly formatted or file is empty.", obWarning);
      return(false);
    }
    if (vs[0] != "@CulgiVersion:")
    {
      obErrorLog.ThrowError(__FUNCTION__, "Unable to read Culgi Object File. First line incorrectly formatted.", obWarning);
      return(false);
    }
    int version = atoi(vs[1].c_str());

    //-Extract the molecule name
    bool bFoundKey = false;
    while(ifs.getline(buffer,BUFF_SIZE))
    {
      tokenize(vs, buffer, delimstr);
      if(vs.size() == 3 && vs[0] == "molecule")
      {
        name = vs[1];
        bFoundKey = true;
        break;
      }
    }
    if (!bFoundKey)
    {
      obErrorLog.ThrowError(__FUNCTION__, "Unable to read Culgi Object File. File does not contain a molecule", obWarning);
      return(false);
    }
    pmol->SetTitle(name);

    //-Extract atoms and bonds.  Ignore other entities like atom groups or connectors
    //-defined in a molecule in culgi object file format.
    pmol->BeginModify();
    pmol->SetAutomaticFormalCharge(true);
    pmol->SetAutomaticPartialCharge(false);
    vector<string> vecAtomNames;
    unsigned int natoms = 0;
    int offset = 0;
    if(version > 7)
      offset = 1;
    while(ifs.getline(buffer,BUFF_SIZE))
    {
      tokenize(vs, buffer, delimstr);
      if(vs.size() == 0)
        continue;
      if(vs[0] == "atom" )
      {
        if(vs.size() < offset+9)
        {
          obErrorLog.ThrowError(__FUNCTION__, "Unable to read Culgi Object File. Atom line appears truncated.", obWarning);
          return(false);
        }
        OBAtom* atom = pmol->NewAtom();
        vecAtomNames.push_back(vs[1+offset]);

        //-set atomic number, or '0' if the atom type is not recognized
        std::string elem = vs[2+offset];
        if (elem == "X")
          elem = "Xx"; // Culgi 'X' is the same as Openbabel 'Xx'.
        int atomicNum = OBElements::GetAtomicNum(elem.c_str());
        atom->SetAtomicNum(atomicNum);
        if(atomicNum == 0)
          atom->SetType(elem);

        //-set charge
        char *endptr;
        double charge = strtod((char*)vs[4+offset].c_str(),&endptr);
        atom->SetPartialCharge(charge);

        //-set coordinates
        double x = strtod((char*)vs[6+offset].c_str(),&endptr);
        double y = strtod((char*)vs[7+offset].c_str(),&endptr);
        double z = strtod((char*)vs[8+offset].c_str(),&endptr);
        atom->SetVector(x,y,z);

        natoms++;
      }
      else if(vs[0] == "bond" )
      {
        if (vs.size() < 4)
        {
          obErrorLog.ThrowError(__FUNCTION__, "Unable to read Culgi Object File. Bond line appears truncated.", obWarning);
          return(false);
        }
        int iniatom = -1;
        int finatom = -1;

        if(version > 7)
        {
          iniatom = atoi(vs[1].c_str());
          finatom = atoi(vs[2].c_str());
        }
        else
        {
          for(int j=0 ; j< natoms ; j++)
          {
            if(vs[1] ==vecAtomNames[j])
              iniatom = j;
            else if(vs[2] == vecAtomNames[j])
              finatom = j;
          }

          if (iniatom < 0 || finatom < 0)
          {
            errorMsg << "Unable to read Culgi Object File. Atoms " << vs[0] << " or "<< vs[1] << "not found in atom list.\n";
            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
            return false;
          }
        }

        char *endptr;
        double order = strtod((char*)vs[3].c_str(),&endptr);
        int iorder = -1;
        unsigned int flags = 0;
        if (fabs(order - 1.0) < 1.e-6)
          iorder = 1;
        else if (fabs(order - 1.5) < 1.e-6)
        {
          // Aromatic bonds get order 1 here. OBKekulize, which is
          // called later, may then leave this 1 or make it 2.
          iorder = 1;
          flags = OB_AROMATIC_BOND;
        }
        else if (fabs(order - 2.0) < 1.e-6)
          iorder = 2;
        else if (fabs(order - 3.0) < 1.e-6)
          iorder = 3;
        else
        {
          errorMsg << "Unable to read Culgi Object File. Bond order " << vs[3] << "not possible.\n";
          obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
          return false;
        }
        pmol->AddBond(iniatom+1, finatom+1, iorder, flags);
      }
      else if(vs[0] == "formalq")
      {
        if(vs.size()!=3)
        {
          obErrorLog.ThrowError(__FUNCTION__,
              "Unable to read Culgi Object File. Formal charge line appears truncated", obWarning);
          return(false);
        }
        int iAtom = atoi(vs[1].c_str());
        int iFormalCharge = atoi(vs[2].c_str());
        pmol->GetAtom(iAtom+1)->SetFormalCharge(iFormalCharge);
      }
    }
    pmol->SetAromaticPerceived();

    //-update neighbour bonds information for each atom.
    vector<OBAtom*>::iterator apos;
    vector<OBBond*>::iterator bpos;
    OBAtom* patom;
    OBBond* pbond;

    for (patom = pmol->BeginAtom(apos); patom; patom = pmol->NextAtom(apos))
    {
      patom->ClearBond();
      for (pbond = pmol->BeginBond(bpos); pbond; pbond = pmol->NextBond(bpos))
      {
        if (patom == pbond->GetBeginAtom() || patom == pbond->GetEndAtom())
          patom->AddBond(pbond);
      }
    }

    // Set atoms on both side of aromatic bond as aromatic, as needed by OBKekulize()
    for (pbond = pmol->BeginBond(bpos); pbond; pbond = pmol->NextBond(bpos))
    {
      if (pbond->IsAromatic())
      {
      pbond->GetBeginAtom()->SetAromatic();
      pbond->GetEndAtom()->SetAromatic();
      }
    }

    bool ok = OBKekulize(pmol);
    if(!ok)
    {
      obErrorLog.ThrowError(__FUNCTION__,
          "Failed to kekulize aromatic bonds in COF file", obWarning);
    }

    OBPairData *dd = new OBPairData;
    dd->SetAttribute("PartialCharges");
    dd->SetValue("USER_CHARGES");
    dd->SetOrigin(external);
    pmol->SetData(dd);
    pmol->SetPartialChargesPerceived();

    pmol->EndModify();
    return true;
  }

  ////////////////////////////////////////////////////////////////

  bool COFFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;
    stringstream errorMsg;

    ostream& ofs = *(pConv->GetOutStream());
    OBMol &mol = *pmol;

    ofs << "@CulgiVersion: 10.0.0\n"
      << "# Culgi Object File\n"
      << "# Generated by Open Babel\n\n"
      << "moleculekeys\n"
      << "\tname\n"
      << "\tid\n"
      << "end_moleculekeys\n\n"
      << "atomkeys\n"
      << "\tindex\n"
      << "\tname\n"
      << "\telement_type\n"
      << "\tforce_field_type\n"
      << "\tcharge\n"
      << "\tfixed\n"
      << "\tx\n"
      << "\ty\n"
      << "\tz\n"
      << "\tvelocity_x\n"
      << "\tvelocity_y\n"
      << "\tvelocity_z\n"
      << "\tvelocity_const_x\n"
      << "\tvelocity_const_y\n"
      << "\tvelocity_const_z\n"
      << "end_atomkeys\n\n"
      << "bondkeys\n"
      << "\tatom1_index\n"
      << "\tatom2_index\n"
      << "\tbond_order\n"
      << "end_bondkeys\n\n";

    FOR_ATOMS_OF_MOL(atom, mol)
    {
      if(atom->GetFormalCharge()!=0)
      {
        ofs << "formalqkeys\n"
          << "\tatom_index\n"
          << "\tformal_charge\n"
          << "end_formalqkeys\n\n";
        break;
      }
    }

    // Get molecule file name.  Some formats add file path
    // etc to the name, so we try to clean that up here
    const char *tit = pmol->GetTitle();
    std::string molname(tit);
    if(molname.empty())
      molname = pConv->GetTitle();
    size_t nPos = molname.find_last_of(".");
    if(nPos != std::string::npos)
      molname = molname.substr(0, nPos);
    nPos = molname.find_last_of("\\");
    if(nPos != std::string::npos)
      molname = molname.substr(nPos+1, molname.size());
    nPos = molname.find_last_of("/");
    if(nPos != std::string::npos)
      molname = molname.substr(nPos+1, molname.size());
    if(molname.empty())
      molname = "mol";
    ofs  << "molecule\t" << molname << "\t0" << endl;

    int i = 0;
    vector<string> names;
    vector<string> elems;
    vector<int> nelem;
    ostringstream sstream;
    FOR_ATOMS_OF_MOL(atom, mol)
    {
      sstream.str("");
      i++;
      string elem = OBElements::GetSymbol(atom->GetAtomicNum());
      // Culgi does not recognize atom type 'Xx' but does know 'X'.
      if(elem == "Xx")
        elem = "X";
      bool found = false;
      string ename;
      for( int j=0; j<elems.size(); j++)
      {
        if(elem == elems[j])
        {
          found = true;
          nelem[j]++;
          sstream.str("");
          sstream << elem << nelem[j];
          ename = sstream.str();
          names.push_back(ename);
          sstream.str("");
          break;
        }
      }
      if (!found)
      {
        elems.push_back(elem);
        nelem.push_back(1);
        sstream << elem << "1";
        ename = sstream.str();
        names.push_back(ename);
        sstream.str("");
      }

      ofs << "atom" << "\t"
        << i - 1 << "\t"
        << ename << "\t"
        << elem << "\t"
        << "X" << "\t"
        << atom->GetPartialCharge() << "\t"
        << "0" << "\t"
        << atom->GetX() << "\t"
        << atom->GetY() << "\t"
        << atom->GetZ() << "\t"
        << "0.00000" << "\t"
        << "0.00000" << "\t"
        << "0.00000" << "\t"
        << "0" << "\t"
        << "0" << "\t"
        << "0" << endl;
    }

    int nbonds = mol.NumBonds();
    OBBond *bond;
    vector<OBBond*>::iterator j;
    int i1, i2;
    for (bond = mol.BeginBond(j);bond;bond = mol.NextBond(j))
    {
      i1 = bond->GetBeginAtomIdx();
      i2 = bond->GetEndAtomIdx();
      string outlabel;
      sstream.str("");
      sstream << "bond\t" << i1-1 << "\t" << i2-1 << "\t";
      if (bond->IsAromatic())
        sstream <<"1.5";
      else
        sstream << bond->GetBondOrder() << ".0";
      outlabel = sstream.str();
      ofs << outlabel << endl;
    }

    int iAt = 0;
    FOR_ATOMS_OF_MOL(atom, mol)
    {
      int iFormalCharge = atom->GetFormalCharge();
      if(iFormalCharge!=0)
        ofs << "formalq\t" << iAt << "\t" << iFormalCharge << endl;
      iAt++;
    }

    ofs << "end_molecule" << endl;
    return true;
  }

} //namespace OpenBabel

