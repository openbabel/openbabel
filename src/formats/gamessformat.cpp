/**********************************************************************
  Copyright (C) 2000 by OpenEye Scientific Software, Inc.
  Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
  Some portions Copyright (C) 2004 by Chris Morley
  Some portions Copyright (C) 2006 by Donald E. Curtis

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
namespace Gamess
{
} 

namespace OpenBabel
{
#define BOHR_TO_ANGSTROM 0.529177249
#define ANGSTROM_TO_BOHR 1.889725989
  class GAMESSOutputFormat : public OBMoleculeFormat
  {

    public:
      //Register this format type ID
      GAMESSOutputFormat()
      {
        OBConversion::RegisterFormat("gam",this);
        OBConversion::RegisterFormat("gamout",this);
      }

      virtual const char* Description() //required
      {
        return
          "GAMESS Output\n \
          Read Options e.g. -as\n\
          s  Output single bonds only\n\
          b  Disable bonding entirely\n\n";
      };

      virtual const char* SpecificationURL()
      { return "http://www.msg.ameslab.gov/GAMESS/doc.menu.html";}; //optional

      virtual const char* GetMIMEType() 
      { return "chemical/x-gamess-output"; };

      //Flags() can return be any the following combined by | or be omitted if none apply
      // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
      virtual unsigned int Flags()
      {
        return READONEONLY | NOTWRITABLE;
      };

      //*** This section identical for most OBMol conversions ***
      ////////////////////////////////////////////////////
      /// The "API" interface functions
      virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);

    private:
      //! \brief Parse GAMESS options section.
      //void ParseSection(char *tag, OBSetData *set, istream &ifs);

  };
  //***

  //Make an instance of the format class
  GAMESSOutputFormat theGAMESSOutputFormat;


  class GAMESSInputFormat : public OBMoleculeFormat
  {
    public:
      //Register this format type ID
      GAMESSInputFormat()
      {
        OBConversion::RegisterFormat("inp",this, "chemical/x-gamess-input");
        OBConversion::RegisterFormat("gamin",this);
        // Command-line keywords
        OBConversion::RegisterOptionParam("k", NULL, 1, OBConversion::OUTOPTIONS);
        // Command-line keyword file
        OBConversion::RegisterOptionParam("f", NULL, 1, OBConversion::OUTOPTIONS);
      }


      virtual const char* Description() //required
      {
        return
          "GAMESS Input\n \
          Write Options e.g. -xk\n\
          k  \"keywords\" Use the specified keywords for input\n\
          f    <file>     Read the file specified for input keywords\n\n";
      };

      virtual const char* SpecificationURL()
      {return "http://www.msg.ameslab.gov/GAMESS/doc.menu.html";}; //optional

      virtual const char* GetMIMEType() 
      { return "chemical/x-gamess-input"; };

      //Flags() can return be any the following combined by | or be omitted if none apply
      // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
      virtual unsigned int Flags()
      {
        return WRITEONEONLY; // | NOTREADABLE;
      };

      ////////////////////////////////////////////////////
      /// The "API" interface functions
      virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
      virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);

  };

  //Make an instance of the format class
  GAMESSInputFormat theGAMESSInputFormat;

  /////////////////////////////////////////////////////////////////
  /* this function is for parsing default options too.  it is decided that
   * we should only parse parameters that the user specified in the input
   * deck and not EVERY option which is defaulted to by GAMESS.
  void GAMESSOutputFormat::ParseSection(char *tag, OBSetData *set, istream &ifs)
  {
    char buffer[BUFF_SIZE];
    OBSetData *curset = (OBSetData *)set->GetData(tag);
    if(!curset)
    {
      curset = new OBSetData();
      curset->SetAttribute(tag);
      set->AddData(curset);
    }

    string attr, value;
    char *ptr;

    for( ; ; )
    {
      ifs.getline(buffer,BUFF_SIZE);
      ptr = buffer;

      // trim initial line whitespace
      while((*ptr == ' ' || *ptr == '\t') && *ptr != '\0') ptr++;
      // If this is it be done
      if(*ptr == '\0') break;

      // parse a line
      while(true)
      {
        attr.clear();
        value.clear();

        // Trim leading whitespace
        while((*ptr == ' ' || *ptr == '\t') && *ptr != '\0') ptr++;

        // Read the attribute name
        while(*ptr != ' ' && *ptr != '=' && *ptr != '\0') attr += toupper(*(ptr++));

        // If this is it, be done
        if(*ptr == '\0') break;

        // Read to next non-whitespace
        while((*ptr == ' ' || *ptr == '\t') && *ptr != '\0') ptr++;

        // Keywords are only one word.  So we must have extra data we don't want.
        // So in this case we just ignore it and go on like we're ready for the
        // next pair.
        if(*ptr != '=') continue;

        // Read to next non-whitespace
        while((*ptr == ' ' || *ptr == '\t') && *ptr != '\0') ptr++;

        while((*ptr == ' ' || *ptr == '\t' || *ptr == '=') && *ptr != '\0') ptr++;

        // Read the attribute value.
        while(*ptr != ' ' && *ptr != '\0') value += toupper(*(ptr++));


        if(attr == "IGAUSS") { attr = "NGAUSS"; }

        // cout << attr << "/" << value << endl;

        OBPairData *data = new OBPairData();
        data = new OBPairData();
        data->SetAttribute(attr);
        data->SetValue(value);

        curset->AddData(data);
      }
    }
  }
*/
  bool GAMESSOutputFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title = pConv->GetTitle();

    char buffer[BUFF_SIZE];
    string str,str1;
    double x,y,z;
    OBAtom *atom;
    vector<string> vs;
    bool hasPartialCharges = false;
    int HOMO = 0;
    vector<double> orbitals;

    // must build generic data while we parse then add at the end.
    OBSetData *gmsset = new OBSetData();
    gmsset->SetAttribute("gamess");

    mol.Clear();
    mol.BeginModify();
    while (ifs.getline(buffer,BUFF_SIZE))
    {
      if(strstr(buffer,"ATOMIC                      COORDINATES (BOHR)") != NULL)
      {
        // mol.EndModify();
        mol.Clear();
        mol.BeginModify();
        ifs.getline(buffer,BUFF_SIZE);	// column headings
        ifs.getline(buffer,BUFF_SIZE);
        tokenize(vs,buffer);
        while (vs.size() == 5)
        {
          atom = mol.NewAtom();
          atom->SetAtomicNum(atoi(vs[1].c_str())); // Parse the current one
          x = atof((char*)vs[2].c_str()) * BOHR_TO_ANGSTROM;
          y = atof((char*)vs[3].c_str()) * BOHR_TO_ANGSTROM;
          z = atof((char*)vs[4].c_str()) * BOHR_TO_ANGSTROM;
          atom->SetVector(x,y,z);
          vs[1].erase(vs[1].size() - 2, 2);

          if (!ifs.getline(buffer,BUFF_SIZE))
            break;
          tokenize(vs,buffer);
        }
      }
      else if(strstr(buffer,"COORDINATES OF ALL ATOMS ARE (ANGS)") != NULL)
      {
        // mol.EndModify();
        mol.Clear();
        mol.BeginModify();
        ifs.getline(buffer,BUFF_SIZE);	// column headings
        ifs.getline(buffer,BUFF_SIZE);	// ---------------
        ifs.getline(buffer,BUFF_SIZE);
        tokenize(vs,buffer);
        while (vs.size() == 5)
        {
          atom = mol.NewAtom();
          atom->SetAtomicNum(atoi(vs[1].c_str())); // Parse the current one
          x = atof((char*)vs[2].c_str());
          y = atof((char*)vs[3].c_str());
          z = atof((char*)vs[4].c_str());
          atom->SetVector(x,y,z);
          vs[1].erase(vs[1].size() - 2, 2);

          if (!ifs.getline(buffer,BUFF_SIZE))
            break;
          tokenize(vs,buffer);
        }
      }
      else if(strstr(buffer,"MOPAC CHARGES") != NULL)
      {
        hasPartialCharges = true;
        ifs.getline(buffer,BUFF_SIZE);	// ---------------
        ifs.getline(buffer,BUFF_SIZE);	// column headings
        ifs.getline(buffer,BUFF_SIZE);
        tokenize(vs,buffer);
        while (vs.size() == 4)
        {
          atom = mol.GetAtom(atoi(vs[0].c_str()));
          atom->SetPartialCharge(atof(vs[2].c_str()));

          if (!ifs.getline(buffer,BUFF_SIZE))
            break;
          tokenize(vs,buffer);
        }
      }
      else if (strstr(buffer,"NUMBER OF OCCUPIED ORBITALS") != NULL)
      {
        tokenize(vs, buffer);
        if (vs.size() == 7) // alpha
          HOMO = atoi(vs[6].c_str());
        else if (vs.size() == 8) //beta
          HOMO = atoi(vs[7].c_str());
      }
      else if (strstr(buffer,"EIGENVECTORS") != NULL ||
          strstr(buffer,"MOLECULAR ORBITALS") != NULL)
      {
        ifs.getline(buffer,BUFF_SIZE); // ------ line
        ifs.getline(buffer,BUFF_SIZE); // blank
        orbitals.clear();

        while (strstr(buffer,"END OF RHF CALCULATION") == NULL &&
            strstr(buffer,"-------") == NULL)
        {
          //loop
          ifs.getline(buffer,BUFF_SIZE); // orbitals!
          ifs.getline(buffer,BUFF_SIZE); // energies in hartree
          tokenize(vs, buffer);
          for (unsigned int i = 0; i < vs.size(); i++)
            orbitals.push_back(27.21 * atof(vs[i].c_str()));

          ifs.getline(buffer,BUFF_SIZE); // symmetries
          // orbital coefficients
          while (ifs.getline(buffer,BUFF_SIZE) && strlen(buffer)
              && strstr(buffer,"END") == NULL 
              && strstr(buffer, "---") == NULL)
          { }
          if (!ifs.good())
            break;
        }
      }
      else if(strstr(buffer, "INPUT CARD> $"))
      {
        string attr, value;
        char *ptr;

        for( ; ; )
        {
          ptr = buffer + 14;
          tokenize(vs, ptr);

          if(vs.size() > 2)
          {
            OBSetData *curset = (OBSetData *)gmsset->GetData(vs[0]);
            if(!curset)
            {
              curset = new OBSetData();
              curset->SetAttribute(vs[0]);
              gmsset->AddData(curset);
            }

            for(int i=1; vs[i] != "$END" && i < vs.size(); i++) {
              string::size_type loc = vs[i].find("=",0);
              if(loc != string::npos)
              {
                // cout << vs[i].substr(0,loc) << " !!!!! " << vs[i].substr(loc+1) << endl;
                OBPairData *data = new OBPairData();
                data = new OBPairData();
                data->SetAttribute(vs[i].substr(0,loc));
                data->SetValue(vs[i].substr(loc+1));
                curset->AddData(data);
              }
            }
          }

          break;

        }
      }
      /*
      else if(strstr(buffer, "$CONTRL OPTIONS"))
      {
        ParseSection("CONTRL", gmsset, ifs);
      }
      else if(strstr(buffer, "$SYSTEM OPTIONS"))
      {
        ParseSection("SYSTEM", gmsset, ifs);
      }
      else if(strstr(buffer, "BASIS OPTIONS"))
      {
        ParseSection("BASIS", gmsset, ifs);
      }
      else if(strstr(buffer, "GUESS OPTIONS"))
      {
        ParseSection("GUESS", gmsset, ifs);
      }
      */
    }
    //    cerr << title << " " << HOMO << " " << orbitals[HOMO - 1] << " " << orbitals[HOMO] << endl;

    // add our gamess set
    mol.SetData(gmsset);

    // if we have basis set data we should set our global pair data
    OBSetData *bset = (OBSetData *) gmsset->GetData("BASIS");
    if(bset)
    {
      OBPairData *pd = NULL;
      string value;

      pd = (OBPairData *) bset->GetData("RUNTYPE");
      if(pd)
      {
        OBPairData *nd = new OBPairData();
        nd->SetAttribute("runtype");
        nd->SetValue(pd->GetValue());
        pmol->SetData(nd);
      }

      pd = (OBPairData *) bset->GetData("DFTTYPE");
      if(pd)
      {
        OBPairData *nd = new OBPairData();
        nd->SetAttribute("dfttype");
        nd->SetValue(pd->GetValue());
        pmol->SetData(nd);
      }

      OBPairData *gbasis = (OBPairData *) gmsset->GetData("GBASIS");
      OBPairData *ngauss = (OBPairData *) gmsset->GetData("NGAUSS");

      value = "";
      if(gbasis)
      {
        value = gbasis->GetValue();

        if(ngauss)
        {
          if(gbasis->GetValue() == "STO")
          {
            value.clear();
            value += "sto-";
            value += ngauss->GetValue();
            value += "g";
          }
          else if(ngauss->GetValue() == "3" || ngauss->GetValue() == "6")
          {
            value.clear();
            value = ngauss->GetValue();
            value += "-";
            value += gbasis->GetValue().substr(1);
            value += "G";
          }
        }

        OBPairData *basis = new OBPairData();
        basis->SetAttribute("basis");
        basis->SetValue(value);
        pmol->SetData(basis);
      }
    }
    if (!pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.ConnectTheDots();
    if (!pConv->IsOption("s",OBConversion::INOPTIONS) && !pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();

    mol.EndModify();
    if (hasPartialCharges)
      mol.SetPartialChargesPerceived();
    mol.SetTitle(title);
    return(true);
  }

  ////////////////////////////////////////////////////////////////
  bool GAMESSInputFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    //const char* title = pConv->GetTitle();

    char buffer[BUFF_SIZE];
    string str,str1;
    double x,y,z;
    OBAtom *atom;
    vector<string> vs;
    bool hasPartialCharges = false;

    mol.BeginModify();
    while	(ifs.getline(buffer,BUFF_SIZE))
    {
      if(strstr(buffer,"$DATA") != NULL)
      {
        // mol.EndModify();
        mol.Clear();
        mol.BeginModify();
        ifs.getline(buffer,BUFF_SIZE);	// title
        tokenize(vs,buffer);
        mol.SetTitle(buffer);
        ifs.getline(buffer,BUFF_SIZE);  // C1
        ifs.getline(buffer,BUFF_SIZE);
        while (strstr(buffer, "$END") == NULL)
        {
          tokenize(vs,buffer);
          if(vs.size() == 5) 
          {
            atom = mol.NewAtom();
            atom->SetAtomicNum(atoi(vs[1].c_str())); // Parse the current one
            x = atof((char*)vs[2].c_str());
            y = atof((char*)vs[3].c_str());
            z = atof((char*)vs[4].c_str());
            atom->SetVector(x,y,z);
            vs[1].erase(vs[1].size() - 2, 2);
          }

          if (!ifs.getline(buffer,BUFF_SIZE))
            break;
        }
      }
    }

    if (!pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.ConnectTheDots();
    if (!pConv->IsOption("s",OBConversion::INOPTIONS) && !pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();

    mol.EndModify();
    if (hasPartialCharges)
      mol.SetPartialChargesPerceived();
    //mol.SetTitle(title);
    return(true);
  }


  bool GAMESSInputFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    unsigned int i;
    char buffer[BUFF_SIZE];

    const char *keywords = pConv->IsOption("k",OBConversion::OUTOPTIONS);
    const char *keywordFile = pConv->IsOption("f",OBConversion::OUTOPTIONS);

    if (!keywords && !keywordFile)
      ofs << " $CONTRL COORD=CART UNITS=ANGS $END" << endl;
    if (keywords)
      ofs << pConv->IsOption("k", OBConversion::OUTOPTIONS) << endl;
    if (keywordFile)
    {
      ifstream kfstream(keywordFile);
      string keyBuffer;
      if (kfstream)
      {
        while (getline(kfstream, keyBuffer))
          ofs << keyBuffer << endl;
      }
    }

    ofs << " $DATA" << endl;
    ofs << mol.GetTitle() << endl;
    if (!mol.HasData(OBGenericDataType::SymmetryData))
      ofs << "C1" << endl;
    else
    {
      // \todo needs to be updated for point group symmetry recognition
      //   particularly for output of the symmetry elements
      //   and any necessary rotation for frame of reference for GAMESS
      ofs << "Put symmetry info here" << endl << endl;
    }

    OBAtom *atom;
    FOR_ATOMS_OF_MOL(atom, mol)
    {
      snprintf(buffer, BUFF_SIZE, "%-3s %4d.0    %8.5f  %8.5f  %8.5f ",
          etab.GetSymbol(atom->GetAtomicNum()),
          atom->GetAtomicNum(),
          atom->GetX(),
          atom->GetY(),
          atom->GetZ());
      ofs << buffer << endl;
    }

    ofs << " $END" << endl << endl << endl;
    return(true);
  }

} //namespace OpenBabel
