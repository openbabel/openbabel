/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
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

  class MOPACFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    MOPACFormat()
    {
      OBConversion::RegisterFormat("mopout",this, "chemical/x-mopac-out");
      OBConversion::RegisterFormat("moo",this, "chemical/x-mopac-out");
    }

    virtual const char* Description() //required
    {
      return
        "MOPAC Output format\n"
        "Read Options e.g. -as\n"
        "  s  Output single bonds only\n"
        "  b  Disable bonding entirely\n\n";
    };

    virtual unsigned int Flags()
    {
      return NOTWRITABLE;
    };

    virtual const char* GetMIMEType() 
    { return "chemical/x-mopac-out"; };

    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    //	virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv); Is Read Only
  };

  //Make an instance of the format class
  MOPACFormat theMOPACFormat;

  /////////////////////////////////////////////////////////////////
  bool MOPACFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title= pConv->GetTitle();

    char buffer[BUFF_SIZE];
    string str,str1;
    double x,y,z;
    OBAtom *atom;
    vector<string> vs;
    vector<double> charges;
    bool hasPartialCharges = false;
    double energy;
    OBVectorData *dipoleMoment = NULL;
    bool readingVibrations = false;
    vector< vector<vector3> > displacements; // vibrational displacements
    vector<double> frequencies, intensities;

    // Translation vectors (if present)
    vector3 translationVectors[3];
    int numTranslationVectors = 0;

    mol.BeginModify();
    while	(ifs.getline(buffer,BUFF_SIZE))
      {
        if(strstr(buffer,"  CARTESIAN COORDINATES") != NULL)
          {
            // mol.EndModify();
            mol.Clear();
            mol.BeginModify();
            ifs.getline(buffer,BUFF_SIZE);	// blank
            ifs.getline(buffer,BUFF_SIZE);	// column headings
            ifs.getline(buffer,BUFF_SIZE);	// blank
            ifs.getline(buffer,BUFF_SIZE);
            tokenize(vs,buffer);
            while (vs.size() == 5)
              {
                atom = mol.NewAtom();
                atom->SetAtomicNum(etab.GetAtomicNum(vs[1].c_str()));
                x = atof((char*)vs[2].c_str());
                y = atof((char*)vs[3].c_str());
                z = atof((char*)vs[4].c_str());
                atom->SetVector(x,y,z);

                if (!ifs.getline(buffer,BUFF_SIZE))
                  break;
                tokenize(vs,buffer);
              }
          }
        else if(strstr(buffer,"UNIT CELL TRANSLATION") != NULL)
          {
            numTranslationVectors = 0; // ignore old translationVectors
            ifs.getline(buffer,BUFF_SIZE);	// blank
            ifs.getline(buffer,BUFF_SIZE);	// column headings
            ifs.getline(buffer,BUFF_SIZE);
            tokenize(vs,buffer);
            while (vs.size() == 5)
              {
                x = atof((char*)vs[2].c_str());
                y = atof((char*)vs[3].c_str());
                z = atof((char*)vs[4].c_str());

                translationVectors[numTranslationVectors++].Set(x, y, z);
                if (!ifs.getline(buffer,BUFF_SIZE))
                  break;
                tokenize(vs,buffer);
              }
          }
        else if(strstr(buffer,"FINAL HEAT") != NULL)
          {
            sscanf(buffer,"%*s%*s%*s%*s%*s%lf",&energy);
            mol.SetEnergy(energy);
          }
        else if(strstr(buffer,"NET ATOMIC CHARGES") != NULL)
          {
            hasPartialCharges = true;
            charges.clear();
            ifs.getline(buffer,BUFF_SIZE);	// blank
            ifs.getline(buffer,BUFF_SIZE);	// column headings
            ifs.getline(buffer,BUFF_SIZE);
            tokenize(vs,buffer);
            if (vs.size() < 1) return false; // timvdm 18/06/2008
            while (strstr(vs[0].c_str(),"DIPOLE") == NULL)
              {
                if (vs.size() < 3) return false; // timvdm 18/06/2008
                atom = mol.GetAtom(atoi(vs[0].c_str()));
                atom->SetPartialCharge(atof(vs[2].c_str()));
                charges.push_back(atof(vs[2].c_str()));

                if (!ifs.getline(buffer,BUFF_SIZE))
                  break;
                tokenize(vs,buffer);
                if (vs.size() < 1) vs.push_back(string()); // timvdm 18/06/2008
              }
            // Now we should be at DIPOLE line
            ifs.getline(buffer,BUFF_SIZE);	// POINT CHARGE
            ifs.getline(buffer,BUFF_SIZE);	// HYBRID
            ifs.getline(buffer,BUFF_SIZE);	// SUM
            tokenize(vs, buffer);
            if (vs.size() == 5) {
              if (dipoleMoment)
                delete dipoleMoment;

              dipoleMoment = new OBVectorData;
              double x, y, z;
              x = atof(vs[1].c_str());
              y = atof(vs[2].c_str());
              z = atof(vs[3].c_str());
              dipoleMoment->SetData(x, y, z);
              dipoleMoment->SetAttribute("Dipole Moment");
              dipoleMoment->SetOrigin(fileformatInput);
            }

            if (!ifs.getline(buffer,BUFF_SIZE))
              break;
          }
        else if(strstr(buffer,"MASS-WEIGHTED COORDINATE ANALYSIS") != NULL)
          { // the correct vibrations -- earlier bits aren't mass-weighted
            readingVibrations = true;
            if (!ifs.getline(buffer,BUFF_SIZE))
              break;
          }
        else if (readingVibrations && strstr(buffer, "Root No.") != NULL)
          {
            ifs.getline(buffer, BUFF_SIZE); // blank line
            ifs.getline(buffer, BUFF_SIZE); // symmetry labels (for OB-2.3)
            ifs.getline(buffer, BUFF_SIZE); // blank
            ifs.getline(buffer, BUFF_SIZE); // frequencies
            tokenize(vs, buffer);
            for (unsigned int i = 0; i < vs.size(); ++i) {
              frequencies.push_back(atof(vs[i].c_str()));
            }
            ifs.getline(buffer, BUFF_SIZE); // blank

            // now real work
            int prevModeCount = displacements.size();
            int newModes = frequencies.size() - displacements.size();
            vector<vector3> displacement;
            for (unsigned int i = 0; i < newModes; ++i) {
              displacements.push_back(displacement);
            } 

            ifs.getline(buffer, BUFF_SIZE);
            tokenize(vs, buffer);
            int modeCount = vs.size();
            vector<double> x, y, z;
            while(modeCount > 1) {
              x.clear();
              for (unsigned int i = 1; i < modeCount; ++i) {
                x.push_back(atof(vs[i].c_str()));
              }
              y.clear();
              ifs.getline(buffer, BUFF_SIZE);
              tokenize(vs, buffer);
              for (unsigned int i = 1; i < modeCount; ++i) {
                y.push_back(atof(vs[i].c_str()));
              }
              
              z.clear();
              ifs.getline(buffer, BUFF_SIZE);
              tokenize(vs, buffer);
              for (unsigned int i = 1; i < modeCount; ++i) {
                z.push_back(atof(vs[i].c_str()));
              }

              // OK, now we have x, y, z for all new modes for one atom
              for (unsigned int i = 0; i < modeCount - 1;  ++i) {
                displacements[prevModeCount + i].push_back(vector3(x[i], y[i], z[i]));
              }

              // Next set of atoms
              ifs.getline(buffer, BUFF_SIZE);
              tokenize(vs, buffer);
              modeCount = vs.size();
            }
          }
        else if (readingVibrations && strstr(buffer, "T-DIPOLE") != NULL)
          {
            unsigned int currentIntensity = intensities.size();
            tokenize(vs, buffer);
            if (vs.size() < 2)
              break;

            double transDipole = atof(vs[1].c_str());
            intensities.push_back(frequencies[currentIntensity] * transDipole * transDipole);
          }
      }

    if (mol.NumAtoms() == 0) { // e.g., if we're at the end of a file PR#1737209
      mol.EndModify();
      return false;
    }

    if (!pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.ConnectTheDots();
    if (!pConv->IsOption("s",OBConversion::INOPTIONS) 
        && !pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();

    mol.EndModify();

    if (hasPartialCharges)
      {
        mol.SetPartialChargesPerceived();
        FOR_ATOMS_OF_MOL(atom, mol) {
          atom->SetPartialCharge(charges[atom->GetIdx()-1]); // atom index issue
        }

        // Annotate that partial charges come from MOPAC Mulliken
        OBPairData *dp = new OBPairData;
        dp->SetAttribute("PartialCharges");
        dp->SetValue("Mulliken");
        dp->SetOrigin(fileformatInput);
        mol.SetData(dp);
      }
    if (dipoleMoment)
      mol.SetData(dipoleMoment);
    if (frequencies.size() != 0) { // we found some vibrations
      OBVibrationData *vd = new OBVibrationData;
      vd->SetData(displacements, frequencies, intensities);
      vd->SetOrigin(fileformatInput);
      mol.SetData(vd);
    }

    // Attach unit cell translation vectors if found
    if (numTranslationVectors > 0) {
      OBUnitCell* uc = new OBUnitCell;
      uc->SetData(translationVectors[0], translationVectors[1], translationVectors[2]);
      uc->SetOrigin(fileformatInput);
      mol.SetData(uc);
    }

    mol.SetTitle(title);

    return(true);
  }

  //************************************************************
  class MOPACCARTFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    MOPACCARTFormat()
    {
      OBConversion::RegisterFormat("mopcrt",this, "chemical/x-mopac-input");
      OBConversion::RegisterFormat("mop",this, "chemical/x-mopac-input");
      OBConversion::RegisterFormat("mpc",this, "chemical/x-mopac-input");
      // Command-line keywords
      OBConversion::RegisterOptionParam("k", NULL, 1, OBConversion::OUTOPTIONS);
      // Command-line keyword file
      OBConversion::RegisterOptionParam("f", NULL, 1, OBConversion::OUTOPTIONS);
    }

    virtual const char* Description() //required
    {
      return
        "MOPAC Cartesian format\n"
        "Read Options e.g. -as\n"
        "  s  Output single bonds only\n"
        "  b  Disable bonding entirely\n"
        "Write Options e.g. -xk\n"
        "  k  \"keywords\" Use the specified keywords for input\n"
        "  f    <file>     Read the file specified for input keywords\n"
        "  u               Write the crystallographic unit cell, if present.\n\n";
    };

    virtual const char* GetMIMEType() 
    { return "chemical/x-mopac-input"; };

    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

    ////////////////////////////////////////////////////
  };

  //Make an instance of the format class
  MOPACCARTFormat theMOPACCARTFormat;

  /////////////////////////////////////////////////////////////////
  bool MOPACCARTFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title= pConv->GetTitle();

    char buffer[BUFF_SIZE];
    string str;
    double x,y,z;
    OBAtom *atom;
    vector<string> vs;

    ifs.getline(buffer,BUFF_SIZE); // keywords
    ifs.getline(buffer,BUFF_SIZE); // filename
    ifs.getline(buffer,BUFF_SIZE); // title (currently ignored)

    mol.BeginModify();

    while	(ifs.getline(buffer,BUFF_SIZE))
      {
        tokenize(vs,buffer);
        if (vs.size() == 0)
          break;
        else if (vs.size() < 7)
          return false;
        atom = mol.NewAtom();
        x = atof((char*)vs[1].c_str());
        y = atof((char*)vs[3].c_str());
        z = atof((char*)vs[5].c_str());
        atom->SetVector(x,y,z); //set coordinates

        //set atomic number
        atom->SetAtomicNum(etab.GetAtomicNum(vs[0].c_str()));
      }

    if (!pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.ConnectTheDots();
    if (!pConv->IsOption("s",OBConversion::INOPTIONS) &&
        !pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();
    mol.SetTitle(title);

    mol.EndModify();

    return(true);
  }

  ////////////////////////////////////////////////////////////////

  bool MOPACCARTFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

//    unsigned int i;
    char buffer[BUFF_SIZE];

    const char *keywords = pConv->IsOption("k",OBConversion::OUTOPTIONS);
    const char *keywordFile = pConv->IsOption("f",OBConversion::OUTOPTIONS);
    bool writeUnitCell = pConv->IsOption("u", OBConversion::OUTOPTIONS);
    string defaultKeywords = "PUT KEYWORDS HERE";

    if(keywords)
      {
        defaultKeywords = keywords;
      }

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
    else
      ofs << defaultKeywords << endl;

    ofs << mol.GetTitle() << endl;
    ofs << endl; // comment

    string str,str1;
    FOR_ATOMS_OF_MOL(atom, mol)
      {
        snprintf(buffer,BUFF_SIZE,"%-3s%8.5f 1 %8.5f 1 %8.5f 1",
                 etab.GetSymbol(atom->GetAtomicNum()),
                 atom->GetX(),
                 atom->GetY(),
                 atom->GetZ());
        ofs << buffer << "\n";
      }

    OBUnitCell *uc = (OBUnitCell*)mol.GetData(OBGenericDataType::UnitCell);
    if (uc && writeUnitCell) {
      uc->FillUnitCell(&mol); // complete the unit cell with symmetry-derived atoms

      vector<vector3> cellVectors = uc->GetCellVectors();
      for (vector<vector3>::iterator i = cellVectors.begin(); i != cellVectors.end(); ++i) {
        snprintf(buffer,BUFF_SIZE,"Tv %8.5f 1 %8.5f 1 %8.5f 1",
                 i->x(),
                 i->y(),
                 i->z());
        ofs << buffer << "\n";
      }
    }

    return(true);
  }

  //************************************************************
  class MOPACINTFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    MOPACINTFormat()
    {
      OBConversion::RegisterFormat("mopin", this, "chemical/x-mopac-input");
      // Command-line keywords
      OBConversion::RegisterOptionParam("k", NULL, 1, OBConversion::OUTOPTIONS);
      // Command-line keyword file
      OBConversion::RegisterOptionParam("f", NULL, 1, OBConversion::OUTOPTIONS);
    }

    virtual const char* Description() //required
    {
      return "MOPAC Internal\n"
        "Write Options e.g. -xk\n"
        "  k  \"keywords\" Use the specified keywords for input\n"
        "  f    <file>     Read the file specified for input keywords\n\n";
    };

    virtual const char* GetMIMEType() 
    { return "chemical/x-mopac-input"; };

    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
  };

  //Make an instance of the format class
  MOPACINTFormat theMOPACINTFormat;

  /////////////////////////////////////////////////////////////////
  bool MOPACINTFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title= pConv->GetTitle();

    char buffer[BUFF_SIZE];
    OBAtom *atom;
    vector<string> vs;
    
    vector<OBInternalCoord*> vic;
    vector<int> indices;
    vic.push_back((OBInternalCoord*)NULL);
    
    ifs.getline(buffer,BUFF_SIZE); // keywords
    ifs.getline(buffer,BUFF_SIZE); // filename
    ifs.getline(buffer,BUFF_SIZE); // title (currently ignored)
    
    mol.BeginModify();

    while (ifs.getline(buffer,BUFF_SIZE)) {
      tokenize(vs,buffer);
      if (vs.size() == 0)
        break;
      else if (vs.size() < 10)
        return false;
      atom = mol.NewAtom();
        
      OBInternalCoord *coord = new OBInternalCoord;
      //vic[atom->GetIdx()]->_dst = atof(vs[1].c_str());
      //vic[atom->GetIdx()]->_ang = atof(vs[3].c_str());
      //vic[atom->GetIdx()]->_tor = atof(vs[5].c_str());
      coord->_dst = atof(vs[1].c_str());
      coord->_ang = atof(vs[3].c_str());
      coord->_tor = atof(vs[5].c_str());
      vic.push_back(coord);

      indices.push_back(atoi(vs[7].c_str()));
      indices.push_back(atoi(vs[8].c_str()));
      indices.push_back(atoi(vs[9].c_str()));
 
      atom->SetAtomicNum(etab.GetAtomicNum(vs[0].c_str()));
    }
    
    int idx = 0;
    FOR_ATOMS_OF_MOL (a, mol) {
      if ((indices[idx] > 0) && (indices[idx] <= mol.NumAtoms()))
        vic[a->GetIdx()]->_a = mol.GetAtom(indices[idx]);
      else
        vic[a->GetIdx()]->_a = NULL;
      
      if ((indices[idx+1] > 0) && (indices[idx+1] <= mol.NumAtoms()))
        vic[a->GetIdx()]->_b = mol.GetAtom(indices[idx+1]);
      else
        vic[a->GetIdx()]->_b = NULL;

      if ((indices[idx+2] > 0) && (indices[idx+2] <= mol.NumAtoms()))
        vic[a->GetIdx()]->_c = mol.GetAtom(indices[idx+2]);
      else
        vic[a->GetIdx()]->_c = NULL;
      
      idx += 3;
    }

    /* 
       vector<OBInternalCoord*>::iterator j;
       for (j = vic.begin(); j != vic.end(); j++) {
       cout << (*j)->_dst << " " << (*j)->_ang << " " << (*j)->_tor << " ";
       if ((*j)->_a)
       cout << (*j)->_a->GetIdx() << " "; 
       if ((*j)->_b)
       cout << (*j)->_b->GetIdx() << " ";
       if ((*j)->_c)
       cout << (*j)->_c->GetIdx() << endl;
       }
    */
    InternalToCartesian(vic,mol);

    if (!pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.ConnectTheDots();
    if (!pConv->IsOption("s",OBConversion::INOPTIONS) && !pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();
    mol.SetTitle(title);

    return(true);
  }

  /////////////////////////////////////////////////////////////////
  bool MOPACINTFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    char type[16], buffer[BUFF_SIZE];
    OBAtom *a,*b,*c;

    vector<OBInternalCoord*> vic;
    vic.push_back((OBInternalCoord*)NULL);
    
    for (unsigned int i = 0; i<mol.NumAtoms(); i++)
      vic.push_back(new OBInternalCoord);

    CartesianToInternal(vic,mol);

    const char *keywords = pConv->IsOption("k",OBConversion::OUTOPTIONS);
    const char *keywordFile = pConv->IsOption("f",OBConversion::OUTOPTIONS);
    string defaultKeywords = "PUT KEYWORDS HERE";

    if(keywords)
      {
        defaultKeywords = keywords;
      }

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
    else
      ofs << defaultKeywords << endl;

    ofs << mol.GetTitle() << endl;
    ofs << endl; // comment

    double r,w,t;
    FOR_ATOMS_OF_MOL (atom, mol) {
      a = vic[atom->GetIdx()]->_a;
      b = vic[atom->GetIdx()]->_b;
      c = vic[atom->GetIdx()]->_c;
      r = vic[atom->GetIdx()]->_dst;
      w = vic[atom->GetIdx()]->_ang;
      t = vic[atom->GetIdx()]->_tor;
	
      strncpy(type, etab.GetSymbol(atom->GetAtomicNum()), 16);
      type[15] = '\0';

      if (t < 0)
        t += 360;
      snprintf(buffer, BUFF_SIZE, "%-2s %10.6f  1  %10.6f  1  %10.6f  1  ", type, r, w, t);
      ofs << buffer;
      if (atom->GetIdx() == 1) 
        snprintf(buffer, BUFF_SIZE, "%4d%4d%4d\n", 0, 0, 0);
      if (atom->GetIdx() == 2) 
        snprintf(buffer, BUFF_SIZE, "%4d%4d%4d\n", a->GetIdx(), 0, 0);
      if (atom->GetIdx() == 3) 
        snprintf(buffer, BUFF_SIZE, "%4d%4d%4d\n", a->GetIdx(), b->GetIdx(), 0);
      if (atom->GetIdx() >= 4) 
        snprintf(buffer, BUFF_SIZE, "%4d%4d%4d\n", a->GetIdx(), b->GetIdx(), c->GetIdx());
      ofs << buffer;
    }

    return(true);
  }


} //namespace OpenBabel
