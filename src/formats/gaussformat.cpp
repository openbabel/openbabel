/**********************************************************************
Copyright (C) 2000 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2010 by Geoffrey R. Hutchison
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

  class GaussianOutputFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    GaussianOutputFormat()
    {
      OBConversion::RegisterFormat("gal",this, "chemical/x-gaussian-log");
      OBConversion::RegisterFormat("g92",this);
      OBConversion::RegisterFormat("g94",this);
      OBConversion::RegisterFormat("g98",this);
      OBConversion::RegisterFormat("g03",this);
      OBConversion::RegisterFormat("g09",this);
    }

    virtual const char* Description() //required
    {
      return
        "Gaussian Output\n"
        "Read Options e.g. -as\n"
        "  s  Output single bonds only\n"
        "  b  Disable bonding entirely\n\n";
    };

    virtual const char* SpecificationURL()
    { return "http://www.gaussian.com/";};

    virtual const char* GetMIMEType()
    { return "chemical/x-gaussian-log"; };

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
      return READONEONLY | NOTWRITABLE;
    };

    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
  };

  //Make an instance of the format class
  GaussianOutputFormat theGaussianOutputFormat;

  class GaussianInputFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    GaussianInputFormat()
    {
      OBConversion::RegisterFormat("com",this, "chemical/x-gaussian-input");
      OBConversion::RegisterFormat("gau",this);
      OBConversion::RegisterFormat("gjc",this);
      OBConversion::RegisterFormat("gjf",this);
      OBConversion::RegisterOptionParam("b", NULL, 0, OBConversion::OUTOPTIONS);
      // Command-line keywords
      OBConversion::RegisterOptionParam("k", NULL, 1, OBConversion::OUTOPTIONS);
      // Command-line keyword file
      OBConversion::RegisterOptionParam("f", NULL, 1, OBConversion::OUTOPTIONS);    }

    virtual const char* Description() //required
    {
      return
        "Gaussian 98/03 Input\n"
        "Write Options e.g. -xk\n"
        "  b               Output includes bonds\n"
        "  k  \"keywords\" Use the specified keywords for input\n"
        "  f    <file>     Read the file specified for input keywords\n"
        "  u               Write the crystallographic unit cell, if present.\n\n";
    };

    virtual const char* SpecificationURL()
    {return "http://www.gaussian.com/g_ur/m_input.htm";};

    virtual const char* GetMIMEType()
    { return "chemical/x-gaussian-input"; };

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
      return NOTREADABLE | WRITEONEONLY;
    };

    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

  };

  //Make an instance of the format class
  GaussianInputFormat theGaussianInputFormat;

  ////////////////////////////////////////////////////////////////

  bool GaussianInputFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    char buffer[BUFF_SIZE];
    const char *keywords = pConv->IsOption("k",OBConversion::OUTOPTIONS);
    const char *keywordsEnable = pConv->IsOption("k",OBConversion::GENOPTIONS);
    const char *keywordFile = pConv->IsOption("f",OBConversion::OUTOPTIONS);
    bool writeUnitCell = (NULL != pConv->IsOption("u", OBConversion::OUTOPTIONS));
    string defaultKeywords = "#Put Keywords Here, check Charge and Multiplicity.";

    if(keywords)
      {
        defaultKeywords = keywords;
      }

    if (keywordsEnable)
      {
        string model;
        string basis;
        string method;

        OBPairData *pd = (OBPairData *) pmol->GetData("model");
        if(pd)
          model = pd->GetValue();

        pd = (OBPairData *) pmol->GetData("basis");
        if(pd)
          basis = pd->GetValue();

        pd = (OBPairData *) pmol->GetData("method");
        if(pd)
          method = pd->GetValue();

        if(method == "optimize")
          {
            method = "opt";
          }

        if(model != "" && basis != "" && method != "")
          {
            ofs << model << "/" << basis << "," << method << endl;
          }
        else
          {
            ofs << "#Unable to translate keywords!" << endl;
            ofs << defaultKeywords << endl;
          }
      }
    else if (keywordFile)
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
      {
        ofs << defaultKeywords << endl;
      }
    ofs << endl; // blank line after keywords
    ofs << " " << mol.GetTitle() << endl << endl;

    snprintf(buffer, BUFF_SIZE, "%d  %d",
             mol.GetTotalCharge(),
             mol.GetTotalSpinMultiplicity());
    ofs << buffer << endl;

    FOR_ATOMS_OF_MOL(atom, mol)
      {
        if (atom->GetIsotope() == 0)
          snprintf(buffer, BUFF_SIZE, "%-3s      %10.5f      %10.5f      %10.5f",
                   etab.GetSymbol(atom->GetAtomicNum()),
                   atom->GetX(), atom->GetY(), atom->GetZ());
        else
          snprintf(buffer, BUFF_SIZE, "%-3s(Iso=%d) %10.5f      %10.5f      %10.5f",
                   etab.GetSymbol(atom->GetAtomicNum()),
                   atom->GetIsotope(),
                   atom->GetX(), atom->GetY(), atom->GetZ());

        ofs << buffer << endl;
      }
    // Translation vectors
    OBUnitCell *uc = (OBUnitCell*)mol.GetData(OBGenericDataType::UnitCell);
    if (uc && writeUnitCell) {
      uc->FillUnitCell(&mol); // complete the unit cell with symmetry-derived atoms

      vector<vector3> cellVectors = uc->GetCellVectors();
      for (vector<vector3>::iterator i = cellVectors.begin(); i != cellVectors.end(); ++i) {
          snprintf(buffer, BUFF_SIZE, "TV       %10.5f      %10.5f      %10.5f",
                   i->x(),
                   i->y(),
                   i->z());
        ofs << buffer << '\n';
      }
    }

    // Bonds, contributed by Daniel Mansfield
    if (pConv->IsOption("b",OBConversion::OUTOPTIONS))
    {
      // first, make begin.GetIdx < end.GetIdx
      OBBond* bond;
      OBAtom *atom;
      vector<OBEdgeBase*>::iterator j;
      vector<OBNodeBase*>::iterator i;
      OBAtom *bgn, *end;
      for (bond = mol.BeginBond(j); bond; bond = mol.NextBond(j))
        {
          if (bond->GetBeginAtomIdx() > bond->GetEndAtomIdx()) {
            bgn = bond->GetBeginAtom();
            end = bond->GetEndAtom();
            bond->SetBegin(end);
            bond->SetEnd(bgn);
          }
        }

      // this seems inefficient -- perhaps using atom neighbor iterators?
      // -GRH
      for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
        {
          ofs << endl << atom->GetIdx() << " ";
          for (bond = mol.BeginBond(j); bond; bond = mol.NextBond(j))
            {
              if (bond->GetBeginAtomIdx() == atom->GetIdx()) {
                snprintf(buffer, BUFF_SIZE, "%d %1.1f ", bond->GetEndAtomIdx(), (float) bond->GetBondOrder());
                ofs << buffer;
              }
            }
        } // iterate through atoms
    } // end writing bonds

    // file should end with a blank line
    ofs << endl;
    return(true);
  }

  // Reading Gaussian output has been tested for G98 and G03 to some degree
  // If you have problems (or examples of older output), please contact
  // the openbabel-discuss@lists.sourceforge.net mailing list and/or post a bug
  bool GaussianOutputFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = pOb->CastAndClear<OBMol>();
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
    int charge = 0;
    unsigned int spin = 1;
    bool hasPartialCharges = false;
    string chargeModel; // descriptor for charges (e.g. "Mulliken")

    // coordinates of all steps
    // Set conformers to all coordinates we adopted
    std::vector<double*> vconf; // index of all frames/conformers
    std::vector<double> coordinates; // coordinates in each frame
    int natoms = 0; // number of atoms -- ensure we don't go to a new job with a different molecule

    // OBConformerData stores information about multiple steps
    // we can change attribute later if needed (e.g., IRC)
    OBConformerData *confData = new OBConformerData();
    confData->SetOrigin(fileformatInput);
    std::vector<unsigned short> confDimensions = confData->GetDimension(); // to be fair, set these all to 3D
    std::vector<double>         confEnergies   = confData->GetEnergies();
    std::vector< std::vector< vector3 > > confForces = confData->GetForces();

    //Vibrational data
    std::vector< std::vector< vector3 > > Lx;
    std::vector<double> Frequencies, Intensities;
    //Rotational data
    std::vector<double> RotConsts(3);
    int RotSymNum=1;
    OBRotationData::RType RotorType;

    // Translation vectors (if present)
    vector3 translationVectors[3];
    int numTranslationVectors = 0;

    //Electronic Excitation data
    std::vector<double> Forces, Wavelengths, EDipole,
      RotatoryStrengthsVelocity, RotatoryStrengthsLength;

    // Orbital data
    std::vector<double> orbitals;
    std::vector<std::string> symmetries;
    int aHOMO, bHOMO, betaStart;

    //Put some metadata into OBCommentData
    string comment("Gaussian ");
    ifs.getline(buffer,BUFF_SIZE);
    if(*buffer)
    {
      comment += strchr(buffer,'=')+2;
      comment += "";
      for(unsigned i=0; i<115, ifs; ++i)
      {
        ifs.getline(buffer,BUFF_SIZE);
        if(buffer[1]=='#')
        {
          //the line describing the method
          comment += buffer;
          OBCommentData *cd = new OBCommentData;
          cd->SetData(comment);
          cd->SetOrigin(fileformatInput);
          mol.SetData(cd);
          break;
        }
      }
    }

    int i=0;
    char coords_type[25];

    //Prescan file to find second instance of "orientation:"
    //This will be the kind of coords used in the chk/fchk file
    while (ifs.getline(buffer,BUFF_SIZE))
      {
        if (strstr(buffer,"orientation:") !=NULL)
           i++;
        if (i==2)
           break;
     }
    if (!ifs.eof())
    {
    tokenize (vs, buffer);
    strcpy (coords_type, vs[0].c_str());
    strcat (coords_type, " orientation:");
    ifs.seekg(0);  //rewind
    //    cout << coords_type << std::endl;
    }

    mol.BeginModify();
    while (ifs.getline(buffer,BUFF_SIZE))
      {
        if (strstr(buffer,"Multiplicity") != NULL)
          {
            tokenize(vs, buffer, " \t\n");
            if (vs.size() == 6)
              {
                charge = atoi(vs[2].c_str());
                spin = atoi(vs[5].c_str());
              }

            ifs.getline(buffer,BUFF_SIZE);
          }
        else if (strstr(buffer, coords_type) != NULL)
          {
            numTranslationVectors = 0; // ignore old translationVectors
            ifs.getline(buffer,BUFF_SIZE);      // ---------------
            ifs.getline(buffer,BUFF_SIZE);      // column headings
            ifs.getline(buffer,BUFF_SIZE);	// column headings
            ifs.getline(buffer,BUFF_SIZE);	// ---------------
            ifs.getline(buffer,BUFF_SIZE);
            tokenize(vs,buffer);
            while (vs.size() == 6)
              {
                x = atof((char*)vs[3].c_str());
                y = atof((char*)vs[4].c_str());
                z = atof((char*)vs[5].c_str());
                int atomicNum = atoi((char*)vs[1].c_str());

                if (atomicNum > 0) // translation vectors are "-2"
                  {
                    if (natoms == 0) { // first time reading the molecule, create each atom
                      atom = mol.NewAtom();
                      atom->SetAtomicNum(atoi((char*)vs[1].c_str()));
                    }
                    coordinates.push_back(x);
                    coordinates.push_back(y);
                    coordinates.push_back(z);
                  }
                else {
                  translationVectors[numTranslationVectors++].Set(x, y, z);
                }

                if (!ifs.getline(buffer,BUFF_SIZE)) {
                  break;
                }
                tokenize(vs,buffer);
              }
            // done with reading atoms
            natoms = mol.NumAtoms();
            // malloc / memcpy
            double *tmpCoords = new double [(natoms)*3];
            memcpy(tmpCoords, &coordinates[0], sizeof(double)*natoms*3);
            vconf.push_back(tmpCoords);
            coordinates.clear();
            confDimensions.push_back(3); // always 3D -- OBConformerData allows mixing 2D and 3D structures
          }
        else if(strstr(buffer,"Dipole moment") != NULL)
            {
              ifs.getline(buffer,BUFF_SIZE); // actual components   X ###  Y #### Z ###
              tokenize(vs,buffer);
              if (vs.size() >= 6)
                {
                  OBVectorData *dipoleMoment = new OBVectorData;
                  dipoleMoment->SetAttribute("Dipole Moment");
                  double x, y, z;
                  x = atof(vs[1].c_str());
                  y = atof(vs[3].c_str());
                  z = atof(vs[5].c_str());
                  dipoleMoment->SetData(x, y, z);
                  dipoleMoment->SetOrigin(fileformatInput);
                  mol.SetData(dipoleMoment);
                }
              if (!ifs.getline(buffer,BUFF_SIZE)) break;
            }
        else if(strstr(buffer,"Total atomic charges") != NULL ||
                strstr(buffer,"Mulliken atomic charges") != NULL)
          {
            hasPartialCharges = true;
            chargeModel = "Mulliken";
            ifs.getline(buffer,BUFF_SIZE);	// column headings
            ifs.getline(buffer,BUFF_SIZE);
            tokenize(vs,buffer);
            while (vs.size() >= 3 &&
                   strstr(buffer,"Sum of ") == NULL)
              {
                atom = mol.GetAtom(atoi(vs[0].c_str()));
                if (!atom)
                  break;
                atom->SetPartialCharge(atof(vs[2].c_str()));

                if (!ifs.getline(buffer,BUFF_SIZE)) break;
                tokenize(vs,buffer);
              }
          }
        else if (strstr(buffer, "Charges from ESP fit") != NULL)
          {
            hasPartialCharges = true;
            chargeModel = "ESP";
            ifs.getline(buffer,BUFF_SIZE);	// Charge / dipole line
            ifs.getline(buffer,BUFF_SIZE); // column header
            ifs.getline(buffer,BUFF_SIZE); // real charges
            tokenize(vs,buffer);
            while (vs.size() >= 3 &&
                   strstr(buffer,"-----") == NULL)
              {
                atom = mol.GetAtom(atoi(vs[0].c_str()));
                if (!atom)
                  break;
                atom->SetPartialCharge(atof(vs[2].c_str()));

                if (!ifs.getline(buffer,BUFF_SIZE)) break;
                tokenize(vs,buffer);
              }
          }
        else if(strstr(buffer, " Frequencies -- ")) //vibrational frequencies
        {
          //The info should appear only once as several blocks starting with this line
          tokenize(vs, buffer);
          for(unsigned int i=2; i<vs.size(); ++i)
            Frequencies.push_back(atof(vs[i].c_str()));
          ifs.getline(buffer,BUFF_SIZE); //Red. masses
          ifs.getline(buffer,BUFF_SIZE); //Frc consts
          ifs.getline(buffer,BUFF_SIZE); //IR Inten
          tokenize(vs, buffer);
          for(unsigned int i=3; i<vs.size(); ++i)
            Intensities.push_back(atof(vs[i].c_str()));

          ifs.getline(buffer, BUFF_SIZE); // column labels or Raman intensity
          if(strstr(buffer, "Raman Activ")) {
            ifs.getline(buffer, BUFF_SIZE); // Depolar (P)
            ifs.getline(buffer, BUFF_SIZE); // Depolar (U)
            ifs.getline(buffer, BUFF_SIZE); // column labels
          }
          ifs.getline(buffer, BUFF_SIZE); // actual displacement data
          tokenize(vs, buffer);
          vector<vector3> vib1, vib2, vib3;
          double x, y, z;
          while(vs.size() > 5) {
            for (unsigned int i = 2; i < vs.size()-2; i += 3) {
              x = atof(vs[i].c_str());
              y = atof(vs[i+1].c_str());
              z = atof(vs[i+2].c_str());

              if (i == 2)
                vib1.push_back(vector3(x, y, z));
              else if (i == 5)
                vib2.push_back(vector3(x, y, z));
              else if (i == 8)
                vib3.push_back(vector3(x, y, z));
            }

            if (!ifs.getline(buffer, BUFF_SIZE))
              break;
            tokenize(vs,buffer);
          }
          Lx.push_back(vib1);
          if (vib2.size())
            Lx.push_back(vib2);
          if (vib3.size())
            Lx.push_back(vib3);
        }

        else if(strstr(buffer, " This molecule is "))//rotational data
        {
          if(strstr(buffer, "asymmetric"))
            RotorType = OBRotationData::ASYMMETRIC;
          else if(strstr(buffer, "symmetric"))
            RotorType = OBRotationData::SYMMETRIC;
          else if(strstr(buffer, "linear"))
            RotorType = OBRotationData::LINEAR;
          else
             RotorType = OBRotationData::UNKNOWN;
          ifs.getline(buffer,BUFF_SIZE); //symmetry number
          tokenize(vs, buffer);
          RotSymNum = atoi(vs[3].c_str());
        }

        else if(strstr(buffer, "Rotational constant"))
        {
          tokenize(vs, buffer);
          RotConsts.clear();
          for(int i=3; i<vs.size(); ++i)
            RotConsts.push_back(atof(vs[i].c_str()));
        }

        else if(strstr(buffer, "alpha electrons")) // # of electrons / orbital
        {
          tokenize(vs, buffer);
          if (vs.size() == 6) {
            // # alpha electrons # beta electrons
            aHOMO = atoi(vs[0].c_str());
            bHOMO = atoi(vs[3].c_str());
          }
        }
        else if(strstr(buffer, "rbital symmetries")) // orbital symmetries
          {
            symmetries.clear();
            std::string label; // used as a temporary to remove "(" and ")" from labels
            int offset = 0;

            while(true) {
              ifs.getline(buffer, BUFF_SIZE);
              tokenize(vs, buffer); // parse first line "Occupied" ...
              for (unsigned int i = 1; i < vs.size(); ++i) {
                label = vs[i].substr(1, vs[i].length() - 2);
                symmetries.push_back(label);
              }
              ifs.getline(buffer, BUFF_SIZE);

              // Parse remaining lines
              while (strstr(buffer, "(")) {
                tokenize(vs, buffer);
                if (strstr(buffer, "Virtual")) {
                  offset = 1; // skip first token
                } else {
                  offset = 0;
                }
                for (unsigned int i = offset; i < vs.size(); ++i) {
                  label = vs[i].substr(1, vs[i].length() - 2);
                  symmetries.push_back(label);
                }
                ifs.getline(buffer, BUFF_SIZE); // get next line
              } // end parsing symmetry labels
              if (!strstr(buffer, "Beta")) // no beta orbitals
                break;
            } // end alpha/beta section
          }
        else if (strstr(buffer, "Alpha") && strstr(buffer, ". eigenvalues --")) {
          orbitals.clear();
          betaStart = 0;
          while (strstr(buffer, ". eigenvalues --")) {
            tokenize(vs, buffer);
            if (vs.size() < 4)
              break;
            if (vs[0].find("Beta") !=string::npos && betaStart == 0) // mark where we switch from alpha to beta
              betaStart = orbitals.size();
            for (unsigned int i = 4; i < vs.size(); ++i) {
              orbitals.push_back(atof(vs[i].c_str()));
            }
            ifs.getline(buffer, BUFF_SIZE);
          }
        }
        else if(strstr(buffer, " Excited State")) // Force and wavelength data
        {
          // The above line appears for each state, so just append the info to the vectors
          tokenize(vs, buffer);
          if (vs.size() == 9) {
            double wavelength = atof(vs[6].c_str());
            double force = atof(vs[8].substr(2).c_str());
            Forces.push_back(force);
            Wavelengths.push_back(wavelength);
          }
        }
        else if(strstr(buffer, " Ground to excited state Transition electric dipole moments (Au):"))
          // Electronic dipole moments
        {
          ifs.getline(buffer, BUFF_SIZE); // Headings
          ifs.getline(buffer, BUFF_SIZE); // First entry
          tokenize(vs, buffer);
          while (vs.size() == 5) {
            double s = atof(vs[4].c_str());
            EDipole.push_back(s);
            ifs.getline(buffer, BUFF_SIZE);
            tokenize(vs, buffer);
          }
        }
        else if(strstr(buffer, "       state          X           Y           Z     R(velocity)")) {
          // Rotatory Strengths
          ifs.getline(buffer, BUFF_SIZE); // First entry
          tokenize(vs, buffer);
          while (vs.size() == 5) {
            double s = atof(vs[4].c_str());
            RotatoryStrengthsVelocity.push_back(s);
            ifs.getline(buffer, BUFF_SIZE);
            tokenize(vs, buffer);
          }
        }
        else if(strstr(buffer, "       state          X           Y           Z     R(length)")) {
          // Rotatory Strengths
          ifs.getline(buffer, BUFF_SIZE); // First entry
          tokenize(vs, buffer);
          while (vs.size() == 5) {
            double s = atof(vs[4].c_str());
            RotatoryStrengthsLength.push_back(s);
            ifs.getline(buffer, BUFF_SIZE);
            tokenize(vs, buffer);
          }
        }

        else if (strstr(buffer, "Forces (Hartrees/Bohr)"))
          {
            ifs.getline(buffer, BUFF_SIZE); // column headers
            ifs.getline(buffer, BUFF_SIZE); // ------
            ifs.getline(buffer, BUFF_SIZE); // real data
          }

        else if (strstr(buffer, "Isotropic = ")) // NMR shifts
          {
            tokenize(vs, buffer);
            if (vs.size() >= 4)
              {
                atom = mol.GetAtom(atoi(vs[0].c_str()));
                OBPairData *nmrShift = new OBPairData();
                nmrShift->SetAttribute("NMR Isotropic Shift");

                string shift = vs[4].c_str();
                nmrShift->SetValue(shift);

                atom->SetData(nmrShift);
              }
          }

        else if(strstr(buffer,"SCF Done:") != NULL)
          {
#define HARTREE_TO_KCAL 627.509469
            tokenize(vs,buffer);
            mol.SetEnergy(atof(vs[4].c_str()) * HARTREE_TO_KCAL);
            confEnergies.push_back(mol.GetEnergy());
          }
/* Temporarily commented out until the handling of energy in OBMol is sorted out
        // MP2 energies also use a different syntax

        // PM3 energies use a different syntax
        else if(strstr(buffer,"E (Thermal)") != NULL)
          {
            ifs.getline(buffer,BUFF_SIZE); //Headers
            ifs.getline(buffer,BUFF_SIZE); //Total energy; what we want
            tokenize(vs,buffer);
            mol.SetEnergy(atof(vs[1].c_str()));
            confEnergies.push_back(mol.GetEnergy());
          }
*/
      } // end while

    if (mol.NumAtoms() == 0) { // e.g., if we're at the end of a file PR#1737209
      mol.EndModify();
      return false;
    }

    mol.EndModify();

    // Set conformers to all coordinates we adopted
    // but remove last geometry -- it's a duplicate
    if (vconf.size() > 1)
      vconf.pop_back();
    mol.SetConformers(vconf);
    mol.SetConformer(mol.NumConformers() - 1);
    // Copy the conformer data too
    confData->SetDimension(confDimensions);
    confData->SetEnergies(confEnergies);
    confData->SetForces(confForces);
    mol.SetData(confData);

    // Attach orbital data, if there is any
    if (orbitals.size() > 0)
      {
        OBOrbitalData *od = new OBOrbitalData;
        if (aHOMO == bHOMO) {
          od->LoadClosedShellOrbitals(orbitals, symmetries, aHOMO);
        } else {
          // we have to separate the alpha and beta vectors
          std::vector<double>      betaOrbitals;
          std::vector<std::string> betaSymmetries;
          int initialSize = orbitals.size();
          for (unsigned int i = betaStart; i < initialSize; ++i) {
            betaOrbitals.push_back(orbitals[i]);
            betaSymmetries.push_back(symmetries[i]);
          }
          // ok, now erase the end elements of orbitals and symmetries
          for (unsigned int i = betaStart; i < initialSize; ++i) {
            orbitals.pop_back();
            symmetries.pop_back();
          }
          // and load the alphas and betas
          od->LoadAlphaOrbitals(orbitals, symmetries, aHOMO);
          od->LoadBetaOrbitals(betaOrbitals, betaSymmetries, bHOMO);
        }
        od->SetOrigin(fileformatInput);
        mol.SetData(od);
      }

    //Attach vibrational data, if there is any, to molecule
    if(Frequencies.size()>0)
    {
      OBVibrationData* vd = new OBVibrationData;
      vd->SetData(Lx, Frequencies, Intensities);
      vd->SetOrigin(fileformatInput);
      mol.SetData(vd);
    }
    //Attach rotational data, if there is any, to molecule
    if(RotConsts[0]!=0.0)
    {
      OBRotationData* rd = new OBRotationData;
      rd->SetData(RotorType, RotConsts, RotSymNum);
      rd->SetOrigin(fileformatInput);
      mol.SetData(rd);
    }
    // Attach unit cell translation vectors if found
    if (numTranslationVectors > 0) {
      OBUnitCell* uc = new OBUnitCell;
      uc->SetData(translationVectors[0], translationVectors[1], translationVectors[2]);
      uc->SetOrigin(fileformatInput);
      mol.SetData(uc);
    }
    //Attach electronic transition data, if there is any, to molecule
    if(Forces.size() > 0 && Forces.size() == Wavelengths.size())
    {
      OBElectronicTransitionData* etd = new OBElectronicTransitionData;
      etd->SetData(Wavelengths, Forces);
      if (EDipole.size() == Forces.size())
        etd->SetEDipole(EDipole);
      if (RotatoryStrengthsLength.size() == Forces.size())
        etd->SetRotatoryStrengthsLength(RotatoryStrengthsLength);
      if (RotatoryStrengthsVelocity.size() == Forces.size())
        etd->SetRotatoryStrengthsVelocity(RotatoryStrengthsVelocity);
      etd->SetOrigin(fileformatInput);
      mol.SetData(etd);
    }

    if (!pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.ConnectTheDots();
    if (!pConv->IsOption("s",OBConversion::INOPTIONS) && !pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();

    if (hasPartialCharges) {
      mol.SetPartialChargesPerceived();

      // Annotate that partial charges come from Mulliken
      OBPairData *dp = new OBPairData;
      dp->SetAttribute("PartialCharges");
      dp->SetValue(chargeModel); // Mulliken, ESP, etc.
      dp->SetOrigin(fileformatInput);
      mol.SetData(dp);
    }
    mol.SetTotalCharge(charge);
    mol.SetTotalSpinMultiplicity(spin);

    mol.SetTitle(title);
    return(true);
  }

} //namespace OpenBabel
