/**********************************************************************
  Copyright (C) 2000 by OpenEye Scientific Software, Inc.
  Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
  Some portions Copyright (C) 2004 by Chris Morley
  Some portions Copyright (C) 2006 by Donald E. Curtis
  Some portions Copyright (C) 2009-2010 by Konstantin L. Tokarev

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

#include <algorithm>

using namespace std;
namespace Gamess {
}

namespace OpenBabel {

#define BOHR_TO_ANGSTROM 0.529177249
#define ANGSTROM_TO_BOHR 1.889725989

  class GAMESSOutputFormat: public OBMoleculeFormat {
    public:
      // Register this format type ID
      GAMESSOutputFormat() {
        OBConversion::RegisterFormat("gam",    this, "chemical/x-gamess-output");
        OBConversion::RegisterFormat("gamout", this);
        OBConversion::RegisterFormat("gamess", this);
      }

      // Required
      virtual const char* Description() {
        return
          "GAMESS Output\n"
          "Read Options e.g. -as\n"
          "  s  Output single bonds only\n"
          "  b  Disable bonding entirely\n"
          "  c  Read multiple conformers\n\n";
      }

      // Optional
      virtual const char* SpecificationURL() {
        return "http://www.msg.ameslab.gov/GAMESS/doc.menu.html";
      }

      virtual const char* GetMIMEType() {
        return "chemical/x-gamess-output";
      }

      // Flags() return can be any of the following combined by | or omitted
      // if none apply
      // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
      virtual unsigned int Flags() {
        return READONEONLY | NOTWRITABLE;
      }

      virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);

    private:
      //! \brief Parse GAMESS options section.
      //void ParseSection(char *tag, OBSetData *set, istream &ifs);
  };

  // Make an instance of the format class
  GAMESSOutputFormat theGAMESSOutputFormat;

  class GAMESSInputFormat: public OBMoleculeFormat {
    public:
      // Register this format type ID
      GAMESSInputFormat() {
        OBConversion::RegisterFormat("inp",   this, "chemical/x-gamess-input");
        OBConversion::RegisterFormat("gamin", this);
        // Command-line keywords
        OBConversion::RegisterOptionParam("k", NULL, 1, OBConversion::OUTOPTIONS);
        // Command-line keyword file
        OBConversion::RegisterOptionParam("f", NULL, 1, OBConversion::OUTOPTIONS);
      }

      // Required
      virtual const char* Description() {
        return
          "GAMESS Input\n"
          "Write Options e.g. -xk\n"
          "  k  \"keywords\" Use the specified keywords for input\n"
          "  f    <file>     Read the file specified for input keywords\n\n";
      }

      // Optional
      virtual const char* SpecificationURL() {
        return "http://www.msg.ameslab.gov/GAMESS/doc.menu.html";
      }

      virtual const char* GetMIMEType() {
        return "chemical/x-gamess-input";
      }

      // Flags() return can be any of the following combined by | or omitted
      // if none apply
      // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
      virtual unsigned int Flags() {
        return WRITEONEONLY; // | NOTREADABLE;
      }

      ////////////////////////////////////////////////////
      /// The "API" interface functions
      virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
      virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
  };

  // Make an instance of the format class
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
   curset->SetOrigin(fileformatInput);
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
   data->SetOrigin(fileformatInput);

   curset->AddData(data);
   }
   }
   }
  */

  bool GAMESSOutputFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv) {
    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if (pmol == NULL)
      return false;

    // Define some references so we can use the old parameter names
    istream& ifs = *pConv->GetInStream();
    OBMol& mol   = *pmol;

    const char* title = pConv->GetTitle();

    char buffer[BUFF_SIZE];
    string str, str1;
    double x, y, z;
    OBAtom* atom;
    vector<string> vs;
    bool hasPartialCharges = false;
    int aHOMO = 0;
    int bHOMO = 0;
    vector<double> orbitals;
    vector<std::string> symmetries;

    // Coordinates of all steps
    // Set conformers to all coordinates we adopted
    std::vector<double*> vconf;      // index of all frames/conformers
    std::vector<double> coordinates; // coordinates in each frame
    int natoms      = 0;
    int ndummyatoms = 0;

    // OBConformerData stores information about multiple steps
    OBConformerData* confData = new OBConformerData();
    confData->SetOrigin(fileformatInput);
    std::vector<unsigned short> confDimensions   = confData->GetDimension(); // to be fair, set these all to 3D
    std::vector<double>         confEnergies     = confData->GetEnergies();
    std::vector< std::vector< vector3 > > confForces = confData->GetForces();

    vector<double> frequencies, intensities, raman_intensities;
    vector< vector<vector3> > displacements;
    int lowFreqModesBegin;           // the number of the first low frequency mode
    int lowFreqModesEnd;             // the number of the last low frequency mode
    int numFreq, numIntens, numDisp; // GAMESS prints rotations & transl., which we ignore
    numFreq = numIntens = numDisp = 0;
    int charge = 0;
    int mult   = 1;

    // Must build generic data while we parse then add at the end.
    OBSetData* gmsset = new OBSetData();
    gmsset->SetAttribute("gamess");
    gmsset->SetOrigin(fileformatInput);

    mol.Clear();
    mol.BeginModify();
    while (ifs.getline(buffer, BUFF_SIZE)) {

      if (strstr(buffer, "ICHARG=")) {
        tokenize(vs, (strstr(buffer, "ICHARG=")));
        charge=atoi(vs[1].c_str());
      }

      if (strstr(buffer, "MULT ")) {
        tokenize(vs, (strstr(buffer, "MULT ")));
        mult=atoi(vs[2].c_str());
      }

      if (strstr(buffer, "ATOMIC                      COORDINATES (BOHR)") != NULL) {
       ifs.getline(buffer, BUFF_SIZE); // column headings
       ifs.getline(buffer, BUFF_SIZE);
       tokenize(vs, buffer);
       while (vs.size() == 5) {
         // Parse the current one
         int atomicNum = atoi(vs[1].c_str());
         x = atof((char*) vs[2].c_str()) * BOHR_TO_ANGSTROM;
         y = atof((char*) vs[3].c_str()) * BOHR_TO_ANGSTROM;
         z = atof((char*) vs[4].c_str()) * BOHR_TO_ANGSTROM;
         // First time reading the molecule, create each atom
         if (natoms == 0) {
           atom = mol.NewAtom();
           atom->SetAtomicNum(atomicNum);
         }
         coordinates.push_back(x);
         coordinates.push_back(y);
         coordinates.push_back(z);

         if (!ifs.getline(buffer, BUFF_SIZE))
           break;
         tokenize(vs, buffer);
       }
       // Done with reading atoms
       natoms = mol.NumAtoms();
       // malloc / memcpy
       double* tmpCoords = new double [(natoms)*3];
       memcpy(tmpCoords, &coordinates[0], sizeof(double)*natoms*3);
       vconf.push_back(tmpCoords);
       coordinates.clear();
       confDimensions.push_back(3); // always 3D -- OBConformerData allows mixing 2D and 3D structures

      } else if (strstr(buffer, "MULTIPOLE COORDINATES, ELECTRONIC AND NUCLEAR CHARGES") != NULL) {
        /*This set of EFP coordinates belongs only to the
         * conformer directly above this (ATOMIC   COORDINATES (BOHR))
         */
        double* tmpCoords = vconf.at(0);
        for (int i=0; i < natoms*3; i++)
          coordinates.push_back(tmpCoords[i]);

        ifs.getline(buffer, BUFF_SIZE); // column headings
        ifs.getline(buffer, BUFF_SIZE);
        ifs.getline(buffer, BUFF_SIZE);
        tokenize(vs, buffer);
        while (vs.size() == 6) {
          int atomicNum;
          /* For the included EFP1 potentials,
           * the atom name may start with "Z"
           * or have a non-zero nuclear charge
           */
          if (atof((char*) vs[5].c_str()) > 0.0) {
            atomicNum = OBElements::GetAtomicNum(vs[0].substr(0, 1).c_str());
            // First time reading the molecule, create each atom
            if (natoms == 0) {
              atom = mol.NewAtom();
              atom->SetAtomicNum(atomicNum);
            }
            x = atof((char*) vs[1].c_str())* BOHR_TO_ANGSTROM;
            y = atof((char*) vs[2].c_str())* BOHR_TO_ANGSTROM;
            z = atof((char*) vs[3].c_str())* BOHR_TO_ANGSTROM;
            coordinates.push_back(x);
            coordinates.push_back(y);
            coordinates.push_back(z);
          } else if (vs[0].substr(0, 1) == "Z") {
            atomicNum = OBElements::GetAtomicNum(vs[0].substr(1, 1).c_str());
            x = atof((char*) vs[1].c_str())* BOHR_TO_ANGSTROM;
            y = atof((char*) vs[2].c_str())* BOHR_TO_ANGSTROM;
            z = atof((char*) vs[3].c_str())* BOHR_TO_ANGSTROM;
            coordinates.push_back(x);
            coordinates.push_back(y);
            coordinates.push_back(z);
          }
          if (!ifs.getline(buffer, BUFF_SIZE))
            break;
          tokenize(vs, buffer);
        }
        // Done with reading atoms
        ndummyatoms = mol.NumAtoms();
        // malloc / memcpy
        memcpy(tmpCoords, &coordinates[0], sizeof(double)*natoms*3);
        vconf[0] = tmpCoords;
        coordinates.clear();

      } else if (strstr(buffer, "COORDINATES OF ALL ATOMS ARE (ANGS)") != NULL) {
        ifs.getline(buffer, BUFF_SIZE); // column headings
        ifs.getline(buffer, BUFF_SIZE); // ---------------
        ifs.getline(buffer, BUFF_SIZE);
        tokenize(vs, buffer);
        while (vs.size() == 5) {
          // Parse the current one
          int atomicNum = atoi(vs[1].c_str());
          x = atof((char*) vs[2].c_str());
          y = atof((char*) vs[3].c_str());
          z = atof((char*) vs[4].c_str());
          // First time reading the molecule, create each atom
          if (natoms == 0) {
            atom = mol.NewAtom();
            atom->SetAtomicNum(atomicNum);
          }
          coordinates.push_back(x);
          coordinates.push_back(y);
          coordinates.push_back(z);

          if (!ifs.getline(buffer, BUFF_SIZE))
            break;
          tokenize(vs, buffer);
        }

        if (strstr(buffer, "COORDINATES OF FRAGMENT") != NULL) {
          ifs.getline(buffer, BUFF_SIZE); // column headings
          ifs.getline(buffer, BUFF_SIZE);
          ifs.getline(buffer, BUFF_SIZE);
          //ifs.getline(buffer, BUFF_SIZE); //FRAGNAME
          tokenize(vs, buffer);
          while (vs.size() > 0) {
            if (vs.size() == 1) {
              vector<string> vs2;
              char delim[] = "=";
              tokenize(vs2, buffer, delim);
            } else {
              /* For the included EFP1 potentials,
               * the atom name may start with "Z"
               */
              int atomicNum;
              if (vs[0].substr(0, 1) == "Z")
                atomicNum = OBElements::GetAtomicNum(vs[0].substr(1, 1).c_str());
              else
                atomicNum = OBElements::GetAtomicNum(vs[0].substr(0, 1).c_str());
              // First time reading the molecule, create each atom
              if (natoms == 0) {
                atom = mol.NewAtom();
                atom->SetAtomicNum(atomicNum);
              }
              x = atof((char*) vs[1].c_str());
              y = atof((char*) vs[2].c_str());
              z = atof((char*) vs[3].c_str());
              coordinates.push_back(x);
              coordinates.push_back(y);
              coordinates.push_back(z);
            }

            if (!ifs.getline(buffer, BUFF_SIZE))
              break;
            tokenize(vs, buffer);
          }
        }

        // Done with reading atoms
        natoms = mol.NumAtoms();
        // malloc / memcpy
        double* tmpCoords = new double [(natoms)*3];
        memcpy(tmpCoords, &coordinates[0], sizeof(double)*natoms*3);
        vconf.push_back(tmpCoords);
        coordinates.clear();
        confDimensions.push_back(3); // always 3D -- OBConformerData allows mixing 2D and 3D structures

      } else if ((strstr(buffer, "NSERCH=") != NULL) && (strstr(buffer, "ENERGY=") != NULL)) {
        char* tok = strtok(buffer, " ="); // my tokenize
        int n = 0;
        while (true) {
          tok = strtok(NULL, " =");
          if (tok == NULL)
            break;
          n++;
          if (n == 3) {
            confEnergies.push_back(atof(tok));
            break;
          }
        }

      } else if (strstr(buffer, "ELECTROSTATIC MOMENTS") != NULL) {
        ifs.getline(buffer, BUFF_SIZE); //-----
        ifs.getline(buffer, BUFF_SIZE); // blank line
        ifs.getline(buffer, BUFF_SIZE); // column headings
        ifs.getline(buffer, BUFF_SIZE); // point charges @todo
        ifs.getline(buffer, BUFF_SIZE); // column headings dipole moment
        ifs.getline(buffer, BUFF_SIZE);

        tokenize(vs, buffer);
        if (vs.size() == 4) {
          OBVectorData* dipoleMoment = new OBVectorData;
          dipoleMoment->SetAttribute("Dipole Moment");
          double x, y, z;
          x = atof(vs[0].c_str());
          y = atof(vs[1].c_str());
          z = atof(vs[2].c_str());
          dipoleMoment->SetData(x, y, z);
          dipoleMoment->SetOrigin(fileformatInput);
          mol.SetData(dipoleMoment);
        }

      } else if (strstr(buffer, "MOPAC CHARGES") != NULL) {
        hasPartialCharges = true;
        ifs.getline(buffer, BUFF_SIZE); // ---------------
        ifs.getline(buffer, BUFF_SIZE); // column headings
        ifs.getline(buffer, BUFF_SIZE);
        tokenize(vs, buffer);
        while (vs.size() == 4) {
          atom = mol.GetAtom(atoi(vs[0].c_str()));
          atom->SetPartialCharge(atof(vs[2].c_str()));

          if (!ifs.getline(buffer, BUFF_SIZE))
            break;
          tokenize(vs, buffer);
        }

      } else if (strstr(buffer, "TOTAL MULLIKEN") != NULL) {
        hasPartialCharges = true;
        ifs.getline(buffer, BUFF_SIZE); // column headings
        ifs.getline(buffer, BUFF_SIZE);
        tokenize(vs, buffer);
        // atom number, atomic symbol, mulliken pop, charge
        while (vs.size() >= 4) {
          int atomNb = atoi(vs[0].c_str());
          if (!atomNb)
            break;
          atom = mol.GetAtom(atomNb);
          atom->SetPartialCharge(atof(vs[3].c_str()));

          if (!ifs.getline(buffer, BUFF_SIZE))
            break;
          tokenize(vs, buffer);
        }

      } else if (strstr(buffer, "NUMBER OF OCCUPIED ORBITALS") != NULL) {
        tokenize(vs, buffer);
        if (vs.size() == 7)      // alpha
          aHOMO = atoi(vs[6].c_str());
        else if (vs.size() == 8) // beta
          bHOMO = atoi(vs[7].c_str());

      } else if (strstr(buffer, "TAKEN AS ROTATIONS AND TRANSLATIONS") != NULL) {
        tokenize(vs, buffer);
        if (vs.size() < 4)
          break;
        lowFreqModesBegin = atoi(vs[1].c_str());
        lowFreqModesEnd   = atoi(vs[3].c_str());

      } else if (strstr(buffer, "TOTAL ENERGY      =") != NULL) {
        tokenize(vs, buffer);
        if (vs.size() == 4)
          mol.SetEnergy(atof(vs[3].c_str()));

      } else if (strstr(buffer, "FREQUENCY:") != NULL) {
        tokenize(vs, buffer);
        for (unsigned int i=1; i < vs.size(); ++i) {
          if (vs[i] == "I") // artifact from previous imaginary frequency
            continue;
          ++numFreq;
          if (numFreq < lowFreqModesBegin) // imaginary frequency
            frequencies.push_back(-atof(vs[i].c_str()));
          if (numFreq > lowFreqModesEnd)
            frequencies.push_back(atof(vs[i].c_str()));
        }
        ifs.getline(buffer, BUFF_SIZE); // possibly symmetry or red. mass
        if (strstr(buffer, "SYMMETRY:") != NULL) {
          // parse the vibrational symmetry
          ifs.getline(buffer, BUFF_SIZE); // reduced mass
        }
        ifs.getline(buffer, BUFF_SIZE); // intensities
        tokenize(vs, buffer);
        for (unsigned int i=2; i < vs.size(); ++i) {
          ++numIntens;
          if (numIntens < lowFreqModesBegin || numIntens > lowFreqModesEnd)
            intensities.push_back(atof(vs[i].c_str()) * 42.255); // conver to km/mol
        }
        ifs.getline(buffer, BUFF_SIZE); // blank or Raman activitie
        if (strstr(buffer, "RAMAN") != NULL) {
          tokenize(vs, buffer);
          for (unsigned int i=2; i < vs.size(); ++i) {
            if (numIntens < lowFreqModesBegin || numIntens > lowFreqModesEnd)
              raman_intensities.push_back(atof(vs[i].c_str()));
          }
          ifs.getline(buffer, BUFF_SIZE); // DEPOLARIZATION
          ifs.getline(buffer, BUFF_SIZE); // blank
        }

        // Now real work -- read displacements
        unsigned int prevModeCount = displacements.size();
        unsigned int newModes = frequencies.size() - displacements.size();
        vector<vector3> displacement;
        for (unsigned int i=0; i < newModes; ++i) {
          displacements.push_back(displacement);
        }

        ifs.getline(buffer, BUFF_SIZE);
        tokenize(vs, buffer);
        int modeCount = vs.size() - 3;
        double massNormalization;
        vector<double> x, y, z;
        while (modeCount >= 1) {
          // 1/sqrt(atomic mass)
          atom = mol.GetAtom(atoi(vs[0].c_str()));
          massNormalization = 1 / sqrt( atom->GetAtomicMass() );

          x.clear();
          // Not a typo -- e.g., atom number, atom label, x, then data
          for (unsigned int i=3; i < vs.size(); ++i) {
            x.push_back(massNormalization * atof(vs[i].c_str()));
          }
          y.clear();
          ifs.getline(buffer, BUFF_SIZE);
          tokenize(vs, buffer);
          for (unsigned int i=1; i < vs.size(); ++i) {
            y.push_back(massNormalization * atof(vs[i].c_str()));
          }

          z.clear();
          ifs.getline(buffer, BUFF_SIZE);
          tokenize(vs, buffer);
          for (unsigned int i=1; i < vs.size(); ++i) {
            z.push_back(massNormalization * atof(vs[i].c_str()));
          }

          // OK, now we have x, y, z for all new modes for one atom
          if (displacements.size()) {
            numDisp = prevModeCount;
            for (unsigned int i=0; i < modeCount;  ++i) {
              if (i >= modeCount - newModes){
                displacements[numDisp++].push_back(vector3(x[i], y[i], z[i]));
              }
            }
          }

          // Next set of atoms
          ifs.getline(buffer, BUFF_SIZE);
          tokenize(vs, buffer);
          modeCount = vs.size() - 3;
        }

      } else if (strstr(buffer, "EIGENVECTORS") != NULL ||
                 strstr(buffer, "MOLECULAR ORBITALS") != NULL) {

        ifs.getline(buffer, BUFF_SIZE); // ------ line
        ifs.getline(buffer, BUFF_SIZE); // blank
        orbitals.clear();
        symmetries.clear();

        while (strstr(buffer, "END OF RHF CALCULATION") == NULL
               && strstr(buffer, "-------") == NULL) {

          ifs.getline(buffer, BUFF_SIZE); // orbitals!
          ifs.getline(buffer, BUFF_SIZE); // energies in hartree
          tokenize(vs, buffer);
          for (unsigned int i=0; i < vs.size(); ++i)
            orbitals.push_back(27.21 * atof(vs[i].c_str()));

          ifs.getline(buffer, BUFF_SIZE); // symmetries
          tokenize(vs, buffer);
          for (unsigned int i=0; i < vs.size(); ++i)
            symmetries.push_back(vs[i]);

          // Orbital coefficients
          while (ifs.getline(buffer, BUFF_SIZE)
                 && strlen(buffer)
                 && strstr(buffer, "END") == NULL
                 && strstr(buffer, "---") == NULL) {
          }
            if (!ifs.good())
              break;
        }

      } else if (strstr(buffer, "INPUT CARD> $")) {
        string attr, value;
        char* ptr;

        for ( ; ; ) {
          ptr = buffer + 14;
          tokenize(vs, ptr);

          if (vs.size() > 2) {
            OBSetData* curset = (OBSetData*) gmsset->GetData(vs[0]);
            if (!curset) {
              curset = new OBSetData();
              curset->SetAttribute(vs[0]);
              curset->SetOrigin(fileformatInput);
              gmsset->AddData(curset);
            }
            for (unsigned int i=1; i < vs.size() && vs[i].substr(0, 4) != "$END"; i++) {
              string::size_type loc = vs[i].find("=", 0);
              if (loc != string::npos) {
                OBPairData *data = new OBPairData();
                data->SetAttribute(vs[i].substr(0, loc));
                data->SetValue(vs[i].substr(loc + 1));
                data->SetOrigin(fileformatInput);
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

    if (mol.NumAtoms() == 0) { // e.g., if we're at the end of a file PR#1737209
      mol.EndModify();
      return false;
    }

    // Add OBOrbitalData
    if (orbitals.size() > 0) {
      OBOrbitalData* od = new OBOrbitalData();

      if (aHOMO == bHOMO) {
        od->LoadClosedShellOrbitals(orbitals, symmetries, aHOMO);
      }
      od->SetOrigin(fileformatInput);
      mol.SetData(od);
    }

    const char* keywordsEnable = pConv->IsOption("k", OBConversion::GENOPTIONS);

    if (keywordsEnable) {
      // Add our gamess set
      pmol->SetData(gmsset);

      // If we have basis set data we should set our global pair data
      OBSetData* cset = (OBSetData*) gmsset->GetData("CONTRL");
      OBSetData* bset = (OBSetData*) gmsset->GetData("BASIS");

      string model = "b3lyp";
      string basis;
      string method;

      if (cset) {
        OBPairData* pd = NULL;

        pd = (OBPairData*) cset->GetData("SCFTYP");
        if (pd) {
          if(pd->GetValue() == "RHF")
            model = "rhf";
        }

        pd = (OBPairData*) cset->GetData("DFTTYP");
        if (pd) {
          if(pd->GetValue() == "BLYP")
            model = "b3lyp";
        }

        pd = (OBPairData*) cset->GetData("MPLEVL");
        if (pd) {
          if(pd->GetValue() == "2")
            model = "mp2";
        }

        pd = (OBPairData*) cset->GetData("CCTYP");
        if (pd) {
          if (pd->GetValue() == "CCSD(T)")
            model = "ccsd(t)";
        }

        pd = (OBPairData*) cset->GetData("RUNTYP");
        if (pd) {
          string value = pd->GetValue();
          if (value == "GRADIENT"
              || value == "HESSIAN"
              || value == "RAMAN"
              || value == "OPTIMIZE"
              || value == "SADPOINT") {

            method = pd->GetValue();
            transform(method.begin(), method.end(), method.begin(), ::tolower);
          }
        }
      }

      if (bset) {
        OBPairData* gbasis = (OBPairData*) bset->GetData("GBASIS");
        OBPairData* ngauss = (OBPairData*) bset->GetData("NGAUSS");

        if (gbasis) {
          string value = gbasis->GetValue();

          if (value == "am1")
            model = "am1";
          else if (value == "pm3")
            model = "pm3";
          else if (ngauss) {
            if (value == "STO") {
              basis.clear();
              basis += "sto-";
              basis += ngauss->GetValue();
              basis += "g";
            } else if (ngauss->GetValue() == "3"
                       || ngauss->GetValue() == "6") {
              basis.clear();
              basis = ngauss->GetValue();
              basis += "-";
              basis += gbasis->GetValue().substr(1);
              basis += "G(d)";
            }
          }
        }
      }

      OBPairData* nd = NULL;
      if (model != "") {
        nd = new OBPairData();
        nd->SetAttribute("model");
        nd->SetValue(model);
        nd->SetOrigin(fileformatInput);
        pmol->SetData(nd);
      }

      if (basis != "") {
        nd = new OBPairData();
        nd->SetAttribute("basis");
        nd->SetValue(basis);
        nd->SetOrigin(fileformatInput);
        pmol->SetData(nd);
      }

      if (method != "") {
        nd = new OBPairData();
        nd->SetAttribute("method");
        nd->SetValue(method);
        nd->SetOrigin(fileformatInput);
        pmol->SetData(nd);
      }
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

    if (!pConv->IsOption("b", OBConversion::INOPTIONS))
      mol.ConnectTheDots();

    if (!pConv->IsOption("s", OBConversion::INOPTIONS)
        && !pConv->IsOption("b", OBConversion::INOPTIONS)) {

      mol.PerceiveBondOrders();
    }

    if (hasPartialCharges) {
      mol.SetPartialChargesPerceived();

      // Annotate that partial charges come from Mulliken
      OBPairData* dp = new OBPairData;
      dp->SetAttribute("PartialCharges");
      dp->SetValue("Mulliken");
      dp->SetOrigin(fileformatInput);
      mol.SetData(dp);
    }

    // Found some vibrations
    if (frequencies.size() != 0) {
      OBVibrationData* vd = new OBVibrationData;
      vd->SetData(displacements, frequencies, intensities, raman_intensities);
      vd->SetOrigin(fileformatInput);
      mol.SetData(vd);
    }

    mol.AssignTotalChargeToAtoms(charge);

    mol.SetTotalSpinMultiplicity(mult);

    mol.SetTitle(title);

    return(true);
  }

  ////////////////////////////////////////////////////////////////
  bool GAMESSInputFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv) {
    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if (pmol == NULL)
      return false;

    // Define some references so we can use the old parameter names
    istream& ifs = *pConv->GetInStream();
    OBMol& mol   = *pmol;

    char buffer[BUFF_SIZE];
    string str, str1;
    double x, y, z;
    OBAtom* atom;
    vector<string> vs;
    bool hasPartialCharges = false;
    string efragName; // used to save identifiers of EFRAG sections

    mol.BeginModify();
    while (ifs.getline(buffer, BUFF_SIZE)) {

      if (strstr(buffer, "$DATA") != NULL) {
        // Title
        ifs.getline(buffer, BUFF_SIZE);
        tokenize(vs, buffer);
        mol.SetTitle(buffer);

        // Symetry
        ifs.getline(buffer, BUFF_SIZE);
        ifs.getline(buffer, BUFF_SIZE);

        while (strstr(buffer, "$END") == NULL) {
          tokenize(vs, buffer);
          if (vs.size() == 5) {
            atom = mol.NewAtom();
            // Parse the current one
            atom->SetAtomicNum(atoi(vs[1].c_str()));
            x = atof((char*) vs[2].c_str());
            y = atof((char*) vs[3].c_str());
            z = atof((char*) vs[4].c_str());
            atom->SetVector(x, y, z);
          }

          if (!ifs.getline(buffer, BUFF_SIZE))
            break;
        }
      }

      if (strstr(buffer, "$FMOXYZ") != NULL) {
        while (strstr(buffer, "$END") == NULL) {
          tokenize(vs, buffer);
          if (vs.size() == 5) {
            atom = mol.NewAtom();
            atom->SetAtomicNum(atoi(vs[1].c_str()));
            x = atof((char*) vs[2].c_str());
            y = atof((char*) vs[3].c_str());
            z = atof((char*) vs[4].c_str());
            atom->SetVector(x, y, z);
          }
          if (!ifs.getline(buffer, BUFF_SIZE))
            break;
        }
      }

      if (strstr(buffer, "$EFRAG") != NULL) {
        while (strstr(buffer, "FRAGNAME") == NULL) {
          // Read $EFRAG parameters
          tokenize(vs, buffer, "=");
          if (vs.size() > 1)
            efragName = vs[1];
          if (!ifs.getline(buffer, BUFF_SIZE))
            break;
        }
        while (strstr(buffer," $END") == NULL) {
          tokenize(vs, buffer);
          if (vs.size() == 4) {
            atom = mol.NewAtom();
            int atomicNum;
            if (vs[0].substr(0, 1) == "Z"
                || vs[0].substr(0, 1) == "z") {

              atomicNum = OBElements::GetAtomicNum(vs[0].substr(1, 1).c_str());
            } else {
              atomicNum = OBElements::GetAtomicNum(vs[0].substr(0, 1).c_str());
            }
            atom->SetAtomicNum(atomicNum);
            x = atof((char*) vs[1].c_str());
            y = atof((char*) vs[2].c_str());
            z = atof((char*) vs[3].c_str());
            atom->SetVector(x, y, z);

            // Tag these atoms as part of a specific EFP fragment
            OBPairData* dp = new OBPairData;
            dp->SetAttribute("EFRAG");
            dp->SetValue(efragName);
            dp->SetOrigin(fileformatInput);
            atom->SetData(dp);
          }
          if (!ifs.getline(buffer, BUFF_SIZE))
            break;
        }
      }
    }

    if (!pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.ConnectTheDots();

    if (!pConv->IsOption("s",OBConversion::INOPTIONS)
        && !pConv->IsOption("b",OBConversion::INOPTIONS)) {

      mol.PerceiveBondOrders();
    }

    mol.EndModify();

    if (hasPartialCharges)
      mol.SetPartialChargesPerceived();

    return(true);
  }


  bool GAMESSInputFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv) {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if (pmol == NULL)
      return false;

    // Define some references so we can use the old parameter names
    ostream& ofs = *pConv->GetOutStream();
    OBMol& mol   = *pmol;

    char buffer[BUFF_SIZE];

    const char* keywords       = pConv->IsOption("k", OBConversion::OUTOPTIONS);
    const char* keywordsEnable = pConv->IsOption("k", OBConversion::GENOPTIONS);
    const char* keywordFile    = pConv->IsOption("f", OBConversion::OUTOPTIONS);

    string defaultKeywords = " $CONTRL COORD=CART UNITS=ANGS $END";

    int a        = 0;
    string s     = "";
    bool wrapped = false;

    vector<int> spacePositions;
    std::vector<int>::reverse_iterator rit;
    std::vector<OBGenericData*>::iterator i, j;

    struct local {
      static const bool cmpfn(const char& a, const char& b) {
        return (a == ' ' && b == ' ');
      }
    };

    if (keywords)
      defaultKeywords = keywords;

    if (keywordsEnable) {
      OBSetData* gmsset = (OBSetData*) pmol->GetData("gamess");

      if (gmsset) {
        for (i=gmsset->GetBegin(); i != gmsset->GetEnd(); ++i) {
          OBSetData* cset = (OBSetData*) (*i);

          if (cset) {
            wrapped = false;
            a = 2 + cset->GetAttribute().length();
            ofs << " $" << cset->GetAttribute();
            for (j=cset->GetBegin(); j != cset->GetEnd(); ++j) {
              OBPairData* pd = (OBPairData*) (*j);

              if (pd) {
                if (a + 2 + pd->GetAttribute().length() + pd->GetValue().length() > 72) {
                  // Reached line end
                  s = pd->GetAttribute();
                  s += "=";
                  s += pd->GetValue();

                  // Remove consecutive spaces
                  s.erase(unique(s.begin(), s.end(), local::cmpfn), s.end());

                  while (s.length() > 0) {
                    if (s.find(' ') != string::npos) {
                      // There are spaces in value
                      // Find space positions
                      spacePositions.clear();
                      for (unsigned int n=0; n < s.length(); ++n) {
                        if (s.at(n) == ' ')
                          spacePositions.push_back(n);
                      }

                      // Try to fit it all
                      wrapped = false;
                      if (a + 1 + s.length() <= 72) {
                        a += 1 + s.length();
                        ofs << " " << s;
                        break;
                      }

                      // Try wrapping
                      for (rit=spacePositions.rbegin(); rit != spacePositions.rend(); ++rit) {
                        if (a + 1 + (*rit) <= 72) {
                          ofs << " " << s.substr(0, *rit)
                              << endl << "   ";
                          a = 3;
                          s = s.substr(*rit);
                          wrapped = true;
                          break;
                        }
                      }

                      if (!wrapped) {
                        // String could not be wrapped
                        // Try putting it on the next line
                        a = 4 + spacePositions.at(0);
                        if (a > 72) {
                          // It exceeds line length
                          ofs << endl
                              << "! Unable to fit " << pd->GetAttribute()
                              << " on the line!" << endl;
                          break;
                        }
                        a = 3;
                        ofs << endl << "   ";
                      }
                    } else {
                      // There are no spaces in the string
                      a = 4 + s.length();
                      if (a > 72) {
                        // It exceeds line length
                        ofs << endl
                            << "! Unable to fit " << pd->GetAttribute()
                            << " on the line!" << endl;
                        break;
                      }

                      if (wrapped) {
                        ofs << s;
                      } else {
                        ofs << endl << "    " << s;
                      }
                      s = "";
                    }
                  }
                } else {
                  // Whole key-value pair fits on the line
                  a += 2 + pd->GetAttribute().length() + pd->GetValue().length();
                  ofs << " " << pd->GetAttribute() << "=" << pd->GetValue();
                }
              }
            }
            ofs << " $END" << endl;
          }
        }
      } else {
        ofs << "! Unable to translate keywords!" << endl;
        ofs << "! Defining default control keywords." << endl;
        ofs << defaultKeywords << endl;
      }

    } else if (keywordFile) {
      ifstream kfstream(keywordFile);
      string keyBuffer;
      if (kfstream) {
        while (getline(kfstream, keyBuffer))
          ofs << keyBuffer << endl;
      }

    } else {
        ofs << defaultKeywords << endl;
    }

    ofs << endl << " $DATA" << endl;
    ofs << mol.GetTitle() << endl;
    if (!mol.HasData(OBGenericDataType::SymmetryData))
      ofs << "C1" << endl;
    else {
      // \todo needs to be updated for point group symmetry recognition
      //   particularly for output of the symmetry elements
      //   and any necessary rotation for frame of reference for GAMESS
      ofs << "Put symmetry info here" << endl << endl;
    }

    FOR_ATOMS_OF_MOL(atom, mol) {
      snprintf(buffer, BUFF_SIZE, "%-3s %4d.0    %14.10f  %14.10f  %14.10f ",
               OBElements::GetSymbol(atom->GetAtomicNum()),
               atom->GetAtomicNum(),
               atom->GetX(),
               atom->GetY(),
               atom->GetZ());
      ofs << buffer << endl;
    }

    ofs << " $END" << endl << endl << endl;

    return(true);
  }

} // namespace OpenBabel
