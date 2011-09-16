/**********************************************************************
Copyright (C) 2004 by Chris Morley for template
Copyright (C) 2009 by David C. Lonie for VASP

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

#include <locale> // For isalpha(int)

#define EV_TO_KCAL_PER_MOL 23.060538

using namespace std;
namespace OpenBabel {
  class VASPFormat : public OBMoleculeFormat
  {
  public:

    VASPFormat()
    {
      // This will actually read the CONTCAR file:
      //      OBConversion::RegisterFormat("vasp",this);
      OBConversion::RegisterFormat("CONTCAR",this);
      OBConversion::RegisterFormat("POSCAR",this);
    }

    virtual const char* Description()
    {
      return
        "VASP format\n"
        "Reads in data from POSCAR and CONTCAR to obtain information from VASP calculations.\n\n"

"Due to limitations in Open Babel's file handling, reading in VASP files can\n"
"be a bit tricky; the client that is using Open Babel must use\n"
"OBConversion::ReadFile() to begin the conversion. This change is usually\n"
"trivial. Also, the complete path to the CONTCAR file must be provided,\n"
"otherwise the other files needed will not be found.\n";
    };

    virtual const char* SpecificationURL(){return "http://cms.mpi.univie.ac.at/vasp/vasp/vasp.html";};

    /* Flags() can return be any of the following combined by |
       or be omitted if none apply
       NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY  DEFAULTFORMAT
       READBINARY  WRITEBINARY  READXML  ZEROATOMSOK */
    virtual unsigned int Flags()
    {
      return READONEONLY | NOTWRITABLE;
    };

    virtual int SkipObjects(int n, OBConversion* pConv)
    {
      return 0;
    };

    ////////////////////////////////////////////////////
    /// Declarations for the "API" interface functions. Definitions are below
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    //    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

  private:
    /* Add declarations for any local function or member variables used.
       Generally only a single instance of a format class is used. Keep this in
       mind if you employ member variables. */
  };
  ////////////////////////////////////////////////////

  //Make an instance of the format class
  VASPFormat theVASPFormat;

  /////////////////////////////////////////////////////////////////

  bool VASPFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    // Move stream to EOF, some apps check ifs position to check for multimolecule files.
    // VASP does not support this, and this parser makes its own streams.
    istream &ifs = *pConv->GetInStream();
    ifs.seekg (0, ios::end);

    char buffer[BUFF_SIZE], tag[BUFF_SIZE];
    double x,y,z, scale;
    unsigned int totalAtoms=0, atomIndex=0, atomCount=0;
    OBAtom *atom;
    bool cartesian;
    string str, path;
    vector<string> vs;
    vector<unsigned int> numAtoms, atomTypes;
    bool hasEnthalpy=false;
    bool needSymbolsInGeometryFile = false;
    double enthalpy_eV, pv_eV;
    vector<vector <vector3> > Lx;
    vector<double> Frequencies, Intensities;


    // Get path of CONTCAR/POSCAR:
    //    ifs_path.getline(buffer,BUFF_SIZE);
    //    path = buffer;
    path = pConv->GetInFilename();
    if (path.empty()) return false; // Should be using ReadFile, not Read!
    size_t found;
    found = path.rfind("/");
    if (found == string::npos) return false; // No "/" in path?
    path = path.substr(0,found);

    // Open files
    string potcar_filename = path + "/POTCAR";
    string outcar_filename = path + "/OUTCAR";
    string doscar_filename = path + "/DOSCAR";
    string contcar_filename = pConv->GetInFilename(); // POSCAR _OR_ CONTCAR
    ifstream ifs_pot (potcar_filename.c_str());
    ifstream ifs_out (outcar_filename.c_str());
    ifstream ifs_dos (doscar_filename.c_str());
    ifstream ifs_cont (contcar_filename.c_str());
    if (!ifs_pot) {
      needSymbolsInGeometryFile = true;
    }
    if (!ifs_cont) {
      return false; // No geometry file?
    }

    pmol->BeginModify();

    // Start working on CONTCAR:
    ifs_cont.getline(buffer,BUFF_SIZE); // Comment line
    ifs_cont.getline(buffer,BUFF_SIZE); // Scale
    scale = atof(buffer);

    ifs_cont.getline(buffer,BUFF_SIZE); // X_Vec vector
    tokenize(vs, buffer);
    x = atof(vs.at(0).c_str()) * scale;
    y = atof(vs.at(1).c_str()) * scale;
    z = atof(vs.at(2).c_str()) * scale;
    vector3 x_vec (x,y,z);

    ifs_cont.getline(buffer,BUFF_SIZE); // Y_Vec vector
    tokenize(vs, buffer);
    x = atof(vs.at(0).c_str()) * scale;
    y = atof(vs.at(1).c_str()) * scale;
    z = atof(vs.at(2).c_str()) * scale;
    vector3 y_vec (x,y,z);

    ifs_cont.getline(buffer,BUFF_SIZE); // Z_Vec vector
    tokenize(vs, buffer);
    x = atof(vs.at(0).c_str()) * scale;
    y = atof(vs.at(1).c_str()) * scale;
    z = atof(vs.at(2).c_str()) * scale;
    vector3 z_vec (x,y,z);

    // Build unit cell
    OBUnitCell *cell = new OBUnitCell;
    cell->SetData(x_vec, y_vec, z_vec);
    pmol->SetData(cell);

    // Next comes either a list of numbers that represent the stoichiometry of
    // the cell. The numbers are the atom counts for each type, in the order
    // listed in the POTCAR file. Since VASP 5.2, a line with a list of atomic
    // element symbols may precede the atom counts. This will be used if the
    // POTCAR file is not present. If both are present, the data in the POSCAR
    // or CONTCAR file will be used.
    ifs_cont.getline(buffer,BUFF_SIZE);
    tokenize(vs, buffer);
    bool symbolsInGeometryFile = false;
    if (vs.size() != 0) {
      if (vs.at(0).size() != 0) {
        if (isalpha(static_cast<int>(vs.at(0).at(0))) != 0) {
          symbolsInGeometryFile = true;
        }
      }
    }

    // If no element data is present anywhere
    if (needSymbolsInGeometryFile && !symbolsInGeometryFile &&
        // and there are atoms specified
        vs.size() != 0) {
      // Abort read
      pmol->EndModify();
      return false;
    }

    if (symbolsInGeometryFile) {
      atomTypes.clear();
      for (size_t i = 0; i < vs.size(); ++i) {
        atomTypes.push_back(OpenBabel::etab.GetAtomicNum(vs.at(i).c_str()));
      }
      // Fetch next line to get stoichiometry
      ifs_cont.getline(buffer,BUFF_SIZE);
      tokenize(vs, buffer);
    }
    else if (ifs_pot) {
      // Get atom types from POTCAR
      while (ifs_pot.getline(buffer,BUFF_SIZE)) {
        if (strstr(buffer,"VRHFIN")) {
          str = buffer;
          size_t start = str.find("=") + 1;
          size_t end = str.find(":");
          str = str.substr(start, end - start);
          // Clean up whitespace:
          for (unsigned int i = 0; i < str.size(); i++)
            if (str.at(i) == ' ') {
              str.erase(i,1);
              --i;
            }
          atomTypes.push_back(OpenBabel::etab.GetAtomicNum(str.c_str()));
        }
      }
      ifs_pot.close();
    }

    // Extract and sum the atom counts. The sum is used to parse the atomic
    // coordinates
    totalAtoms = 0;
    for (unsigned int i = 0; i < vs.size(); i++) {
      int currentCount = atoi(vs.at(i).c_str());
      numAtoms.push_back(currentCount);
      totalAtoms += currentCount;
    }

    // Do the number of atom types match the number of atom counts?
    if (numAtoms.size() != atomTypes.size()) {
      // If not, abort read
      pmol->EndModify();
      return false;
    }

    // Cartesian or fractional?
    ifs_cont.getline(buffer,BUFF_SIZE);
    // Skip selective dynamics line if present.
    if (buffer[0] == 'S' || buffer[0] == 's') {
      ifs_cont.getline(buffer,BUFF_SIZE);
    }
    // [C|c|K|k] indicates cartesian coordinates, anything else (including
    // an empty string, buffer[0] == 0) indicates fractional coordinates
    if ( buffer[0] == 'C' || buffer[0] == 'c' ||
         buffer[0] == 'K' || buffer[0] == 'k' ) {
      cartesian = true;
    }
    else {
      cartesian = false;
    }

    atomCount = 0;
    for (unsigned int i = 0; i < totalAtoms; i++) {
      // Things get a little nasty here. VASP just prints all the coordinates with no
      // identifying information one after another here. So in the above sections we've
      // parsed out which atom types and how many of each are present in atomTypes and
      // numAtoms, respectively. The counters atomIndex and atomCount have the following
      // roles: atomIndex keeps track of where we are in atomTypes so that we know the
      // atomic species we're inserting. atomCount tracks how many of the current
      // atomTypes.at(atomIndex) species have been inserted, so that when we reach
      // (atomCount >= numAtoms.at(atomIndex) ) we should stop. Phew.
      ifs_cont.getline(buffer,BUFF_SIZE); // atom location

      // Let's first figure out the atomic number we are dealing with:
      while (atomCount >= numAtoms.at(atomIndex)) {
        atomIndex++;
        atomCount = 0;
      }

      // If we made it past that check, we have atomic number = atomTypes.at(atomIndex)
      // Parse the buffer now.
      tokenize(vs, buffer);
      atom = pmol->NewAtom();
      atom->SetAtomicNum(atomTypes.at(atomIndex));
      x = atof((char*)vs[0].c_str());
      y = atof((char*)vs[1].c_str());
      z = atof((char*)vs[2].c_str());
      vector3 coords (x,y,z);
      if (!cartesian)
        coords = cell->FractionalToCartesian( coords );
      atom->SetVector(coords);
      atomCount++;
    };

    // There is some trailing garbage, but AFAIK it's not useful for anything.
    ifs_cont.close();

    // Read density of states info from DOSCAR, if available
    if (ifs_dos) {
      // Create DOS object
      OBDOSData *dos = new OBDOSData();

      // skip header
      ifs_dos.getline(buffer,BUFF_SIZE); // Junk
      ifs_dos.getline(buffer,BUFF_SIZE); // Junk
      ifs_dos.getline(buffer,BUFF_SIZE); // Junk
      ifs_dos.getline(buffer,BUFF_SIZE); // Junk
      ifs_dos.getline(buffer,BUFF_SIZE); // Junk

      // Get fermi level
      double fermi;
      if (ifs_dos.getline(buffer,BUFF_SIZE)) { // startE endE res fermi ???
        tokenize(vs, buffer);
        fermi = atof(vs[3].c_str());
      }

      // Start pulling out energies and densities
      std::vector<double> energies;
      std::vector<double> densities;
      std::vector<double> integration;
      while (ifs_dos.getline(buffer,BUFF_SIZE)) {
        tokenize(vs, buffer);
        energies.push_back(atof(vs[0].c_str()));
        densities.push_back(atof(vs[1].c_str()));
        integration.push_back(atof(vs[2].c_str()));
      }

      if (energies.size() != 0) {
        dos->SetData(fermi, energies, densities, integration);
        pmol->SetData(dos);
      }
    }

    ifs_dos.close();

    // Read in optional information from outcar
    if (ifs_out) {
      while (ifs_out.getline(buffer,BUFF_SIZE)) {
        // Enthalphy
        if (strstr(buffer, "enthalpy is")) {
          hasEnthalpy = true;
          tokenize(vs, buffer);
          enthalpy_eV = atof(vs[4].c_str());
          pv_eV = atof(vs[8].c_str());
        }

        // Free energy
        if (strstr(buffer, "free  energy")) {
          tokenize(vs, buffer);
          pmol->SetEnergy(atof(vs[4].c_str()) * EV_TO_KCAL_PER_MOL);
        }
        // Frequencies
        if (strstr(buffer, "Eigenvectors") && Frequencies.size() == 0) {
          double x, y, z;
          ifs_out.getline(buffer,BUFF_SIZE);  // dash line
          ifs_out.getline(buffer,BUFF_SIZE);  // blank line
          ifs_out.getline(buffer,BUFF_SIZE);  // blank line
          ifs_out.getline(buffer,BUFF_SIZE);  // first frequency line
          while (!strstr(buffer, "Eigenvectors")) {
            vector<vector3> vib;
            tokenize(vs, buffer);
            if (vs.size() < 2) {
              // No more frequencies
              break;
            }
            int freqnum = atoi(vs[0].c_str());
            if (strstr(vs[1].c_str(), "f/i=")) {
              // Imaginary frequency
              Frequencies.push_back(-atof(vs[6].c_str()));
            } else {
              Frequencies.push_back(atof(vs[7].c_str()));
            }
            // TODO: Intensities not parsed yet
            Intensities.push_back(0.0);
            ifs_out.getline(buffer,BUFF_SIZE);  // header line
            ifs_out.getline(buffer,BUFF_SIZE);  // first displacement line
            tokenize(vs, buffer);
            while (vs.size() == 6) {
              x = atof(vs[3].c_str());
              y = atof(vs[4].c_str());
              z = atof(vs[5].c_str());
              vib.push_back(vector3(x, y, z));
              ifs_out.getline(buffer,BUFF_SIZE);  // next displacement line
              tokenize(vs, buffer);
            }
            Lx.push_back(vib);
            ifs_out.getline(buffer,BUFF_SIZE);  // next frequency line
          }
          OBVibrationData* vd = new OBVibrationData;
          vd->SetData(Lx, Frequencies, Intensities);
          pmol->SetData(vd);
        }
      }
    }

    ifs_out.close();

    // Set enthalpy
    if (hasEnthalpy) {
      OBPairData *enthalpyPD = new OBPairData();
      OBPairData *enthalpyPD_pv = new OBPairData();
      OBPairData *enthalpyPD_eV = new OBPairData();
      OBPairData *enthalpyPD_pv_eV = new OBPairData();
      enthalpyPD->SetAttribute("Enthalpy (kcal/mol)");
      enthalpyPD_pv->SetAttribute("Enthalpy PV term (kcal/mol)");
      enthalpyPD_eV->SetAttribute("Enthalpy (eV)");
      enthalpyPD_pv_eV->SetAttribute("Enthalpy PV term (eV)");
      double en_kcal_per_mole = enthalpy_eV * EV_TO_KCAL_PER_MOL;
      double pv_kcal_per_mole = pv_eV * EV_TO_KCAL_PER_MOL;
      snprintf(tag, BUFF_SIZE, "%f", en_kcal_per_mole);
      enthalpyPD->SetValue(tag);
      snprintf(tag, BUFF_SIZE, "%f", pv_kcal_per_mole);
      enthalpyPD_pv->SetValue(tag);
      snprintf(tag, BUFF_SIZE, "%f", enthalpy_eV);
      enthalpyPD_eV->SetValue(tag);
      snprintf(tag, BUFF_SIZE, "%f", pv_eV);
      enthalpyPD_pv_eV->SetValue(tag);
      pmol->SetData(enthalpyPD);
      pmol->SetData(enthalpyPD_pv);
      pmol->SetData(enthalpyPD_eV);
      pmol->SetData(enthalpyPD_pv_eV);
    }

    pmol->EndModify();

    return true;
  }

} //namespace OpenBabel
