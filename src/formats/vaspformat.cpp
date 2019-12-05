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
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/elements.h>
#include <openbabel/generic.h>


#include <limits>
#include <locale> // For isalpha(int)
#include <functional>
#include <iostream>
#include <algorithm>
#include <cstdlib>

#define EV_TO_KCAL_PER_MOL 23.060538

using namespace std;
namespace OpenBabel {
  class VASPFormat : public OBMoleculeFormat
  {
  protected:
    class compare_sort_items
    {
      std::vector<int> csm;
      bool num_sort;
    public:
      compare_sort_items(const std::vector<int> &_custom_sort_nums, bool _num_sort):
                         csm(_custom_sort_nums), num_sort(_num_sort) {};
      bool operator()(const OBAtom *a, const OBAtom *b)
      {
        int a_num = a->GetAtomicNum();
        int b_num = b->GetAtomicNum();
        int dist = std::distance(std::find(csm.begin(), csm.end(), b_num),
                                 std::find(csm.begin(), csm.end(), a_num));
        
        if ( dist != 0)
          return dist < 0;

        if( (num_sort) && ( a_num - b_num != 0 ) )
          return a_num < b_num;
        
        return false;
      }
    };
  public:

    VASPFormat()
    {
      // This will actually read the CONTCAR file:
      OBConversion::RegisterFormat("CONTCAR",this);
      OBConversion::RegisterFormat("POSCAR",this);
      OBConversion::RegisterFormat("VASP",this);
      OBConversion::RegisterOptionParam("s", this, 0, OBConversion::INOPTIONS);
      OBConversion::RegisterOptionParam("b", this, 0, OBConversion::INOPTIONS);
      OBConversion::RegisterOptionParam("w", this, 0, OBConversion::OUTOPTIONS);
      OBConversion::RegisterOptionParam("z", this, 0, OBConversion::OUTOPTIONS);
      OBConversion::RegisterOptionParam("4", this, 0, OBConversion::OUTOPTIONS);
    }

    virtual const char* Description()
    {
      return
        "VASP format\n"
        "Reads in data from POSCAR and CONTCAR to obtain information from "
        "VASP calculations.\n\n"

        "Due to limitations in Open Babel's file handling, reading in VASP\n"
        "files can be a bit tricky; the client that is using Open Babel must\n"
        "use OBConversion::ReadFile() to begin the conversion. This change is\n"
        "usually trivial. Also, the complete path to the CONTCAR/POSCAR file\n"
        "must be provided, otherwise the other files needed will not be\n"
        "found.\n\n"

        "Both VASP 4.x and 5.x POSCAR formats are supported.\n\n"

	"By default, atoms are written out in the order they are present in the input\n"
	"molecule. To sort by atomic number specify ``-xw``. To specify the sort\n"
	"order, use the ``-xz`` option.\n\n"

        "Read Options e.g. -as\n"
        "  s Output single bonds only\n"
        "  b Disable bonding entirely\n\n"

        "Write Options e.g. -x4\n"
        " w  Sort atoms by atomic number\n"
        " z <list of atoms>  Specify the order to write out atoms\n"
	"       'atom1 atom2 ...': atom1 first, atom2 second, etc. The remaining\n"
	"       atoms are written in the default order or (if ``-xw`` is specified)\n"
	"       in order of atomic number.\n"
        "  4 Write a POSCAR using the VASP 4.x specification.\n"
        "    The default is to use the VASP 5.x specification.\n\n"
        ;

    };

    virtual const char* SpecificationURL(){return "http://cms.mpi.univie.ac.at/vasp/vasp/vasp.html";};

    /* Flags() can return be any of the following combined by |
       or be omitted if none apply
       NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY  DEFAULTFORMAT
       READBINARY  WRITEBINARY  READXML  ZEROATOMSOK */
    virtual unsigned int Flags()
    {
      return READONEONLY;
    };

    virtual int SkipObjects(int n, OBConversion* pConv)
    {
      return 0;
    };

    ////////////////////////////////////////////////////
    /// Declarations for the "API" interface functions. Definitions are below
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

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
    bool selective;    // is selective dynamics used?
    string key, value; // store the info about constraints
    OBPairData *cp;    // in this PairData
    bool hasEnthalpy=false;
    bool hasVibrations=false;
    bool needSymbolsInGeometryFile = false;
    double enthalpy_eV, pv_eV;
    vector<vector <vector3> > Lx;
    vector<double> Frequencies;
    vector<matrix3x3> dipGrad;

    // Get path of CONTCAR/POSCAR:
    //    ifs_path.getline(buffer,BUFF_SIZE);
    //    path = buffer;
    path = pConv->GetInFilename();
    if (path.empty()) return false; // Should be using ReadFile, not Read!
    size_t found;
    found = path.rfind("/");
    path = path.substr(0, found);
    if (found == string::npos) path = "./"; // No "/" in path?

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
    pmol->SetTitle(buffer);
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
    cell->SetSpaceGroup(1);
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
        atomTypes.push_back(OpenBabel::OBElements::GetAtomicNum(vs.at(i).c_str()));
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
          atomTypes.push_back(OpenBabel::OBElements::GetAtomicNum(str.c_str()));
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
    selective = false;
    // Set the variable selective accordingly
    if (buffer[0] == 'S' || buffer[0] == 's') {
      selective = true;
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
      // If we have Cartesian coordinates, we need to apply the scaling factor
      else coords *= scale;
      atom->SetVector(coords);
      //if the selective dynamics info is present then read it into OBPairData
      //this needs to be kept somehow to be able to write out the same as input
      //it's string so it wastes memory :(
      if (selective && vs.size() >= 6) {
        key = "move";
        value  = " "; value += vs[3].c_str();
        value += " "; value += vs[4].c_str();
        value += " "; value += vs[5].c_str();
        cp = new OBPairData;
        cp->SetAttribute(key);
        cp->SetValue(value);
        cp->SetOrigin(fileformatInput);
        atom->SetData(cp);
      }

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

    // Vibration intensities
    vector3 prevDm;
    vector<vector3> prevXyz;
    vector3 currDm;
    vector<vector3> currXyz;

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
          hasVibrations = true;
          double x, y, z;
          ifs_out.getline(buffer,BUFF_SIZE);  // dash line
          ifs_out.getline(buffer,BUFF_SIZE);  // blank line
          ifs_out.getline(buffer,BUFF_SIZE);  // blank line
          ifs_out.getline(buffer,BUFF_SIZE);  // first frequency line
          while (!strstr(buffer, "Eigenvectors")) {
            vector<vector3> vib;
            tokenize(vs, buffer);
            int freqnum = atoi(vs[0].c_str());
            if (vs[1].size() == 1 and vs[1].compare("f") == 0) {
              // Real frequency
              Frequencies.push_back(atof(vs[7].c_str()));
            } else if (strstr(vs[1].c_str(), "f/i=")) {
              // Imaginary frequency
              Frequencies.push_back(-atof(vs[6].c_str()));
            } else {
              // No more frequencies
              break;
            }
            ifs_out.getline(buffer,BUFF_SIZE);  // header line
            ifs_out.getline(buffer,BUFF_SIZE);  // first displacement line
            tokenize(vs, buffer);
            // normal modes
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
        }

        if (strstr(buffer, "dipolmoment")) {
          tokenize(vs, buffer);
          x = atof(vs[1].c_str());
          y = atof(vs[2].c_str());
          z = atof(vs[3].c_str());
          currDm.Set(x, y, z);
        }
        if (strstr(buffer, "TOTAL-FORCE")) {
          currXyz.clear();
          ifs_out.getline(buffer, BUFF_SIZE);  // header line
          ifs_out.getline(buffer, BUFF_SIZE);
          tokenize(vs, buffer);
          while (vs.size() == 6) {
            x = atof(vs[0].c_str());
            y = atof(vs[1].c_str());
            z = atof(vs[2].c_str());
            currXyz.push_back(vector3(x, y, z));
            ifs_out.getline(buffer, BUFF_SIZE);  // next line
            tokenize(vs, buffer);
          }
        }
        if (strstr(buffer, "BORN EFFECTIVE CHARGES")) {
          // IBRION = 7; IBRION = 8
          dipGrad.clear();
          ifs_out.getline(buffer, BUFF_SIZE);  // header line
          ifs_out.getline(buffer, BUFF_SIZE);  // `ion    #`
          tokenize(vs, buffer);
          while (vs.size() == 2) {
            matrix3x3 dmudq;
            for (int row = 0; row < 3; ++row) {
              ifs_out.getline(buffer, BUFF_SIZE);
              tokenize(vs, buffer);
              x = atof(vs[1].c_str());
              y = atof(vs[2].c_str());
              z = atof(vs[3].c_str());
              dmudq.SetRow(row, vector3(x, y, z));
            }
            dipGrad.push_back(dmudq);
            ifs_out.getline(buffer, BUFF_SIZE);  // next line
            tokenize(vs, buffer);
          }
        } else if (strstr(buffer, "free  energy")) {
          // IBRION = 5
          // reached the end of an iteration, use the values
          if (dipGrad.empty()) {
            // first iteration: nondisplaced ions
            dipGrad.resize(pmol->NumAtoms());
          } else if (prevXyz.empty()) {
            // even iteration: store values
            prevXyz = currXyz;
            prevDm = currDm;
          } else {
            // odd iteration: compute dipGrad = dmu / dxyz for moved ion
            for (size_t natom = 0; natom < pmol->NumAtoms(); ++natom) {
              const vector3 dxyz = currXyz[natom] - prevXyz[natom];
              vector3::const_iterator iter = std::find_if(dxyz.begin(), dxyz.end(),
                      std::bind2nd(std::not_equal_to<double>(), 0.0));
              if (iter != dxyz.end()) dipGrad[natom].SetRow(iter - dxyz.begin(),
                                                            (currDm - prevDm) / *iter);
            }
            prevXyz.clear();
          }
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

    // Set vibrations
    if (hasVibrations) {
      // compute dDip/dQ
      vector<double> Intensities;
      for (vector<vector<vector3> >::const_iterator
           lxIter = Lx.begin(); lxIter != Lx.end(); ++lxIter) {
        vector3 intensity;
        for (size_t natom = 0; natom < dipGrad.size(); ++natom) {
          intensity += dipGrad[natom].transpose() * lxIter->at(natom)
              / sqrt(pmol->GetAtomById(natom)->GetAtomicMass());
        }
        Intensities.push_back(dot(intensity, intensity));
      }
      const double max = *max_element(Intensities.begin(), Intensities.end());
      if (max != 0.0) {
        // Normalize
        std::transform(Intensities.begin(), Intensities.end(), Intensities.begin(),
                       std::bind2nd(std::divides<double>(), max / 100.0));
      } else {
        Intensities.clear();
      }
      OBVibrationData* vd = new OBVibrationData;
      vd->SetData(Lx, Frequencies, Intensities);
      pmol->SetData(vd);
    }

    pmol->EndModify();

    const char *noBonding  = pConv->IsOption("b", OBConversion::INOPTIONS);
    const char *singleOnly = pConv->IsOption("s", OBConversion::INOPTIONS);

    if (noBonding == NULL) {
      pmol->ConnectTheDots();
      if (singleOnly == NULL) {
        pmol->PerceiveBondOrders();
      }
    }

    return true;
  }

  bool VASPFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    //No surprises in this routine, cartesian coordinates are written out
    //and if at least a single atom has information about constraints,
    //then selective dynamics is used and the info is written out.
    //The atoms are ordered according to their atomic number so that the
    //output looks nice, this can be reversed by using command line flag "-xw".
    //
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if (pmol == NULL) {
      return false;
    }

    ostream& ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    char buffer[BUFF_SIZE];
    OBUnitCell *uc = NULL;
    vector<vector3> cell;

    const char * sortAtomsNum = pConv->IsOption("w", OBConversion::OUTOPTIONS);
    const char * sortAtomsCustom = pConv->IsOption("z", OBConversion::OUTOPTIONS);

    // Create a list of ids. These may be sorted by atomic number depending
    // on the value of keepOrder.
    std::vector<OBAtom *> atoms_sorted;
    atoms_sorted.reserve(mol.NumAtoms());

    FOR_ATOMS_OF_MOL(atom, mol) {
      atoms_sorted.push_back(&(*atom));
    }

    std::vector<int> custom_sort_nums;
    
    if (sortAtomsCustom != NULL)
    {
      vector<string> vs;
      tokenize(vs, sortAtomsCustom);
      for(size_t i = 0; i < vs.size(); ++i)
        custom_sort_nums.push_back(OBElements::GetAtomicNum(vs[i].c_str()));
    }

    compare_sort_items csi(custom_sort_nums, sortAtomsNum != NULL);
    std::stable_sort(atoms_sorted.begin(), atoms_sorted.end(), csi);

    // Use the atomicNums vector to determine the composition line.
    // atomicNumsCondensed and atomCounts contain the same data as atomicNums:
    // if:
    //   atoms_sorted[i]->GetAtomicNum() = [ 3 3 3 2 2 8 2 6 6 ]
    // then:
    //   atomicNums =  [(3 3) (2 2) (8 1) (2 1) (6 2)] 
    
    std::vector<std::pair<int, int> > atomicNums;    
    
    int prev_anum = -20; //not a periodic table number
    for(int i = 0; i < atoms_sorted.size(); i++)
    {
      const int anum = atoms_sorted[i]->GetAtomicNum();
      
      if( prev_anum != anum )
      {
        std::pair<int, int> x(anum, 1);
        atomicNums.push_back(x);
      }
      else
      {    
        if(atomicNums.size() > 0)  
          atomicNums.rbegin()->second++;
      }  
      
      prev_anum = anum;
    }

    // write title
    ofs << mol.GetTitle() << endl;
    // write the multiplication factor, set this to one
    // and write the cell using the 3x3 cell matrix
    ofs << "1.000 " << endl;

    if (!mol.HasData(OBGenericDataType::UnitCell)) {
      // the unit cell has not been defined. Leave as all zeros so the user
      // can fill it in themselves
      for (int ii = 0; ii < 3; ii++) {
        snprintf(buffer, BUFF_SIZE, "0.0  0.0  0.0");
        ofs << buffer << endl;
      }
    }
    else
    {
      // there is a unit cell, write it out
      uc = static_cast<OBUnitCell*>(mol.GetData(OBGenericDataType::UnitCell));
      cell = uc->GetCellVectors();
      for (vector<vector3>::const_iterator i = cell.begin();
           i != cell.end(); ++i) {
        snprintf(buffer, BUFF_SIZE, "%20.15f%20.15f%20.15f",
                 i->x(), i->y(), i->z());
        ofs << buffer << endl;
      }
    }

    // go through the atoms first to write out the element names if using
    // VASP 5 format
    const char *vasp4Format = pConv->IsOption("4", OBConversion::OUTOPTIONS);
    if (!vasp4Format) {
      for (vector< std::pair<int, int> >::const_iterator
           it = atomicNums.begin(),
           it_end = atomicNums.end(); it != it_end; ++it) {
        snprintf(buffer, BUFF_SIZE, "%-3s ", OBElements::GetSymbol(it->first));
        ofs << buffer ;
      }
      ofs << endl;
    }

    // then do the same to write out the number of ions of each element
    for (vector< std::pair<int, int> >::const_iterator
           it = atomicNums.begin(),
           it_end = atomicNums.end(); it != it_end; ++it) {
      snprintf(buffer, BUFF_SIZE, "%-3u ", it->second);
      ofs << buffer ;
    }
    ofs << endl;

    // assume that there are no constraints on the atoms
    bool selective = false;
    // and test if any of the atoms has constraints
    FOR_ATOMS_OF_MOL(atom, mol) {
      if (atom->HasData("move")){
        selective = true;
        break;
      }
    }
    if (selective) {
      ofs << "SelectiveDyn" << endl;
    }

    // print the atomic coordinates in \AA
    ofs << "Cartesian" << endl;

    for (std::vector<OBAtom *>::const_iterator it = atoms_sorted.begin();
         it != atoms_sorted.end(); ++it) 
    {
      // Print coordinates
      snprintf(buffer,BUFF_SIZE, "%26.19f %26.19f %26.19f",
               (*it)->GetX(), (*it)->GetY(), (*it)->GetZ());
      ofs << buffer;

      // if at least one atom has info about constraints
      if (selective) {
        // if this guy has, write it out
        if ((*it)->HasData("move")) {
          OBPairData *cp = (OBPairData*)(*it)->GetData("move");
          // seemingly ridiculous number of digits is written out
          // but sometimes you just don't want to change them
          ofs << " " << cp->GetValue().c_str();
        }
        else {
          // the atom has been created and the info has not been copied
          ofs << "  T T T";
        }
      }
      ofs << endl;
    }

    return true;
  }

} //namespace OpenBabel

