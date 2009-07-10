/**********************************************************************
Copyright (C) 2004 by Chris Morley for template
Copyright (C) 2004 by David C. Lonie for VASP
 
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
        "VASP format\n \
	Reads in data from POTCAR and CONTCAR to obtain information from VASP calculations.\n \
	Due to limitations in OB's file handling, reading in VASP files can be a bit tricky:\n \
	\tThe client that is using OpenBabel must use OBConversion::ReadFile() to begin the conversion.\n \
	\tThis change is usually trivial, ask the package mantainer to look into it. Also, the complete\n \
	\tpath to the CONTCAR file must be provided, otherwise the other files need won't be found.\n \
		";
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
    vector<int> numAtoms, atomTypes;

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
    string contcar_filename = pConv->GetInFilename(); // POSCAR _OR_ CONTCAR 
    ifstream ifs_pot (potcar_filename.c_str());
    ifstream ifs_out (outcar_filename.c_str());
    ifstream ifs_cont (contcar_filename.c_str());
     if (!ifs_pot || !ifs_cont)
       return false; // Missing some files

     pmol->BeginModify();

     // Read in optional information from outcar
     if (ifs_out) {
       while (ifs_out.getline(buffer,BUFF_SIZE)) {
         // Enthalphy
         if (strstr(buffer, "enthalpy is")) {
           tokenize(vs, buffer);
           OBPairData *enthalpy = new OBPairData();
           OBPairData *enthalpy_pv = new OBPairData();
           OBPairData *enthalpy_eV = new OBPairData();
           OBPairData *enthalpy_pv_eV = new OBPairData();
           enthalpy->SetAttribute("Enthalpy (kcal/mol)");
           enthalpy_pv->SetAttribute("Enthalpy PV term (kcal/mol)");
           enthalpy_eV->SetAttribute("Enthalpy (eV)");
           enthalpy_pv_eV->SetAttribute("Enthalpy PV term (eV)");
           float en = atof(vs[4].c_str()) * EV_TO_KCAL_PER_MOL;
           float pv = atof(vs[8].c_str()) * EV_TO_KCAL_PER_MOL;
           float en_eV = atof(vs[4].c_str());
           float pv_eV = atof(vs[8].c_str());
           snprintf(tag, BUFF_SIZE, "%f", en);
           enthalpy->SetValue(tag);
           snprintf(tag, BUFF_SIZE, "%f", pv);
           enthalpy_pv->SetValue(tag);
           snprintf(tag, BUFF_SIZE, "%f", en_eV);
           enthalpy_eV->SetValue(tag);
           snprintf(tag, BUFF_SIZE, "%f", pv_eV);
           enthalpy_pv_eV->SetValue(tag);
           pmol->SetData(enthalpy);
           pmol->SetData(enthalpy_pv);
           pmol->SetData(enthalpy_eV);
           pmol->SetData(enthalpy_pv_eV);
         }

         // Free energy
         if (strstr(buffer, "free  energy")) {
           tokenize(vs, buffer);
           pmol->SetEnergy(atof(vs[4].c_str()) * EV_TO_KCAL_PER_MOL);
         }
       }
     }

     ifs_out.close();

     // Get atom types from POTCAR
     while (ifs_pot.getline(buffer,BUFF_SIZE)) {
       if (strstr(buffer,"VRHFIN")) {
         str = buffer;
         size_t start = str.find("=") + 1;
         size_t end = str.find(":");
         str = str.substr(start, end - start);
         // Clean up whitespace:
         for (unsigned int i = 0; i < str.size(); i++)
           if (str.at(i) == ' ')
             str.erase(i,1);
         atomTypes.push_back(OpenBabel::etab.GetAtomicNum(str.c_str()));
       };
     };

     ifs_pot.close();

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

     ifs_cont.getline(buffer,BUFF_SIZE); // Number of atoms per atom type. This is in the same order as atomTypes.
     tokenize(vs, buffer);
     for (unsigned int i = 0; i < atomTypes.size(); i++) {
       numAtoms.push_back(atoi(vs.at(i).c_str()));
       totalAtoms += atoi(vs.at(i).c_str());
     }

     ifs_cont.getline(buffer,BUFF_SIZE); // Cartesian or fractional?
     if (buffer[0] == 'S' || buffer[0] == 's') // Skip selective dynamics line if present.
       ifs_cont.getline(buffer,BUFF_SIZE);
     if ( buffer[0] == 'C' || buffer[0] == 'c' ||
          buffer[0] == 'K' || buffer[0] == 'k' )
       cartesian = true;
     else
       cartesian = false;

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
         coords = cell->GetOrthoMatrix() * coords;
       atom->SetVector(coords);
       atomCount++;
     };

     // There is some trailing garbage, but AFAIK it's not useful for anything.
     ifs_cont.close();

     pmol->EndModify();

     return true;
   }

 } //namespace OpenBabel

