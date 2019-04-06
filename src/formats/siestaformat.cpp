/**********************************************************************
Copyright (C) 2004 by Chris Morley for template
Copyright (C) 2015 by Patrick S. Avery for SIESTA

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
#include <cstdlib>


#define EV_TO_KCAL_PER_MOL 23.060538

using namespace std;
namespace OpenBabel {
  class SIESTAFormat : public OBMoleculeFormat
  {
  public:

    SIESTAFormat()
    {
      OBConversion::RegisterFormat("siesta",this);
    }

    virtual const char* Description()
    {
      return
        "SIESTA format\n"
        "The format used by SIESTA (Spanish Initiative for Electronic Simulations with Thousands of Atoms).\n\n";
    };

    virtual const char* SpecificationURL(){return "http://departments.icmab.es/leem/siesta/";};

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
  SIESTAFormat theSIESTAFormat;

  /////////////////////////////////////////////////////////////////

  bool SIESTAFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();

    std::stringstream errorMsg;

    bool struct_out_found = false;
    bool coordsAreFractional = false;
    bool coordsAreAngstrom = false;
    bool cellDataWasSet = false;
    bool EOFReached = false;
    char buffer[BUFF_SIZE];
    double x,y,z,a,b,c,alpha,beta,gamma;
    double latticeConstant = 1.0;
    size_t size;
    vector<string> vs;
    // The int is the species label, and the string is the atomic symbol
    map<int, string> atomTypeLabels;
    int atomicNum, numAtoms, numSpecies;
    matrix3x3 cellMatrix;
    OBUnitCell *cell = new OBUnitCell();
    pmol->BeginModify();

    // We're gonna try to find "<name>.STRUCT_OUT"
    // If we find it and read it, we will skip all the reading steps
    // for atomic positions and cell coordinates.
    string filePath = pConv->GetInFilename();
    if (filePath.empty()) {
      delete cell;
      errorMsg << "Invalid path specified for siesta output file.\n";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      return false;
    }

    size_t found;
    found = filePath.rfind("/");
    string path = filePath.substr(0, found) + "/";
    if (found == string::npos) path = "./";
    string fileName = filePath.substr(found + 1, filePath.size());

    // We want to chop off the last four characters: ".out"
    found = fileName.find(".out");
    if (found != std::string::npos && found == fileName.size() - 4)
      fileName.erase(fileName.size() - 4, fileName.size() - 1);

    string struct_name;
    ifstream ifs_struct_out;
    vector<string> directories;
    // Append more directories to this list in the future if there are
    // other common directories to find the .STRUCT_OUT file in.
    directories.push_back(".");
    directories.push_back("work");
    for (size_t i = 0; i < directories.size(); i++) {
      struct_name = path + directories.at(i) + "/" + fileName + ".STRUCT_OUT";
      ifs_struct_out.open(struct_name.c_str());
      if (ifs_struct_out) break;
      else ifs_struct_out.close();
    }
    if (ifs_struct_out) struct_out_found = true;

    // If these failed, we'll just look for the coordinates in the .out file
    // Send a message to the user if the .STRUCT_OUT was not found
    if (!struct_out_found) {
      errorMsg << "Could not find " << fileName
               << ".STRUCT_OUT\nAttempting to read " << "coordinates from "
               << fileName << ".out " << "instead.\n";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      errorMsg.clear();
    }

    // Read the .STRUCT_OUT file if it was found
    else if (struct_out_found) {
      // First three lines are cell vectors
      // They are in Angstroms by default
      for (size_t i = 0; i < 3; i++) {
        ifs_struct_out.getline(buffer, BUFF_SIZE);
        tokenize(vs, buffer);
        for (size_t j = 0; j < 3; j++) {
          cellMatrix(i,j) = atof(vs.at(j).c_str());
        }
      }

      // Build unit cell
      cell->SetData(cellMatrix);
      cellDataWasSet = true;

      // .STRUCT_OUT gives fractional coordinates for atom positions
      coordsAreFractional = true;

      // Get number of atoms
      ifs_struct_out.getline(buffer, BUFF_SIZE);
      tokenize(vs, buffer);
      numAtoms = atoi(vs.at(0).c_str());

      // Clear old atoms from pmol in case they were set for some reason
      vector<OBAtom*> toDelete;
      FOR_ATOMS_OF_MOL(a, *pmol)
        toDelete.push_back(&*a);
      for (size_t i = 0; i < toDelete.size(); i++)
        pmol->DeleteAtom(toDelete.at(i));

      // Store atom info
      int atomsIterated = 0;
      while (atomsIterated < numAtoms &&
             ifs_struct_out.getline(buffer, BUFF_SIZE)) {

        tokenize(vs, buffer);

        // save atomic number
        atomicNum = atoi(vs.at(1).c_str());

        x = atof(vs.at(2).c_str());
        y = atof(vs.at(3).c_str());
        z = atof(vs.at(4).c_str());

        // Add atom
        OBAtom *atom = pmol->NewAtom();
        atom->SetAtomicNum(atomicNum);
        vector3 coords (x,y,z);
        atom->SetVector(coords);
        atomsIterated++;
      }
      if (atomsIterated != numAtoms) {
        delete cell;
        errorMsg << "Error reading the .STRUCT_OUT file. Make sure it was "
                 << "saved correctly\n";
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
        return false;
      }
    } // Done reading .STRUCT_OUT !

    // Will attempt to read coordinates from the .out file if the .STRUCT_OUT
    // file was not found.
    while (ifs.getline(buffer, BUFF_SIZE)) {
/*
      // These are currently unused, but may be used in the future...
      if(strstr(buffer, "NumberOfSpecies")) {
        tokenize(vs, buffer);
        numSpecies = atoi(vs.at(1).c_str());
      }

      if (strstr(buffer, "NumberOfAtoms")) {
        tokenize(vs, buffer);
        numAtoms = atoi(vs.at(1).c_str());
      }
*/
      // ChemicalSpeciesLabels are stored in atomTypeLabels
      if (strstr(buffer, "%block ChemicalSpeciesLabel")) {
        ifs.getline(buffer, BUFF_SIZE);
        tokenize(vs, buffer);
        while (!strstr(buffer, "%endblock ChemicalSpeciesLabel")) {
          atomTypeLabels[atoi(vs.at(0).c_str())] = vs.at(2).c_str();
          ifs.getline(buffer, BUFF_SIZE);
          tokenize(vs, buffer);
        }
      }

      // Input coordinates are fractional!
      if (strstr(buffer, "AtomicCoordinatesFormat") && !struct_out_found) {
        if (strstr(buffer, "Fractional")) coordsAreFractional = true;
      }

      // Input atom info
      // This will be overwritten by the output atom info if the output atom
      // info exists (i.e., if a relaxation was performed).
      if (strstr(buffer, "%block AtomicCoordinatesAndAtomicSpecies") &&
          !struct_out_found) {

        // Clear old atoms from pmol in case they were set for some reason
        vector<OBAtom*> toDelete;
        FOR_ATOMS_OF_MOL(a, *pmol)
          toDelete.push_back(&*a);
        for (size_t i = 0; i < toDelete.size(); i++)
          pmol->DeleteAtom(toDelete.at(i));

        ifs.getline(buffer, BUFF_SIZE);
        tokenize(vs, buffer);

        while (!strstr(buffer, "%endblock AtomicCoordinatesAndAtomicSpecies")) {
          // Find the correct atomic number
          std::map<int, string>::iterator it;
          it = atomTypeLabels.find(atoi(vs.at(3).c_str()));
          // Just a basic find() error check
          if(it == atomTypeLabels.end()) {
             delete cell;
             errorMsg << "Error reading AtomicSpecies\n";
             obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
             pmol->EndModify();
             return false;
          }

          atomicNum = OBElements::GetAtomicNum(it->second.c_str());

          x = atof(vs.at(0).c_str());
          y = atof(vs.at(1).c_str());
          z = atof(vs.at(2).c_str());

          // Add atom
          OBAtom *atom = pmol->NewAtom();
          atom->SetAtomicNum(atomicNum);
          vector3 coords (x,y,z);
          atom->SetVector(coords);

          // Reset vars
          ifs.getline(buffer,BUFF_SIZE);
          tokenize(vs, buffer);
          size = vs.size();
        }
      }

      // Output atom info
      // This will overwrite the input atom info
      if (strstr(buffer, "outcoor: Relaxed atomic coordinates") &&
          !struct_out_found) {
        // Check to see if they are fractional coordinates
        coordsAreFractional = false;

        if (strstr(buffer, "fractional")) coordsAreFractional = true;

        // Clear old atoms from pmol. It was previously set from the input
        // section
        vector<OBAtom*> toDelete;
        FOR_ATOMS_OF_MOL(a, *pmol)
          toDelete.push_back(&*a);
        for (size_t i = 0; i < toDelete.size(); i++)
          pmol->DeleteAtom(toDelete.at(i));

        ifs.getline(buffer, BUFF_SIZE);
        tokenize(vs, buffer);
        size_t size = vs.size();
        // All of these are 6 words in length!
        while (size == 6) {
          atomicNum = OBElements::GetAtomicNum(vs.at(5).c_str());

          x = atof(vs.at(0).c_str());
          y = atof(vs.at(1).c_str());
          z = atof(vs.at(2).c_str());

          // Add atom
          OBAtom *atom = pmol->NewAtom();
          atom->SetAtomicNum(atomicNum);
          vector3 coords (x,y,z);
          atom->SetVector(coords);

          // Reset vars
          ifs.getline(buffer,BUFF_SIZE);
          tokenize(vs, buffer);
          size = vs.size();
        }
      }

      // Lattice constant for the input lattice parameters
      if (strstr(buffer, "LatticeConstant") && !struct_out_found) {
        tokenize(vs, buffer);
        latticeConstant = atof(vs.at(1).c_str());

        // TODO: Maybe check for other length units in the future?
        if (strstr(buffer, "Ang")) coordsAreAngstrom = true;
      }

      // input lattice parameters. This will be overwritten if there are
      // relaxed output lattice parameters elsewhere in the file
      // (i.e., if a cell relaxation was performed...)
      if (strstr(buffer, "%block LatticeParameters") && !struct_out_found) {
        ifs.getline(buffer, BUFF_SIZE);
        tokenize(vs, buffer);
        a = atof(vs.at(0).c_str()) * latticeConstant;
        b = atof(vs.at(1).c_str()) * latticeConstant;
        c = atof(vs.at(2).c_str()) * latticeConstant;

        alpha = atof(vs.at(3).c_str());
        beta  = atof(vs.at(4).c_str());
        gamma = atof(vs.at(5).c_str());
        cell->SetData(a, b, c, alpha, beta, gamma);
        cellDataWasSet = true;
      }

      // output lattice parameters
      // This will overwrite the input lattice parameters
      // TODO: there are many "outcell: Unit cell vectors" that appear
      // before the final one is given. Perhaps we should come up with a way
      // to identify the final one so this code block doesn't get called over
      // and over again?
      if (strstr(buffer, "outcell: Unit cell vectors") && !struct_out_found) {
        for (size_t i = 0; i < 3; i++) {
          ifs.getline(buffer, BUFF_SIZE);
          tokenize(vs, buffer);
          for (size_t j = 0; j < 3; j++) {
            cellMatrix(i,j) = atof(vs.at(j).c_str());
          }
        }

        // Build unit cell
        cell->SetData(cellMatrix);
        cellDataWasSet = true;
      }

      // Final energy
      if (strstr(buffer, "Final energy (eV):")) {
        // Loop through the pieces until we get to the total energy
        // We may need to add more energies to record in the future...
        while (!strstr(buffer, "Total")) ifs.getline(buffer, BUFF_SIZE);
        tokenize(vs, buffer);
        pmol->SetEnergy(atof(vs[3].c_str()) * EV_TO_KCAL_PER_MOL);
      }

      // We've reached the end!
      if (strstr(buffer, "End of run:")) EOFReached = true;
    }

    if (!EOFReached) {
      delete cell;
      errorMsg << "Error! The EOF for siesta was not reached. Check the file "
               << "to see if it was saved properly.\n";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      pmol->EndModify();
      return false;
    }

    // Delete the dynamically allocated cell if it was not set
    if (!cellDataWasSet) delete cell;

    // Convert coords to cartesian if needed
    if (coordsAreFractional) {
      FOR_ATOMS_OF_MOL(atom, pmol) {
        atom->SetVector(cell->FractionalToCartesian(atom->GetVector()));
      }
    }

    // set final unit cell
    if (cellDataWasSet) pmol->SetData(cell);

    pmol->EndModify();

    return true;
  }

} //namespace OpenBabel
