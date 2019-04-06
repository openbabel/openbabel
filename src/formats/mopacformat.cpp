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
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/elements.h>
#include <openbabel/internalcoord.h>
#include <openbabel/generic.h>
#include <cstdlib>

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
    vector<double> orbitalEnergies;
    vector<string> orbitalSymmetries; // left empty for now

    // Translation vectors (if present)
    vector3 translationVectors[3];
    int numTranslationVectors = 0;
    int alphaHOMO = 0;

    mol.BeginModify();
    while	(ifs.getline(buffer,BUFF_SIZE))
      {
        // Avoid "FORCE CONSTANT IN CARTESIAN COORDINATES" (PR#3417992)
        if(strstr(buffer,"  CARTESIAN COORDINATES") != NULL)
          {
            // mol.EndModify();
            mol.Clear();
            ifs.getline(buffer,BUFF_SIZE);	// blank

            // could either be columns or real data
            ifs.getline(buffer,BUFF_SIZE);
            tokenize(vs, buffer);
            if (vs.size() != 5 || vs[0][0] != '1') { // first character should be atom 1
              // those were column headings
              ifs.getline(buffer,BUFF_SIZE);	// blank
              ifs.getline(buffer,BUFF_SIZE);
              tokenize(vs,buffer);
            }
            // now we're at real data
            while (vs.size() == 5)
              {
                atom = mol.NewAtom();
                atom->SetAtomicNum(OBElements::GetAtomicNum(vs[1].c_str()));
                x = atof((char*)vs[2].c_str());
                y = atof((char*)vs[3].c_str());
                z = atof((char*)vs[4].c_str());
                atom->SetVector(x,y,z);

                if (!ifs.getline(buffer,BUFF_SIZE))
                  break;
                tokenize(vs,buffer);
              }
          }
        // ANGSTROMS but not DEGREES (cartesians, not angles)
        else if(strstr(buffer,"(ANGSTROMS)") != NULL && strstr(buffer,"(DEGREES)") == NULL)
          { // newer versions don't print CARTESIAN for final geometry
            mol.Clear();
            ifs.getline(buffer,BUFF_SIZE);	// blank
            ifs.getline(buffer,BUFF_SIZE);
            tokenize(vs,buffer);
            while (vs.size() == 8)
              {
                if (strcmp(vs[1].c_str(), "Tv") != 0)
                  {
                    atom = mol.NewAtom();
                    atom->SetAtomicNum(OBElements::GetAtomicNum(vs[1].c_str()));
                    x = atof((char*)vs[2].c_str());
                    y = atof((char*)vs[4].c_str());
                    z = atof((char*)vs[6].c_str());
                    atom->SetVector(x,y,z);
                  }

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
        // Optimized translation vectors:
        else if (strstr(buffer, "FINAL  POINT  AND  DERIVATIVES") != NULL)
          {
            numTranslationVectors = 0; // Reset
            ifs.getline(buffer,BUFF_SIZE);	// blank
            ifs.getline(buffer,BUFF_SIZE);	// column headings
            ifs.getline(buffer,BUFF_SIZE);
            tokenize(vs,buffer);
            while (vs.size() == 8)
              {
                // Skip coords -- these would be overwritten by the later
                // CARTESIAN COORDINATES block anyway
                if (strcmp(vs.at(2).c_str(), "Tv") != 0)
                  {
                    if (!ifs.getline(buffer,BUFF_SIZE))
                      break;
                    tokenize(vs,buffer);
                    continue;
                  }
                const char coord = vs[4].at(0);
                double val = atof(vs[5].c_str());
                bool isZ = false;
                switch (coord) {
                case 'X':
                  x = val;
                  break;
                case 'Y':
                  y = val;
                  break;
                case 'Z':
                  z = val;
                  isZ = true;
                  break;
                default:
                  cerr << "Reading MOPAC Tv values: unknown coordinate '"
                       << coord << "', value: " << val << endl;
                  break;
                }

                if (isZ)
                  translationVectors[numTranslationVectors++].Set(x, y, z);

                if (!ifs.getline(buffer,BUFF_SIZE))
                  break;
                tokenize(vs,buffer);
              }
          }
        else if(strstr(buffer,"NC:NB:NA:I") != NULL) // z-matrix
          {
            mol.Clear();
            vector<OBInternalCoord*> vic;
            vector<unsigned int> indices;
            vic.push_back((OBInternalCoord*)NULL);

            while (ifs.getline(buffer,BUFF_SIZE)) {
              tokenize(vs,buffer);
              if (vs.size() == 0)
                break;
              else if (vs.size() < 11)
                break;

              atom = mol.NewAtom();

              OBInternalCoord *coord = new OBInternalCoord;
              coord->_dst = atof(vs[2].c_str());
              coord->_ang = atof(vs[4].c_str());
              coord->_tor = atof(vs[6].c_str());
              vic.push_back(coord);

              indices.push_back(atoi(vs[8].c_str()));
              indices.push_back(atoi(vs[9].c_str()));
              indices.push_back(atoi(vs[10].c_str()));

              // symbol in column 1
              atom->SetAtomicNum(OBElements::GetAtomicNum(vs[1].c_str()));
            }
            // read the z-matrix

            // now fill in the atom ids into the internal coords
            unsigned int idx = 0;
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
            InternalToCartesian(vic,mol);
            // coordinates should be set
          }
        else if(strstr(buffer,"DOUBLY OCCUPIED LEVELS") != NULL)
          {
            tokenize(vs, buffer);
            if (vs.size() < 9)
              continue;
            alphaHOMO = atoi(vs[8].c_str());
          }
        else if(strstr(buffer,"EIGENVALUES") != NULL)
          {
            ifs.getline(buffer, BUFF_SIZE); // real data
            tokenize(vs, buffer);
            while(vs.size() > 0) { // ends with a blank line
              for (unsigned int orbital = 0; orbital < vs.size(); ++orbital) {
                // orbitals are listed in eV already, no conversion needed
                orbitalEnergies.push_back(atof(vs[orbital].c_str()));
              }
              ifs.getline(buffer, BUFF_SIZE);
              tokenize(vs, buffer);
            }
          }
        else if(strstr(buffer,"FINAL HEAT") != NULL)
          {
            sscanf(buffer,"%*s%*s%*s%*s%*s%lf",&energy);
            mol.SetEnergy(energy);
          }
        else if(strstr(buffer,"ELECTROSTATIC POTENTIAL CHARGES") != NULL)
          {
            hasPartialCharges = true;
            charges.clear(); // Mulliken Charges
            ifs.getline(buffer,BUFF_SIZE);	// blank
            ifs.getline(buffer,BUFF_SIZE);	// column headings
            ifs.getline(buffer,BUFF_SIZE);
            tokenize(vs,buffer);
            if (vs.size() < 1) return false; // timvdm 18/06/2008
            while (vs.size() > 0 && strstr(vs[0].c_str(),"DIPOLE") == NULL)
              {
                if (vs.size() < 3) break;
                atom = mol.GetAtom(atoi(vs[0].c_str()));
                if (atom != NULL)
                  atom->SetPartialCharge(atof(vs[2].c_str()));
                charges.push_back(atof(vs[2].c_str()));

                if (!ifs.getline(buffer,BUFF_SIZE))
                  break;
                tokenize(vs,buffer);
              }
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
            while (vs.size() > 0 && strstr(vs[0].c_str(),"DIPOLE") == NULL)
              {
                if (vs.size() < 3) break;
                atom = mol.GetAtom(atoi(vs[0].c_str()));
                if (atom != NULL)
                  atom->SetPartialCharge(atof(vs[2].c_str()));
                charges.push_back(atof(vs[2].c_str()));

                if (!ifs.getline(buffer,BUFF_SIZE))
                  break;
                tokenize(vs,buffer);
              }
            // Now we should be at DIPOLE line. If missing, break out of block
            // and continue parsing file.
            if (vs.size() == 0 || strstr(vs[0].c_str(), "DIPOLE") != NULL)
              continue;
            if (!ifs.getline(buffer,BUFF_SIZE))	// POINT CHARGE
              continue; // let the outer loop handle this
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
            unsigned int prevModeCount = displacements.size();
            unsigned int newModes = frequencies.size() - displacements.size();
            vector<vector3> displacement;
            for (unsigned int i = 0; i < newModes; ++i) {
              displacements.push_back(displacement);
            }

            ifs.getline(buffer, BUFF_SIZE);
            tokenize(vs, buffer);
            unsigned int modeCount = vs.size();
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
                if (displacements.size() < prevModeCount + i + 1)
                  displacements.push_back(displacement);
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
    if (numTranslationVectors == 3) {
      OBUnitCell* uc = new OBUnitCell;
      uc->SetData(translationVectors[0], translationVectors[1], translationVectors[2]);
      uc->SetOrigin(fileformatInput);
      mol.SetData(uc);
    }

    // Attach orbitals if found
    if (alphaHOMO > 0) {
      OBOrbitalData *od = new OBOrbitalData();
      od->LoadClosedShellOrbitals(orbitalEnergies, orbitalSymmetries, alphaHOMO);
      od->SetOrigin(fileformatInput);
      mol.SetData(od);
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
        "  b  Disable bonding entirely\n\n"
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
  // Here is a to-do list for a more complete MOPAC input reader
  // - cjh 2011-07-02
  //
  // A. Comment lines
  //
  // A comment line begins with * and may be specified anywhere.
  //
  // Status: implemented in the geometry block, not in header
  //
  //
  // B. Header
  //
  // MOPAC supports line continuation for keywords using the special keywords & and +
  // if & is present, keywords continue on the next line
  // & may be specified on lines 1 and 2 only
  // the total length of the header remains fixed at three lines; the number of lines
  // available for description is reduced accordingly
  // If + is present, keywords continue on the next line
  // AND the total length of the header is extended by one line
  // Up to two + may be used
  //
  // References
  // ----------
  // MOPAC 7.1:
  // MOPAC 2009: http://openmopac.net/manual/allkeys.html
  //
  // Status: not implemented
  //
  //
  // C. Processing atom name
  //
  // 1. MOPAC offers some unique atom names
  // In MOPAC 7.1:
  //
  // XX dummy atom - OB already understands this
  //
  // sparkles
  // +	A 100% ionic alkali metal
  // ++	A 100% ionic alkaline earth metal
  // -	A 100% ionic halogen-like atom
  // --	A 100% ionic group VI-like atom.
  // (Section 6.12 of MOPAC 7 Manual)
  //
  // Cb	(Capped bond) A special type of monovalent atom
  //    existing purely to satisfy valence
  // (Section 3.5 of MOPAC 7 Manual)
  //
  // Tv - Translation vector defining 1-D periodicity for polymers
  //
  // In MOPAC 2009:
  // 2. All of the above, plus:
  // +3 - A +3 sparkle
  // -3 - A -3 sparkle
  // Fr - A sparkle with charge  1/2 (NOT Francium!)
  // At - A sparkle with charge -1/2 (NOT Actinium)
  // X  - also a dummy atom
  // D  - Deuterium - OB already understands this
  // T  - Tritium - OB already understands this
  // Tv - up to 3 translation vectors can be specified for periodic cells
  //      in 1D, 2D and 3D
  //
  // Isotopes can be specified with isotopic mass
  // e.g. C13.0034
  //
  // 3. optional atom labels can be specified with ()
  // e.g. "Mg(At center of porphyrin ring)"
  // label is text in () and can be up to 38 characters long
  // it CAN include spaces
  // if it is "+" or "-", this specifies atomic charges in MOPAC 2009
  // Both mass and label can be specified, e.g. C1(on C5)34.96885
  //
  // In MOPAC 7.x only the Z-matrix format is documented to support labels
  // but in MOPAC 2009 labels are officially supported in all formats
  //
  // References
  // ----------
  // MOPAC 7.x: http://nova.colombo58.unimi.it/manual/pdf/Mopac7.pdf
  // MOPAC 2009: http://openmopac.net/manual/Labels.html
  //
  // Status: atom labels and isotopes recognized but thrown away
  //
  //
  // D. Mixed coordinate format
  //
  // MOPAC2009 supports mixed internal and Cartesian coordinate specification
  // but the code as it stands will currently fail to process the coordinates
  // correctly in this forma
  //
  // Status: throws error if this format is encountered
  //
  // -cjh 2011-07-02
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
    string str, atomLabel, elementSymbol;
    double x,y,z, isotopeMass;
    OBAtom *atom;
    vector<string> vs;

    // Translation vectors (if present)
    vector3 translationVectors[3];
    int numTranslationVectors = 0;
    //

    ifs.getline(buffer,BUFF_SIZE); // keywords
    ifs.getline(buffer,BUFF_SIZE); // filename
    ifs.getline(buffer,BUFF_SIZE); // title (currently ignored)

    mol.BeginModify();

    while (ifs.getline(buffer,BUFF_SIZE))
      {
        isotopeMass = 0;
        elementSymbol = "";

        //First see if this is a comment line - skip comment lines
        if (buffer[0] == '*') continue;

        //First see if there is a label defined
        tokenize(vs,buffer,"()");
        if (vs.size() > 3) //Only one label allowed per line
          {
            //TODO Replace with correct OBError.ThrowError() call
            cerr << "Invalid format in geometry specification: There appears to be more than one atom label specified!\n";
            return false;
          }
        else if (1 < vs.size() && vs.size() <= 3) //There is a label
          {
            elementSymbol = vs[0];
            atomLabel = vs[1];
            strcpy(buffer,vs[2].c_str());
          }
        else //no label, reset buffer
          strcpy(buffer,vs[0].c_str());

        //Now parse the rest of the line
        //There should be three cases:
        //1. There are 7 tokens and the first token is a number specifying the isotope mass
        //2. There are 7 tokens and the first token is a string containing the element symbol
        //3. There are 6 tokens and the first token is a number specifying the Cartesian x coordinate
        tokenize(vs,buffer);
        if (vs.size() == 0)
          break;
        else if (vs.size() < 6)
          {
            //TODO Replace with correct OBError.ThrowError() call
            cerr << "Invalid format in geometry specification.\n";
            return false;
          }
        else if (vs.size() > 7) //cjh 2011-07-02
          {
            //TODO Replace with correct OBError.ThrowError() call
            cerr << "Mixed Cartesian and internal coordinates are currently not supported.\n";
            return false;
          }
        else if (vs.size() == 7)
          {
            if (elementSymbol == "")
              elementSymbol = vs[0];
            else
              isotopeMass = atof((char*)vs[0].c_str());

            x = atof((char*)vs[1].c_str());
            y = atof((char*)vs[3].c_str());
            z = atof((char*)vs[5].c_str());
          }
        else //vs.size() == 6
          {
            x = atof((char*)vs[0].c_str());
            y = atof((char*)vs[2].c_str());
            z = atof((char*)vs[4].c_str());
          }

        if (elementSymbol == "Tv") //MOPAC translation vector
          {
            translationVectors[numTranslationVectors++].Set(x, y, z);
          }
        else
          {
            atom = mol.NewAtom();
            atom->SetVector(x,y,z); //set coordinates
            //set atomic number
            atom->SetAtomicNum(OBElements::GetAtomicNum(elementSymbol.c_str()));
          }
      }

    // Attach unit cell translation vectors if found
    if (numTranslationVectors > 0) {
      OBUnitCell* uc = new OBUnitCell;
      uc->SetData(translationVectors[0], translationVectors[1], translationVectors[2]);
      uc->SetOrigin(fileformatInput);
      mol.SetData(uc);
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
    bool writeUnitCell = (NULL != pConv->IsOption("u", OBConversion::OUTOPTIONS));
    string defaultKeywords = "PUT KEYWORDS HERE";

    if(keywords)
      defaultKeywords = keywords;

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
    else {
      ofs << defaultKeywords;
      if (mol.GetTotalCharge() != 0)
        ofs << " CHARGE=" << mol.GetTotalCharge();
      
      // should handle GetTotalSpinMultiplicity() too
      ofs << endl;
    }

    ofs << mol.GetTitle() << endl;
    ofs << endl; // comment

    string str,str1;
    FOR_ATOMS_OF_MOL(atom, mol)
      {
        snprintf(buffer,BUFF_SIZE,"%-3s%8.5f 1 %8.5f 1 %8.5f 1",
                 OBElements::GetSymbol(atom->GetAtomicNum()),
                 atom->GetX(),
                 atom->GetY(),
                 atom->GetZ());
        ofs << buffer << "\n";
      }

    OBUnitCell *uc = (OBUnitCell*)mol.GetData(OBGenericDataType::UnitCell);
    if (uc && writeUnitCell) {
      //      uc->FillUnitCell(&mol); // complete the unit cell with symmetry-derived atoms

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
    vector<unsigned int> indices;
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

      atom->SetAtomicNum(OBElements::GetAtomicNum(vs[0].c_str()));
    }

    unsigned int idx = 0;
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

      strncpy(type, OBElements::GetSymbol(atom->GetAtomicNum()), 16);
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
