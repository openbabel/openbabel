/**********************************************************************
Copyright (C) 2001-2006 by Geoffrey R. Hutchison
Some portions Copyright (C) 2004 by Chris Morley
Some portions Copyright (C) 2009 by Michael Banck

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


// Required for imaginary frequencies detection
#include <cmath>
// Required for TS detection in ZTS calculation
#include <algorithm>
#define HARTREE_TO_KCAL 627.509469
#define AU_TO_ANGSTROM 0.529177249
#define EV_TO_NM(x) 1239.84193/x

using namespace std;
namespace OpenBabel
{

  class NWChemOutputFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    NWChemOutputFormat()
    {
      OBConversion::RegisterFormat("nwo",this);
    }

    virtual const char* Description() //required
    {
      return
        "NWChem output format\n"
        "Read Options e.g. -as\n"
        " s  Output single bonds only\n"
        " f  Overwrite molecule if more than one\n"
        "    calculation with different molecules\n"
        "    is present in the output file\n"
        "    (last calculation will be prefered)\n"
        " b  Disable bonding entirely\n\n";
    };

    virtual const char* SpecificationURL()
    {return "http://www.emsl.pnl.gov/docs/nwchem/";}; //optional

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
      return READONEONLY | NOTWRITABLE;
    };

    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);

  private:
    void ReadCoordinates(istream* ifs, OBMol* molecule);
    void ReadPartialCharges(istream* ifs, OBMol* molecule);
    void ReadOrbitals(istream* ifs, OBMol* molecule);
    void ReadMultipoleMoment(istream* ifs, OBMol* molecule);

    void ReadFrequencyCalculation(istream* ifs, OBMol* molecule);
    void ReadGeometryOptimizationCalculation(istream* ifs, OBMol* molecule);
    void ReadSinglePointCalculation(istream* ifs, OBMol* molecule);
    void ReadZTSCalculation(istream* ifs, OBMol* molecule);
    void ReadTDDFTCalculation(istream* ifs, OBMol* molecule);
    void ReadMEPCalculation(istream* ifs, OBMol* molecule);
    void ReadNEBCalculation(istream* ifs, OBMol* molecule);
  };

static const char* COORDINATES_PATTERN = "Output coordinates";
static const char* GEOMETRY_OPTIMIZATION_PATTERN = "NWChem Geometry Optimization";
static const char* PROPERTY_CALCULATION_PATTERN = "NWChem Property Module";
static const char* ZTS_CALCULATION_PATTERN = " String method.";
static const char* NEB_CALCULATION_PATTERN = "NWChem Minimum Energy Pathway Program (NEB)";
static const char* PYTHON_CALCULATION_PATTERN = "NWChem Python program";
static const char* ESP_CALCULATION_PATTERN = "NWChem Electrostatic Potential Fit Module";
static const char* SCF_CALCULATION_PATTERN = "SCF Module";
static const char* DFT_CALCULATION_PATTERN = "DFT Module";
static const char* TDDFT_CALCULATION_PATTERN = "TDDFT Module";
static const char* MEP_CALCULATION_PATTERN = "Gonzalez & Schlegel IRC Optimization";
static const char* SCF_ENERGY_PATTERN = "SCF energy =";
static const char* DFT_ENERGY_PATTERN = "DFT energy =";
static const char* FREQUENCY_PATTERN = "NWChem Nuclear Hessian and Frequency Analysis";
static const char* OPTIMIZATION_STEP_PATTERN = "Step       Energy";
static const char* VIBRATIONS_TABLE_PATTERN = "P.Frequency";
static const char* INTENSITIES_TABLE_PATTERN = "Projected Infra Red Intensities";
static const char* DIGITS = "1234567890";
static const char* END_OF_CALCULATION_PATTERN = "times  cpu";
static const char* ORBITAL_START_PATTERN = "Vector";
static const char* ORBITAL_SECTION_PATTERN_1 = "Analysis";
static const char* ORBITAL_SECTION_PATTERN_2 = "rbital";
static const char* BETA_ORBITAL_PATTERN = "Beta";
static const char* MULLIKEN_CHARGES_PATTERN = "Mulliken analysis of the total density";
static const char* GEOMETRY_PATTERN = "Geometry \"geometry\"";
static const char* ZTS_CONVERGED_PATTERN = " The string calculation ";
static const char* NBEADS_PATTERN = " Number of replicas";
static const char* ROOT_PATTERN = "Root";
static const char* OSCILATOR_STRENGTH_PATTERN = "Oscillator Strength";
static const char* SPIN_FORBIDDEN_PATTERN = "Spin forbidden";
static const char* MULTIPOLE_MOMENT_PATTERN = "Multipole analysis of the density";
static const char* MEP_STEP_END_PATTERN = "&  Point";
static const char* NEB_BEAD_START_PATTERN = "neb: running bead";
static const char* NEB_BEAD_ENERGY_PATTERN = "neb: final energy";
static const char* NEB_NBEADS_PATTERN = "number of images in path";
static const char* GRADIENT_PATTERN = "ENERGY GRADIENTS";
// Two spaces are nessesary to avoid matching "IRC Optimization converged"
static const char* OPTIMIZATION_END_PATTERN = "  Optimization converged";

  //Make an instance of the format class
  NWChemOutputFormat theNWChemOutputFormat;

  class NWChemInputFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    NWChemInputFormat()
    {
      OBConversion::RegisterFormat("nw",this);
    }

    virtual const char* Description() //required
    {
      return
        "NWChem input format\n"
        "No comments yet\n";
    };

    virtual const char* SpecificationURL()
    {return "http://www.emsl.pnl.gov/docs/nwchem/";}; //optional

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
  NWChemInputFormat theNWChemInputFormat;

  /////////////////////////////////////////////////////////////////
  /**
  Moves stream (ifs) position to end of calculation.
  */
  static void GotoCalculationEnd(istream* ifs)
  {
  char buffer[BUFF_SIZE];
    while ( (strstr(buffer,END_OF_CALCULATION_PATTERN) == NULL))
        if (!ifs->getline(buffer,BUFF_SIZE))
            break;
  }


  //////////////////////////////////////////////////////
  /**
  Method reads coordinates from input stream (ifs) and
  writes it into supplied OBMol object (molecule).
  Input stream must be set to begining of coordinates
  table in nwo file. (Line after "Output coordinates...")
  Stream will be set at next line after geometry table.
  If one of input arguments is NULL method returns without
  any changes
  If "molecule" already contain geometry - method will write
  new geometry as conformer.
  If "molecule" contain geometry which incompatible with read
  data method returns without changes.
  */
  void NWChemOutputFormat::ReadCoordinates(istream* ifs, OBMol* molecule)
  {
    if ((molecule == NULL) || (ifs == NULL))
        return;
    vector<string> vs;
    char buffer[BUFF_SIZE];
    double x, y, z;
    unsigned int natoms = molecule->NumAtoms();
    bool from_scratch = false;
    double* coordinates;
    if (natoms == 0)
        from_scratch = true;
    else
        coordinates = new double[3*natoms];
    ifs->getline(buffer,BUFF_SIZE);	// blank
    ifs->getline(buffer,BUFF_SIZE);	// column headings
    ifs->getline(buffer,BUFF_SIZE);	// ---- ----- ----
    ifs->getline(buffer,BUFF_SIZE);
    tokenize(vs,buffer);
    unsigned int i=0;
    while (vs.size() == 6)
    {
        x = atof((char*)vs[3].c_str());
        y = atof((char*)vs[4].c_str());
        z = atof((char*)vs[5].c_str());
        if (from_scratch)
        {
            // set atomic number
            OBAtom* atom = molecule->NewAtom();
            atom->SetAtomicNum(atoi(vs[2].c_str()));
            atom->SetVector(x,y,z);
        }
        else
        {
            // check atomic number
            if ((i>=natoms) || (molecule->GetAtom(i+1)->GetAtomicNum() != atoi(vs[2].c_str())))
            {
                delete[] coordinates;
                return;
            }
            coordinates[i*3] = x;
            coordinates[i*3+1] = y;
            coordinates[i*3+2] = z;
            i++;
        }
        if (!ifs->getline(buffer,BUFF_SIZE))
          break;
        tokenize(vs,buffer);
    }
    if (from_scratch) 
    {
        return;
    }
    if (i != natoms) {
        delete[] coordinates;
        return;
    }
    molecule->AddConformer(coordinates);
  }

//////////////////////////////////////////////////////
  /**
  Method reads charge, dipole and quadrupole moment from input stream (ifs)
  and writes them to supplied OBMol object (molecule)
  Input stream must be set to begining of Multipole moment
  section in nwo file. (Line after "Multipole analysis of the density")
  Stream will be set to the end of multipole moment section.
  */
  void NWChemOutputFormat::ReadMultipoleMoment(istream* ifs, OBMol* molecule)
  {
    if ((ifs == NULL) || (molecule == NULL))
        return;

    char buffer[BUFF_SIZE];
    vector<string> vs;
    matrix3x3 quadrupole;
    double dipole[3];
    int charge;
    bool blank_line = false;

    ifs->getline(buffer, BUFF_SIZE); // -------
    ifs->getline(buffer, BUFF_SIZE); // blank
    ifs->getline(buffer, BUFF_SIZE); // Header
    ifs->getline(buffer, BUFF_SIZE); // -------

    while (ifs->getline(buffer, BUFF_SIZE))
    {
        tokenize(vs, buffer);
        // L   x y z        total         alpha         beta         nuclear
        // L   x y z        total         open         nuclear
        // 0   1 2 3          4             5            6             7
        if (vs.size() < 7)
        {
            if (blank_line)
            {
                molecule->SetTotalCharge(charge);
                OBVectorData* dipole_moment = new OBVectorData;
                dipole_moment->SetData(vector3(dipole));
                dipole_moment->SetAttribute("Dipole Moment");
                molecule->SetData(dipole_moment);
                OBMatrixData* quadrupole_moment = new OBMatrixData;
                quadrupole_moment->SetData(quadrupole);
                quadrupole_moment->SetAttribute("Quadrupole Moment");
                molecule->SetData(quadrupole_moment);
                return;
            }
            // Second blank line means end of multipole section
            blank_line = true;
            continue;
        }
        blank_line = false;
        if (vs[0][0] == '0')
            charge = atoi(vs[4].c_str());
        else if (vs[0][0] == '1')
            for (unsigned int i = 0; i < 3; i++)
                if (vs[i+1][0] == '1')
                    dipole[i] = atof(vs[4].c_str());
        else if (vs[0][0] == '2')
        {
            double value = atof(vs[4].c_str());
            unsigned int i[2], j = 0;
            for (unsigned int k = 0 ; k<3; k++)
            {
                if (vs[k+1][0] == '2')
                    i[0] = i[1] = k; // Diagonal elements
                else if (vs[k+1][0] == '1')
                    i[j++] = k;
            }
            quadrupole.Set(i[0], i[1], value);
            quadrupole.Set(i[1], i[0], value);
        }
        else
            return;
    }
  }

  //////////////////////////////////////////////////////
  /**
  Method reads UV Spectra from input stream (ifs)
  and writes them to supplied OBMol object (molecule)
  Input stream must be set to begining of TDDFT
  calculation in nwo file. (Line after "NWChem TDDFT Module")
  Stream will be set to the end of calculation.
  */
  void NWChemOutputFormat::ReadTDDFTCalculation(istream* ifs, OBMol* molecule)
  {
    if ((ifs == NULL) || (molecule == NULL))
        return;

    char buffer[BUFF_SIZE];
    vector<string> vs;
    vector<double> wavelengths;
    vector<double> oscilator_strengths;

    while (ifs->getline(buffer, BUFF_SIZE))
    {
        if (strstr(buffer, ROOT_PATTERN) != NULL)
        {
            tokenize(vs, buffer);
            //  Root   1 singlet b2             0.294221372 a.u.                8.0062 eV
            //   0     1    2    3                  4        5                    6    7
            if (vs.size() < 8)
                break;
            wavelengths.push_back(EV_TO_NM(atof(vs[6].c_str())));
        }
        else if (strstr(buffer, OSCILATOR_STRENGTH_PATTERN) != NULL)
        {
            if (strstr(buffer, SPIN_FORBIDDEN_PATTERN) != NULL)
                oscilator_strengths.push_back(0);
            else
            {
                tokenize(vs, buffer);
                // Dipole Oscillator Strength                         0.01418
                //   0        1         2                                3
                if (vs.size() < 4)
                    break;
                oscilator_strengths.push_back(atof(vs[3].c_str()));
            }
        }
        else if (strstr(buffer, END_OF_CALCULATION_PATTERN) != NULL)
            break;
    }
    if (wavelengths.size() != oscilator_strengths.size())
        return;
    OBElectronicTransitionData* et_data = new OBElectronicTransitionData;
    et_data->SetData(wavelengths, oscilator_strengths);
    molecule->SetData(et_data);
  }

  //////////////////////////////////////////////////////
  /**
  Method reads partial charges from input stream (ifs)
  and writes them to supplied OBMol object (molecule)
  Input stream must be set to begining of charges
  table in nwo file. (Line after "Mulliken analysis of the total density")
  Stream will be set at next line after charges table.
  If reading charges failed or "molecule" contains
  data incompatible with read charges then "molecule"
  wont be changed.
  */
  void NWChemOutputFormat::ReadPartialCharges(istream* ifs, OBMol* molecule)
  {
    if ((molecule == NULL) || (ifs == NULL))
        return;
    vector<string> vs;
    char buffer[BUFF_SIZE];
    bool from_scratch = false;
    vector<int> charges;
    vector<double> partial_charges;
    unsigned int natoms = molecule->NumAtoms();

    if (natoms == 0)
        from_scratch = true;
    ifs->getline(buffer,BUFF_SIZE); // ---- ----- ----
    ifs->getline(buffer,BUFF_SIZE);	// blank
    ifs->getline(buffer,BUFF_SIZE);	// column headings
    ifs->getline(buffer,BUFF_SIZE);	// ---- ----- ----
    ifs->getline(buffer,BUFF_SIZE);
    tokenize(vs, buffer);

    // N Symbol    Charge     PartialCharge+Charge   ShellCharges
    // 0   1          2                3                4,etc
    unsigned int i = 1;
    while (vs.size() >= 4)
    {
        int charge = atoi(vs[2].c_str());
        if (!from_scratch)
        {
            if (i > natoms)
                return;
            if (molecule->GetAtom(i++)->GetAtomicNum() != charge)
                return;
        }
        else
            charges.push_back(charge);
        partial_charges.push_back(atof(vs[3].c_str()) - charge);
        ifs->getline(buffer,BUFF_SIZE);
        tokenize(vs, buffer);
    }

    if (from_scratch)
        molecule->ReserveAtoms(partial_charges.size());
    else if (partial_charges.size() != natoms)
        return;
    for(unsigned int j=0;j<partial_charges.size();j++)
    {
        OBAtom* atom;
        if (from_scratch)
        {
            atom = molecule->NewAtom();
            atom->SetAtomicNum(charges[j]);
        }
        else
        {
            atom = molecule->GetAtom(j+1);
        }
        atom->SetPartialCharge(partial_charges[j]);
    }
  }


  //////////////////////////////////////////////////////
  /**
  Method reads orbital information from input stream (ifs)
  and writes them to supplied OBMol object (molecule).
  Input stream must be set to begining of orbital data
  section in nwo file. (Line after "... Molecular Orbital Analysis")
  Stream will be set at next line after end of orbital section.
  */
  void NWChemOutputFormat::ReadOrbitals(istream* ifs, OBMol* molecule)
  {
    if ((ifs == NULL) || (molecule == NULL))
        return;
    vector<string> vs;
    char buffer[BUFF_SIZE];
    vector<OBOrbital> orbitals;
    OBOrbitalData* orbital_data = new OBOrbitalData;
    ifs->getline(buffer, BUFF_SIZE); // ---------
    ifs->getline(buffer, BUFF_SIZE); // blank line

    while (ifs->getline(buffer,BUFF_SIZE))
    {
        if (strstr(buffer, ORBITAL_START_PATTERN))
        {
            tokenize(vs, buffer);
            // Vector   N  Occ=X  E= Y  Symmetry=a'
            //   0      1    2    3  4  5(optional)
            if (vs.size() < 5)
                break; // Orbital data is broken

            double energy = atof(vs[4].c_str()) * HARTREE_TO_KCAL;
            double occupation = atof(vs[2].c_str()+4); // Start from symbol after '='
            string symbol;
            if (vs.size() > 5)
                symbol = vs[5].substr(9, string::npos);
            else
                symbol = " "; // Symmetry is unknown
            OBOrbital orbital;
            orbital.SetData(energy, occupation, symbol);
            orbitals.push_back(orbital);

            ifs->getline(buffer, BUFF_SIZE); // MO Center ...
            ifs->getline(buffer, BUFF_SIZE); // Table header
            ifs->getline(buffer,BUFF_SIZE); // ----------
            while (ifs->getline(buffer,BUFF_SIZE))
                if (strlen(buffer) < 2) // If blank line detected
                    break;
        }// if Vector ...
        else if ((strstr(buffer, ORBITAL_SECTION_PATTERN_2) != NULL)&&(strstr(buffer, ORBITAL_SECTION_PATTERN_1) != NULL))
        {
            orbital_data->SetAlphaOrbitals(orbitals);
            orbital_data->SetOpenShell(true);
            orbitals.clear();
            ifs->getline(buffer, BUFF_SIZE); // ---------
            ifs->getline(buffer, BUFF_SIZE); // blank line
        }// if beta orbital section found
        else
        {
            if (orbital_data->IsOpenShell())
                orbital_data->SetBetaOrbitals(orbitals);
            else
                orbital_data->SetAlphaOrbitals(orbitals);
            molecule->SetData(orbital_data);
            return;
        }
    }
  delete orbital_data;
  }


  //////////////////////////////////////////////////////
  /**
  Method reads IRC steps from input stream (ifs)
  and writes it to supplied OBMol object (molecule).
  Input stream must be set to begining of Minimal Energy
  Path IRC calculation in nwo file.
  (Line after "Gonzalez & Schlegel IRC Optimization")
  Method wont work if "molecule" already contains data
  about conformers.
  After all stream will be set at the end of calculation.
  */
  void NWChemOutputFormat::ReadMEPCalculation(istream* ifs, OBMol* molecule)
  {
    if ((molecule == NULL) || (ifs == NULL))
        return;
    if (molecule->NumConformers() > 0)
        return;

    vector<string> vs;
    char buffer[BUFF_SIZE];
    vector<double> energies;

    while (ifs->getline(buffer, BUFF_SIZE))
    {
        if(strstr(buffer, OPTIMIZATION_END_PATTERN) != NULL)
        {
            while(ifs->getline(buffer, BUFF_SIZE))
            {
                if (strstr(buffer, COORDINATES_PATTERN))
                    ReadCoordinates(ifs, molecule);
                else if (strstr(buffer, OPTIMIZATION_STEP_PATTERN))
                {
                    ifs->getline(buffer, BUFF_SIZE); // ------
                    ifs->getline(buffer, BUFF_SIZE);
                    tokenize(vs, buffer);
                    molecule->SetConformer(molecule->NumConformers() - 1);
                    if (vs.size() > 2) // @ NStep   Energy...
                        energies.push_back(atof(vs[2].c_str()) * HARTREE_TO_KCAL);
                }
                else if (strstr(buffer, MULTIPOLE_MOMENT_PATTERN) != NULL)
                    ReadMultipoleMoment(ifs, molecule);
                else if (strstr(buffer, MEP_STEP_END_PATTERN) != NULL)
                    break;
            }
        }
        else if(strstr(buffer, END_OF_CALCULATION_PATTERN) != NULL)
            break;
    }
    if (energies.size() != molecule->NumConformers())
    {
        cerr << "Number of read energies (" << energies.size();
        cerr << ") does not match number of read conformers (";
        cerr << molecule->NumConformers() << ")!" << endl;
        return;
    }
    molecule->SetEnergies(energies);
  }


  //////////////////////////////////////////////////////
  /**
  Method reads optimization steps from input stream (ifs)
  and writes it to supplied OBMol object (molecule).
  Input stream must be set to begining of geometry optimization
  calculation in nwo file. (Line after "NWChem Geometry Optimization")
  If no geometry data found then "molecule" wont be changed.
  After all stream will be set at the end of calculation.
  */
  void NWChemOutputFormat::ReadGeometryOptimizationCalculation(istream* ifs, OBMol* molecule)
  {
    if ((molecule == NULL) || (ifs == NULL))
        return;
    vector<string> vs;
    char buffer[BUFF_SIZE];
    vector<double> energies;

    while (ifs->getline(buffer, BUFF_SIZE))
    {
        if(strstr(buffer,COORDINATES_PATTERN) != NULL)
        {
            ReadCoordinates(ifs, molecule);
            molecule->SetConformer(molecule->NumConformers() - 1);
        }
        else if ((strstr(buffer, ORBITAL_SECTION_PATTERN_2) != NULL)&&(strstr(buffer, ORBITAL_SECTION_PATTERN_1) != NULL))
            ReadOrbitals(ifs, molecule);
        else if(strstr(buffer, OPTIMIZATION_STEP_PATTERN) != NULL)
        {
            // Extract energy
            ifs->getline(buffer, BUFF_SIZE); // ------
            ifs->getline(buffer, BUFF_SIZE);
            tokenize(vs, buffer);
            molecule->SetConformer(molecule->NumConformers() - 1);
            if (vs.size() > 2) // @ NStep   Energy...
                energies.push_back(atof(vs[2].c_str()) * HARTREE_TO_KCAL);
        }
        else if(strstr(buffer, MULTIPOLE_MOMENT_PATTERN) != NULL)
            ReadMultipoleMoment(ifs, molecule);
        else if(strstr(buffer, MULLIKEN_CHARGES_PATTERN) != NULL)
            ReadPartialCharges(ifs, molecule);
        else if(strstr(buffer, END_OF_CALCULATION_PATTERN) != NULL)
            break;
    }
    vector<double> old_energies = molecule->GetEnergies();
    old_energies.reserve(old_energies.size() + energies.size());
    old_energies.insert(old_energies.end(), energies.begin(), energies.end());
    molecule->SetEnergies(old_energies);
  }

  //////////////////////////////////////////////////////
  /**
  Method reads vibration data and all other avalible data
  from input stream (ifs) and writes it to supplied OBMol
  object (molecule).
  If any of arguments are NULL method will quit without changes.
  If molecule does not contain geometry data method quits
  without changes.
  Input stream must be set to begining of frequency
  calculation in nwo file.
  (Line after "NWChem Nuclear Hessian and Frequency Analysis")
  If vibration data not found then only avalible data will be
  attached.
  Input stream will be set at the end of calculation.
  */
  void NWChemOutputFormat::ReadFrequencyCalculation(istream* ifs, OBMol* molecule)
  {
    if ((ifs == NULL) || (molecule == NULL))
        return;
    if (molecule->NumAtoms() == 0)
        return;
    OBVibrationData* vibration_data = NULL;
    vector<double>  Frequencies, Intensities;
    vector<vector<vector3> > Lx;
    vector<string> vs;
    char buffer[BUFF_SIZE];

    while (ifs->getline(buffer, BUFF_SIZE))
    {
        if (strstr(buffer, VIBRATIONS_TABLE_PATTERN) != NULL)
        {
            vector<double> freq;
            vector<vector<vector3> > vib;
            // freq and vib are auxiliary vectors which hold the data for
            // every block of 6 vibrations.
            tokenize(vs,buffer);
            for(unsigned int i=1; i<vs.size(); ++i)
            {
                vib.push_back(vector<vector3>());
                freq.push_back(atof(vs[i].c_str()));
            }
            ifs->getline(buffer,BUFF_SIZE);     // blank line
            ifs->getline(buffer,BUFF_SIZE);
            tokenize(vs,buffer);
            while(vs.size() > 2)
            {
                vector<double> x, y, z;
                for (unsigned int i = 1; i < vs.size(); i++)
                    x.push_back(atof(vs[i].c_str()));
                ifs->getline(buffer, BUFF_SIZE);
                tokenize(vs,buffer);
                for (unsigned int i = 1; i < vs.size(); i++)
                    y.push_back(atof(vs[i].c_str()));
                ifs->getline(buffer, BUFF_SIZE);
                tokenize(vs,buffer);
                for (unsigned int i = 1; i < vs.size(); i++)
                    z.push_back(atof(vs[i].c_str()));
                ifs->getline(buffer, BUFF_SIZE);
                tokenize(vs,buffer);
                if (x.size() == y.size() && y.size() == z.size()) {
                  // make sure the arrays are equal or we'll crash
                  // not sure how to recover if it's not true
                  for (unsigned int i = 0; i < freq.size(); i++)
                  {
                    vib[i].push_back(vector3(x[i], y[i], z[i]));
                  }
                }
            }// while vs.size() > 2
            for (unsigned int i = 0; i < freq.size(); i++)
            {
              if (abs(freq[i]) > 10.0) {
                Frequencies.push_back(freq[i]);
                Lx.push_back(vib[i]);
              }
            }// for (unsigned int i = 0; i < freq.size(); i++)
        }// if P.Frequency
        else if(strstr(buffer, INTENSITIES_TABLE_PATTERN) != NULL)
        {
            ifs->getline(buffer, BUFF_SIZE); // table header
            ifs->getline(buffer, BUFF_SIZE); // table delimiter
            ifs->getline(buffer, BUFF_SIZE);
            tokenize(vs,buffer);
            while (vs.size() == 7)
            {
                if (abs(atof(vs[1].c_str())) > 10.0)
                    Intensities.push_back(atof(vs[5].c_str()));
                ifs->getline(buffer, BUFF_SIZE);
                tokenize(vs,buffer);
            }
        } // if "Projected Infra Red Intensities"
        else if(strstr(buffer, MULLIKEN_CHARGES_PATTERN) != NULL)
            ReadPartialCharges(ifs, molecule);
        else if(strstr(buffer, MULTIPOLE_MOMENT_PATTERN) != NULL)
            ReadMultipoleMoment(ifs, molecule);
        else if ((strstr(buffer, ORBITAL_SECTION_PATTERN_2) != NULL)&&(strstr(buffer, ORBITAL_SECTION_PATTERN_1) != NULL))
            ReadOrbitals(ifs, molecule);
        else if(strstr(buffer, END_OF_CALCULATION_PATTERN) != NULL) // End of task
            break;
    }
    if (Frequencies.size() == 0)
        return;

    vibration_data = new OBVibrationData;
    vibration_data->SetData(Lx, Frequencies, Intensities);
    molecule->SetData(vibration_data);
  }

  /////////////////////////////////////////////////////////////////
  /**
  Method reads single point energy and all avalible data from input
  stream (ifs) and writes it to supplied OBMol object (molecule)
  Input stream must be set to begining of energy calculation
  in nwo file. (Line after "NWChem <theory> Module")
  If energy not found then "molecule" wont be changed.
  */
  void NWChemOutputFormat::ReadSinglePointCalculation(istream* ifs, OBMol* molecule)
  {
    if ((molecule == NULL) || (ifs == NULL))
        return;
    double energy;
    vector<string> vs;
    char buffer[BUFF_SIZE];

    while (ifs->getline(buffer, BUFF_SIZE))
    {
        if ((strstr(buffer, DFT_ENERGY_PATTERN) != NULL) || (strstr(buffer, SCF_ENERGY_PATTERN) != NULL))
        {
            tokenize(vs, buffer);
            energy = atof(vs[4].c_str()) * HARTREE_TO_KCAL;
        }
        else if ((strstr(buffer, ORBITAL_SECTION_PATTERN_2) != NULL)&&(strstr(buffer, ORBITAL_SECTION_PATTERN_1) != NULL))
            ReadOrbitals(ifs, molecule);
        else if(strstr(buffer, MULTIPOLE_MOMENT_PATTERN) != NULL)
            ReadMultipoleMoment(ifs, molecule);
        else if (strstr(buffer, MULLIKEN_CHARGES_PATTERN) != NULL)
            ReadPartialCharges(ifs, molecule);
        else if (strstr(buffer, TDDFT_CALCULATION_PATTERN) != NULL)
            ReadTDDFTCalculation(ifs, molecule);
        else if (strstr(buffer, END_OF_CALCULATION_PATTERN) != NULL)
            break;
    }
    if (energy == 0)
        return;
    molecule->SetEnergy(energy);
  }

  /**
  Method reads beads and their energies from NEB calculation from
  input stream (ifs) and writes them to supplied OBMol object (molecule)
  Input stream must be set to begining of NEB calculation
  in nwo file. (Line after "NWChem Minimum Energy Pathway Program (NEB)")
  If method failed then "molecule" wont be changed.
  */
  void NWChemOutputFormat::ReadNEBCalculation(istream* ifs, OBMol* molecule)
  {
    if ((ifs == NULL) || (molecule == NULL))
        return;
    unsigned int natoms = molecule->NumAtoms();
    // Inital geometry must be supplied
    if (natoms == 0)
        return;
    char buffer[BUFF_SIZE];
    vector<string> vs;
    vector<double*> beads;
    vector<double> energies;
    unsigned int nbeads = 0;
    unsigned int current_bead = UINT_MAX;

    while(ifs->getline(buffer, BUFF_SIZE))
    {
        if (strstr(buffer, NEB_BEAD_START_PATTERN) != NULL)
        {
            tokenize(vs, buffer);
            // neb: running bead                    N
            //  0      1      2                     3
            if (vs.size() < 4)
                break;
            current_bead = atoi(vs[3].c_str()) - 1;
            // Bead index in array starts from 0
            // but in log it starts from 1
        }
        else if (strstr(buffer, NEB_BEAD_ENERGY_PATTERN) != NULL)
        {
            tokenize(vs, buffer);
            // neb: final energy  N
            //  0     1      2    3
            if (vs.size() < 4)
                break;
            if (current_bead >= nbeads)
            {
                cerr << "Current bead out of range: " << current_bead << " of " << nbeads << endl;
                break;
            }
            energies[current_bead] = atof(vs[3].c_str());
        }
        else if (strstr(buffer, GRADIENT_PATTERN) != NULL)
        {
            ifs->getline(buffer, BUFF_SIZE); // blank line
            ifs->getline(buffer, BUFF_SIZE); // 1st level header
            ifs->getline(buffer, BUFF_SIZE); // 2nd level header
            for (unsigned int i = 0; i<natoms; i++)
            {
                ifs->getline(buffer, BUFF_SIZE);
                tokenize(vs, buffer);
                // N Symbol     x   y  z    x_grad  y_grad  z_grad
                // 0   1        2   3  4       5      6       7
                if (vs.size() < 8)
                    break;
                unsigned int end_of_symbol = vs[1].find_last_not_of(DIGITS) + 1;
                if (OBElements::GetAtomicNum(vs[1].substr(0, end_of_symbol).c_str()) != molecule->GetAtom(i+1)->GetAtomicNum())
                    break;
                if (current_bead >= nbeads)
                {
                    cerr << "Current bead out of range: " << current_bead << " of " << nbeads << endl;
                    break;
                }
                beads[current_bead][i*3] = atof(vs[2].c_str())*AU_TO_ANGSTROM;
                beads[current_bead][1+i*3] = atof(vs[3].c_str())*AU_TO_ANGSTROM;
                beads[current_bead][2+i*3] = atof(vs[4].c_str())*AU_TO_ANGSTROM;
            }
        }
        else if (strstr(buffer, NEB_NBEADS_PATTERN) != NULL)
        {
            tokenize(vs, buffer);
            // number of images in path         (nbeads) =   N
            //   0    1     2    3   4             5     6   7
            if (vs.size() < 8)
                break;
            nbeads = atoi(vs[7].c_str());
            beads.reserve(nbeads);
            energies.reserve(nbeads);
            for (unsigned int i = 0;i<nbeads;i++)
            {
                beads.push_back(new double[natoms*3]);
                energies.push_back(0.0);
            }
        }
        else if (strstr(buffer, END_OF_CALCULATION_PATTERN) != NULL)
        {
            molecule->SetConformers(beads);
            molecule->SetEnergies(energies);
            return;
        }
    }
    cerr << "Failed to read NEB calculation!" << endl;
    for(unsigned int i = 0; i < beads.size();i++)
        delete beads[i];
  }

  /////////////////////////////////////////////////////////////////
  /**
  Method reads beads and their energies from ZTS calculation from
  input stream (ifs) and writes them to supplied OBMol object (molecule)
  Input stream must be set to begining of ZTS calculation
  in nwo file. (Line after "@ String method.")
  If method failed then "molecule" wont be changed.
  */
  void NWChemOutputFormat::ReadZTSCalculation(istream* ifs, OBMol* molecule)
  {
    if ((ifs == NULL) || (molecule == NULL))
        return;
    unsigned int natoms = molecule->NumAtoms();
    // Inital geometry must be supplied
    if (natoms == 0)
        return;
    char buffer[BUFF_SIZE];
    vector<string> vs;
    vector<double*> beads;
    vector<double> energies;
    unsigned int nbeads;
    while(ifs->getline(buffer, BUFF_SIZE))
    {
        if (strstr(buffer, NBEADS_PATTERN) != NULL)
        {
            tokenize(vs, buffer);
            // @ Number of replicas   =        24
            // 0   1     2    3       4        5
            if (vs.size() < 6)
                break; // Line with number of beads is incomplete
            nbeads = atoi(vs[5].c_str());
            beads.reserve(nbeads);
        }// @ Number of replicas
        else if (strstr(buffer, ZTS_CONVERGED_PATTERN) != NULL)
        {
            // NWChem does not mark end in this type of calculation,
            // so end will be there, where all nessesary data have
            // obtained
            ifs->getline(buffer, BUFF_SIZE); // blank line
            ifs->getline(buffer, BUFF_SIZE);
            // @ Bead number =     <N>  Potential Energy =     <Energy>
            // 0  1     2    3      4       5       6    7        8
            tokenize(vs, buffer);
            // Thanks to the commit jeffhammond/nwchem@76d2b8c the beads
            // output was broken (in nwchem 6.6+ there is no equal sign after
            // 'number'. So all indicies will be counted from the end.
            unsigned int vsize = vs.size();
            while (vsize > 7)
            {
                unsigned int bead_number = atoi(vs[vsize-5].c_str());
                double bead_energy = atof(vs[vsize-1].c_str()) * HARTREE_TO_KCAL;
                ifs->getline(buffer, BUFF_SIZE); // natoms
                if (atoi(buffer) != natoms)
                    break; // table contains geometry of different molecule
                ifs->getline(buffer, BUFF_SIZE); // comment
                double* bead = new double[natoms*3];
                for(unsigned int i = 0; i<natoms; i++)
                {
                    ifs->getline(buffer, BUFF_SIZE);
                    tokenize(vs, buffer);
                    //  Symbol              X     Y     Z
                    //    0                 1     2     3
                    if ((vs.size() < 4) || (molecule->GetAtom(i+1)->GetAtomicNum() != OBElements::GetAtomicNum(vs[0].c_str())))
                        break; // molecule has no such atom or table row incomplete

                    unsigned int atom_idx = i*3;
                    bead[atom_idx] = atof(vs[1].c_str()); // X
                    bead[atom_idx+1] = atof(vs[2].c_str()); // Y
                    bead[atom_idx+2] = atof(vs[3].c_str()); // Z
                }
                beads.push_back(bead);
                energies.push_back(bead_energy);
                ifs->getline(buffer, BUFF_SIZE);
                tokenize(vs, buffer);
                if (vs.size() <= 1) // blank line
                {
                    // Looks like it's end of calculation.
                    if (bead_number != nbeads)
                        break;
                    molecule->SetEnergies(energies);
                    molecule->SetConformers(beads);
                    unsigned int ts_position = distance(energies.begin(), max_element(energies.begin(), energies.end()));
                    molecule->SetConformer(ts_position);
                    return;
                }
            }
            break;// It is the end of calculation anyway
        }//@ Bead number
        else if (strstr(buffer, END_OF_CALCULATION_PATTERN) != NULL)
        {
            // End of all calculations still required to handle
            molecule->SetEnergies(energies);
            molecule->SetConformers(beads);
            unsigned int ts_position = distance(energies.begin(), max_element(energies.begin(), energies.end()));
            molecule->SetConformer(ts_position);
            return;
        }
    }
    // Something went wrong. Do some cleanup and exit
    for(unsigned int i = 0; i < beads.size();i++)
        delete beads[i];
  }

  /////////////////////////////////////////////////////////////////
  bool NWChemOutputFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title = pConv->GetTitle();

    char buffer[BUFF_SIZE];

    mol.BeginModify();

    // Reading inital parameters of calculation such as
    // used theory and calculation type for better
    // recognition futher output
    while	(ifs.getline(buffer,BUFF_SIZE))
    {
        if(strstr(buffer,GEOMETRY_PATTERN) != NULL)
        {
            // Input coordinates for calculation
            if ((mol.NumAtoms() == 0) || (pConv->IsOption("f",OBConversion::INOPTIONS) != NULL))
            {
                // If coordinates had redefined while calculation
                // in input file and "f" option had supplied then overwrite
                // all previous calculations. Otherwise calculations for
                // new geometry will be considered as new molecule.
                mol.Clear();
                mol.BeginModify();
                ifs.getline(buffer,BUFF_SIZE);// -------------------------
                ifs.getline(buffer,BUFF_SIZE);// blank
                ifs.getline(buffer,BUFF_SIZE);// Output coordinates...
                ReadCoordinates(&ifs, &mol);
            }
            else
            {
                int i;
                for(i=0; buffer[i] != '\0';i++);
                ifs.seekg(-i, ios_base::cur);
                break;
            }
        }
        else if(strstr(buffer, GEOMETRY_OPTIMIZATION_PATTERN) != NULL)
            ReadGeometryOptimizationCalculation(&ifs, &mol);
        else if(strstr(buffer, FREQUENCY_PATTERN) != NULL)
            ReadFrequencyCalculation(&ifs, &mol);
        else if(strstr(buffer, SCF_CALCULATION_PATTERN) != strstr(buffer, DFT_CALCULATION_PATTERN))
            ReadSinglePointCalculation(&ifs, &mol);
        else if(strstr(buffer, ZTS_CALCULATION_PATTERN) != NULL)
            ReadZTSCalculation(&ifs, &mol);
        else if(strstr(buffer, MEP_CALCULATION_PATTERN) != NULL)
            ReadMEPCalculation(&ifs, &mol);
        else if(strstr(buffer, NEB_CALCULATION_PATTERN) != NULL)
            ReadNEBCalculation(&ifs, &mol);
        // These calculation handlers still not implemented
        // so we just skip them
        else if(strstr(buffer, PROPERTY_CALCULATION_PATTERN) != NULL)
            GotoCalculationEnd(&ifs);
        else if (strstr(buffer, ESP_CALCULATION_PATTERN) != NULL)
            GotoCalculationEnd(&ifs);
        else if (strstr(buffer, PYTHON_CALCULATION_PATTERN) != NULL)
            GotoCalculationEnd(&ifs);
    }//while

    if (mol.NumAtoms() == 0) { // e.g., if we're at the end of a file PR#1737209
      mol.EndModify();
      return false;
    }

    if (!pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.ConnectTheDots();
    if (!pConv->IsOption("s",OBConversion::INOPTIONS) && !pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();

    mol.EndModify();
    // EndModify adds new conformer equals to current
    // molecule geometry so we will just delete it
    unsigned int nconformers = mol.NumConformers();
    if (nconformers > 1)
        mol.DeleteConformer(nconformers - 1);

    mol.SetTitle(title);
    return(true);
  }

  ////////////////////////////////////////////////////////////////

  bool NWChemInputFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    char buffer[BUFF_SIZE];

    ofs << "start molecule" << "\n\n";
    ofs << "title " << endl << " " << mol.GetTitle() << "\n\n";

    ofs << "geometry units angstroms print xyz autosym\n";

    FOR_ATOMS_OF_MOL(atom, mol)
      {
        snprintf(buffer, BUFF_SIZE, "%3s%15.5f%15.5f%15.5f\n",
                OBElements::GetSymbol(atom->GetAtomicNum()),
                atom->GetX(),
                atom->GetY(),
                atom->GetZ());
        ofs << buffer;
      }

    ofs << "end\n";

    return(true);
  }

} //namespace OpenBabel
