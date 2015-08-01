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

//Required for double abs(double)
#include <cmath> 

#define HARTREE_TO_KCAL 627.509469

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
        "calculations with different molecules\n"
        "input in one output file detected\n"
        "(last calculation will be prefered)\n"
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
    vector<OBOrbital> ReadOrbitals(istream* ifs);

    OBVibrationData* ReadFrequencyCalculation(istream* ifs);
    void ReadGeometryOptimizationCalculation(istream* ifs, OBMol* molecule);
    void ReadSinglePointCalculation(istream* ifs, OBMol* molecule);

  };

static const char* COORDINATES_PATTERN = "Output coordinates";
static const char* GEOMETRY_OPTIMIZATION_PATTERN = "NWChem Geometry Optimization";
static const char* PROPERTY_CALCULATION_PATTERN = "NWChem Property Module";
static const char* ZTS_CALCULATION_PATTERN = "@ String method.";
static const char* SCF_CALCULATION_PATTERN = "SCF Module";
static const char* DFT_CALCULATION_PATTERN = "DFT Module";
static const char* SCF_ENERGY_PATTERN = "SCF energy =";
static const char* DFT_ENERGY_PATTERN = "DFT energy =";
static const char* FREQUENCY_PATTERN = "NWChem Nuclear Hessian and Frequency Analysis";
static const char* OPTIMIZATION_STEP_PATTERN = "Step       Energy";
static const char* VIBRATIONS_TABLE_PATTERN = "P.Frequency";
static const char* INTENSITIES_TABLE_PATTERN = "Projected Infra Red Intensities";
static const char* DIGITS = "1234567890";
static const char* END_OF_CALCULATION_PATTERN = "Task  times  cpu";
static const char* ORBITAL_START_PATTERN = "Vector";
static const char* ORBITAL_SECTION_PATTERN_1 = "Analysis";
static const char* ORBITAL_SECTION_PATTERN_2 = "rbital";
static const char* BETA_ORBITAL_PATTERN = "Beta";
static const char* MULLIKEN_CHARGES_PATTERN = "Mulliken analysis of the total density";

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
  Returns true if molecules contains same atoms regardless of
  their coordinates.
  */
  static bool CheckMoleculesEqual(OBMol* mol1, OBMol* mol2)
  {
    if (mol1->NumAtoms() != mol2->NumAtoms())
        return false;
    for(unsigned int i = 0; i < mol1->NumAtoms();i++)
        if(mol1->GetAtomById(i)->GetAtomicNum() != mol2->GetAtomById(i)->GetAtomicNum())
            return false;
    return true;
  }
  /////////////////////////////////////////////////////////////////
  /**
  Moves stream (ifs) position to end of calculation. 
  */
  static void GotoCalculationEnd(istream* ifs)
  {
  char buffer[BUFF_SIZE];
    while ( (strstr(buffer,END_OF_CALCULATION_PATTERN) == NULL))
        if (ifs->getline(buffer,BUFF_SIZE) == NULL)
            break;
  }

  /////////////////////////////////////////////////////////////////
  /**
  Adds non existent data from "source" OBMol object to "destination"
  and returns it then.
  If "source" object has incompatible data, NULL will be returned
  */
  static OBMol* AddNonExistentData(OBMol* destination, OBMol* source)
  {
    if ((destination->NumAtoms() == 0))
    {
        *destination += *source;
        return destination;
    }

    if (source->NumAtoms() != 0)
        if (!CheckMoleculesEqual(destination, source))
            return NULL;

    // Datas
    vector<OBGenericData* > datas = source->GetData();
    for(unsigned int i = 0;i < datas.size();i++)
    {
        destination->CloneData(datas[i]);
    }
    return destination;
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
                delete coordinates;
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
    if ((from_scratch)||(i != natoms))
        return;
    molecule->AddConformer(coordinates);
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
  and returns vector of OBOrbital objects.
  Input stream must be set to begining of orbital data
  section in nwo file. (Line after "... Molecular Orbital Analysis")
  Stream will be set at next line after end of orbital section.
  */
  vector<OBOrbital> NWChemOutputFormat::ReadOrbitals(istream* ifs)
  {
    vector<string> vs;
    char buffer[BUFF_SIZE];
    vector<OBOrbital> orbitals;

    ifs->getline(buffer, BUFF_SIZE); // ---------
    ifs->getline(buffer, BUFF_SIZE); // blank line

    ifs->getline(buffer,BUFF_SIZE);

    while (strstr(buffer, ORBITAL_START_PATTERN))
    {
        tokenize(vs, buffer);
        // Vector   N  Occ=X  E= Y  Symmetry=a'
        //   0      1    2    3  4  5(optional)
        if (vs.size() < 5)
            return orbitals;

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
        while (ifs->getline(buffer,BUFF_SIZE) != NULL)
            if (strlen(buffer) < 2) // If blank line detected
                break;
        ifs->getline(buffer,BUFF_SIZE); // New line that can contain "Vector"
    }

    return orbitals;
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

    while (ifs->getline(buffer, BUFF_SIZE) != NULL)
    {
        if(strstr(buffer,COORDINATES_PATTERN) != NULL)
        {
            ReadCoordinates(ifs, molecule);
            molecule->SetConformer(molecule->NumConformers() - 1);
        }
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
  Method reads vibration data from input stream (ifs)
  and writes it to new OBVibrationData object, then
  returns it.
  Input stream must be set to begining of frequency
  calculation in nwo file. (Line after "NWChem <theory> Module")
  If vibration data not found then NULL will be returned.
  Stream will be set at the end of calculation.
  */
  OBVibrationData* NWChemOutputFormat::ReadFrequencyCalculation(istream* ifs)
  {
    OBVibrationData* vibration_data = NULL;
    vector<double>  Frequencies, Intensities;
    vector<vector<vector3> > Lx;
    vector<string> vs;
    char buffer[BUFF_SIZE];

    while (ifs->getline(buffer, BUFF_SIZE) != NULL)
    {
        if (strstr(buffer, VIBRATIONS_TABLE_PATTERN) != NULL)
        {
            vector<double> freq;
            vector<vector<vector3> > vib;
            // freq and vib are auxiliary vectors which hold the data for
            // every block of 6 vibrations.
            tokenize(vs,buffer);
            for(unsigned int i=1; i<vs.size(); ++i)
                freq.push_back(atof(vs[i].c_str()));
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
              for (unsigned int i = 0; i < freq.size(); i++)
              {
                vib.push_back(vector<vector3>());
                vib[i].push_back(vector3(x[i], y[i], z[i]));
              }
              ifs->getline(buffer, BUFF_SIZE);
              tokenize(vs,buffer);
            }// while vs.size() > 2
            for (unsigned int i = 0; i < freq.size(); i++)
              {
                if (abs(freq[i]) > 10.0)
                {
                   // skip rotational and translational modes
                   Frequencies.push_back(freq[i]);
                   Lx.push_back(vib[i]);
                }// if abs(freq[i]) > 10.0
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
        else if(strstr(buffer, END_OF_CALCULATION_PATTERN) != NULL)
        {
            // End of task
            break;
        } // if "Task  times  cpu"
    }


    if (Frequencies.size() == 0)
        return NULL;

    vibration_data = new OBVibrationData;
    vibration_data->SetData(Lx, Frequencies, Intensities);

    return vibration_data;
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
    OBOrbitalData *orbital_data = NULL;

    while (ifs->getline(buffer, BUFF_SIZE))
    {
        if ((strstr(buffer, DFT_ENERGY_PATTERN) != NULL) || (strstr(buffer, SCF_ENERGY_PATTERN) != NULL))
        {
            tokenize(vs, buffer);
            energy = atof(vs[4].c_str()) * HARTREE_TO_KCAL;
        }
        else if ((strstr(buffer, ORBITAL_SECTION_PATTERN_2) != NULL)&&(strstr(buffer, ORBITAL_SECTION_PATTERN_1) != NULL))
        {
            if (orbital_data == NULL)
                orbital_data = new OBOrbitalData;

            vector<OBOrbital> orbitals = ReadOrbitals(ifs);

            if (strstr(buffer, BETA_ORBITAL_PATTERN))
            {
                orbital_data->SetOpenShell(true);
                orbital_data->SetBetaOrbitals(orbitals);
            }
            else
                orbital_data->SetAlphaOrbitals(orbitals);
        }
        else if(strstr(buffer, MULLIKEN_CHARGES_PATTERN) != NULL)
            ReadPartialCharges(ifs, molecule);
        else if (strstr(buffer, END_OF_CALCULATION_PATTERN) != NULL)
            break;
    }
    if (energy == 0)
    {
        if (orbital_data != NULL)
            delete orbital_data;
        return;
    }
    molecule->SetEnergy(energy);
    if (orbital_data != NULL)
        molecule->SetData(orbital_data);
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
    while	(ifs.getline(buffer,BUFF_SIZE) != NULL)
    {
        if(strstr(buffer,COORDINATES_PATTERN) != NULL)
        {
            // Input coordinates for calculation
            if ((mol.NumAtoms() == 0) | (pConv->IsOption("f",OBConversion::INOPTIONS) != NULL))
            {
                // If coordinates had redefined while calculation
                // in input file and "f" option had supplied then overwrite
                // all previous calculations. Otherwise calculations for
                // new geometry will be considered as new molecule.
                mol.Clear();
                mol.BeginModify();
                ReadCoordinates(&ifs, &mol);
            }
            else
            {
                unsigned int i;
                for(i=0; buffer[i] != '\0';i++);
                ifs.seekg(-i, ios_base::cur);
                break;
            }
        }
        else if(strstr(buffer, GEOMETRY_OPTIMIZATION_PATTERN) != NULL)
            ReadGeometryOptimizationCalculation(&ifs, &mol);
        else if(strstr(buffer, FREQUENCY_PATTERN) != NULL)
        {
            OBVibrationData* result = ReadFrequencyCalculation(&ifs);
            if (result != NULL)
                mol.SetData(result);
        }// if "Frequency Analysis"
        else if(strstr(buffer, SCF_CALCULATION_PATTERN) != strstr(buffer, DFT_CALCULATION_PATTERN))
            ReadSinglePointCalculation(&ifs, &mol);
        // These calculation handlers still not implemented
        // so we just skip them
        else if(strstr(buffer, ZTS_CALCULATION_PATTERN) != NULL)
            GotoCalculationEnd(&ifs); 
        else if(strstr(buffer, PROPERTY_CALCULATION_PATTERN) != NULL)
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
    mol.SetConformer(nconformers - 1);

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
                etab.GetSymbol(atom->GetAtomicNum()),
                atom->GetX(),
                atom->GetY(),
                atom->GetZ());
        ofs << buffer;
      }

    ofs << "end\n";

    return(true);
  }

} //namespace OpenBabel
