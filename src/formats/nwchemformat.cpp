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
    OBMol* ReadCoordinates(istream *ifs);
    OBVibrationData* ReadFrequencyCalculation(istream* ifs);
    OBMol* ReadGeometryOptimizationCalculation(istream* ifs);

    enum CalculationType
    {
        SinglePoint,
        GeometryOptimization,
        VibrationalAnalysis,
        ZTS,
        MEP,
        Property,
        Unknown
    };

  };

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

  //////////////////////////////////////////////////////
  /*
  Method reads coordinates from input stream (ifs) and creates
  new OBMol object then return reference to it.
  Input stream must be set to begining of coordinates
  table in nwo file. (Line after "Output coordinates...")
  */
  OBMol* NWChemOutputFormat::ReadCoordinates(istream *ifs)
  {
    OBMol* atoms;
    atoms = new OBMol;
    vector<string> vs;
    char buffer[BUFF_SIZE];
    double x, y, z;
    atoms->BeginModify();
    ifs->getline(buffer,BUFF_SIZE);	// blank
    ifs->getline(buffer,BUFF_SIZE);	// column headings
    ifs->getline(buffer,BUFF_SIZE);	// ---- ----- ----
    ifs->getline(buffer,BUFF_SIZE);
    tokenize(vs,buffer);
    for (unsigned int i=0; vs.size() == 6;i++)
      {
        OBAtom* atom = atoms->NewAtom();

        x = atof((char*)vs[3].c_str());
        y = atof((char*)vs[4].c_str());
        z = atof((char*)vs[5].c_str());
        atom->SetVector(x,y,z); //set coordinates

        //set atomic number
        size_t end_of_atom_symbol = vs[1].find_last_not_of("1234567890") + 1;
        atom->SetAtomicNum(etab.GetAtomicNum(vs[1].substr(0,end_of_atom_symbol).c_str()));
        if (!ifs->getline(buffer,BUFF_SIZE))
          break;
        tokenize(vs,buffer);
      }
    atoms->EndModify();
    return atoms;
  }

  //////////////////////////////////////////////////////
  /**
  Method reads optimization steps from input stream (ifs)
  and writes it to new OBMol object, then returns it.
  Input stream must be set to begining of geometry optimization
  calculation in nwo file. (Line after "NWChem Geometry Optimization")
  If no geometry data found then NULL will be returned.
  */
  OBMol* NWChemOutputFormat::ReadGeometryOptimizationCalculation(istream* ifs)
  {
    vector<double> energies;
    vector<string> vs;
    char buffer[BUFF_SIZE];

    OBMol *pmol = NULL;

    while (ifs->getline(buffer, BUFF_SIZE) != NULL)
    {
        if(strstr(buffer,"Output coordinates") != NULL)
        {
            OBMol* geometry = ReadCoordinates(ifs);
            unsigned int natoms = geometry->NumAtoms();
            double *conformer = new double[3*natoms];
            double *coordinates = geometry->GetCoordinates();
            memcpy(conformer, coordinates, sizeof(double)*3*natoms);

            if (pmol == NULL)
            {
                pmol = geometry;
                delete [] conformer;
            }
            else
            {
                //for(unsigned int i = 0;i < natoms;i++)
                //    pmol->GetAtomById(i)->SetVector(geometry->GetAtomById(i)->GetVector());
                delete geometry;
                pmol->AddConformer(conformer);
            }
        }
        else if(strstr(buffer,"Step       Energy") != NULL)
        {
            // Extract energy
            ifs->getline(buffer, BUFF_SIZE); // ------
            ifs->getline(buffer, BUFF_SIZE);
            tokenize(vs, buffer);
            energies.push_back(atof(vs[2].c_str()));
        }
        else if(strstr(buffer, "Task  times  cpu") != NULL)
        {
            // End of task
            break;
        } // if "Task  times  cpu"
    }
    for(unsigned int i = 0;i < pmol->NumAtoms();i++)
    {
        double *geometry = pmol->GetConformer(pmol->NumConformers() - 1);
        pmol->GetAtomById(i)->SetVector(geometry[i*3+0],
                                        geometry[i*3+1],
                                        geometry[i*3+2]);
    }
    pmol->SetEnergies(energies);
    return pmol;
  }

  //////////////////////////////////////////////////////
  /**
  Method reads vibration data from input stream (ifs)
  and writes it to new OBVibrationData object, then
  returns it.
  Input stream must be set to begining of frequency
  calculation in nwo file. (Line after "NWChem <theory> Module")
  If vibration data not found then NULL will be returned.
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
        if (strstr(buffer,"P.Frequency") != NULL)
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
        else if(strstr(buffer,"Projected Infra Red Intensities") != NULL)
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
        else if(strstr(buffer, "Task  times  cpu") != NULL)
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
  static bool CheckMoleculesEqual(OBMol* mol1, OBMol* mol2)
  {
    if (mol1->NumAtoms() != mol2->NumAtoms())
        return false;
    for(unsigned int i = 0; i < mol1->NumAtoms();i++)
        if(mol1->GetAtomById(i)->GetAtomicNum() != mol2->GetAtomById(i)->GetAtomicNum())
            return false;
    return true;
  }

  static void goto_calculation_end(istream *ifs)
  {
  char buffer[BUFF_SIZE];
    while ( (strstr(buffer,"Task  times  cpu") == NULL))
        if (ifs->getline(buffer,BUFF_SIZE) == NULL)
            break;
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
    string str;
    double x,y,z;
    OBAtom *atom;
    vector<string> vs;
    double last_calculated_energy;

    mol.BeginModify();

    // Reading inital parameters of calculation such as
    // used theory and calculation type for better
    // recognition futher output
    while	(ifs.getline(buffer,BUFF_SIZE) != NULL)
    {
        if(strstr(buffer,"Output coordinates") != NULL)
        {
            // Input coordinates for calculation
            if ((mol.NumAtoms() == 0) | (pConv->IsOption("f",OBConversion::INOPTIONS) != NULL))
            {
                OBMol* geometry = ReadCoordinates(&ifs);
                // If coordinates had redefined while calculation
                // in input file and "f" option had supplied then overwrite
                // all previous calculations. Otherwise calculations for
                // new geometry will be considered as new molecule.
                mol.Clear();
                mol += *geometry;
                delete geometry;
            }
            else
            {
                unsigned int i;
                for(i=0; buffer[i] != '\0';i++);
                ifs.seekg(-i, ios_base::cur);
                break;
            }
        }
        else if(strstr(buffer,"NWChem Geometry Optimization") != NULL)
        {
            OBMol* result = ReadGeometryOptimizationCalculation(&ifs);
            if (result != NULL)
            {
                if (mol.NumAtoms() == 0)
                {
                    mol += *result;
                    delete result;
                    continue;
                }
                else if (!CheckMoleculesEqual(&mol, result))
                {
                    delete result;
                    continue;
                }
                unsigned int natoms = result->NumAtoms();
                vector<double*> conformers;
                conformers.reserve(result->NumConformers());
                for(unsigned int i = 0;i < result->NumConformers();i++)
                {
                    double * conformer = new double[natoms*3];
                    memcpy(conformer, result->GetConformer(i),sizeof(double)*3*natoms);
                    mol.AddConformer(conformer);
                }
                vector<double> energies = result->GetEnergies();
                mol.SetEnergies(energies);
                last_calculated_energy = energies[energies.size() -1];
                mol.SetConformer(mol.NumConformers() - 1);
                for(unsigned int i = 0;i < natoms;i++)
                    mol.GetAtomById(i)->SetVector(result->GetAtomById(i)->GetVector());
                delete result;
            }
        }// if "Geometry Optimization"
        else if(strstr(buffer,"NWChem Nuclear Hessian and Frequency Analysis") != NULL)
        {
            OBVibrationData* result = ReadFrequencyCalculation(&ifs);
            if (result != NULL)
                mol.SetData(result);
        }// if "Frequency Analysis"
        // These calculation handlers still not implemented
        // so we just skip them
        else if(strstr(buffer,"@ String method.") != NULL)
            goto_calculation_end(&ifs); 
        else if(strstr(buffer,"NWChem Property Module") != NULL)
            goto_calculation_end(&ifs);
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
    // EndModify adds new conformer equals to current molecule geometry
    // so we will just assign last calculated energy to this conformer
    mol.SetConformer(mol.NumConformers() - 1);
    mol.SetEnergy(last_calculated_energy);

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
