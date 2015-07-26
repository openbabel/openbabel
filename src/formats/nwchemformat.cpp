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
    OBVibrationData* ReadFrequencies(istream* ifs, OBVibrationData* vibration_data);

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
  /*
  Method reads vibration data from input stream (ifs) and writes it to existent
  OBVibrationData object or creates new OBVibrationData object when it does
  not exist.
  Input stream must be set to begining of frequency table in nwo file. (Line after
  "(Projected Frequencies expressed in cm-1)")
  Method returns true when vibration data read and written to OBVibrationData successful
  or returns false when supplied OBVibrationData object contain frequencies not related to
  read ones.
  */
  OBVibrationData* NWChemOutputFormat::ReadFrequencies(istream* ifs, OBVibrationData* vibration_data)
  {

    vector<double>  Frequencies;
    vector<vector<vector3> > Lx;
    vector<string> vs;
    char buffer[BUFF_SIZE];

    ifs->getline(buffer, BUFF_SIZE); // blank line
    ifs->getline(buffer, BUFF_SIZE); // frequency counting
    ifs->getline(buffer, BUFF_SIZE); // blank line
    ifs->getline(buffer, BUFF_SIZE); // P.Frequency  etc.
    while (strstr(buffer,"P.Frequency") != NULL)
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
        while(vs.size() > 2) {
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
          for (unsigned int i = 0; i < freq.size(); i++) {
            vib.push_back(vector<vector3>());
            vib[i].push_back(vector3(x[i], y[i], z[i]));
          }
          ifs->getline(buffer, BUFF_SIZE);
          tokenize(vs,buffer);
        } // while

        for (unsigned int i = 0; i < freq.size(); i++)
          {
            if (abs(freq[i]) > 10.0)
             {
               // skip rotational and translational modes
               Frequencies.push_back(freq[i]);
               Lx.push_back(vib[i]);
             }// if abs(freq[i]) > 10.0
          }// for

        ifs->getline(buffer, BUFF_SIZE); // frequency counting
        ifs->getline(buffer, BUFF_SIZE); // blank line
        ifs->getline(buffer, BUFF_SIZE); // P.Frequency  etc.
      }//while P.Frequency

    if (vibration_data == NULL)
        vibration_data = new OBVibrationData;
    else
      {
        //Check if read frequencies is not related to supplied OBVibrationData
        if (vibration_data->GetNumberOfFrequencies() != Frequencies.size())
            return NULL;
        vector<double> supplied_frequencies = vibration_data->GetFrequencies();
        for (unsigned int i = 0;i < supplied_frequencies.size();i++)
          {
            if (supplied_frequencies[i] != Frequencies[i])
                return NULL;
          }//for
      } // if !vibration_data

    vibration_data->SetData(Lx, Frequencies, vibration_data->GetIntensities());

    return vibration_data;
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

    //Vibrational data
    OBVibrationData* vibration_data = NULL;
    std::vector<double> Intensities;

    char buffer[BUFF_SIZE];
    string str;
    double x,y,z;
    OBAtom *atom;
    vector<string> vs;
    
    
    
    // Reading inital parameters of calculation such as
    // used theory and calculation type for better
    // recognition futher output
    string theory;
    string energy_pattern;
    CalculationType current_calculation = SinglePoint;
    while	(ifs.getline(buffer,BUFF_SIZE))
      {
        if(strstr(buffer,"NWChem SCF Module") != NULL)
          {
            theory = "SCF";
            break;
          }// if "SCF Module"
        else if(strstr(buffer,"NWChem DFT Module") != NULL)
          {
            theory = "DFT";
            break;
          }// if "DFT Module"
        else if(strstr(buffer,"NWChem Geometry Optimization") != NULL)
          {
            current_calculation = GeometryOptimization;
          }// if "Geometry Optimization"
        else if(strstr(buffer,"NWChem Nuclear Hessian and Frequency Analysis") != NULL)
          {
            current_calculation = VibrationalAnalysis;
          }// if "Frequency Analysis"
        else if(strstr(buffer,"@ String method.") != NULL)
          {
            current_calculation = ZTS;
          }// if "String method"
        else if(strstr(buffer,"NWChem Property Module") != NULL)
          {
            current_calculation = Property;
          }// if "Property"
          
      }//while
    
    if (!theory.empty())
    {
        stringstream ss;
        ss << "Total " << theory << " energy";
        energy_pattern = ss.str();
    }
    // "Task  times  cpu" End marker
    mol.BeginModify();
    while	(ifs.getline(buffer,BUFF_SIZE))
      {
        if(strstr(buffer,"Output coordinates") != NULL)
          {
            // mol.EndModify();
            OBMol* atoms = ReadCoordinates(&ifs);
            mol.Clear();
            mol.BeginModify();
            mol.ReserveAtoms(atoms->NumAtoms());
            for(unsigned int i = 0;i < atoms->NumAtoms();i++)
              {
                OBAtom* atom = atoms->GetAtomById(i);
                OBAtom* writeatom = mol.NewAtom();
                writeatom->SetVector(atom->GetVector());
                writeatom->SetAtomicNum(atom->GetAtomicNum());
              }
            delete atoms;
          } // if "output coordinates"
        if(strstr(buffer,"(Projected Frequencies expressed in cm-1)") != NULL)
          {
            vibration_data = ReadFrequencies(&ifs, vibration_data);
          } // if "P.Frequency"
        if(strstr(buffer,"Projected Infra Red Intensities") != NULL)
          {
           ifs.getline(buffer, BUFF_SIZE); // table header
           ifs.getline(buffer, BUFF_SIZE); // table delimiter
           ifs.getline(buffer, BUFF_SIZE);
           tokenize(vs,buffer);
           while (vs.size() == 7) {
             if (abs(atof(vs[1].c_str())) > 10.0)
                Intensities.push_back(atof(vs[5].c_str()));
             ifs.getline(buffer, BUFF_SIZE);
             tokenize(vs,buffer);
           }
         } // if "Projected Infra Red Intensities"
       if((!energy_pattern.empty()) && (strstr(buffer, energy_pattern.c_str()) != NULL))
         {
            tokenize(vs,buffer);
            mol.SetEnergy(atof(vs[3].c_str()) * HARTREE_TO_KCAL);
         }// if "Total energy"
      } // while

    if (mol.NumAtoms() == 0) { // e.g., if we're at the end of a file PR#1737209
      mol.EndModify();
      return false;
    }

    //Attach vibrational data, if there is any, to molecule
    if(vibration_data != NULL)
    {
      if ((Intensities.size() > 0) && (vibration_data->GetNumberOfFrequencies() == Intensities.size()))
        vibration_data->SetData(vibration_data->GetLx(),vibration_data->GetFrequencies(),Intensities);
      mol.SetData(vibration_data);
    }

    if (!pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.ConnectTheDots();
    if (!pConv->IsOption("s",OBConversion::INOPTIONS) && !pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();

    mol.EndModify();
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
