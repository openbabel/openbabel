/**********************************************************************
Copyright (C) 2009 by Michael Banck

Based on nwchemformat.cpp,
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
#include <cstdlib>


using namespace std;
namespace OpenBabel
{
#define BOHR_TO_ANGSTROM 0.529177249

  class MolproOutputFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    MolproOutputFormat()
    {
      OBConversion::RegisterFormat("mpo",this);
    }

    virtual const char* Description() //required
    {
      return
        "Molpro output format\n"
        "Read Options e.g. -as\n"
        " s  Output single bonds only\n"
        " b  Disable bonding entirely\n\n";
    };

//TODO    virtual const char* SpecificationURL()
//    {return "http://www.example.com";}; //optional

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
      return READONEONLY | NOTWRITABLE;
    };

    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
  };

  //Make an instance of the format class
  MolproOutputFormat theMolproOutputFormat;

  class MolproInputFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    MolproInputFormat()
    {
      OBConversion::RegisterFormat("mp",this);
    }

    virtual const char* Description() //required
    {
      return
        "Molpro input format\n"
        "No comments yet\n";
    };

//TODO    virtual const char* SpecificationURL()
//    {return "http://www.example.com";}; //optional

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
  MolproInputFormat theMolproInputFormat;


  /////////////////////////////////////////////////////////////////
  bool MolproOutputFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title = pConv->GetTitle();

    //Vibrational data
    std::vector< std::vector< vector3 > > Lx;
    std::vector<double> Frequencies, Intensities;

    // 0 for no vibrations, 1 for regular vibrations, 2 for imaginary
    // vibrations, 3 for low/zero vibrations
    int vibration_state = 0;

    char buffer[BUFF_SIZE];
    string str;
    double x,y,z;
    OBAtom *atom;
    vector<string> vs;

    mol.BeginModify();
    while	(ifs.getline(buffer,BUFF_SIZE))
      {
        if(strstr(buffer,"ATOMIC COORDINATES") != NULL)
          {
            // mol.EndModify();
            mol.Clear();
            mol.BeginModify();
            ifs.getline(buffer,BUFF_SIZE);	// blank
            ifs.getline(buffer,BUFF_SIZE);	// column headings
            ifs.getline(buffer,BUFF_SIZE);	// blank
            ifs.getline(buffer,BUFF_SIZE);
            tokenize(vs,buffer);
            while (vs.size() == 6)
              {
                atom = mol.NewAtom();
                x = atof((char*)vs[3].c_str())*BOHR_TO_ANGSTROM;
                y = atof((char*)vs[4].c_str())*BOHR_TO_ANGSTROM;
                z = atof((char*)vs[5].c_str())*BOHR_TO_ANGSTROM;
                atom->SetVector(x,y,z); //set coordinates

                //set atomic number
                int n;
                atom->SetAtomicNum(0);
                while (vs[1].length()!=0) { // recognize name with number
                    n = OBElements::GetAtomicNum(vs[1].c_str());
                  if (n!=0) {
                    atom->SetAtomicNum(n);
                    break;
                  } else {
                    vs[1].erase(vs[1].end()-1,vs[1].end());
                  }
                }

                if (!ifs.getline(buffer,BUFF_SIZE))
                  break;
                tokenize(vs,buffer);
              }
          } // if "ATOMIC COORDINATES"
        if(strstr(buffer,"Normal Modes") != NULL && strstr(buffer,"of") == NULL) {
          vibration_state = 1;
          continue;
        }
        if(strstr(buffer,"Normal Modes of imag") != NULL) {
          vibration_state = 2;
          continue;
        }
        if(strstr(buffer,"Normal Modes of low") != NULL) {
          vibration_state = 3;
          continue;
        }
        if(strstr(buffer,"Wavenumbers [cm-1]") != NULL && vibration_state < 3)
          {
            // freq, intens and vib are auxiliary vectors which hold the data
            // for every block of 5 vibrations.
            vector<double> freq;
            vector<double> intens;
            vector<vector<vector3> > vib;

            tokenize(vs,buffer);
            for(unsigned int i=2; i<vs.size(); ++i) {
                if (vibration_state == 2) {
                  // imaginary frequency
                  freq.push_back(-atof(vs[i].c_str()));
                } else {
                  freq.push_back(atof(vs[i].c_str()));
                }
	    }

            ifs.getline(buffer,BUFF_SIZE);
            tokenize(vs,buffer);
            for(unsigned int i=2; i<vs.size(); ++i) {
                intens.push_back(atof(vs[i].c_str()));
	    }
            ifs.getline(buffer,BUFF_SIZE); // relative intensities
            ifs.getline(buffer,BUFF_SIZE);
            tokenize(vs,buffer);
	    while(vs.size() > 1) {
              vector<double> x, y, z;
              for (unsigned int i = 1; i < vs.size(); i++)
                x.push_back(atof(vs[i].c_str()));
              ifs.getline(buffer, BUFF_SIZE);
              tokenize(vs,buffer);
              for (unsigned int i = 1; i < vs.size(); i++)
                y.push_back(atof(vs[i].c_str()));
              ifs.getline(buffer, BUFF_SIZE);
              tokenize(vs,buffer);
              for (unsigned int i = 1; i < vs.size(); i++)
                z.push_back(atof(vs[i].c_str()));
              for (unsigned int i = 0; i < freq.size(); i++) {
                vib.push_back(vector<vector3>());
                vib[i].push_back(vector3(x[i], y[i], z[i]));
              }
              ifs.getline(buffer, BUFF_SIZE);
              tokenize(vs,buffer);
            } // while
            for (unsigned int i = 0; i < freq.size(); i++) {
	      Frequencies.push_back(freq[i]);
              Intensities.push_back(intens[i]);
              Lx.push_back(vib[i]);
            }
          } // if "Normal Modes"
        if(strstr(buffer,"STATE") != NULL && strstr(buffer,"DIPOLE MOMENT") != NULL)
          {
            tokenize(vs,buffer);
            if (vs.size() == 8) {
              OBVectorData *dipoleMoment = new OBVectorData;
              dipoleMoment->SetAttribute("Dipole Moment");
              double x, y, z;
              x = atof(vs[5].c_str());
              y = atof(vs[6].c_str());
              z = atof(vs[7].c_str());
              dipoleMoment->SetData(x, y, z);
              dipoleMoment->SetOrigin(fileformatInput);
              mol.SetData(dipoleMoment);
            }
          } // "STATE" && "DIPOLE MOMENT"

        if(strstr(buffer, "   HF-SCF") || strstr(buffer, "   RHF-SCF") || strstr(buffer, "   OPTG(RHF)"))
          {
            ifs.getline(buffer,BUFF_SIZE);
            const double myHARTREE_TO_KCAL = 627.509;
            mol.SetEnergy(atof(buffer) * myHARTREE_TO_KCAL); //Uses first value on line
          }
      } // while

    if (mol.NumAtoms() == 0) { // e.g., if we're at the end of a file PR#1737209
      mol.EndModify();
      return false;
    }

    //Attach vibrational data, if there is any, to molecule
    if(Frequencies.size()>0)
    {
      OBVibrationData* vd = new OBVibrationData;
      vd->SetData(Lx, Frequencies, Intensities);
      mol.SetData(vd);
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

  bool MolproInputFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    char buffer[BUFF_SIZE];

    ofs << "*** " << mol.GetTitle() << endl;
    ofs << "!file,2,INSERT WAVEFUNCTION FILE LOCATION HERE" << endl;
    ofs << "!memory,INSERT MEMORY HERE" << endl;
    ofs << "!basis,INSERT BASIS SET HERE" << endl;
    ofs << "\n";
    ofs << "geomtyp=xyz" << endl;
    ofs << "geometry={" << endl;
    ofs << mol.NumAtoms() << endl;
    ofs << "Geometry specification:" << endl;

    FOR_ATOMS_OF_MOL(atom, mol)
      {
        snprintf(buffer, BUFF_SIZE, "%3s,%15.5f,%15.5f,%15.5f\n",
                OBElements::GetSymbol(atom->GetAtomicNum()),
                atom->GetX(),
                atom->GetY(),
                atom->GetZ());
        ofs << buffer;
      }

    ofs << "}\n\n";

    ofs << "!INSERT QM METHODS HERE" << endl;
    ofs << "!hf" << endl;
    ofs << "---" << endl;
    return(true);
  }

} //namespace OpenBabel
