/**********************************************************************
Copyright (C) 2011 by Michael Banck

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

  class AcesOutputFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    AcesOutputFormat()
    {
      OBConversion::RegisterFormat("acesout",this);
    }

    virtual const char* Description() //required
    {
      return
        "ACES output format\n"
        "ACES is a set of programs that performs ab initio quantum chemistry calculations.\n"
        "Read Options e.g. -as\n"
        " s  Output single bonds only\n"
        " b  Disable bonding entirely\n\n";
    };

    virtual const char* SpecificationURL()
    {return "http://www.qtp.ufl.edu/ACES/";}; //optional

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
  AcesOutputFormat theAcesOutputFormat;

  class AcesInputFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    AcesInputFormat()
    {
      OBConversion::RegisterFormat("acesin",this);
    }

    virtual const char* Description() //required
    {
      return
        "ACES input format\n"
        "ACES is a set of programs that performs ab initio quantum chemistry calculations.\n"
        ;
    };

    virtual const char* SpecificationURL()
    {return "http://www.qtp.ufl.edu/ACES/";}; //optional

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
  AcesInputFormat theAcesInputFormat;


  /////////////////////////////////////////////////////////////////
  bool AcesOutputFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
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

    char buffer[BUFF_SIZE];
    string str;
    double x,y,z;
    OBAtom *atom;
    vector<string> vs;

    mol.BeginModify();
    while	(ifs.getline(buffer,BUFF_SIZE))
      {
        if(strstr(buffer,"Cartesian Coordinates") != NULL)
          {
            // mol.EndModify();
            mol.Clear();
            mol.BeginModify();
            ifs.getline(buffer,BUFF_SIZE);	// ---------------
            ifs.getline(buffer,BUFF_SIZE);	// blank
            ifs.getline(buffer,BUFF_SIZE);	// Total number of coordinates:
            ifs.getline(buffer,BUFF_SIZE);	// blank
            ifs.getline(buffer,BUFF_SIZE);	// blank
            ifs.getline(buffer,BUFF_SIZE);
            tokenize(vs,buffer);
            while (vs.size() == 5 || vs.size() == 6)
              {
                atom = mol.NewAtom();

                //set atomic number
                atom->SetAtomicNum(OBElements::GetAtomicNum(vs[1].c_str()));

		if (vs.size() == 6 ) {
                	x = atof((char*)vs[5].c_str());
		} else {
                	x = atof((char*)vs[4].c_str());
		}

            	ifs.getline(buffer,BUFF_SIZE);
            	tokenize(vs,buffer);
                y = atof((char*)vs[2].c_str());

		ifs.getline(buffer,BUFF_SIZE);
		tokenize(vs,buffer);
                z = atof((char*)vs[2].c_str());

                atom->SetVector(x,y,z); //set coordinates

		ifs.getline(buffer,BUFF_SIZE);	// blank

                if (!ifs.getline(buffer,BUFF_SIZE))
                  break;
                tokenize(vs,buffer);
              }
          } // if "Cartesian Coordinates"
        if(strstr(buffer,"Normal Coordinates") != NULL)
          {
            ifs.getline(buffer,BUFF_SIZE);     // [Dimensions are Mass**-1/2 Distance]
            ifs.getline(buffer,BUFF_SIZE);     // blank line
            ifs.getline(buffer,BUFF_SIZE);     // Symmetries
	    while((strstr(buffer,"Normal modes in internal coordinates") == NULL) && 
		  (strstr(buffer,"Dipole Moment Function") == NULL)) 
	      {
	        vector<vector3> vib1, vib2, vib3;
                ifs.getline(buffer,BUFF_SIZE);     // Frequencies
                tokenize(vs,buffer);
                for(unsigned int i = 0; i < vs.size(); i++) {
		    if (vs[i].find("i") != string::npos) {
		      // imaginary frequency
	              Frequencies.push_back(-atof(vs[i].c_str()));
		    } else {
	              Frequencies.push_back(atof(vs[i].c_str()));
		    }
	        }		
                ifs.getline(buffer,BUFF_SIZE);     // VIBRATION
                ifs.getline(buffer,BUFF_SIZE);     // X Y Z
                ifs.getline(buffer,BUFF_SIZE);	
                tokenize(vs,buffer);
	        while(vs.size() > 2) {
                  for (unsigned int i = 1; i < vs.size(); i += 3) {
	            double x, y, z;
                    x = atof(vs[i+0].c_str());
                    y = atof(vs[i+1].c_str());
                    z = atof(vs[i+2].c_str());
                    if (i == 1)
                      vib1.push_back(vector3(x, y, z));
                    else if (i == 4)
                      vib2.push_back(vector3(x, y, z));
                    else if (i == 7)
                      vib3.push_back(vector3(x, y, z));
	          }
                  ifs.getline(buffer, BUFF_SIZE);
                  tokenize(vs,buffer);
                } // while
                Lx.push_back(vib1);
	        if (vib2.size())
                  Lx.push_back(vib2);
	        if (vib3.size())
                  Lx.push_back(vib3);
	        ifs.getline(buffer,BUFF_SIZE);     // Symmetries (or end of frequencies)
	      } // while
          } // if "Normal Coordinates"
        if(strstr(buffer,"Infrared") != NULL)
          {
            ifs.getline(buffer, BUFF_SIZE); // table header
            ifs.getline(buffer, BUFF_SIZE); // table delimiter
            ifs.getline(buffer, BUFF_SIZE); // table header
            ifs.getline(buffer, BUFF_SIZE); // table delimiter
            ifs.getline(buffer, BUFF_SIZE);
            tokenize(vs,buffer);
            while (vs.size() == 4) {
              if (strstr(buffer,"VIBRATION") != NULL) {
                 Intensities.push_back(atof(vs[2].c_str()));
	      }	
              ifs.getline(buffer, BUFF_SIZE);
              tokenize(vs,buffer);
            }
          } // "Infrared"
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

  bool AcesInputFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    char buffer[BUFF_SIZE];

    ofs << mol.GetTitle() << "\n";

    FOR_ATOMS_OF_MOL(atom, mol)
      {
        snprintf(buffer, BUFF_SIZE, "%3s%15.5f%15.5f%15.5f\n",
                OBElements::GetSymbol(atom->GetAtomicNum()),
                atom->GetX(),
                atom->GetY(),
                atom->GetZ());
        ofs << buffer;
      }

    ofs << "\n*ACES2(__ADD_SETUP_HERE__)\n\n";

    return(true);
  }

} //namespace OpenBabel
