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
        if(strstr(buffer,"Output coordinates") != NULL)
          {
            // mol.EndModify();
            mol.Clear();
            mol.BeginModify();
            ifs.getline(buffer,BUFF_SIZE);	// blank
            ifs.getline(buffer,BUFF_SIZE);	// column headings
            ifs.getline(buffer,BUFF_SIZE);	// ---- ----- ----
            ifs.getline(buffer,BUFF_SIZE);
            tokenize(vs,buffer);
            while (vs.size() == 6)
              {
                atom = mol.NewAtom();
                x = atof((char*)vs[3].c_str());
                y = atof((char*)vs[4].c_str());
                z = atof((char*)vs[5].c_str());
                atom->SetVector(x,y,z); //set coordinates

                //set atomic number
                atom->SetAtomicNum(etab.GetAtomicNum(vs[1].c_str()));

                if (!ifs.getline(buffer,BUFF_SIZE))
                  break;
                tokenize(vs,buffer);
              }
          } // if "output coordinates"
        if(strstr(buffer,"P.Frequency") != NULL)
          {
            // freq and vib are auxiliary vectors which hold the data for
            // every block of 6 vibrations.
            vector<double> freq;
            vector<vector<vector3> > vib;
            tokenize(vs,buffer);
            for(unsigned int i=1; i<vs.size(); ++i)
                freq.push_back(atof(vs[i].c_str()));
            ifs.getline(buffer,BUFF_SIZE);     // blank line
            ifs.getline(buffer,BUFF_SIZE);
            tokenize(vs,buffer);
	    while(vs.size() > 2) {
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
              if (freq[i] > 10.0) {
                // skip rotational and translational modes
	        Frequencies.push_back(freq[i]);
                Lx.push_back(vib[i]);
              }
            }
         } // if "P.Frequency"
       if(strstr(buffer,"Projected Infra Red Intensities") != NULL)
         {
           ifs.getline(buffer, BUFF_SIZE); // table header
           ifs.getline(buffer, BUFF_SIZE); // table delimiter
           ifs.getline(buffer, BUFF_SIZE);
           tokenize(vs,buffer);
           while (vs.size() == 7) {
             if (atof(vs[1].c_str()) > 10.0)
                Intensities.push_back(atof(vs[5].c_str()));
             ifs.getline(buffer, BUFF_SIZE);
             tokenize(vs,buffer);
           }
         } // if "Projected Infra Red Intensities"
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
