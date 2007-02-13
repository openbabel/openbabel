/**********************************************************************
Copyright (C) 2000 by OpenEye Scientific Software, Inc.
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

using namespace std;
namespace OpenBabel
{

  class GaussianOutputFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    GaussianOutputFormat()
    {
      OBConversion::RegisterFormat("gal",this, "chemical/x-gaussian-log");
      OBConversion::RegisterFormat("g98",this);
      OBConversion::RegisterFormat("g03",this);
    }

    virtual const char* Description() //required
    {
      return
        "Gaussian98/03 Output\n \
       Read Options e.g. -as\n\
        s  Output single bonds only\n\
        b  Disable bonding entirely\n\n";
    };

    virtual const char* SpecificationURL()
    { return "http://www.gaussian.com/";};

    virtual const char* GetMIMEType() 
    { return "chemical/x-gaussian-log"; };

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
      return READONEONLY | NOTWRITABLE;
    };

    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
  };

  //Make an instance of the format class
  GaussianOutputFormat theGaussianOutputFormat;

  class GaussianInputFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    GaussianInputFormat()
    {
      OBConversion::RegisterFormat("com",this, "chemical/x-gaussian-input");
      OBConversion::RegisterFormat("gau",this);
      OBConversion::RegisterFormat("gjc",this);
      OBConversion::RegisterFormat("gjf",this);
      // Command-line keywords
      OBConversion::RegisterOptionParam("k", NULL, 1, OBConversion::OUTOPTIONS);
      // Command-line keyword file
      OBConversion::RegisterOptionParam("f", NULL, 1, OBConversion::OUTOPTIONS);    }

    virtual const char* Description() //required
    {
      return
        "Gaussian 98/03 Cartesian Input\n \
       Write Options e.g. -xk\n\
        k  \"keywords\" Use the specified keywords for input\n\
        f    <file>     Read the file specified for input keywords\n\n";
    };

    virtual const char* SpecificationURL()
    {return "http://www.gaussian.com/g_ur/m_input.htm";};

    virtual const char* GetMIMEType() 
    { return "chemical/x-gaussian-input"; };

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
  GaussianInputFormat theGaussianInputFormat;

  ////////////////////////////////////////////////////////////////

  bool GaussianInputFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    char buffer[BUFF_SIZE];
    const char *keywords = pConv->IsOption("k",OBConversion::OUTOPTIONS);
    const char *keywordsEnable = pConv->IsOption("k",OBConversion::GENOPTIONS);
    const char *keywordFile = pConv->IsOption("f",OBConversion::OUTOPTIONS);
    string defaultKeywords = "#Put Keywords Here, check Charge and Multiplicity.";

    if(keywords)
    {
      defaultKeywords = keywords;
    }
    
    if (keywordsEnable)
    {
      string model;
      string basis;
      string method;

      OBPairData *pd = (OBPairData *) pmol->GetData("model");
      if(pd)
        model = pd->GetValue();

      pd = (OBPairData *) pmol->GetData("basis");
      if(pd)
        basis = pd->GetValue();

      pd = (OBPairData *) pmol->GetData("method");
      if(pd)
        method = pd->GetValue();

      if(method == "optimize")
      {
        method = "opt";
      }

      if(model != "" && basis != "" && method != "")
      {
        ofs << model << "/" << basis << "," << method << endl;
      }
      else
      {
        ofs << "#Unable to translate keywords!" << endl;
        ofs << defaultKeywords << endl;
      }
    }
    else if (keywordFile)
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
    {
      ofs << defaultKeywords << endl;
    }
    ofs << endl; // blank line after keywords
    ofs << " " << mol.GetTitle() << endl << endl;

    snprintf(buffer, BUFF_SIZE, "%d  %d", mol.GetTotalCharge(), mol.GetTotalSpinMultiplicity());
    ofs << buffer << endl;
    FOR_ATOMS_OF_MOL(atom, mol)
      {
        if (atom->GetIsotope() == 0)
          snprintf(buffer, BUFF_SIZE, "%-3s      %10.5f      %10.5f      %10.5f ",
                  etab.GetSymbol(atom->GetAtomicNum()),
                  atom->GetX(), atom->GetY(), atom->GetZ());
        else
          snprintf(buffer, BUFF_SIZE, "%-3s(Iso=%d) %10.5f      %10.5f      %10.5f ",
                  etab.GetSymbol(atom->GetAtomicNum()),
                  atom->GetIsotope(),
                  atom->GetX(), atom->GetY(), atom->GetZ());
	
        ofs << buffer << endl;
      }
    // file should end with a blank line
    ofs << endl;
    return(true);
  }

  // Reading Gaussian output has been tested for G98 and G03 to some degree
  // If you have problems (or examples of older output), please contact
  // the openbabel-discuss@lists.sourceforge.net mailing list and/or post a bug
  bool GaussianOutputFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title = pConv->GetTitle();

    char buffer[BUFF_SIZE];
    string str,str1;
    double x,y,z;
    OBAtom *atom;
    vector<string> vs;
    int charge = 0;
    unsigned int spin = 1;
    bool hasPartialCharges = false;
    
    mol.BeginModify();
    
    while (ifs.getline(buffer,BUFF_SIZE))
      {
        if (strstr(buffer,"Symbolic Z-matrix:") != NULL)
          {
            ifs.getline(buffer,BUFF_SIZE); // Charge = # Multiplicty = #
            tokenize(vs, buffer, " \t\n");
            if (vs.size() == 6)
              {
                charge = atoi(vs[2].c_str());
                spin = atoi(vs[5].c_str());
              }
	    
            ifs.getline(buffer,BUFF_SIZE);
          }
        else if(strstr(buffer,"Coordinates (Angstroms)") != NULL)
          {
            // mol.EndModify();
            mol.Clear();
            mol.BeginModify();
            ifs.getline(buffer,BUFF_SIZE);	// column headings
            ifs.getline(buffer,BUFF_SIZE);	// ---------------
            ifs.getline(buffer,BUFF_SIZE);
            tokenize(vs,buffer);
            while (vs.size() == 6)
              {
                int atomicNum = atoi((char*)vs[1].c_str());
                if (atomicNum > 0) // translation vectors are "-2"
                  {
                    atom = mol.NewAtom();
                    atom->SetAtomicNum(atoi((char*)vs[1].c_str()));
                    x = atof((char*)vs[3].c_str());
                    y = atof((char*)vs[4].c_str());
                    z = atof((char*)vs[5].c_str());
                    atom->SetVector(x,y,z);
                  }
		
                if (!ifs.getline(buffer,BUFF_SIZE)) break;
                tokenize(vs,buffer);
              }
          }
        else if(strstr(buffer,"Total atomic charges") != NULL ||
                strstr(buffer,"Mulliken atomic charges") != NULL)
          {
            hasPartialCharges = true;
            ifs.getline(buffer,BUFF_SIZE);	// column headings
            ifs.getline(buffer,BUFF_SIZE);
            tokenize(vs,buffer);
            while (vs.size() >= 3 && 
                   strstr(buffer,"Sum of ") == NULL)
              {
                atom = mol.GetAtom(atoi(vs[0].c_str()));
                atom->SetPartialCharge(atof(vs[2].c_str()));
		
                if (!ifs.getline(buffer,BUFF_SIZE)) break;
                tokenize(vs,buffer);
              }
          }
        else if(strstr(buffer,"SCF Done:") != NULL)
          {
#define HARTREE_TO_KCAL 627.509
            tokenize(vs,buffer);
            mol.SetEnergy(atof(vs[4].c_str()) * HARTREE_TO_KCAL);
          }
      } // end while

    if (!pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.ConnectTheDots();
    if (!pConv->IsOption("s",OBConversion::INOPTIONS) && !pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();

    mol.EndModify();
    if (hasPartialCharges)
      mol.SetPartialChargesPerceived();
    mol.SetTotalCharge(charge);
    mol.SetTotalSpinMultiplicity(spin);
    
    mol.SetTitle(title);
    return(true);
  }
    
} //namespace OpenBabel
