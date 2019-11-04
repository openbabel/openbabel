/**********************************************************************
Copyright (C) 2011 by Geoffrey Hutchison

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

  class XSFFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    XSFFormat()
    {
      OBConversion::RegisterFormat("xsf",this);
      // animation variant
      OBConversion::RegisterFormat("axsf",this);
    }

    virtual const char* Description() //required
    {
      return
        "XCrySDen Structure Format\n"
        "Read Options e.g. -as\n"
        "  s  Output single bonds only\n"
        "  b  Disable bonding entirely\n\n";
    };

    virtual const char* SpecificationURL()
    { return "http://www.xcrysden.org/doc/XSF.html/" ;}; //optional

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
      return READONEONLY | NOTWRITABLE;
    };

    //*** This section identical for most OBMol conversions ***
    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
  };
  //***

  //Make an instance of the format class
  XSFFormat theXSFFormat;

  /////////////////////////////////////////////////////////////////
  bool XSFFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
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
    vector3 translationVectors[3];
    int numTranslationVectors = 0;
    vector<string> vs;
    vector<vector3> atomPositions;
    bool createdAtoms = false;
    int atomicNum;

    mol.BeginModify();

    while (ifs.getline(buffer, BUFF_SIZE))
      {
        if (buffer[0] == '#')
          continue; // comment
        if (strstr(buffer, "ATOMS") != NULL) {
          // Minimum of 4 columns -- AtNum, x, y, z (forces)
          // where AtNum stands for atomic number (or symbol), while X Y Z are
          ifs.getline(buffer, BUFF_SIZE);
          tokenize(vs, buffer);
          while (vs.size() >= 4) {
            if (!createdAtoms) {
              atom = mol.NewAtom();
              //set atomic number
              atomicNum = OBElements::GetAtomicNum(vs[0].c_str());
              if (atomicNum == 0) {
                atomicNum = atoi(vs[0].c_str());
              }
              atom->SetAtomicNum(atomicNum);
            }
            x = atof((char*)vs[1].c_str());
            y = atof((char*)vs[2].c_str());
            z = atof((char*)vs[3].c_str());
            atomPositions.push_back(vector3(x, y, z)); // we may have a movie or animation

            ifs.getline(buffer, BUFF_SIZE);
            tokenize(vs, buffer);
          }
          createdAtoms = true; // don't run NewAtom() anymore
        }
        else if ( strstr(buffer, "PRIMVEC")
                 || strstr(buffer, "CONVVEC") ) {
          // translation vectors
          numTranslationVectors = 0; // if we have an animation
          while (numTranslationVectors < 3 && ifs.getline(buffer,BUFF_SIZE)) {
            tokenize(vs,buffer); // we really need to check that it's 3 entries only
            if (vs.size() < 3) return false; // timvdm 18/06/2008
            x = atof((char*)vs[0].c_str());
            y = atof((char*)vs[1].c_str());
            z = atof((char*)vs[2].c_str());
            translationVectors[numTranslationVectors++].Set(x, y, z);
          }
        }
        else if (strstr(buffer, "PRIMCOORD") != NULL) {
          // read the coordinates
          ifs.getline(buffer, BUFF_SIZE);
          tokenize(vs, buffer);
          if (vs.size() < 2) return false;
          int numAtoms = atoi(vs[0].c_str());
          for (int a = 0; a < numAtoms; ++a) {
            if (!ifs.getline(buffer,BUFF_SIZE))
              break;
            tokenize(vs,buffer);
            if (vs.size() < 4)
              break;

            if (!createdAtoms) {
              atom = mol.NewAtom();
              //set atomic number
              atomicNum = OBElements::GetAtomicNum(vs[0].c_str());
              if (atomicNum == 0) {
                atomicNum = atoi(vs[0].c_str());
              }
              atom->SetAtomicNum(atomicNum);
            }
            x = atof((char*)vs[1].c_str());
            y = atof((char*)vs[2].c_str());
            z = atof((char*)vs[3].c_str());
            atomPositions.push_back(vector3(x, y, z));
          }
        }
      }

    mol.EndModify();

    int natom = mol.NumAtoms();
    int numConformers = atomPositions.size() / natom;
    for (int i = 0; i < numConformers; ++i) {
      double *coordinates = new double[natom * 3];
      for (int j = 0; j < natom; ++j) {
        vector3 currentPosition = atomPositions[i*natom + j];
        coordinates[j*3] = currentPosition.x();
        coordinates[j*3 + 1] = currentPosition.y();
        coordinates[j*3 + 2] = currentPosition.z();
      }
      mol.AddConformer(coordinates);
    }
    // Delete first conformer, created by EndModify, bunch of 0s
    mol.DeleteConformer(0);
    // Set geometry to last one
    mol.SetConformer(mol.NumConformers() - 1);

    if (!pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.ConnectTheDots();
    if (!pConv->IsOption("s",OBConversion::INOPTIONS) && !pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();

    // Add final properties
    mol.SetTitle(title);
    if (numTranslationVectors == 3) {
      OBUnitCell *uc = new OBUnitCell;
      uc->SetOrigin(fileformatInput);
      uc->SetData(translationVectors[0],
                  translationVectors[1],
                  translationVectors[2]);
      mol.SetData(uc);
    }

    return(true);
  }

} //namespace OpenBabel
