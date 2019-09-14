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
#include <openbabel/elements.h>


using namespace std;
namespace OpenBabel
{

  class CSRFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    CSRFormat()
    {
      OBConversion::RegisterFormat("csr",this);
    }

    virtual const char* Description() //required
    {
      return
        "Accelrys/MSI Quanta CSR format\n"
        "No comments yet\n";
    };

    virtual const char* SpecificationURL()
    {return "";}; //optional

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
      return NOTREADABLE;
    };

    //*** This section identical for most OBMol conversions ***
    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

  private:
    //	static bool FirstTime = true; Use new framework functions
    int MolCount; //was = 1;

    void WriteSize(int,ostream&);
    char *PadString(char*,int);
    void WriteCSRHeader(ostream&,OBMol&);
    void WriteCSRCoords(ostream&,OBMol&);


  };
  //***

  //Make an instance of the format class
  CSRFormat theCSRFormat;

  ////////////////////////////////////////////////////////////////

  bool CSRFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    //  if (FirstTime)
    if(pConv->GetOutputIndex()==1)
      {
        WriteCSRHeader(ofs,mol);
        //FirstTime = false;
        MolCount=1;
      }

    WriteCSRCoords(ofs,mol);
    MolCount++;

    return(true);
  }

  void CSRFormat::WriteCSRHeader(ostream &ofs,OBMol &mol)
  {
    char *molnames;
    int nmol, natom;

    molnames = PadString((char*)mol.GetTitle(),100);

    nmol = 1;
    natom = mol.NumAtoms();

    WriteSize(4*sizeof(char),ofs);
    ofs.write("V33 ",strlen("V33 ")*sizeof(char));
    WriteSize(4*sizeof(char),ofs);

    WriteSize(2*sizeof(int),ofs);
    ofs.write((char*)&natom,sizeof(int));
    ofs.write((char*)&nmol,sizeof(int));
    WriteSize(2*sizeof(int),ofs);

    WriteSize(100*sizeof(char),ofs);
    ofs.write(molnames,100*sizeof(char));
    WriteSize(100*sizeof(char),ofs);

    WriteSize(sizeof(int),ofs);
    ofs.write((char*)&natom,sizeof(int));
    WriteSize(sizeof(int),ofs);

    delete [] molnames;
  }

  void CSRFormat::WriteCSRCoords(ostream &ofs,OBMol &mol)
  {
    int the_size,jconf;
    double x,y,z,energy;
    char title[100];
    char *tag;

    the_size = sizeof(int) + sizeof(double) + (80 * sizeof(char));

    jconf = 1;
    energy = -2.584565;

    snprintf(title, 80, "%s:%d",mol.GetTitle(),MolCount);
    tag = PadString(title,80);

    WriteSize(the_size,ofs);
    ofs.write((char*)&jconf,sizeof(int));
    ofs.write((char*)&energy,sizeof(double));
    ofs.write(tag,80*sizeof(char));
    WriteSize(the_size,ofs);

    WriteSize(mol.NumAtoms()*sizeof(double),ofs);

    OBAtom *atom;
    vector<OBAtom*>::iterator i;
    for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
      {
        x = atom->x();
        ofs.write((char*)&x,sizeof(double));
      }
    WriteSize(mol.NumAtoms()*sizeof(double),ofs);

    WriteSize(mol.NumAtoms()*sizeof(double),ofs);
    for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
      {
        y = atom->y();
        ofs.write((char*)&y,sizeof(double));
      }
    WriteSize(mol.NumAtoms()*sizeof(double),ofs);

    WriteSize(mol.NumAtoms()*sizeof(double),ofs);
    for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
      {
        z = atom->z();
        ofs.write((char*)&z,sizeof(double));
      }
    WriteSize(mol.NumAtoms()*sizeof(double),ofs);

    delete [] tag;
  }

  void CSRFormat::WriteSize(int size,ostream &ofs)
  {
    ofs.write((char*)&size,sizeof(int));
  }

  char* CSRFormat::PadString(char *input, int size)
  {
    char *output;

    output = new char[size];
    memset(output, ' ', size);
    strncpy(output, input, strlen(input));
    output[ size - 1] = '\0';
    return(output);
  }

} //namespace OpenBabel
