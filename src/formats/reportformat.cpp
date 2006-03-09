/**********************************************************************
reportformat.cpp - Report information about the molecule: charge, distance 
             matrix angle, chiral info, etc.
 
Copyright (C) 2000 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison
Some portions Copyright (C) 2004 by Chris Morley
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#include "babelconfig.h"

#include "obmolecformat.h"
#include "obconversion.h"

using namespace std;
namespace OpenBabel
{

class ReportFormat : public OBMoleculeFormat
{
public:
    //Register this format type ID
    ReportFormat()
    {
        OBConversion::RegisterFormat("report",this);
    }

    virtual const char* Description() //required
    {
        return
            "Open Babel report format\n \
            No comments yet\n";
    };

  virtual const char* SpecificationURL()
  {return "";}; //optional

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
        return NOTREADABLE;
    };

    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

    static void WriteCharges(ostream &ofs,OBMol &mol);
    static void WriteDistanceMatrix(ostream &ofs,OBMol &mol);
    static void WriteTorsions(ostream &ofs,OBMol &mol);
    static void WriteAngles(ostream &ofs,OBMol &mol);
    static void WriteChiral(ostream &ofs,OBMol &mol);

};

//Make an instance of the format class
ReportFormat theReportFormat;

////////////////////////////////////////////////////////////////

bool ReportFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
        return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    char buffer[BUFF_SIZE];
    ofs << "FILENAME: " << mol.GetTitle() << endl;
    ofs << "FORMULA: " << mol.GetFormula() << endl;
    ofs << "MASS: ";
    sprintf(buffer, "%5.4f", mol.GetMolWt());
    ofs << buffer << endl;
    ofs << "EXACT MASS: ";
    sprintf(buffer, "%5.7f", mol.GetExactMass());
    ofs << buffer << endl;
    if (mol.GetTotalCharge() != 0)
    {
        ofs << "TOTAL CHARGE: ";
        sprintf(buffer, "%d", mol.GetTotalCharge());
        ofs << buffer << endl;
    }
    if (mol.GetTotalSpinMultiplicity() != 1)
    {
        ofs << "TOTAL SPIN: ";
        sprintf(buffer, "%d", mol.GetTotalSpinMultiplicity());
        ofs << buffer << endl;
    }
    ofs << "INTERATOMIC DISTANCES" << endl;
    WriteDistanceMatrix(ofs, mol);
    ofs << endl << endl << "ATOMIC CHARGES" << endl;
    WriteCharges(ofs, mol);
    ofs << endl << endl << "BOND ANGLES" << endl;
    WriteAngles(ofs, mol);
    ofs << endl << endl << "TORSION ANGLES" << endl;
    WriteTorsions(ofs, mol);
    if (mol.IsChiral())
    {
        ofs << endl << endl << "CHIRAL ATOMS" << endl;
        WriteChiral(ofs, mol);
    }
    if (mol.HasData(OBGenericDataType::CommentData))
    {
        ofs << endl << endl << "COMMENTS" << endl;
        OBCommentData *cd = (OBCommentData*)mol.GetData(OBGenericDataType::CommentData);
        ofs << cd->GetData() << endl;
    }
    ofs << endl << endl;
    return(true);
}

/////////////////////////////////////////////////////////////
void ReportFormat::WriteCharges(ostream &ofs,OBMol &mol)
{
    unsigned int i;
    OBAtom *atom;
    char buffer[BUFF_SIZE];

    for(i = 1;i <= mol.NumAtoms(); i++)
    {
        atom = mol.GetAtom(i);
        sprintf(buffer,"%4s%4d   % 2.10f",
                etab.GetSymbol(atom->GetAtomicNum()),
                i,
                atom->GetPartialCharge());

        ofs << buffer << endl;
    }
}

void ReportFormat::WriteDistanceMatrix(ostream &ofs,OBMol &mol)
{
    int columns = 7;
    unsigned int max, min = 1;
    unsigned int i,j;
    string type;
    OBAtom *atom, *atom2;
    char buffer[BUFF_SIZE];
    double dst;

    max = columns;
    while (max <= mol.NumAtoms() + columns)
    {
        ofs << endl;
        if (min > mol.NumAtoms())
            break;
        atom = mol.GetAtom(min);

        sprintf(buffer,"%15s%4d",
                etab.GetSymbol(atom->GetAtomicNum()),
                min);
        ofs << buffer;

        for (i = min + 1; ((i < max) && (i <= mol.NumAtoms())); i++)
            if (i <= mol.NumAtoms())
            {
                atom = mol.GetAtom(i);
                sprintf(buffer,"%7s%4d",
                        etab.GetSymbol(atom->GetAtomicNum()),
                        i);
                ofs << buffer;
            }
        ofs << endl;

        sprintf(buffer,"%14s","");
        ofs << buffer;
        for (i = min; i < max; i++)
            if (i <= mol.NumAtoms())
            {
                sprintf(buffer,"%11s","-----------");
                ofs << buffer;
            }

        ofs << endl;
        for (i = min; i <= mol.NumAtoms(); i++)
        {
            atom = mol.GetAtom(i);
            sprintf(buffer,"%4s%4d",
                    etab.GetSymbol(atom->GetAtomicNum()),
                    i);
            ofs << buffer;
            for (j = min; j < max; j++)
                if (j <= i)
                {
                    atom2 = mol.GetAtom(j);
                    dst = SQUARE(atom->GetX() - atom2->GetX());
                    dst += SQUARE(atom->GetY() - atom2->GetY());
                    dst += SQUARE(atom->GetZ() - atom2->GetZ());
                    dst = sqrt(dst);
                    sprintf(buffer,"%10.4f ",dst);
                    ofs << buffer;
                }
            ofs << endl;
        }
        max += columns - 1;
        min += columns - 1;
    }
    ofs << endl;
}

void ReportFormat::WriteTorsions(ostream &ofs,OBMol &mol)
{
    vector<OBEdgeBase*>::iterator bi1,bi2,bi3;
    OBBond* bond;
    OBAtom *a,*b,*c,*d;
    char buffer[BUFF_SIZE];

    //loop through all bonds generating torsions
    for(bond = mol.BeginBond(bi1); bond; bond = mol.NextBond(bi1))
    {
        b = bond->GetBeginAtom();
        c = bond->GetEndAtom();

        for(a = b->BeginNbrAtom(bi2);a;a = b->NextNbrAtom(bi2))
        {
            if(a == c)
                continue;

            for(d = c->BeginNbrAtom(bi3);d;d = c->NextNbrAtom(bi3))
            {
                if(d == b)
                    continue;

                sprintf(buffer,"%4d %4d %4d %4d %10.3f",
                        a->GetIdx(), b->GetIdx(),c->GetIdx(),d->GetIdx(),
                        CalcTorsionAngle(a->GetVector(), b->GetVector(),
                                         c->GetVector(), d->GetVector()));
                ofs << buffer << endl;
            }
        }
    }
}

void ReportFormat::WriteAngles(ostream &ofs,OBMol &mol)
{
    // Alas, we still need to sort these to only list unique entries...
    vector3 v1, v2;
    OBAtom *a, *b, *c, *d;
    OBBond *bond1, *bond2, *bond3;
    vector<OBEdgeBase*>::iterator i, j, k;
    char buffer[BUFF_SIZE];

    for (bond1 = mol.BeginBond(i); bond1; bond1 = mol.NextBond(i))
    {
        b = bond1->GetBeginAtom();
        c = bond1->GetEndAtom();

        for (bond2 = b->BeginBond(j); bond2; bond2 = b->NextBond(j))
        {
            if (bond2->GetEndAtomIdx() != c->GetIdx()
                    && bond2->GetEndAtomIdx() != b->GetIdx())
            {
                a = bond2->GetEndAtom();

                v1 = a->GetVector() - b->GetVector();
                v2 = c->GetVector() - b->GetVector();

                sprintf(buffer,"%4d %4d %4d %4s %4s %4s %10.3f",
                        a->GetIdx(),b->GetIdx(),c->GetIdx(),
                        a->GetType(),b->GetType(),c->GetType(),
                        vectorAngle(v1, v2));
                ofs << buffer << endl;

                for (bond3 = c->BeginBond(k); bond3; bond3 = c->NextBond(k))
                    if (bond3->GetEndAtomIdx() != b->GetIdx()
                            && bond3->GetEndAtomIdx() != c->GetIdx())
                    {
                        d = bond3->GetEndAtom();

                        v1 = b->GetVector() - c->GetVector();
                        v2 = d->GetVector() - c->GetVector();

                        sprintf(buffer,"%4d %4d %4d %4s %4s %4s %10.3f",
                                b->GetIdx(),c->GetIdx(),d->GetIdx(),
                                b->GetType(),c->GetType(),d->GetType(),
                                vectorAngle(v1, v2));
                        ofs << buffer << endl;
                    }
            }
        }
    }
}

void ReportFormat::WriteChiral(ostream &ofs,OBMol &mol)
{
    OBAtom *atom;
    vector<OBNodeBase*>::iterator i;
    char buffer[BUFF_SIZE];

    for (atom = mol.BeginAtom(i);atom;atom = mol.NextAtom(i))
    {
        if (atom->IsChiral())
        {
            sprintf(buffer,"%4s %5d is chiral: %s",
                    etab.GetSymbol(atom->GetAtomicNum()),
                    atom->GetIdx(),
                    (atom->IsClockwise() ? "clockwise" : "counterclockwise"));

            ofs << buffer << endl;
        }
    }
}

} //namespace OpenBabel
