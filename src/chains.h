/**********************************************************************
chains.h - Parse for macromolecule chains and residues
 
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison
 
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

#ifndef OB_CHAINS_H
#define OB_CHAINS_H

#define MaxMonoAtom 20
#define MaxMonoBond 20

namespace OpenBabel
{

class OBAtom;
class OBMol;

//! Structure for atomic patterns (templates) in residues for OBChainsParser
typedef struct
{
    int flag;
    short elem, count;
    int n1, n2, n3, n4;
}
Template;

//! \brief Perceives peptide or nucleotide chains and residues in an OBMol
class OBAPI OBChainsParser
{
public:

    OBChainsParser(void);
    ~OBChainsParser(void);

    bool PerceiveChains(OBMol &, bool nukeSingleResidue = false);

private: // methods

    bool  DetermineHetAtoms(OBMol &);
    bool  DetermineConnectedChains(OBMol &);
    bool  DeterminePeptideBackbone(OBMol &);
    bool  DeterminePeptideSidechains(OBMol &);
    bool  DetermineNucleicBackbone(OBMol &);
    bool  DetermineNucleicSidechains(OBMol &);
    bool  DetermineHydrogens(OBMol &);

    void  SetupMol(OBMol &);
    void  SetResidueInformation(OBMol &, bool nukeSingleResidue);
    void  ClearResidueInformation(OBMol &);
    void  CleanupMol(void);

    void  AssignResidue(OBMol &, int, int, int);
    int   IdentifyResidue(void *, OBMol &, int, int); // ByteCode *

    void  DefineMonomer(void **, int, char *); // ByteCode **
    int   IdentifyElement(char *);

    bool  MatchConstraint(OBAtom *, int);
    bool  Match2Constraints(Template *, OBAtom *, OBAtom *);
    bool  Match3Constraints(Template *, OBAtom *, OBAtom *, OBAtom *);
    bool  Match4Constraints(Template *, OBAtom *, OBAtom *, OBAtom *, OBAtom *);

    void  ConstrainBackbone(OBMol &, Template *, int);

    int   RecurseChain(OBMol &, int, int);
    void  TraceNucleicChain(OBMol &, int, int);
    void  TracePeptideChain(OBMol &, int, int);

    char *ParseSmiles(char *, int);

private: // members

    void *PDecisionTree; // ByteCode *
    void *NDecisionTree; // ByteCode *

    int   ResMonoAtom[MaxMonoAtom];
    int   ResMonoBond[MaxMonoBond];

    unsigned short *bitmasks;
    unsigned char  *resids;
    unsigned char  *flags;
    bool           *hetflags;
    short          *atomids;
    short          *resnos;
    short          *sernos;
    char           *hcounts;
    char           *chains;
};

}
#endif // OB_CHAINS_H

//! \file chains.h
//! \brief Parse for macromolecule chains and residues.
