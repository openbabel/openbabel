/**********************************************************************
chains.h - Parse for macromolecule chains and residues
 
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
 
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

  // Structure template for atomic patterns in residues for OBChainsParser
  // implementation in chains.cpp
  typedef struct Template Template;

  //! \class OBChainsParser chains.h <openbabel/chains.h>
  //! \brief Perceives peptide or nucleotide chains and residues in an OBMol
  //
  //! Perceive peptide or nucleotide chains and residues from atom connectivity.
  //! Based on original RasMol code by Roger Sayle and modified by Joe Corkery.
  //! For more on Roger's original talk, see:
  //!  http://www.daylight.com/meetings/mug96/sayle/sayle.html
  class OBAPI OBChainsParser
    {
    public:

      OBChainsParser(void);
      ~OBChainsParser(void);

      //! Perceive macromolecular (peptide and nucleotide) residues and chains
      //! \param mol - the molecule to parse and update
      //! \param nukeSingleResidue - if only one residue is found, clear information
      //!   default = false  -- single residue files should still be recognized.
      bool PerceiveChains(OBMol &mol, bool nukeSingleResidue = false);

    private: // internal methods

      //! Determine HETATOM records for ligands vs. possible residues
      bool  DetermineHetAtoms(OBMol &);
      //! Determine connected chains (e.g., subunits)
      bool  DetermineConnectedChains(OBMol &);
      //! Walk a peptide "backbone" atom sequence, from one residue to the next
      bool  DeterminePeptideBackbone(OBMol &);
      //! Identify peptide residues based on side chains
      bool  DeterminePeptideSidechains(OBMol &);
      //! Walk a nucleic acid "backone" from one residue to the next
      bool  DetermineNucleicBackbone(OBMol &);
      //! Identify nucleic acids based on "side chains"
      bool  DetermineNucleicSidechains(OBMol &);
      //! Determine any bonded hydrogens based on peptide/nucleic residue info
      bool  DetermineHydrogens(OBMol &);

      //! Set up the chain perception to operate on the supplied molecule
      void  SetupMol(OBMol &);
      void  SetResidueInformation(OBMol &, bool nukeSingleResidue);
      void  ClearResidueInformation(OBMol &);
      void  CleanupMol(void);

      void  AssignResidue(OBMol &, int, int, int);
      int   IdentifyResidue(void *, OBMol &, int, int); // ByteCode *

      void  DefineMonomer(void **, int, const char *); // ByteCode **
      int   IdentifyElement(const char *);

      bool  MatchConstraint(OBAtom *, int);
      bool  Match2Constraints(Template *, OBAtom *, OBAtom *);
      bool  Match3Constraints(Template *, OBAtom *, OBAtom *, OBAtom *);
      bool  Match4Constraints(Template *, OBAtom *, OBAtom *, OBAtom *, OBAtom *);

      void  ConstrainBackbone(OBMol &, Template *, int);

      int   RecurseChain(OBMol &, int, int);
      void  TraceNucleicChain(OBMol &, int, int);
      void  TracePeptideChain(OBMol &, int, int);

      const char *ParseSmiles(const char *, int);

    private: // members

      void *PDecisionTree; // ByteCode *
      void *NDecisionTree; // ByteCode *

      int   ResMonoAtom[MaxMonoAtom];
      int   ResMonoBond[MaxMonoBond];

      unsigned short *bitmasks; 
      bool           *visits;   //!< mark visits to prevent looping
      unsigned char  *resids;
      unsigned char  *flags;
      bool           *hetflags;
      int            *atomids;
      short          *resnos;
      short          *sernos;   //!< array of residue serial numbers
      char           *hcounts;
      char           *chains;
    };

}
#endif // OB_CHAINS_H

//! \file chains.h
//! \brief Parse for macromolecule chains and residues.
