/**********************************************************************
chains.h - Parse for macromolecule chains and residues

Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
Some portions Copyright (C) 2008 by Tim Vandermeersch

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

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

#include <openbabel/babelconfig.h>
#include <vector>

namespace OpenBabel
{

  class OBAtom;
  class OBMol;

  //! Structure template for atomic patterns in residues for OBChainsParser
  // implementation in chains.cpp
  struct Template;
  typedef struct Template Template;

  /** @class OBChainsParser chains.h <openbabel/chains.h>
      @brief Perceives peptide or nucleotide chains and residues in an OBMol

      Perceive peptide or nucleotide chains and residues from atom connectivity.
      Based on original RasMol code by Roger Sayle and modified by Joe Corkery.
      For more on Roger's original talk, see:
      http://www.daylight.com/meetings/mug96/sayle/sayle.html
   */
  class OBAPI OBChainsParser
  {
    public:

      OBChainsParser(void);
      ~OBChainsParser(void);

      /**
       * Perceive macromolecular (peptide and nucleotide) residues and chains
       * @param mol The molecule to parse and update
       * @param nukeSingleResidue If only one residue is found, clear information
       * default = false  -- single residue files should still be recognized.
       */
      bool PerceiveChains(OBMol &mol, bool nukeSingleResidue = false);

    private: // internal methods

      //! @name Step 1: Determine hetero atoms
      //@{
      /**
       * Determine HETATOM records for all atoms with a heavy valance of 0.
       * This includes HOH, Cl, Fe, ...
       *
       * Sets resids[i] & hetflags[i] for these atoms.
       * @todo add ions (Cl, Fe, ...)
       */
      bool DetermineHetAtoms(OBMol &);
      //@}

      //! @name Step 2: Determine connected chains
      //@{
      /**
       * Determine connected chains (e.g., subunits). Chains will be labeled A, B, C, ...
       * Ligands also get assigned a chain label. The chain for ligands will later be
       * replaced by ' '. The residue numbers will also be updated in this process to
       * make sure all ligands, HOH, ions, etc. have a unique residue number in the ' '
       * chain.
       *
       * Sets chains[i] for all atoms. (through RecurseChain())
       */
      bool DetermineConnectedChains(OBMol &);
      /**
       * Perform the actual work for DetermineConnectedChains(). Set chains[i]
       * to @p c for all atoms of the recursed chain.
       * @param mol The molecule.
       * @param i Index for the current atom. (RecurseChain() will be called for all neighbours)
       * @param c The chain which we are recusring. ('A' + count)
       * @return The number of heavy atoms in the recursed chain.
       */
      unsigned int RecurseChain(OBMol &mol, unsigned int i, int c);
      //@}

      //! @name Step 3: Determine peptide backbone
      //@{
      /**
       * Walk a peptide "backbone" atom sequence, from one residue to the next. This
       * function will look for N-CA-C-O sequences and mark these atoms.
       *
       * Sets bitmaks[i] for these atoms. (through ConstrainBackbone())
       * Sets resnos[i] for these atoms. (through TracePeptideChain())
       */
      bool DeterminePeptideBackbone(OBMol &);
      /**
       * First the bitmasks[i] will be OR-ed with Template::flag for all atoms based on
       * on Template::element and Template::count.
       *
       * Next, the bitmasks[i] are iteratively resolved by matching the
       * constraints in OpenBabel::Peptide or OpenBabel::Nucleotide.
       * @param mol The molecule.
       * @param templ OpenBabel::Peptide or OpenBabel::Nucleotide
       * @param tmax Number of entries in @p templ
       */
      void ConstrainBackbone(OBMol &mol, Template *templ, int tmax);
      /**
       * @return True if the bitmasks[i] for @p atom matches @p mask.
       */
      bool MatchConstraint(OBAtom *atom, int mask);
      /**
       * @return True if atom @p na and @p nb match the Template::n1 and
       * Template::n2.
       */
      bool Match2Constraints(Template *templ, OBAtom *na, OBAtom *nb);
      /**
       * @return True if atom @p na, @p nb and @p nc match the Template::n1,
       * Template::n2 and Template::n3.
       */
      bool Match3Constraints(Template *templ, OBAtom *na, OBAtom *nb, OBAtom *nc);
      /**
       * @return True if atom @p na, @p nb, @p nc and @p nd match the Template::n1,
       * Template::n2, Template::n3 and Template::n4.
       */
      bool Match4Constraints(Template *templ, OBAtom *na, OBAtom *nb, OBAtom *nc, OBAtom *nd);
      /**
       * Now we have the constrained bitmaks[i], trace N-CA-C-O-... and set
       * resnos[i] and atomids[i] for each N-CA-C-O sequence.
       *
       * Also adds BF_DOUBLE to flags[b] for< each carbonyl bond in N-CA-C=O.
       * @param mol The molecule.
       * @param i Index for the current atom. (TracePeptideChain() will be called for all neighbours)
       * @param r The residue number which we are tracing.
       */
      void TracePeptideChain(OBMol &mol, unsigned int i, int r);
      //@}

      //! @name Step 4: Determine peptide side chains
      //@{
      /**
       * Look for atoms with atomids[i] CA and identify their side chain.
       *
       * Sets resnos[i] and resids[i] for all identified residues (including the N-CA-C-O).
       * (through IdentifyResidue() and AssignResidue())
       */
      bool  DeterminePeptideSidechains(OBMol &);
      /**
       * Identify a residue based on the @p tree ByteCode.
       *
       * Sets resnos[i] for all sidechain atoms to the residue number of
       * the seed CA/C1 atom.
       * @param tree Bytecode for the residues. (OBChainsParser::PDecisionTree or OBChainsParser::NDecisionTree)
       * @param mol The molecule.
       * @param seed Atom index for the CA (peptides) or C1 (nucleotides) atom.
       * @param resno The residue number for this residue.
       * @return The resids[i] for the identified residue.
       */
      int IdentifyResidue(void *tree, OBMol &mol, unsigned int seed, int resno); // ByteCode *
      /**
       * Set resids[i] for all atoms where resids[i] = @p r and chains[i] = @p c.
       * @param mol The molecule.
       * @param r The residue number.
       * @param c The chain number.
       * @param i The residue id (resids[i] returned by IdentifyResidue())
       */
      void  AssignResidue(OBMol &mol, int r, int c, int i);
      //@}

      //! @name Step 5: Assign hydrogens
      //@{
      /**
       * Assign the resids[i], resnos[i], ... for all hydrogens based on the
       * atom they are bound to.
       */
      bool  DetermineHydrogens(OBMol &);
      //@}

      //! @name Step 6: Set the residue information
      //@{
      /**
       * Convert the private data vectors to OBResidue objects and add them to @p mol.
       * @param mol The molecule to parse and update
       * @param nukeSingleResidue If only one residue is found, clear information
       * default = false  -- single residue files should still be recognized.
       */
      void  SetResidueInformation(OBMol &, bool nukeSingleResidue);
      //@}

      //! @name Nucleic acids (analog to peptides)
      //@{
      /**
       * Walk a nucleic "backbone" atom sequence, from one residue to the next. This
       * function will look for ribose-5-P sequences and mark these atoms.
       *
       * Sets bitmaks[i] for these atoms. (through ConstrainBackbone())
       * Sets resnos[i] for these atoms. (through TraceNucleicChain())
       */
      bool  DetermineNucleicBackbone(OBMol &);
      /**
       * Now we have the constrained bitmaks[i], trace nucleic backbone and set
       * resnos[i] and atomids[i] for each ribose-5-P sequence.
       * @param mol The molecule.
       * @param i Index for the current atom. (TraceNucleicChain() will be called for all neighbours)
       * @param r The residue number which we are tracing.
       */
      void  TraceNucleicChain(OBMol &, unsigned int i, int r);
      /**
       * Look for atoms with atomids[i] C1 and identify their side chain.
       *
       * Sets resnos[i] and resids[i] for all identified residues.
       * (through IdentifyResidue() and AssignResidue())
       */
      bool  DetermineNucleicSidechains(OBMol &);
      //@}

      /**
       * Set up the chain perception to operate on the supplied molecule
       * by resizing and initializing the private data vectors.
       */
      void  SetupMol(OBMol &);
      /**
       * Delete all residues in @p mol
       */
      void  ClearResidueInformation(OBMol &mol);
      /**
       * Clear all private data vectors
       */
      void CleanupMol();
      /**
       * Construct and add ByteCode to the @p tree for a single residue.
       * @param tree Bytecode for the residues. (OBChainsParser::PDecisionTree or OBChainsParser::NDecisionTree)
       * @param resid The residue id.
       * @param smiles The pseudo-smiles string (from OpenBabel::AminoAcids or OpenBabel::Nucleotides)
       */
      void  DefineMonomer(void **tree, int resid, const char *smiles); // ByteCode **
      /**
       * @param ptr Element id (from OpenBabel::ChainsAtomName)
       * @return The element number.
       */
      int   IdentifyElement(char *ptr);
      /**
       * Parse a pseudo smiles from OpenBabel::AminoAcids or OpenBabel::Nucleotides.
       * @param smiles The pseudo-smiles string.
       * @param prev The previous position (used for recursing, use -1 to start).
       */
      const char *ParseSmiles(const char *smiles, int prev);
      /**
       * Debugging function.
       */
      void DumpState();

      void *PDecisionTree; //!< ByteCode decision tree for peptides
      void *NDecisionTree; //!< ByteCode decision tree for nucleotides

      int   ResMonoAtom[MaxMonoAtom];
      int   ResMonoBond[MaxMonoBond];

      std::vector<unsigned short> bitmasks;
      std::vector<bool>           visits;   //!< mark visits to prevent looping
      std::vector<unsigned char>  resids;
      std::vector<unsigned char>  flags;
      std::vector<bool>           hetflags;
      std::vector<int>            atomids;
      std::vector<short>          resnos;
      std::vector<short>          sernos;   //!< array of residue serial numbers
      std::vector<char>           hcounts;
      std::vector<char>           chains;
    };

    //! Global OBChainsParser for detecting macromolecular chains and residues
    OB_EXTERN  OBChainsParser   chainsparser;

}
#endif // OB_CHAINS_H

//! \file chains.h
//! \brief Parse for macromolecule chains and residues.
