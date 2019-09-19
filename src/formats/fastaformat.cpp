/**********************************************************************
Copyright (C) 2006 by Sangwoo Shim

Read sequence functions, helix generation and completed IUPAC coverage.
Copyright (C) Scarlet Line 2007-9

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

#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/elements.h>
#include <openbabel/obiter.h>
#include <openbabel/bond.h>

#include <map>
#include <string>
#include <cstdlib>

using namespace std;
namespace OpenBabel
{

  class FASTAFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID in the constructor
    FASTAFormat() {
      OBConversion::RegisterFormat("fasta", this, "chemical/x-fasta");
      OBConversion::RegisterFormat("fa", this);
      OBConversion::RegisterFormat("fsa", this);

      OBConversion::RegisterOptionParam("s", this);
      OBConversion::RegisterOptionParam("b", this);
      OBConversion::RegisterOptionParam("n", this);
      OBConversion::RegisterOptionParam("1", this);
      OBConversion::RegisterOptionParam("t", NULL, 1, OBConversion::INOPTIONS);
    }

    virtual const char* Description() //required
    {
      return
        "FASTA format\n"
        "A file format used to exchange information between genetic sequence databases\n\n"
        "Read Options e.g. -as\n"
        "  1  Output single-stranded DNA\n"
        "  t <turns>  Use the specified number of base pairs per turn (e.g., 10)\n"
        "  s  Output single bonds only\n"
        "  b  Disable bonding entirely\n\n"
        "Write Options e.g. -xn \n"
        "  n  Omit title and comments\n";
    };

    virtual const char* SpecificationURL() {
      return "http://www.ebi.ac.uk/help/formats_frame.html";
    };
    // Additionally http://www.ncbi.nlm.nih.gov/blast/fasta.shtml

    virtual const char* GetMIMEType()
    { return "chemical/x-fasta"; }

    virtual unsigned int Flags() {
      return READONEONLY | WRITEONEONLY;
    };

    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

    enum SequenceType
      {
        UnknownSequence,
        ProteinSequence,
        DNASequence,
        RNASequence,
        MAXSequence
      };
  private:
    char conv_3to1(const string & three) const;
  };

  FASTAFormat theFASTAFormat;

  bool FASTAFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    string seq;
    OBMol* pmol;
    //   OBResidue *res;

    pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol == NULL)
      return false;
    ostream &ofs = *pConv->GetOutStream();

    int seq_count = 0;
    FOR_RESIDUES_OF_MOL(res,pmol) {
      if (res->GetAtoms().size() < 3)
        continue;
      seq.append(1, conv_3to1(res->GetName()));
      ++ seq_count;
      if (seq_count >= 60) {
        seq_count = 0;
        seq.append("\n");
      }
    }
    if(!pConv->IsOption("n")) {
      if (strlen(pmol->GetTitle()) > 0)
        ofs << ">" << pmol->GetTitle();
      else
        ofs << ">Unknown molecule";
      ofs << " " << pmol->NumResidues() << " bp";
      ofs << "; generated with OpenBabel " << BABEL_VERSION << endl;
    }
    ofs << seq << endl;
    return true;
  }
  struct residue_lookup_record
  {
    char TLA[4];
    char code;
  };
  residue_lookup_record residue_lookup_table[] =
    {
      { "ALA", 'A' }, // Alanine
      { "ARG", 'R' }, // Arginine
      { "ASN", 'N' }, // Asparagine
      { "ASP", 'D' }, // Aspartic acid
      { "CYS", 'C' }, // Cysteine
      { "GLN", 'Q' }, // Glutamine
      { "GLU", 'E' }, // Glutamic acid
      { "GLY", 'G' }, // Glycine
      { "HIS", 'H' }, // Histidine
      { "ILE", 'I' }, // Isoleucine
      { "LEU", 'L' }, // Leucine
      { "LYS", 'K' }, // Lysine
      { "MET", 'M' }, // Methionine
      { "PHE", 'F' }, // Phenylalanine
      { "PRO", 'P' }, // Proline
      { "SER", 'S' }, // Serine
      { "THR", 'T' }, // Threonine
      { "TRP", 'W' }, // Tryptophan
      { "TYR", 'Y' }, // Tyrosine
      { "VAL", 'V' }, // Valine
      { "SEC", 'U' }, // Selenocysteine
      { "PYL", 'O' }, // Pyrrolysine
      { "ASX", 'B' }, // Asparagine or Aspartic acid
      { "GLX", 'Z' }, // Glutamine or Glutamic acid
      { "XLE", 'J' }, // Leucine or Isoleucine
      { "XAA", 'X' }, // any amino acid

      { "DG", 'G' }, // Guanine
      { "DC", 'C' }, // Cytosine
      { "DA", 'A' }, // Adenine
      { "DT", 'T' }, // Thymine
      { "GAX", 'R' }, // Purine (G or A)
      { "TCX", 'Y' }, // Pyrimidine (T or C)
      { "UCX", 'Y' }, // Pyrimidine (U or C)
      { "GTX", 'K' }, // Keto (G or T)
      { "GUX", 'K' }, // Keto (G or U)
      { "ACX", 'M' }, // Amino (A or C)
      { "GCX", 'S' }, // Strong (G or C)
      { "ATX", 'W' }, // Weak (A or T)
      { "AUX", 'W' }, // Weak (A or U)
      { "GTC", 'B' }, // (G T C)
      { "GUC", 'B' }, // (G U C)
      { "GAT", 'D' }, // (G A T)
      { "GAU", 'D' }, // (G A U)
      { "ACT", 'H' }, // (A C T)
      { "ACU", 'H' }, // (A C U)
      { "GCA", 'V' }, // (G C A)
      { "XNA", 'N' }, // Unknown Nucleic Acid
      { "", '\0' }
    };
  typedef std::map< std::string, char > residue_lookup_map;
  residue_lookup_map residue_lookup;

  char
  FASTAFormat::conv_3to1(const std::string & three) const
  {
    char return_code = 'X';
    if (residue_lookup.empty())
      {
        for (residue_lookup_record * rl = residue_lookup_table; rl->code; ++ rl)
          {
            residue_lookup.insert(residue_lookup_map::value_type(std::string(rl->TLA), rl->code));
          }
      }
    char tla_buf[4];
    strncpy(tla_buf, three.c_str(), 3);
    tla_buf[3] = 0;
    for (int idx = 0; idx < 3; ++ idx)
      tla_buf[idx] = (char)toupper(tla_buf[idx]);

    residue_lookup_map::const_iterator mx = residue_lookup.find(std::string(tla_buf));
    if (mx != residue_lookup.end())
      {
        return_code = (* mx).second;
      }
    else if (strlen(tla_buf) == 1)
      {
        return_code = tla_buf[0];
      }
    return return_code;
  }
#define IUPAC_Start 0
#define IUPAC_End 1
#define IUPAC_Unknown 2
  const char * IUPAC_DNA_codes = "01NACGTRYKMSWBDHV";
  enum IUPAC_DNA_code
    {
      IUPAC_DNA_Start = 0,
      IUPAC_DNA_End = 1,
      IUPAC_DNA_N, // N XNA Unknown Nucleic Acid N -> N
      IUPAC_DNA_A, // A DA  Adenine A -> T
      IUPAC_DNA_C, // C DC  Cytosine C -> G
      IUPAC_DNA_G, // G DG  Guanine G -> C
      IUPAC_DNA_T, // T DT  Thymine T -> A
      IUPAC_DNA_R, // R GAX Purine (G or A) R -> Y
      IUPAC_DNA_Y, // Y TCX Pyrimidine (T or C) Y -> R
      IUPAC_DNA_K, // K GTX Keto (G or T) K -> M
      IUPAC_DNA_M, // M ACX Amino (A or C) M -> K
      IUPAC_DNA_S, // S GCX Strong (G or C) S -> W
      IUPAC_DNA_W, // W ATX Weak (A or T) W -> S
      IUPAC_DNA_B, // B GTC (G T C) B -> V
      IUPAC_DNA_D, // D GAT (G A T) D -> H
      IUPAC_DNA_H, // H ACT (A C T) H -> D
      IUPAC_DNA_V, // V GCA (G C A) V -> B
      IUPAC_DNA_max
      // -  gap of indeterminate length
    };
  const char * IUPAC_RNA_codes = "01NACGURYKMSWBDHV";
  enum IUPAC_RNA_code
    {
      IUPAC_RNA_Start = 0,
      IUPAC_RNA_End = 1,
      IUPAC_RNA_N, // N XNA Unknown Nucleic Acid N -> N
      IUPAC_RNA_A, // A A   Adenine A -> U
      IUPAC_RNA_C, // C C   Cytosine C -> G
      IUPAC_RNA_G, // G G   Guanine G -> C
      IUPAC_RNA_U, // U U   Uridine U -> A
      IUPAC_RNA_R, // R GAX Purine (G or A) R -> Y
      IUPAC_RNA_Y, // Y UCX Pyrimidine (U or C) Y -> R
      IUPAC_RNA_K, // K GUX Keto (G or U) K -> M
      IUPAC_RNA_M, // M ACX Amino (A or C) M -> K
      IUPAC_RNA_S, // S GCX Strong (G or C) S -> W
      IUPAC_RNA_W, // W AUX Weak (A or U) W -> S
      IUPAC_RNA_B, // B GUC (G U C) B -> V
      IUPAC_RNA_D, // D GAU (G A U) D -> H
      IUPAC_RNA_H, // H ACU (A C U) H -> D
      IUPAC_RNA_V, // V GCA (G C A) V -> B
      IUPAC_RNA_max
      // -  gap of indeterminate length
    };
  /*
  // Additional RNA Not yet implemented:
  T --> ribosylthymine (not thymidine)
  I --> inosine
  X --> xanthosine
  Q --> pseudouridine
  */
  const char * IUPAC_Protein_codes = "01XABCDEFGHIJKLMNOPQRSTUVWYZ";
  enum IUPAC_Protein_code
    {
      IUPAC_Protein_Start = 0,
      IUPAC_Protein_End = 1,
      IUPAC_Protein_X, // X XAA any (unknown) amino acid
      IUPAC_Protein_A, // A ALA Alanine
      IUPAC_Protein_B, // B ASX Asparagine or Aspartic acid i.e. N or D
      IUPAC_Protein_C, // C CYS Cysteine
      IUPAC_Protein_D, // D ASP Aspartic acid
      IUPAC_Protein_E, // E GLU Glutamic acid
      IUPAC_Protein_F, // F PHE Phenylalanine
      IUPAC_Protein_G, // G GLY Glycine
      IUPAC_Protein_H, // H HIS Histidine
      IUPAC_Protein_I, // I ILE Isoleucine
      IUPAC_Protein_J, // J XLE Leucine or Isoleucine i.e. L or I
      IUPAC_Protein_K, // K LYS Lysine
      IUPAC_Protein_L, // L LEU Leucine
      IUPAC_Protein_M, // M MET Methionine
      IUPAC_Protein_N, // N ASN Asparagine
      IUPAC_Protein_O, // O PYL Pyrrolysine
      IUPAC_Protein_P, // P PRO Proline
      IUPAC_Protein_Q, // Q GLN Glutamine
      IUPAC_Protein_R, // R ARG Arginine
      IUPAC_Protein_S, // S SER Serine
      IUPAC_Protein_T, // T THR Threonine
      IUPAC_Protein_U, // U SEC Selenocysteine
      IUPAC_Protein_V, // V VAL Valine
      IUPAC_Protein_W, // W TRP Tryptophan
      IUPAC_Protein_Y, // Y TYR Tyrosine
      IUPAC_Protein_Z, // Z GLX Glutamine or Glutamic acid i.e. Q or E
      IUPAC_Protein_max
      // *  translation stop
      // -  gap of indeterminate length
    };

  struct ResidueAtomRecord
  {
    char label[6]; // e.g. CA, CB, CC
    char symbol[4]; // e.g. C, H, N, O, P
    double x, r, Theta; // cylindrical polar co-ordinates
  };
  struct ResidueBondRecord
  {
    size_t from_idx, to_idx;
    int bond_order;
  };
  struct ResidueRecord
  {
    char IUPACcode;
    char Name[6]; // usually three-letter code
    ResidueAtomRecord  atom[48];
    ResidueBondRecord  bond[48];
  };
  extern ResidueRecord RNAResidues[IUPAC_RNA_max];
  extern ResidueRecord DNAResidues[IUPAC_DNA_max];
  extern ResidueRecord DNAPairResidues[IUPAC_DNA_max];
  extern ResidueRecord ProteinResidues[IUPAC_Protein_max];
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795028841971694
#endif
  struct HelixParameters
  {
    double unit_X; // distance along the x-axis added for each unit
    double unit_Theta; // angle in radians between each unit
    size_t  bond_connect; // the zero-based index of the residue atom to bond from
  };
  HelixParameters protein_helix = { 1.4781, 1.74474, 2 }; // HNCA-.N-CA-C-.N-CA-C-.OXT-HOCA
  HelixParameters RNA_helix = { 2.95057, 2*M_PI/11.0, 8 }; // HTER-OXT-.P-O5'-C5'-C4'-C3'-O3'-.P-O5'-C5'-C4'-C3'-O3'-.HCAP
  HelixParameters DNA_helix = { 3.37998, 2*M_PI/10.0, 8 }; // HTER-OXT-.P-O5'-C5'-C4'-C3'-O3'-.P-O5'-C5'-C4'-C3'-O3'-.HCAP
  HelixParameters DNA_pair_helix = { -3.37998, -2*M_PI/10.0, 8 }; // Simply the negative of the above
  // HelixParameters DNA_strandAtoB = { 4.1600027464119, 3.3230513172128164 };
  typedef OBAtom * ptrAtom;
  void add_bond(OBMol * pmol, OBAtom * from, OBAtom * to, int bond_order)
  {
    OBBond * bond = pmol->NewBond();
    bond->SetBegin(from);
    bond->SetEnd(to);
    bond->SetBondOrder(bond_order);
  }
  // #define GEN_BONDS 1
  void add_residue(OBMol * pmol, OBResidue * res, double offset_x, double offset_Theta, unsigned long & serial_no, ResidueRecord * res_rec, int resBondFromOffset, ptrAtom & resBondFrom, bool create_bonds, bool bond_orders)
  {
    typedef std::vector<OBAtom *> atom_list;
    atom_list bond_refs;
    for (ResidueAtomRecord * atom_rec = res_rec->atom; atom_rec->symbol[0]; ++ atom_rec)
      {
        OBAtom * atom = pmol->NewAtom();
        atom->SetAtomicNum(OBElements::GetAtomicNum(atom_rec->symbol));
        atom->SetType(atom_rec->symbol);
        double theta = offset_Theta + atom_rec->Theta;
        atom->SetVector(offset_x + atom_rec->x, atom_rec->r * cos(theta), atom_rec->r * sin(theta));
        res->AddAtom(atom);
        res->SetAtomID(atom, atom_rec->label);
        res->SetSerialNum(atom, serial_no);
        ++ serial_no;
        bond_refs.push_back(atom);
      }

    if (create_bonds)
      {
        size_t atom_count = bond_refs.size();
        if (resBondFrom && atom_count)
          { // insert the bond from the previous residue
            add_bond(pmol, resBondFrom, bond_refs[0], 1);
          }
        resBondFrom = 0;
        for (ResidueBondRecord * bond_rec = res_rec->bond; bond_rec->bond_order; ++ bond_rec)
          {
            size_t from = bond_rec->from_idx - 1;
            size_t to = bond_rec->to_idx - 1;
            if (from < atom_count && to < atom_count)
              {
                add_bond(pmol, bond_refs[from], bond_refs[to], bond_orders?bond_rec->bond_order:1);
              }
          }
        if (atom_count && resBondFromOffset != -2)
          {
            if (resBondFromOffset == -1)
              { // start
                resBondFrom = bond_refs[atom_count - 1];
              }
            else if (atom_count > resBondFromOffset)
              {
                resBondFrom = bond_refs[resBondFromOffset];
              }
          }
      }

    bond_refs.clear();
  }


  void generate_sequence(const std::string & sequence, OBMol * pmol, unsigned long chain_no, const HelixParameters & helix, const char * IUPAC_codes, ResidueRecord * Residues, double & offset_x, double & offset_Theta, unsigned long & serial_no, bool create_bonds, bool bond_orders)
  {
    unsigned long residue_num = 1;
    OBResidue * res = 0;
    ptrAtom resBondFrom = 0;
    for (std::string::const_iterator sx = sequence.begin(), sy = sequence.end(); sx != sy; ++ sx, ++ residue_num)
      {
        bool is_gap = (((* sx) == '-') || ((* sx) == '*'));
        if (is_gap)
          {
            offset_x += helix.unit_X * 2;
            resBondFrom = 0;
          }
        else
          {
            const char * idx = strchr(IUPAC_codes, (* sx)); // e.g. "01NACGURYKMSWBDHV"
            size_t unit_code = (size_t)( idx ? (idx - IUPAC_codes) : IUPAC_Unknown );
            ResidueRecord * res_rec = & Residues[unit_code];
            if (res_rec->IUPACcode)
              {
                res = pmol->NewResidue();
                res->SetChainNum(chain_no);
                res->SetNum(residue_num);
                res->SetName(res_rec->Name);

                if (residue_num == 1)
                  { // Add the start terminal
                    add_residue(pmol, res, offset_x, offset_Theta, serial_no, & Residues[IUPAC_Start], -1, resBondFrom, create_bonds, bond_orders);
                  }

                add_residue(pmol, res, offset_x, offset_Theta, serial_no, res_rec, helix.bond_connect, resBondFrom, create_bonds, bond_orders);
              }
            offset_x += helix.unit_X;
            offset_Theta += helix.unit_Theta;
          }
      }

    if (res != 0)
      { // Add the end terminal
        add_residue(pmol, res, offset_x - helix.unit_X, offset_Theta - helix.unit_Theta, serial_no, & Residues[IUPAC_End], -2, resBondFrom, create_bonds, bond_orders);
      }

  }

  bool ReadFASTASequence(OBMol * pmol, int seq_type, std::istream * in, bool create_bonds, bool bond_orders,
                         bool singleStrand, const char *turns = 0)
  {
    /*
      Sequence is from 5' to 3' left -> right.

    */
    std::string line, sequence;

    FASTAFormat::SequenceType sequence_type = (FASTAFormat::SequenceType)seq_type, sequence_na = FASTAFormat::UnknownSequence;
    while (!in->eof())
      { // We have a new line to fetch into 'line'
        getline( * in, line);
        if (line[0] == '>')
          { // comment data
            if (pmol->GetTitle()[0] == 0)
              {
                pmol->SetTitle( & (line.c_str()[1]) );
              }
            if (sequence_type == FASTAFormat::UnknownSequence)
              { // attempt to determine the sequence type
                if (line.find("RNA") != std::string::npos)
                  sequence_type = FASTAFormat::RNASequence;
                else if (line.find("DNA") != std::string::npos)
                  sequence_type = FASTAFormat::DNASequence;
                else if (line.find("gene") != std::string::npos)
                  sequence_type = FASTAFormat::DNASequence;
                else if (line.find("protein") != std::string::npos)
                  sequence_type = FASTAFormat::ProteinSequence;
                else if (line.find("peptide") != std::string::npos)
                  sequence_type = FASTAFormat::ProteinSequence;
                else if (line.find("Protein") != std::string::npos)
                  sequence_type = FASTAFormat::ProteinSequence;
                else if (line.find("Peptide") != std::string::npos)
                  sequence_type = FASTAFormat::ProteinSequence;
              }
          }
        else
          { // sequence data
            for (std::string::size_type pos = 0, endpos = line.size(); pos < endpos; ++ pos)
              {
                char current = (char)toupper(line[pos]);
                if (isupper((unsigned char)current) || strchr("*-", current))
                  {
                    sequence.append(1, current);
                    if (sequence_type == FASTAFormat::UnknownSequence)
                      { // attempt to determine the sequence type
                        // DNA/RNA == ABCDGHKMNRSTUVWY
                        if (strchr("EFIJLOPQXZ*", current)) // Protein == ABCDEFGHIJKLMNOPQRSTUVWXYZ*
                          sequence_type = FASTAFormat::ProteinSequence;
                        else if (current == 'U')
                          sequence_na = FASTAFormat::RNASequence;
                        else if (current == 'T')
                          sequence_na = FASTAFormat::DNASequence;
                      }
                  }
              }
          }
      }
    if (sequence_type == FASTAFormat::UnknownSequence)
      sequence_type = sequence_na;
    if (sequence_type == FASTAFormat::UnknownSequence)
      sequence_type = FASTAFormat::DNASequence;

    // We now have the sequence and know it's type
    double offset_x = 0, offset_Theta = 0;
    unsigned long serial_no = 1;

    if (turns) {
      // set everything even if we don't use it
      double unit_Theta = 2*M_PI/atof(turns);
      protein_helix.unit_Theta = unit_Theta;
      RNA_helix.unit_Theta = unit_Theta;
      DNA_helix.unit_Theta = unit_Theta;
      DNA_pair_helix.unit_Theta = -unit_Theta;
    }

    switch (sequence_type)
      {
      case FASTAFormat::RNASequence:
        generate_sequence(sequence, pmol, 1, RNA_helix, IUPAC_RNA_codes, RNAResidues, offset_x, offset_Theta, serial_no, create_bonds, bond_orders);
        break;
      case FASTAFormat::ProteinSequence:
        generate_sequence(sequence, pmol, 1, protein_helix, IUPAC_Protein_codes, ProteinResidues, offset_x, offset_Theta, serial_no, create_bonds, bond_orders);
        break;
      case FASTAFormat::DNASequence:
        {
          generate_sequence(sequence, pmol, 1, DNA_helix, IUPAC_DNA_codes, DNAResidues, offset_x, offset_Theta, serial_no, create_bonds, bond_orders);
          if (!singleStrand) {
            offset_x -= DNA_helix.unit_X;
            offset_Theta -= DNA_helix.unit_Theta;
            std::string rsequence;
            for (std::string::const_reverse_iterator sx = sequence.rbegin(), sy = sequence.rend(); sx != sy; ++ sx)
              rsequence.append(1, * sx);
            generate_sequence(rsequence, pmol, 2, DNA_pair_helix, IUPAC_DNA_codes, DNAPairResidues, offset_x, offset_Theta, serial_no, create_bonds, bond_orders);
          }
        }
        break;
      default:
        break;
      }
    pmol->SetChainsPerceived();
    return (pmol->NumAtoms() > 0 ? true : false);
  }

  /////////////////////////////////////////////////////////////////
  bool FASTAFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if (pmol == 0)
      return false;
    pmol->BeginModify();
    bool rv = ReadFASTASequence(pmol, UnknownSequence, pConv->GetInStream(),
                                !pConv->IsOption("b",OBConversion::INOPTIONS),
                                !pConv->IsOption("s",OBConversion::INOPTIONS),
                                pConv->IsOption("1",OBConversion::INOPTIONS),
                                pConv->IsOption("t",OBConversion::INOPTIONS));
    pmol->EndModify();
    return rv;
  }

  ResidueRecord DNAResidues[IUPAC_DNA_max] =
    {
      { 0, "", // DNA Start
        {
          { "HTER",  "H",  -1.547,  9.456,  -0.129},
          { "OXT",  "O",  -1.440,  8.905,  -0.045},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 2, 1, 1},
          { 0, 0, 0}
        }
      },
      { 0, "", // DNA End
        {
          { "HCAP",  "H",  1.699,  8.181,  0.680},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 0, 0, 0}
        }
      },
      { 'N', "XNA", // N XNA Unknown Nucleic Acid N -> N
        {
          { "P",  "P",  0.000,  8.910,  0.000},
          { "O1P",  "O",  0.220,  10.199,  0.072},
          { "O2P",  "O",  0.790,  8.819,  -0.141},
          { "O5'",  "O",  0.250,  7.730,  0.126},
          { "C5'",  "C",  -0.690,  7.700,  0.269},
          { "C4'",  "C",  0.040,  7.590,  0.442},
          { "O4'",  "O",  0.250,  6.220,  0.510},
          { "C3'",  "C",  1.440,  8.200,  0.442},
          { "O3'",  "O",  1.830,  8.750,  0.590},
          { "C2'",  "C",  2.320,  7.040,  0.384},
          { "C1'",  "C",  1.610,  5.860,  0.485},
          { "N9",  "N",  1.660,  4.630,  0.325},
          { "H5'1",  "H",  -1.179,  8.685,  0.271},
          { "H5'2",  "H",  -1.449,  6.917,  0.250},
          { "H4'",  "H",  -0.582,  8.159,  0.531},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 3, 1, 2 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 4, 5, 1 },
          { 12, 11, 1 },
          { 14, 5, 1 },
          { 5, 13, 1 },
          { 5, 6, 1 },
          { 10, 11, 1 },
          { 10, 8, 1 },
          { 11, 7, 1 },
          { 7, 6, 1 },
          { 6, 8, 1 },
          { 6, 15, 1 },
          { 8, 9, 1 },
          { 0, 0, 0}
        }
      },
      { 'A', "DA", // A DA  Adenine A -> T
        {
          { "P",  "P",  0.000,  8.910,  0.000},
          { "O1P",  "O",  0.220,  10.201,  0.072},
          { "O2P",  "O",  0.790,  8.821,  -0.141},
          { "O5'",  "O",  0.250,  7.730,  0.126},
          { "C5'",  "C",  -0.690,  7.701,  0.269},
          { "C4'",  "C",  0.040,  7.590,  0.442},
          { "O4'",  "O",  0.250,  6.221,  0.510},
          { "C3'",  "C",  1.440,  8.201,  0.442},
          { "O3'",  "O",  1.830,  8.750,  0.590},
          { "C2'",  "C",  2.320,  7.040,  0.384},
          { "C1'",  "C",  1.610,  5.860,  0.485},
          { "N9",  "N",  1.660,  4.630,  0.325},
          { "C8",  "C",  1.580,  4.840,  0.038},
          { "N7",  "N",  1.650,  3.950,  -0.178},
          { "C5",  "C",  1.800,  2.741,  0.021},
          { "C6",  "C",  1.930,  1.410,  -0.210},
          { "N6",  "N",  1.940,  1.810,  -1.026},
          { "N1",  "N",  2.050,  0.860,  0.962},
          { "C2",  "C",  2.040,  2.170,  1.127},
          { "N3",  "N",  1.920,  3.240,  0.841},
          { "C4",  "C",  1.800,  3.330,  0.431},
          { "H5'1",  "H",  -1.178,  8.686,  0.271},
          { "H5'2",  "H",  -1.450,  6.919,  0.250},
          { "H4'",  "H",  -0.582,  8.160,  0.531},
          { "H1'",  "H",  2.073,  5.804,  0.656},
          { "H2'1",  "H",  2.353,  7.073,  0.228},
          { "H2'2",  "H",  3.344,  7.146,  0.439},
          { "H3'",  "H",  1.510,  9.057,  0.362},
          { "H8",  "H",  1.469,  5.889,  -0.013},
          { "H61",  "H",  1.850,  2.827,  -0.968},
          { "H62",  "H",  2.039,  1.740,  -1.611},
          { "H2",  "H",  2.147,  2.724,  1.510},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 30, 17, 1 },
          { 31, 17, 1 },
          { 17, 16, 1 },
          { 3, 1, 2 },
          { 14, 15, 1 },
          { 14, 13, 2 },
          { 16, 15, 2 },
          { 16, 18, 1 },
          { 29, 13, 1 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 15, 21, 1 },
          { 13, 12, 1 },
          { 18, 19, 2 },
          { 4, 5, 1 },
          { 21, 12, 1 },
          { 21, 20, 2 },
          { 12, 11, 1 },
          { 26, 10, 1 },
          { 23, 5, 1 },
          { 19, 20, 1 },
          { 19, 32, 1 },
          { 5, 22, 1 },
          { 5, 6, 1 },
          { 10, 11, 1 },
          { 10, 27, 1 },
          { 10, 8, 1 },
          { 11, 7, 1 },
          { 11, 25, 1 },
          { 7, 6, 1 },
          { 28, 8, 1 },
          { 6, 8, 1 },
          { 6, 24, 1 },
          { 8, 9, 1 },
          { 0, 0, 0}
        }
      },
      { 'C', "DC", // C DC  Cytosine C -> G
        {
          { "P",  "P",  0.000,  8.910,  0.000},
          { "O1P",  "O",  0.220,  10.200,  0.072},
          { "O2P",  "O",  0.790,  8.820,  -0.141},
          { "O5'",  "O",  0.250,  7.730,  0.126},
          { "C5'",  "C",  -0.690,  7.700,  0.269},
          { "C4'",  "C",  0.040,  7.590,  0.442},
          { "O4'",  "O",  0.250,  6.220,  0.510},
          { "C3'",  "C",  1.440,  8.200,  0.442},
          { "O3'",  "O",  1.830,  8.750,  0.590},
          { "C2'",  "C",  2.320,  7.040,  0.384},
          { "C1'",  "C",  1.610,  5.860,  0.485},
          { "N1",  "N",  1.660,  4.629,  0.325},
          { "C2",  "C",  1.810,  3.400,  0.485},
          { "O2",  "O",  1.900,  3.690,  0.826},
          { "N3",  "N",  1.860,  2.311,  0.110},
          { "C4",  "C",  1.760,  2.940,  -0.258},
          { "N4",  "N",  1.810,  2.760,  -0.723},
          { "C5",  "C",  1.610,  4.350,  -0.206},
          { "C6",  "C",  1.560,  4.990,  0.052},
          { "H5'1",  "H",  -1.180,  8.684,  0.271},
          { "H5'2",  "H",  -1.449,  6.916,  0.250},
          { "H4'",  "H",  -0.582,  8.160,  0.531},
          { "H1'",  "H",  2.073,  5.803,  0.656},
          { "H2'1",  "H",  2.352,  7.073,  0.228},
          { "H2'2",  "H",  3.344,  7.145,  0.439},
          { "H3'",  "H",  1.510,  9.056,  0.362},
          { "H6",  "H",  1.438,  6.067,  0.034},
          { "H5",  "H",  1.540,  5.087,  -0.376},
          { "H41",  "H",  1.928,  2.018,  -1.023},
          { "H42",  "H",  1.731,  3.718,  -0.838},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 30, 17, 1 },
          { 28, 18, 1 },
          { 17, 29, 1 },
          { 17, 16, 1 },
          { 3, 1, 2 },
          { 18, 16, 1 },
          { 18, 19, 2 },
          { 16, 15, 2 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 27, 19, 1 },
          { 15, 13, 1 },
          { 19, 12, 1 },
          { 4, 5, 1 },
          { 12, 13, 1 },
          { 12, 11, 1 },
          { 13, 14, 2 },
          { 24, 10, 1 },
          { 21, 5, 1 },
          { 5, 20, 1 },
          { 5, 6, 1 },
          { 10, 11, 1 },
          { 10, 25, 1 },
          { 10, 8, 1 },
          { 11, 7, 1 },
          { 11, 23, 1 },
          { 7, 6, 1 },
          { 26, 8, 1 },
          { 6, 8, 1 },
          { 6, 22, 1 },
          { 8, 9, 1 },
          { 0, 0, 0}
        }
      },
      { 'G', "DG", // G DG  Guanine G -> C
        {
          { "P",  "P",  0.000,  8.910,  0.000},
          { "O1P",  "O",  0.220,  10.199,  0.072},
          { "O2P",  "O",  0.790,  8.819,  -0.141},
          { "O5'",  "O",  0.250,  7.730,  0.126},
          { "C5'",  "C",  -0.690,  7.700,  0.269},
          { "C4'",  "C",  0.040,  7.590,  0.442},
          { "O4'",  "O",  0.250,  6.220,  0.510},
          { "C3'",  "C",  1.440,  8.200,  0.442},
          { "O3'",  "O",  1.830,  8.750,  0.590},
          { "C2'",  "C",  2.320,  7.040,  0.384},
          { "C1'",  "C",  1.610,  5.860,  0.485},
          { "N9",  "N",  1.660,  4.630,  0.325},
          { "C8",  "C",  1.580,  4.819,  0.035},
          { "N7",  "N",  1.660,  3.920,  -0.183},
          { "C5",  "C",  1.800,  2.700,  0.021},
          { "C6",  "C",  1.930,  1.390,  -0.246},
          { "O6",  "O",  1.950,  1.709,  -1.037},
          { "N1",  "N",  2.050,  0.920,  1.001},
          { "C2",  "C",  2.050,  2.280,  1.161},
          { "N2",  "N",  2.180,  3.010,  1.588},
          { "N3",  "N",  1.920,  3.290,  0.847},
          { "C4",  "C",  1.800,  3.330,  0.435},
          { "H5'1",  "H",  -1.179,  8.685,  0.271},
          { "H5'2",  "H",  -1.449,  6.917,  0.250},
          { "H4'",  "H",  -0.582,  8.159,  0.531},
          { "H1'",  "H",  2.073,  5.803,  0.657},
          { "H2'1",  "H",  2.352,  7.072,  0.228},
          { "H2'2",  "H",  3.344,  7.145,  0.439},
          { "H3'",  "H",  1.510,  9.056,  0.362},
          { "H8",  "H",  1.461,  5.865,  -0.018},
          { "H1",  "H",  2.143,  0.769,  2.251},
          { "H21",  "H",  2.180,  4.016,  1.524},
          { "H22",  "H",  2.281,  2.908,  1.935},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 17, 16, 2 },
          { 3, 1, 2 },
          { 14, 15, 1 },
          { 14, 13, 2 },
          { 16, 15, 1 },
          { 16, 18, 1 },
          { 30, 13, 1 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 15, 22, 2 },
          { 13, 12, 1 },
          { 31, 18, 1 },
          { 18, 19, 1 },
          { 4, 5, 1 },
          { 22, 12, 1 },
          { 22, 21, 1 },
          { 12, 11, 1 },
          { 27, 10, 1 },
          { 24, 5, 1 },
          { 5, 23, 1 },
          { 5, 6, 1 },
          { 19, 21, 2 },
          { 19, 20, 1 },
          { 10, 11, 1 },
          { 10, 28, 1 },
          { 10, 8, 1 },
          { 33, 20, 1 },
          { 11, 7, 1 },
          { 11, 26, 1 },
          { 20, 32, 1 },
          { 7, 6, 1 },
          { 29, 8, 1 },
          { 6, 8, 1 },
          { 6, 25, 1 },
          { 8, 9, 1 },
          { 0, 0, 0}
        }
      },
      { 'T', "DT", // T DT  Thymine T -> A
        {
          { "P",  "P",  0.000,  8.910,  0.000},
          { "O1P",  "O",  0.220,  10.200,  0.072},
          { "O2P",  "O",  0.790,  8.820,  -0.141},
          { "O5'",  "O",  0.250,  7.730,  0.126},
          { "C5'",  "C",  -0.690,  7.699,  0.269},
          { "C4'",  "C",  0.040,  7.590,  0.442},
          { "O4'",  "O",  0.250,  6.220,  0.510},
          { "C3'",  "C",  1.440,  8.200,  0.442},
          { "O3'",  "O",  1.830,  8.750,  0.590},
          { "C2'",  "C",  2.320,  7.040,  0.384},
          { "C1'",  "C",  1.610,  5.860,  0.485},
          { "N1",  "N",  1.660,  4.629,  0.325},
          { "C2",  "C",  1.810,  3.420,  0.487},
          { "O2",  "O",  1.900,  3.640,  0.827},
          { "N3",  "N",  1.850,  2.360,  0.174},
          { "C4",  "C",  1.760,  2.980,  -0.291},
          { "O4",  "O",  1.810,  2.820,  -0.717},
          { "C5",  "C",  1.610,  4.380,  -0.204},
          { "C5M",  "C",  1.500,  5.400,  -0.429},
          { "C6",  "C",  1.560,  5.010,  0.051},
          { "H5'1",  "H",  -1.180,  8.685,  0.271},
          { "H5'2",  "H",  -1.449,  6.916,  0.250},
          { "H4'",  "H",  -0.582,  8.159,  0.531},
          { "H1'",  "H",  2.073,  5.802,  0.656},
          { "H2'1",  "H",  2.352,  7.073,  0.228},
          { "H2'2",  "H",  3.344,  7.145,  0.439},
          { "H3'",  "H",  1.510,  9.056,  0.362},
          { "H6",  "H",  1.439,  6.089,  0.034},
          { "H51",  "H",  2.469,  5.917,  -0.442},
          { "H52",  "H",  0.732,  6.117,  -0.373},
          { "H53",  "H",  1.220,  5.101,  -0.624},
          { "H3",  "H",  1.956,  1.414,  0.342},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 31, 19, 1 },
          { 29, 19, 1 },
          { 19, 30, 1 },
          { 19, 18, 1 },
          { 17, 16, 2 },
          { 3, 1, 2 },
          { 18, 16, 1 },
          { 18, 20, 2 },
          { 16, 15, 1 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 28, 20, 1 },
          { 20, 12, 1 },
          { 15, 32, 1 },
          { 15, 13, 1 },
          { 4, 5, 1 },
          { 12, 13, 1 },
          { 12, 11, 1 },
          { 25, 10, 1 },
          { 13, 14, 2 },
          { 22, 5, 1 },
          { 5, 21, 1 },
          { 5, 6, 1 },
          { 10, 11, 1 },
          { 10, 26, 1 },
          { 10, 8, 1 },
          { 11, 7, 1 },
          { 11, 24, 1 },
          { 7, 6, 1 },
          { 27, 8, 1 },
          { 6, 8, 1 },
          { 6, 23, 1 },
          { 8, 9, 1 },
          { 0, 0, 0}
        }
      },
      { 'R', "GAX", // R GAX Purine (G or A) R -> Y
        {
          { "P",  "P",  0.000,  8.910,  0.000},
          { "O1P",  "O",  0.220,  10.199,  0.072},
          { "O2P",  "O",  0.790,  8.819,  -0.141},
          { "O5'",  "O",  0.250,  7.730,  0.126},
          { "C5'",  "C",  -0.690,  7.700,  0.269},
          { "C4'",  "C",  0.040,  7.590,  0.442},
          { "O4'",  "O",  0.250,  6.220,  0.510},
          { "C3'",  "C",  1.440,  8.200,  0.442},
          { "O3'",  "O",  1.830,  8.750,  0.590},
          { "C2'",  "C",  2.320,  7.040,  0.384},
          { "C1'",  "C",  1.610,  5.860,  0.485},
          { "N9",  "N",  1.660,  4.630,  0.325},
          { "C8",  "C",  1.580,  4.819,  0.035},
          { "N7",  "N",  1.660,  3.920,  -0.183},
          { "C5",  "C",  1.800,  2.700,  0.021},
          { "C6",  "C",  1.930,  1.390,  -0.246},
          { "N1",  "N",  2.050,  0.920,  1.001},
          { "C2",  "C",  2.050,  2.280,  1.161},
          { "N2",  "N",  2.180,  3.010,  1.588},
          { "N3",  "N",  1.920,  3.290,  0.847},
          { "C4",  "C",  1.800,  3.330,  0.435},
          { "H5'1",  "H",  -1.179,  8.685,  0.271},
          { "H5'2",  "H",  -1.449,  6.917,  0.250},
          { "H4'",  "H",  -0.582,  8.159,  0.531},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 3, 1, 2 },
          { 14, 15, 1 },
          { 14, 13, 2 },
          { 16, 15, 2 },
          { 16, 17, 1 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 15, 21, 1 },
          { 13, 12, 1 },
          { 17, 18, 2 },
          { 4, 5, 1 },
          { 21, 12, 1 },
          { 21, 20, 2 },
          { 12, 11, 1 },
          { 23, 5, 1 },
          { 5, 22, 1 },
          { 5, 6, 1 },
          { 18, 20, 1 },
          { 18, 19, 1 },
          { 10, 11, 1 },
          { 10, 8, 1 },
          { 11, 7, 1 },
          { 7, 6, 1 },
          { 6, 8, 1 },
          { 6, 24, 1 },
          { 8, 9, 1 },
          { 0, 0, 0}
        }
      },
      { 'Y', "TCX", // Y TCX Pyrimidine (T or C) Y -> R
        {
          { "P",  "P",  0.000,  8.910,  0.000},
          { "O1P",  "O",  0.220,  10.200,  0.072},
          { "O2P",  "O",  0.790,  8.820,  -0.141},
          { "O5'",  "O",  0.250,  7.730,  0.126},
          { "C5'",  "C",  -0.690,  7.699,  0.269},
          { "C4'",  "C",  0.040,  7.590,  0.442},
          { "O4'",  "O",  0.250,  6.220,  0.510},
          { "C3'",  "C",  1.440,  8.200,  0.442},
          { "O3'",  "O",  1.830,  8.750,  0.590},
          { "C2'",  "C",  2.320,  7.040,  0.384},
          { "C1'",  "C",  1.610,  5.860,  0.485},
          { "N1",  "N",  1.660,  4.629,  0.325},
          { "C2",  "C",  1.810,  3.420,  0.487},
          { "O2",  "O",  1.900,  3.640,  0.827},
          { "N3",  "N",  1.850,  2.360,  0.174},
          { "C4",  "C",  1.760,  2.980,  -0.291},
          { "C5",  "C",  1.610,  4.380,  -0.204},
          { "C6",  "C",  1.560,  5.010,  0.051},
          { "H5'1",  "H",  -1.180,  8.685,  0.271},
          { "H5'2",  "H",  -1.449,  6.916,  0.250},
          { "H4'",  "H",  -0.582,  8.159,  0.531},
          { "H1'",  "H",  2.073,  5.802,  0.656},
          { "H2'1",  "H",  2.352,  7.073,  0.228},
          { "H2'2",  "H",  3.344,  7.145,  0.439},
          { "H3'",  "H",  1.510,  9.056,  0.362},
          { "H6",  "H",  1.439,  6.089,  0.034},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 3, 1, 2 },
          { 17, 16, 1 },
          { 17, 18, 2 },
          { 16, 15, 2 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 26, 18, 1 },
          { 18, 12, 1 },
          { 15, 13, 1 },
          { 4, 5, 1 },
          { 12, 13, 1 },
          { 12, 11, 1 },
          { 23, 10, 1 },
          { 13, 14, 2 },
          { 20, 5, 1 },
          { 5, 19, 1 },
          { 5, 6, 1 },
          { 10, 11, 1 },
          { 10, 24, 1 },
          { 10, 8, 1 },
          { 11, 7, 1 },
          { 11, 22, 1 },
          { 7, 6, 1 },
          { 25, 8, 1 },
          { 6, 8, 1 },
          { 6, 21, 1 },
          { 8, 9, 1 },
          { 0, 0, 0}
        }
      },
      { 'K', "GTX", // K GTX Keto (G or T) K -> M
        {
          { "P",  "P",  0.000,  8.910,  0.000},
          { "O1P",  "O",  0.220,  10.199,  0.072},
          { "O2P",  "O",  0.790,  8.819,  -0.141},
          { "O5'",  "O",  0.250,  7.730,  0.126},
          { "C5'",  "C",  -0.690,  7.700,  0.269},
          { "C4'",  "C",  0.040,  7.590,  0.442},
          { "O4'",  "O",  0.250,  6.220,  0.510},
          { "C3'",  "C",  1.440,  8.200,  0.442},
          { "O3'",  "O",  1.830,  8.750,  0.590},
          { "C2'",  "C",  2.320,  7.040,  0.384},
          { "C1'",  "C",  1.610,  5.860,  0.485},
          { "N9",  "N",  1.660,  4.630,  0.325},
          { "H5'1",  "H",  -1.179,  8.685,  0.271},
          { "H5'2",  "H",  -1.449,  6.917,  0.250},
          { "H4'",  "H",  -0.582,  8.159,  0.531},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 3, 1, 2 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 4, 5, 1 },
          { 12, 11, 1 },
          { 14, 5, 1 },
          { 5, 13, 1 },
          { 5, 6, 1 },
          { 10, 11, 1 },
          { 10, 8, 1 },
          { 11, 7, 1 },
          { 7, 6, 1 },
          { 6, 8, 1 },
          { 6, 15, 1 },
          { 8, 9, 1 },
          { 0, 0, 0}
        }
      },
      { 'M', "ACX", // M ACX Amino (A or C) M -> K
        {
          { "P",  "P",  0.000,  8.910,  0.000},
          { "O1P",  "O",  0.220,  10.201,  0.072},
          { "O2P",  "O",  0.790,  8.821,  -0.141},
          { "O5'",  "O",  0.250,  7.730,  0.126},
          { "C5'",  "C",  -0.690,  7.701,  0.269},
          { "C4'",  "C",  0.040,  7.590,  0.442},
          { "O4'",  "O",  0.250,  6.221,  0.510},
          { "C3'",  "C",  1.440,  8.201,  0.442},
          { "O3'",  "O",  1.830,  8.750,  0.590},
          { "C2'",  "C",  2.320,  7.040,  0.384},
          { "C1'",  "C",  1.610,  5.860,  0.485},
          { "N9",  "N",  1.660,  4.630,  0.325},
          { "H5'1",  "H",  -1.178,  8.686,  0.271},
          { "H5'2",  "H",  -1.450,  6.919,  0.250},
          { "H4'",  "H",  -0.582,  8.160,  0.531},
          { "H1'",  "H",  2.073,  5.804,  0.656},
          { "H2'1",  "H",  2.353,  7.073,  0.228},
          { "H2'2",  "H",  3.344,  7.146,  0.439},
          { "H3'",  "H",  1.510,  9.057,  0.362},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 3, 1, 2 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 4, 5, 1 },
          { 12, 11, 1 },
          { 17, 10, 1 },
          { 14, 5, 1 },
          { 5, 13, 1 },
          { 5, 6, 1 },
          { 10, 11, 1 },
          { 10, 18, 1 },
          { 10, 8, 1 },
          { 11, 7, 1 },
          { 11, 16, 1 },
          { 7, 6, 1 },
          { 19, 8, 1 },
          { 6, 8, 1 },
          { 6, 15, 1 },
          { 8, 9, 1 },
          { 0, 0, 0}
        }
      },
      { 'S', "GCX", // S GCX Strong (G or C) S -> W
        {
          { "P",  "P",  0.000,  8.910,  0.000},
          { "O1P",  "O",  0.220,  10.199,  0.072},
          { "O2P",  "O",  0.790,  8.819,  -0.141},
          { "O5'",  "O",  0.250,  7.730,  0.126},
          { "C5'",  "C",  -0.690,  7.700,  0.269},
          { "C4'",  "C",  0.040,  7.590,  0.442},
          { "O4'",  "O",  0.250,  6.220,  0.510},
          { "C3'",  "C",  1.440,  8.200,  0.442},
          { "O3'",  "O",  1.830,  8.750,  0.590},
          { "C2'",  "C",  2.320,  7.040,  0.384},
          { "C1'",  "C",  1.610,  5.860,  0.485},
          { "N9",  "N",  1.660,  4.630,  0.325},
          { "H5'1",  "H",  -1.179,  8.685,  0.271},
          { "H5'2",  "H",  -1.449,  6.917,  0.250},
          { "H4'",  "H",  -0.582,  8.159,  0.531},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 3, 1, 2 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 4, 5, 1 },
          { 12, 11, 1 },
          { 14, 5, 1 },
          { 5, 13, 1 },
          { 5, 6, 1 },
          { 10, 11, 1 },
          { 10, 8, 1 },
          { 11, 7, 1 },
          { 7, 6, 1 },
          { 6, 8, 1 },
          { 6, 15, 1 },
          { 8, 9, 1 },
          { 0, 0, 0}
        }
      },
      { 'W', "ATX", // W ATX Weak (A or T) W -> S
        {
          { "P",  "P",  0.000,  8.910,  0.000},
          { "O1P",  "O",  0.220,  10.201,  0.072},
          { "O2P",  "O",  0.790,  8.821,  -0.141},
          { "O5'",  "O",  0.250,  7.730,  0.126},
          { "C5'",  "C",  -0.690,  7.701,  0.269},
          { "C4'",  "C",  0.040,  7.590,  0.442},
          { "O4'",  "O",  0.250,  6.221,  0.510},
          { "C3'",  "C",  1.440,  8.201,  0.442},
          { "O3'",  "O",  1.830,  8.750,  0.590},
          { "C2'",  "C",  2.320,  7.040,  0.384},
          { "C1'",  "C",  1.610,  5.860,  0.485},
          { "N9",  "N",  1.660,  4.630,  0.325},
          { "H5'1",  "H",  -1.178,  8.686,  0.271},
          { "H5'2",  "H",  -1.450,  6.919,  0.250},
          { "H4'",  "H",  -0.582,  8.160,  0.531},
          { "H1'",  "H",  2.073,  5.804,  0.656},
          { "H2'1",  "H",  2.353,  7.073,  0.228},
          { "H2'2",  "H",  3.344,  7.146,  0.439},
          { "H3'",  "H",  1.510,  9.057,  0.362},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 3, 1, 2 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 4, 5, 1 },
          { 12, 11, 1 },
          { 17, 10, 1 },
          { 14, 5, 1 },
          { 5, 13, 1 },
          { 5, 6, 1 },
          { 10, 11, 1 },
          { 10, 18, 1 },
          { 10, 8, 1 },
          { 11, 7, 1 },
          { 11, 16, 1 },
          { 7, 6, 1 },
          { 19, 8, 1 },
          { 6, 8, 1 },
          { 6, 15, 1 },
          { 8, 9, 1 },
          { 0, 0, 0}
        }
      },
      { 'B', "GTC", // B GTC (G T C) B -> V
        {
          { "P",  "P",  0.000,  8.910,  0.000},
          { "O1P",  "O",  0.220,  10.199,  0.072},
          { "O2P",  "O",  0.790,  8.819,  -0.141},
          { "O5'",  "O",  0.250,  7.730,  0.126},
          { "C5'",  "C",  -0.690,  7.700,  0.269},
          { "C4'",  "C",  0.040,  7.590,  0.442},
          { "O4'",  "O",  0.250,  6.220,  0.510},
          { "C3'",  "C",  1.440,  8.200,  0.442},
          { "O3'",  "O",  1.830,  8.750,  0.590},
          { "C2'",  "C",  2.320,  7.040,  0.384},
          { "C1'",  "C",  1.610,  5.860,  0.485},
          { "N9",  "N",  1.660,  4.630,  0.325},
          { "H5'1",  "H",  -1.179,  8.685,  0.271},
          { "H5'2",  "H",  -1.449,  6.917,  0.250},
          { "H4'",  "H",  -0.582,  8.159,  0.531},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 3, 1, 2 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 4, 5, 1 },
          { 12, 11, 1 },
          { 14, 5, 1 },
          { 5, 13, 1 },
          { 5, 6, 1 },
          { 10, 11, 1 },
          { 10, 8, 1 },
          { 11, 7, 1 },
          { 7, 6, 1 },
          { 6, 8, 1 },
          { 6, 15, 1 },
          { 8, 9, 1 },
          { 0, 0, 0}
        }
      },
      { 'D', "GAT", // D GAT (G A T) D -> H
        {
          { "P",  "P",  0.000,  8.910,  0.000},
          { "O1P",  "O",  0.220,  10.199,  0.072},
          { "O2P",  "O",  0.790,  8.819,  -0.141},
          { "O5'",  "O",  0.250,  7.730,  0.126},
          { "C5'",  "C",  -0.690,  7.700,  0.269},
          { "C4'",  "C",  0.040,  7.590,  0.442},
          { "O4'",  "O",  0.250,  6.220,  0.510},
          { "C3'",  "C",  1.440,  8.200,  0.442},
          { "O3'",  "O",  1.830,  8.750,  0.590},
          { "C2'",  "C",  2.320,  7.040,  0.384},
          { "C1'",  "C",  1.610,  5.860,  0.485},
          { "N9",  "N",  1.660,  4.630,  0.325},
          { "H5'1",  "H",  -1.179,  8.685,  0.271},
          { "H5'2",  "H",  -1.449,  6.917,  0.250},
          { "H4'",  "H",  -0.582,  8.159,  0.531},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 3, 1, 2 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 4, 5, 1 },
          { 12, 11, 1 },
          { 14, 5, 1 },
          { 5, 13, 1 },
          { 5, 6, 1 },
          { 10, 11, 1 },
          { 10, 8, 1 },
          { 11, 7, 1 },
          { 7, 6, 1 },
          { 6, 8, 1 },
          { 6, 15, 1 },
          { 8, 9, 1 },
          { 0, 0, 0}
        }
      },
      { 'H', "ACT", // H ACT (A C T) H -> D
        {
          { "P",  "P",  0.000,  8.910,  0.000},
          { "O1P",  "O",  0.220,  10.201,  0.072},
          { "O2P",  "O",  0.790,  8.821,  -0.141},
          { "O5'",  "O",  0.250,  7.730,  0.126},
          { "C5'",  "C",  -0.690,  7.701,  0.269},
          { "C4'",  "C",  0.040,  7.590,  0.442},
          { "O4'",  "O",  0.250,  6.221,  0.510},
          { "C3'",  "C",  1.440,  8.201,  0.442},
          { "O3'",  "O",  1.830,  8.750,  0.590},
          { "C2'",  "C",  2.320,  7.040,  0.384},
          { "C1'",  "C",  1.610,  5.860,  0.485},
          { "N9",  "N",  1.660,  4.630,  0.325},
          { "H5'1",  "H",  -1.178,  8.686,  0.271},
          { "H5'2",  "H",  -1.450,  6.919,  0.250},
          { "H4'",  "H",  -0.582,  8.160,  0.531},
          { "H1'",  "H",  2.073,  5.804,  0.656},
          { "H2'1",  "H",  2.353,  7.073,  0.228},
          { "H2'2",  "H",  3.344,  7.146,  0.439},
          { "H3'",  "H",  1.510,  9.057,  0.362},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 3, 1, 2 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 4, 5, 1 },
          { 12, 11, 1 },
          { 17, 10, 1 },
          { 14, 5, 1 },
          { 5, 13, 1 },
          { 5, 6, 1 },
          { 10, 11, 1 },
          { 10, 18, 1 },
          { 10, 8, 1 },
          { 11, 7, 1 },
          { 11, 16, 1 },
          { 7, 6, 1 },
          { 19, 8, 1 },
          { 6, 8, 1 },
          { 6, 15, 1 },
          { 8, 9, 1 },
          { 0, 0, 0}
        }
      },
      { 'V', "GCA", // V GCA (G C A) V -> B
        {
          { "P",  "P",  0.000,  8.910,  0.000},
          { "O1P",  "O",  0.220,  10.199,  0.072},
          { "O2P",  "O",  0.790,  8.819,  -0.141},
          { "O5'",  "O",  0.250,  7.730,  0.126},
          { "C5'",  "C",  -0.690,  7.700,  0.269},
          { "C4'",  "C",  0.040,  7.590,  0.442},
          { "O4'",  "O",  0.250,  6.220,  0.510},
          { "C3'",  "C",  1.440,  8.200,  0.442},
          { "O3'",  "O",  1.830,  8.750,  0.590},
          { "C2'",  "C",  2.320,  7.040,  0.384},
          { "C1'",  "C",  1.610,  5.860,  0.485},
          { "N9",  "N",  1.660,  4.630,  0.325},
          { "H5'1",  "H",  -1.179,  8.685,  0.271},
          { "H5'2",  "H",  -1.449,  6.917,  0.250},
          { "H4'",  "H",  -0.582,  8.159,  0.531},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 3, 1, 2 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 4, 5, 1 },
          { 12, 11, 1 },
          { 14, 5, 1 },
          { 5, 13, 1 },
          { 5, 6, 1 },
          { 10, 11, 1 },
          { 10, 8, 1 },
          { 11, 7, 1 },
          { 7, 6, 1 },
          { 6, 8, 1 },
          { 6, 15, 1 },
          { 8, 9, 1 },
          { 0, 0, 0}
        }
      }
    };
  ResidueRecord DNAPairResidues[IUPAC_DNA_max] =
    {
      { 0, "", // DNA Pair Start
        {
          { "HTER",  "H",  5.721,  9.511,  -2.837},
          { "OXT",  "O",  5.601,  8.904,  -2.915},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 2, 1, 1},
          { 0, 0, 0}
        }
      },
      { 0, "", // DNA Pair End
        {
          { "HCAP",  "H",  1.613,  9.365,  2.752},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 0, 0, 0}
        }
      },
      { 'N', "XNA", // N XNA Pair Unknown Nucleic Acid N -> N
        {
          { "P",  "P",  4.160,  8.910,  -2.960},
          { "O1P",  "O",  3.940,  10.200,  -3.032},
          { "O2P",  "O",  3.370,  8.820,  -2.819},
          { "O5'",  "O",  3.910,  7.730,  -3.086},
          { "C5'",  "C",  4.850,  7.701,  3.054},
          { "C4'",  "C",  4.120,  7.590,  2.881},
          { "O4'",  "O",  3.910,  6.220,  2.813},
          { "C3'",  "C",  2.720,  8.200,  2.882},
          { "O3'",  "O",  2.330,  8.750,  2.733},
          { "C2'",  "C",  1.840,  7.040,  2.939},
          { "C1'",  "C",  2.550,  5.860,  2.838},
          { "N9",  "N",  2.500,  4.630,  2.998},
          { "H5'1",  "H",  5.339,  8.686,  3.052},
          { "H5'2",  "H",  5.611,  6.918,  3.073},
          { "H4'",  "H",  4.743,  8.159,  2.792},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 3, 1, 2 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 4, 5, 1 },
          { 12, 11, 1 },
          { 14, 5, 1 },
          { 5, 13, 1 },
          { 5, 6, 1 },
          { 10, 11, 1 },
          { 10, 8, 1 },
          { 11, 7, 1 },
          { 7, 6, 1 },
          { 6, 8, 1 },
          { 6, 15, 1 },
          { 8, 9, 1 },
          { 0, 0, 0}
        }
      },
      { 'T', "DT", // T DT  Pair Thymine T -> A
        {
          { "P",  "P",  4.160,  8.910,  -2.960},
          { "O1P",  "O",  3.940,  10.200,  -3.032},
          { "O2P",  "O",  3.370,  8.820,  -2.819},
          { "O5'",  "O",  3.910,  7.730,  -3.086},
          { "C5'",  "C",  4.850,  7.699,  3.054},
          { "C4'",  "C",  4.120,  7.590,  2.881},
          { "O4'",  "O",  3.910,  6.220,  2.813},
          { "C3'",  "C",  2.720,  8.199,  2.881},
          { "O3'",  "O",  2.330,  8.750,  2.733},
          { "C2'",  "C",  1.840,  7.039,  2.939},
          { "C1'",  "C",  2.550,  5.859,  2.838},
          { "N1",  "N",  2.500,  4.629,  2.998},
          { "C2",  "C",  2.350,  3.419,  2.836},
          { "O2",  "O",  2.260,  3.640,  2.496},
          { "N3",  "N",  2.310,  2.360,  -3.135},
          { "C4",  "C",  2.400,  2.979,  -2.669},
          { "O4",  "O",  2.350,  2.820,  -2.243},
          { "C5",  "C",  2.550,  4.379,  -2.756},
          { "C5M",  "C",  2.660,  5.399,  -2.531},
          { "C6",  "C",  2.600,  5.010,  -3.011},
          { "H5'1",  "H",  5.341,  8.684,  3.052},
          { "H5'2",  "H",  5.610,  6.916,  3.073},
          { "H4'",  "H",  4.743,  8.159,  2.792},
          { "H1'",  "H",  2.088,  5.802,  2.667},
          { "H2'1",  "H",  1.809,  7.072,  3.095},
          { "H2'2",  "H",  0.817,  7.144,  2.884},
          { "H3'",  "H",  2.651,  9.055,  2.961},
          { "H6",  "H",  2.722,  6.088,  -2.994},
          { "H51",  "H",  1.692,  5.916,  -2.518},
          { "H52",  "H",  3.429,  6.117,  -2.587},
          { "H53",  "H",  2.941,  5.100,  -2.336},
          { "H3",  "H",  2.205,  1.414,  2.980},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 31, 19, 1 },
          { 29, 19, 1 },
          { 19, 30, 1 },
          { 19, 18, 1 },
          { 17, 16, 2 },
          { 3, 1, 2 },
          { 18, 16, 1 },
          { 18, 20, 2 },
          { 16, 15, 1 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 28, 20, 1 },
          { 20, 12, 1 },
          { 15, 32, 1 },
          { 15, 13, 1 },
          { 4, 5, 1 },
          { 12, 13, 1 },
          { 12, 11, 1 },
          { 25, 10, 1 },
          { 13, 14, 2 },
          { 22, 5, 1 },
          { 5, 21, 1 },
          { 5, 6, 1 },
          { 10, 11, 1 },
          { 10, 26, 1 },
          { 10, 8, 1 },
          { 11, 7, 1 },
          { 11, 24, 1 },
          { 7, 6, 1 },
          { 27, 8, 1 },
          { 6, 8, 1 },
          { 6, 23, 1 },
          { 8, 9, 1 },
          { 0, 0, 0}
        }
      },
      { 'G', "DG", // G DG  Pair Guanine G -> C
        {
          { "P",  "P",  4.160,  8.910,  -2.960},
          { "O1P",  "O",  3.940,  10.200,  -3.032},
          { "O2P",  "O",  3.370,  8.820,  -2.819},
          { "O5'",  "O",  3.910,  7.730,  -3.086},
          { "C5'",  "C",  4.850,  7.701,  3.054},
          { "C4'",  "C",  4.120,  7.590,  2.881},
          { "O4'",  "O",  3.910,  6.220,  2.813},
          { "C3'",  "C",  2.720,  8.200,  2.882},
          { "O3'",  "O",  2.330,  8.750,  2.733},
          { "C2'",  "C",  1.840,  7.040,  2.939},
          { "C1'",  "C",  2.550,  5.860,  2.838},
          { "N9",  "N",  2.500,  4.630,  2.998},
          { "C8",  "C",  2.580,  4.820,  -2.995},
          { "N7",  "N",  2.500,  3.920,  -2.777},
          { "C5",  "C",  2.360,  2.700,  -2.981},
          { "C6",  "C",  2.230,  1.389,  -2.714},
          { "O6",  "O",  2.210,  1.710,  -1.923},
          { "N1",  "N",  2.110,  0.920,  2.323},
          { "C2",  "C",  2.110,  2.281,  2.162},
          { "N2",  "N",  1.980,  3.010,  1.735},
          { "N3",  "N",  2.240,  3.290,  2.477},
          { "C4",  "C",  2.360,  3.330,  2.889},
          { "H5'1",  "H",  5.339,  8.686,  3.052},
          { "H5'2",  "H",  5.611,  6.918,  3.073},
          { "H4'",  "H",  4.743,  8.159,  2.792},
          { "H1'",  "H",  2.088,  5.804,  2.667},
          { "H2'1",  "H",  1.808,  7.073,  3.095},
          { "H2'2",  "H",  0.817,  7.146,  2.884},
          { "H3'",  "H",  2.651,  9.057,  2.961},
          { "H8",  "H",  2.700,  5.866,  -2.942},
          { "H1",  "H",  2.018,  0.768,  1.073},
          { "H21",  "H",  1.981,  4.016,  1.799},
          { "H22",  "H",  1.880,  2.906,  1.388},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 3, 1, 2 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 17, 16, 2 },
          { 14, 13, 2 },
          { 14, 15, 1 },
          { 30, 13, 1 },
          { 13, 12, 1 },
          { 16, 15, 1 },
          { 16, 18, 1 },
          { 15, 22, 2 },
          { 4, 5, 1 },
          { 27, 10, 1 },
          { 24, 5, 1 },
          { 12, 22, 1 },
          { 12, 11, 1 },
          { 18, 31, 1 },
          { 18, 19, 1 },
          { 5, 23, 1 },
          { 5, 6, 1 },
          { 22, 21, 1 },
          { 10, 11, 1 },
          { 10, 28, 1 },
          { 10, 8, 1 },
          { 29, 8, 1 },
          { 11, 7, 1 },
          { 11, 26, 1 },
          { 19, 21, 2 },
          { 19, 20, 1 },
          { 6, 7, 1 },
          { 6, 8, 1 },
          { 6, 25, 1 },
          { 8, 9, 1 },
          { 33, 20, 1 },
          { 20, 32, 1 },
          { 0, 0, 0}
        }
      },
      { 'C', "DC", // C DC  Pair Cytosine C -> G
        {
          { "P",  "P",  4.160,  8.910,  -2.960},
          { "O1P",  "O",  3.940,  10.200,  -3.032},
          { "O2P",  "O",  3.370,  8.820,  -2.819},
          { "O5'",  "O",  3.910,  7.731,  -3.086},
          { "C5'",  "C",  4.850,  7.701,  3.054},
          { "C4'",  "C",  4.120,  7.591,  2.882},
          { "O4'",  "O",  3.910,  6.221,  2.813},
          { "C3'",  "C",  2.720,  8.200,  2.882},
          { "O3'",  "O",  2.330,  8.751,  2.733},
          { "C2'",  "C",  1.840,  7.040,  2.939},
          { "C1'",  "C",  2.550,  5.860,  2.838},
          { "N1",  "N",  2.500,  4.630,  2.998},
          { "C2",  "C",  2.350,  3.400,  2.838},
          { "O2",  "O",  2.260,  3.690,  2.497},
          { "N3",  "N",  2.300,  2.310,  -3.070},
          { "C4",  "C",  2.400,  2.940,  -2.702},
          { "N4",  "N",  2.350,  2.760,  -2.238},
          { "C5",  "C",  2.550,  4.350,  -2.754},
          { "C6",  "C",  2.600,  4.990,  -3.013},
          { "H5'1",  "H",  5.340,  8.685,  3.052},
          { "H5'2",  "H",  5.610,  6.918,  3.073},
          { "H4'",  "H",  4.743,  8.160,  2.792},
          { "H1'",  "H",  2.088,  5.804,  2.667},
          { "H2'1",  "H",  1.808,  7.073,  3.095},
          { "H2'2",  "H",  0.817,  7.146,  2.884},
          { "H3'",  "H",  2.651,  9.056,  2.961},
          { "H6",  "H",  2.723,  6.067,  -2.994},
          { "H5",  "H",  2.621,  5.088,  -2.584},
          { "H41",  "H",  2.233,  2.019,  -1.938},
          { "H42",  "H",  2.429,  3.718,  -2.122},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 30, 17, 1 },
          { 28, 18, 1 },
          { 17, 29, 1 },
          { 17, 16, 1 },
          { 3, 1, 2 },
          { 18, 16, 1 },
          { 18, 19, 2 },
          { 16, 15, 2 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 27, 19, 1 },
          { 15, 13, 1 },
          { 19, 12, 1 },
          { 4, 5, 1 },
          { 12, 13, 1 },
          { 12, 11, 1 },
          { 13, 14, 2 },
          { 24, 10, 1 },
          { 21, 5, 1 },
          { 5, 20, 1 },
          { 5, 6, 1 },
          { 10, 11, 1 },
          { 10, 25, 1 },
          { 10, 8, 1 },
          { 11, 7, 1 },
          { 11, 23, 1 },
          { 7, 6, 1 },
          { 26, 8, 1 },
          { 6, 8, 1 },
          { 6, 22, 1 },
          { 8, 9, 1 },
          { 0, 0, 0}
        }
      },
      { 'A', "DA", // A DA  Pair Adenine A -> T
        {
          { "P",  "P",  4.160,  8.910,  -2.960},
          { "O1P",  "O",  3.940,  10.200,  -3.032},
          { "O2P",  "O",  3.370,  8.820,  -2.819},
          { "O5'",  "O",  3.910,  7.730,  -3.086},
          { "C5'",  "C",  4.850,  7.701,  3.054},
          { "C4'",  "C",  4.120,  7.590,  2.881},
          { "O4'",  "O",  3.910,  6.220,  2.813},
          { "C3'",  "C",  2.720,  8.200,  2.882},
          { "O3'",  "O",  2.330,  8.750,  2.733},
          { "C2'",  "C",  1.840,  7.040,  2.939},
          { "C1'",  "C",  2.550,  5.860,  2.838},
          { "N9",  "N",  2.500,  4.630,  2.999},
          { "C8",  "C",  2.580,  4.840,  -2.998},
          { "N7",  "N",  2.510,  3.950,  -2.782},
          { "C5",  "C",  2.360,  2.740,  -2.981},
          { "C6",  "C",  2.230,  1.410,  -2.750},
          { "N6",  "N",  2.220,  1.810,  -1.934},
          { "N1",  "N",  2.110,  0.860,  2.361},
          { "C2",  "C",  2.120,  2.170,  2.196},
          { "N3",  "N",  2.240,  3.240,  2.482},
          { "C4",  "C",  2.360,  3.330,  2.892},
          { "H5'1",  "H",  5.339,  8.686,  3.052},
          { "H5'2",  "H",  5.611,  6.918,  3.074},
          { "H4'",  "H",  4.743,  8.159,  2.792},
          { "H1'",  "H",  2.088,  5.803,  2.667},
          { "H2'1",  "H",  1.808,  7.072,  3.095},
          { "H2'2",  "H",  0.817,  7.145,  2.884},
          { "H3'",  "H",  2.651,  9.056,  2.961},
          { "H8",  "H",  2.692,  5.889,  -2.947},
          { "H61",  "H",  2.311,  2.827,  -1.992},
          { "H62",  "H",  2.122,  1.741,  -1.349},
          { "H2",  "H",  2.014,  2.724,  1.813},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 30, 17, 1 },
          { 31, 17, 1 },
          { 17, 16, 1 },
          { 3, 1, 2 },
          { 14, 15, 1 },
          { 14, 13, 2 },
          { 16, 15, 2 },
          { 16, 18, 1 },
          { 29, 13, 1 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 15, 21, 1 },
          { 13, 12, 1 },
          { 18, 19, 2 },
          { 4, 5, 1 },
          { 21, 12, 1 },
          { 21, 20, 2 },
          { 12, 11, 1 },
          { 26, 10, 1 },
          { 23, 5, 1 },
          { 19, 20, 1 },
          { 19, 32, 1 },
          { 5, 22, 1 },
          { 5, 6, 1 },
          { 10, 11, 1 },
          { 10, 27, 1 },
          { 10, 8, 1 },
          { 11, 7, 1 },
          { 11, 25, 1 },
          { 7, 6, 1 },
          { 28, 8, 1 },
          { 6, 8, 1 },
          { 6, 24, 1 },
          { 8, 9, 1 },
          { 0, 0, 0}
        }
      },
      { 'Y', "TCX", // Y TCX Pair Pyrimidine (T or C) Y -> R
        {
          { "P",  "P",  4.160,  8.910,  -2.960},
          { "O1P",  "O",  3.940,  10.200,  -3.032},
          { "O2P",  "O",  3.370,  8.820,  -2.819},
          { "O5'",  "O",  3.910,  7.730,  -3.086},
          { "C5'",  "C",  4.850,  7.699,  3.054},
          { "C4'",  "C",  4.120,  7.590,  2.881},
          { "O4'",  "O",  3.910,  6.220,  2.813},
          { "C3'",  "C",  2.720,  8.199,  2.881},
          { "O3'",  "O",  2.330,  8.750,  2.733},
          { "C2'",  "C",  1.840,  7.039,  2.939},
          { "C1'",  "C",  2.550,  5.859,  2.838},
          { "N1",  "N",  2.500,  4.629,  2.998},
          { "C2",  "C",  2.350,  3.419,  2.836},
          { "O2",  "O",  2.260,  3.640,  2.496},
          { "N3",  "N",  2.310,  2.360,  -3.135},
          { "C4",  "C",  2.400,  2.979,  -2.669},
          { "C5",  "C",  2.550,  4.379,  -2.756},
          { "C6",  "C",  2.600,  5.010,  -3.011},
          { "H5'1",  "H",  5.341,  8.684,  3.052},
          { "H5'2",  "H",  5.610,  6.916,  3.073},
          { "H4'",  "H",  4.743,  8.159,  2.792},
          { "H1'",  "H",  2.088,  5.802,  2.667},
          { "H2'1",  "H",  1.809,  7.072,  3.095},
          { "H2'2",  "H",  0.817,  7.144,  2.884},
          { "H3'",  "H",  2.651,  9.055,  2.961},
          { "H6",  "H",  2.722,  6.088,  -2.994},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 3, 1, 2 },
          { 17, 16, 1 },
          { 17, 18, 2 },
          { 16, 15, 2 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 26, 18, 1 },
          { 18, 12, 1 },
          { 15, 13, 1 },
          { 4, 5, 1 },
          { 12, 13, 1 },
          { 12, 11, 1 },
          { 23, 10, 1 },
          { 13, 14, 2 },
          { 20, 5, 1 },
          { 5, 19, 1 },
          { 5, 6, 1 },
          { 10, 11, 1 },
          { 10, 24, 1 },
          { 10, 8, 1 },
          { 11, 7, 1 },
          { 11, 22, 1 },
          { 7, 6, 1 },
          { 25, 8, 1 },
          { 6, 8, 1 },
          { 6, 21, 1 },
          { 8, 9, 1 },
          { 0, 0, 0}
        }
      },
      { 'R', "GAX", // R GAX Pair Purine (G or A) R -> Y
        {
          { "P",  "P",  4.160,  8.910,  -2.960},
          { "O1P",  "O",  3.940,  10.200,  -3.032},
          { "O2P",  "O",  3.370,  8.820,  -2.819},
          { "O5'",  "O",  3.910,  7.730,  -3.086},
          { "C5'",  "C",  4.850,  7.701,  3.054},
          { "C4'",  "C",  4.120,  7.590,  2.881},
          { "O4'",  "O",  3.910,  6.220,  2.813},
          { "C3'",  "C",  2.720,  8.200,  2.882},
          { "O3'",  "O",  2.330,  8.750,  2.733},
          { "C2'",  "C",  1.840,  7.040,  2.939},
          { "C1'",  "C",  2.550,  5.860,  2.838},
          { "N9",  "N",  2.500,  4.630,  2.998},
          { "C8",  "C",  2.580,  4.820,  -2.995},
          { "N7",  "N",  2.500,  3.920,  -2.777},
          { "C5",  "C",  2.360,  2.700,  -2.981},
          { "C6",  "C",  2.230,  1.389,  -2.714},
          { "N1",  "N",  2.110,  0.920,  2.323},
          { "C2",  "C",  2.110,  2.281,  2.162},
          { "N2",  "N",  1.980,  3.010,  1.735},
          { "N3",  "N",  2.240,  3.290,  2.477},
          { "C4",  "C",  2.360,  3.330,  2.889},
          { "H5'1",  "H",  5.339,  8.686,  3.052},
          { "H5'2",  "H",  5.611,  6.918,  3.073},
          { "H4'",  "H",  4.743,  8.159,  2.792},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 3, 1, 2 },
          { 14, 15, 1 },
          { 14, 13, 2 },
          { 16, 15, 2 },
          { 16, 17, 1 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 15, 21, 1 },
          { 13, 12, 1 },
          { 17, 18, 2 },
          { 4, 5, 1 },
          { 21, 12, 1 },
          { 21, 20, 2 },
          { 12, 11, 1 },
          { 23, 5, 1 },
          { 5, 22, 1 },
          { 5, 6, 1 },
          { 18, 20, 1 },
          { 18, 19, 1 },
          { 10, 11, 1 },
          { 10, 8, 1 },
          { 11, 7, 1 },
          { 7, 6, 1 },
          { 6, 8, 1 },
          { 6, 24, 1 },
          { 8, 9, 1 },
          { 0, 0, 0}
        }
      },
      { 'M', "ACX", // M ACX Pair Amino (A or C) M -> K
        {
          { "P",  "P",  4.160,  8.910,  -2.960},
          { "O1P",  "O",  3.940,  10.200,  -3.032},
          { "O2P",  "O",  3.370,  8.820,  -2.819},
          { "O5'",  "O",  3.910,  7.730,  -3.086},
          { "C5'",  "C",  4.850,  7.701,  3.054},
          { "C4'",  "C",  4.120,  7.590,  2.881},
          { "O4'",  "O",  3.910,  6.220,  2.813},
          { "C3'",  "C",  2.720,  8.200,  2.882},
          { "O3'",  "O",  2.330,  8.750,  2.733},
          { "C2'",  "C",  1.840,  7.040,  2.939},
          { "C1'",  "C",  2.550,  5.860,  2.838},
          { "N9",  "N",  2.500,  4.630,  2.999},
          { "H5'1",  "H",  5.339,  8.686,  3.052},
          { "H5'2",  "H",  5.611,  6.918,  3.074},
          { "H4'",  "H",  4.743,  8.159,  2.792},
          { "H1'",  "H",  2.088,  5.803,  2.667},
          { "H2'1",  "H",  1.808,  7.072,  3.095},
          { "H2'2",  "H",  0.817,  7.145,  2.884},
          { "H3'",  "H",  2.651,  9.056,  2.961},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 3, 1, 2 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 4, 5, 1 },
          { 12, 11, 1 },
          { 17, 10, 1 },
          { 14, 5, 1 },
          { 5, 13, 1 },
          { 5, 6, 1 },
          { 10, 11, 1 },
          { 10, 18, 1 },
          { 10, 8, 1 },
          { 11, 7, 1 },
          { 11, 16, 1 },
          { 7, 6, 1 },
          { 19, 8, 1 },
          { 6, 8, 1 },
          { 6, 15, 1 },
          { 8, 9, 1 },
          { 0, 0, 0}
        }
      },
      { 'K', "GTX", // K GTX Pair Keto (G or T) K -> M
        {
          { "P",  "P",  4.160,  8.910,  -2.960},
          { "O1P",  "O",  3.940,  10.200,  -3.032},
          { "O2P",  "O",  3.370,  8.820,  -2.819},
          { "O5'",  "O",  3.910,  7.730,  -3.086},
          { "C5'",  "C",  4.850,  7.701,  3.054},
          { "C4'",  "C",  4.120,  7.590,  2.881},
          { "O4'",  "O",  3.910,  6.220,  2.813},
          { "C3'",  "C",  2.720,  8.200,  2.882},
          { "O3'",  "O",  2.330,  8.750,  2.733},
          { "C2'",  "C",  1.840,  7.040,  2.939},
          { "C1'",  "C",  2.550,  5.860,  2.838},
          { "N9",  "N",  2.500,  4.630,  2.998},
          { "H5'1",  "H",  5.339,  8.686,  3.052},
          { "H5'2",  "H",  5.611,  6.918,  3.073},
          { "H4'",  "H",  4.743,  8.159,  2.792},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 3, 1, 2 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 4, 5, 1 },
          { 12, 11, 1 },
          { 14, 5, 1 },
          { 5, 13, 1 },
          { 5, 6, 1 },
          { 10, 11, 1 },
          { 10, 8, 1 },
          { 11, 7, 1 },
          { 7, 6, 1 },
          { 6, 8, 1 },
          { 6, 15, 1 },
          { 8, 9, 1 },
          { 0, 0, 0}
        }
      },
      { 'W', "ATX", // W ATX Pair Weak (A or T) W -> S
        {
          { "P",  "P",  4.160,  8.910,  -2.960},
          { "O1P",  "O",  3.940,  10.200,  -3.032},
          { "O2P",  "O",  3.370,  8.820,  -2.819},
          { "O5'",  "O",  3.910,  7.730,  -3.086},
          { "C5'",  "C",  4.850,  7.701,  3.054},
          { "C4'",  "C",  4.120,  7.590,  2.881},
          { "O4'",  "O",  3.910,  6.220,  2.813},
          { "C3'",  "C",  2.720,  8.200,  2.882},
          { "O3'",  "O",  2.330,  8.750,  2.733},
          { "C2'",  "C",  1.840,  7.040,  2.939},
          { "C1'",  "C",  2.550,  5.860,  2.838},
          { "N9",  "N",  2.500,  4.630,  2.999},
          { "H5'1",  "H",  5.339,  8.686,  3.052},
          { "H5'2",  "H",  5.611,  6.918,  3.074},
          { "H4'",  "H",  4.743,  8.159,  2.792},
          { "H1'",  "H",  2.088,  5.803,  2.667},
          { "H2'1",  "H",  1.808,  7.072,  3.095},
          { "H2'2",  "H",  0.817,  7.145,  2.884},
          { "H3'",  "H",  2.651,  9.056,  2.961},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 3, 1, 2 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 4, 5, 1 },
          { 12, 11, 1 },
          { 17, 10, 1 },
          { 14, 5, 1 },
          { 5, 13, 1 },
          { 5, 6, 1 },
          { 10, 11, 1 },
          { 10, 18, 1 },
          { 10, 8, 1 },
          { 11, 7, 1 },
          { 11, 16, 1 },
          { 7, 6, 1 },
          { 19, 8, 1 },
          { 6, 8, 1 },
          { 6, 15, 1 },
          { 8, 9, 1 },
          { 0, 0, 0}
        }
      },
      { 'S', "GCX", // S GCX Pair Strong (G or C) S -> W
        {
          { "P",  "P",  4.160,  8.910,  -2.960},
          { "O1P",  "O",  3.940,  10.200,  -3.032},
          { "O2P",  "O",  3.370,  8.820,  -2.819},
          { "O5'",  "O",  3.910,  7.730,  -3.086},
          { "C5'",  "C",  4.850,  7.701,  3.054},
          { "C4'",  "C",  4.120,  7.590,  2.881},
          { "O4'",  "O",  3.910,  6.220,  2.813},
          { "C3'",  "C",  2.720,  8.200,  2.882},
          { "O3'",  "O",  2.330,  8.750,  2.733},
          { "C2'",  "C",  1.840,  7.040,  2.939},
          { "C1'",  "C",  2.550,  5.860,  2.838},
          { "N9",  "N",  2.500,  4.630,  2.998},
          { "H5'1",  "H",  5.339,  8.686,  3.052},
          { "H5'2",  "H",  5.611,  6.918,  3.073},
          { "H4'",  "H",  4.743,  8.159,  2.792},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 3, 1, 2 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 4, 5, 1 },
          { 12, 11, 1 },
          { 14, 5, 1 },
          { 5, 13, 1 },
          { 5, 6, 1 },
          { 10, 11, 1 },
          { 10, 8, 1 },
          { 11, 7, 1 },
          { 7, 6, 1 },
          { 6, 8, 1 },
          { 6, 15, 1 },
          { 8, 9, 1 },
          { 0, 0, 0}
        }
      },
      { 'V', "GCA", // V GCA Pair (G C A) V -> B
        {
          { "P",  "P",  4.160,  8.910,  -2.960},
          { "O1P",  "O",  3.940,  10.200,  -3.032},
          { "O2P",  "O",  3.370,  8.820,  -2.819},
          { "O5'",  "O",  3.910,  7.730,  -3.086},
          { "C5'",  "C",  4.850,  7.701,  3.054},
          { "C4'",  "C",  4.120,  7.590,  2.881},
          { "O4'",  "O",  3.910,  6.220,  2.813},
          { "C3'",  "C",  2.720,  8.200,  2.882},
          { "O3'",  "O",  2.330,  8.750,  2.733},
          { "C2'",  "C",  1.840,  7.040,  2.939},
          { "C1'",  "C",  2.550,  5.860,  2.838},
          { "N9",  "N",  2.500,  4.630,  2.998},
          { "H5'1",  "H",  5.339,  8.686,  3.052},
          { "H5'2",  "H",  5.611,  6.918,  3.073},
          { "H4'",  "H",  4.743,  8.159,  2.792},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 3, 1, 2 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 4, 5, 1 },
          { 12, 11, 1 },
          { 14, 5, 1 },
          { 5, 13, 1 },
          { 5, 6, 1 },
          { 10, 11, 1 },
          { 10, 8, 1 },
          { 11, 7, 1 },
          { 7, 6, 1 },
          { 6, 8, 1 },
          { 6, 15, 1 },
          { 8, 9, 1 },
          { 0, 0, 0}
        }
      },
      { 'H', "ACT", // H ACT Pair (A C T) H -> D
        {
          { "P",  "P",  4.160,  8.910,  -2.960},
          { "O1P",  "O",  3.940,  10.200,  -3.032},
          { "O2P",  "O",  3.370,  8.820,  -2.819},
          { "O5'",  "O",  3.910,  7.730,  -3.086},
          { "C5'",  "C",  4.850,  7.701,  3.054},
          { "C4'",  "C",  4.120,  7.590,  2.881},
          { "O4'",  "O",  3.910,  6.220,  2.813},
          { "C3'",  "C",  2.720,  8.200,  2.882},
          { "O3'",  "O",  2.330,  8.750,  2.733},
          { "C2'",  "C",  1.840,  7.040,  2.939},
          { "C1'",  "C",  2.550,  5.860,  2.838},
          { "N9",  "N",  2.500,  4.630,  2.999},
          { "H5'1",  "H",  5.339,  8.686,  3.052},
          { "H5'2",  "H",  5.611,  6.918,  3.074},
          { "H4'",  "H",  4.743,  8.159,  2.792},
          { "H1'",  "H",  2.088,  5.803,  2.667},
          { "H2'1",  "H",  1.808,  7.072,  3.095},
          { "H2'2",  "H",  0.817,  7.145,  2.884},
          { "H3'",  "H",  2.651,  9.056,  2.961},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 3, 1, 2 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 4, 5, 1 },
          { 12, 11, 1 },
          { 17, 10, 1 },
          { 14, 5, 1 },
          { 5, 13, 1 },
          { 5, 6, 1 },
          { 10, 11, 1 },
          { 10, 18, 1 },
          { 10, 8, 1 },
          { 11, 7, 1 },
          { 11, 16, 1 },
          { 7, 6, 1 },
          { 19, 8, 1 },
          { 6, 8, 1 },
          { 6, 15, 1 },
          { 8, 9, 1 },
          { 0, 0, 0}
        }
      },
      { 'D', "GAT", // D GAT Pair (G A T) D -> H
        {
          { "P",  "P",  4.160,  8.910,  -2.960},
          { "O1P",  "O",  3.940,  10.200,  -3.032},
          { "O2P",  "O",  3.370,  8.820,  -2.819},
          { "O5'",  "O",  3.910,  7.730,  -3.086},
          { "C5'",  "C",  4.850,  7.701,  3.054},
          { "C4'",  "C",  4.120,  7.590,  2.881},
          { "O4'",  "O",  3.910,  6.220,  2.813},
          { "C3'",  "C",  2.720,  8.200,  2.882},
          { "O3'",  "O",  2.330,  8.750,  2.733},
          { "C2'",  "C",  1.840,  7.040,  2.939},
          { "C1'",  "C",  2.550,  5.860,  2.838},
          { "N9",  "N",  2.500,  4.630,  2.998},
          { "H5'1",  "H",  5.339,  8.686,  3.052},
          { "H5'2",  "H",  5.611,  6.918,  3.073},
          { "H4'",  "H",  4.743,  8.159,  2.792},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 3, 1, 2 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 4, 5, 1 },
          { 12, 11, 1 },
          { 14, 5, 1 },
          { 5, 13, 1 },
          { 5, 6, 1 },
          { 10, 11, 1 },
          { 10, 8, 1 },
          { 11, 7, 1 },
          { 7, 6, 1 },
          { 6, 8, 1 },
          { 6, 15, 1 },
          { 8, 9, 1 },
          { 0, 0, 0}
        }
      },
      { 'B', "GTC", // B GTC Pair (G T C) B -> V
        {
          { "P",  "P",  4.160,  8.910,  -2.960},
          { "O1P",  "O",  3.940,  10.200,  -3.032},
          { "O2P",  "O",  3.370,  8.820,  -2.819},
          { "O5'",  "O",  3.910,  7.730,  -3.086},
          { "C5'",  "C",  4.850,  7.701,  3.054},
          { "C4'",  "C",  4.120,  7.590,  2.881},
          { "O4'",  "O",  3.910,  6.220,  2.813},
          { "C3'",  "C",  2.720,  8.200,  2.882},
          { "O3'",  "O",  2.330,  8.750,  2.733},
          { "C2'",  "C",  1.840,  7.040,  2.939},
          { "C1'",  "C",  2.550,  5.860,  2.838},
          { "N9",  "N",  2.500,  4.630,  2.998},
          { "H5'1",  "H",  5.339,  8.686,  3.052},
          { "H5'2",  "H",  5.611,  6.918,  3.073},
          { "H4'",  "H",  4.743,  8.159,  2.792},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 3, 1, 2 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 4, 5, 1 },
          { 12, 11, 1 },
          { 14, 5, 1 },
          { 5, 13, 1 },
          { 5, 6, 1 },
          { 10, 11, 1 },
          { 10, 8, 1 },
          { 11, 7, 1 },
          { 7, 6, 1 },
          { 6, 8, 1 },
          { 6, 15, 1 },
          { 8, 9, 1 },
          { 0, 0, 0}
        }
      }
    };

  ResidueRecord RNAResidues[IUPAC_RNA_max] =
    {
      { 0, "", // RNA Start
        {
          { "HTER",  "H",  -2.588,  8.491,  -0.008},
          { "OXT",  "O",  -1.492,  9.373,  0.009},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 2, 1, 1},
          { 0, 0, 0}
        }
      },
      { 0, "", // RNA End
        {
          { "HCAP",  "H",  3.646,  8.242,  0.567},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 0, 0, 0}
        }
      },
      { 'N', "XNA", // N XNA Unknown Nucleic Acid N -> N [As K Below]
        {
          { "P",  "P",  -0.000,  8.791,  -0.000},
          { "O1P",  "O",  1.422,  9.311,  0.003},
          { "O2P",  "O",  -0.360,  7.848,  -0.133},
          { "O5'",  "O",  -0.323,  8.218,  0.173},
          { "C5'",  "C",  -0.483,  9.385,  0.274},
          { "C4'",  "C",  0.057,  9.172,  0.426},
          { "O4'",  "O",  -0.999,  8.642,  0.501},
          { "C3'",  "C",  1.275,  8.307,  0.458},
          { "O3'",  "O",  2.013,  8.747,  0.588},
          { "C2'",  "C",  0.606,  6.917,  0.468},
          { "O2'",  "O",  1.405,  6.122,  0.603},
          { "C1'",  "C",  -0.731,  7.281,  0.552},
          { "H1'",  "H",  -0.699,  7.317,  0.703},
          { "H2'",  "H",  0.484,  6.539,  0.315},
          { "H3'",  "H",  1.968,  8.387,  0.356},
          { "H4'",  "H",  0.274,  10.189,  0.463},
          { "H5'1",  "H",  0.084,  10.229,  0.231},
          { "H5'2",  "H",  -1.546,  9.658,  0.281},
          { "HO2'",  "H",  1.093,  6.345,  0.742},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 3, 1, 2 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 4, 5, 1 },
          { 14, 10, 1 },
          { 17, 5, 1 },
          { 5, 18, 1 },
          { 5, 6, 1 },
          { 15, 8, 1 },
          { 10, 11, 1 },
          { 10, 8, 1 },
          { 10, 12, 1 },
          { 11, 19, 1 },
          { 8, 6, 1 },
          { 8, 9, 1 },
          { 6, 7, 1 },
          { 6, 16, 1 },
          { 12, 7, 1 },
          { 12, 13, 1 },
          { 0, 0, 0}
        }
      },
      { 'A', "A", // A A   Adenine A -> U
        {
          { "P",  "P",  -0.000,  8.791,  -0.000},
          { "O1P",  "O",  1.421,  9.312,  0.003},
          { "O2P",  "O",  -0.360,  7.848,  -0.133},
          { "O5'",  "O",  -0.323,  8.218,  0.173},
          { "C5'",  "C",  -0.484,  9.386,  0.274},
          { "C4'",  "C",  0.057,  9.172,  0.426},
          { "O4'",  "O",  -0.998,  8.641,  0.501},
          { "C3'",  "C",  1.275,  8.307,  0.458},
          { "O3'",  "O",  2.013,  8.747,  0.588},
          { "C2'",  "C",  0.607,  6.917,  0.468},
          { "O2'",  "O",  1.405,  6.121,  0.602},
          { "C1'",  "C",  -0.731,  7.281,  0.552},
          { "N9",  "N",  -1.774,  6.424,  0.474},
          { "C8",  "C",  -1.605,  5.989,  0.269},
          { "N7",  "N",  -2.632,  5.390,  0.170},
          { "C5",  "C",  -3.517,  5.270,  0.368},
          { "C6",  "C",  -4.784,  4.675,  0.397},
          { "N6",  "N",  -5.383,  4.048,  0.157},
          { "N1",  "N",  -5.349,  5.072,  0.640},
          { "C2",  "C",  -4.747,  6.066,  0.771},
          { "N3",  "N",  -3.570,  6.563,  0.726},
          { "C4",  "C",  -2.995,  6.046,  0.547},
          { "H1'",  "H",  -0.699,  7.316,  0.703},
          { "H2'",  "H",  0.484,  6.538,  0.316},
          { "H3'",  "H",  1.967,  8.387,  0.356},
          { "H4'",  "H",  0.273,  10.189,  0.464},
          { "H5'1",  "H",  0.083,  10.229,  0.231},
          { "H5'2",  "H",  -1.547,  9.658,  0.281},
          { "H2",  "H",  -5.283,  6.564,  0.898},
          { "H8",  "H",  -0.705,  6.267,  0.180},
          { "H61",  "H",  -6.313,  3.610,  0.176},
          { "H62",  "H",  -4.909,  4.211,  -0.061},
          { "HO2'",  "H",  1.093,  6.344,  0.742},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 3, 1, 2 },
          { 32, 18, 1 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 31, 18, 1 },
          { 18, 17, 1 },
          { 15, 14, 2 },
          { 15, 16, 1 },
          { 30, 14, 1 },
          { 4, 5, 1 },
          { 14, 13, 1 },
          { 17, 16, 2 },
          { 17, 19, 1 },
          { 16, 22, 1 },
          { 24, 10, 1 },
          { 27, 5, 1 },
          { 5, 28, 1 },
          { 5, 6, 1 },
          { 25, 8, 1 },
          { 13, 22, 1 },
          { 13, 12, 1 },
          { 19, 20, 2 },
          { 10, 11, 1 },
          { 10, 8, 1 },
          { 10, 12, 1 },
          { 22, 21, 2 },
          { 11, 33, 1 },
          { 8, 6, 1 },
          { 8, 9, 1 },
          { 6, 7, 1 },
          { 6, 26, 1 },
          { 12, 7, 1 },
          { 12, 23, 1 },
          { 20, 21, 1 },
          { 20, 29, 1 },
          { 0, 0, 0}
        }
      },
      { 'C', "C", // C C   Cytosine C -> G
        {
          { "P",  "P",  -0.000,  8.791,  -0.000},
          { "O1P",  "O",  1.401,  9.367,  -0.002},
          { "O2P",  "O",  -0.378,  7.899,  -0.137},
          { "O5'",  "O",  -0.323,  8.217,  0.173},
          { "C5'",  "C",  -0.483,  9.385,  0.274},
          { "C4'",  "C",  0.057,  9.171,  0.426},
          { "O4'",  "O",  -0.999,  8.641,  0.501},
          { "C3'",  "C",  1.274,  8.307,  0.458},
          { "O3'",  "O",  2.013,  8.747,  0.588},
          { "C2'",  "C",  0.607,  6.917,  0.468},
          { "O2'",  "O",  1.405,  6.121,  0.603},
          { "C1'",  "C",  -0.731,  7.281,  0.552},
          { "N1",  "N",  -1.774,  6.424,  0.474},
          { "C2",  "C",  -2.809,  6.085,  0.608},
          { "O2",  "O",  -2.770,  6.659,  0.783},
          { "N3",  "N",  -3.821,  5.311,  0.529},
          { "C4",  "C",  -3.773,  5.045,  0.275},
          { "N4",  "N",  -4.783,  4.492,  0.136},
          { "C5",  "C",  -2.677,  5.611,  0.141},
          { "C6",  "C",  -1.711,  6.187,  0.263},
          { "H1'",  "H",  -0.699,  7.316,  0.703},
          { "H2'",  "H",  0.484,  6.538,  0.316},
          { "H3'",  "H",  1.968,  8.388,  0.357},
          { "H4'",  "H",  0.274,  10.189,  0.464},
          { "H5'1",  "H",  0.083,  10.229,  0.231},
          { "H5'2",  "H",  -1.547,  9.657,  0.281},
          { "H5",  "H",  -2.629,  5.749,  -0.049},
          { "H6",  "H",  -0.865,  6.641,  0.182},
          { "H41",  "H",  -5.573,  4.074,  0.256},
          { "H42",  "H",  -4.778,  4.661,  -0.086},
          { "HO2'",  "H",  1.189,  6.451,  0.740},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 3, 1, 2 },
          { 30, 18, 1 },
          { 27, 19, 1 },
          { 2, 1, 1 },
          { 1, 4, 1 },
          { 18, 29, 1 },
          { 18, 17, 1 },
          { 19, 17, 1 },
          { 19, 20, 2 },
          { 28, 20, 1 },
          { 17, 16, 2 },
          { 4, 5, 1 },
          { 20, 13, 1 },
          { 22, 10, 1 },
          { 25, 5, 1 },
          { 5, 26, 1 },
          { 5, 6, 1 },
          { 16, 14, 1 },
          { 23, 8, 1 },
          { 13, 14, 1 },
          { 13, 12, 1 },
          { 10, 11, 1 },
          { 10, 8, 1 },
          { 10, 12, 1 },
          { 11, 31, 1 },
          { 14, 15, 2 },
          { 8, 6, 1 },
          { 8, 9, 1 },
          { 6, 7, 1 },
          { 6, 24, 1 },
          { 12, 7, 1 },
          { 12, 21, 1 },
          { 0, 0, 0}
        }
      },
      { 'G', "G", // G G   Guanine G -> C
        {
          { "P",  "P",  -0.000,  8.791,  -0.000},
          { "O1P",  "O",  1.422,  9.313,  0.003},
          { "O2P",  "O",  -0.360,  7.848,  -0.133},
          { "O5'",  "O",  -0.323,  8.218,  0.173},
          { "C5'",  "C",  -0.483,  9.385,  0.274},
          { "C4'",  "C",  0.057,  9.171,  0.426},
          { "O4'",  "O",  -0.998,  8.642,  0.501},
          { "C3'",  "C",  1.275,  8.307,  0.458},
          { "O3'",  "O",  2.013,  8.747,  0.588},
          { "C2'",  "C",  0.607,  6.917,  0.468},
          { "O2'",  "O",  1.405,  6.121,  0.603},
          { "C1'",  "C",  -0.731,  7.281,  0.552},
          { "N9",  "N",  -1.774,  6.424,  0.474},
          { "C8",  "C",  -1.712,  6.189,  0.266},
          { "N7",  "N",  -2.757,  5.529,  0.179},
          { "C5",  "C",  -3.491,  5.243,  0.379},
          { "C6",  "C",  -4.797,  4.658,  0.414},
          { "O6",  "O",  -5.530,  4.203,  0.212},
          { "N1",  "N",  -5.227,  4.911,  0.678},
          { "C2",  "C",  -4.478,  5.823,  0.825},
          { "N2",  "N",  -5.075,  6.378,  1.007},
          { "N3",  "N",  -3.325,  6.311,  0.765},
          { "C4",  "C",  -2.868,  5.934,  0.566},
          { "H1'",  "H",  -0.699,  7.317,  0.703},
          { "H2'",  "H",  0.484,  6.538,  0.315},
          { "H3'",  "H",  1.969,  8.387,  0.357},
          { "H4'",  "H",  0.274,  10.189,  0.463},
          { "H5'1",  "H",  0.084,  10.229,  0.231},
          { "H5'2",  "H",  -1.547,  9.658,  0.281},
          { "H1",  "H",  -6.123,  4.611,  0.747},
          { "H8",  "H",  -0.916,  6.638,  0.173},
          { "H21",  "H",  -5.986,  6.065,  1.065},
          { "H22",  "H",  -4.610,  7.201,  1.067},
          { "HO2'",  "H",  1.094,  6.344,  0.742},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 3, 1, 2 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 18, 17, 2 },
          { 15, 14, 2 },
          { 15, 16, 1 },
          { 31, 14, 1 },
          { 4, 5, 1 },
          { 14, 13, 1 },
          { 17, 16, 1 },
          { 17, 19, 1 },
          { 16, 23, 2 },
          { 25, 10, 1 },
          { 28, 5, 1 },
          { 5, 29, 1 },
          { 5, 6, 1 },
          { 26, 8, 1 },
          { 13, 23, 1 },
          { 13, 12, 1 },
          { 19, 30, 1 },
          { 19, 20, 1 },
          { 10, 11, 1 },
          { 10, 8, 1 },
          { 10, 12, 1 },
          { 23, 22, 1 },
          { 11, 34, 1 },
          { 8, 6, 1 },
          { 8, 9, 1 },
          { 6, 7, 1 },
          { 6, 27, 1 },
          { 12, 7, 1 },
          { 12, 24, 1 },
          { 20, 22, 2 },
          { 20, 21, 1 },
          { 32, 21, 1 },
          { 21, 33, 1 },
          { 0, 0, 0}
        }
      },
      { 'U', "U", // U U   Uridine U -> A
        {
          { "P",  "P",  -0.000,  8.791,  -0.000},
          { "O1P",  "O",  1.422,  9.311,  0.003},
          { "O2P",  "O",  -0.360,  7.848,  -0.133},
          { "O5'",  "O",  -0.323,  8.218,  0.173},
          { "C5'",  "C",  -0.483,  9.385,  0.274},
          { "C4'",  "C",  0.057,  9.172,  0.426},
          { "O4'",  "O",  -0.999,  8.642,  0.501},
          { "C3'",  "C",  1.275,  8.307,  0.458},
          { "O3'",  "O",  2.013,  8.747,  0.588},
          { "C2'",  "C",  0.606,  6.917,  0.468},
          { "O2'",  "O",  1.405,  6.122,  0.603},
          { "C1'",  "C",  -0.731,  7.281,  0.552},
          { "N1",  "N",  -1.774,  6.424,  0.474},
          { "C2",  "C",  -2.902,  6.238,  0.596},
          { "O2",  "O",  -3.049,  6.885,  0.754},
          { "N3",  "N",  -3.840,  5.395,  0.501},
          { "C4",  "C",  -3.761,  4.935,  0.245},
          { "O4",  "O",  -4.675,  4.400,  0.123},
          { "C5",  "C",  -2.555,  5.445,  0.125},
          { "C6",  "C",  -1.622,  6.024,  0.264},
          { "H1'",  "H",  -0.699,  7.317,  0.703},
          { "H2'",  "H",  0.484,  6.539,  0.315},
          { "H3'",  "H",  1.968,  8.387,  0.356},
          { "H4'",  "H",  0.274,  10.189,  0.463},
          { "H5'1",  "H",  0.084,  10.229,  0.231},
          { "H5'2",  "H",  -1.546,  9.658,  0.281},
          { "H3",  "H",  -4.674,  5.253,  0.601},
          { "H5",  "H",  -2.401,  5.543,  -0.071},
          { "H6",  "H",  -0.700,  6.362,  0.188},
          { "HO2'",  "H",  1.093,  6.345,  0.742},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 3, 1, 2 },
          { 28, 19, 1 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 18, 17, 2 },
          { 19, 17, 1 },
          { 19, 20, 2 },
          { 29, 20, 1 },
          { 17, 16, 1 },
          { 4, 5, 1 },
          { 20, 13, 1 },
          { 22, 10, 1 },
          { 25, 5, 1 },
          { 5, 26, 1 },
          { 5, 6, 1 },
          { 16, 27, 1 },
          { 16, 14, 1 },
          { 23, 8, 1 },
          { 13, 14, 1 },
          { 13, 12, 1 },
          { 10, 11, 1 },
          { 10, 8, 1 },
          { 10, 12, 1 },
          { 11, 30, 1 },
          { 14, 15, 2 },
          { 8, 6, 1 },
          { 8, 9, 1 },
          { 6, 7, 1 },
          { 6, 24, 1 },
          { 12, 7, 1 },
          { 12, 21, 1 },
          { 0, 0, 0}
        }
      },
      { 'R', "GAX", // R GAX Purine (G or A) R -> Y [A - N6,H61,H62,H2]
        {
          { "P",  "P",  -0.000,  8.791,  -0.000},
          { "O1P",  "O",  1.421,  9.312,  0.003},
          { "O2P",  "O",  -0.360,  7.848,  -0.133},
          { "O5'",  "O",  -0.323,  8.218,  0.173},
          { "C5'",  "C",  -0.484,  9.386,  0.274},
          { "C4'",  "C",  0.057,  9.172,  0.426},
          { "O4'",  "O",  -0.998,  8.641,  0.501},
          { "C3'",  "C",  1.275,  8.307,  0.458},
          { "O3'",  "O",  2.013,  8.747,  0.588},
          { "C2'",  "C",  0.607,  6.917,  0.468},
          { "O2'",  "O",  1.405,  6.121,  0.602},
          { "C1'",  "C",  -0.731,  7.281,  0.552},
          { "N9",  "N",  -1.774,  6.424,  0.474},
          { "C8",  "C",  -1.605,  5.989,  0.269},
          { "N7",  "N",  -2.632,  5.390,  0.170},
          { "C5",  "C",  -3.517,  5.270,  0.368},
          { "C6",  "C",  -4.784,  4.675,  0.397},
          { "N1",  "N",  -5.349,  5.072,  0.640},
          { "C2",  "C",  -4.747,  6.066,  0.771},
          { "N3",  "N",  -3.570,  6.563,  0.726},
          { "C4",  "C",  -2.995,  6.046,  0.547},
          { "H2'",  "H",  0.484,  6.538,  0.316},
          { "H3'",  "H",  1.967,  8.387,  0.356},
          { "H4'",  "H",  0.273,  10.189,  0.464},
          { "H5'1",  "H",  0.083,  10.229,  0.231},
          { "H5'2",  "H",  -1.547,  9.658,  0.281},
          { "H8",  "H",  -0.705,  6.267,  0.180},
          { "HO2'",  "H",  1.093,  6.344,  0.742},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 3, 1, 2 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 15, 14, 2 },
          { 15, 16, 1 },
          { 27, 14, 1 },
          { 4, 5, 1 },
          { 14, 13, 1 },
          { 17, 16, 2 },
          { 17, 18, 1 },
          { 16, 21, 1 },
          { 22, 10, 1 },
          { 25, 5, 1 },
          { 5, 26, 1 },
          { 5, 6, 1 },
          { 23, 8, 1 },
          { 13, 21, 1 },
          { 13, 12, 1 },
          { 18, 19, 2 },
          { 10, 11, 1 },
          { 10, 8, 1 },
          { 10, 12, 1 },
          { 21, 20, 2 },
          { 11, 28, 1 },
          { 8, 6, 1 },
          { 8, 9, 1 },
          { 6, 7, 1 },
          { 6, 24, 1 },
          { 12, 7, 1 },
          { 19, 20, 1 },
          { 0, 0, 0}
        }
      },
      { 'Y', "UCX", // Y UCX Pyrimidine (U or C) Y -> R [C - N4,H41,H42]
        {
          { "P",  "P",  -0.000,  8.791,  -0.000},
          { "O1P",  "O",  1.401,  9.367,  -0.002},
          { "O2P",  "O",  -0.378,  7.899,  -0.137},
          { "O5'",  "O",  -0.323,  8.217,  0.173},
          { "C5'",  "C",  -0.483,  9.385,  0.274},
          { "C4'",  "C",  0.057,  9.171,  0.426},
          { "O4'",  "O",  -0.999,  8.641,  0.501},
          { "C3'",  "C",  1.274,  8.307,  0.458},
          { "O3'",  "O",  2.013,  8.747,  0.588},
          { "C2'",  "C",  0.607,  6.917,  0.468},
          { "O2'",  "O",  1.405,  6.121,  0.603},
          { "C1'",  "C",  -0.731,  7.281,  0.552},
          { "N1",  "N",  -1.774,  6.424,  0.474},
          { "C2",  "C",  -2.809,  6.085,  0.608},
          { "O2",  "O",  -2.770,  6.659,  0.783},
          { "N3",  "N",  -3.821,  5.311,  0.529},
          { "C4",  "C",  -3.773,  5.045,  0.275},
          { "C5",  "C",  -2.677,  5.611,  0.141},
          { "C6",  "C",  -1.711,  6.187,  0.263},
          { "H1'",  "H",  -0.699,  7.316,  0.703},
          { "H2'",  "H",  0.484,  6.538,  0.316},
          { "H3'",  "H",  1.968,  8.388,  0.357},
          { "H4'",  "H",  0.274,  10.189,  0.464},
          { "H5'1",  "H",  0.083,  10.229,  0.231},
          { "H5'2",  "H",  -1.547,  9.657,  0.281},
          { "H5",  "H",  -2.629,  5.749,  -0.049},
          { "H6",  "H",  -0.865,  6.641,  0.182},
          { "HO2'",  "H",  1.189,  6.451,  0.740},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 3, 1, 2 },
          { 26, 18, 1 },
          { 2, 1, 1 },
          { 1, 4, 1 },
          { 18, 17, 1 },
          { 18, 19, 2 },
          { 27, 19, 1 },
          { 17, 16, 2 },
          { 4, 5, 1 },
          { 19, 13, 1 },
          { 21, 10, 1 },
          { 24, 5, 1 },
          { 5, 25, 1 },
          { 5, 6, 1 },
          { 16, 14, 1 },
          { 22, 8, 1 },
          { 13, 14, 1 },
          { 13, 12, 1 },
          { 10, 11, 1 },
          { 10, 8, 1 },
          { 10, 12, 1 },
          { 11, 28, 1 },
          { 14, 15, 2 },
          { 8, 6, 1 },
          { 8, 9, 1 },
          { 6, 7, 1 },
          { 6, 23, 1 },
          { 12, 7, 1 },
          { 12, 20, 1 },
          { 0, 0, 0}
        }
      },
      { 'K', "GUX", // K GUX Keto (G or U) K -> M [U - N1,C6,H6,C5,H5,C4,O4,N3,H3,C2,O2]
        {
          { "P",  "P",  -0.000,  8.791,  -0.000},
          { "O1P",  "O",  1.422,  9.311,  0.003},
          { "O2P",  "O",  -0.360,  7.848,  -0.133},
          { "O5'",  "O",  -0.323,  8.218,  0.173},
          { "C5'",  "C",  -0.483,  9.385,  0.274},
          { "C4'",  "C",  0.057,  9.172,  0.426},
          { "O4'",  "O",  -0.999,  8.642,  0.501},
          { "C3'",  "C",  1.275,  8.307,  0.458},
          { "O3'",  "O",  2.013,  8.747,  0.588},
          { "C2'",  "C",  0.606,  6.917,  0.468},
          { "O2'",  "O",  1.405,  6.122,  0.603},
          { "C1'",  "C",  -0.731,  7.281,  0.552},
          { "H1'",  "H",  -0.699,  7.317,  0.703},
          { "H2'",  "H",  0.484,  6.539,  0.315},
          { "H3'",  "H",  1.968,  8.387,  0.356},
          { "H4'",  "H",  0.274,  10.189,  0.463},
          { "H5'1",  "H",  0.084,  10.229,  0.231},
          { "H5'2",  "H",  -1.546,  9.658,  0.281},
          { "HO2'",  "H",  1.093,  6.345,  0.742},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 3, 1, 2 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 4, 5, 1 },
          { 14, 10, 1 },
          { 17, 5, 1 },
          { 5, 18, 1 },
          { 5, 6, 1 },
          { 15, 8, 1 },
          { 10, 11, 1 },
          { 10, 8, 1 },
          { 10, 12, 1 },
          { 11, 19, 1 },
          { 8, 6, 1 },
          { 8, 9, 1 },
          { 6, 7, 1 },
          { 6, 16, 1 },
          { 12, 7, 1 },
          { 12, 13, 1 },
          { 0, 0, 0}
        }
      },
      { 'M', "ACX", // M ACX Amino (A or C) M -> K [As K Above]
        {
          { "P",  "P",  -0.000,  8.791,  -0.000},
          { "O1P",  "O",  1.422,  9.311,  0.003},
          { "O2P",  "O",  -0.360,  7.848,  -0.133},
          { "O5'",  "O",  -0.323,  8.218,  0.173},
          { "C5'",  "C",  -0.483,  9.385,  0.274},
          { "C4'",  "C",  0.057,  9.172,  0.426},
          { "O4'",  "O",  -0.999,  8.642,  0.501},
          { "C3'",  "C",  1.275,  8.307,  0.458},
          { "O3'",  "O",  2.013,  8.747,  0.588},
          { "C2'",  "C",  0.606,  6.917,  0.468},
          { "O2'",  "O",  1.405,  6.122,  0.603},
          { "C1'",  "C",  -0.731,  7.281,  0.552},
          { "H1'",  "H",  -0.699,  7.317,  0.703},
          { "H2'",  "H",  0.484,  6.539,  0.315},
          { "H3'",  "H",  1.968,  8.387,  0.356},
          { "H4'",  "H",  0.274,  10.189,  0.463},
          { "H5'1",  "H",  0.084,  10.229,  0.231},
          { "H5'2",  "H",  -1.546,  9.658,  0.281},
          { "HO2'",  "H",  1.093,  6.345,  0.742},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 3, 1, 2 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 4, 5, 1 },
          { 14, 10, 1 },
          { 17, 5, 1 },
          { 5, 18, 1 },
          { 5, 6, 1 },
          { 15, 8, 1 },
          { 10, 11, 1 },
          { 10, 8, 1 },
          { 10, 12, 1 },
          { 11, 19, 1 },
          { 8, 6, 1 },
          { 8, 9, 1 },
          { 6, 7, 1 },
          { 6, 16, 1 },
          { 12, 7, 1 },
          { 12, 13, 1 },
          { 0, 0, 0}
        }
      },
      { 'S', "GCX", // S GCX Strong (G or C) S -> W [As K Above]
        {
          { "P",  "P",  -0.000,  8.791,  -0.000},
          { "O1P",  "O",  1.422,  9.311,  0.003},
          { "O2P",  "O",  -0.360,  7.848,  -0.133},
          { "O5'",  "O",  -0.323,  8.218,  0.173},
          { "C5'",  "C",  -0.483,  9.385,  0.274},
          { "C4'",  "C",  0.057,  9.172,  0.426},
          { "O4'",  "O",  -0.999,  8.642,  0.501},
          { "C3'",  "C",  1.275,  8.307,  0.458},
          { "O3'",  "O",  2.013,  8.747,  0.588},
          { "C2'",  "C",  0.606,  6.917,  0.468},
          { "O2'",  "O",  1.405,  6.122,  0.603},
          { "C1'",  "C",  -0.731,  7.281,  0.552},
          { "H1'",  "H",  -0.699,  7.317,  0.703},
          { "H2'",  "H",  0.484,  6.539,  0.315},
          { "H3'",  "H",  1.968,  8.387,  0.356},
          { "H4'",  "H",  0.274,  10.189,  0.463},
          { "H5'1",  "H",  0.084,  10.229,  0.231},
          { "H5'2",  "H",  -1.546,  9.658,  0.281},
          { "HO2'",  "H",  1.093,  6.345,  0.742},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 3, 1, 2 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 4, 5, 1 },
          { 14, 10, 1 },
          { 17, 5, 1 },
          { 5, 18, 1 },
          { 5, 6, 1 },
          { 15, 8, 1 },
          { 10, 11, 1 },
          { 10, 8, 1 },
          { 10, 12, 1 },
          { 11, 19, 1 },
          { 8, 6, 1 },
          { 8, 9, 1 },
          { 6, 7, 1 },
          { 6, 16, 1 },
          { 12, 7, 1 },
          { 12, 13, 1 },
          { 0, 0, 0}
        }
      },
      { 'W', "AUX", // W AUX Weak (A or U) W -> S [As K Above]
        {
          { "P",  "P",  -0.000,  8.791,  -0.000},
          { "O1P",  "O",  1.422,  9.311,  0.003},
          { "O2P",  "O",  -0.360,  7.848,  -0.133},
          { "O5'",  "O",  -0.323,  8.218,  0.173},
          { "C5'",  "C",  -0.483,  9.385,  0.274},
          { "C4'",  "C",  0.057,  9.172,  0.426},
          { "O4'",  "O",  -0.999,  8.642,  0.501},
          { "C3'",  "C",  1.275,  8.307,  0.458},
          { "O3'",  "O",  2.013,  8.747,  0.588},
          { "C2'",  "C",  0.606,  6.917,  0.468},
          { "O2'",  "O",  1.405,  6.122,  0.603},
          { "C1'",  "C",  -0.731,  7.281,  0.552},
          { "H1'",  "H",  -0.699,  7.317,  0.703},
          { "H2'",  "H",  0.484,  6.539,  0.315},
          { "H3'",  "H",  1.968,  8.387,  0.356},
          { "H4'",  "H",  0.274,  10.189,  0.463},
          { "H5'1",  "H",  0.084,  10.229,  0.231},
          { "H5'2",  "H",  -1.546,  9.658,  0.281},
          { "HO2'",  "H",  1.093,  6.345,  0.742},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 3, 1, 2 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 4, 5, 1 },
          { 14, 10, 1 },
          { 17, 5, 1 },
          { 5, 18, 1 },
          { 5, 6, 1 },
          { 15, 8, 1 },
          { 10, 11, 1 },
          { 10, 8, 1 },
          { 10, 12, 1 },
          { 11, 19, 1 },
          { 8, 6, 1 },
          { 8, 9, 1 },
          { 6, 7, 1 },
          { 6, 16, 1 },
          { 12, 7, 1 },
          { 12, 13, 1 },
          { 0, 0, 0}
        }
      },
      { 'B', "GUC", // B GUC (G U C) B -> V [As K Above]
        {
          { "P",  "P",  -0.000,  8.791,  -0.000},
          { "O1P",  "O",  1.422,  9.311,  0.003},
          { "O2P",  "O",  -0.360,  7.848,  -0.133},
          { "O5'",  "O",  -0.323,  8.218,  0.173},
          { "C5'",  "C",  -0.483,  9.385,  0.274},
          { "C4'",  "C",  0.057,  9.172,  0.426},
          { "O4'",  "O",  -0.999,  8.642,  0.501},
          { "C3'",  "C",  1.275,  8.307,  0.458},
          { "O3'",  "O",  2.013,  8.747,  0.588},
          { "C2'",  "C",  0.606,  6.917,  0.468},
          { "O2'",  "O",  1.405,  6.122,  0.603},
          { "C1'",  "C",  -0.731,  7.281,  0.552},
          { "H1'",  "H",  -0.699,  7.317,  0.703},
          { "H2'",  "H",  0.484,  6.539,  0.315},
          { "H3'",  "H",  1.968,  8.387,  0.356},
          { "H4'",  "H",  0.274,  10.189,  0.463},
          { "H5'1",  "H",  0.084,  10.229,  0.231},
          { "H5'2",  "H",  -1.546,  9.658,  0.281},
          { "HO2'",  "H",  1.093,  6.345,  0.742},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 3, 1, 2 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 4, 5, 1 },
          { 14, 10, 1 },
          { 17, 5, 1 },
          { 5, 18, 1 },
          { 5, 6, 1 },
          { 15, 8, 1 },
          { 10, 11, 1 },
          { 10, 8, 1 },
          { 10, 12, 1 },
          { 11, 19, 1 },
          { 8, 6, 1 },
          { 8, 9, 1 },
          { 6, 7, 1 },
          { 6, 16, 1 },
          { 12, 7, 1 },
          { 12, 13, 1 },
          { 0, 0, 0}
        }
      },
      { 'D', "GAU", // D GAU (G A U) D -> H [As K Above]
        {
          { "P",  "P",  -0.000,  8.791,  -0.000},
          { "O1P",  "O",  1.422,  9.311,  0.003},
          { "O2P",  "O",  -0.360,  7.848,  -0.133},
          { "O5'",  "O",  -0.323,  8.218,  0.173},
          { "C5'",  "C",  -0.483,  9.385,  0.274},
          { "C4'",  "C",  0.057,  9.172,  0.426},
          { "O4'",  "O",  -0.999,  8.642,  0.501},
          { "C3'",  "C",  1.275,  8.307,  0.458},
          { "O3'",  "O",  2.013,  8.747,  0.588},
          { "C2'",  "C",  0.606,  6.917,  0.468},
          { "O2'",  "O",  1.405,  6.122,  0.603},
          { "C1'",  "C",  -0.731,  7.281,  0.552},
          { "H1'",  "H",  -0.699,  7.317,  0.703},
          { "H2'",  "H",  0.484,  6.539,  0.315},
          { "H3'",  "H",  1.968,  8.387,  0.356},
          { "H4'",  "H",  0.274,  10.189,  0.463},
          { "H5'1",  "H",  0.084,  10.229,  0.231},
          { "H5'2",  "H",  -1.546,  9.658,  0.281},
          { "HO2'",  "H",  1.093,  6.345,  0.742},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 3, 1, 2 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 4, 5, 1 },
          { 14, 10, 1 },
          { 17, 5, 1 },
          { 5, 18, 1 },
          { 5, 6, 1 },
          { 15, 8, 1 },
          { 10, 11, 1 },
          { 10, 8, 1 },
          { 10, 12, 1 },
          { 11, 19, 1 },
          { 8, 6, 1 },
          { 8, 9, 1 },
          { 6, 7, 1 },
          { 6, 16, 1 },
          { 12, 7, 1 },
          { 12, 13, 1 },
          { 0, 0, 0}
        }
      },
      { 'H', "ACU", // H ACU (A C U) H -> D [As K Above]
        {
          { "P",  "P",  -0.000,  8.791,  -0.000},
          { "O1P",  "O",  1.422,  9.311,  0.003},
          { "O2P",  "O",  -0.360,  7.848,  -0.133},
          { "O5'",  "O",  -0.323,  8.218,  0.173},
          { "C5'",  "C",  -0.483,  9.385,  0.274},
          { "C4'",  "C",  0.057,  9.172,  0.426},
          { "O4'",  "O",  -0.999,  8.642,  0.501},
          { "C3'",  "C",  1.275,  8.307,  0.458},
          { "O3'",  "O",  2.013,  8.747,  0.588},
          { "C2'",  "C",  0.606,  6.917,  0.468},
          { "O2'",  "O",  1.405,  6.122,  0.603},
          { "C1'",  "C",  -0.731,  7.281,  0.552},
          { "H1'",  "H",  -0.699,  7.317,  0.703},
          { "H2'",  "H",  0.484,  6.539,  0.315},
          { "H3'",  "H",  1.968,  8.387,  0.356},
          { "H4'",  "H",  0.274,  10.189,  0.463},
          { "H5'1",  "H",  0.084,  10.229,  0.231},
          { "H5'2",  "H",  -1.546,  9.658,  0.281},
          { "HO2'",  "H",  1.093,  6.345,  0.742},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 3, 1, 2 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 4, 5, 1 },
          { 14, 10, 1 },
          { 17, 5, 1 },
          { 5, 18, 1 },
          { 5, 6, 1 },
          { 15, 8, 1 },
          { 10, 11, 1 },
          { 10, 8, 1 },
          { 10, 12, 1 },
          { 11, 19, 1 },
          { 8, 6, 1 },
          { 8, 9, 1 },
          { 6, 7, 1 },
          { 6, 16, 1 },
          { 12, 7, 1 },
          { 12, 13, 1 },
          { 0, 0, 0}
        }
      },
      { 'V', "GCA", // V GCA (G C A) V -> B [As K Above]
        {
          { "P",  "P",  -0.000,  8.791,  -0.000},
          { "O1P",  "O",  1.422,  9.311,  0.003},
          { "O2P",  "O",  -0.360,  7.848,  -0.133},
          { "O5'",  "O",  -0.323,  8.218,  0.173},
          { "C5'",  "C",  -0.483,  9.385,  0.274},
          { "C4'",  "C",  0.057,  9.172,  0.426},
          { "O4'",  "O",  -0.999,  8.642,  0.501},
          { "C3'",  "C",  1.275,  8.307,  0.458},
          { "O3'",  "O",  2.013,  8.747,  0.588},
          { "C2'",  "C",  0.606,  6.917,  0.468},
          { "O2'",  "O",  1.405,  6.122,  0.603},
          { "C1'",  "C",  -0.731,  7.281,  0.552},
          { "H1'",  "H",  -0.699,  7.317,  0.703},
          { "H2'",  "H",  0.484,  6.539,  0.315},
          { "H3'",  "H",  1.968,  8.387,  0.356},
          { "H4'",  "H",  0.274,  10.189,  0.463},
          { "H5'1",  "H",  0.084,  10.229,  0.231},
          { "H5'2",  "H",  -1.546,  9.658,  0.281},
          { "HO2'",  "H",  1.093,  6.345,  0.742},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 3, 1, 2 },
          { 1, 2, 1 },
          { 1, 4, 1 },
          { 4, 5, 1 },
          { 14, 10, 1 },
          { 17, 5, 1 },
          { 5, 18, 1 },
          { 5, 6, 1 },
          { 15, 8, 1 },
          { 10, 11, 1 },
          { 10, 8, 1 },
          { 10, 12, 1 },
          { 11, 19, 1 },
          { 8, 6, 1 },
          { 8, 9, 1 },
          { 6, 7, 1 },
          { 6, 16, 1 },
          { 12, 7, 1 },
          { 12, 13, 1 },
          { 0, 0, 0}
        }
      }
    };

  ResidueRecord ProteinResidues[IUPAC_Protein_max] =
    {
      { 0, "", // Protein Start
        {
          { "HNCA",  "H",  -0.908,  1.986,  -0.059},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 0, 0, 0}
        }
      },
      { 0, "", // Protein End
        {
          { "OXT",  "O",  1.517,  1.615,  1.659},
          { "HOCA",  "H",  2.208,  1.937,  1.983},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 1, 2, 1},
          { 0, 0, 0}
        }
      },
      { 'X', "XAA", // X XAA any (unknown) amino acid
        {
          { "N",  "N",  0.000,  1.576,  0.000},
          { "CA",  "C",  0.963,  2.292,  0.440},
          { "C",  "C",  1.963,  1.608,  0.913},
          { "O",  "O",  3.153,  1.742,  0.774},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 1, 2, 1 },
          { 2, 3, 1 },
          { 4, 3, 2 },
          { 0, 0, 0}
        }
      },
      { 'A', "ALA", // A ALA Alanine
        {
          { "N",  "N",  0.000,  1.576,  0.000},
          { "CA",  "C",  0.963,  2.292,  0.440},
          { "C",  "C",  1.963,  1.608,  0.913},
          { "O",  "O",  3.153,  1.742,  0.774},
          { "CB",  "C",  0.285,  3.409,  0.741},
          { "HA",  "H",  1.478,  2.999,  0.183},
          { "H",  "H",  -0.969,  1.563,  0.159},
          { "HB3",  "H",  -0.234,  3.243,  1.029},
          { "HB2",  "H",  -0.449,  3.921,  0.565},
          { "HB1",  "H",  1.041,  4.176,  0.802},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 1, 7, 1 },
          { 1, 2, 1 },
          { 6, 2, 1 },
          { 2, 3, 1 },
          { 2, 5, 1 },
          { 4, 3, 2 },
          { 9, 5, 1 },
          { 5, 8, 1 },
          { 5, 10, 1 },
          { 0, 0, 0}
        }
      },
      { 'B', "ASX", // B ASX Asparagine or Aspartic acid i.e. N or D
        {
          { "N",  "N",  0.000,  1.576,  0.000},
          { "CA",  "C",  1.277,  1.956,  0.378},
          { "C",  "C",  1.657,  2.111,  1.142},
          { "O",  "O",  1.005,  2.285,  1.603},
          { "H",  "H",  -0.460,  0.827,  0.439},
          { "HA",  "H",  2.086,  1.697,  -0.015},
          { "CB",  "C",  1.215,  3.490,  0.329},
          { "CG",  "C",  2.563,  4.246,  0.369},
          { "OD2",  "O",  2.446,  5.465,  0.396},
          { "HB1",  "H",  0.858,  3.913,  0.070},
          { "HB2",  "H",  0.481,  3.976,  0.508},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 6, 2, 1 },
          { 1, 5, 1 },
          { 1, 2, 1 },
          { 10, 7, 1 },
          { 2, 7, 1 },
          { 2, 3, 1 },
          { 7, 8, 1 },
          { 7, 11, 1 },
          { 8, 9, 1 },
          { 3, 4, 2 },
          { 0, 0, 0}
        }
      },
      { 'C', "CYS", // C CYS Cysteine
        {
          { "N",  "N",  0.000,  1.576,  0.000},
          { "CA",  "C",  1.240,  2.044,  0.344},
          { "C",  "C",  1.604,  1.920,  1.098},
          { "O",  "O",  2.805,  2.688,  1.146},
          { "CB",  "C",  1.076,  3.574,  0.377},
          { "SG",  "S",  2.635,  4.388,  0.481},
          { "H",  "H",  -0.547,  0.711,  0.326},
          { "HA",  "H",  2.055,  1.955,  -0.029},
          { "HB1",  "H",  0.759,  4.104,  0.137},
          { "HB2",  "H",  0.290,  3.934,  0.563},
          { "HG",  "H",  2.305,  4.962,  0.729},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 8, 2, 1 },
          { 1, 7, 1 },
          { 1, 2, 1 },
          { 9, 5, 1 },
          { 2, 5, 1 },
          { 2, 3, 1 },
          { 5, 6, 1 },
          { 5, 10, 1 },
          { 3, 4, 2 },
          { 6, 11, 1 },
          { 0, 0, 0}
        }
      },
      { 'D', "ASP", // D ASP Aspartic acid
        {
          { "N",  "N",  0.000,  1.576,  0.000},
          { "CA",  "C",  1.277,  1.956,  0.378},
          { "C",  "C",  1.657,  2.111,  1.142},
          { "O",  "O",  1.005,  2.285,  1.603},
          { "CB",  "C",  1.215,  3.490,  0.329},
          { "CG",  "C",  2.563,  4.246,  0.369},
          { "OD1",  "O",  3.587,  3.562,  0.350},
          { "OD2",  "O",  2.446,  5.465,  0.396},
          { "H",  "H",  -0.460,  0.827,  0.439},
          { "HA",  "H",  2.086,  1.697,  -0.015},
          { "HB1",  "H",  0.858,  3.913,  0.070},
          { "HB2",  "H",  0.481,  3.976,  0.508},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 10, 2, 1 },
          { 1, 9, 1 },
          { 1, 2, 1 },
          { 11, 5, 1 },
          { 2, 5, 1 },
          { 2, 3, 1 },
          { 5, 6, 1 },
          { 5, 12, 1 },
          { 7, 6, 2 },
          { 6, 8, 1 },
          { 3, 4, 2 },
          { 0, 0, 0}
        }
      },
      { 'E', "GLU", // E GLU Glutamic acid
        {
          { "N",  "N",  0.000,  1.576,  0.000},
          { "CA",  "C",  0.910,  2.346,  0.444},
          { "C",  "C",  1.926,  1.740,  0.913},
          { "O",  "O",  3.102,  2.014,  0.814},
          { "CB",  "C",  0.167,  3.467,  0.718},
          { "CG",  "C",  -0.772,  4.293,  0.478},
          { "CD",  "C",  -1.539,  5.337,  0.640},
          { "OE1",  "O",  -1.316,  5.704,  0.855},
          { "OE2",  "O",  -2.382,  6.019,  0.527},
          { "HA",  "H",  1.444,  3.037,  0.191},
          { "H",  "H",  -0.971,  1.494,  0.147},
          { "HB2",  "H",  0.911,  4.214,  0.802},
          { "HB1",  "H",  -0.414,  3.237,  0.989},
          { "HG2",  "H",  -1.507,  3.741,  0.327},
          { "HG1",  "H",  -0.174,  4.922,  0.331},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 1, 11, 1 },
          { 1, 2, 1 },
          { 10, 2, 1 },
          { 2, 3, 1 },
          { 2, 5, 1 },
          { 14, 6, 1 },
          { 3, 4, 2 },
          { 15, 6, 1 },
          { 6, 5, 1 },
          { 6, 7, 1 },
          { 5, 13, 1 },
          { 5, 12, 1 },
          { 9, 7, 2 },
          { 7, 8, 1 },
          { 0, 0, 0}
        }
      },
      { 'F', "PHE", // F PHE Phenylalanine
        {
          { "N",  "N",  0.000,  1.576,  0.000},
          { "CA",  "C",  1.247,  2.025,  0.345},
          { "C",  "C",  1.575,  1.945,  1.105},
          { "O",  "O",  0.830,  1.957,  1.633},
          { "CB",  "C",  1.129,  3.574,  0.392},
          { "CG",  "C",  0.795,  4.509,  0.099},
          { "CD1",  "C",  -0.543,  4.812,  0.040},
          { "CD2",  "C",  1.809,  5.212,  -0.036},
          { "CE1",  "C",  -0.863,  5.867,  -0.120},
          { "CE2",  "C",  1.489,  6.326,  -0.170},
          { "CZ",  "C",  0.153,  6.675,  -0.201},
          { "H",  "H",  -0.577,  0.741,  0.345},
          { "HA",  "H",  2.070,  1.915,  -0.029},
          { "HB1",  "H",  0.375,  3.918,  0.590},
          { "HB2",  "H",  2.072,  4.010,  0.494},
          { "HD1",  "H",  -1.338,  4.286,  0.153},
          { "HD2",  "H",  2.846,  4.963,  -0.001},
          { "HE1",  "H",  -1.898,  6.114,  -0.150},
          { "HE2",  "H",  2.276,  6.964,  -0.227},
          { "HZ",  "H",  -0.097,  7.602,  -0.270},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 20, 11, 1 },
          { 19, 10, 1 },
          { 11, 10, 2 },
          { 11, 9, 1 },
          { 10, 8, 1 },
          { 18, 9, 1 },
          { 9, 7, 2 },
          { 8, 17, 1 },
          { 8, 6, 2 },
          { 13, 2, 1 },
          { 1, 12, 1 },
          { 1, 2, 1 },
          { 7, 6, 1 },
          { 7, 16, 1 },
          { 6, 5, 1 },
          { 2, 5, 1 },
          { 2, 3, 1 },
          { 5, 15, 1 },
          { 5, 14, 1 },
          { 3, 4, 2 },
          { 0, 0, 0}
        }
      },
      { 'G', "GLY", // G GLY Glycine
        {
          { "N",  "N",  0.000,  1.576,  0.000},
          { "CA",  "C",  1.248,  2.021,  0.345},
          { "C",  "C",  1.569,  1.950,  1.106},
          { "O",  "O",  0.839,  1.959,  1.616},
          { "H",  "H",  -0.581,  0.747,  0.349},
          { "HA1",  "H",  2.081,  1.949,  -0.021},
          { "HA2",  "H",  1.177,  3.121,  0.379},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 6, 2, 1 },
          { 1, 5, 1 },
          { 1, 2, 1 },
          { 2, 7, 1 },
          { 2, 3, 1 },
          { 3, 4, 2 },
          { 0, 0, 0}
        }
      },
      { 'H', "HIS", // H HIS Histidine
        {
          { "N",  "N",  0.000,  1.576,  0.000},
          { "CA",  "C",  0.950,  2.280,  0.454},
          { "C",  "C",  1.934,  1.609,  0.944},
          { "O",  "O",  3.126,  1.770,  0.820},
          { "CB",  "C",  0.255,  3.433,  0.731},
          { "CG",  "C",  -0.659,  4.241,  0.496},
          { "ND1",  "N",  -0.268,  5.058,  0.290},
          { "CD2",  "C",  -1.970,  4.460,  0.547},
          { "CE1",  "C",  -1.269,  5.817,  0.232},
          { "NE2",  "N",  -2.293,  5.436,  0.342},
          { "HA",  "H",  1.509,  2.939,  0.189},
          { "H",  "H",  -0.977,  1.624,  0.127},
          { "HB2",  "H",  1.020,  4.169,  0.809},
          { "HB1",  "H",  -0.320,  3.247,  1.008},
          { "HD1",  "H",  0.691,  5.179,  0.219},
          { "HD2",  "H",  -2.632,  4.175,  0.736},
          { "HE1",  "H",  -1.242,  6.719,  0.134},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 1, 12, 1 },
          { 1, 2, 1 },
          { 11, 2, 1 },
          { 17, 9, 1 },
          { 2, 3, 1 },
          { 2, 5, 1 },
          { 15, 7, 1 },
          { 4, 3, 2 },
          { 9, 7, 1 },
          { 9, 10, 2 },
          { 7, 6, 1 },
          { 10, 8, 1 },
          { 6, 5, 1 },
          { 6, 8, 2 },
          { 5, 14, 1 },
          { 5, 13, 1 },
          { 8, 16, 1 },
          { 0, 0, 0}
        }
      },
      { 'I', "ILE", // I ILE Isoleucine
        {
          { "N",  "N",  0.000,  1.576,  0.000},
          { "CA",  "C",  1.249,  2.017,  0.346},
          { "C",  "C",  1.562,  1.956,  1.107},
          { "O",  "O",  2.768,  2.715,  1.159},
          { "CB",  "C",  1.224,  3.599,  0.360},
          { "CG1",  "C",  0.832,  4.480,  0.053},
          { "CG2",  "C",  2.576,  4.245,  0.472},
          { "CD1",  "C",  0.483,  5.914,  0.141},
          { "H",  "H",  -0.587,  0.754,  0.351},
          { "HA",  "H",  2.064,  1.835,  -0.033},
          { "HB",  "H",  0.446,  3.960,  0.548},
          { "HG11",  "H",  1.620,  4.608,  -0.116},
          { "HG12",  "H",  -0.065,  4.164,  -0.078},
          { "HG21",  "H",  2.532,  5.347,  0.462},
          { "HG22",  "H",  2.875,  4.163,  0.725},
          { "HG23",  "H",  3.408,  3.978,  0.307},
          { "HD11",  "H",  0.147,  6.555,  0.006},
          { "HD12",  "H",  -0.333,  5.998,  0.265},
          { "HD13",  "H",  1.346,  6.470,  0.207},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 12, 6, 1 },
          { 13, 6, 1 },
          { 10, 2, 1 },
          { 1, 9, 1 },
          { 1, 2, 1 },
          { 17, 8, 1 },
          { 6, 8, 1 },
          { 6, 5, 1 },
          { 2, 5, 1 },
          { 2, 3, 1 },
          { 8, 19, 1 },
          { 8, 18, 1 },
          { 16, 7, 1 },
          { 5, 7, 1 },
          { 5, 11, 1 },
          { 3, 4, 2 },
          { 7, 14, 1 },
          { 7, 15, 1 },
          { 0, 0, 0}
        }
      },
      { 'J', "XLE", // J XLE Leucine or Isoleucine i.e. L or I
        {
          { "N",  "N",  0.000,  1.576,  0.000},
          { "CA",  "C",  1.249,  2.017,  0.346},
          { "C",  "C",  1.562,  1.956,  1.107},
          { "O",  "O",  2.768,  2.715,  1.159},
          { "CB",  "C",  1.224,  3.599,  0.360},
          { "CG1",  "C",  0.832,  4.480,  0.053},
          { "CD1",  "C",  0.483,  5.914,  0.141},
          { "H",  "H",  -0.587,  0.754,  0.351},
          { "HA",  "H",  2.064,  1.835,  -0.033},
          { "HB",  "H",  0.446,  3.960,  0.548},
          { "HG11",  "H",  1.620,  4.608,  -0.116},
          { "HG12",  "H",  -0.065,  4.164,  -0.078},
          { "HD11",  "H",  0.147,  6.555,  0.006},
          { "HD12",  "H",  -0.333,  5.998,  0.265},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 11, 6, 1 },
          { 12, 6, 1 },
          { 9, 2, 1 },
          { 1, 8, 1 },
          { 1, 2, 1 },
          { 13, 7, 1 },
          { 6, 7, 1 },
          { 6, 5, 1 },
          { 2, 5, 1 },
          { 2, 3, 1 },
          { 7, 14, 1 },
          { 5, 10, 1 },
          { 3, 4, 2 },
          { 0, 0, 0}
        }
      },
      { 'K', "LYS", // K LYS Lysine
        {
          { "N",  "N",  0.000,  1.576,  0.000},
          { "CA",  "C",  0.947,  2.308,  0.444},
          { "C",  "C",  1.917,  1.661,  0.940},
          { "O",  "O",  3.102,  1.908,  0.855},
          { "CB",  "C",  0.242,  3.451,  0.721},
          { "CG",  "C",  -0.645,  4.320,  0.478},
          { "CD",  "C",  -1.510,  5.312,  0.650},
          { "CE",  "C",  -2.471,  6.121,  0.490},
          { "NZ",  "N",  -3.334,  7.021,  0.612},
          { "HA",  "H",  1.506,  2.980,  0.186},
          { "H",  "H",  -0.957,  1.471,  0.178},
          { "HB2",  "H",  1.001,  4.172,  0.813},
          { "HB1",  "H",  -0.369,  3.228,  0.988},
          { "HG2",  "H",  -1.317,  3.804,  0.304},
          { "HG1",  "H",  -0.005,  4.988,  0.348},
          { "HD2",  "H",  -0.851,  6.067,  0.731},
          { "HD1",  "H",  -2.112,  4.877,  0.810},
          { "HE2",  "H",  -3.107,  5.498,  0.378},
          { "HE1",  "H",  -1.882,  6.791,  0.389},
          { "HZ3",  "H",  -3.936,  6.560,  0.725},
          { "HZ2",  "H",  -4.019,  7.560,  0.524},
          { "HZ1",  "H",  -2.751,  7.795,  0.676},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 1, 11, 1 },
          { 1, 2, 1 },
          { 10, 2, 1 },
          { 2, 3, 1 },
          { 2, 5, 1 },
          { 14, 6, 1 },
          { 3, 4, 2 },
          { 15, 6, 1 },
          { 6, 5, 1 },
          { 6, 7, 1 },
          { 18, 8, 1 },
          { 5, 13, 1 },
          { 5, 12, 1 },
          { 19, 8, 1 },
          { 8, 7, 1 },
          { 8, 9, 1 },
          { 7, 17, 1 },
          { 7, 16, 1 },
          { 21, 9, 1 },
          { 9, 20, 1 },
          { 9, 22, 1 },
          { 0, 0, 0}
        }
      },
      { 'L', "LEU", // L LEU Leucine
        {
          { "N",  "N",  0.000,  1.576,  0.000},
          { "CA",  "C",  0.930,  2.335,  0.438},
          { "C",  "C",  1.911,  1.701,  0.925},
          { "O",  "O",  3.096,  1.942,  0.837},
          { "CB",  "C",  0.226,  3.473,  0.717},
          { "CG",  "C",  -0.820,  4.280,  0.501},
          { "CD1",  "C",  -1.494,  5.347,  0.691},
          { "CD2",  "C",  -0.194,  5.223,  0.275},
          { "HA",  "H",  1.500,  2.994,  0.182},
          { "H",  "H",  -0.954,  1.440,  0.181},
          { "HB2",  "H",  0.999,  4.214,  0.785},
          { "HB1",  "H",  -0.277,  3.266,  1.002},
          { "HG",  "H",  -1.618,  3.666,  0.387},
          { "HD13",  "H",  -1.970,  5.030,  0.873},
          { "HD12",  "H",  -2.270,  5.876,  0.589},
          { "HD11",  "H",  -0.746,  6.106,  0.740},
          { "HD23",  "H",  0.668,  5.733,  0.359},
          { "HD22",  "H",  -0.935,  6.001,  0.232},
          { "HD21",  "H",  0.135,  4.849,  0.080},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 1, 10, 1 },
          { 1, 2, 1 },
          { 19, 8, 1 },
          { 9, 2, 1 },
          { 2, 3, 1 },
          { 2, 5, 1 },
          { 3, 4, 2 },
          { 18, 8, 1 },
          { 13, 6, 1 },
          { 8, 17, 1 },
          { 8, 6, 1 },
          { 6, 5, 1 },
          { 6, 7, 1 },
          { 5, 12, 1 },
          { 5, 11, 1 },
          { 15, 7, 1 },
          { 7, 14, 1 },
          { 7, 16, 1 },
          { 0, 0, 0}
        }
      },
      { 'M', "MET", // M MET Methionine
        {
          { "N",  "N",  0.000,  1.576,  0.000},
          { "CA",  "C",  0.933,  2.282,  0.466},
          { "C",  "C",  1.911,  1.626,  0.963},
          { "O",  "O",  3.095,  1.885,  0.883},
          { "CB",  "C",  0.222,  3.387,  0.768},
          { "CG",  "C",  -0.732,  4.265,  0.542},
          { "SD",  "S",  0.101,  5.355,  0.289},
          { "CE",  "C",  1.257,  6.259,  0.476},
          { "HA",  "H",  1.490,  2.973,  0.214},
          { "H",  "H",  -0.977,  1.544,  0.137},
          { "HB2",  "H",  0.981,  4.105,  0.862},
          { "HB1",  "H",  -0.358,  3.141,  1.046},
          { "HG2",  "H",  -1.165,  5.040,  0.683},
          { "HG1",  "H",  -1.563,  3.740,  0.418},
          { "HE3",  "H",  0.729,  6.775,  0.601},
          { "HE2",  "H",  1.672,  7.081,  0.385},
          { "HE1",  "H",  2.080,  5.645,  0.541},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 1, 10, 1 },
          { 1, 2, 1 },
          { 9, 2, 1 },
          { 2, 3, 1 },
          { 2, 5, 1 },
          { 3, 4, 2 },
          { 14, 6, 1 },
          { 7, 6, 1 },
          { 7, 8, 1 },
          { 6, 5, 1 },
          { 6, 13, 1 },
          { 5, 12, 1 },
          { 5, 11, 1 },
          { 16, 8, 1 },
          { 8, 17, 1 },
          { 8, 15, 1 },
          { 0, 0, 0}
        }
      },
      { 'N', "ASN", // N ASN Asparagine
        {
          { "N",  "N",  0.000,  1.576,  0.000},
          { "CA",  "C",  0.948,  2.326,  0.449},
          { "C",  "C",  1.942,  1.694,  0.928},
          { "O",  "O",  3.148,  1.800,  0.841},
          { "H",  "H",  -0.477,  0.752,  0.474},
          { "HA",  "H",  1.494,  2.984,  0.184},
          { "CB",  "C",  0.243,  3.477,  0.720},
          { "HB1",  "H",  -0.434,  3.998,  0.534},
          { "HB2",  "H",  1.006,  4.240,  0.780},
          { "CG",  "C",  -0.535,  3.369,  1.100},
          { "OD2",  "O",  -0.005,  4.017,  1.342},
          { "ND1",  "N",  -1.799,  2.906,  1.124},
          { "HD1",  "H",  -2.252,  2.576,  0.819},
          { "HD2",  "H",  -2.294,  3.144,  1.401},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 1, 5, 1 },
          { 1, 2, 1 },
          { 6, 2, 1 },
          { 2, 3, 1 },
          { 2, 7, 1 },
          { 4, 3, 2 },
          { 13, 12, 1 },
          { 8, 7, 1 },
          { 7, 9, 1 },
          { 7, 10, 1 },
          { 12, 10, 1 },
          { 12, 14, 1 },
          { 10, 11, 2 },
          { 0, 0, 0}
        }
      },
      { 'O', "PYL", // O PYL Pyrrolysine (K Lysine derivative)
        {
          { "N",  "N",  0.000,  1.576,  0.000},
          { "CA",  "C",  0.947,  2.308,  0.444},
          { "C",  "C",  1.917,  1.661,  0.940},
          { "O",  "O",  3.102,  1.908,  0.855},
          { "CB",  "C",  0.242,  3.451,  0.721},
          { "CG",  "C",  -0.645,  4.320,  0.478},
          { "CD",  "C",  -1.510,  5.312,  0.650},
          { "CE",  "C",  -2.471,  6.121,  0.490},
          { "NZ",  "N",  -3.334,  7.021,  0.612},
          { "HA",  "H",  1.506,  2.980,  0.186},
          { "H",  "H",  -0.957,  1.471,  0.178},
          { "HB2",  "H",  1.001,  4.172,  0.813},
          { "HB1",  "H",  -0.369,  3.228,  0.988},
          { "HG2",  "H",  -1.317,  3.804,  0.304},
          { "HG1",  "H",  -0.005,  4.988,  0.348},
          { "HD2",  "H",  -0.851,  6.067,  0.731},
          { "HD1",  "H",  -2.112,  4.877,  0.810},
          { "HE2",  "H",  -3.107,  5.498,  0.378},
          { "HE1",  "H",  -1.882,  6.791,  0.389},
          { "HZ3",  "H",  -3.936,  6.560,  0.725},
          { "HZ2",  "H",  -4.019,  7.560,  0.524},
          { "HZ1",  "H",  -2.751,  7.795,  0.676},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 1, 11, 1 },
          { 1, 2, 1 },
          { 10, 2, 1 },
          { 2, 3, 1 },
          { 2, 5, 1 },
          { 14, 6, 1 },
          { 3, 4, 2 },
          { 15, 6, 1 },
          { 6, 5, 1 },
          { 6, 7, 1 },
          { 18, 8, 1 },
          { 5, 13, 1 },
          { 5, 12, 1 },
          { 19, 8, 1 },
          { 8, 7, 1 },
          { 8, 9, 1 },
          { 7, 17, 1 },
          { 7, 16, 1 },
          { 21, 9, 1 },
          { 9, 20, 1 },
          { 9, 22, 1 },
          { 0, 0, 0}
        }
      },
      { 'P', "PRO", // P PRO Proline
        {
          { "N",  "N",  0.000,  1.576,  0.000},
          { "CA",  "C",  1.411,  1.867,  0.295},
          { "C",  "C",  1.762,  1.962,  1.090},
          { "O",  "O",  0.887,  2.164,  1.520},
          { "CB",  "C",  1.482,  3.388,  0.297},
          { "CG",  "C",  0.041,  3.764,  0.404},
          { "CD",  "C",  -0.812,  2.773,  0.156},
          { "HA",  "H",  2.139,  1.587,  -0.143},
          { "HB1",  "H",  2.160,  3.805,  0.510},
          { "HB2",  "H",  1.672,  3.881,  0.031},
          { "HG1",  "H",  -0.227,  3.726,  0.690},
          { "HG2",  "H",  -0.114,  4.791,  0.320},
          { "HD1",  "H",  -1.750,  2.547,  0.358},
          { "HD2",  "H",  -0.967,  3.433,  -0.126},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 14, 7, 1 },
          { 8, 2, 1 },
          { 1, 7, 1 },
          { 1, 2, 1 },
          { 10, 5, 1 },
          { 7, 13, 1 },
          { 7, 6, 1 },
          { 2, 5, 1 },
          { 2, 3, 1 },
          { 5, 6, 1 },
          { 5, 9, 1 },
          { 6, 12, 1 },
          { 6, 11, 1 },
          { 3, 4, 2 },
          { 0, 0, 0}
        }
      },
      { 'Q', "GLN", // Q GLN Glutamine
        {
          { "N",  "N",  0.000,  1.576,  0.000},
          { "CA",  "C",  0.932,  2.292,  0.461},
          { "C",  "C",  1.945,  1.628,  0.932},
          { "O",  "O",  3.134,  1.785,  0.797},
          { "CB",  "C",  0.234,  3.389,  0.768},
          { "CG",  "C",  -0.593,  4.374,  0.542},
          { "CD",  "C",  -1.534,  5.262,  0.707},
          { "OE1",  "O",  -2.508,  5.692,  0.597},
          { "NE2",  "N",  -1.306,  5.805,  0.928},
          { "HA",  "H",  1.459,  2.997,  0.206},
          { "H",  "H",  -0.980,  1.619,  0.116},
          { "HB2",  "H",  1.012,  4.047,  0.880},
          { "HB1",  "H",  -0.418,  3.126,  1.029},
          { "HG2",  "H",  -1.198,  3.942,  0.347},
          { "HG1",  "H",  0.090,  5.088,  0.439},
          { "HE21",  "H",  -0.512,  5.657,  1.031},
          { "HE22",  "H",  -1.955,  6.518,  0.972},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 1, 11, 1 },
          { 1, 2, 1 },
          { 10, 2, 1 },
          { 2, 3, 1 },
          { 2, 5, 1 },
          { 4, 3, 2 },
          { 14, 6, 1 },
          { 15, 6, 1 },
          { 6, 5, 1 },
          { 6, 7, 1 },
          { 5, 13, 1 },
          { 5, 12, 1 },
          { 8, 7, 2 },
          { 7, 9, 1 },
          { 9, 16, 1 },
          { 9, 17, 1 },
          { 0, 0, 0}
        }
      },
      { 'R', "ARG", // R ARG Arginine
        {
          { "N",  "N",  0.000,  1.576,  0.000},
          { "CA",  "C",  1.416,  1.596,  0.225},
          { "C",  "C",  1.616,  2.015,  1.050},
          { "O",  "O",  0.672,  2.525,  1.316},
          { "CB",  "C",  1.985,  3.007,  0.077},
          { "CG",  "C",  3.533,  3.130,  0.098},
          { "CD",  "C",  4.008,  4.580,  0.046},
          { "NE",  "N",  5.499,  4.632,  0.058},
          { "CZ",  "C",  6.234,  5.734,  0.038},
          { "NH1",  "N",  5.730,  6.928,  0.017},
          { "NH2",  "N",  7.523,  5.621,  0.048},
          { "H",  "H",  -0.326,  2.124,  -0.443},
          { "HA",  "H",  1.937,  0.948,  -0.377},
          { "HB1",  "H",  1.643,  3.517,  -0.208},
          { "HB2",  "H",  1.517,  3.772,  0.272},
          { "HG1",  "H",  3.904,  2.951,  0.439},
          { "HG2",  "H",  3.985,  2.558,  -0.199},
          { "HD1",  "H",  3.629,  5.079,  -0.145},
          { "HD2",  "H",  3.559,  5.266,  0.199},
          { "HE",  "H",  6.050,  3.779,  0.092},
          { "HH11",  "H",  6.363,  7.728,  0.011},
          { "HH12",  "H",  4.708,  6.950,  0.012},
          { "HH21",  "H",  7.899,  4.681,  0.074},
          { "HH22",  "H",  8.072,  6.480,  0.035},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 12, 1, 1 },
          { 18, 7, 1 },
          { 14, 5, 1 },
          { 17, 6, 1 },
          { 13, 2, 1 },
          { 1, 2, 1 },
          { 22, 10, 1 },
          { 21, 10, 1 },
          { 10, 9, 2 },
          { 7, 8, 1 },
          { 7, 6, 1 },
          { 7, 19, 1 },
          { 9, 8, 1 },
          { 9, 11, 1 },
          { 24, 11, 1 },
          { 5, 6, 1 },
          { 5, 2, 1 },
          { 5, 15, 1 },
          { 8, 20, 1 },
          { 11, 23, 1 },
          { 6, 16, 1 },
          { 2, 3, 1 },
          { 3, 4, 2 },
          { 0, 0, 0}
        }
      },
      { 'S', "SER", // S SER Serine
        {
          { "N",  "N",  0.000,  1.576,  0.000},
          { "CA",  "C",  1.225,  2.074,  0.344},
          { "C",  "C",  1.650,  1.881,  1.087},
          { "O",  "O",  1.008,  1.813,  1.663},
          { "CB",  "C",  1.056,  3.606,  0.401},
          { "OG",  "O",  0.160,  4.124,  0.658},
          { "H",  "H",  -0.499,  0.664,  0.284},
          { "HA",  "H",  2.031,  2.046,  -0.026},
          { "HB1",  "H",  2.041,  4.080,  0.449},
          { "HB2",  "H",  0.716,  4.159,  0.169},
          { "HG",  "H",  -0.707,  3.729,  0.625},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 8, 2, 1 },
          { 1, 7, 1 },
          { 1, 2, 1 },
          { 2, 5, 1 },
          { 2, 3, 1 },
          { 10, 5, 1 },
          { 5, 9, 1 },
          { 5, 6, 1 },
          { 3, 4, 2 },
          { 11, 6, 1 },
          { 0, 0, 0}
        }
      },
      { 'T', "THR", // T THR Threonine
        {
          { "N",  "N",  0.000,  1.576,  0.000},
          { "CA",  "C",  1.228,  2.069,  0.343},
          { "C",  "C",  1.643,  1.888,  1.089},
          { "O",  "O",  1.000,  1.825,  1.651},
          { "CB",  "C",  1.051,  3.615,  0.406},
          { "OG1",  "O",  0.772,  4.381,  0.112},
          { "CG2",  "C",  2.273,  4.433,  0.524},
          { "H",  "H",  -0.507,  0.671,  0.293},
          { "HA",  "H",  2.037,  2.045,  -0.025},
          { "HB",  "H",  0.192,  3.874,  0.581},
          { "HG1",  "H",  0.606,  5.245,  0.193},
          { "HG21",  "H",  2.080,  5.520,  0.508},
          { "HG22",  "H",  2.557,  4.392,  0.766},
          { "HG23",  "H",  3.163,  4.253,  0.378},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 9, 2, 1 },
          { 1, 8, 1 },
          { 1, 2, 1 },
          { 6, 11, 1 },
          { 6, 5, 1 },
          { 2, 5, 1 },
          { 2, 3, 1 },
          { 5, 10, 1 },
          { 5, 7, 1 },
          { 14, 7, 1 },
          { 3, 4, 2 },
          { 7, 12, 1 },
          { 7, 13, 1 },
          { 0, 0, 0}
        }
      },
      { 'U', "SEC", // U SEC Selenocysteine (CYS with Se in place of S)
        {
          { "N",  "N",  0.000,  1.576,  0.000},
          { "CA",  "C",  1.240,  2.044,  0.344},
          { "C",  "C",  1.604,  1.920,  1.098},
          { "O",  "O",  2.805,  2.688,  1.146},
          { "CB",  "C",  1.076,  3.574,  0.377},
          { "SEG",  "Se",  2.635,  4.388,  0.481},
          { "H",  "H",  -0.547,  0.711,  0.326},
          { "HA",  "H",  2.055,  1.955,  -0.029},
          { "HB1",  "H",  0.759,  4.104,  0.137},
          { "HB2",  "H",  0.290,  3.934,  0.563},
          { "HG",  "H",  2.305,  4.962,  0.729},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 8, 2, 1 },
          { 1, 7, 1 },
          { 1, 2, 1 },
          { 9, 5, 1 },
          { 2, 5, 1 },
          { 2, 3, 1 },
          { 5, 6, 1 },
          { 5, 10, 1 },
          { 3, 4, 2 },
          { 6, 11, 1 },
          { 0, 0, 0}
        }
      },
      { 'V', "VAL", // V VAL Valine
        {
          { "N",  "N",  0.000,  1.576,  0.000},
          { "CA",  "C",  1.246,  2.028,  0.345},
          { "C",  "C",  1.579,  1.942,  1.104},
          { "O",  "O",  0.854,  1.945,  1.626},
          { "CB",  "C",  1.137,  3.601,  0.378},
          { "CG1",  "C",  0.990,  4.498,  0.065},
          { "CG2",  "C",  2.309,  4.401,  0.538},
          { "H",  "H",  -0.572,  0.737,  0.343},
          { "HA",  "H",  2.070,  1.900,  -0.027},
          { "HB",  "H",  0.220,  3.885,  0.528},
          { "HG11",  "H",  1.892,  4.613,  -0.073},
          { "HG12",  "H",  0.135,  4.338,  -0.089},
          { "HG13",  "H",  0.806,  5.520,  0.142},
          { "HG21",  "H",  3.282,  4.161,  0.428},
          { "HG22",  "H",  2.158,  5.489,  0.509},
          { "HG23",  "H",  2.419,  4.421,  0.788},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 12, 6, 1 },
          { 11, 6, 1 },
          { 9, 2, 1 },
          { 1, 8, 1 },
          { 1, 2, 1 },
          { 6, 13, 1 },
          { 6, 5, 1 },
          { 2, 5, 1 },
          { 2, 3, 1 },
          { 5, 10, 1 },
          { 5, 7, 1 },
          { 14, 7, 1 },
          { 3, 4, 2 },
          { 7, 15, 1 },
          { 7, 16, 1 },
          { 0, 0, 0}
        }
      },
      { 'W', "TRP", // W TRP Tryptophan
        {
          { "N",  "N",  0.000,  1.576,  0.000},
          { "CA",  "C",  1.246,  2.034,  0.341},
          { "C",  "C",  1.601,  1.924,  1.099},
          { "O",  "O",  0.894,  1.916,  1.630},
          { "CB",  "C",  1.088,  3.574,  0.392},
          { "CG",  "C",  0.917,  4.497,  0.089},
          { "CD1",  "C",  -0.313,  5.053,  0.002},
          { "CD2",  "C",  1.885,  5.143,  -0.065},
          { "NE1",  "N",  -0.133,  6.139,  -0.156},
          { "CE2",  "C",  1.238,  6.200,  -0.187},
          { "CE3",  "C",  3.293,  4.979,  -0.080},
          { "CZ2",  "C",  1.996,  7.193,  -0.286},
          { "CZ3",  "C",  4.022,  5.914,  -0.215},
          { "CH2",  "C",  3.384,  7.044,  -0.298},
          { "H",  "H",  -0.517,  0.782,  0.288},
          { "HA",  "H",  2.067,  1.947,  -0.031},
          { "HB1",  "H",  0.238,  3.884,  0.564},
          { "HB2",  "H",  1.970,  4.047,  0.517},
          { "HD1",  "H",  -1.275,  4.736,  0.080},
          { "HE1",  "H",  -0.851,  6.770,  -0.213},
          { "HE3",  "H",  3.791,  4.228,  0.050},
          { "HZ2",  "H",  1.509,  8.081,  -0.334},
          { "HZ3",  "H",  5.095,  5.818,  -0.228},
          { "HH2",  "H",  3.975,  7.842,  -0.355},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 24, 14, 1 },
          { 22, 12, 1 },
          { 14, 12, 1 },
          { 14, 13, 2 },
          { 12, 10, 2 },
          { 20, 9, 1 },
          { 23, 13, 1 },
          { 13, 11, 1 },
          { 10, 9, 1 },
          { 10, 8, 1 },
          { 9, 7, 1 },
          { 11, 8, 2 },
          { 11, 21, 1 },
          { 8, 6, 1 },
          { 16, 2, 1 },
          { 1, 15, 1 },
          { 1, 2, 1 },
          { 7, 19, 1 },
          { 7, 6, 2 },
          { 6, 5, 1 },
          { 2, 5, 1 },
          { 2, 3, 1 },
          { 5, 18, 1 },
          { 5, 17, 1 },
          { 3, 4, 2 },
          { 0, 0, 0}
        }
      },
      { 'Y', "TYR", // Y TYR Tyrosine
        {
          { "N",  "N",  0.000,  1.576,  0.000},
          { "CA",  "C",  1.230,  2.066,  0.343},
          { "C",  "C",  1.638,  1.892,  1.091},
          { "O",  "O",  0.981,  1.837,  1.656},
          { "CB",  "C",  1.043,  3.595,  0.425},
          { "CG",  "C",  0.707,  4.570,  0.145},
          { "CD1",  "C",  1.704,  5.199,  -0.010},
          { "CD2",  "C",  -0.611,  5.023,  0.121},
          { "CE1",  "C",  1.383,  6.374,  -0.127},
          { "CE2",  "C",  -0.930,  6.111,  -0.025},
          { "CZ",  "C",  0.068,  6.830,  -0.126},
          { "OH",  "O",  -0.244,  8.039,  -0.201},
          { "H",  "H",  -0.512,  0.676,  0.298},
          { "HA",  "H",  2.042,  2.054,  -0.023},
          { "HB1",  "H",  0.272,  3.867,  0.627},
          { "HB2",  "H",  1.961,  4.059,  0.533},
          { "HD1",  "H",  2.727,  4.850,  -0.004},
          { "HD2",  "H",  -1.394,  4.610,  0.251},
          { "HE1",  "H",  2.150,  6.976,  -0.196},
          { "HE2",  "H",  -1.953,  6.455,  -0.029},
          { "HH",  "H",  -1.185,  8.185,  -0.187},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 12, 21, 1 },
          { 12, 11, 1 },
          { 19, 9, 1 },
          { 11, 9, 2 },
          { 11, 10, 1 },
          { 9, 7, 1 },
          { 20, 10, 1 },
          { 10, 8, 2 },
          { 7, 17, 1 },
          { 7, 6, 2 },
          { 14, 2, 1 },
          { 1, 13, 1 },
          { 1, 2, 1 },
          { 8, 6, 1 },
          { 8, 18, 1 },
          { 6, 5, 1 },
          { 2, 5, 1 },
          { 2, 3, 1 },
          { 5, 16, 1 },
          { 5, 15, 1 },
          { 3, 4, 2 },
          { 0, 0, 0}
        }
      },
      { 'Z', "GLX", // Z GLX Glutamine or Glutamic acid i.e. Q or E
        {
          { "N",  "N",  0.000,  1.576,  0.000},
          { "CA",  "C",  0.932,  2.292,  0.461},
          { "C",  "C",  1.945,  1.628,  0.932},
          { "O",  "O",  3.134,  1.785,  0.797},
          { "CB",  "C",  0.234,  3.389,  0.768},
          { "CG",  "C",  -0.593,  4.374,  0.542},
          { "CD",  "C",  -1.534,  5.262,  0.707},
          { "OE1",  "O",  -2.508,  5.692,  0.597},
          { "HA",  "H",  1.459,  2.997,  0.206},
          { "H",  "H",  -0.980,  1.619,  0.116},
          { "HB2",  "H",  1.012,  4.047,  0.880},
          { "HB1",  "H",  -0.418,  3.126,  1.029},
          { "HG2",  "H",  -1.198,  3.942,  0.347},
          { "HG1",  "H",  0.090,  5.088,  0.439},
          { "", "", 0.0, 0.0, 0.0 }
        },
        {
          { 1, 10, 1 },
          { 1, 2, 1 },
          { 9, 2, 1 },
          { 2, 3, 1 },
          { 2, 5, 1 },
          { 4, 3, 2 },
          { 13, 6, 1 },
          { 14, 6, 1 },
          { 6, 5, 1 },
          { 6, 7, 1 },
          { 5, 12, 1 },
          { 5, 11, 1 },
          { 8, 7, 2 },
          { 0, 0, 0}
        }
      }
    };

} //namespace OpenBabel
