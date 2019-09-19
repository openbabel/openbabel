/**********************************************************************
mmcifformat.cpp - Conversion to and from mmCIF format.
Copyright (C) Scarlet Line 2007

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
#include <openbabel/generic.h>

#include <openbabel/op.h>

#include <iostream>
#include <algorithm>
#include <ctype.h>

using namespace std;
namespace OpenBabel
{
 static const string UNKNOWN_VALUE = "?";

 class mmCIFFormat : public OBMoleculeFormat
 {
 public:
   //Register this format type ID
   mmCIFFormat()
   { // Copied from the Chemical MIME Page at http://www.ch.ic.ac.uk/chemime/
     OBConversion::RegisterFormat("mcif", this, "chemical/x-mmcif");
     OBConversion::RegisterFormat("mmcif", this, "chemical/x-mmcif");
     // Uncomment the following line, and this file will handle all CIF formats
     // OBConversion::RegisterFormat("cif", this, "chemical/x-cif");

     OBConversion::RegisterOptionParam("s", this);
     OBConversion::RegisterOptionParam("b", this);
     OBConversion::RegisterOptionParam("w", this);
   }

   virtual const char* Description() //required
   {
     return
       "Macromolecular Crystallographic Info\n "
       "Read Options e.g. -as\n"
       "  s  Output single bonds only\n"
       "  b  Disable bonding entirely\n"
       "  w  Wrap atomic coordinates into unit cell box\n\n";
   };

   virtual const char* SpecificationURL()
   { return "http://mmcif.pdb.org/";}; //optional
   // CIF itself is at http://www.iucr.org/iucr-top/cif/index.html

   virtual const char* GetMIMEType()
   { return "chemical/x-mmcif"; };

   //*** This section identical for most OBMol conversions ***
   ////////////////////////////////////////////////////
   /// The "API" interface functions
   virtual int SkipObjects(int n, OBConversion* pConv);
   virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
   virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
 };

 //Make an instance of the format class
 mmCIFFormat themmCIFFormat;

 struct CIFTagID
   {
   enum CIFCatName
     {
     unread_CIFCatName,
     atom_site,
     cell,
     chemical,
     chemical_formula,
     symmetry,
     symmetry_equiv,
     space_group,
     atom_type,
     MAX_CIFCatName
     };
   enum CIFDataName
     {
     unread_CIFDataName,
     _atom_site_fract_x, // The x coordinate specified as a fraction of _cell_length_a
     _atom_site_fract_y, // The y coordinate specified as a fraction of _cell_length_b
     _atom_site_fract_z, // The z coordinate specified as a fraction of _cell_length_c
     _atom_site_Cartn_x, // The x coordinate in angstroms
     _atom_site_Cartn_y, // The y coordinate in angstroms
     _atom_site_Cartn_z, // The z coordinate in angstroms
     _atom_site_label, // The atomic label if more detailed label info unavailable
     _atom_site_label_atom_id, // The atomic label within the residue
     _atom_site_label_comp_id, // The residue abbreviation, e.g. ILE
     _atom_site_label_entity_id, // The chain entity number of the residue, e.g. 2
     _atom_site_label_asym_id, // The unique chain id
     _atom_site_label_seq_id, // The sequence number of the residue, within the chain, e.g. 12
     _atom_site_type_symbol, // Atomic symbol, e.g. C
     _atom_site_occupancy,
     MAX_atom_site,
     _cell_length_a, // Unit-cell length a in Angstroms
     _cell_length_b, // Unit-cell length b in Angstroms
     _cell_length_c, // Unit-cell length c in Angstroms
     _cell_angle_alpha, // Unit-cell angle alpha in degrees
     _cell_angle_beta, // Unit-cell angle beta in degrees
     _cell_angle_gamma, // Unit-cell angle gamma in degrees
     MAX_cell,
     // chemical name organized by increasing desirability
     _chemical_name_common,
     _chemical_name_structure_type,
     _chemical_name_mineral,
     _chemical_name_systematic,
     MAX_chemical,
     // chemical formulae organized by increasing desirability
     _chemical_formula_moiety,
     _chemical_formula_iupac,
     _chemical_formula_structural,
     _chemical_formula_analytical,
     MAX_chemical_formula,
     _symmetry_Int_Tables_number,
     _symmetry_space_group_name_Hall,
     _symmetry_space_group_name_H_M,
     MAX_symmetry,
     _symmetry_equiv_pos_as_xyz,
     MAX_symmetry_equiv,
     _space_group_IT_number,
     _space_group_name_Hall,
     _space_group_name_H_M_alt,
     MAX_space_group,
     _atom_type_symbol,
     _atom_type_oxidation_number,
     MAX_atom_type,
     MAX_CIFDataName
     };
   char  tagname[76];
   CIFDataName tagid;
   };
 typedef vector<CIFTagID::CIFDataName> CIFColumnList;
 typedef map<string, CIFTagID::CIFDataName> CIFtagmap;
 struct CIFResidueID
   {
   unsigned long ChainNum; // The number of the chain
   unsigned long ResNum;  // The number of the residue within the chain
   CIFResidueID()
     {}
   CIFResidueID(unsigned long c, unsigned long r)
   :ChainNum(c), ResNum(r)
     {}
   CIFResidueID(const CIFResidueID & other)
   :ChainNum(other.ChainNum), ResNum(other.ResNum)
     {}
   CIFResidueID & operator=(const CIFResidueID & other)
     {
     ChainNum = other.ChainNum;
     ResNum = other.ResNum;
     return (* this);
     }
   bool operator< (const CIFResidueID & other) const
     {
     return ( ChainNum < other.ChainNum ? true : ( other.ChainNum < ChainNum ? false : ResNum < other.ResNum ) );
     }
   };
 typedef map<CIFResidueID, int> CIFResidueMap;
 CIFtagmap CIFtagLookupTable;

 CIFTagID CIFTagsRead[] =
   {
   { "_atom_site_fract_x", CIFTagID::_atom_site_fract_x },
   { "_atom_site_fract_y", CIFTagID::_atom_site_fract_y },
   { "_atom_site_fract_z", CIFTagID::_atom_site_fract_z },
   { "_atom_site_cartn_x", CIFTagID::_atom_site_Cartn_x },
   { "_atom_site_cartn_y", CIFTagID::_atom_site_Cartn_y },
   { "_atom_site_cartn_z", CIFTagID::_atom_site_Cartn_z },
   { "_atom_site_type_symbol", CIFTagID::_atom_site_type_symbol },
   { "_atom_site_occupancy", CIFTagID::_atom_site_occupancy},
   { "_atom_site_id", CIFTagID::_atom_site_label },
   { "_atom_site_label", CIFTagID::_atom_site_label },
   { "_atom_site_label_atom_id", CIFTagID::_atom_site_label_atom_id },
   { "_atom_site_label_comp_id", CIFTagID::_atom_site_label_comp_id },
   { "_atom_site_label_entity_id", CIFTagID::_atom_site_label_entity_id },
   { "_atom_site_label_seq_id", CIFTagID::_atom_site_label_seq_id },
   { "_atom_site_label_asym_id", CIFTagID::_atom_site_label_asym_id },
   { "_cell_length_a", CIFTagID::_cell_length_a },
   { "_cell_length_b", CIFTagID::_cell_length_b },
   { "_cell_length_c", CIFTagID::_cell_length_c },
   { "_cell_angle_alpha", CIFTagID::_cell_angle_alpha },
   { "_cell_angle_beta", CIFTagID::_cell_angle_beta },
   { "_cell_angle_gamma", CIFTagID::_cell_angle_gamma },
   { "_chemical_name_systematic", CIFTagID::_chemical_name_systematic },
   { "_chemical_name_mineral", CIFTagID::_chemical_name_mineral },
   { "_chemical_name_structure_type", CIFTagID::_chemical_name_structure_type },
   { "_chemical_name_common", CIFTagID::_chemical_name_common },
   { "_chemical_formula_analytical", CIFTagID::_chemical_formula_analytical },
   { "_chemical_formula_structural", CIFTagID::_chemical_formula_structural },
   { "_chemical_formula_iupac", CIFTagID::_chemical_formula_iupac },
   { "_chemical_formula_moiety", CIFTagID::_chemical_formula_moiety },
   { "_space_group_it_number", CIFTagID::_space_group_IT_number },
   { "_space_group_name_hall", CIFTagID::_space_group_name_Hall },
   { "_space_group_name_h-m_alt", CIFTagID::_space_group_name_H_M_alt },
   { "_symmetry_int_tables_number", CIFTagID::_symmetry_Int_Tables_number },
   { "_symmetry_space_group_name_hall", CIFTagID::_symmetry_space_group_name_Hall },
   { "_symmetry_space_group_name_h-m", CIFTagID::_symmetry_space_group_name_H_M },
   { "_symmetry_equiv_pos_as_xyz", CIFTagID::_symmetry_equiv_pos_as_xyz },
   { "_space_group_symop_operation_xyz", CIFTagID::_symmetry_equiv_pos_as_xyz },
   { "_atom_type_symbol", CIFTagID::_atom_type_symbol },
   { "_atom_type_oxidation_number",CIFTagID::_atom_type_oxidation_number },
   { "", CIFTagID::unread_CIFDataName }
   };

 class CIFLexer
 {
 public:
   enum TokenType
     {
     UnknownToken,
     KeyDataToken,
     KeyLoopToken,
     KeySaveToken,
     KeySaveEndToken,
     KeyStopToken,
     KeyGlobalToken,
     TagToken,
     ValueToken,
     ValueOrKeyToken,
     MAXTokenType
     };
   struct Token
     {
     TokenType type;
     string as_text;
     double  as_number() const
       { return strtod(as_text.c_str(), 0); }
     unsigned long  as_unsigned() const
       { return strtoul(as_text.c_str(), 0, 10); }
     };
   CIFLexer(std::istream * in)
   :input(in)
     {
     last_char = 0;
     next_char = input->get();
     }
   bool next_token(CIFLexer::Token & token);
   static CIFTagID::CIFDataName lookup_tag(const string & tag_name);
   static CIFTagID::CIFCatName lookup_cat(CIFTagID::CIFDataName tagid);
   void advance()
     {
     last_char = next_char;
     next_char = input->get();
     }
   void backup(size_t count)
     {
     for ( ++ count; count; -- count )
       input->unget();
     last_char = 0;
     next_char = input->get();
     }
   void backup(size_t count, char next)
     {
     for ( ; count; -- count )
       input->unget();
     last_char = 0;
     next_char = next;
     }
   bool good() const
     { return input->good(); }
 private:
   istream  * input;
   int  last_char, next_char;
 };
 CIFTagID::CIFDataName CIFLexer::lookup_tag(const string & tag_name)
 {
    if (CIFtagLookupTable.empty())
      {
      for (size_t idx = 0; CIFTagsRead[idx].tagid != CIFTagID::unread_CIFDataName; ++ idx)
        {
        CIFtagLookupTable.insert(CIFtagmap::value_type(string(CIFTagsRead[idx].tagname), CIFTagsRead[idx].tagid ));
        }
      }
   CIFTagID::CIFDataName rtn = CIFTagID::unread_CIFDataName;
    CIFtagmap::const_iterator found = CIFtagLookupTable.find(tag_name);
    if (found != CIFtagLookupTable.end())
      rtn = (* found).second;
    return rtn;
 }
 CIFTagID::CIFCatName CIFLexer::lookup_cat(CIFTagID::CIFDataName tagid)
 {
   CIFTagID::CIFCatName catid = CIFTagID::unread_CIFCatName;
   if (tagid > CIFTagID::unread_CIFDataName)
     {
     if (tagid < CIFTagID::MAX_atom_site)
       catid = CIFTagID::atom_site;
     else if (tagid < CIFTagID::MAX_cell)
       catid = CIFTagID::cell;
     else if (tagid < CIFTagID::MAX_chemical)
       catid = CIFTagID::chemical;
     else if (tagid < CIFTagID::MAX_chemical_formula)
       catid = CIFTagID::chemical_formula;
     else if (tagid < CIFTagID::MAX_symmetry)
       catid = CIFTagID::symmetry;
     else if (tagid < CIFTagID::MAX_symmetry_equiv)
       catid = CIFTagID::symmetry_equiv;
     else if (tagid < CIFTagID::MAX_space_group)
       catid = CIFTagID::space_group;
     else if (tagid < CIFTagID::MAX_atom_type)
       catid = CIFTagID::atom_type;
     }
   return catid;
 }

 bool CIFLexer::next_token(CIFLexer::Token & token)
 {
 token.type = CIFLexer::UnknownToken;
 token.as_text.clear();
 while (token.type == CIFLexer::UnknownToken && input->good())
   {
   if (next_char <= ' ')
     { // whitespace
     advance();
     }
   else
     { // i.e. not WhiteSpace
     switch(next_char)
       {
     // Comment handling
     case '#':
       do // eat comment to the end of the line
         {
         advance();
         } while (next_char != '\n' && input->good());
       // We are now pointing at EOL or EOF
       break;
     // Tag handling
     case '_':
       do // read name to the next whitespace
         {
         if (next_char == '.') // combines DDL1 and DDL2 tag names
           next_char = '_';
         else
           next_char = tolower(next_char);
         token.as_text.push_back((char)next_char);
         advance();
         } while (next_char > ' ' && input->good());
       // We are now pointing at the next whitespace
       token.type = CIFLexer::TagToken;
       break;
     // Quoted data handling
     case '"':
       do // read name to the next quote-whitespace
         {
         advance();
         if (next_char == '"')
           {
           while (next_char == '"')
             {
             advance();
             if (next_char <= ' ') // whitespace
               break;
             token.as_text.push_back((char)last_char);
             }
           if (next_char <= ' ') // whitespace
             break;
           }
         token.as_text.push_back((char)next_char);
         } while (input->good());
       // We are now pointing at the next whitespace
       token.type = CIFLexer::ValueToken;
       break;
     case '\'':
       do // read name to the next quote-whitespace
         {
         advance();
         if (next_char == '\'')
           {
           while (next_char == '\'')
             {
             advance();
             if (next_char <= ' ') // whitespace
               break;
             token.as_text.push_back((char)last_char);
             }
           if (next_char <= ' ') // whitespace
             break;
           }
         token.as_text.push_back((char)next_char);
         } while (input->good());
       // We are now pointing at the next whitespace
       token.type = CIFLexer::ValueToken;
       break;
     case ';':
       if (last_char == '\n')
         {
         do // read name to the next <eol>-;
           {
           advance();
           if (next_char == '\n')
             {
             while (next_char == '\n')
               {
               advance();
               if (next_char == ';') // end
                 break;
               token.as_text.push_back((char)last_char);
               }
             if (next_char == ';') // end
               {
               advance(); // go past the end
               break;
               }
             }
           token.as_text.push_back((char)next_char);
           } while (input->good());
         // We are now pointing at the next whitespace
         token.type = CIFLexer::ValueToken;
         break;
         }
       // drop through to the default case
     default: // reading an un-quoted text string
       do // read text to the next whitespace
         {
         token.as_text.push_back((char)next_char);
         advance();
         } while (next_char > ' ' && input->good());
       token.type = CIFLexer::ValueOrKeyToken;
       // We are now pointing at the next whitespace
       break;
       }
     }
   }
 if (token.type == CIFLexer::ValueOrKeyToken)
   {
   string::size_type len = token.as_text.size();
   if (len == 1 && token.as_text[0] == '.')
     token.type = CIFLexer::ValueToken;
   else if (!strncasecmp(token.as_text.c_str(), "data_", 5))
     {
     token.type = CIFLexer::KeyDataToken;
     token.as_text.erase(0, 5);
     }
   else if (!strcasecmp(token.as_text.c_str(), "loop_"))
     token.type = CIFLexer::KeyLoopToken;
   else if (!strncasecmp(token.as_text.c_str(), "save_", 5))
     {
     if (len == 5)
       {
       token.type = CIFLexer::KeySaveEndToken;
       }
     else
       {
       token.type = CIFLexer::KeySaveToken;
       token.as_text.erase(0, 5);
       }
     }
   else if (!strcasecmp(token.as_text.c_str(), "stop_"))
     token.type = CIFLexer::KeyStopToken;
   else if (!strcasecmp(token.as_text.c_str(), "global_"))
     token.type = CIFLexer::KeyGlobalToken;
   else
     token.type = CIFLexer::ValueToken;
   }
 return token.type != CIFLexer::UnknownToken;
 }
 /////////////////////////////////////////////////////////////////
  int mmCIFFormat::SkipObjects(int n, OBConversion* pConv)
 {
   if (n == 0)
     ++ n;
   CIFLexer lexer(pConv->GetInStream());
   CIFLexer::Token token;
   while (n && lexer.good())
     {
     while ( lexer.next_token(token) && token.type != CIFLexer::KeyDataToken);
     -- n;
     }
   if (lexer.good())
     lexer.backup(5 + token.as_text.size(), 'd'); // length of "data_<name>"

   return lexer.good() ? 1 : -1;
 }
 bool mmCIFFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
 {
   OBMol* pmol = pOb->CastAndClear<OBMol>();
   if(pmol==NULL)
     return false;

   CIFLexer lexer(pConv->GetInStream());
   CIFLexer::Token token;

   typedef map<string, unsigned> CIFasymmap;
   CIFasymmap asym_map;
   string last_asym_id = "";
   unsigned next_asym_no = 0;
   bool has_residue_information = false;

   pmol->SetChainsPerceived(); // avoid perception if we are setting residues

   bool wrap_coords = pConv->IsOption("w",OBConversion::INOPTIONS);

   // move to the next data block (i.e. molecule, we hope )
   while (lexer.next_token(token) && token.type != CIFLexer::KeyDataToken);
   if (token.type == CIFLexer::KeyDataToken)
     { // we have found the next data block:
     pmol->BeginModify();
     pmol->SetTitle(token.as_text);
     bool finished = false, token_peeked = false;
     double cell_a = 1.0, cell_b = 1.0, cell_c = 1.0;
     double cell_alpha = 90.0, cell_beta = 90.0, cell_gamma = 90.0;
     CIFTagID::CIFDataName name_tag = CIFTagID::unread_CIFDataName;
     CIFTagID::CIFDataName formula_tag = CIFTagID::unread_CIFDataName;
     int use_cell = 0, use_fract = 0;
     string space_group_name("P1");
     SpaceGroup space_group;
     bool space_group_failed = false;
     std::map<string, double> atomic_charges;
     while (!finished && (token_peeked || lexer.next_token(token)))
       {
       token_peeked = false;
       switch (token.type)
         {
       case CIFLexer::KeyGlobalToken:
         // We have come to the next block:
         if (pmol->NumAtoms() > 0)
           { // Found a molecule, so finished
           finished = true;
           // move back to the start of the global block:
           lexer.backup(token.as_text.size(), 'g'); // length of "global_"
           }
         else // not yet found a molecule, so go to the next data block
           {
           while (lexer.next_token(token) && token.type != CIFLexer::KeyDataToken);
           if (token.type == CIFLexer::KeyDataToken)
             { // we have found the next data block:
             pmol->SetTitle(token.as_text);
             }
           }
         break;
       case CIFLexer::KeyDataToken:
         // We have come to the next data block:
         if (pmol->NumAtoms() > 0)
           { // Found a molecule, so finished
           finished = true;
           // move back to the start of the data block:
           lexer.backup(5 + token.as_text.size(), 'd'); // length of "data_<name>"
           }
         else // not yet found a molecule, so try again
           pmol->SetTitle(token.as_text);
         break;
       case CIFLexer::KeySaveToken:
         { // Simply eat tokens until the save_ ending token
         while (lexer.next_token(token) && token.type != CIFLexer::KeySaveEndToken);
         }
         break;
       case CIFLexer::KeyLoopToken:
         {
         CIFColumnList  columns;
         CIFTagID::CIFCatName catid = CIFTagID::unread_CIFCatName;
         while ( (token_peeked = lexer.next_token(token)) == true && token.type == CIFLexer::TagToken)
           { // Read in the tags
           CIFTagID::CIFDataName tagid = lexer.lookup_tag(token.as_text);
           columns.push_back(tagid);
           if (catid == CIFTagID::unread_CIFCatName && tagid != CIFTagID::unread_CIFDataName)
             catid = lexer.lookup_cat(tagid);
           }
         size_t column_count = columns.size();
         switch (catid)
           {
         case CIFTagID::atom_site:
           {
           int use_cartn = 0, use_residue = 0;
           use_fract = 0;
           CIFTagID::CIFDataName atom_type_tag = CIFTagID::unread_CIFDataName;
           for (CIFColumnList::const_iterator colx = columns.begin(), coly = columns.end(); colx != coly; ++ colx)
             {
             switch (* colx)
               {
             case CIFTagID::_atom_site_Cartn_x:
             case CIFTagID::_atom_site_Cartn_y:
             case CIFTagID::_atom_site_Cartn_z:
               ++ use_cartn;
               break;
             case CIFTagID::_atom_site_fract_x:
             case CIFTagID::_atom_site_fract_y:
             case CIFTagID::_atom_site_fract_z:
               ++ use_fract;
               break;
             case CIFTagID::_atom_site_label_comp_id:
             case CIFTagID::_atom_site_label_seq_id:
               ++ use_residue;
               break;
             case CIFTagID::_atom_site_type_symbol:
             case CIFTagID::_atom_site_label_atom_id:
             case CIFTagID::_atom_site_label:
               if (atom_type_tag < (* colx))
                 atom_type_tag = (* colx);
               break;
             default:
               break;
               }
             }
           if (use_cartn)
             {
             for (CIFColumnList::iterator colx = columns.begin(), coly = columns.end(); colx != coly; ++ colx)
               if ( (* colx) >= CIFTagID::_atom_site_fract_x && (* colx) <= CIFTagID::_atom_site_fract_z)
                 (* colx) = CIFTagID::unread_CIFDataName;
             use_fract = 0;
             }
           size_t column_idx = 0;
           OBAtom * atom = 0;
           double x = 0.0, y = 0.0, z = 0.0;
           CIFResidueMap ResidueMap;
           unsigned long chain_num = 1, residue_num = 1;
           unsigned int nbc=0;
           string residue_name, atom_label, atom_mol_label, tmpSymbol;
           int atomicNum;
           OBPairData *label;
           OBPairFloatingPoint * occup;
           double occupancy = 1.0;
           while (token.type == CIFLexer::ValueToken) // Read in the Fields
             {
             if (column_idx == 0)
               {
               atom  = pmol->NewAtom();
               x = y = z = 0.0;
               }
             switch (columns[column_idx])
               {
             case CIFTagID::_atom_site_label: // The atomic label within the molecule
               label = new OBPairData;
               label->SetAttribute("_atom_site_label");
               label->SetValue(token.as_text);
               label->SetOrigin(fileformatInput);
               atom->SetData(label);
               atom_mol_label.assign(token.as_text);

               if (atom_type_tag != CIFTagID::_atom_site_label)
                 break;
               // Else remove everything starting from the first digit
               // and drop through to _atom_site_type_symbol
               if(string::npos != token.as_text.find_first_of("0123456789"))
                 {token.as_text.erase(token.as_text.find_first_of("0123456789"), token.as_text.size());}
             case CIFTagID::_atom_site_type_symbol:
               // Problem: posat->mSymbol is not guaranteed to actually be a
               // symbol see http://www.iucr.org/iucr-top/cif/cifdic_html/1/cif_core.dic/Iatom_type_symbol.html
               // Try to strip the string to have a better chance to have a
               // valid symbol
               // This is not guaranteed to work still, as the CIF standard
               // allows about any string...
               tmpSymbol=token.as_text.c_str();
               if ((tmpSymbol.size()==1) && isalpha(tmpSymbol[0]))
                 {
                 nbc=1;
                 }
               else if (tmpSymbol.size()>=2)
                 {
                 if (isalpha(tmpSymbol[0]) && isalpha(tmpSymbol[1]))
                   {
                   nbc=2;
                   }
                 else if (isalpha(tmpSymbol[0]))
                   {
                   nbc=1;
                   }
                 }
               else
                 {
                 nbc = 0;
                 }
               if (tmpSymbol.size()>nbc)
                 {// Try to find a formal charge in the symbol
                 int charge=0;
                 int sign=0;
                 for(unsigned int i=nbc;i<tmpSymbol.size();++i)
                   {// Use first number found as formal charge
                   if (isdigit(tmpSymbol[i]) && (charge==0))
                     {
                     charge=atoi(tmpSymbol.substr(i,1).c_str());
                     }
                   if ('-'==tmpSymbol[i])
                     {
                     sign-=1;
                     }
                   if ('+'==tmpSymbol[i])
                     {
                     sign+=1;
                     }
                   }
                   if (0!=sign) // no sign, no charge
                     {
                     if (charge==0)
                       {
                       charge=1;
                       }
                     stringstream ss;
                     ss<< tmpSymbol <<" / symbol="<<tmpSymbol.substr(0,nbc)
                       <<" charge= "<<sign*charge;
                     obErrorLog.ThrowError(__FUNCTION__, ss.str(), obDebug);
                     atom->SetFormalCharge(sign*charge);
                     }
                 }
               if (nbc>0)
                 {
                 tmpSymbol=tmpSymbol.substr(0,nbc);
                 }
               else
                 {
                 stringstream ss;
                 ss<< tmpSymbol <<" / could not derive a symbol"
                   <<" for atomic number. Setting it to default "
                   <<" Xx(atomic number 0)";
                 obErrorLog.ThrowError(__FUNCTION__, ss.str(), obDebug);
                 tmpSymbol="Xx";//Something went wrong, no symbol ! Default to Xx
                 }
               atomicNum = OBElements::GetAtomicNum(tmpSymbol.c_str());
               // Test for some oxygens with subscripts
               if (atomicNum == 0 && tmpSymbol[0] == 'O')
                 {
                 atomicNum = 8; // e.g. Ob, OH, etc.
                 }

               atom->SetAtomicNum(atomicNum); //set atomic number, or '0' if the atom type is not recognized
               atom->SetType(tmpSymbol);
               break;
             case CIFTagID::_atom_site_fract_x:
             case CIFTagID::_atom_site_Cartn_x:
               x = token.as_number();
               break;
             case CIFTagID::_atom_site_fract_y:
             case CIFTagID::_atom_site_Cartn_y:
               y = token.as_number();
               break;
             case CIFTagID::_atom_site_fract_z:
             case CIFTagID::_atom_site_Cartn_z:
               z = token.as_number();
               break;
             case CIFTagID::_atom_site_label_atom_id: // The atomic label within the residue
               atom_label.assign(token.as_text);
               if (atom_type_tag == CIFTagID::_atom_site_label_atom_id)
                 {
                 for (string::iterator posx = token.as_text.begin(), posy = token.as_text.end(); posx != posy; ++ posx)
                   {
                   char c = (char)toupper(* posx);
                   if ( c < 'A' || c > 'Z' )
                     {
                     token.as_text.erase(posx, posy);
                     break;
                     }
                   }
                 atom->SetAtomicNum(OBElements::GetAtomicNum(token.as_text.c_str()));
                 atom->SetType(token.as_text);
                 }
               break;
             case CIFTagID::_atom_site_label_comp_id: // The residue abbreviation, e.g. ILE
               residue_name.assign(token.as_text);
               break;
             case CIFTagID::_atom_site_label_entity_id: // The chain entity number of the residue, e.g. 2
    // ignored and replaced by unique id for label_asym_id
               break;
             case CIFTagID::_atom_site_label_asym_id: // The strand number of the residue
                   if (token.as_text != last_asym_id) {
                       CIFasymmap::const_iterator asym_it = asym_map.find(token.as_text);
                          if (asym_it == asym_map.end()) {
                              ++next_asym_no;
                              asym_it =
                                  asym_map.insert(CIFasymmap::value_type(token.as_text,
                                                                         next_asym_no)).first;
                          }
                          chain_num = asym_it->second;
                          last_asym_id = token.as_text;
               }
               break;
             case CIFTagID::_atom_site_label_seq_id: // The sequence number of the residue, within the chain, e.g. 12
               residue_num = token.as_unsigned();
               break;
             case CIFTagID::_atom_site_occupancy: // The occupancy of the site.
               occup = new OBPairFloatingPoint;
               occup->SetAttribute("_atom_site_occupancy");
               occupancy = token.as_number();
               if (occupancy <= 0.0 || occupancy > 1.0){
                 occupancy = 1.0;
               }
               occup->SetValue(occupancy);
               occup->SetOrigin(fileformatInput);
               atom->SetData(occup);
               break;
             case CIFTagID::unread_CIFDataName:
             default:
               break;
               }
             ++ column_idx;
             if (column_idx == column_count)
               {
               atom->SetVector(x, y, z);
               if (use_residue == 2)
                 {
                 has_residue_information = true;
                 CIFResidueID res_id(chain_num, residue_num);
                 CIFResidueMap::const_iterator resx = ResidueMap.find(res_id);
                 OBResidue * res;
                 if (resx == ResidueMap.end())
                   {
                   ResidueMap[res_id] = pmol->NumResidues();
                   res  = pmol->NewResidue();
                   res->SetChainNum(chain_num);
                   res->SetNum(residue_num);
                   res->SetName(residue_name);
                   }
                 else
                   res = pmol->GetResidue( (* resx).second );
                 res->AddAtom(atom);
                 if (!atom_label.empty())
                   res->SetAtomID(atom, atom_label);
                 unsigned long serial_no = strtoul(atom_mol_label.c_str(), 0, 10);
                 if (serial_no > 0)
                   res->SetSerialNum(atom, serial_no);
                 }
               column_idx = 0;
               }
             token_peeked = lexer.next_token(token);
             }
           }
           break;
         case CIFTagID::symmetry_equiv:
           {
           size_t column_idx = 0;
           while (token.type == CIFLexer::ValueToken) // Read in the Fields
             {
             if ((columns[column_idx] == CIFTagID::_symmetry_equiv_pos_as_xyz)
               && token.as_text.find(UNKNOWN_VALUE) == string::npos)
               space_group.AddTransform(token.as_text);
             ++ column_idx;
             if (column_idx == column_count)
               column_idx = 0;
             token_peeked = lexer.next_token(token);
             }
           }
           break;

         case CIFTagID::atom_type: //Atoms oxidations
           {
           size_t column_idx = 0;
           string atom_label = "";
           double charge = 0;
           while (token.type == CIFLexer::ValueToken) // Read in the Fields
             {
             if (columns[column_idx] == CIFTagID::_atom_type_symbol)
               atom_label = token.as_text;
             if (columns[column_idx] == CIFTagID::_atom_type_oxidation_number)
               charge = token.as_number();
             ++ column_idx;
             if (column_idx == column_count)
             {
               atomic_charges[atom_label] = charge;
               column_idx = 0;
             }
             token_peeked = lexer.next_token(token);
             }
           }
           break;


         case CIFTagID::unread_CIFCatName:
         default:
           while (token.type == CIFLexer::ValueToken) // Eat the values, we don't want them
             token_peeked = lexer.next_token(token);
           break;
           }
         }
         break;
       case CIFLexer::TagToken:
         {
         CIFTagID::CIFDataName tag_id = lexer.lookup_tag(token.as_text);
         // get the value
         lexer.next_token(token);
         switch (tag_id)
           {
         case CIFTagID::_cell_length_a:
           cell_a = token.as_number();
           ++ use_cell;
           break;
         case CIFTagID::_cell_length_b:
           cell_b = token.as_number();
           ++ use_cell;
           break;
         case CIFTagID::_cell_length_c:
           cell_c = token.as_number();
           ++ use_cell;
           break;
         case CIFTagID::_cell_angle_alpha:
           cell_alpha = token.as_number();
           ++ use_cell;
           break;
         case CIFTagID::_cell_angle_beta:
           cell_beta = token.as_number();
           ++ use_cell;
           break;
         case CIFTagID::_cell_angle_gamma:
           cell_gamma = token.as_number();
           ++ use_cell;
           break;
         case CIFTagID::_chemical_name_systematic:
         case CIFTagID::_chemical_name_mineral:
         case CIFTagID::_chemical_name_structure_type:
         case CIFTagID::_chemical_name_common:
           if (tag_id > name_tag)
             {
             name_tag = tag_id;
             pmol->SetTitle(token.as_text);
             }
           break;
         case CIFTagID::_chemical_formula_analytical:
         case CIFTagID::_chemical_formula_structural:
         case CIFTagID::_chemical_formula_iupac:
         case CIFTagID::_chemical_formula_moiety:
           if (tag_id > formula_tag)
             {
             formula_tag = tag_id;
             pmol->SetFormula(token.as_text);
             }
           break;
         case CIFTagID::_space_group_IT_number:
         case CIFTagID::_symmetry_Int_Tables_number:
           space_group_name.assign(token.as_text);
           space_group.SetId(atoi(space_group_name.c_str()));
           break;
         case CIFTagID::_space_group_name_Hall:
         case CIFTagID::_symmetry_space_group_name_Hall:
           space_group_name.assign(token.as_text);
           space_group.SetHallName(space_group_name.c_str());
           break;
         case CIFTagID::_space_group_name_H_M_alt:
         case CIFTagID::_symmetry_space_group_name_H_M:
           space_group_name.assign(token.as_text);
           space_group.SetHMName(space_group_name.c_str());
           break;
         case CIFTagID::_symmetry_equiv_pos_as_xyz:
           space_group.AddTransform(token.as_text);
           break;
         default: // eat the value for this tag
           break;
           }
         }
         break;
       case CIFLexer::KeyStopToken:
       case CIFLexer::ValueToken:
       default:
         break;
         }
       }
     if (pmol->NumAtoms() > 0)
       {
       if (use_cell >= 6)
         {
         OBUnitCell * pCell = new OBUnitCell;
         pCell->SetOrigin(fileformatInput);
         pCell->SetData(cell_a, cell_b, cell_c,
                        cell_alpha,
                        cell_beta,
                        cell_gamma
                        );
         pCell->SetSpaceGroup(space_group_name);
         const SpaceGroup * pSpaceGroup = SpaceGroup::Find( & space_group);
         if (pSpaceGroup)
           pCell->SetSpaceGroup(pSpaceGroup);
         else
           space_group_failed = true;
         pmol->SetData(pCell);
         if (use_fract)
           {
           for (OBAtomIterator atom_x = pmol->BeginAtoms(), atom_y = pmol->EndAtoms(); atom_x != atom_y; ++ atom_x)
             {
             OBAtom * atom = (* atom_x);
             if (wrap_coords)
               atom->SetVector(pCell->FractionalToCartesian(
                               pCell->WrapFractionalCoordinate(atom->GetVector())));
             else
               atom->SetVector(pCell->FractionalToCartesian(atom->GetVector()));
             }
           }
         }
       for (OBAtomIterator atom_x = pmol->BeginAtoms(), atom_y = pmol->EndAtoms(); atom_x != atom_y; ++atom_x )
       {
         OBAtom * atom = (* atom_x);
         OBPairData * pd = dynamic_cast<OBPairData *>( atom->GetData( "_atom_site_label" ) );
         if ( pd != NULL )
         {
           if( atomic_charges.count( pd->GetValue() ) > 0 )
           {
               OBPairFloatingPoint * charge_obd = new OBPairFloatingPoint;
               charge_obd->SetAttribute("input_charge");
               charge_obd->SetValue(atomic_charges[pd->GetValue()] );
               charge_obd->SetOrigin(fileformatInput);
               atom->SetData(charge_obd);
           }
         }
       }

       if (!pConv->IsOption("b",OBConversion::INOPTIONS))
         {
         pmol->ConnectTheDots();
         if (!pConv->IsOption("s",OBConversion::INOPTIONS))
           pmol->PerceiveBondOrders();
         }
       }

       if (space_group_failed)
       {
         string transformations;
         transform3dIterator ti;
         const transform3d *t = space_group.BeginTransform(ti);
         while(t){
           transformations += t->DescribeAsString() + " ";
           t = space_group.NextTransform(ti);
         }

         OBOp* pOp = OBOp::FindType("fillUC");
         if (pOp && transformations.length())
         {
           map<string, string> m;
           m.insert(pair<string, string>("transformations", transformations));
           pOp->Do(pmol, "strict", &m);
         }
       }

     pmol->EndModify();
     }
   if (has_residue_information)
     pmol->SetChainsPerceived();
   return (pmol->NumAtoms() > 0 ? true : false);
 }

 ////////////////////////////////////////////////////////////////

 bool mmCIFFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
 {
   OBMol* pmol = dynamic_cast<OBMol*>(pOb);
   if(pmol==NULL)
     return false;

   //Define some references so we can use the old parameter names
   ostream & ofs = * pConv->GetOutStream();

   char buffer[BUFF_SIZE];

   string id;
   for (const char * p = pmol->GetTitle(); * p; ++ p)
     if ( (* p) > ' ' && (* p) <= '~' )
       id.append(1, (char)toupper(* p));
   if (id.empty())
     {
     snprintf(buffer, BUFF_SIZE, "T%lu", (unsigned long)time(0));
     id.assign(buffer);
     }
   ofs << "# --------------------------------------------------------------------------" << endl;
   ofs << "#" << endl;
   ofs << "# CIF file generated by openbabel " << BABEL_VERSION << " http://openbabel.org/" << endl;
   ofs << "# to comply with the Macromolecular CIF Dictionary  (cif_mm.dic) version  2.0.11 http://mmcif.pdb.org/" << endl;
   ofs << "# The contents of this file were derived from " << pConv->GetInFilename() << endl;
   ofs << "#" << endl;
   ofs << "#---------------------------------------------------------------------------" << endl;
   ofs << endl;
   ofs << "data_" << id << endl;
   ofs << endl;
   ofs << "###########" << endl;
   ofs << "## ENTRY ##" << endl;
   ofs << "###########" << endl;
   ofs << endl;
   ofs << "_entry.id\t" << id << endl;
   ofs << endl;
   if (* (pmol->GetTitle()))
     {
     ofs << "##############" << endl;
     ofs << "## CHEMICAL ##" << endl;
     ofs << "##############" << endl;
     ofs << endl;
     ofs << "_chemical.entry_id\t" << id << endl;
     ofs << "_chemical.name_common\t'" << pmol->GetTitle() << "'" << endl;
     ofs << endl;
     }
   if (! (pmol->GetSpacedFormula().empty()))
     {
     ofs << "######################" << endl;
     ofs << "## CHEMICAL FORMULA ##" << endl;
     ofs << "######################" << endl;
     ofs << endl;
     ofs << "_chemical_formula.entry_id\t" << id << endl;
     ofs << "_chemical_formula.structural\t'" << pmol->GetFormula() << "'" << endl;
     ofs << endl;
     }
   ofs << "###############" << endl;
   ofs << "## ATOM_SITE ##" << endl;
   ofs << "###############" << endl;
   ofs << endl;
   ofs << "loop_" << endl;
   ofs << "_atom_site.id" << endl;
   ofs << "_atom_site.type_symbol" << endl;
   bool has_residues = (pmol->NumResidues() > 0);
   if (has_residues)
     {
     ofs << "_atom_site.label_atom_id" << endl;
     ofs << "_atom_site.label_comp_id" << endl;
     ofs << "_atom_site.label_entity_id" << endl;
     ofs << "_atom_site.label_seq_id" << endl;
     }
   ofs << "_atom_site.Cartn_x" << endl;
   ofs << "_atom_site.Cartn_y" << endl;
   ofs << "_atom_site.Cartn_z" << endl;
   size_t site_id = 1;
   for (OBAtomIterator atom_x = pmol->BeginAtoms(), atom_y = pmol->EndAtoms(); atom_x != atom_y; ++ atom_x, ++ site_id)
     {
     OBAtom * atom = (* atom_x);
     ofs << '\t' << site_id << '\t' << OBElements::GetSymbol(atom->GetAtomicNum());
     if (has_residues)
       {
       OBResidue * pRes = atom->GetResidue();
       string resname(pRes->GetName()), atomname(pRes->GetAtomID(atom));
       if (atomname.empty())
         {
         snprintf(buffer, BUFF_SIZE, "%s%lu", OBElements::GetSymbol(atom->GetAtomicNum()), (unsigned long)site_id);
         atomname.assign(buffer);
         }
       if (resname.empty())
         resname.assign("UNK");
       ofs << '\t' << atomname << '\t' << resname << '\t' << pRes->GetChainNum() << '\t' << pRes->GetNum() << endl;
       }
     ofs << '\t' << atom->GetX() << '\t' << atom->GetY() << '\t' << atom->GetZ() << endl;
     }
   ofs << endl;
   if (pmol->HasData(OBGenericDataType::UnitCell))
     {
     OBUnitCell * pCell = (OBUnitCell * )pmol->GetData(OBGenericDataType::UnitCell);
     ofs << "##########" << endl;
     ofs << "## CELL ##" << endl;
     ofs << "##########" << endl;
     ofs << endl;
     ofs << "_cell.entry_id\t" << id << endl;
     ofs << "_cell.length_a\t" << pCell->GetA() << endl;
     ofs << "_cell.length_b\t" << pCell->GetB() << endl;
     ofs << "_cell.length_c\t" << pCell->GetC() << endl;
     ofs << "_cell.angle_alpha\t" << pCell->GetAlpha() << endl;
     ofs << "_cell.angle_beta\t"  << pCell->GetBeta() << endl;
     ofs << "_cell.angle_gamma\t" << pCell->GetGamma() << endl;
     ofs << endl;
     const SpaceGroup * pSG = pCell->GetSpaceGroup();
     if (pSG)
       {
       ofs << "#################" << endl;
       ofs << "## SPACE GROUP ##" << endl;
       ofs << "#################" << endl;
       ofs << endl;
       ofs << "_space_group.id\t" << id << endl;
       bool loop_transforms = false;
       if (pSG->GetId())
         {
         ofs << "_space_group.IT_number\t" << pSG->GetId() << endl;
         }
       if (! (pSG->GetHallName().empty()))
         {
         loop_transforms = true;
         ofs << "_space_group.name_Hall\t'" << pSG->GetHallName() << "'" << endl;
         }
       if (! (pSG->GetHMName().empty()))
         {
         loop_transforms = true;
         // Do we have an extended HM symbol, with origin choice as ":1" or ":2" ? If so, remove it.
         size_t n=pSG->GetHMName().find(":");
         if(n==string::npos)
           ofs << "_space_group_name_H-M_alt '" << pSG->GetHMName() << "'" << endl;
         else
           ofs << "_space_group_name_H-M_alt '" << pSG->GetHMName().substr(0,n) << "'" << endl;
         }
       ofs << endl;

       transform3dIterator transx;
       const transform3d * trans = pSG->BeginTransform(transx);
       if (trans)
         {
         ofs << "####################" << endl;
         ofs << "## SYMMETRY EQUIV ##" << endl;
         ofs << "####################" << endl;
         ofs << endl;
         size_t symid = 1;
         if (!loop_transforms)
           {
           ofs << "_symmetry_equiv.id\t" << symid << endl;
           ofs << "_symmetry_equiv.pos_as_xyz\t'" << trans->DescribeAsString() << "'" << endl;
           }
         else
           {
           ofs << "loop_" << endl;
           ofs << "_symmetry_equiv.id" << endl;
           ofs << "_symmetry_equiv.pos_as_xyz" << endl;
           while (trans)
             {
             ofs << '\t' << symid << "\t'" << trans->DescribeAsString() << "'" << endl;
             trans = pSG->NextTransform(transx);
             }
           }
         ofs << endl;
         }
       }
     }

   return true;
 }

} //namespace OpenBabel
