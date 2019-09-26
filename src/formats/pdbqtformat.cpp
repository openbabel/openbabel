/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2003-2006 Geoffrey R. Hutchison
Some portions Copyright (C) 2004 by Chris Morley

Original Copyright refers to the pdbformat.cpp file, for reading and
writing pdb format files.
Extensively modified 2010 Stuart Armstrong (Source Science/InhibOx)
for the purpose of reading and writing pdbqt format files.
Some portions Copyright (C) 2010 by Stuart Armstrong of Source Science/
InhibOx

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
#include <openbabel/bond.h>
#include <openbabel/data.h>
#include <openbabel/obiter.h>
#include <openbabel/typer.h>

#include <algorithm>
#include <cstdlib>
#include <vector>
#include <map>
#include <set>

#include <sstream>

using namespace std;
namespace OpenBabel
{
  class branch
  {
    public:
    vector <int> atoms;
    bool done;
    unsigned int index;
    set <unsigned int> children;
    vector <unsigned int> parents;
    unsigned int depth;
    unsigned int connecting_atom_parent;
    unsigned int connecting_atom_branch;
    unsigned int how_many_atoms_moved;

    set <unsigned int> rigid_with; //the other branches that move rigidly with this one

    void clear() {done=false; index=0; depth=0; connecting_atom_parent=0; connecting_atom_branch=0;
      how_many_atoms_moved=0; children.clear(); parents.clear(); atoms.clear(); rigid_with.clear(); parents.push_back(0);}
    unsigned int UpOne() {if (parents.size()>=2) {return parents.at(parents.size()-2);} return 0;}
    branch() {clear();}
    void all_atoms(OBMol& mol) {clear(); rigid_with.insert(0); for (unsigned int i=1; i <= mol.NumAtoms(); i++) {atoms.push_back(i);}}
  };

  class PDBQTFormat : public OBMoleculeFormat
  {
    public:
    //Register this format type ID
    PDBQTFormat()
    {
      OBConversion::RegisterFormat("pdbqt",this, "chemical/x-pdbqt");
    }

    virtual const char* Description() //required
    {
      return

      "AutoDock PDBQT format\n"
      "Reads and writes AutoDock PDBQT (Protein Data Bank, Partial Charge (Q), & Atom Type (T)) format\n"
      "Note that the torsion tree is by default. Use the ``r`` write option\n"
      "to prevent this.\n\n"

      "Read Options, e.g. -ab\n"
      "  b  Disable automatic bonding\n"
      "  d  Input file is in dlg (AutoDock docking log) format\n\n"

      "Write Options, e.g. -xr\n"
      "  b  Enable automatic bonding\n"
      "  r  Output as a rigid molecule (i.e. no branches or torsion tree)\n"
      "  c  Combine separate molecular pieces of input into a single rigid molecule (requires \"r\" option or will have no effect)\n"
      "  s  Output as a flexible residue\n"
      "  p  Preserve atom indices from input file (default is to renumber atoms sequentially)\n"
      "  h  Preserve hydrogens\n"
			"  n  Preserve atom names\n\n";
    };

    virtual const char* SpecificationURL()
      {return "http://autodock.scripps.edu/faqs-help/faq/what-is-the-format-of-a-pdbqt-file";};

    virtual const char* GetMIMEType()
      {return "chemical/x-pdbqt";};

    //*** This section identical for most OBMol conversions ***
    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual int SkipObjects(int n, OBConversion* pConv);
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

  };
  //***
  //Make an instance of the format class
  PDBQTFormat thePDBQTFormat;

  ////////////////////////////////////////////////////
  /// Utility functions
  static bool parseAtomRecord(char *buffer, OBMol & mol, int chainNum);
  static bool IsRotBond_PDBQT(OBBond * the_bond);
  static bool IsIn(const vector<int>& vec, const int num);
  static void OutputAtom(OBAtom* atom, ostream& ofs, unsigned int index);
  static void OutputGroup(OBMol& mol, ostream& ofs, const vector <int>& group, map <unsigned int, unsigned int> new_indexes, bool use_new_indexes);
  static unsigned int AtomsSoFar(const map <unsigned int, branch >& tree, unsigned int depth);
  static bool FindBondedPiece(const vector<int>& root, const vector<int>& branch, unsigned int& root_atom, unsigned int& branch_atom,
                unsigned int& root_atom_rank, unsigned int& branch_atom_rank, const OBMol& mol, unsigned int & atoms_moved);
  static bool OutputTree(OBConversion *pConv, OBMol& mol, ostream& ofs, map <unsigned int, branch >& tree, unsigned int depth, bool moves_many, bool preserve_original_index);
  static void ConstructTree (map <unsigned int, branch >& tree, vector <vector <int> > rigid_fragments, unsigned int root_piece, const OBMol& mol, bool flexible);
  static bool DeleteHydrogens(OBMol & mol);
  static bool Separate_preserve_charges(OBMol & mol, vector<OBMol> & result);
  static unsigned int FindFragments(OBMol mol, vector <vector <int> >& rigid_fragments);
  static unsigned int RotBond_count(OBMol & mol);

  /////////////////////////////////////////////////////////////////
  int PDBQTFormat::SkipObjects(int n, OBConversion* pConv)
  {
    if (n == 0)
    ++ n;
    istream &ifs = *pConv->GetInStream();
    char buffer[BUFF_SIZE];
    while (n && ifs.getline(buffer,BUFF_SIZE))
    {
      if (EQn(buffer,"ENDMDL",6)) {-- n;}
    }

    return ifs.good() ? 1 : -1;
  }
  /////////////////////////////////////////////////////////////////
  bool PDBQTFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
    return false;


    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title = pConv->GetTitle();

    bool dlg=false;
    if (pConv->IsOption("d",OBConversion::INOPTIONS)) {dlg=true;} //check whether we have a file in dlg format

    int chainNum = 1;
    char buffer[BUFF_SIZE];
    string line, key, value;
    OBPairData *dp;

    mol.SetTitle(title);
    mol.SetChainsPerceived(); // It's a PDBQT file, we read all chain/res info.

    mol.BeginModify();
    while (ifs.good() && ifs.getline(buffer,BUFF_SIZE))
    {
      line = buffer;
      if (dlg) //if we have a dlg file, only care about the lines starting with "DOCKED: "
      {
        if (!EQn(buffer,"DOCKED: ",8)) {continue;}
        else
        {
          for (unsigned int i=0; i<BUFF_SIZE-8; i++)
          {
            buffer[i]=buffer[i+8];
            if (buffer[i]=='\0') {break;}
          }
        }
      }
      if (line.length() == 0)
      {
        stringstream errorMsg;
        errorMsg << "Warning: empty line, ignoring." << endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
        continue;
      }
      if (line.length() < 3)
      {
        stringstream errorMsg;
        errorMsg << "ERROR: not a valid PDBQT file" << endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obError);
        return false;
      }
      if (EQn(buffer,"ROOT",4)) {continue;}
      if (EQn(buffer,"ENDROOT",7)) {continue;}
      if (EQn(buffer,"BRANCH",6)) {continue;}
      if (EQn(buffer,"ENDBRANCH",9)) {continue;}

      if (EQn(buffer,"ENDMDL",6)) {break;}
      if (EQn(buffer,"END_RES",7)) {break;}

      if (EQn(buffer,"END",3))
      {
        // eat anything until the next ENDMDL
        while (ifs.getline(buffer,BUFF_SIZE) && !EQn(buffer,"ENDMDL",6));
        break;
      }

      if (EQn(buffer,"ATOM",4) || EQn(buffer,"HETATM",6))
      {
        if( ! parseAtomRecord(buffer,mol,chainNum))
        {
          stringstream errorMsg;
          errorMsg << "WARNING: Problems reading a PDBQT file\n"
            << "  Problems reading a ATOM/HETATM record.\n"
            << endl << buffer << endl;
          obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obError);
        }
        continue;
      }
      if ((EQn(buffer,"REMARK",6)) || (EQn(buffer,"USER",4)))
      {
        stringstream sst;
        string buffstring=buffer;
        sst.str(buffstring);
        sst >> buffstring;
        sst >> buffstring;
        if (buffstring == "Name")
        {
          sst >> buffstring;
          if (sst.good())
          {
            sst >> buffstring;
            mol.SetTitle(buffstring);
          }
        }
        else if (buffstring == "NEWDPF")
        {
          sst >> buffstring;
          if (buffstring == "move")
          {
            if (sst.good())
                                          {
                                                  sst >> buffstring;
                                                  mol.SetTitle(buffstring);
                                          }
          }
        }
      }

      if (line.length() <= 6)
      {
        continue;
      }


      key = line.substr(0,6); // the first 6 characters are the record name
      Trim(key);
      value = line.substr(6);

      // We haven't found this record yet
      if (!mol.HasData(key))
      {
        dp = new OBPairData;
        dp->SetAttribute(key);
        dp->SetValue(value);
        dp->SetOrigin(fileformatInput);
        mol.SetData(dp);
      }
      // Add on additional lines
      else
      {
        dp = static_cast<OBPairData*>(mol.GetData(key));
        line = dp->GetValue();
        line += '\n';
        line += value;
        dp->SetValue(line);
      }
    }

    if (!mol.NumAtoms())
    { // skip the rest of this processing
      mol.EndModify();
      return(false);
    }

    resdat.AssignBonds(mol);
    /*assign hetatm bonds based on distance*/

    mol.EndModify();
    // Clear all virtual bond data
    vector<OBGenericData*> vbonds = mol.GetAllData(OBGenericDataType::VirtualBondData);
    mol.DeleteData(vbonds);

    if (!pConv->IsOption("b",OBConversion::INOPTIONS)) {mol.ConnectTheDots(); mol.PerceiveBondOrders();}

    mol.SetChainsPerceived();

    // clean out remaining blank lines
    std::streampos ipos;
    do
    {
      ipos = ifs.tellg();
      ifs.getline(buffer,BUFF_SIZE);
    }
    while(strlen(buffer) == 0 && !ifs.eof() );
    ifs.seekg(ipos);

    mol.SetPartialChargesPerceived();
    return(true);
  }

  /////////////////////////////////////////////////////////////////////////
  void OutputAtom(OBAtom* atom, ostream& ofs, const unsigned int index)
  {
    char buffer[BUFF_SIZE];
    char type_name[10], padded_name[10];
    char the_res[10];
    char the_chain = ' ';
    char the_icode = ' ';
    const char *element_name;
    string element_name_string;
    int res_num;
    bool het=false;

    OBResidue *res;
    strncpy(type_name, OBElements::GetSymbol(atom->GetAtomicNum()), sizeof(type_name));
    type_name[sizeof(type_name) - 1] = '\0';
    //two char. elements are on position 13 and 14 one char. start at 14

    if (strlen(type_name) > 1)
      type_name[1] = toupper(type_name[1]);
    else
    {
      char tmp[10];
      strncpy(tmp, type_name, 9); // Make sure to null-terminate
      snprintf(type_name, sizeof(type_name), " %-3s", tmp);
    }

    if ( (res = atom->GetResidue()) != 0 )
    {
      snprintf(the_res,4,"%s",(char*)res->GetName().c_str());
      snprintf(type_name,5,"%s",(char*)res->GetAtomID(atom).c_str());
      the_chain = res->GetChain();
      the_icode = res->GetInsertionCode();
      if(the_icode == 0) the_icode = ' ';

      //two char. elements are on position 13 and 14 one char. start at 14
      if (strlen(OBElements::GetSymbol(atom->GetAtomicNum())) == 1)
      {
        if (strlen(type_name) < 4)
        {
          char tmp[10];
          strncpy(tmp, type_name, 9); // make sure to null-terminate
          snprintf(padded_name, sizeof(padded_name), " %-3s", tmp);
          strncpy(type_name,padded_name,4);
          type_name[4] = '\0';
        }
        else
        {
          type_name[4] = '\0';
        }
      }
      res_num = res->GetNum();
    }
    else
    {
      strcpy(the_res,"UNK");
      snprintf(padded_name,sizeof(padded_name), "%s",type_name);
      strncpy(type_name,padded_name,4);
      type_name[4] = '\0';
      res_num = 1;
    }

    element_name = OBElements::GetSymbol(atom->GetAtomicNum());
    char element_name_final[3];
    element_name_final[2] = '\0';

    if (atom->GetAtomicNum() == OBElements::Hydrogen) {element_name_final[0]='H'; element_name_final[1]='D';}
    else if ((atom->GetAtomicNum() == OBElements::Carbon) && (atom->IsAromatic())) {element_name_final[0]='A'; element_name_final[1]=' ';}
    else if (atom->GetAtomicNum() == OBElements::Oxygen)  {element_name_final[0]='O'; element_name_final[1]='A';}
    else if ((atom->GetAtomicNum() == OBElements::Nitrogen) && (atom->IsHbondAcceptor())) {element_name_final[0]='N'; element_name_final[1]='A';}
    else if ((atom->GetAtomicNum() == OBElements::Sulfur) && (atom->IsHbondAcceptor())) {element_name_final[0]='S'; element_name_final[1]='A';}
    else
    {
      if (!isalnum(element_name[0])) {element_name_final[0]=' ';}
      else element_name_final[0]=element_name[0];
      if (!isalnum(element_name[1])) {element_name_final[1]=' ';}
      else element_name_final[1]=element_name[1];
    }

    double charge = atom->GetPartialCharge();
    snprintf(buffer, BUFF_SIZE, "%s%5d %-4s %-3s %c%4d%c   %8.3f%8.3f%8.3f  0.00  0.00    %+5.3f %.2s",
      het?"HETATM":"ATOM  ",
      index,
      type_name,
      the_res,
      the_chain,
      res_num,
      the_icode,
      atom->GetX(),
      atom->GetY(),
      atom->GetZ(),
      charge,
      element_name_final);
    ofs << buffer;
    ofs << endl;
  }

  void OutputGroup(OBMol& mol, ostream& ofs, const vector <int>& group, map <unsigned int, unsigned int> new_indexes, bool use_new_indexes)
  {
    for (vector <int>::const_iterator it = group.begin(); it != group.end(); ++it)
    {
      if (use_new_indexes) {OutputAtom(mol.GetAtom((*it)), ofs, new_indexes.find(*it)->second);}
      else {OutputAtom(mol.GetAtom((*it)), ofs, (*it));}
    }
  }


  unsigned int FindFragments(OBMol mol, vector <vector <int> >& rigid_fragments)
  {
    unsigned int best_root_atom=1;
    unsigned int shortest_maximal_remaining_subgraph=mol.NumAtoms();
    for (unsigned int i=1; i <= mol.NumAtoms(); i++)
    //finds the root atom by copying the molecule, deleting each atom in turn, and finding the sizes of the resulting pieces
    {
      OBMol mol_pieces = mol;
      OBAtom * atom_to_del = mol_pieces.GetAtom(i);
      vector <vector <int> > frag_list;

      mol_pieces.DeleteAtom(atom_to_del, true);
      mol_pieces.ContigFragList(frag_list);
      unsigned int smrsi=0;
      for (unsigned int j = 0; j < frag_list.size(); j++)
      {
        smrsi= smrsi > frag_list.at(j).size() ? smrsi : frag_list.at(j).size();
      }
      if (smrsi < shortest_maximal_remaining_subgraph)
      {
        shortest_maximal_remaining_subgraph=smrsi;
        best_root_atom=i;
      }
    }

    vector <OBBond*> bonds_to_delete;
    OBMol mol_pieces = mol;
    for (OBBondIterator it=mol_pieces.BeginBonds(); it != mol_pieces.EndBonds(); it++)
    {
      if (IsRotBond_PDBQT((*it)))
      {
        bonds_to_delete.push_back(*it);
      }
    }
    for (vector<OBBond*>::iterator bit = bonds_to_delete.begin(); bit != bonds_to_delete.end(); ++bit)
    {
      mol_pieces.DeleteBond(*bit, true);
    }
    mol_pieces.ContigFragList(rigid_fragments);

    return best_root_atom;
  }

  bool DeleteHydrogens(OBMol & mol)
  {
    for (OBAtomIterator it=mol.BeginAtoms(); it != mol.EndAtoms(); it++)
    {
      if ( (*it)->IsNonPolarHydrogen() )
      {
        OBBondIterator voider;
        double charger = (*it)->GetPartialCharge();
        charger += ((*it)->BeginNbrAtom(voider))->GetPartialCharge();
        ((*it)->BeginNbrAtom(voider))->SetPartialCharge(charger);
      }
    }
    return mol.DeleteNonPolarHydrogens();
  }

  bool OutputTree(OBConversion* pConv, OBMol& mol, ostream& ofs, map <unsigned int, branch> & tree, unsigned int depth, bool moves_many, bool preserve_original_index)
  {
    if (tree.empty()) {return false;}
    if (depth>= tree.size()-1) {depth=tree.size()-1;}

    set <unsigned int> free_bonds; //this section is to allow the code to be generalised when using obabel rather than babel, which accepts numerical arguments as to how many bonds to fix. As laid out, it will prioritise those bonds that move the fewest atoms, unless moves_many is true, where it will prioritise those that move the most. This is moot for the moment, as either all rotatable bonds are free, or they are all rigid.

    free_bonds.insert(0);
    multimap <unsigned int, unsigned int> how_many_atoms_move;
    for (unsigned int i=1; i<tree.size(); i++)
    {
      how_many_atoms_move.insert(pair<unsigned int, unsigned int>( (*tree.find(i)).second.how_many_atoms_moved,i));
    }

    multimap <unsigned int, unsigned int>::iterator it=how_many_atoms_move.begin();
    if ((!moves_many) && !how_many_atoms_move.empty()) {
      it=how_many_atoms_move.end();
      if (it!=how_many_atoms_move.begin()) // don't move past begin
        --it;
    }
    for (unsigned int i = 1; i <= depth; i++)
    {
      free_bonds.insert((*it).second);
      if (!moves_many) {
        if (it!=how_many_atoms_move.begin())
          --it;
      }
      else{
        ++it;
      }
    }


    for (unsigned int i=tree.size()-1; i > 0; i--)
    {
      if (!free_bonds.count(i)) //adds the index of any branch with rigid rotations to its parent
      {
        unsigned int parent=(*tree.find(i)).second.UpOne();
        (*tree.find(parent)).second.rigid_with.insert(
          (*tree.find(i)).second.rigid_with.begin(),(*tree.find(i)).second.rigid_with.end() );
      }
    }

    map <unsigned int, unsigned int> new_order; //gives the new ordering of the indexes of atoms, so that they are in increasing order from 1 in the output


    if (!preserve_original_index) //generates the new ordering
    {
      unsigned int current_atom_index=1; //the index of the current atom
      for (unsigned int i=0; i < tree.size(); i++)
      {
        if (free_bonds.count(i))
        {
          for (set <unsigned int>::iterator it= (*tree.find(i)).second.rigid_with.begin() ; it != (*tree.find(i)).second.rigid_with.end(); ++it)
                                  {
            vector <int> atoms=(*tree.find(*it)).second.atoms;
            for (unsigned int j=0; j < atoms.size(); j++)
            {
              new_order.insert(pair<unsigned int, unsigned int> (atoms.at(j), current_atom_index));
              current_atom_index++;
            }
          }
        }
      }
    }

    if (!(pConv->IsOption("r",OBConversion::OUTOPTIONS)))
      ofs << "ROOT" << endl;
    for (set <unsigned int>::iterator it= (*tree.find(0)).second.rigid_with.begin() ; it != (*tree.find(0)).second.rigid_with.end(); ++it)
    {
      OutputGroup(mol, ofs, (*tree.find(*it)).second.atoms, new_order, !preserve_original_index);
    }

   if (!(pConv->IsOption("r",OBConversion::OUTOPTIONS)))
     ofs << "ENDROOT" << endl;

    for (unsigned int i=1; i < tree.size(); i++)
    {
      if (free_bonds.count(i))
      {
        ofs << "BRANCH";
        ofs.width(4);
        unsigned int parent_atom=(*tree.find(i)).second.connecting_atom_parent;
        unsigned int child_atom=(*tree.find(i)).second.connecting_atom_branch;
        if (!preserve_original_index) { ofs << (new_order.find(parent_atom))-> second;}
        else {ofs << parent_atom;}
        ofs.width(4);
        if (!preserve_original_index) {ofs << (new_order.find(child_atom))-> second;}
        else {ofs << child_atom;}
        ofs << endl;
        for (set <unsigned int>::iterator it= (*tree.find(i)).second.rigid_with.begin() ; it != (*tree.find(i)).second.rigid_with.end(); ++it)
        {
          OutputGroup(mol, ofs, (*tree.find(*it)).second.atoms, new_order, !preserve_original_index);
        }
      }
      unsigned int child=i;
      for (vector <unsigned int>::iterator it=(*tree.find(i)).second.parents.end(); it != (*tree.find(i)).second.parents.begin(); )
      {
        --it;
        if ((*it)==0) {break;} //do not close the main root; that is closed seperately
        vector <unsigned int>::iterator it_parent=it;
        --it_parent;
        if ((*tree.find(*it)).second.children.size() == 0)
        {
          if (free_bonds.count(*it))
          {
            ofs << "ENDBRANCH";
            ofs.width(4);
            unsigned int parent_atom=(*tree.find(*it)).second.connecting_atom_parent;
            unsigned int child_atom=(*tree.find(*it)).second.connecting_atom_branch;
            if (!preserve_original_index) { ofs << (new_order.find(parent_atom))-> second;}
                                    else {ofs << parent_atom;}
            ofs.width(4);
            if (!preserve_original_index) {ofs << (new_order.find(child_atom))-> second;}
                                    else {ofs << child_atom;}
            ofs << endl;
          }
          (*tree.find(*it_parent)).second.children.erase(*it);
        }
      }
    }
    return true;
  }

  void ConstructTree (map <unsigned int, branch>& tree, vector <vector <int> > rigid_fragments, unsigned int root_piece, const OBMol& mol, bool flexible)
  {
    unsigned int first_atom = 0;
    unsigned int second_atom = 0;
    unsigned int first_atom_rank = 0;
    unsigned int second_atom_rank = 0;


    branch sprog;

    sprog.atoms=rigid_fragments.at(root_piece);
    sprog.rigid_with.insert(0);

    tree.insert(pair<unsigned int, branch> (0, sprog));

    rigid_fragments.erase(rigid_fragments.begin() + root_piece);

    unsigned int position=0;
    unsigned int atoms_moved=0;
    bool fecund;
    while (!((*tree.find(0)).second.done))
    {
      fecund=!((*tree.find(position)).second.done);
      if (fecund)
      {
        bool sterile=true;
        for (unsigned int i = 0; i < rigid_fragments.size(); i++)
        {
          if (FindBondedPiece( (*tree.find(position)).second.atoms, rigid_fragments.at(i),
            first_atom, second_atom, first_atom_rank, second_atom_rank, mol, atoms_moved))
          {
            sprog.connecting_atom_parent = first_atom;
            sprog.connecting_atom_branch = second_atom;
            sprog.how_many_atoms_moved = atoms_moved;
            sprog.atoms = rigid_fragments.at(i);

            sprog.depth=(*tree.find(position)).second.depth+1;
            sprog.parents=(*tree.find(position)).second.parents; //all parents of the parent are parents too
            sprog.parents.push_back(tree.size()); //a branch is its own parent
            sprog.index=tree.size(); //the index is simply the number of precursors
            sprog.rigid_with.clear();
            sprog.rigid_with.insert(sprog.index);

            (*tree.find(position)).second.children.insert(tree.size()); //tells the current parent it has an extra child
                        tree.insert(pair<unsigned int, branch> (tree.size(), sprog)); //adds the current branch to the tree

            rigid_fragments.erase(rigid_fragments.begin() + i);
            sterile=false;
            position = tree.size()-1;
            break;
          }
        }
        if (sterile)
        {
          (*tree.find(position)).second.done=true;
        }
      }
      else {position--;}
    }
  }
  unsigned int RotBond_count(OBMol & mol)
  {
    unsigned int count=0;
    for (OBBondIterator it=mol.BeginBonds(); it!=mol.EndBonds(); it++)
    {
      if (IsRotBond_PDBQT((*it))) {count++;}
    }
    return count;
  }

  static bool IsImide(OBBond* querybond)
  {
    if (querybond->GetBondOrder() != 2)
      return(false);

    OBAtom* bgn = querybond->GetBeginAtom();
    OBAtom* end = querybond->GetEndAtom();
    if ((bgn->GetAtomicNum() == 6 && end->GetAtomicNum() == 7) ||
      (bgn->GetAtomicNum() == 7 && end->GetAtomicNum() == 6))
      return(true);

    return(false);
  }

  static bool IsAmidine(OBBond* querybond)
  {
    OBAtom *c, *n;
    c = n = NULL;

    // Look for C-N bond
    OBAtom* bgn = querybond->GetBeginAtom();
    OBAtom* end = querybond->GetEndAtom();
    if (bgn->GetAtomicNum() == 6 && end->GetAtomicNum() == 7)
    {
      c = bgn;
      n = end;
    }
    if (bgn->GetAtomicNum() == 7 && end->GetAtomicNum() == 6)
    {
      c = end;
      n =bgn;
    }
    if (!c || !n) return(false);
    if (querybond->GetBondOrder() != 1) return(false);
    if (n->GetTotalDegree() != 3) return false; // must be a degree 3 nitrogen

    // Make sure C is attached to =N
    OBBond *bond;
    vector<OBBond*>::iterator i;
    for (bond = c->BeginBond(i); bond; bond = c->NextBond(i))
    {
      if (IsImide(bond)) return(true);
    }

    // Return
    return(false);
  }


  /////////////////////////////////////////////////////////////////////////
  bool IsRotBond_PDBQT(OBBond * the_bond)
  //identifies a bond as rotatable if it is a single bond, not amide, not in a ring,
  //and if both atoms it connects have at least one other atom bounded to them
  {
    if ( the_bond->GetBondOrder() != 1 || the_bond->IsAromatic() || 
         the_bond->IsAmide() || IsAmidine(the_bond) || the_bond->IsInRing() )
      return false;
    if ( ((the_bond->GetBeginAtom())->GetExplicitDegree() == 1) || ((the_bond->GetEndAtom())->GetExplicitDegree() == 1) ) {return false;}
    return true;
  }

  bool IsIn(const vector<int>& vec, const int num) //checks whether a vector of int contains a specific int
  {
    for (vector<int>::const_iterator itv=vec.begin(); itv != vec.end(); ++itv)
    {
      if ((*itv) == num ) {return true;}
    }
    return false;
  }

  unsigned int AtomsSoFar(const map <unsigned int, branch> & tree, unsigned int depth)
  {
    if (depth > tree.size()) {return 0;}
    unsigned int numberAtoms=0;
//		for (unsigned int i = 0; i < depth; i++) {numberAtoms+= tree.at(i).second.first.size();}
    return numberAtoms;
  }

  bool FindBondedPiece(const vector<int>& root, const vector<int>& branched, unsigned int& root_atom, unsigned int& branch_atom,
    unsigned int& root_atom_rank, unsigned int& branch_atom_rank, const OBMol& mol, unsigned int & atoms_moved)
  {
    OBBond* the_bond;
    for (unsigned int i=0; i < root.size(); i++)
    {
      for (unsigned int j=0; j < branched.size(); j++)
      {
        the_bond=mol.GetBond(mol.GetAtom(root.at(i)), mol.GetAtom(branched.at(j)));
        if (the_bond != NULL)
        {
          root_atom=root.at(i);
          branch_atom=branched.at(j);
          root_atom_rank=i;
          branch_atom_rank=j;
          OBMol mol_copy=mol;
          the_bond=mol_copy.GetBond(mol_copy.GetAtom(root.at(i)), mol_copy.GetAtom(branched.at(j)));
          mol_copy.DeleteBond(the_bond, true);

          vector <vector <int> > two_pieces;
          mol_copy.ContigFragList(two_pieces);
          atoms_moved = two_pieces.at(1).size();
          return true;
        }
      }
    }
    return false;
  }
  /////////////////////////////////////////////////////////////////////////

  bool Separate_preserve_charges(OBMol & mol, vector <OBMol> & result)
  {
//    vector<OBMol> result;
    if( mol.NumAtoms() == 0 )
      return false; // nothing to do, but let's prevent a crash

    OBMolAtomDFSIter iter( mol, 1);
    OBMol newMol;
    newMol.SetAutomaticPartialCharge(false);
    int fragments = 0;
    while( mol.GetNextFragment( iter, newMol ) )
    {
      result.push_back( newMol );
      newMol.Clear();
      newMol.SetAutomaticPartialCharge(false);
    }

    return true;
  }
  /////////////////////////////////////////////////////////////////////////
  int CompareBondAtoms(const void *a, const void *b)
  {
    const OBAtom **da = (const OBAtom **)a;
    const OBAtom **db = (const OBAtom **)b;
    unsigned int aIdx = (*da)->GetIdx();
    unsigned int bIdx = (*db)->GetIdx();

    return ((aIdx > bIdx) - (aIdx < bIdx));
  }
  /////////////////////////////////////////////////////////////////////////
  int CompareBonds(const void *a, const void *b)
  {
    const OBAtom ***da = (const OBAtom ***)a;
    const OBAtom ***db = (const OBAtom ***)b;
    unsigned int aIdx[2] = { (*da)[0]->GetIdx(), (*da)[1]->GetIdx() };
    unsigned int bIdx[2] = { (*db)[0]->GetIdx(), (*db)[1]->GetIdx() };
    int cmp1;


    cmp1 = ((aIdx[0] > bIdx[0]) - (aIdx[0] < bIdx[0]));
    return (cmp1 ? cmp1 : ((aIdx[1] > bIdx[1]) - (aIdx[1] < bIdx[1])));
  }
  /////////////////////////////////////////////////////////////////////////
  bool PDBQTFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol & mol = *pmol;

    if(!mol.HasAromaticPerceived()) { //need aromaticity for correct atom typing
      aromtyper.AssignAromaticFlags(mol);
    }

    if (pConv->IsOption("b",OBConversion::OUTOPTIONS)) {mol.ConnectTheDots(); mol.PerceiveBondOrders();}
    vector <OBMol> all_pieces;
    if ( ((pConv->IsOption("c",OBConversion::OUTOPTIONS)!=NULL) && (pConv->IsOption("r",OBConversion::OUTOPTIONS)!=NULL))
      || (pConv->IsOption("n",OBConversion::OUTOPTIONS))
    )
    {
      mol.SetAutomaticPartialCharge(false);
      all_pieces.push_back(mol);
    }
    else
    {
      Separate_preserve_charges(mol, all_pieces);
    }

    for (unsigned int i = 0; i < all_pieces.size(); i++)
    {
      bool residue=false;
      string res_name="";
      string res_chain="";
      int res_num=1;
      if (pConv->IsOption("s",OBConversion::OUTOPTIONS))
      {
        residue=true;
        OBResidue* res=mol.GetResidue(0);
        res_name=res->GetName();
        res_name.resize(3);
        res_chain=res->GetChain();
        res_num=res->GetNum();
      }

      all_pieces.at(i).SetAutomaticPartialCharge(false);
      all_pieces.at(i).SetAromaticPerceived(); //retain aromatic flags in fragments
      if (!(pConv->IsOption("h",OBConversion::OUTOPTIONS))) {
        DeleteHydrogens(all_pieces.at(i));
      }

      int model_num = 0;
      char buffer[BUFF_SIZE];
      if (!residue)
      {
        if (!pConv->IsLast() || pConv->GetOutputIndex() > 1)
        { // More than one molecule record
          model_num = pConv->GetOutputIndex(); // MODEL 1-based index
          snprintf(buffer, BUFF_SIZE, "MODEL %8d", model_num);
          ofs << buffer << endl;
        }
        ofs << "REMARK  Name = " << mol.GetTitle(true) << endl;
//        ofs << "USER    Name = " << mol.GetTitle(true) << endl;
        if (!(pConv->IsOption("r",OBConversion::OUTOPTIONS)))
        {
          char type_name[10];
          int nRotBond=RotBond_count(mol);
          OBAtom ***rotBondTable = new OBAtom **[nRotBond];
          int rotBondId=0;
          int bondAtomNum;
          unsigned int end;
          OBResidue *res;
          for (OBBondIterator it=mol.BeginBonds(); it != mol.EndBonds(); it++)
          {
            if (IsRotBond_PDBQT((*it)))
            {
              rotBondTable[rotBondId] = new OBAtom *[2];
              rotBondTable[rotBondId][0] = (*it)->GetBeginAtom();
              rotBondTable[rotBondId][1] = (*it)->GetEndAtom();
              qsort(rotBondTable[rotBondId], 2, sizeof(OBAtom *), CompareBondAtoms);
              rotBondId++;
            }
          }
          qsort(rotBondTable, nRotBond, sizeof(OBAtom **), CompareBonds);
          ofs << "REMARK  " << nRotBond << " active torsions:" << endl;
          ofs << "REMARK  status: ('A' for Active; 'I' for Inactive)" << endl;
          for (rotBondId=0; rotBondId < nRotBond; rotBondId++)
          {
            snprintf(buffer, BUFF_SIZE, "REMARK  %3d  A    between atoms: ", rotBondId + 1);
            ofs << buffer;
            for (bondAtomNum=0; bondAtomNum < 2; bondAtomNum++)
            {
              memset(type_name, 0, sizeof(type_name));
              strncpy(type_name, OBElements::GetSymbol(rotBondTable[rotBondId][bondAtomNum]->GetAtomicNum()), sizeof(type_name));
              if (strlen(type_name) > 1)
                type_name[1] = toupper(type_name[1]);
              if ((res = rotBondTable[rotBondId][bondAtomNum]->GetResidue()) != 0)
              {
                snprintf(type_name,5,"%s",(char*)res->GetAtomID(rotBondTable[rotBondId][bondAtomNum]).c_str());
                // AtomIDs may start with space if read from a PDB file (rather than perceived)
                end = isspace(type_name[0]) ? 1 : 0;
                // Use sizeof() - 1 to ensure there's room for the NULL termination!
                while (end < sizeof(type_name) - 1 && type_name[end] && !isspace(type_name[end]))
                  end++;
                type_name[end] = '\0';
              }
              snprintf(buffer, BUFF_SIZE, "%s_%d", type_name + (isspace(type_name[0]) ? 1 : 0),
                rotBondTable[rotBondId][bondAtomNum]->GetIdx());
              ofs << buffer;
              if (bondAtomNum == 0)
                ofs << "  and  ";
            }
            delete [] rotBondTable[rotBondId];
            ofs << endl;
          }
          delete [] rotBondTable;
        }
        ofs << "REMARK                            x       y       z     vdW  Elec       q    Type" << endl;
//        ofs << "USER                              x       y       z     vdW  Elec       q    Type" << endl;
        ofs << "REMARK                         _______ _______ _______ _____ _____    ______ ____" << endl;
//        ofs << "USER                           _______ _______ _______ _____ _____    ______ ____" << endl;
      }
      else
      {
        ofs << "BEGIN_RES" << " " << res_name << " " << res_chain << " ";
        ofs.width(3);
        ofs << right << res_num << endl;
      }

      // before we write any records, we should check to see if any coord < -1000
      // which will cause errors in the formatting
      double minX, minY, minZ;
      minX = minY = minZ = -999.0f;
      FOR_ATOMS_OF_MOL(a, all_pieces.at(i))
      {
        if (a->GetX() < minX)
          minX = a->GetX();
        if (a->GetY() < minY)
          minY = a->GetY();
        if (a->GetZ() < minZ)
          minZ = a->GetZ();
      }
      vector3 transV = VZero;
      if (minX < -999.0)
        transV.SetX(-1.0*minX - 900.0);
      if (minY < -999.0)
        transV.SetY(-1.0*minY - 900.0);
      if (minZ < -999.0)
        transV.SetZ(-1.0*minZ - 900.0);

      // if minX, minY, or minZ was never changed, shift will be 0.0f
      // otherwise, move enough so that smallest coord is > -999.0f
      all_pieces.at(i).Translate(transV);

      bool flexible=!pConv->IsOption("r",OBConversion::OUTOPTIONS);

      vector <vector <int> > rigid_fragments; //the vector of all the rigid molecule fragments, using atom indexes
      unsigned int best_root_atom=1;
      map <unsigned int, branch> tree;
      unsigned int torsdof=0;
      unsigned int root_piece=0;
      unsigned int rotatable_bonds=0;

      if (flexible)
      {

        best_root_atom=FindFragments(all_pieces.at(i), rigid_fragments); //the return value is the root atom index

        if (residue) {best_root_atom=1;} //if this is a residue, uses the first atom as the root

        torsdof=rigid_fragments.size()-1;

        for (unsigned int j = 0; j < rigid_fragments.size(); j++)
        {
          if (IsIn((rigid_fragments.at(j)), best_root_atom)) {root_piece=j; break;} //this is the root rigid molecule fragment
        }
        ConstructTree(tree, rigid_fragments, root_piece, all_pieces.at(i), true);
        rotatable_bonds=torsdof;
      }
      else //if no rotatable bonds are selected, then won't construct a tree, instead get a whole branch directly from the OBMol
      {
        branch all_molecule_branch;
        all_molecule_branch.all_atoms(all_pieces.at(i));
        tree.insert(pair<unsigned int, branch> (0, all_molecule_branch));
        torsdof=RotBond_count(all_pieces.at(i));
      }

      bool preserve_original_index = (pConv->IsOption("p",OBConversion::OUTOPTIONS));
      if (!flexible) {preserve_original_index=false;} //no need to relabel if we are preserving the original order anyway

      if (!OutputTree(pConv, all_pieces.at(i), ofs, tree, rotatable_bonds, false, preserve_original_index) )
      {
        stringstream errorMsg;
        errorMsg << "WARNING: Problems writing a PDBQT file\n"
          <<  "  The torsion tree is wrong.\n";
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
        return false;
      }


      if (!residue)
      {
        if (!(pConv->IsOption("r",OBConversion::OUTOPTIONS)))
          ofs << "TORSDOF " << torsdof << endl;
        else
          ofs << "TER " << endl;
//        ofs << "TER" << endl;
        if (model_num)
        {
          ofs << "ENDMDL" << endl;
        }
      }
      else
      {
        ofs << "END_RES" << " " << res_name << " " << res_chain << " ";
        ofs.width(3);
        ofs << right << res_num << endl;
      }
    }
    return true;
  }

  /////////////////////////////////////////////////////////////////////////

  static bool parseAtomRecord(char *buffer, OBMol &mol,int chainNum)
  /* ATOMFORMAT "(i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3,2f6.2,a2,a2)" */
  {
    string sbuf = &buffer[6];
    if (sbuf.size() < 48)
      return(false);

    bool hetatm = (EQn(buffer,"HETATM",6)) ? true : false;
    bool elementFound = false; // true if correct element found in col 77-78

    /* serial number */
    string serno = sbuf.substr(0,5);

    /* atom name */
    string atmid = sbuf.substr(6,4);

    /* chain */
    char chain = sbuf.substr(15,1)[0];

    /* element */
    string element = "  ";
    string pdbqt_element = "  "; //the literal pdbqt element
    if (sbuf.size() == 72) {sbuf+=" ";}
    if (sbuf.size() > 72)
    {
      pdbqt_element = sbuf.substr(71,2);
      if ( (pdbqt_element == "A ") || (pdbqt_element == " A") || (pdbqt_element == "Z ") || (pdbqt_element == " Z") ||
        (pdbqt_element == "G ") || (pdbqt_element == " G") || (pdbqt_element == "GA") || (pdbqt_element == "J ") ||
        (pdbqt_element == " J") || (pdbqt_element == "Q ") || (pdbqt_element == " Q") )
        {element = "C"; element += " ";} //all these are carbons
      else if (pdbqt_element == "HD") {element = "H"; element += " ";} //HD are hydrogens
      else if (pdbqt_element == "HS") {element = "H"; element += " ";} //HS are hydrogens
      else if (pdbqt_element == "NA") {element = "N"; element += " ";} //NA are nitrogens
      else if (pdbqt_element == "NS") {element = "N"; element += " ";} //NS are nitrogens
      else if (pdbqt_element == "OA") {element = "O"; element += " ";} //OA are oxygens
      else if (pdbqt_element == "OS") {element = "O"; element += " ";} //OS are oxygens
      else if (pdbqt_element == "SA") {element = "S"; element += " ";} //SA are sulphurs
      else {element = pdbqt_element;}

      if (isalpha(element[0])) //trims the name if needed
      {
        if (element[1] == ' ')
        {
          element.erase(1, 1);
          elementFound = true;
        }
        else if (isalpha(element[1]))
        {
          elementFound = true;
        }
      }
      else if (isalpha(element[1]))
      {
        element.erase(0,1);
        elementFound = true;
      }
    }

    if (!elementFound)
    {
      stringstream errorMsg;
      errorMsg << "WARNING: Problems reading a PDBQT file\n"
        << "  Problems reading a HETATM or ATOM record.\n"
        << "  According to the PDBQT specification,\n"
        << "  columns 78-79 should contain the element symbol of an atom.\n"
        << "  but OpenBabel found '" << element << "' (atom " << mol.NumAtoms()+1 << ")";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      return false;
    }

    // charge - optional
    string scharge;
    if (sbuf.size() > 69)
    {
      scharge = sbuf.substr(64,6);
    }

    //trim spaces on the right and left sides
    while (!atmid.empty() && atmid[0] == ' ')
      atmid = atmid.erase(0, 1);

    while (!atmid.empty() && atmid[atmid.size()-1] == ' ')
      atmid = atmid.substr(0,atmid.size()-1);

    /* residue name */
    string resname = sbuf.substr(11,3);
    if (resname == "   ") resname = "UNK";
    else
    {
      while (!resname.empty() && resname[0] == ' ')
        resname = resname.substr(1,resname.size()-1);
      while (!resname.empty() && resname[resname.size()-1] == ' ')
        resname = resname.substr(0,resname.size()-1);
    }

    OBAtom atom;
    /* X, Y, Z */
    string xstr = sbuf.substr(24,8);
    string ystr = sbuf.substr(32,8);
    string zstr = sbuf.substr(40,8);
    vector3 v(atof(xstr.c_str()),atof(ystr.c_str()),atof(zstr.c_str()));
    atom.SetVector(v);

    // useful for debugging unknown atom types (e.g., PR#1577238)
    //    cout << mol.NumAtoms() + 1  << " : '" << element << "'" << " " << OBElements::GetAtomicNum(element.c_str()) << endl;


    atom.SetAtomicNum(OBElements::GetAtomicNum(element.c_str()));

    if ( (! scharge.empty()) && "     " != scharge )
    {
      stringstream sst;
      sst.str(scharge);
      double charge;
      sst >> charge;
      if ( !sst.fail() )
      {
        atom.SetPartialCharge(charge);
      }
      else
      {
        stringstream errorMsg;
        errorMsg << "WARNING: Problems reading a PDBQT file\n"
          << "  Problems reading a HETATM or ATOM record.\n"
          << "  According to the PDBQT specification,\n"
          << "  columns 72-76 should contain charge of the atom\n"
          << "  but OpenBabel found '" << scharge << "' (atom " << mol.NumAtoms()+1 << ").";
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      }
    }
    else
    {
      atom.SetPartialCharge(0.);
    }

    /* residue sequence number */
    string resnum = sbuf.substr(16,4);
    char icode = sbuf[20];
    if(icode == ' ') icode = 0;

    OBResidue *res  = (mol.NumResidues() > 0) ? mol.GetResidue(mol.NumResidues()-1) : NULL;
    if (res == NULL || res->GetName() != resname
      || res->GetNumString() != resnum || res->GetInsertionCode() != icode)
    {
      vector<OBResidue*>::iterator ri;
      for (res = mol.BeginResidue(ri) ; res ; res = mol.NextResidue(ri))
      if (res->GetName() == resname
        && res->GetNumString() == resnum
        && res->GetInsertionCode() == icode
        && static_cast<int>(res->GetChain()) == chain)
        break;

      if (res == NULL)
      {
        res = mol.NewResidue();
        res->SetChain(chain);
        res->SetName(resname);
        res->SetNum(resnum);
        res->SetInsertionCode(icode);
      }
    }

    if (!mol.AddAtom(atom))
      return(false);
    else
    {
      OBAtom *atom = mol.GetAtom(mol.NumAtoms());

      res->AddAtom(atom);
      res->SetSerialNum(atom, atoi(serno.c_str()));
      res->SetAtomID(atom, sbuf.substr(6,4));
      res->SetHetAtom(atom, hetatm);
      return(true);
    }
  } // end reading atom records

} //namespace OpenBabel
