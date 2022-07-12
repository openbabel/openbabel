/**********************************************************************
Contribution by the University of Southampton

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

#include <vector>
#include <iostream>
#include <cmath>
#include <stack>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include <openbabel/ring.h>
#include <iterator>

#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>

////////////////////////////////////////////////////////////////////////
// Utility Classes

class BondContainer{
public:
    unsigned int atom_1_index{};
    unsigned int atom_2_index{};

    unsigned int atom_1_degree{};
    unsigned int atom_2_degree{};

    unsigned int atom_1_label{};
    unsigned int atom_2_label{};

    unsigned int bond_order{};
    bool aromatic{};
    bool ring_bond{};
    unsigned int bonded_connections{};

    unsigned int state{0};
    unsigned int stack_index{};
public:
    BondContainer();
    BondContainer(  unsigned int stack_pos,
                    unsigned int start, unsigned int end,
                    unsigned int deg_1, unsigned deg_2,
                    unsigned int type_1,unsigned int type_2,
                    unsigned int order_val,
                    bool is_aromatic,
                    bool is_ring);

    BondContainer(const BondContainer &source);

    void Display();
};

BondContainer::BondContainer() {
    atom_1_index = atom_2_index = NAN;
    atom_1_degree = atom_2_degree = NAN;
    atom_1_label = atom_2_label = NAN;
    bond_order = NAN;
    aromatic = {false};
    state = 0;
    stack_index = NAN;
    bonded_connections = 0;
}

BondContainer::BondContainer(   unsigned int stack_pos,
                                unsigned int start, unsigned int end,
                                unsigned int deg_1, unsigned int deg_2,
                                unsigned int type_1, unsigned int type_2,
                                unsigned int order_val,
                                bool is_aromatic, bool is_ring)
        :   stack_index{stack_pos},
            atom_1_index{start}, atom_2_index{end},
            atom_1_degree{deg_1}, atom_2_degree{deg_2},
            atom_1_label{type_1}, atom_2_label{type_2},
            bond_order{order_val},
            aromatic{is_aromatic},
            ring_bond{is_ring}{state = 1;};


BondContainer::BondContainer(const BondContainer &source)
        :   stack_index{source.stack_index},
            atom_1_index{source.atom_1_index}, atom_2_index{source.atom_2_index},
            atom_1_degree{source.atom_1_degree}, atom_2_degree{source.atom_2_degree},
            atom_1_label{source.atom_1_label}, atom_2_label{source.atom_2_label},
            bond_order{source.bond_order},
            aromatic{source.aromatic},bonded_connections{source.bonded_connections},
            ring_bond{source.ring_bond}{state=1;};


void BondContainer::Display() {
    std::cout << "Bond: Atoms(" << atom_1_index << "," << atom_2_index <<")";
    std::cout << " Degrees(" << atom_1_degree << ',' << atom_2_degree << ')';
    std::cout << " Connections(" << bonded_connections <<")";
    std::cout << " Type(" << atom_1_label << ',' << atom_2_label <<")";
    std::cout << " Order(" << bond_order << ")";
    std::cout << " Aromatic(" << aromatic << ")";
    std::cout << " Ring(" << ring_bond << ")";
    std::cout << " Stack Index(" << stack_index << ")";
    std::cout << std::endl;
}



class AtomContainer{
public:
    unsigned int atom_index;
    unsigned int atom_degree;
    unsigned int atom_label;
    bool isinring;
    std::vector<BondContainer> connected_bonds{};
    std::vector<unsigned int> inbound_atoms{};
    std::vector<unsigned int> outbound_atoms{};
    unsigned int atom_state =0;

public:
    AtomContainer();
    AtomContainer(unsigned int atom_index_val, unsigned int atom_degree_val, unsigned int atom_label_val,
                  std::vector<BondContainer> connected_bonds_val, bool isinring_val);
    AtomContainer(const AtomContainer &source);
    void Display();
};

// Constructors etc.
AtomContainer::AtomContainer()
        :atom_index{0},atom_degree{0},atom_label{0},atom_state{0}{};

AtomContainer::AtomContainer(unsigned int atom_index_val, unsigned int atom_degree_val, unsigned int atom_label_val,
                             std::vector<BondContainer> connected_bonds_val, bool isinring_val)
        :atom_index{atom_index_val},atom_degree{atom_degree_val},atom_label{atom_label_val},connected_bonds{std::move(connected_bonds_val)}
        ,atom_state{1}, isinring(isinring_val){};

AtomContainer::AtomContainer(const AtomContainer &source)
        :atom_index{source.atom_index},atom_degree{source.atom_degree},atom_label{source.atom_label},
         connected_bonds{source.connected_bonds}, inbound_atoms{source.inbound_atoms}, outbound_atoms{source.outbound_atoms},
         atom_state{source.atom_state}, isinring{source.isinring}{};

void AtomContainer::Display() {
    std::cout << "Atom: Index(" << atom_index << ")    ";
    std::cout << "Label(" << atom_label << ")    ";
    std::cout << "Degree(" << atom_degree << ")    ";
    std::cout << "InRing(" << isinring << ")    ";

    if (!connected_bonds.empty()){
        std::cout << "   Connected Bonds: (";
        std::cout << connected_bonds.size();
        std::cout << ")";
    };

    std::cout << "   Inbounds: (";
    std::copy(inbound_atoms.begin(), inbound_atoms.end(), std::ostream_iterator<unsigned int>(std::cout, ","));
    std::cout << ")";


    std::cout << "   Outbounds: (";
    std::copy(outbound_atoms.begin(), outbound_atoms.end(), std::ostream_iterator<unsigned int>(std::cout, ","));
    std::cout << ")";

    std::cout << "\n";
}


class RingContainer{
public:
    std::vector<AtomContainer> ring_atoms{};
    std::vector<BondContainer> ring_bonds{};
    std::vector<std::tuple< AtomContainer,unsigned int>> active_atoms{};
    AtomContainer start_atom;
    AtomContainer end_atom;
    unsigned int stack_val;
    char hetero = 'L';
    unsigned int size;
    bool aromatic{};
    std::string type;

public:
    RingContainer();
    RingContainer(const RingContainer &source);

    void Display();
};


RingContainer::RingContainer()
        :ring_atoms{}{};


RingContainer::RingContainer(const RingContainer &source)
        :ring_atoms{source.ring_atoms},
         ring_bonds{source.ring_bonds},
         active_atoms{source.active_atoms},
         start_atom{source.start_atom},
         end_atom{source.end_atom},
         hetero{source.hetero},
         size{source.size},
         aromatic{source.aromatic},
         type{source.type},
         stack_val{source.stack_val}{};


void RingContainer::Display() {
    std::cout << stack_val << " ";
    std::cout << "Ring: Size(" << size << ")   ";
    std::cout << "Type(" << type << ")   ";
    std::cout << "Aromatic(" << aromatic << ")    ";
    std::cout << "Atom Start(" << start_atom.atom_index << ")   ";
    std::cout << "Atom End(" << end_atom.atom_index << ")   ";
    std::cout << "Hetero: " << hetero << "   ";
    std::cout << "Active Atoms[" << active_atoms.size() << "]   ";
    std::cout << "\n";
}



class Graph{
public:

    int size;
    std::vector<std::vector<int>> graph_2d;
    std::vector<int> contributions;
    bool HasCycle = false;
    std::vector<std::vector<int>> indexes;


    bool CreateGraph();
    void DisplayGraph();

    bool AddEdge(int src, int trg);
    bool DestroyEdge(int src, int trg);


    bool IsSafe(int n,std::vector<int> &path, int pos);
    bool CycleCheck(std::vector<int>&path, std::vector<bool> visited ,int pos);
    bool GeneratePath();
    bool ReorderPaths(int start_point);

    Graph(int size_val);
};

Graph::Graph(int size_val)
        :size{size_val}
{
    std::vector<int>blank (size, 0);
    for (int i=0; i<size;i++)
        graph_2d.push_back(blank);

}

void Graph::DisplayGraph() {

    for (int i=0; i<size; i++){
        std::vector<int> row = graph_2d.at(i);
        std::cout << "[ ";
        for (auto term: row){
            std::cout << term << " , ";

        }
        std::cout << "]\n";
    }
}
bool Graph::AddEdge(int src, int trg) {
    if (src > size || trg > size)
        return false;

    graph_2d.at(src).at(trg) = 1;
    graph_2d.at(trg).at(src) = 1;
    return true;
}
bool Graph::DestroyEdge(int src, int trg) {

    if (graph_2d.at(src).at(trg) == 1)
        graph_2d.at(src).at(trg) = 0;
    else
        return false;

    if (graph_2d.at(trg).at(src) == 1)
        graph_2d.at(trg).at(src) = 0;
    else
        return false;

    return true;
}


// Hamiltonian Cycle Checks,
bool Graph::GeneratePath() {

    HasCycle = false;
    int pos = 1;

    std::vector<int> path;
    path.push_back(0);

    std::vector<bool> visited(size,false);
    visited.at(0) = true;

    if (!CycleCheck(path, visited,pos)){
        return false;
    }
    return true;
}
bool Graph::IsSafe(int n,std::vector<int> &path, int pos){
    if( graph_2d.at(path.at(pos-1)).at(n) == 0)
        return false;
    for (int i=0; i<pos;i++){
        if (path.at(i)==n)
            return false;
    }
    return true;
}
bool Graph::CycleCheck(std::vector<int>&path, std::vector<bool> visited ,int pos){
    if (pos == size){
        if(graph_2d.at(path.at(path.size()-1)).at(path.at(0)) !=0){
            //path.push_back(0);
            //path.pop_back();
            HasCycle = true;
            indexes.push_back(path);
        }

        if(graph_2d.at(path.at(path.size()-1)).at(path.at(0)) ==0){
            indexes.push_back(path);
            HasCycle = false;
        }

        return HasCycle;
    }

    for (int n=1; n < size;n++){ // +1 indexed as we are starting at zero <-- "start point" etc..
        if (IsSafe(n, path, pos) && !visited.at(n)){

            path.push_back(n);
            visited.at(n) = true;

            CycleCheck(path, visited, pos+1);
            visited.at(n) = false;
            path.pop_back();
        }
    }
    return HasCycle;
}
bool Graph::ReorderPaths(int start_point){
    for (unsigned int i=0;i<indexes.size();i++){
        std::vector<int> &path = indexes.at(i);

        unsigned int shift=0;
        std::vector<int> stored;
        for (unsigned int k=0; k<path.size();k++){
            if (path.at(k) == start_point){
                shift = k;
                break;
            }
            stored.push_back(path.at(k));
        }
        path.erase(path.begin(), path.begin() + shift);
        path.insert(path.end(),stored.begin(),stored.end());
    }
    return true;
}

////////////////////////////////////////////////////////////////////////
// WLN Symbol Class

class WLNSymbol{
public:
    char symbol{};
    std::vector<char> multi_symbol{};
    bool terminator{};
    unsigned int serve_ring{0};
    unsigned int wln_valence{};
    int start_atom{-1};
    int end_atom{-1};
    unsigned int pos{};
    bool active{false};

    std::vector<BondContainer> inbound_connections{};
    std::vector<BondContainer> outbound_connections{};
    std::vector<AtomContainer> atom_container{};
    std::vector<unsigned char> poly_chars{};
    std::vector<unsigned char> peri_chars{};

    std::vector<AtomContainer> peri_atoms{};
    std::vector<AtomContainer> poly_atoms{};

    RingContainer contained_ring{};
    std::vector<RingContainer> SSRS{};
    std::string ring_type{"NULL"};

public:
    WLNSymbol();
    WLNSymbol(char mode_val);
    WLNSymbol(std::vector<AtomContainer> container_val);
    WLNSymbol(const WLNSymbol &source);

    // Operators
    void Display();
    void BuildSymbol();
};

WLNSymbol::WLNSymbol()
        :inbound_connections{},outbound_connections{},wln_valence{0}{};

WLNSymbol::WLNSymbol(char mode_val)
        :inbound_connections{},outbound_connections{},wln_valence{0}{};


WLNSymbol::WLNSymbol(std::vector<AtomContainer> container_val)
        :atom_container{std::move(container_val)}{};


WLNSymbol::WLNSymbol(const WLNSymbol &source)
        :symbol{source.symbol},terminator{source.terminator},wln_valence{source.wln_valence},
         start_atom{source.start_atom}, end_atom{source.end_atom}, inbound_connections{source.inbound_connections},
         outbound_connections{source.outbound_connections}, pos{source.pos}
        ,atom_container{source.atom_container},
         serve_ring{source.serve_ring}, contained_ring{source.contained_ring},
         active{source.active}, SSRS{source.SSRS}, ring_type{source.ring_type}
        ,poly_chars{source.poly_chars},peri_chars{source.peri_chars},
         multi_symbol{source.multi_symbol},peri_atoms{source.peri_atoms},poly_atoms{source.poly_atoms} {};


void WLNSymbol::Display() {
    std::cout << "Symbol: " << symbol << "\t";
    std::cout << "Contained Atoms[" << start_atom << "->" << end_atom << "]";

    std::cout << "  InCount[" << inbound_connections.size() << "]";
    std::cout << "  OutCount[" << outbound_connections.size() << "]";

    std::cout << "   WLN Valence(" << wln_valence << ")";
    std::cout << "    Position Given(" << pos << ")";
    std::cout << "   Terminator(" << terminator << ")";
    std::cout << "   ServeRing?(" << serve_ring << ")";

    if (!ring_type.empty())
        std::cout << "   Ring Type(" << ring_type << ")";
    std::cout << std::endl;
}
void WLNSymbol::BuildSymbol(){
    // Build WLN symbol from atom_container
    symbol = (atom_container.size()) + '0';
    wln_valence = 0;
    terminator = false;
    start_atom = atom_container.front().atom_index;
    end_atom = atom_container.back().atom_index;
}

////////////////////////////////////////////////////////////////////////
// WLN Writer Class

class WriteWLN{
public:
    char* smiles;
    OpenBabel::OBMol *Mol;
    std::vector<BondContainer> global_bond_stack{};
    std::vector<AtomContainer> global_atom_stack{};
    std::vector<RingContainer> global_ring_stack{};

public:
    // Constructors/Destructors
    WriteWLN();
    WriteWLN(OpenBabel::OBMol *Mol);

    // Methods
    bool BuildGlobalStacks(OpenBabel::OBMol *Mol);
    bool SegmentAtomStack(std::vector<WLNSymbol> &wln_vector);
    bool ParseAtoms(std::vector<WLNSymbol> &wln_vector);
    bool FormSystem(std::vector<WLNSymbol> &wln_vector, std::vector<char> &wln_string);
    bool SanitiseString(std::vector<char> &wln_string);

    bool CreateString(std::vector<char> &wln_string);
    bool Error(const char* type);

    // Branch Functions
    bool ChainAtom(AtomContainer atom);
    bool Joined(AtomContainer atom_1, AtomContainer atom_2, int order);
    bool RemovePos(unsigned int symbol_id, std::vector<WLNSymbol> &wln_vector);
    std::vector<BondContainer> GetConnections(unsigned int atom_index);
    std::vector<BondContainer> ForwardConnections(unsigned int atom_index);
    std::vector<BondContainer> BackwardConnections(unsigned int atom_index);
    std::vector<std::tuple<WLNSymbol,int>> NextSymbols(WLNSymbol check_symbol, std::vector<WLNSymbol> wln_vector);
    std::vector<std::tuple<WLNSymbol,int>> PreviousSymbols(WLNSymbol check_symbol, std::vector<WLNSymbol> wln_vector);
    AtomContainer FetchAtom(unsigned int atom_index);
    void CheckOrder(unsigned int order, std::vector<char> &wln_string);
    char PeriodicCharacter(AtomContainer atom);
    std::vector<unsigned int> SymbolAtomIndexes(WLNSymbol subject);
    unsigned int ExternalBond(AtomContainer atom, std::vector<AtomContainer> atom_list);
    bool MultiRing(WLNSymbol subject, WLNSymbol target);

    void CheckDuplicates(std::vector<AtomContainer> &vec);
    void CheckDuplicates(std::vector<BondContainer> &vec);

    // Ring Functions
    char FetchRingPosition(unsigned int atom_index, RingContainer ring, char mode);
    bool LinkRings(std::vector<WLNSymbol> &wln_vector);
    bool BuildRing(WLNSymbol virtual_ring, std::vector<char> &wln_string);
    bool CreateVirtualRings(std::vector<WLNSymbol> &wln_vector);
    void PushJunction(WLNSymbol &hold_junction, std::stack<WLNSymbol> &junctions);
    void StartInline(WLNSymbol &next_symbol ,std::vector<char> &wln_string);
    bool CheckPolyBinding(WLNSymbol &subject, WLNSymbol &target, std::vector<unsigned int> &atom_intersection);
    WLNSymbol MergeRingSymbol(WLNSymbol subject, WLNSymbol target, std::vector<unsigned int> atom_intersection);
    bool CreateGraph(std::vector<AtomContainer> atom_list, Graph &ring_graph);

    template<class T>
    std::vector<T> VectorReorder(std::vector<T> vA, std::vector<int> vOrder );
    bool CloneRingData(WLNSymbol &merged_ring, WLNSymbol &subject ,WLNSymbol &target, std::string type);
    unsigned int ScorePaths(std::vector<AtomContainer> atom_list, Graph &ring_graph);
    unsigned int ScoreLists(std::vector<AtomContainer> list_1, std::vector<AtomContainer> list_2);

};
std::string TranslateToWLN(const char* smiles_string);




WriteWLN::WriteWLN()
        :smiles{"NULL"}{};

WriteWLN::WriteWLN(OpenBabel::OBMol *Mol_val)
        :Mol{Mol_val}{};


// --- Methods ---
bool WriteWLN::BuildGlobalStacks(OpenBabel::OBMol *Mol) {
    FOR_BONDS_OF_MOL(bond, Mol){
        OpenBabel::OBAtom *start = bond->GetBeginAtom();
        OpenBabel::OBAtom *end = bond->GetEndAtom();

        unsigned int start_index = start->GetIndex();
        unsigned int end_index = end->GetIndex();

        unsigned int start_degree = start->GetExplicitValence();
        unsigned int end_degree = end->GetExplicitValence();

        unsigned int start_label =  start->GetAtomicNum();
        unsigned int end_label =  end->GetAtomicNum();

        unsigned int order = bond->GetBondOrder();
        bool aromatic = bond->IsAromatic();
        bool ring_bond = bond->IsInRing();

        BondContainer bond_info = BondContainer(global_bond_stack.size(),start_index,end_index,
                                                start_degree,end_degree,
                                                start_label,end_label,
                                                order, aromatic, ring_bond);
        global_bond_stack.push_back(bond_info);
    }
    FOR_ATOMS_OF_MOL(atom, Mol){
        unsigned int atom_index = atom->GetIndex();
        bool ring_check = atom->IsInRing();
        unsigned int atom_label = atom->GetAtomicNum();
        unsigned int atom_degree = atom->GetExplicitValence();
        std::vector<BondContainer> connections = GetConnections(atom_index);
        AtomContainer atom_info = AtomContainer(atom_index,atom_degree,atom_label,connections, ring_check);

        for (unsigned int i=0;i<connections.size();i++){
            BondContainer bond = connections.at(i);
            if (bond.atom_2_index != atom_index)
                atom_info.outbound_atoms.push_back(bond.atom_2_index);
            if (bond.atom_1_index != atom_index)
                atom_info.inbound_atoms.push_back(bond.atom_1_index);
        }
        global_atom_stack.push_back(atom_info);
    }
    for (unsigned int i=0; i< global_bond_stack.size();i++){
        BondContainer bond = global_bond_stack.at(i);
        global_bond_stack.at(i).bonded_connections =  GetConnections(bond.atom_1_index).size();
    }

    unsigned int iter=0;
    FOR_RINGS_OF_MOL(ring, Mol){
        RingContainer ring_container {};
        ring_container.aromatic = ring->IsAromatic();
        std::vector<int> ring_indexes = ring->_path;
        std::for_each(ring_indexes.begin(), ring_indexes.end(), [](int &n){ n+=-1; }); // add -1 to each element. zero index the atoms, complies with ISSC OB rules
        for (unsigned int i=0; i<global_atom_stack.size();i++){
            if ( std::find(ring_indexes.begin(), ring_indexes.end(), global_atom_stack.at(i).atom_index) != ring_indexes.end() ) {
                ring_container.ring_atoms.push_back(global_atom_stack.at(i));
                if (global_atom_stack.at(i).atom_label != 6)
                    ring_container.hetero = 'T';
                if (global_atom_stack.at(i).connected_bonds.size()>2){
                    std::tuple<AtomContainer, unsigned int> active_info{global_atom_stack.at(i), global_atom_stack.at(i).connected_bonds.size()-2};
                    ring_container.active_atoms.push_back(active_info);
                }
            }
        }

        for (unsigned int i=0; i<global_bond_stack.size();i++){
            BondContainer bond = global_bond_stack.at(i);
            if ( std::find(ring_indexes.begin(), ring_indexes.end(), global_bond_stack.at(i).atom_1_index) != ring_indexes.end() ) {
                ring_container.ring_bonds.push_back(global_bond_stack.at(i));
            }
        }

        ring_container.size = ring_container.ring_atoms.size();
        ring_container.start_atom = ring_container.ring_atoms.front();
        ring_container.end_atom = ring_container.ring_atoms.back();
        ring_container.stack_val = iter;

        iter++;
        global_ring_stack.push_back(ring_container);
    }

    // Mol deletion checked
    delete Mol;
    return true;
}
bool WriteWLN::SegmentAtomStack(std::vector<WLNSymbol> &wln_vector){
    if (global_atom_stack.empty())
        return false;

    // These HAVE TO BE IN INDEX ORDER <-- lost in rings notation
    std::sort(global_atom_stack.begin(),
              global_atom_stack.end(),
              [](const AtomContainer& lhs, const AtomContainer& rhs)
              {
                  return lhs.atom_index < rhs.atom_index;
              });

    bool chain_state{false};
    AtomContainer previous;
    std::vector<AtomContainer> chain{};
    for (unsigned int i=0; i< global_atom_stack.size();i++){
        AtomContainer atom = global_atom_stack.at(i);
        //handle the null start state
        if (previous.atom_state == 0){
            previous = atom;
            if (ChainAtom(atom)){
                chain.push_back(atom);
                chain_state = true;}

            if (i==global_atom_stack.size()-1){
                WLNSymbol symbol{chain};
                symbol.BuildSymbol();
                symbol.inbound_connections = BackwardConnections(symbol.start_atom);
                symbol.outbound_connections = ForwardConnections(symbol.end_atom);
                symbol.active = true;
                wln_vector.push_back(symbol);
                chain.clear();
            }
            continue;
        }

        if (chain_state){
            if (Joined(previous, atom,1)){
                if (ChainAtom(atom)){
                    previous = atom;
                    if (i == global_atom_stack.size() -1){
                        chain_state = false;
                        chain.push_back(atom);
                        WLNSymbol symbol{chain};
                        symbol.BuildSymbol();
                        symbol.inbound_connections = BackwardConnections(symbol.start_atom);
                        symbol.outbound_connections = ForwardConnections(symbol.end_atom);
                        symbol.active = true;
                        wln_vector.push_back(symbol);
                        chain.clear();
                        continue;
                    }
                    chain.push_back(atom);
                }else{
                    chain_state = false;
                    //chain.push_back(atom);
                    WLNSymbol symbol{chain};
                    symbol.BuildSymbol();
                    symbol.inbound_connections = BackwardConnections(symbol.start_atom);
                    symbol.outbound_connections = ForwardConnections(symbol.end_atom);
                    symbol.active = true;
                    wln_vector.push_back(symbol);
                    chain.clear();
                    continue;
                }

            }else{
                chain_state = false;
                WLNSymbol symbol{chain};
                symbol.BuildSymbol();
                symbol.inbound_connections = BackwardConnections(symbol.start_atom);
                symbol.outbound_connections = ForwardConnections(symbol.end_atom);
                symbol.active = true;
                wln_vector.push_back(symbol);
                chain.clear();
            }
        }

        // This needs a slight fix... some chain is being double counted.
        if (!chain_state){
            if (ChainAtom(atom)){
                chain_state = true;
                chain.push_back(atom);
                previous = atom;
                if (i == global_atom_stack.size() -1){
                    chain_state = false;
                    WLNSymbol symbol{chain};
                    symbol.BuildSymbol();
                    symbol.inbound_connections = BackwardConnections(symbol.start_atom);
                    symbol.outbound_connections = ForwardConnections(symbol.end_atom);
                    symbol.active = true;
                    wln_vector.push_back(symbol);
                    chain.clear();
                }
            }

        }
    }
    return true;
}
bool WriteWLN::ParseAtoms(std::vector<WLNSymbol> &wln_vector) {

    for (unsigned int i=0; i< global_atom_stack.size();i++){
        AtomContainer atom = global_atom_stack.at(i);


        // Special Dual Cases

        //li
        if (atom.atom_label == 3 && !atom.isinring){
            WLNSymbol lithium;
            lithium.symbol = '$'; // this can be my multi denotion.
            lithium.multi_symbol  = std::vector<char> {'-','L','I','-'};
            lithium.terminator = false;
            lithium.wln_valence = 0;
            lithium.start_atom = lithium.end_atom = atom.atom_index;
            lithium.inbound_connections = BackwardConnections(lithium.start_atom);
            lithium.outbound_connections = ForwardConnections(lithium.end_atom);
            lithium.active = true;
            wln_vector.push_back(lithium);
            continue;
        }

        //Be
        if (atom.atom_label == 4 && !atom.isinring){
            WLNSymbol bery;
            bery.symbol = '$'; // this can be my multi denotion.
            bery.multi_symbol  = std::vector<char> {'-','B','E','-'};
            bery.terminator = false;
            bery.wln_valence = 0;
            bery.start_atom = bery.end_atom = atom.atom_index;
            bery.inbound_connections = BackwardConnections(bery.start_atom);
            bery.outbound_connections = ForwardConnections(bery.end_atom);
            bery.active = true;
            wln_vector.push_back(bery);
            continue;
        }

        //Na
        if (atom.atom_label == 11 && !atom.isinring){
            WLNSymbol sodium;
            sodium.symbol = '$'; // this can be my multi denotion.
            sodium.multi_symbol  = std::vector<char> {'-','N','A','-'};
            sodium.terminator = false;
            sodium.wln_valence = 0;
            sodium.start_atom = sodium.end_atom = atom.atom_index;
            sodium.inbound_connections = BackwardConnections(sodium.start_atom);
            sodium.outbound_connections = ForwardConnections(sodium.end_atom);
            sodium.active = true;
            wln_vector.push_back(sodium);
            continue;
        }

        //Mg
        if (atom.atom_label == 12 && !atom.isinring){
            WLNSymbol magnesium;
            magnesium.symbol = '$'; // this can be my multi denotion.
            magnesium.multi_symbol  = std::vector<char> {'-','M','G','-'};
            magnesium.terminator = false;
            magnesium.wln_valence = 0;
            magnesium.start_atom = magnesium.end_atom = atom.atom_index;
            magnesium.inbound_connections = BackwardConnections(magnesium.start_atom);
            magnesium.outbound_connections = ForwardConnections(magnesium.end_atom);
            magnesium.active = true;
            wln_vector.push_back(magnesium);
            continue;
        }

        //Al
        if (atom.atom_label == 13 && !atom.isinring){
            WLNSymbol alum;
            alum.symbol = '$'; // this can be my multi denotion.
            alum.multi_symbol  = std::vector<char> {'-','A','L','-'};
            alum.terminator = false;
            alum.wln_valence = 0;
            alum.start_atom = alum.end_atom = atom.atom_index;
            alum.inbound_connections = BackwardConnections(alum.start_atom);
            alum.outbound_connections = ForwardConnections(alum.end_atom);
            alum.active = true;
            wln_vector.push_back(alum);
            continue;
        }

        //Si
        if (atom.atom_label == 14 && !atom.isinring){
            WLNSymbol silicon;
            silicon.symbol = '$'; // this can be my multi denotion.
            silicon.multi_symbol  = std::vector<char> {'-','S','I','-'};
            silicon.terminator = false;
            silicon.wln_valence = 0;
            silicon.start_atom = silicon.end_atom = atom.atom_index;
            silicon.inbound_connections = BackwardConnections(silicon.start_atom);
            silicon.outbound_connections = ForwardConnections(silicon.end_atom);
            silicon.active = true;
            wln_vector.push_back(silicon);
            continue;
        }

        //K
        if (atom.atom_label == 19 && !atom.isinring){
            WLNSymbol potatassium;
            potatassium.symbol = '$'; // this can be my multi denotion.
            potatassium.multi_symbol  = std::vector<char> {'-','K','A','-'};
            potatassium.terminator = false;
            potatassium.wln_valence = 0;
            potatassium.start_atom = potatassium.end_atom = atom.atom_index;
            potatassium.inbound_connections = BackwardConnections(potatassium.start_atom);
            potatassium.outbound_connections = ForwardConnections(potatassium.end_atom);
            potatassium.active = true;
            wln_vector.push_back(potatassium);
            continue;
        }

        //CAL
        if (atom.atom_label == 20 && !atom.isinring){
            WLNSymbol calcium;
            calcium.symbol = '$'; // this can be my multi denotion.
            calcium.multi_symbol  = std::vector<char> {'-','C','A','-'};
            calcium.terminator = false;
            calcium.wln_valence = 0;
            calcium.start_atom = calcium.end_atom = atom.atom_index;
            calcium.inbound_connections = BackwardConnections(calcium.start_atom);
            calcium.outbound_connections = ForwardConnections(calcium.end_atom);
            calcium.active = true;
            wln_vector.push_back(calcium);
            continue;
        }

        //Scand
        if (atom.atom_label == 21 && !atom.isinring){
            WLNSymbol scand;
            scand.symbol = '$'; // this can be my multi denotion.
            scand.multi_symbol  = std::vector<char> {'-','S','C','-'};
            scand.terminator = false;
            scand.wln_valence = 0;
            scand.start_atom = scand.end_atom = atom.atom_index;
            scand.inbound_connections = BackwardConnections(scand.start_atom);
            scand.outbound_connections = ForwardConnections(scand.end_atom);
            scand.active = true;
            wln_vector.push_back(scand);
            continue;
        }

        //Titanium
        if (atom.atom_label == 22 && !atom.isinring){
            WLNSymbol titanium;
            titanium.symbol = '$'; // this can be my multi denotion.
            titanium.multi_symbol  = std::vector<char> {'-','T','I','-'};
            titanium.terminator = false;
            titanium.wln_valence = 0;
            titanium.start_atom = titanium.end_atom = atom.atom_index;
            titanium.inbound_connections = BackwardConnections(titanium.start_atom);
            titanium.outbound_connections = ForwardConnections(titanium.end_atom);
            titanium.active = true;
            wln_vector.push_back(titanium);
            continue;
        }

        //Vanadium
        if (atom.atom_label == 23 && !atom.isinring){
            WLNSymbol vanadium;
            vanadium.symbol = '$'; // this can be my multi denotion.
            vanadium.multi_symbol  = std::vector<char> {'-','V','A','-'};
            vanadium.terminator = false;
            vanadium.wln_valence = 0;
            vanadium.start_atom = vanadium.end_atom = atom.atom_index;
            vanadium.inbound_connections = BackwardConnections(vanadium.start_atom);
            vanadium.outbound_connections = ForwardConnections(vanadium.end_atom);
            vanadium.active = true;
            wln_vector.push_back(vanadium);
            continue;
        }

        //Chromium
        if (atom.atom_label == 24 && !atom.isinring){
            WLNSymbol chromium;
            chromium.symbol = '$'; // this can be my multi denotion.
            chromium.multi_symbol  = std::vector<char> {'-','C','R','-'};
            chromium.terminator = false;
            chromium.wln_valence = 0;
            chromium.start_atom = chromium.end_atom = atom.atom_index;
            chromium.inbound_connections = BackwardConnections(chromium.start_atom);
            chromium.outbound_connections = ForwardConnections(chromium.end_atom);
            chromium.active = true;
            wln_vector.push_back(chromium);
            continue;
        }

        //manganese
        if (atom.atom_label == 25 && !atom.isinring){
            WLNSymbol manganese;
            manganese.symbol = '$'; // this can be my multi denotion.
            manganese.multi_symbol  = std::vector<char> {'-','M','N','-'};
            manganese.terminator = false;
            manganese.wln_valence = 0;
            manganese.start_atom = manganese.end_atom = atom.atom_index;
            manganese.inbound_connections = BackwardConnections(manganese.start_atom);
            manganese.outbound_connections = ForwardConnections(manganese.end_atom);
            manganese.active = true;
            wln_vector.push_back(manganese);
            continue;
        }

        //iron
        if (atom.atom_label == 26 && !atom.isinring){
            WLNSymbol iron;
            iron.symbol = '$'; // this can be my multi denotion.
            iron.multi_symbol  = std::vector<char> {'-','F','E','-'};
            iron.terminator = false;
            iron.wln_valence = 0;
            iron.start_atom = iron.end_atom = atom.atom_index;
            iron.inbound_connections = BackwardConnections(iron.start_atom);
            iron.outbound_connections = ForwardConnections(iron.end_atom);
            iron.active = true;
            wln_vector.push_back(iron);
            continue;
        }

        //cobalt
        if (atom.atom_label == 27 && !atom.isinring){
            WLNSymbol cobalt;
            cobalt.symbol = '$'; // this can be my multi denotion.
            cobalt.multi_symbol  = std::vector<char> {'-','C','O','-'};
            cobalt.terminator = false;
            cobalt.wln_valence = 0;
            cobalt.start_atom = cobalt.end_atom = atom.atom_index;
            cobalt.inbound_connections = BackwardConnections(cobalt.start_atom);
            cobalt.outbound_connections = ForwardConnections(cobalt.end_atom);
            cobalt.active = true;
            wln_vector.push_back(cobalt);
            continue;
        }

        //nickel
        if (atom.atom_label == 28 && !atom.isinring){
            WLNSymbol nickel;
            nickel.symbol = '$'; // this can be my multi denotion.
            nickel.multi_symbol  = std::vector<char> {'-','N','I','-'};
            nickel.terminator = false;
            nickel.wln_valence = 0;
            nickel.start_atom = nickel.end_atom = atom.atom_index;
            nickel.inbound_connections = BackwardConnections(nickel.start_atom);
            nickel.outbound_connections = ForwardConnections(nickel.end_atom);
            nickel.active = true;
            wln_vector.push_back(nickel);
            continue;
        }

        //copper
        if (atom.atom_label == 29 && !atom.isinring){
            WLNSymbol copper;
            copper.symbol = '$'; // this can be my multi denotion.
            copper.multi_symbol  = std::vector<char> {'-','C','U','-'};
            copper.terminator = false;
            copper.wln_valence = 0;
            copper.start_atom = copper.end_atom = atom.atom_index;
            copper.inbound_connections = BackwardConnections(copper.start_atom);
            copper.outbound_connections = ForwardConnections(copper.end_atom);
            copper.active = true;
            wln_vector.push_back(copper);
            continue;
        }

        //Zinc
        if (atom.atom_label == 30 && !atom.isinring){
            WLNSymbol zinc;
            zinc.symbol = '$'; // this can be my multi denotion.
            zinc.multi_symbol  = std::vector<char> {'-','Z','N','-'};
            zinc.terminator = false;
            zinc.wln_valence = 0;
            zinc.start_atom = zinc.end_atom = atom.atom_index;
            zinc.inbound_connections = BackwardConnections(zinc.start_atom);
            zinc.outbound_connections = ForwardConnections(zinc.end_atom);
            zinc.active = true;
            wln_vector.push_back(zinc);
            continue;
        }


        // Standard SYM atoms

        // Carbons
        if (atom.atom_label == 6){
            if (atom.connected_bonds.size() == 4 && !atom.isinring) {
                WLNSymbol x_char;
                x_char.symbol = 'X';
                x_char.wln_valence = 4;
                x_char.terminator = false;
                x_char.start_atom = x_char.end_atom = atom.atom_index;
                x_char.inbound_connections = BackwardConnections(x_char.start_atom);
                x_char.outbound_connections = ForwardConnections(x_char.end_atom);
                x_char.active = true;
                wln_vector.push_back(x_char);
                continue;
            }
            if (atom.connected_bonds.size() == 3 && !atom.isinring) {
                WLNSymbol y_char;
                y_char.symbol = 'Y';
                y_char.wln_valence = 3;
                y_char.terminator = false;
                y_char.start_atom = y_char.end_atom = atom.atom_index;
                y_char.inbound_connections = BackwardConnections(y_char.start_atom);
                y_char.outbound_connections = ForwardConnections(y_char.end_atom);
                y_char.active= true;
                wln_vector.push_back(y_char);
                continue;
            }
        }

        // Boron
        if (atom.atom_label ==5 && !atom.isinring){
            WLNSymbol boron;
            boron.symbol = 'B';
            boron.wln_valence  = 3;
            boron.terminator = false;
            boron.start_atom = boron.end_atom = atom.atom_index;
            boron.inbound_connections = BackwardConnections(boron.start_atom);
            boron.outbound_connections = ForwardConnections(boron.end_atom);
            boron.active = true;
            wln_vector.push_back(boron);
            continue;
        }

        // Nitrogen
        if (atom.atom_label ==7 && !atom.isinring){
            WLNSymbol nitrogen;
            if (atom.connected_bonds.size() == 4) {
                nitrogen.symbol = 'K';
                nitrogen.wln_valence  = 4;
                nitrogen.terminator = false;
            }
            if (atom.connected_bonds.size() == 3) {
                nitrogen.symbol = 'N';
                nitrogen.wln_valence  = 3;
                nitrogen.terminator = false;
            }

            if (atom.connected_bonds.size() == 2) {
                nitrogen.symbol = 'M';
                nitrogen.wln_valence = 0;
                nitrogen.terminator = false;
            }

            if (atom.connected_bonds.size() == 1) {
                if (atom.connected_bonds.front().bond_order == 2){
                    nitrogen.symbol = 'M';
                    nitrogen.wln_valence  = 0;
                    nitrogen.terminator = false;
                }
                else{
                    nitrogen.symbol = 'Z';
                    nitrogen.wln_valence  = 0;
                    nitrogen.terminator = true;
                }
            }
            nitrogen.start_atom = nitrogen.end_atom = atom.atom_index;
            nitrogen.inbound_connections = BackwardConnections(nitrogen.start_atom);
            nitrogen.outbound_connections = ForwardConnections(nitrogen.end_atom);
            nitrogen.active=true;
            wln_vector.push_back(nitrogen);
            continue;
        }

        // Oxygen
        if (atom.atom_label ==8 && !atom.isinring){
            WLNSymbol oxygen;

            if (atom.connected_bonds.size() > 1) {
                oxygen.symbol = 'O';
                oxygen.wln_valence  = 0;
                oxygen.terminator = false;
            }

            if (atom.connected_bonds.size() == 1) {
                if (atom.connected_bonds.front().bond_order == 2){
                    oxygen.symbol = 'O';
                    oxygen.wln_valence  = 0;
                    oxygen.terminator = false;
                }
                else{
                    oxygen.symbol = 'Q';
                    oxygen.wln_valence  = 0;
                    oxygen.terminator = true;
                }
            }
            oxygen.start_atom = oxygen.end_atom = atom.atom_index;
            oxygen.inbound_connections = BackwardConnections(oxygen.start_atom);
            oxygen.outbound_connections = ForwardConnections(oxygen.end_atom);
            oxygen.active=true;
            wln_vector.push_back(oxygen);
            continue;
        }

        // Florine
        if (atom.atom_label ==9 && !atom.isinring){
            WLNSymbol fluorine;
            fluorine.symbol = 'F';
            fluorine.wln_valence = 0;
            fluorine.terminator = true;
            fluorine.start_atom = fluorine.end_atom = atom.atom_index;
            fluorine.inbound_connections = BackwardConnections(fluorine.start_atom);
            fluorine.outbound_connections = ForwardConnections(fluorine.end_atom);
            fluorine.active = true;
            wln_vector.push_back(fluorine);
            continue;
        }

        // Phosphorus
        if (atom.atom_label ==15 && !atom.isinring){
            WLNSymbol phos;
            phos.symbol = 'P';
            phos.wln_valence = 3;
            phos.terminator = false;
            phos.start_atom = phos.end_atom = atom.atom_index;
            phos.inbound_connections = BackwardConnections(phos.start_atom);
            phos.outbound_connections = ForwardConnections(phos.end_atom);
            phos.active = true;
            wln_vector.push_back(phos);
            continue;
        }

        // Sulphur
        if (atom.atom_label ==16 && !atom.isinring){
            WLNSymbol sulphur;
            sulphur.symbol = 'S';
            sulphur.wln_valence = 3;
            sulphur.terminator = false;
            sulphur.start_atom = sulphur.end_atom = atom.atom_index;
            sulphur.inbound_connections = BackwardConnections(sulphur.start_atom);
            sulphur.outbound_connections = ForwardConnections(sulphur.end_atom);
            sulphur.active = true;
            wln_vector.push_back(sulphur);
            continue;
        }

        // Chlorine
        if (atom.atom_label ==17 && !atom.isinring){
            WLNSymbol chlorine;
            chlorine.symbol = 'G';
            chlorine.wln_valence  = 0;
            chlorine.terminator = true;
            chlorine.start_atom = chlorine.end_atom = atom.atom_index;
            chlorine.inbound_connections = BackwardConnections(chlorine.start_atom);
            chlorine.outbound_connections = ForwardConnections(chlorine.end_atom);
            chlorine.active = true;
            wln_vector.push_back(chlorine);
            continue;
        }

        // Bromine
        if (atom.atom_label ==35 && !atom.isinring){
            WLNSymbol bromine;
            bromine.symbol = 'E';
            bromine.wln_valence  = 0;
            bromine.terminator = true;
            bromine.start_atom = bromine.end_atom = atom.atom_index;
            bromine.inbound_connections = BackwardConnections(bromine.start_atom);
            bromine.outbound_connections = ForwardConnections(bromine.end_atom);
            bromine.active = true;
            wln_vector.push_back(bromine);
            continue;
        }

        // Iodine
        if (atom.atom_label ==53 && !atom.isinring){
            WLNSymbol Iodine;
            Iodine.symbol = 'I';
            Iodine.wln_valence = 0;
            Iodine.terminator = true;
            Iodine.start_atom = Iodine.end_atom = atom.atom_index;
            Iodine.inbound_connections = BackwardConnections(Iodine.start_atom);
            Iodine.outbound_connections = ForwardConnections(Iodine.end_atom);
            Iodine.active = true;
            wln_vector.push_back(Iodine);
            continue;
        }
    }

    return true;
}
bool WriteWLN::FormSystem(std::vector<WLNSymbol> &wln_vector, std::vector<char> &wln_string){
    std::sort(wln_vector.begin(),
              wln_vector.end(),
              [](const WLNSymbol& lhs, const WLNSymbol& rhs)
              {
                  return lhs.inbound_connections.size() < rhs.inbound_connections.size();
              });


    // Generate positions for each Symbol plus remove 1 values that are connected to X/Y
    for (unsigned int i=0; i <wln_vector.size();i++){
        // Positional
        wln_vector.at(i).pos = i;

        // The lone carbon
        if (wln_vector.at(i).symbol=='1' && wln_vector.at(i).outbound_connections.empty() && wln_vector.at(i).inbound_connections.empty()){
            continue;
        }


        // The bonds we're looking for should only have one term, therefore only have a front...
        if (wln_vector.at(i).symbol=='1' && wln_vector.at(i).outbound_connections.empty()){
            std::tuple<WLNSymbol,int> previous = PreviousSymbols(wln_vector.at(i), wln_vector).front();
            int order;
            WLNSymbol prev_symbol;
            std::tie(prev_symbol,order) = previous;
            if (order==1 && (prev_symbol.symbol=='X' || prev_symbol.symbol=='Y') ){
                wln_vector.erase(wln_vector.begin() + i);
                i+= -1;
            }
        }
    }


    // Chain Stack and Decompose Block
    std::stack<WLNSymbol> branch_stack;
    std::stack<WLNSymbol> ring_stack;
    std::stack<WLNSymbol> junctions;
    std::vector <std::tuple<WLNSymbol,int>> next_symbols_and_bonds;
    WLNSymbol current{};
    WLNSymbol hold_junction{};
    WLNSymbol next_symbol;
    int order;
    for(;;){
        // Inbound connection 0 should be a unique string start case.
        if(current.inbound_connections.empty()){
            if (!current.active)
                current = WLNSymbol(wln_vector.front());

            if (current.symbol == '*'){
                if (!ring_stack.empty()){
                    if(current.pos != ring_stack.top().pos){
                        if(!BuildRing(current,wln_string))
                            return false;
                        ring_stack.push(current);
                        RemovePos(current.pos, wln_vector);
                    }
                }
                if (ring_stack.empty()){
                    if(!BuildRing(current,wln_string))
                        return false;
                    ring_stack.push(current);
                    RemovePos(current.pos, wln_vector);
                }
                if (!current.outbound_connections.empty()){
                    wln_string.push_back(' '); // opens the ring assigment notation
                    next_symbols_and_bonds = NextSymbols(current,wln_vector);
                    std::tie(next_symbol,order) = next_symbols_and_bonds.front();
                    char ring_pos = FetchRingPosition(current.outbound_connections.front().atom_2_index, current.contained_ring, 'B');
                    wln_string.push_back(ring_pos); // opens the ring assigment notation
                    CheckOrder(order, wln_string);
                    current.outbound_connections.erase(current.outbound_connections.begin()); // Takes the bond away.
                    ring_stack.top() = current; // overwrite happens here;

                    // Double ring clause with single binding bond.
                    if (next_symbol.symbol == '*'){
                        StartInline(next_symbol, wln_string);
                        current = WLNSymbol(next_symbol);
                    }else
                        current = WLNSymbol(next_symbol);

                    continue;
                }

                // This is where the rings end, so all terminator characters MUST be handled here.
                if (current.outbound_connections.empty()){
                    ring_stack.pop();
                    if (ring_stack.empty() && junctions.empty())
                        return true;

                    if (ring_stack.empty() && !junctions.empty()){
                        current = WLNSymbol(junctions.top());
                        next_symbols_and_bonds = NextSymbols(current, wln_vector);
                        if (next_symbols_and_bonds.empty())
                            return true;
                        wln_string.push_back('&'); // This is required.
                        continue;
                    }

                    if (!ring_stack.empty() && junctions.empty()){
                        wln_string.push_back('&');
                        current = WLNSymbol(ring_stack.top());
                        continue;
                    }

                    if (!ring_stack.empty() && !junctions.empty()){
                        current = WLNSymbol(junctions.top());
                        next_symbols_and_bonds = NextSymbols(current, wln_vector); // isn't assigned it's mainly a check for return value
                        if (next_symbols_and_bonds.empty()){

                            // place on ring
                            current = WLNSymbol(ring_stack.top());
                            next_symbols_and_bonds = NextSymbols(current, wln_vector); // isn

                            if (next_symbols_and_bonds.empty())
                                return true;

                            continue;
                        }

                        wln_string.push_back('&'); // This is required.
                        continue;
                    }
                }
            }
            else{
                if (current.symbol=='$'){
                    wln_string.insert(wln_string.end(), current.multi_symbol.begin(),current.multi_symbol.end());
                }else
                    wln_string.push_back(current.symbol);
                next_symbols_and_bonds = NextSymbols(current,wln_vector);
                if (next_symbols_and_bonds.empty()){
                    if (ring_stack.empty())
                        return true;
                    else{
                        wln_string.push_back('&');
                        current = WLNSymbol(ring_stack.top());
                        branch_stack = std::stack<WLNSymbol>(); // clears the branch stack so the loop can iterate as normal.
                        continue;
                    }
                }

                std::tie(next_symbol,order) = next_symbols_and_bonds.front();
                CheckOrder(order, wln_string);
                current = WLNSymbol(next_symbol);

            }
        }


        if(branch_stack.empty()){
            if (current.wln_valence>1){
                // Starting a Chain split character
                current.wln_valence += -1;
                branch_stack.push(current);
                hold_junction = current; // hold this junction for rings.
                if (!junctions.empty()){
                    if (current.pos != junctions.top().pos){
                        if (current.symbol=='$'){
                            wln_string.insert(wln_string.end(), current.multi_symbol.begin(),current.multi_symbol.end());
                        }else
                            wln_string.push_back(current.symbol);
                    }
                }
                else{
                    if (current.symbol=='$'){
                        wln_string.insert(wln_string.end(), current.multi_symbol.begin(),current.multi_symbol.end());
                    }else
                        wln_string.push_back(current.symbol);
                }

                next_symbols_and_bonds = NextSymbols(current,wln_vector);
                if (next_symbols_and_bonds.empty()){
                    if (ring_stack.empty())
                        return true;
                    else{
                        wln_string.push_back('&');
                        current = WLNSymbol(ring_stack.top());
                        branch_stack = std::stack<WLNSymbol>(); // clears the branch stack so the loop can iterate as normal.
                        continue;
                    }
                }

                std::tie(next_symbol,order) = next_symbols_and_bonds.front();
                CheckOrder(order, wln_string);

                if (next_symbol.symbol == '*'){
                    StartInline(next_symbol, wln_string);
                    current = WLNSymbol(next_symbol);
                    if (hold_junction.active)
                        PushJunction(hold_junction, junctions);
                }else
                    current = WLNSymbol(next_symbol);
                continue;
            }
            if (current.wln_valence < 1){
                // Flows into a chain character
                if (current.symbol=='$'){
                    wln_string.insert(wln_string.end(), current.multi_symbol.begin(),current.multi_symbol.end());
                }else
                    wln_string.push_back(current.symbol);
                next_symbols_and_bonds = NextSymbols(current,wln_vector);

                if (next_symbols_and_bonds.empty()){
                    if (ring_stack.empty())
                        return true;
                    else{
                        current = WLNSymbol(ring_stack.top());
                        branch_stack = std::stack<WLNSymbol>(); // clears the branch stack so the loop can iterate as normal.
                        continue;
                    }
                }

                std::tie(next_symbol,order) = next_symbols_and_bonds.front();
                CheckOrder(order, wln_string);

                if (next_symbol.symbol == '*'){
                    StartInline(next_symbol, wln_string);
                    current = WLNSymbol(next_symbol);
                }else
                    current = WLNSymbol(next_symbol);


                continue;
            }
        }

        if(!branch_stack.empty()){
            if (current.wln_valence >= 1){
                if (current.pos == branch_stack.top().pos){ // Check if we've looped back, wln_valence should of been handled
                    next_symbols_and_bonds = NextSymbols(current,wln_vector);
                    if (next_symbols_and_bonds.empty()){
                        branch_stack.pop();
                        RemovePos(current.pos, wln_vector);
                        if (!branch_stack.empty()){
                            branch_stack.top().wln_valence += -1;
                            current = WLNSymbol(branch_stack.top());
                            continue;
                        }else{
                            if (ring_stack.empty()){
                                return true;
                            }
                            else{
                                //SanitiseString(wln_string);
                                current = WLNSymbol(ring_stack.top());
                                branch_stack = std::stack<WLNSymbol>(); // clears the branch stack so the loop can iterate as normal.
                                continue;
                            }
                        }

                    }else{

                        std::tie(next_symbol,order) = next_symbols_and_bonds.front();
                        CheckOrder(order, wln_string);


                        if (next_symbol.symbol == '*'){
                            StartInline(next_symbol, wln_string);
                            current = WLNSymbol(next_symbol);
                            branch_stack = std::stack<WLNSymbol>();
                        }else
                            current = WLNSymbol(next_symbol);
                        continue;
                    }
                }

                if (current.pos != branch_stack.top().pos){
                    current.wln_valence += -1;
                    branch_stack.push(current);
                    if (current.symbol=='$'){
                        wln_string.insert(wln_string.end(), current.multi_symbol.begin(),current.multi_symbol.end());
                    }else
                        wln_string.push_back(current.symbol);
                    next_symbols_and_bonds = NextSymbols(current,wln_vector);

                    if (next_symbols_and_bonds.empty())
                        return true;
                    if (next_symbols_and_bonds.empty()){
                        if (ring_stack.empty())
                            return true;
                        else{
                            wln_string.push_back('&');
                            current = WLNSymbol(ring_stack.top());
                            branch_stack = std::stack<WLNSymbol>(); // clears the branch stack so the loop can iterate as normal.
                            continue;
                        }
                    }

                    std::tie(next_symbol,order) = next_symbols_and_bonds.front();
                    CheckOrder(order, wln_string);
                    current = WLNSymbol(next_symbol);
                    continue;
                }
            }
            if (current.wln_valence < 1){
                // Where a chain character has maxed out its bonds
                if (current.pos == branch_stack.top().pos){
                    branch_stack.pop();
                    RemovePos(current.pos, wln_vector);
                    if (branch_stack.empty()){
                        if (ring_stack.empty())
                            return true;
                        else{
                            current = WLNSymbol(ring_stack.top());
                            branch_stack = std::stack<WLNSymbol>(); // clears the branch stack so the loop can iterate as normal.
                            continue;
                        }
                    }
                    if (!branch_stack.empty()){
                        branch_stack.top().wln_valence += -1;
                        current = WLNSymbol(branch_stack.top());
                        continue;
                    }
                }
                next_symbols_and_bonds = NextSymbols(current,wln_vector);
                if (current.symbol=='$'){
                    wln_string.insert(wln_string.end(), current.multi_symbol.begin(),current.multi_symbol.end());
                }else
                    wln_string.push_back(current.symbol);
                if (next_symbols_and_bonds.empty()){
                    if(!current.terminator) // terminator not needed for chain ending characters
                        wln_string.push_back('&'); // splits chain
                    RemovePos(current.pos, wln_vector);
                    branch_stack.top().wln_valence += -1;  // minus the degree from the branch stack, then copy.
                    current = WLNSymbol(branch_stack.top());
                    continue;
                }
                if (!next_symbols_and_bonds.empty()){
                    // If the chain continues more than one symbol
                    RemovePos(current.pos, wln_vector);
                    std::tie(next_symbol,order) = next_symbols_and_bonds.front();
                    CheckOrder(order, wln_string);

                    if (next_symbol.symbol == '*'){
                        StartInline(next_symbol, wln_string);
                        current = WLNSymbol(next_symbol);
                        if (hold_junction.active)
                            PushJunction(hold_junction, junctions);
                    }else
                        current = WLNSymbol(next_symbol);

                    continue;
                }
            }
        }
    }
}
bool WriteWLN::SanitiseString(std::vector<char> &wln_string){
    if (wln_string.back() == '&')
        wln_string.pop_back();


    std::vector<char> contracted_string{};
    // Contractions into shorthand characters
    unsigned int iter=0;
    unsigned int size =wln_string.size()-1;

    do{
        char focus = wln_string.at(iter);

        switch(focus){

            case 'A':
            case 'B':
            case 'C':
            case 'D':
            case 'E':
            case 'F':
            case 'G':
            case 'H':
            case 'I':
            case 'J':
            case 'K':
            case 'M':
            case 'N':
            case 'O':
            case 'P':
            case 'Q':
            case 'R':
            case 'S':
            case 'T':
            case 'V':
            case 'W':
            case 'X':
            case 'Z':
            case '&':
            case '-':
            case ' ':
            case '/':
            case '0':
                contracted_string.push_back(wln_string.at(iter));
                iter++;
                continue;

            case '1':
            case '2':
            case '3':
            case '4':
            case '5':
            case '6':
            case '7':
            case '8':
            case '9':
                if (size-iter >= 2){
                    if(wln_string.at(iter+1) == 'U' &&
                       wln_string.at(iter+2) == 'O'){
                        contracted_string.push_back('V');
                        iter+=3;
                        continue;

                    }
                }

                contracted_string.push_back(wln_string.at(iter));
                iter++;
                continue;

            case 'L':
                if (size-iter >= 2){
                    if(wln_string.at(iter+1) == '6' &&
                       wln_string.at(iter+2) == 'J'){
                        contracted_string.push_back('R');
                        iter+=3;
                        continue;

                    }
                }

                contracted_string.push_back(wln_string.at(iter));
                iter++;
                continue;


            case 'U':
                if (size-iter >= 4){
                    if(wln_string.at(iter+1) == 'O' &&
                       wln_string.at(iter+2) == '&' &&
                       wln_string.at(iter+3) == 'U' &&
                       wln_string.at(iter+4) == 'O'){

                        contracted_string.push_back('W');
                        iter+=5;
                        continue;

                    }
                }
                contracted_string.push_back(wln_string.at(iter));
                iter++;
                continue;


            case 'Y':
                // highest order must be at the top of each case
                if (size-iter >= 3){
                    if(wln_string.at(iter+1) == 'U' &&
                       wln_string.at(iter+2) == 'O' &&
                       wln_string.at(iter+3) == '&'){

                        contracted_string.push_back('V');
                        iter+=4;
                        continue;
                    }
                }

                if (size-iter >= 2) {
                    if (wln_string.at(iter + 1) == 'U' &&
                        wln_string.at(iter + 2) == 'O') {
                        contracted_string.push_back('V');
                        iter += 3;
                        continue;
                    }
                }

                contracted_string.push_back(wln_string.at(iter));
                iter++;
                continue;

        }




    }while(iter<wln_string.size());


    wln_string = contracted_string;
    return true;
}


//  --- Output Function ---
bool WriteWLN::CreateString(std::vector<char> &wln_string){

    std::vector<WLNSymbol> wln_vector;
    if(!BuildGlobalStacks(Mol))
        return Error("stacks");

    if (global_ring_stack.size()==0){
        if (!SegmentAtomStack(wln_vector))
            return Error("segment");
        if (!ParseAtoms(wln_vector))
            return Error("parser");
        if (!FormSystem(wln_vector, wln_string))
            return Error("system");
        if (!SanitiseString(wln_string))
            return Error("sanitise");
    }

    if (global_ring_stack.size()>0){
        if (!SegmentAtomStack(wln_vector))
            return Error("segment");
        if (!ParseAtoms(wln_vector))
            return Error("parser");
        if (!LinkRings(wln_vector))
            return Error("ring-linker");
        if (!CreateVirtualRings(wln_vector))
            return Error("ring-merger");
        if (!FormSystem(wln_vector, wln_string))
            return Error("system");
        if (!SanitiseString(wln_string))
            return Error("sanitise");
    }

    else
        return Error("unknown");
}
bool WriteWLN::Error(const char* type){
    fprintf(stderr, "Error in wln write\n");
    fprintf(stderr, "Type: %s \n", type);
    return false;
}


// --- Functions ---
bool WriteWLN::RemovePos(unsigned int symbol_id, std::vector<WLNSymbol> &wln_vector){
    for (unsigned int i=0; i<wln_vector.size();i++){
        if (wln_vector.at(i).pos == symbol_id){
            wln_vector.erase(wln_vector.begin()+i);
            return true;
        }
    }
    return false;
}
bool WriteWLN::ChainAtom(AtomContainer atom) {
    if (atom.atom_label==6 && atom.connected_bonds.size() <=2 && !atom.isinring) // Now this is global.
        return true;
    return false;
};
bool WriteWLN::Joined(AtomContainer atom_1, AtomContainer atom_2, int order=-1){

    if (order != -1){
        for (unsigned int i=0;i<atom_1.connected_bonds.size();i++){
            BondContainer bond = atom_1.connected_bonds.at(i);
            if (atom_2.atom_index == bond.atom_1_index || atom_2.atom_index == bond.atom_2_index  && bond.bond_order==order)
                return true;
        }
        return false;
    }else{
        for (unsigned int i=0;i<atom_1.connected_bonds.size();i++){
            BondContainer bond = atom_1.connected_bonds.at(i);
            if (atom_2.atom_index == bond.atom_1_index || atom_2.atom_index == bond.atom_2_index)
                return true;
        }
        return false;
    }
};
std::vector<BondContainer> WriteWLN::GetConnections(unsigned int atom_index){
    std::vector<BondContainer> connected_bonds{};
    for (unsigned int i=0; i< global_bond_stack.size();i++){
        BondContainer bond = global_bond_stack.at(i);
        if (bond.atom_1_index == atom_index || bond.atom_2_index == atom_index){
            connected_bonds.push_back(bond);
        }
    }
    return connected_bonds;
}
std::vector<BondContainer> WriteWLN::ForwardConnections(unsigned int atom_index){
    std::vector<BondContainer> forward_bonds{};
    for (unsigned int i=0; i< global_bond_stack.size();i++){
        BondContainer bond = global_bond_stack.at(i);
        if (bond.atom_1_index == atom_index){
            forward_bonds.push_back(bond);
        }
    }
    return forward_bonds;
}
std::vector<BondContainer> WriteWLN::BackwardConnections(unsigned int atom_index){
    std::vector<BondContainer> backward_bonds{};
    for (unsigned int i=0; i< global_bond_stack.size();i++){
        BondContainer bond = global_bond_stack.at(i);
        if (bond.atom_2_index == atom_index){
            backward_bonds.push_back(bond);
        }
    }
    return backward_bonds;
}
std::vector<std::tuple<WLNSymbol,int>> WriteWLN::NextSymbols(WLNSymbol check_symbol, std::vector<WLNSymbol> wln_vector){
    std::vector< std::tuple<WLNSymbol,int>> symbol_and_bonds;
    if(check_symbol.outbound_connections.empty())
        return symbol_and_bonds;

    for (unsigned int i=0;i<check_symbol.outbound_connections.size();i++){
        BondContainer bond = check_symbol.outbound_connections.at(i);
        for (unsigned int j=0;j<wln_vector.size();j++){
            WLNSymbol search_wln = wln_vector.at(j);
            if (search_wln.inbound_connections.empty())
                continue;
            for (unsigned int k=0; k<search_wln.inbound_connections.size();k++){
                BondContainer symbol_bond = search_wln.inbound_connections.at(k);
                if (symbol_bond.stack_index == bond.stack_index){
                    std::tuple<WLNSymbol,int> push{search_wln,bond.bond_order};
                    symbol_and_bonds.push_back(push);
                }
            }
        }
    }

    return symbol_and_bonds;
}
std::vector<std::tuple<WLNSymbol,int>> WriteWLN::PreviousSymbols(WLNSymbol check_symbol, std::vector<WLNSymbol> wln_vector){
    std::vector< std::tuple<WLNSymbol,int>> symbol_and_bonds;
    if(check_symbol.inbound_connections.empty())
        return symbol_and_bonds;

    for (unsigned int i=0;i<check_symbol.inbound_connections.size();i++){
        BondContainer bond = check_symbol.inbound_connections.at(i);
        for (unsigned int j=0;j<wln_vector.size();j++){
            WLNSymbol search_wln = wln_vector.at(j);
            if (search_wln.outbound_connections.empty())
                continue;
            for (unsigned int k=0; k<search_wln.outbound_connections.size();k++){
                BondContainer symbol_bond = search_wln.outbound_connections.at(k);
                if (symbol_bond.stack_index == bond.stack_index){
                    std::tuple<WLNSymbol,int> push{search_wln,bond.bond_order};
                    symbol_and_bonds.push_back(push);
                }
            }
        }
    }
    return symbol_and_bonds;
}
AtomContainer WriteWLN::FetchAtom(unsigned int atom_index){
    AtomContainer blank;
    for (unsigned int i=0;i<global_atom_stack.size();i++){
        AtomContainer atom = global_atom_stack.at(i);
        if(atom.atom_index == atom_index)
            return atom;
    }
    std::cout << "Blank Hit" << std::endl;
    return blank;
}
void WriteWLN::CheckOrder(unsigned int order, std::vector<char> &wln_string){
    if (order==2)
        wln_string.push_back('U');
    if (order==3){
        wln_string.push_back('U');
        wln_string.push_back('U');
    }
}
char WriteWLN::PeriodicCharacter(AtomContainer atom){

    unsigned int total_valence = 0;
    for (unsigned int i=0;i<atom.connected_bonds.size();i++){
        total_valence += atom.connected_bonds.at(i).bond_order;
    }

    // Boron
    if (atom.atom_label ==5)
        return 'B';

    // Nitrogen
    if (atom.atom_label ==7){
        if (total_valence == 4)
            return 'K';
        if (total_valence == 3)
            return 'N';
        if (total_valence == 2)
            return 'M';

        if (total_valence == 1)
            return 'Z';
    }

    // Oxygen
    if (atom.atom_label ==8){
        if (total_valence > 1)
            return 'O';
        if (total_valence > 1)
            return 'Q';
    }

    // Florine
    if (atom.atom_label ==9){
        return 'F';
    }

    // Phosphorus
    if (atom.atom_label ==15){
        return 'P';
    }

    // Sulphur
    if (atom.atom_label ==16){
        return 'S';
    }

    // Chlorine
    if (atom.atom_label ==17){
        return 'G';
    }

    // Bromine
    if (atom.atom_label ==35){
        return 'E';
    }

    // Iodine
    if (atom.atom_label ==53){
        return 'I';
    }

    else
        return '#';
}
std::vector<unsigned int> WriteWLN::SymbolAtomIndexes(WLNSymbol subject){

    std::sort(subject.atom_container.begin(),
              subject.atom_container.end(),
              [](const AtomContainer& lhs, const AtomContainer& rhs)
              {
                  return lhs.atom_index < rhs.atom_index;
              });

    std::vector<unsigned int> subject_indexes{};
    for (unsigned int sub_atom=0;sub_atom < subject.atom_container.size();sub_atom++)
        subject_indexes.push_back(subject.atom_container.at(sub_atom).atom_index);

    return subject_indexes;
}
unsigned int WriteWLN::ExternalBond(AtomContainer atom, std::vector<AtomContainer> atom_list){

    // Changing to hold a priority value for ring bonds over the R groups, although it should return at least something.
    // to mimising will work as expected.

    // returns 2 for joint ring bond not held in the current list
    // return 1 for R group
    // 0 for non bonding atom
    std::vector<unsigned int> indexes;
    for (unsigned int i=0;i<atom_list.size();i++)
        indexes.push_back(atom_list.at(i).atom_index);

    for (unsigned int i=0;i<atom.connected_bonds.size();i++){
        BondContainer bond = atom.connected_bonds.at(i);

        if (atom.atom_index != bond.atom_1_index){
            if( std::find(indexes.begin(), indexes.end(), bond.atom_1_index) != indexes.end())
                continue;
            else{
                if (bond.ring_bond)
                    return 2;
                else
                    return 1;
            }
        }

        if (atom.atom_index != bond.atom_2_index){
            if( std::find(indexes.begin(), indexes.end(), bond.atom_2_index) != indexes.end())
                continue;
            else{
                if (bond.ring_bond)
                    return 2;
                else
                    return 1;
            }
        }

    }
    return 0;
}

template<class T>
std::vector<T> WriteWLN::VectorReorder(std::vector<T> vA, std::vector<int> vOrder )  {

    assert(vA.size() == vOrder.size());
    for( size_t vv = 0; vv < vA.size() - 1; ++vv )
    {
        if (vOrder[vv] == vv)
        {
            continue;
        }
        size_t oo;
        for(oo = vv + 1; oo < vOrder.size(); ++oo)
        {
            if (vOrder[oo] == vv)
            {
                break;
            }
        }
        std::swap( vA[vv], vA[vOrder[vv]] );
        std::swap( vOrder[vv], vOrder[oo] );
    }
    return vA;
}
void WriteWLN::CheckDuplicates(std::vector<AtomContainer> &vec){

    std::vector<unsigned int> remove_indexes;
    std::vector<unsigned int> atom_positions;

    for (unsigned int i=0; i< vec.size(); i++){
        unsigned int atom_pos = vec.at(i).atom_index;
        if (std::find(atom_positions.begin(), atom_positions.end(), atom_pos) != atom_positions.end()) {
            remove_indexes.push_back(i);
        }
        else {
            atom_positions.push_back(atom_pos);

        }
    }
    std::sort(remove_indexes.begin(), remove_indexes.end(), std::greater<unsigned int>());
    for ( int i : remove_indexes ) vec.erase( std::next( vec.begin(), i ) );
}
void WriteWLN::CheckDuplicates(std::vector<BondContainer> &vec){

    std::vector<unsigned int> remove_indexes;
    std::vector<unsigned int> bond_positions;

    for (unsigned int i=0; i< vec.size(); i++){
        unsigned int bond_pos = vec.at(i).stack_index;
        if (std::find(bond_positions.begin(), bond_positions.end(), bond_pos) != bond_positions.end()) {
            remove_indexes.push_back(i);
        }
        else {
            bond_positions.push_back(bond_pos);

        }
    }
    std::sort(remove_indexes.begin(), remove_indexes.end(), std::greater<unsigned int>());
    for ( int i : remove_indexes ) vec.erase( std::next( vec.begin(), i ) );
}


// --- Ring Functions ---
char WriteWLN::FetchRingPosition(unsigned int atom_index, RingContainer ring, char mode){
    if (mode=='B'){
        std::vector<BondContainer> bonds = BackwardConnections(atom_index);
        unsigned int ring_label = bonds.front().atom_1_index;
        // now need to find what relation the ring label has to the rest of the ring.
        for (unsigned int i=0; i<ring.ring_atoms.size();i++){
            if (ring.ring_atoms.at(i).atom_index == ring_label){
                ring_label = i;
                break;
            }
        }

        return (ring_label)  + 'A';
    }

    if (mode=='F'){
        unsigned int ring_label;
        for (unsigned int i=0; i<ring.ring_atoms.size();i++){
            if (ring.ring_atoms.at(i).atom_index == atom_index){
                ring_label = i;
                break;
            }
        }
        return (ring_label  + 'A');
    }
};
bool WriteWLN::LinkRings(std::vector<WLNSymbol> &wln_vector){
    for (unsigned int i=0; i<wln_vector.size();i++){
        for (unsigned int j=0; j<wln_vector.at(i).outbound_connections.size();j++){
            BondContainer bond = wln_vector.at(i).outbound_connections.at(j);
            AtomContainer query = FetchAtom(bond.atom_2_index);
            if (query.isinring)
                wln_vector.at(i).serve_ring += 1;
        }
    }
    return true;
}

bool WriteWLN::BuildRing(WLNSymbol virtual_ring, std::vector<char> &wln_string){

    std::cout << virtual_ring.ring_type << std::endl;

    wln_string.push_back(virtual_ring.contained_ring.hetero);

    if (virtual_ring.ring_type == "NULL")
        return false;

    if (virtual_ring.ring_type == "SINGLE")
        wln_string.push_back(virtual_ring.contained_ring.size + '0');

    if(virtual_ring.ring_type == "BICYCLIC"){
        for (unsigned int i=0; i<virtual_ring.SSRS.size();i++){
            wln_string.push_back(virtual_ring.SSRS.at(i).size + '0');
        }
    }

    if (virtual_ring.ring_type == "POLYCYCLIC"){
        wln_string.push_back(' ');
        for (int i=virtual_ring.poly_chars.size()-1;i>=0;i--){
            wln_string.push_back(virtual_ring.poly_chars.at(i)); // pushes B

            if (virtual_ring.SSRS.size() > 3){
                wln_string.push_back(virtual_ring.SSRS.back().size + '0');
                virtual_ring.SSRS.pop_back();
                wln_string.push_back(' ');
            }else{
                for (unsigned int k = 0; k < 3; k++) { // three is the max in can have in a pure peri cyclic foundation
                    wln_string.push_back(virtual_ring.SSRS.at(k).size + '0');
                    if (k == 2)
                        wln_string.push_back(' ');
                }
            }

        }
    }

    if (virtual_ring.ring_type == "PERICYCLIC") {

        for (unsigned int i = 0; i < 3; i++) { // three is the max in can have in a pure peri cyclic foundation
            wln_string.push_back(virtual_ring.SSRS.at(i).size + '0');
            if (i == 2)
                wln_string.push_back(' ');
        }

        if (virtual_ring.SSRS.size()>3){
            for (unsigned int i=3;i<virtual_ring.SSRS.size();i++){
                wln_string.push_back(virtual_ring.peri_chars.at(i-2));
                wln_string.push_back(virtual_ring.SSRS.at(i).size + '0');
                wln_string.push_back(' ');
            }
        }

        // defines that first peri position.
        wln_string.push_back(virtual_ring.peri_chars.size() + '0');
        for (unsigned int i = 0; i < virtual_ring.peri_chars.size(); i++) {
            wln_string.push_back(virtual_ring.peri_chars.at(i));
            if (i == virtual_ring.peri_chars.size() - 1)
                wln_string.push_back(' ');
        }

        // ends the perinotation on the end atom
        wln_string.push_back(virtual_ring.atom_container.size()-1 + 'A');
    }


    if (virtual_ring.ring_type == "MULTICYCLIC") {
        wln_string.push_back(' ');
        wln_string.push_back(virtual_ring.poly_chars.front());
        wln_string.push_back(virtual_ring.SSRS.back().size + '0');
        virtual_ring.SSRS.pop_back();

        // handle the poly bind pos
        if (virtual_ring.poly_chars.size() == 1){
            wln_string.push_back(virtual_ring.poly_chars.front());
            wln_string.push_back(virtual_ring.SSRS.back().size + '0');
            virtual_ring.SSRS.pop_back();
            for (unsigned int i = 0; i < 3; i++) { // three is the max in can have in a pure peri cyclic foundation
                wln_string.push_back(virtual_ring.SSRS.at(i).size + '0');
                if (i == 2)
                    wln_string.push_back(' ');
            }

            wln_string.push_back(virtual_ring.peri_chars.size() + '0');
            for (unsigned int i = 0; i < virtual_ring.peri_chars.size(); i++) {
                wln_string.push_back(virtual_ring.peri_chars.at(i));
                if (i == virtual_ring.peri_chars.size() - 1)
                    wln_string.push_back(' ');
            }

        }

    }

    if (virtual_ring.ring_type == "BRIDGED"){
        for (unsigned int i=0;i<virtual_ring.SSRS.size();i++){
            wln_string.push_back(virtual_ring.SSRS.at(i).size + '0');
        }
        wln_string.push_back(' ');
        for (unsigned int i=0; i<virtual_ring.poly_chars.size();i++)
            wln_string.push_back(virtual_ring.poly_chars.at(i));
        wln_string.push_back(' ');
    }


    // Internal Bond assignments
    for (unsigned int i=0; i<virtual_ring.contained_ring.ring_bonds.size();i++){
        BondContainer ring_bond = virtual_ring.contained_ring.ring_bonds.at(i);
        if (ring_bond.bond_order == 2 && !ring_bond.aromatic){
            if (ring_bond.atom_1_index - ring_bond.atom_2_index == 1 || ring_bond.atom_2_index - ring_bond.atom_1_index ==1){  // this small condition will help with outside double characters
                wln_string.push_back(' ');
                char internal_pos;
                for (unsigned int k=0; k<virtual_ring.atom_container.size(); k++){
                    AtomContainer atom = virtual_ring.atom_container.at(k);
                    if (atom.atom_index == ring_bond.atom_1_index || atom.atom_index == ring_bond.atom_2_index){
                        internal_pos = k + 'A';
                        break;
                    }
                }
                wln_string.push_back(internal_pos);
                wln_string.push_back('U');
            }
        }
    }


    // Handles hetero ring cases
    if (virtual_ring.contained_ring.hetero=='T'){
        for (unsigned int i=0; i<virtual_ring.contained_ring.ring_atoms.size();i++){
            AtomContainer ring_atom = virtual_ring.contained_ring.ring_atoms.at(i);
            if (ring_atom.atom_label != 6){
                wln_string.push_back(' ');
                char internal_pos = (ring_atom.atom_index - virtual_ring.contained_ring.start_atom.atom_index)  + 'A';
                wln_string.push_back(internal_pos);
                wln_string.push_back(PeriodicCharacter(ring_atom));

            }
        }
    }

    if (virtual_ring.SSRS.size() > 1){
        std::vector<char> unsaturated;
        for (unsigned int i=0; i<virtual_ring.SSRS.size();i++){
            RingContainer ring = virtual_ring.SSRS.at(i);
            if (!ring.aromatic){
                unsaturated.push_back('&');
            }
        }

        // handle unsaturated bonded cycles
        if (!unsaturated.empty()){
            if (unsaturated.size() != virtual_ring.SSRS.size())
                wln_string.insert(wln_string.end(), unsaturated.begin(), unsaturated.end());
        }
        wln_string.push_back('T');

    }
    else{
        if (!virtual_ring.contained_ring.aromatic)
            wln_string.push_back('T');
    }

    wln_string.push_back('J');
    return true;
}

bool WriteWLN::CreateVirtualRings(std::vector<WLNSymbol> &wln_vector){

    std::vector<WLNSymbol> rings{};
    for (unsigned int i=0;i<global_ring_stack.size();i++){
        RingContainer ring = global_ring_stack.at(i);
        WLNSymbol virtual_ring{};
        virtual_ring.symbol = '*';
        virtual_ring.wln_valence =0;
        virtual_ring.active = true;
        virtual_ring.start_atom = ring.ring_atoms.front().atom_index;
        virtual_ring.end_atom = ring.ring_atoms.back().atom_index;
        virtual_ring.atom_container = ring.ring_atoms;
        virtual_ring.contained_ring = ring;
        virtual_ring.ring_type = "SINGLE";

        for (unsigned int i=0;i<ring.active_atoms.size();i++){
            std::vector<BondContainer> inbound_check = BackwardConnections(std::get<0>(ring.active_atoms.at(i)).atom_index);
            for (unsigned int k=0;k<inbound_check.size();k++){
                if (!inbound_check.at(k).ring_bond)
                    virtual_ring.inbound_connections.push_back(inbound_check.at(k));
            }
        }

        for (unsigned int j=0;j<ring.active_atoms.size();j++){
            std::vector<BondContainer> outbound_check = ForwardConnections(std::get<0>(ring.active_atoms.at(j)).atom_index);
            for (unsigned int k=0;k<outbound_check.size();k++){
                if (!outbound_check.at(k).ring_bond)
                    virtual_ring.outbound_connections.push_back(outbound_check.at(k));
            }
        }
        rings.push_back(virtual_ring);
    }


    // Look for rings that are bonded together for bi/poly/peri fused systems.
    // If we do this approach with atoms, it's possible to tease out spiro systems naturally. Even though bonds
    // on the surface make more intuitive sense.


    // give it some positional help
    for (unsigned int i = 0; i < rings.size(); i++) {
        rings.at(i).pos = i;
    }
    std::vector<WLNSymbol> ringvec_1{rings};
    std::vector<WLNSymbol> ringvec_2{rings};
    // new merger that goes based off index sum of the subject ring.
    // should take care of priority paths, needed for peri fusing of multiple rings.
    for(;;){
        unsigned int insertion;
        std::vector<WLNSymbol> condensed{};
        for (unsigned int i = 0; i < ringvec_1.size()-1; i++) {
            WLNSymbol src = ringvec_1.at(i);
            std::vector<unsigned int> bound_indexes;
            std::vector<std::vector<unsigned int>> bound_intersections;
            std::vector<int> bound_score;
            bool bind_available=false;
            for (unsigned int j=i+1; j<ringvec_1.size();j++){
                std::vector<unsigned int> atom_intersection;
                WLNSymbol trg = ringvec_1.at(j);
                if (CheckPolyBinding(src, trg, atom_intersection)) {
                    bind_available = true;
                    int pos_score=0;
                    for (int inter_index=0; inter_index<atom_intersection.size();inter_index++){
                        for (int atom_index=0; atom_index<src.atom_container.size();atom_index++){
                            if (atom_intersection.at(inter_index) == src.atom_container.at(atom_index).atom_index)
                                pos_score += atom_index;
                        }
                    }
                    bound_score.push_back(pos_score);
                    bound_indexes.push_back(j);
                    bound_intersections.push_back(atom_intersection);
                }
            }
            if (bind_available){
                std::vector<int>::iterator it = std::min_element(std::begin(bound_score), std::end(bound_score));
                int lowest_index = std::distance(std::begin(bound_score), it);
                unsigned int j_index = bound_indexes.at(lowest_index);
                std::vector<unsigned int> selected_intersection = bound_intersections.at(lowest_index);
                WLNSymbol selected = ringvec_1.at(j_index);
                WLNSymbol new_ring = MergeRingSymbol(src, selected, selected_intersection); // created the new ring
                condensed.push_back(new_ring);
                RemovePos(src.pos, ringvec_2);
                RemovePos(selected.pos, ringvec_2);
                insertion = i;
                break;
            }
        }

        if (!condensed.empty()){

            // ring build condition
            ringvec_2.insert(ringvec_2.begin() + insertion, condensed.begin(), condensed.end());
            for (unsigned int i = 0; i < ringvec_2.size(); i++) {
                ringvec_2.at(i).pos = i;
            }
            ringvec_1 = ringvec_2;

        }

        if (condensed.empty()){

            // loop end condition
            break;
        }

        // infinite loop ending
    }

    wln_vector.insert(wln_vector.end(), ringvec_1.begin(), ringvec_1.end());
    return true;
}
void WriteWLN::PushJunction(WLNSymbol &hold_junction, std::stack<WLNSymbol> &junctions){

    if (!junctions.empty()){
        if (junctions.top().pos != hold_junction.pos){
            hold_junction.wln_valence +=-1;
            junctions.push(hold_junction);
        }
        else
            junctions.top().wln_valence += -1;
    }
    else{
        hold_junction.wln_valence +=-1;
        junctions.push(hold_junction);
    }
}
void WriteWLN::StartInline(WLNSymbol &next_symbol ,std::vector<char> &wln_string){
    wln_string.push_back('-'); wln_string.push_back(' ');
    char ring_pos = FetchRingPosition(next_symbol.inbound_connections.front().atom_2_index, next_symbol.contained_ring, 'F');
    wln_string.push_back(ring_pos); // opens the ring assigment notation
    next_symbol.inbound_connections.erase(next_symbol.inbound_connections.begin()); // Takes the bond away.
}
bool WriteWLN::CheckPolyBinding(WLNSymbol &subject, WLNSymbol &target, std::vector<unsigned int> &atom_intersection){
    std::vector<unsigned int> subject_indexes = SymbolAtomIndexes(subject);
    std::vector<unsigned int> target_indexes = SymbolAtomIndexes(target);
    std::set_intersection(
            subject_indexes.begin(), subject_indexes.end(),
            target_indexes.begin(), target_indexes.end(),
            std::back_inserter(atom_intersection)
    );

    // if some of the atoms match, this is a merging symbol pair.
    if (!atom_intersection.empty()){
        return true;
    }else
        return false;
}
bool WriteWLN::MultiRing(WLNSymbol subject, WLNSymbol target){
    // careful here, it's a look ahead function
    if ( (subject.SSRS.size() + 1) == (subject.peri_chars.size()+1)+2 )
        return false;
    else
        return true;
}

WLNSymbol WriteWLN::MergeRingSymbol(WLNSymbol subject, WLNSymbol target, std::vector<unsigned int> atom_intersection){

    WLNSymbol merged_ring{};
    merged_ring.symbol = '*';
    merged_ring.wln_valence=0;
    merged_ring.active = true;

    unsigned char binding_char;
    unsigned int binding_pos;
    unsigned int binding_atom_index;

    bool peri = false;
    bool multi = false;
    bool poly = false;

    std::cout << subject.ring_type <<":" << subject.atom_container.size();
    std::cout << "   ";
    std::cout << target.ring_type << ":" << target.atom_container.size();
    std::cout << "   " <<atom_intersection.size() << std::endl;

    // Main path block - always consistent - Merges the subject and target rings
    std::vector<AtomContainer> AllRingAtoms;
    AllRingAtoms.insert(AllRingAtoms.end(), subject.atom_container.begin(), subject.atom_container.end());
    AllRingAtoms.insert(AllRingAtoms.end(), target.atom_container.begin(), target.atom_container.end());
    CheckDuplicates(AllRingAtoms); // Checks the intersection duplicates;

    // do the transfers early for each to keep as standard
    merged_ring.poly_chars.insert(merged_ring.poly_chars.begin(),subject.poly_chars.begin(),subject.poly_chars.end());
    merged_ring.peri_chars.insert(merged_ring.peri_chars.begin(),subject.peri_chars.begin(),subject.peri_chars.end());
    merged_ring.poly_atoms.insert(merged_ring.poly_atoms.begin(),subject.poly_atoms.begin(),subject.poly_atoms.end());
    merged_ring.peri_atoms.insert(merged_ring.peri_atoms.begin(),subject.peri_atoms.begin(),subject.peri_atoms.end());

    // This function is perfect, works very nicely for single-single merge rings.
    if (subject.ring_type == "SINGLE"){

        // "Locant" Path algorithm using Hamiltonian Circuits.
        // Start point = the lowest index out of A and B from intersection.

        // two potential starting points, here, need to generate the paths, reduce, and minimise to always get the lowest potential start.

        int start_index_a;
        int start_index_b;
        for (int i =0; i<AllRingAtoms.size();i++){
            if (AllRingAtoms.at(i).atom_index == atom_intersection.front()){
                start_index_a = i;
            }
            if (AllRingAtoms.at(i).atom_index == atom_intersection.back()){
                start_index_b = i;
            }
        }

        // This keep position really well, now to check individual paths in the graph object.
        // based off set intersection, the indexes should be preserved, so we just need to shift atoms based on intersection index.

        // Build Connection Matrix from AllRingAtoms vector.
        Graph ring_graph(AllRingAtoms.size()); // creates graph of the right size;
        CreateGraph(AllRingAtoms, ring_graph); // Adds all the bonds in relation to position in AllRingAtoms vector
        ring_graph.GeneratePath();                      // Generates the Locant Paths into the graph object. // first one is usually the lowest index.

        // this should all be consistent between the two start indexes.
        ring_graph.ReorderPaths(start_index_a); // reorder based off the start index.

        int first_index = ScorePaths(AllRingAtoms, ring_graph);
        std::vector<AtomContainer> Ordered_1 = VectorReorder(AllRingAtoms, ring_graph.indexes.at(first_index));

        ring_graph.ReorderPaths(start_index_b); // reorder based off the alternative start index.
        int second_index = ScorePaths(AllRingAtoms, ring_graph);
        std::vector<AtomContainer> Ordered_2 = VectorReorder(AllRingAtoms, ring_graph.indexes.at(second_index));

        if(ScoreLists(Ordered_1,Ordered_2) == 0){
            merged_ring.atom_container = Ordered_1;
        }
        if(ScoreLists(Ordered_1,Ordered_2) == 1){
            merged_ring.atom_container = Ordered_2;
        }

        // Allringatoms now contain the fully reordered set for the new symbol
        CloneRingData(merged_ring, subject, target, "BICYCLIC");
        return merged_ring;
    }

    if (subject.ring_type == "BICYCLIC" || subject.ring_type == "POLYCYCLIC"){
        if (atom_intersection.size() == 2){
            poly = true;
            for (int i=0;i<subject.atom_container.size();i++){
                if (subject.atom_container.at(i).atom_index == atom_intersection.front() ||
                    subject.atom_container.at(i).atom_index == atom_intersection.back()){
                    binding_char = i + 'A';
                    binding_pos = i;
                    binding_atom_index = subject.atom_container.at(i).atom_index;
                    break;
                }
            }
        } // standard fuse
        if (atom_intersection.size() == 3){
            peri = true;
            // path must be consistent with the peri going first
            // since these are now split, there is no way for a BICYCLIC or POLYCYCLIC to have peri atoms
            binding_char = 'A';
            binding_pos = 0;
            binding_atom_index = subject.atom_container.front().atom_index;
            merged_ring.peri_atoms.push_back(subject.atom_container.front());
            merged_ring.peri_chars.push_back(binding_char);

        } // standard peri fuse
        if (atom_intersection.size() > 3){
            // now we know whats going on, lets start generalising here. need the ring size to be independent
            // of the method

            // checks multi or peri
            if (MultiRing(subject, target))
                multi = true;
            else
                peri = true;

            for (unsigned int i=0; i<atom_intersection.size()-2;i++){ // peri atoms is always -2 the size coming from poly
                merged_ring.peri_chars.push_back( char(i + 'A') );
                merged_ring.peri_atoms.push_back(subject.atom_container.at(i));
            }

            binding_pos = 0;
            binding_atom_index = subject.atom_container.front().atom_index;
        } // complex peri fuse - all cases


        unsigned int start_index; // Keep the path definition
        // the position of the start atom.
        for (unsigned int i =0; i<AllRingAtoms.size();i++){
            if (AllRingAtoms.at(i).atom_index == subject.start_atom){
                start_index = i;
                break;
            }
        }

        // Build Connection Matrix from AllRingAtoms vector.
        Graph ring_graph(AllRingAtoms.size()); // creates graph of the right size;
        CreateGraph(AllRingAtoms, ring_graph); // Adds all the bonds in relation to position in AllRingAtoms vector
        ring_graph.GeneratePath();                      // Generates the Locant Paths into the graph object. // first one is usually the lowest index.
        ring_graph.ReorderPaths(start_index); // reorder based off the start index.

        // we'll come back to minimising once all paths are done properly.
        std::vector<std::vector<AtomContainer>> potential_paths;
        for (unsigned int i=0;i<ring_graph.indexes.size();i++){
            std::vector<AtomContainer> Ordered = VectorReorder(AllRingAtoms, ring_graph.indexes.at(i));
            if (Ordered.at(binding_pos).atom_index == binding_atom_index)
                potential_paths.push_back(Ordered);
        }

        merged_ring.atom_container = potential_paths.front();

        // type building
        if (poly){

            CloneRingData(merged_ring, subject, target, "POLYCYCLIC");
            merged_ring.poly_chars.push_back(binding_char);
            return merged_ring;}
        if (peri){
            // Allringatoms now contain the fully reordered set for the new symbol
            CloneRingData(merged_ring, subject, target, "PERICYCLIC");
            return merged_ring;}
        if (multi){
            CloneRingData(merged_ring, subject, target, "MULTICYCLIC");
            merged_ring.poly_chars.push_back(binding_char);
            return merged_ring;
        }

    }

    if (subject.ring_type == "PERICYCLIC"){

        if (atom_intersection.size() == 2){
            multi = true;
            for (int i=0;i<subject.atom_container.size();i++){
                if (subject.atom_container.at(i).atom_index == atom_intersection.front() ||
                    subject.atom_container.at(i).atom_index == atom_intersection.back()){
                    binding_char = i + 'A';
                    binding_pos = i;
                    binding_atom_index = subject.atom_container.at(i).atom_index; // always keep the path here.
                    merged_ring.poly_chars.push_back(binding_char);
                    break;
                }
            }
        }
        if (atom_intersection.size() == 3){
            // path must be consistent with the peri going first
            peri = true;
            binding_char = subject.peri_chars.size() + 'A';
            binding_pos = subject.peri_chars.size();
            merged_ring.peri_chars.push_back(binding_char);
            binding_atom_index = atom_intersection.front();
        }

        if (atom_intersection.size() > 4){
            peri = true;
            binding_char = subject.peri_chars.size() + 'A';
            binding_pos = subject.peri_chars.size();
            merged_ring.peri_chars.push_back(binding_char);
            binding_atom_index = atom_intersection.front();
        }

        AtomContainer start_atom = subject.peri_atoms.back();
        unsigned int start_index; // Keep the path definition
        // the position of the start atom.
        for (unsigned int i =0; i<AllRingAtoms.size();i++){
            if (AllRingAtoms.at(i).atom_index == start_atom.atom_index)
                start_index = i;
        }

        // Build Connection Matrix from AllRingAtoms vector.
        Graph ring_graph(AllRingAtoms.size()); // creates graph of the right size;
        CreateGraph(AllRingAtoms, ring_graph); // Adds all the bonds in relation to position in AllRingAtoms vector
        ring_graph.GeneratePath();                      // Generates the Locant Paths into the graph object. // first one is usually the lowest index.
        ring_graph.ReorderPaths(start_index); // reorder based off the start index.

        // we'll come back to minimising once all paths are done properly.
        std::vector<std::vector<AtomContainer>> potential_paths;
        for (unsigned int i=0;i<ring_graph.indexes.size();i++){
            std::vector<AtomContainer> Ordered = VectorReorder(AllRingAtoms, ring_graph.indexes.at(i));
            if (Ordered.at(binding_pos).atom_index == binding_atom_index)
                potential_paths.push_back(Ordered);
        }

        merged_ring.atom_container = potential_paths.front();

        // type building
        if (multi){
            CloneRingData(merged_ring, subject, target, "MULTICYCLIC");
            return merged_ring;}
        if (peri){
            // Allringatoms now contain the fully reordered set for the new symbol
            CloneRingData(merged_ring, subject, target, "PERICYCLIC");
            return merged_ring;}

    }

    if (subject.ring_type == "MULTICYCLIC"){

        if (atom_intersection.size() == 2){
            // standard poly cyclic binding to multi system
            multi = true;
            for (int i=0;i<subject.atom_container.size();i++){
                if (subject.atom_container.at(i).atom_index == atom_intersection.front() ||
                    subject.atom_container.at(i).atom_index == atom_intersection.back()){
                    binding_char = i + 'A';
                    binding_pos = i;
                    binding_atom_index = subject.atom_container.at(i).atom_index;
                    break;
                }
            }
        }
        if (atom_intersection.size() == 3){

            for (auto v: subject.peri_chars)
                std::cout << v << " ";
            std::cout << "\n";
            // path must be consistent with the peri going first
            // this is where the formula i wrote for peri check comes in handy.
            if (MultiRing(subject,target))
                multi = true;
            else
                peri = true;

        }
        if (atom_intersection.size() > 3){
            if (MultiRing(subject, target))
                multi = true;
            else
                peri = true;

            for (unsigned int i=0; i<atom_intersection.size()-2;i++){ // peri atoms is always -2 the size coming from poly
                merged_ring.peri_chars.push_back(i + 'A');
                merged_ring.peri_atoms.push_back(subject.atom_container.at(i));
            }

            binding_pos = 0;
            binding_atom_index = subject.atom_container.front().atom_index;

        }

        unsigned int start_index; // Keep the path definition
        // the position of the start atom.
        for (unsigned int i =0; i<AllRingAtoms.size();i++){
            if (AllRingAtoms.at(i).atom_index == subject.start_atom){
                start_index = i;
                break;
            }
        }

        // Build Connection Matrix from AllRingAtoms vector.
        Graph ring_graph(AllRingAtoms.size()); // creates graph of the right size;
        CreateGraph(AllRingAtoms, ring_graph); // Adds all the bonds in relation to position in AllRingAtoms vector
        ring_graph.GeneratePath();                      // Generates the Locant Paths into the graph object. // first one is usually the lowest index.
        ring_graph.ReorderPaths(start_index); // reorder based off the start index.

        // we'll come back to minimising once all paths are done properly.
        std::vector<std::vector<AtomContainer>> potential_paths;
        for (unsigned int i=0;i<ring_graph.indexes.size();i++){
            std::vector<AtomContainer> Ordered = VectorReorder(AllRingAtoms, ring_graph.indexes.at(i));
            if (Ordered.at(binding_pos).atom_index == binding_atom_index)
                potential_paths.push_back(Ordered);
        }

        merged_ring.atom_container = potential_paths.front();

        // type building
        if (multi){
            CloneRingData(merged_ring, subject, target, "MULTICYCLIC");
            merged_ring.poly_chars.push_back(binding_char);
            if (!subject.poly_chars.empty())
                merged_ring.poly_chars.insert(merged_ring.poly_chars.end(),subject.poly_chars.begin(),subject.poly_chars.end());
            return merged_ring;
        }
        if (peri){
            // Allringatoms now contain the fully reordered set for the new symbol
            CloneRingData(merged_ring, subject, target, "PERICYCLIC");
            return merged_ring;}

    }

}

bool WriteWLN::CreateGraph(std::vector<AtomContainer> atom_list, Graph &ring_graph){
    for (int i=0; i<atom_list.size(); i++){
        AtomContainer src = atom_list.at(i);
        for (int k=i+1;k<atom_list.size();k++){
            AtomContainer trg = atom_list.at(k);
            if(Joined(src, trg)){
                ring_graph.AddEdge(i,k);
            }
        }
    }

    return true;
}
unsigned int WriteWLN::ScorePaths(std::vector<AtomContainer> atom_list, Graph &ring_graph){
    // Check for the minimal path available based on bonding order.

    bool rings = false;
    bool groups = false;

    std::vector<int> earliest_ring_pos;
    std::vector<int> earliest_group_pos;
    for (unsigned int i=0;i<ring_graph.indexes.size();i++){
        std::vector<int> index_vec = ring_graph.indexes.at(i);
        std::vector<int> priority_check;
        for (unsigned int k=1;k< index_vec.size();k++){ // Miminising on the first term causes big issues!
            int index = index_vec.at(k);
            AtomContainer atom = atom_list.at(index);
            int priority = ExternalBond(atom, atom_list);
            if (priority == 2)
                rings = true;
            if (priority == 1)
                groups = true;

            priority_check.push_back(priority);
        }

        if (rings){
            int ring_index;
            auto itr=std::find(priority_check.begin(), priority_check.end(), 2);
            if(itr!=priority_check.end())
                ring_index =std::distance(priority_check.begin(), itr);

            earliest_ring_pos.push_back(ring_index);
        }

        if (groups){
            int group_index;
            auto itr=std::find(priority_check.begin(), priority_check.end(), 1);
            if(itr!=priority_check.end())
                group_index =std::distance(priority_check.begin(), itr);

            earliest_group_pos.push_back(group_index);
        }
    }


    if (rings){
        std::vector<int>::iterator it = std::min_element(std::begin(earliest_ring_pos), std::end(earliest_ring_pos));
        return std::distance(std::begin(earliest_ring_pos), it);
    }

    if (groups){
        std::vector<int>::iterator it = std::min_element(std::begin(earliest_group_pos), std::end(earliest_group_pos));
        return std::distance(std::begin(earliest_group_pos), it);
    }

}
unsigned int WriteWLN::ScoreLists(std::vector<AtomContainer> list_1, std::vector<AtomContainer> list_2){
    bool rings;
    bool groups;
    std::vector<int> score_A;
    for (unsigned int i=0;i<list_1.size();i++){
        AtomContainer atom = list_1.at(i);
        int priority = ExternalBond(atom, list_1);
        if (priority == 2)
            rings = true;
        if (priority == 1)
            groups = true;

        score_A.push_back(priority);
    }

    std::vector<int> score_B;
    for (unsigned int i=0;i<list_2.size();i++){
        AtomContainer atom = list_2.at(i);
        int priority = ExternalBond(atom, list_2);
        if (priority == 2)
            rings = true;
        if (priority == 1)
            groups = true;

        score_B.push_back(priority);
    }

    if (rings){
        int ring_index_A;
        auto itr_a=std::find(score_A.begin(), score_A.end(), 2);
        if(itr_a!=score_A.end())
            ring_index_A =std::distance(score_A.begin(), itr_a);

        int ring_index_B;
        auto itr_b=std::find(score_B.begin(), score_B.end(), 2);
        if(itr_b!=score_B.end())
            ring_index_B =std::distance(score_B.begin(), itr_b);

        if (ring_index_A<ring_index_B)
            return 0;
        else
            return 1;

    }

    if (groups){
        int group_index_A;
        auto itr_a=std::find(score_A.begin(), score_A.end(), 1);
        if(itr_a!=score_A.end())
            group_index_A =std::distance(score_A.begin(), itr_a);

        int group_index_B;
        auto itr_b=std::find(score_B.begin(), score_B.end(), 1);
        if(itr_b != score_B.end())
            group_index_B =std::distance(score_B.begin(), itr_b);

        if (group_index_A<group_index_B)
            return 0;
        else
            return 1;

    }


}
bool WriteWLN::CloneRingData(WLNSymbol &merged_ring, WLNSymbol &subject ,WLNSymbol &target, std::string type){

    // Add the non inter ring OUTbounds
    for (unsigned int i=0;i<target.outbound_connections.size();i++){
        BondContainer bond = target.outbound_connections.at(i);
        if (!bond.ring_bond)
            merged_ring.outbound_connections.push_back(bond);
    }

    for (unsigned int i=0;i<subject.outbound_connections.size();i++){
        BondContainer bond = subject.outbound_connections.at(i);
        if (!bond.ring_bond)
            merged_ring.outbound_connections.push_back(bond);
    }


    // Add the non inter ring INbounds
    for (unsigned int i=0;i<target.inbound_connections.size();i++){
        BondContainer bond = target.inbound_connections.at(i);
        if (!bond.ring_bond)
            merged_ring.inbound_connections.push_back(bond);
    }

    for (unsigned int i=0;i<subject.inbound_connections.size();i++){
        BondContainer bond = subject.inbound_connections.at(i);
        if (!bond.ring_bond)
            merged_ring.inbound_connections.push_back(bond);
    }

    merged_ring.ring_type = merged_ring.contained_ring.type = std::move(type);

    if (subject.SSRS.empty())
        merged_ring.SSRS.push_back(subject.contained_ring);
    else
        merged_ring.SSRS.insert(merged_ring.SSRS.end(), subject.SSRS.begin(), subject.SSRS.end());

    merged_ring.SSRS.push_back(target.contained_ring);
    merged_ring.start_atom = merged_ring.atom_container.front().atom_index;
    merged_ring.end_atom = merged_ring.atom_container.back().atom_index;
    merged_ring.contained_ring.size = merged_ring.atom_container.size();

    // Idiot for missing this really...
    merged_ring.contained_ring.ring_atoms = merged_ring.atom_container;

    merged_ring.contained_ring.ring_bonds.insert(merged_ring.contained_ring.ring_bonds.end(),
                                                 subject.contained_ring.ring_bonds.begin(),subject.contained_ring.ring_bonds.end());
    merged_ring.contained_ring.ring_bonds.insert(merged_ring.contained_ring.ring_bonds.end(),
                                                 target.contained_ring.ring_bonds.begin(),target.contained_ring.ring_bonds.end());


    CheckDuplicates(merged_ring.contained_ring.ring_bonds);
    CheckDuplicates(merged_ring.outbound_connections);
    CheckDuplicates(merged_ring.inbound_connections);
    return true;
}


// Working Function for executable testing
std::string TranslateToWLN(const char* smiles_string){
    OpenBabel::OBMol *Mol = new OpenBabel::OBMol;
    OpenBabel::OBConversion conv;
    conv.SetInFormat("smi");
    if(!conv.ReadString(Mol, smiles_string))
        return "NULL - SMILES READ ERROR";
    Mol->DeleteHydrogens();

    WriteWLN bond_engine{Mol};
    std::vector<char> wln_out;
    bond_engine.CreateString(wln_out);
    std::string res(begin(wln_out), end(wln_out));
    return res;
};


bool MBWriterWLN(OpenBabel::OBMol *mol, std::string &buffer){
    // converts wln to a given buffer
    mol->DeleteHydrogens();
    WriteWLN wr{mol};
    std::vector<char> wln_string;
    if (!wr.CreateString(wln_string))
        return false;

    // this should return the wln string
    std::string res(begin(wln_string), end(wln_string));
    buffer = res;
    return true;
}


