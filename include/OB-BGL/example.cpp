/**********************************************************************
example.cpp - example of interfacing of OpenBabel with Boost Graph Library
 
Copyright (C) 2007 by Gerde Menche
 
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

#include "babelconfig.h"
#include "openbabel/mol.h"
#include "openbabel/obconversion.h"

#include <iostream>
//#include <unistd.h>

#include "mol_graph.h"
#include "verbose_visitors.h"

#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/undirected_dfs.hpp>
#include <boost/pending/queue.hpp>

using namespace std;
using OpenBabel::OBMol;
using OpenBabel::OBAtom;
using OpenBabel::OBBond;
using OpenBabel::OBConversion;
using OpenBabel::OBFormat;

int main(int argc,char **argv)
{
    int c;
    char *FileIn = NULL;

    if (argc != 2)
    {
        cout << "give file name\n";
        return 1;
    }

    FileIn  = argv[1];

    // Find Input filetype
    OBConversion conv;
    OBFormat *format = conv.FormatFromExt(FileIn);
    
    if (!format || !conv.SetInAndOutFormats(format, format))
    {
        cerr << argv[0] << ": cannot read input format!" << endl;
        exit (-1);
    }

    ifstream ifs;

    // Read the file
    ifs.open(FileIn);
    if (!ifs)
    {
        cerr << argv[0] << ": cannot read input file!" << endl;
        exit (-1);
    }

    OBMol mol;

    for (c=0;;)
    {
        mol.Clear();
	conv.Read(&mol, &ifs);
        if (mol.Empty())
            break;
        OBAtom* atom = *(mol.BeginAtoms()); 

        // Examples using breadth first search

        std::cout << "\nBREADTH FIRST SEARCH\n";

        boost::queue< OBAtom* > buffer;
        boost::breadth_first_search( mol,atom,buffer, OpenBabel::OBVerboseBFSVisitor(cout),
            OBAtomMap< boost::default_color_type >(mol) );

        OBAtomMap< OBAtom* >predecessor_map_bfs(mol,(OBAtom*)0);
        OBAtomMap< int >distance_map_bfs(mol, 0);

        boost::breadth_first_search( mol,atom,buffer,
            boost::make_bfs_visitor( std::make_pair ( 
                boost::record_predecessors(predecessor_map_bfs,boost::on_tree_edge()),
                boost::record_distances(distance_map_bfs,boost::on_tree_edge()) )),
            OBAtomMap< boost::default_color_type >(mol) );

        for ( OpenBabel::OBAtomIterator beg = mol.BeginAtoms(); 
             beg != mol.EndAtoms(); ++beg )
        {
            OBAtom* parent = predecessor_map_bfs[*beg]; 
            if ( parent )
            {
                std::cout << (*beg)->GetIdx() - 1 << " has parent : " 
                          << parent->GetIdx() - 1 << "\t";
            }
            else
            {
                std::cout << (*beg)->GetIdx() - 1 << " has no parent\t\t";
            }
            std::cout << " graph theoretical distance : " << distance_map_bfs[*beg] << "\n"; 
        }

        std::vector< boost::graph_traits<OpenBabel::OBMol>::edge_descriptor > bv;
        boost::breadth_first_search( mol,atom,buffer,
            boost::make_bfs_visitor( 
                boost::record( std::back_inserter(bv), boost::on_gray_target() )),
            OBAtomMap< boost::default_color_type >(mol) );

        std::cout << "\nBACK EDGES\n";
        for ( std::vector< boost::graph_traits<OpenBabel::OBMol>::edge_descriptor >::iterator 
              beg = bv.begin(), end = bv.end();
              beg != end; ++beg )
        {
            std::cout << '(' << (*beg)->GetBeginAtom()->GetIdx() << ")--";
            std::cout << (*beg)->GetIdx() << "--(";
            std::cout << (*beg)->GetEndAtom()->GetIdx() << ")\n";
        }

        // Examples using depth first search

        std::cout << "\nDEPTH FIRST SEARCH\n";

        boost::undirected_dfs( mol, OpenBabel::OBVerboseDFSVisitor(cout),
            OBAtomMap< boost::default_color_type >(mol),
            OBBondMap< boost::default_color_type >(mol),
            atom );

        OBAtomMap< OBAtom* >predecessor_map_dfs(mol,(OBAtom*)0);
        OBAtomMap< int >distance_map_dfs(mol,0);
            
        boost::undirected_dfs( mol, 
            boost::make_dfs_visitor( std::make_pair ( 
                boost::record_predecessors(predecessor_map_dfs,boost::on_tree_edge()),
                boost::record_distances(distance_map_dfs,boost::on_tree_edge()) )),
            OBAtomMap< boost::default_color_type >(mol),
            OBBondMap< boost::default_color_type >(mol),
            atom );

        for ( OpenBabel::OBAtomIterator beg = mol.BeginAtoms(); 
              beg != mol.EndAtoms(); ++beg )
        {
            OBAtom* parent = predecessor_map_dfs[*beg]; 
            if ( parent )
            {
                std::cout << (*beg)->GetIdx() - 1 << " has parent : " 
                          << parent->GetIdx() - 1 << "\t";
            }
            else
            {
                std::cout << (*beg)->GetIdx() - 1 << " has no parent\t\t";
            }
            std::cout << " graph theoretical distance : " << distance_map_dfs[*beg] << "\n"; 
        }

    } // end for loop

    return(1);
}
