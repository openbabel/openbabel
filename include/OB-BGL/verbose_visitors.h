/**********************************************************************
verbose_visitors.h - visitor classes dumping a short text message for every visitor function.
 
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

#ifndef VERBOSE_VISITORS
#define VERBOSE_VISITORS

#include<ostream>

namespace OpenBabel {

class OBVerboseDFSVisitor {

    std::ostream& os;

    public:

    OBVerboseDFSVisitor( std::ostream& o ) : os(o) {}

    template <class Vertex, class Graph>
    void initialize_vertex(Vertex u, const Graph& g)
    {
        os << "initialize_vertex() called on:     " << u->GetIdx() << "\n";   
    }

    template <class Vertex, class Graph>
    void start_vertex(Vertex u, const Graph& g) 
    {
        os << "start_vertex() called on:          " << u->GetIdx() << "\n";   
    }
    template <class Vertex, class Graph>
    void discover_vertex(Vertex u, const Graph& g)
    {
        os << "discover_vertex() called on:       " << u->GetIdx() << "\n";   
    }
    template <class Edge, class Graph>
    void examine_edge(Edge u, const Graph& g)
    {
        os << "examine_edge() called on:          " << u->GetIdx() << "\n";   
    }
    template <class Edge, class Graph>
    void tree_edge(Edge u, const Graph& g)
    {
        os << "tree_edge() called on:             " << u->GetIdx() << "\n";   
    }
    template <class Edge, class Graph>
    void back_edge(Edge u, const Graph& g)
    {
        os << "back_edge() called on:             " << u->GetIdx() << "\n";   
    }
    template <class Edge, class Graph>
    void forward_or_cross_edge(Edge u, const Graph& g)
    {
        os << "forward_or_cross_edge() called on: " << u->GetIdx() << "\n";   
    }
    template <class Vertex, class Graph>
    void finish_vertex(Vertex u, const Graph& g)
    {
        os << "finish_vertex() called on:         " << u->GetIdx() << "\n";   
    }
};

class OBVerboseBFSVisitor
{
    std::ostream& os;
 
  public :

    OBVerboseBFSVisitor( std::ostream& o ) : os(o) {}

    template< typename vertex, typename graph >
    void initialize_vertex( vertex e, const graph& g)
    {
        os << "initialize_vertex() called on: " << e->GetIdx() << "\n" ;
    }
    template< typename vertex, typename graph >
    void discover_vertex( vertex e, const graph& g)
    {
        os << "discover_vertex() called on:   " << e->GetIdx() << "\n" ;
    }
    template< typename vertex, typename graph >
    void examine_vertex( vertex e, const graph& g)
    {
        os << "examine_vertex() called on:    " << e->GetIdx() << "\n" ;
    }
    template< typename edge, typename graph >
    void examine_edge( edge e, const graph& g)
    {
        os << "examine_edge() called on:      " << e->GetIdx() << "\n" ;
    }
    template< typename edge, typename graph >
    void tree_edge( edge e, const graph& g)
    {
        os << "tree_edge() called on:         " << e->GetIdx() << "\n" ;
    }
    template< typename edge, typename graph >
    void non_tree_edge( edge e, const graph& g)
    {
        os << "non_tree_edge() called on:     " << e->GetIdx() << "\n" ;
    }
    template< typename edge, typename graph >
    void cycle_edge( edge e, const graph& g)
    {
        os << "cycle_edge() called on:        " << e->GetIdx() << "\n" ;
    }
    template< typename edge, typename graph >
    void gray_target( edge e, const graph& g)
    {
        os << "gray_target() called on:       " << e->GetIdx() << "\n" ;
    }
    template< typename edge, typename graph >
    void black_target( edge e, const graph& g)
    {
        os << "black_target() called on:      " << e->GetIdx() << "\n" ;
    }
    template< typename vertex, typename graph >
    void finish_vertex( vertex e, const graph& g)
    {
        os << "finish_vertex() called on:     " << e->GetIdx() << "\n" ;
    }
};

} // end namespace

#endif // VERBOSE_VISITORS
