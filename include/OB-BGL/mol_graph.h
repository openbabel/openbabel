/**********************************************************************
mol_graph.h - interface with Boost Graph Library.
 
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

#ifndef MOL_GRAPH_H
#define MOL_GRAPH_H

#include <boost/graph/graph_traits.hpp>
#include <boost/shared_ptr.hpp> // OBAtomMap, OBBondMap
#include <boost/property_map.hpp> // OBAtomMap, OBBondMap
#include <boost/graph/visitors.hpp> // functor_caller, recorder

#include"openbabel/mol.h"

namespace boost // cannot specialize boost::graph_traits outside namespace
{
   //  There is a bunch of functions in the BGL modifying the structure 
   //  of the graph -- like inserting or deleting vertices or edges. None of 
   //  these will be implemented here.

  struct obmol_traversal_tag :

   //  For the purpose of the "tag classes" and the ideas behind design 
   //  and implementation of the Boost Graph Library see 
   //  The Boost Graph Library User and Reference Manual (TBGL) 
   //  by Siek, Lee and Lumsdaine or
   //  www.boost.org/libs/graph/doc/table_of_contents.html 

    public virtual incidence_graph_tag, 
    public virtual edge_list_graph_tag,
    public virtual vertex_list_graph_tag {};

  template<>
  struct graph_traits< OpenBabel::OBMol > 
  {

    typedef obmol_traversal_tag  traversal_category;
    typedef undirected_tag       directed_category;
    typedef disallow_parallel_edge_tag edge_parallel_category;
    typedef OpenBabel::OBAtom*   vertex_descriptor;

    typedef int   vertices_size_type;
    typedef int   edges_size_type;
    typedef int   degree_size_type;
    typedef OpenBabel::OBAtomIterator  vertex_iterator;
    typedef OpenBabel::OBBondIterator  edge_iterator;

    class edge_descriptor
    {
        OpenBabel::OBBond* bond;
        bool begin_is_source; // a dereferenced out_edge_iterator
                              // is required to know its source and target
                              // vertex, so we keep track in a bool variable 
      public:
        edge_descriptor( OpenBabel::OBBond* a = 0, bool b = true )
        :  bond(a), begin_is_source(b) {}
        OpenBabel::OBBond* operator->() { return bond; }
        OpenBabel::OBBond& operator*()  { return *bond; }
        vertex_descriptor source()     
        { 
            if ( begin_is_source ) return bond->GetBeginAtom();
            return bond->GetEndAtom();
        }
        vertex_descriptor target()     
        { 
            if ( begin_is_source ) return bond->GetEndAtom();
            return bond->GetBeginAtom();
        }
        // this function is to hide the edge_descriptor class
        // from users of the obmol BGL interface
        operator OpenBabel::OBBond*()  { return bond; } 
    };

    class adjacent_iterator
    {
        OpenBabel::OBBondIterator ip_bond;
        OpenBabel::OBAtom* p_atom;
      public :
        typedef std::forward_iterator_tag iterator_category;
        typedef OpenBabel::OBAtom*  value_type;
        typedef OpenBabel::OBAtom** pointer;
        typedef OpenBabel::OBBondIterator::difference_type difference_type;
        typedef OpenBabel::OBBondIterator::reference       reference;

        adjacent_iterator () {};
        adjacent_iterator( OpenBabel::OBBondIterator i, OpenBabel::OBAtom* a ) : 
            ip_bond(i), p_atom(a) {}
        adjacent_iterator( const adjacent_iterator& i ) : 
            ip_bond(i.ip_bond), p_atom(i.p_atom) {}
        adjacent_iterator& operator++() { ++ip_bond; return *this; }
        adjacent_iterator  operator++(int) 
               { adjacent_iterator tmp(ip_bond,p_atom); ++ip_bond; return tmp; }
        bool operator==( const adjacent_iterator& rhs ) 
            { return ip_bond == rhs.ip_bond; }
        bool operator!=( const adjacent_iterator& rhs ) { return !(*this == rhs); }
        vertex_descriptor operator*() { return (*ip_bond)->GetNbrAtom(p_atom); }
    };

    class out_edge_iterator
    {
        OpenBabel::OBBondIterator ip_bond;
        OpenBabel::OBAtom* p_atom; // source OBAtom*
      public :
        typedef std::forward_iterator_tag iterator_category;
        typedef edge_descriptor value_type;
        typedef edge_descriptor* pointer;
        typedef OpenBabel::OBBondIterator::difference_type difference_type;
        typedef OpenBabel::OBBondIterator::reference       reference;

        out_edge_iterator () {};
        out_edge_iterator( const out_edge_iterator& i ) : 
            ip_bond(i.ip_bond), p_atom(i.p_atom) {}
        out_edge_iterator( OpenBabel::OBBondIterator i, OpenBabel::OBAtom* a ) : 
            ip_bond(i), p_atom(a) {}
        out_edge_iterator& operator++() { ++ip_bond; return *this; }
        out_edge_iterator  operator++(int) 
                   { out_edge_iterator tmp(ip_bond,p_atom); ++ip_bond; return tmp; }
        bool operator==( const out_edge_iterator& rhs ) 
            { return ip_bond == rhs.ip_bond; }
        bool operator!=( const out_edge_iterator& rhs ) { return !(*this == rhs); }
        edge_descriptor operator*() 
        {
            // this is why we need to keep track of the source atom 
            return edge_descriptor(*ip_bond,(*ip_bond)->GetBeginAtom() == p_atom);
        }
    };    

  };

  std::pair < graph_traits< OpenBabel::OBMol >::edge_iterator,
              graph_traits< OpenBabel::OBMol >::edge_iterator >
  edges( const OpenBabel::OBMol& m )
  { 
      // Now that cast is ugly: but m.BeginBonds() violates
      // the const qualifier of the argument
      // and this is the signature imposed by BGL.
      // A better solution would be very welcome
      // propably a ConstOBBondIterator
      OpenBabel::OBMol* mol = const_cast< OpenBabel::OBMol* >(&m);
      return std::make_pair(mol->BeginBonds(),mol->EndBonds());
  }

  std::pair < graph_traits< OpenBabel::OBMol >::vertex_iterator,
              graph_traits< OpenBabel::OBMol >::vertex_iterator >
  vertices( const OpenBabel::OBMol& m )
  { 
      // Now that cast is ugly: but m.BeginAtoms() violates
      // the const qualifier of the argument
      // and this is the signature imposed by BGL.
      // A better solution would be very welcome
      // propably a ConstOBAtomIterator
      OpenBabel::OBMol* mol = const_cast< OpenBabel::OBMol* >(&m);
      return std::make_pair(mol->BeginAtoms(),mol->EndAtoms()); 
  }

  std::pair < graph_traits< OpenBabel::OBMol >::adjacent_iterator,
              graph_traits< OpenBabel::OBMol >::adjacent_iterator >
  adjacent_vertices ( graph_traits< OpenBabel::OBMol >::vertex_descriptor& v,
                      const OpenBabel::OBMol& )
  {
      return std::make_pair(
      graph_traits< OpenBabel::OBMol >::adjacent_iterator(v->BeginBonds(),v),
      graph_traits< OpenBabel::OBMol >::adjacent_iterator(v->EndBonds()  ,v) );
  }

  graph_traits< OpenBabel::OBMol >::vertex_descriptor
  source( graph_traits< OpenBabel::OBMol >::edge_descriptor& e,
          const OpenBabel::OBMol& )
  {
      return e.source();
  }

  graph_traits< OpenBabel::OBMol >::vertex_descriptor
  target( graph_traits< OpenBabel::OBMol >::edge_descriptor e,
          const OpenBabel::OBMol& )
  {
      return e.target();
  }

  std::pair < graph_traits< OpenBabel::OBMol >::out_edge_iterator,
              graph_traits< OpenBabel::OBMol >::out_edge_iterator >
  out_edges ( graph_traits< OpenBabel::OBMol >::vertex_descriptor v,
              const OpenBabel::OBMol& )
  {
      return std::make_pair(
      graph_traits< OpenBabel::OBMol >::out_edge_iterator(v->BeginBonds(),v),
      graph_traits< OpenBabel::OBMol >::out_edge_iterator(v->EndBonds()  ,v) );
  }

  int out_degree( graph_traits< OpenBabel::OBMol >::vertex_descriptor v,
                  const OpenBabel::OBMol& )
  { 
      return v->EndBonds() - v->BeginBonds(); 
  }

  int num_vertices( const OpenBabel::OBMol& m )
  { 
      return m.NumAtoms(); 
  }

  int num_edges( const OpenBabel::OBMol& m )
  { 
      return m.NumBonds(); 
  }

} // namespace

// When putting these classes inside the OpenBabel namespace 
// the put and get functions will not work with the BGL algorithms.
// Don't know why.

template< class value >
class OBAtomMap 
{
  public :

    typedef value value_type;
    typedef OpenBabel::OBAtom* key_type;
    typedef value_type& reference;
    typedef boost::read_write_property_map_tag category;
 
    OBAtomMap( const OpenBabel::OBMol& m ) :
        _values( new std::vector< value >( m.NumAtoms() ) )
    {}
    OBAtomMap( const OpenBabel::OBMol& m, const value v ) :
        _values( new std::vector< value >( m.NumAtoms(),v ) )
    {}

    reference operator[]( OpenBabel::OBAtom* a ) 
        { return (*_values)[ a->GetIdx()-1 ]; }

  private :

    boost::shared_ptr< std::vector< value_type > > _values;

};

template< class value >
class OBBondMap 
{
  public :

    typedef value value_type;
    typedef OpenBabel::OBBond* key_type;
    typedef value_type& reference;
    typedef boost::read_write_property_map_tag category;
 
    OBBondMap( const OpenBabel::OBMol& m ) :
        _values( new std::vector< value >( m.NumBonds() ) )
    {}
    OBBondMap( const OpenBabel::OBMol& m, const value v ) :
        _values( new std::vector< value >( m.NumBonds(),v ) )
    {}

    reference operator[]( OpenBabel::OBBond* b ) 
        { return (*_values)[ b->GetIdx() ]; }

  private :

    boost::shared_ptr< std::vector< value_type > > _values;
};

template< class map_t >
void put( map_t& map, typename map_t::key_type key,  typename map_t::value_type value )
    { map[key] = value; }
 
template< class map_t >
typename map_t::value_type get( map_t& map, typename map_t::key_type key )
    { return map[key]; }

namespace boost {

template <class functor, class Tag >
struct functor_caller
    : public base_visitor< functor_caller < functor, Tag> >
{
    typedef Tag event_filter;
    functor_caller(functor f) : m_func(f) { }
    template <class T, class Graph>
    void operator()(T t, const Graph& g) {
      m_func.operator()(t);
    }
    functor m_func;
};

template < class functor, class Tag >
functor_caller< functor, Tag >
call_functor(functor func, Tag) {
    return functor_caller< functor, Tag > (func);
}

template <class output_iterator, class Tag >
struct recorder
    : public base_visitor< recorder < output_iterator, Tag> >
{
    typedef Tag event_filter;
    recorder(output_iterator o) : m_output(o) { }
    template <class T, class Graph>
    void operator()(T t, const Graph& g) {
      m_output = t;
    }
    output_iterator m_output;
};

template <class output_iterator, class Tag >
recorder< output_iterator, Tag > record(output_iterator out, Tag) {
    return recorder< output_iterator, Tag > (out);
}

} // namespace

#endif // MOL_GRAPH_H
