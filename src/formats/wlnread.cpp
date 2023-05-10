/**********************************************************************
 
Author : Michael Blakey

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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <set>
#include <deque>
#include <vector>
#include <stack>
#include <map>

#include <utility> // std::pair
#include <iterator>
#include <sstream>

#include <openbabel/mol.h>
#include <openbabel/plugin.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include <openbabel/kekulize.h>
#include <openbabel/ring.h>
#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>


#define REASONABLE 1024

const char *cli_inp;

// --- options ---
static bool opt_wln2dot = false;
static bool opt_debug = false;


const char *wln_string;
struct WLNSymbol;
struct WLNEdge; 
struct WLNRing;
struct WLNGraph;
struct ObjectStack;


enum WLNTYPE
{
  STANDARD = 0,
  RING = 1,     
  SPECIAL = 2  // for now this is only going to be the pi bond
};

unsigned char static int_to_locant(unsigned int i){
  return i + 64;
}

unsigned int static locant_to_int(unsigned char loc){
  return loc - 64;
}


std::string get_notation(unsigned int s, unsigned int e)
{
  std::string res; 
  for (unsigned int i = s; i <= e; i++)
  {
    res.push_back(wln_string[i]);
  }
  return res; 
}

void Fatal(unsigned int pos)
{
  fprintf(stderr, "Fatal: %s\n", wln_string);
  fprintf(stderr, "       ");
  for (unsigned int i = 0; i < pos; i++)
    fprintf(stderr, " ");

  fprintf(stderr, "^\n");

  exit(1);
}


/**********************************************************************
                          STRUCT DEFINTIONS
**********************************************************************/
 

struct WLNEdge{
  WLNSymbol *parent;
  WLNSymbol *child;
  WLNEdge *nxt;

  bool aromatic;
  unsigned int order;

  WLNEdge(){
    parent = 0;
    child = 0;
    aromatic = 0;
    order = 0;
    nxt = 0;
  }
  ~WLNEdge(){};
};


struct WLNSymbol
{
  unsigned char ch;
  std::string special; // string for element, or ring, if value = '*'
  
  unsigned int type;
  unsigned int allowed_edges;
  unsigned int num_edges;

  WLNSymbol *previous;
  WLNEdge   *bonds; // array of bonds

  // if default needed
  WLNSymbol()
  {
    ch = '\0';
    allowed_edges = 0;
    num_edges = 0;
    type = 0;
    previous = 0;
    bonds = 0;
  }
  ~WLNSymbol(){};

  void set_edge_and_type(unsigned int e, unsigned int t=STANDARD){
    allowed_edges = e;
    type = t;
  }

  void add_special(unsigned int s, unsigned int e)
  {
    for (unsigned int i = s; i <= e; i++)
      special.push_back(wln_string[i]);
  }

};

struct WLNRing
{
  std::vector<unsigned int> rings;
  std::map<unsigned char, WLNSymbol *> locants; 
  std::map<WLNSymbol*,unsigned char> locants_ch;
  std::vector<std::pair<unsigned char,int>> post_charges; 
  
  WLNRing(){}
  ~WLNRing(){};
};


// handles all memory and 'global' vars
struct WLNGraph
{
  
  WLNSymbol *root;

  unsigned int edge_count;
  unsigned int symbol_count;
  unsigned int ring_count;

  WLNSymbol *SYMBOLS[REASONABLE];
  WLNEdge   *EDGES  [REASONABLE];
  WLNRing   *RINGS  [REASONABLE];

  std::map<WLNSymbol *, unsigned int> index_lookup;
  std::map<unsigned int, WLNSymbol *> symbol_lookup;

  unsigned int glob_index; // babel starts from 1, keep consistent  

    // ionic parsing
  std::map<unsigned int,WLNSymbol*> string_positions; 
  std::map<WLNSymbol*,int> charge_additions;

  WLNGraph(){
    edge_count   = 0;
    symbol_count = 0;
    ring_count   = 0;
    glob_index   = 1; // babel atoms are +1 indexed

    // pointer safety
    root = 0;
    for (unsigned int i = 0; i < REASONABLE;i++){
      SYMBOLS[i] = 0;
      EDGES[i] = 0;
      RINGS[i] = 0;
    }
  };

  ~WLNGraph(){
    for (unsigned int i = 0; i < REASONABLE;i++){
      delete SYMBOLS[i];
      delete EDGES[i];
      delete RINGS[i];
    }
  }
};

// some bridge notation is index dependent
struct indexed_pair{
  unsigned char bind_1 = '\0';
  unsigned char bind_2 = '\0';
  unsigned int index   = 0;

  void set(unsigned char a, unsigned char b, unsigned int p){
    bind_1 = a;
    bind_2 = b;
    index = p;
  }

};

// needs to be able to hold both a WLNSymbol and WLNRing for branch returns
struct ObjectStack{  
  std::vector<std::pair<WLNRing*,WLNSymbol*>> stack; // vector so i can iterate down and instant access
  WLNRing   *ring;
  WLNSymbol *branch;
  unsigned int size;

  ObjectStack(){
    ring = 0;
    branch = 0;
    size = 0;
  }

  void reserve(unsigned int n){
    stack.reserve(n);
  }

  bool peek(){
    if(!size){
      fprintf(stderr,"Error: peeking empty ring stack\n");
      return false;
    }
    else{
      fprintf(stderr,"top: ring: %p   branch: %p\n",stack.back().first,stack.back().second);
      return true;
    }
     
  }

  bool pop(){
    stack.pop_back();
    size--;

    ring = 0;
    branch = 0;

    if(stack.empty()){
      fprintf(stderr,"Error: popping empty ring stack\n");
      return false;
    }
     
    for (int i=size-1;i>-1;i--){
      
      if(!ring && stack[i].first)
        ring = stack[i].first;
      
      if(!branch && stack[i].second)
        branch = stack[i].second;        
    }

    return true; 
  }

  void push(std::pair<WLNRing*,WLNSymbol*> pair,bool verbose = false){
    stack.push_back(pair);
    if(pair.first)
      ring = pair.first;

    if(pair.second)
      branch = pair.second;

    if(verbose){
      fprintf(stderr,"pushed: ring: %p    branch: %p\n",pair.first,pair.second);
    }

    size++;
  }


  bool empty(){
    if (stack.empty())
      return true;
    else 
      return false;
  }

  void clear_all(){
    ring = 0;
    branch = 0;
    stack.clear();
  }

  std::pair<WLNRing*,WLNSymbol*> & top(){
    return stack.back();
  }

  // cleans branches
  bool branch_avaliable(){
    if(branch && branch->num_edges < branch->allowed_edges)
      return true;
    else
      return false;
  }

};



/**********************************************************************
                         WLNSYMBOL Functions
**********************************************************************/


WLNSymbol *AllocateWLNSymbol(unsigned char ch, WLNGraph &graph)
{

  graph.symbol_count++;
  if(graph.symbol_count > REASONABLE){
    fprintf(stderr,"Error: creating more than 1024 wln symbols - is this reasonable?\n");
    return 0;
  }

  if(!ch){
    fprintf(stderr,"Error: null char used to symbol creation\n");
    return 0;
  }

  WLNSymbol *wln = new WLNSymbol;
  graph.SYMBOLS[graph.symbol_count] = wln;
  
  wln->ch = ch;
  graph.index_lookup[wln] = graph.glob_index;
  graph.symbol_lookup[graph.glob_index] = wln;
  graph.glob_index++;
  return wln;
}

WLNSymbol* define_hypervalent_element(unsigned char sym, WLNGraph &graph){

  if(!sym){
    fprintf(stderr,"Error: null char used for hypervalent element allocation\n");
    return 0;
  }

  WLNSymbol *new_symbol = 0;
  
  switch(sym){
    
    case 'P':
    case 'S':
    case 'G':
    case 'E':
    case 'I':
    case 'F':
      new_symbol = AllocateWLNSymbol(sym,graph);
      new_symbol->set_edge_and_type(6);            // allows FCl6
      break;

    default:
      fprintf(stderr,"Error: character %c does not need - notation for valence expansion, please remove -\n",sym);
      break;
  }
  
  return new_symbol;
}

/* allocate new or override exisiting node*/
WLNSymbol* define_element(std::string special, WLNGraph &graph){
    
  WLNSymbol *created_wln = 0;
  
  switch (special[0]){

    case 'A':
      switch(special[1]){
        case 'C':
        case 'G':
        case 'L':
        case 'M':
        case 'R':
        case 'S':
        case 'T':
        case 'U':
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }
      break;

    case 'B':
      switch(special[1]){
        case 'A':
        case 'E':
        case 'H':
        case 'I':
        case 'K':
        case 'R':
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }
      break;
      

    case 'C':
      switch(special[1]){
        case 'A':
        case 'D':
        case 'E':
        case 'F':
        case 'M':
        case 'N':
        case 'O':
        case 'R':
        case 'S':
        case 'U':
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }
      break;
      
    case 'D':
      switch(special[1]){
        case 'B':
        case 'S':
        case 'Y':
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }
      break;

    case 'E':
      switch(special[1]){
        case 'R':
        case 'S':
        case 'U':
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }
      break;

    case 'F':
      switch(special[1]){
        case 'E':
        case 'L':
        case 'M':
        case 'R':
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }
      break;

    case 'G':
      switch(special[1]){
        case 'A':
        case 'D':
        case 'E':
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }
      break;

    case 'H':
      switch(special[1]){
        case 'E':
        case 'F':
        case 'G':
        case 'O':
        case 'S':
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }
      break;

    case 'I':
      switch(special[1]){
        case 'N':
        case 'R':
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }
      break;

    case 'K':
      switch(special[1]){
        case 'R':
        case 'A':
          created_wln = AllocateWLNSymbol('*',graph);
          break;

        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }
      break;
      

    case 'L':
      switch(special[1]){
        case 'A':
        case 'I':
        case 'R':
        case 'U':
        case 'V':
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }
      break;

    case 'M':
      switch(special[1]){
        case 'C':
        case 'D':
        case 'G':
        case 'N':
        case 'O':
        case 'T':
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }
      break;

    case 'N':
      switch(special[1]){
        case 'A':
        case 'B':
        case 'D':
        case 'E':
        case 'H':
        case 'I':
        case 'O':
        case 'P':
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }
      break;


    case 'O':
      switch(special[1]){
        case 'O':
        case 'G':
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }
      break;

    case 'P':
      switch(special[1]){
        case 'A':
        case 'B':
        case 'D':
        case 'M':
        case 'O':
        case 'R':
        case 'T':
        case 'U':
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }
      break;

    case 'R':
      switch(special[1]){
        case 'A':
        case 'B':
        case 'E':
        case 'F':
        case 'G':
        case 'H':
        case 'N':
        case 'U':
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }
      break;
     

    case 'S':
      switch(special[1]){
        case 'B':
        case 'C':
        case 'E':
        case 'G':
        case 'I':
        case 'M':
        case 'N':
        case 'R':
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }
      break;


    case 'T':
      switch(special[1]){
        case 'A':
        case 'B':
        case 'C':
        case 'E':
        case 'H':
        case 'I':
        case 'L':
        case 'M':
        case 'S':
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }
      break;

    case 'U':
      if(special[1] == 'R')
        created_wln = AllocateWLNSymbol('*',graph);
      else
      {
        fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
        return (WLNSymbol *)0;
      }
      break;

    case 'V':
      if (special[1] == 'A')
        created_wln = AllocateWLNSymbol('*',graph);
      else
      {
        fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
        return (WLNSymbol *)0;
      }
      break;
    
    case 'W':
      if(special[1] == 'T')
        created_wln = AllocateWLNSymbol('*',graph);
      else
      {
        fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
        return (WLNSymbol *)0;
      }
      break;
    

    case 'X':
      if (special[1] == 'E')
        created_wln = AllocateWLNSymbol('*',graph);
      else
      {
        fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
        return (WLNSymbol *)0;
      }
      break;

    case 'Y':
      switch(special[1]){
        case 'B':
        case 'T':
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }
      break;

    case 'Z':
      switch(special[1]){
        case 'N':
        case 'R':
          created_wln = AllocateWLNSymbol('*',graph);
          break;
           
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }
      break;

    default:
      fprintf(stderr, "Error: invalid character in special definition switch\n");
      return (WLNSymbol *)0;
  }

  created_wln->special = special;
  created_wln->allowed_edges = 8; // allow anything for now;
  return created_wln;
}


/* checks are already made, this should just return*/
unsigned int special_element_atm(std::string &special){
  
  switch (special[0]){

    case 'A':
      if (special[1] == 'C')
        return 89;
      else if (special[1] == 'G')
        return 47;
      else if (special[1] == 'L')
        return 13;
      else if (special[1] == 'M')
        return 95;
      else if (special[1] == 'R')
        return 18;
      else if (special[1] == 'S')
        return 33;
      else if (special[1] == 'T')
        return 85;
      else if (special[1] == 'U')
        return 79;
      break;

    case 'B':
      if (special[1] == 'A')
        return 56;
      else if (special[1] == 'E')
        return 4;
      else if (special[1] == 'H')
        return 107;
      else if (special[1] == 'I')
        return 83;
      else if (special[1] == 'K')
        return 97;
      else if (special[1] == 'R')
        return 35;
      break;

    case 'C':
      if (special[1] == 'A')
        return 20;
      else if (special[1] == 'D')
        return 48;
      else if (special[1] == 'E')
        return 58;
      else if (special[1] == 'F')
        return 98;
      else if (special[1] == 'M')
        return 96;
      else if (special[1] == 'N')
        return 112;
      else if (special[1] == 'O')
        return 27;
      else if (special[1] == 'R')
        return 24;
      else if (special[1] == 'S')
        return 55;
      else if (special[1] == 'U')
        return 29;
      break;

    case 'D':
      if (special[1] == 'B')
        return 105;
      else if (special[1] == 'S')
        return 110;
      else if (special[1] == 'Y')
        return 66;
      break;

    case 'E':
      if (special[1] == 'R')
        return 68;
      else if (special[1] == 'S')
        return 99; 
      else if (special[1] == 'U')
        return 63;
      break;

    case 'F':
      if (special[1] == 'E')
        return 26;
      else if (special[1] == 'L')
        return 114;
      else if (special[1] == 'M')
        return 100;
      else if (special[1] == 'R')
        return 97;
      break;

    case 'G':
      if (special[1] == 'A')
        return 31;
      else if (special[1] == 'D')
        return 64;
      else if (special[1] == 'E')
        return 32;
      break;

    case 'H':
      if (special[1] == 'E')
        return 2;
      else if (special[1] == 'F')
        return 72;
      else if (special[1] == 'G')
        return 80;
      else if (special[1] == 'O')
        return 67;
      else if (special[1] == 'S')
        return 108;

      break;

    case 'I':
      if (special[1] == 'N')
        return 49;
      else if (special[1] == 'R')
        return 77;
      break;

    case 'K':
      if (special[1] == 'R')
        return 39;
      else if(special[1] == 'A')
        return 19;
      break;

    case 'L':
      if (special[1] == 'A')
        return 57;
      else if (special[1] == 'I')
        return 3;
      else if (special[1] == 'R')
        return 103;
      else if (special[1] == 'U')
        return 71;
      else if (special[1] == 'V')
        return 116;
      break;

    case 'M':
      if (special[1] == 'C')
        return 115;
      else if (special[1] == 'D')
        return 101;
      else if (special[1] == 'G')
        return 12;
      else if (special[1] == 'N')
        return 25;
      else if (special[1] == 'O')
        return 42;
      else if (special[1] == 'T')
        return 109;
      break;

    case 'N':
      if (special[1] == 'A')
       return 11;
      else if (special[1] == 'B')
        return 41;
      else if (special[1] == 'D')
        return 60;
      else if (special[1] == 'E')
        return 10;
      else if (special[1] == 'H')
        return 113;
      else if (special[1] == 'I')
        return 28;
      else if (special[1] == 'O')
        return 102;
      else if (special[1] == 'P')
        return 93;
      break;

    case 'O':
      if (special[1] == 'G')
        return 118;
      else if (special[1] == 'S')
        return 76;
      break;

    case 'P':
      if (special[1] == 'A')
        return 91;       
      else if (special[1] == 'B')
        return 82;
      else if (special[1] == 'D')
        return 46;
      else if (special[1] == 'M')
        return 61;
      else if (special[1] == 'O')
        return 84;
      else if (special[1] == 'R')
        return 59;
      else if (special[1] == 'T')
        return 78;
      else if (special[1] == 'U')
        return 94;
      
      break;

    case 'R':
      if (special[1] == 'A')
        return 88;
      else if (special[1] == 'B')
        return 37;
      else if (special[1] == 'E')
        return 75;
      else if (special[1] == 'F')
        return 104;
      else if (special[1] == 'G')
        return 111;
      else if (special[1] == 'H')
        return 45;
      else if (special[1] == 'N')
        return 86;
      else if (special[1] == 'U')
        return 44;
      break;

    case 'S':
      if (special[1] == 'B')
        return 51;
      else if (special[1] == 'C')
        return 21;
      else if (special[1] == 'E')
        return 34;
      else if (special[1] == 'G')
        return 106;
      else if (special[1] == 'I')
        return 14;
      else if (special[1] == 'M')
        return 62;
      else if (special[1] == 'N')
        return 50;
      else if (special[1] == 'R')
        return 38;
      
      break;

    case 'T':
      if (special[1] == 'A')
        return 73;
      else if (special[1] == 'B')
        return 65;
      else if (special[1] == 'C')
        return 43;
      else if (special[1] == 'E')
        return 52;
      else if (special[1] == 'H')
        return 90;
      else if (special[1] == 'I')
        return 22;
      else if (special[1] == 'L')
        return 81;
      else if (special[1] == 'M')
        return 69;
      else if (special[1] == 'S')
        return 117;

      break;

    case 'U':
      if(special[1] == 'R')
        return 92;
      break;

    case 'V':
      if(special[1] == 'A')
        return 23;
      break;

    case 'X':
      if (special[1] == 'E')
        return 54;
      break;

    case 'Y':
      if(special[1] == 'T')
        return 39;
      else if (special[1] == 'B')
        return 70;
      break;

    case 'Z':
      if (special[1] == 'N')
        return 30;
      else if (special[1] == 'R')
        return 40;
  
      break;

    default:
      fprintf(stderr, "Error: invalid character in special definition switch\n");
      return 0;
  }

  return 0;
}


// this one pops based on bond numbers
WLNSymbol *return_object_symbol(ObjectStack &branch_stack){
    // only used for characters that can 'act' like a '&' symbol

  WLNSymbol *top = 0;
  while(!branch_stack.empty()){
    top = branch_stack.top().second;
    if(!top)
      return top; // only iterate to the next
    else if(top->num_edges == top->allowed_edges)
      branch_stack.pop();
    else
      return top;
  }

  return top;
}

/**********************************************************************
                          WLNEdge Functions
**********************************************************************/


WLNEdge *AllocateWLNEdge(WLNSymbol *child, WLNSymbol *parent,WLNGraph &graph){

  if(!child || !parent){
    fprintf(stderr,"Error: attempting bond of non-existent symbols - %s|%s is dead\n",child ? "":"child",parent ? "":"parent");
    return 0;
  }

  graph.edge_count++;
  if(graph.edge_count > REASONABLE){
    fprintf(stderr,"Error: creating more than 1024 wln symbols - is this reasonable?\n");
    return 0;
  }
  
  if ((child->num_edges + 1) > child->allowed_edges){
    fprintf(stderr, "Error: wln character[%c] is exceeding allowed connections %d/%d\n", child->ch,child->num_edges+1, child->allowed_edges);
    return 0;
  }
  
  if ((parent->num_edges + 1) > parent->allowed_edges){
    fprintf(stderr, "Error: wln character[%c] is exceeding allowed connections %d/%d\n", parent->ch,parent->num_edges+1, parent->allowed_edges);
    return 0;
  }

  WLNEdge *edge = new WLNEdge;
  graph.EDGES[graph.edge_count] = edge;

  // use a linked list to store the bond, can also check if it already exists

  WLNEdge *curr = parent->bonds;
  if(curr){
    
    while(curr->nxt){
      if(curr->child == child){
        fprintf(stderr,"Error: trying to bond already bonded symbols\n");
        return 0;
      }
      curr = curr->nxt;
    }
      
    curr->nxt = edge;
  }
  else
    parent->bonds = edge; 

  // set the previous for look back
  child->previous = parent; 

  child->num_edges++;
  parent->num_edges++;

  edge->parent = parent; 
  edge->child = child;
  edge->order = 1;
  return edge;
}


WLNEdge *search_edge(WLNSymbol *child, WLNSymbol*parent, bool verbose=true){
  if(!child || !parent){
    fprintf(stderr,"Error: searching edge on nullptrs\n");
    return 0;
  }
  
  WLNEdge *edge = 0;
  for (edge=parent->bonds;edge;edge = edge->nxt){
    if(edge->child == child)
      return edge;
  }
  if(verbose)
    fprintf(stderr,"Error: could not find edge in search\n");
  return 0;
}

WLNEdge *unsaturate_edge(WLNEdge *edge,unsigned int n){
  if(!edge){
    fprintf(stderr,"Error: unsaturating non-existent edge\n");
    return 0;
  }

  edge->order += n; 
  edge->parent->num_edges += n;
  edge->child->num_edges+= n;

  if(edge->parent->num_edges > edge->parent->allowed_edges){
    fprintf(stderr, "Error: wln character[%c] is exceeding allowed connections %d/%d\n", edge->parent->ch,edge->parent->num_edges, edge->parent->allowed_edges);
    return 0;
  }

  if(edge->child->num_edges > edge->child->allowed_edges){
    fprintf(stderr, "Error: wln character[%c] is exceeding allowed connections %d/%d\n", edge->child->ch,edge->child->num_edges, edge->child->allowed_edges);
    return 0;
  }

  return edge;
}

WLNEdge *saturate_edge(WLNEdge *edge,unsigned int n){
  if(!edge){
    fprintf(stderr,"Error: saturating non-existent edge\n");
    return 0;
  }

  if(edge->order < 2)
    return edge;
  
  edge->order -= n; 
  edge->parent->num_edges -= n;
  edge->child->num_edges -= n;

  return edge;
}


bool remove_edge(WLNSymbol *head,WLNEdge *edge){
  if(!head || !edge){
    fprintf(stderr,"Error: removing bond of non-existent symbols\n");
    return false;
  }
  
  head->num_edges--;
  edge->child->num_edges--;

  if(head->bonds == edge){
    head->bonds = 0;
    return true;
  }

  bool found = false;
  WLNEdge *search = head->bonds;

  WLNEdge *prev = 0;
  while(search){
    if(search == edge){ 
      found = true;
      break;
    }
    prev = search; 
    search = search->nxt;
  }

  if(!found){
    fprintf(stderr,"Error: trying to remove bond from wln character[%c] - bond not found\n",head->ch);
    return false;
  }
  else{
    WLNEdge *tmp = edge->nxt;
    prev->nxt = tmp;
    // dont null the edges as we use the mempool to release them
  }

  return true;
}


WLNEdge* add_methyl(WLNSymbol *head, WLNGraph &graph){

  WLNSymbol *carbon = AllocateWLNSymbol('C',graph);
  WLNSymbol *hydrogen = 0;
  WLNEdge   *edge = 0;
  carbon->set_edge_and_type(4); // used for hydrogens
  
  for(unsigned int i=0;i<3;i++){
    hydrogen = AllocateWLNSymbol('H',graph);
    hydrogen->set_edge_and_type(1);
    edge = AllocateWLNEdge(hydrogen,carbon,graph);
    if(!edge)
      return 0;
  }

  WLNEdge *bond = AllocateWLNEdge(carbon,head,graph);
  return bond; 
}


WLNSymbol* create_carbon_chain(WLNSymbol *head,unsigned int size, WLNGraph &graph){

  if (size > REASONABLE){
    fprintf(stderr,"Error: making carbon chain over 1024 long, reasonable molecule?\n");
    return 0;
  }
    
  head->ch = '1';
  head->set_edge_and_type(4);

  if(size == 1)
    return head;
  
  WLNEdge *edge = 0;
  WLNSymbol *prev = head;
  for(unsigned int i=0;i<size-1;i++){
    WLNSymbol* carbon = AllocateWLNSymbol('1',graph);
    carbon->set_edge_and_type(4); // allows hydrogen resolve
    edge = AllocateWLNEdge(carbon,prev,graph);
    if(!edge)
      return 0;
    prev = carbon;
  } 

  return prev;
}

bool add_diazo(WLNSymbol *head,WLNGraph &graph){

  WLNEdge *edge = 0;
  WLNSymbol *oxygen = 0;

  oxygen = AllocateWLNSymbol('O',graph);
  oxygen->set_edge_and_type(2,head->type);
  graph.charge_additions[oxygen] = -1;

  edge = AllocateWLNEdge(oxygen,head,graph);
  
  oxygen = AllocateWLNSymbol('O',graph);
  oxygen->set_edge_and_type(2,head->type);
  edge = AllocateWLNEdge(oxygen,head,graph);
  
  edge = unsaturate_edge(edge,1);
  if(!edge)
    return false;
  
  return true;
}


/* resolve carbon methyl assumptions */
bool resolve_methyls(WLNSymbol *target, WLNGraph &graph){

  switch(target->ch){

    case 'Y':
    case 'X':
    case 'K':
      while(target->num_edges < target->allowed_edges){
        if(!add_methyl(target,graph))
          return false;
      }
      target->num_edges = target->allowed_edges;
      break;

    default:
      fprintf(stderr,"Error: resolving methyls performed on invalid symbol: %c\n",target->ch);
      return false;
  }

  return true;
}



/**********************************************************************
                          WLNRing Functions
**********************************************************************/

WLNRing *AllocateWLNRing(WLNGraph &graph)
{
  graph.ring_count++;
  if(graph.ring_count > REASONABLE){
    fprintf(stderr,"Error: creating more than 1024 wln rings - is this reasonable?\n");
    return 0;
  }

  WLNRing *wln_ring = new WLNRing;
  graph.RINGS[graph.ring_count] = wln_ring;
  return wln_ring;
}


// both lookups needed for QOL in ring building
WLNSymbol* assign_locant(unsigned char loc,WLNSymbol *locant, WLNRing *ring){
    
  if(!locant)
    return 0;
  
  ring->locants[loc] = locant; 
  ring->locants_ch[locant] = loc;
  locant->type = RING;
  return locant; 
}  


bool assign_aromatics(std::deque<unsigned char> &ring_path, WLNRing *ring){
  std::map<WLNSymbol*, bool> assigned; 

  for (unsigned int i=1;i<ring_path.size();i++){
  
    WLNSymbol *par = ring->locants[ring_path[i-1]];
    WLNSymbol *chi = ring->locants[ring_path[i]];

    if(i == 1){
      switch(par->ch){
        case 'N':
          par->allowed_edges = 4;
        break;
        default:
          break;
      }
    }

    switch(chi->ch){
      case 'N':
        chi->allowed_edges = 4;
      break;

      default:
        break;
    }

    // cannot have double bonds following each other
    if(assigned[par] || assigned[chi])
      continue;

    // must have double bonds outside of the ring
    if(par->ch == 'Y' || chi->ch == 'Y')
      continue;

    if(par->num_edges < par->allowed_edges && chi->num_edges < chi->allowed_edges){
      assigned[par]  = true;
      assigned[chi]  = true;
      WLNEdge * edge = search_edge(chi,par,false);
      if(!edge){
        edge = search_edge(par,chi);
      }
      if(!unsaturate_edge(edge,1)){
        //fprintf(stderr,"par: %c  chi: %c\n",locants_ch[par],locants_ch[chi]);
        return false;
      } 
        
    }
  }
  

  return true;
}


/* creates poly rings, aromaticity is defined in reverse due to the nature of notation build
return the size of the ring, or zero if failure */
unsigned int CreatePolyCyclic(std::vector<std::pair<unsigned int,
                              unsigned char>> &ring_assignments, 
                              std::vector<bool> &aromaticity,
                              std::map<unsigned char,bool> &bridge_locants,
                              WLNRing *ring,
                              WLNGraph &graph)
  {
    
  unsigned int local_size = 0; 
  for (unsigned int i=0;i<ring_assignments.size();i++){
    std::pair<unsigned int, unsigned char> component = ring_assignments[i]; 
    if(local_size)
      local_size += component.first - 2;
    else
      local_size = component.first;
  }

  for (unsigned int i=0;i<252;i++){
    if(bridge_locants[i])
      local_size+= -1; 
  }

  // create all the nodes in a large straight chain
  
  WLNSymbol *curr= 0; 
  WLNSymbol *prev = 0; 
  for (unsigned int i=1;i<=local_size;i++){
    unsigned char loc = int_to_locant(i);
    if(!ring->locants[loc]){
      curr = AllocateWLNSymbol('C',graph);
      curr->set_edge_and_type(4,RING);
      curr = assign_locant(loc,curr,ring);
    }
    else
      curr = ring->locants[loc];

    if(prev){
      WLNEdge *edge = AllocateWLNEdge(curr,prev,graph);
      if(!edge)
        return false;
    }
    prev = curr;
  }


  // calculate bindings and then traversals round the loops
  unsigned int comp_size = 0;
  unsigned char bind_1 = '\0';
  unsigned char bind_2 = '\0';
  unsigned int fuses = 0; 
  bool aromatic = false; 


  for (unsigned int i=0;i<ring_assignments.size();i++){
    std::pair<unsigned int, unsigned char> component = ring_assignments[i];
    comp_size = component.first;
    bind_1 = component.second;
    aromatic = aromaticity[i];
    WLNSymbol *path = ring->locants[bind_1];

    while(bridge_locants[bind_1] && ring->locants[bind_1]->num_edges >= 2){
      bind_1++;
    }

    std::deque<unsigned char> ring_path;
    // first pair can be calculated directly without a path travel
    if(!fuses){
      bind_2 = bind_1 + comp_size - 1; // includes start atom
      for (unsigned int i=0; i<comp_size;i++)
        ring_path.push_back(bind_1+i);
    }
    else{
      //there needs to be a graph travel here taking the longest locant

      // 1. starting on bind_1, travel n-1 places through the maximal locant path, to calculate fuse

      unsigned char highest_loc = '\0';
      for (unsigned int i=0;i<comp_size - 1; i++){
        ring_path.push_back(ring->locants_ch[path]);

        WLNEdge *lc = 0;
        for(lc = path->bonds;lc;lc = lc->nxt){
          WLNSymbol *child = lc->child;
          unsigned char child_loc = ring->locants_ch[child];
          if(child_loc > highest_loc)
            highest_loc = child_loc;
        }    
        path = ring->locants[highest_loc];
      }

      ring_path.push_back(ring->locants_ch[path]); // add the last symbol
      bind_2 = highest_loc;
    }


    if(opt_debug){
      fprintf(stderr,"  %d  fusing: %c <-- %c   [",fuses,bind_2,bind_1);
      for (unsigned char ch : ring_path){
        fprintf(stderr," %c(%d)",ch,ch);
      }
      fprintf(stderr," ]\n");
    }


    WLNEdge *edge = AllocateWLNEdge(ring->locants[bind_2],ring->locants[bind_1],graph);
    if(!edge)
      return false;

    if(aromatic){
      if(!assign_aromatics(ring_path,ring))
        return false;
    }
      
    fuses++;
  }

  return local_size; 
}


/* interesting here that the multicyclic points are not explicitly used */
unsigned int CreateMultiCyclic( std::vector<std::pair<unsigned int,unsigned char>> &ring_assignments, 
                                std::vector<bool> &aromaticity,
                                std::vector<unsigned char> &multicyclic_locants,
                                std::vector<indexed_pair> &pseudo_locants,
                                std::set<unsigned char> &broken_locants,
                                std::map<unsigned char,bool> &bridge_locants,
                                unsigned char size_designator,
                                WLNRing *ring,
                                WLNGraph &graph) 
{
  
  // create a chain size of ring designator
  unsigned int local_size = locant_to_int(size_designator);

  std::map<unsigned int,std::vector<indexed_pair>>  pseudo_lookup;
  std::map<unsigned char,std::deque<unsigned char>> broken_lookup; // want to pop front
  std::map<unsigned char, bool>                     spawned_broken;
  std::map<unsigned char,unsigned int>              shared_rings; 

  // create all the nodes in a large straight chain
  WLNSymbol *curr = 0; 
  WLNSymbol *prev = 0; 
  for (unsigned int i=1;i<=local_size;i++){
    unsigned char loc = int_to_locant(i);
    if(!ring->locants[loc]){
      curr = AllocateWLNSymbol('C',graph);
      curr->set_edge_and_type(4,RING);
      curr = assign_locant(loc,curr,ring);
    }
    else{
      curr = ring->locants[loc];
      if(!ring->locants_ch[curr])
        ring->locants_ch[curr] = loc;
    }
      
    if(prev){
      WLNEdge *edge = AllocateWLNEdge(curr,prev,graph);
      if(!edge)
        return false;
    }
    prev = curr;
  }

  // shared rings cannot be done purely on ring paths annoyingly, leads to missteps in 
  // bridged ring system where the ring is spawned outside the path. 

  // pseudo pairs
  for (indexed_pair psd_pair : pseudo_locants)
    pseudo_lookup[psd_pair.index].push_back(psd_pair);
  
  // broken locants
  if(!broken_locants.empty()){
    // create the atoms, 
    for (unsigned char loc_broken : broken_locants){
      unsigned char calculate_origin = loc_broken;
      unsigned int pos = 0;
      while( (calculate_origin - 23) > 128){
        calculate_origin += -23;
        pos++;
      }
      // position here decodes where to link them
      unsigned char parent = '\0';
      parent = int_to_locant(128 + calculate_origin); // relative positioning
      if(pos == 2 || pos == 3)
        parent = locant_to_int(parent) + 128;
      else if(pos > 3){
        fprintf(stderr,"Error: non-locant links past a two-level tree are unsuitable for this parser\n");
        return false;
      }

      if(opt_debug)
        fprintf(stderr,"  ghost linking %d to parent %c\n",loc_broken,parent);
      
      if(!ring->locants[loc_broken]){
        // bond them in straight away
        WLNSymbol *broken = AllocateWLNSymbol('C',graph);
        broken->set_edge_and_type(4,RING);
        broken = assign_locant(loc_broken,broken,ring);
        broken_lookup[parent].push_back(loc_broken);
        WLNEdge *edge = AllocateWLNEdge(ring->locants[loc_broken],ring->locants[parent],graph);
        if(!edge)
          return false; 
      }
      else{
        fprintf(stderr,"Error: branching locants are overlapping created elements already in the locant path\n");
        return false;
      }
    }
  }

  // calculate bindings and then traversals round the loops
  unsigned int comp_size = 0;
  unsigned char bind_1 = '\0';
  unsigned char bind_2 = '\0';
  unsigned int fuses = 0; 
  bool aromatic = false;

  for (unsigned int i=0;i<ring_assignments.size();i++){
    std::pair<unsigned int, unsigned char> component = ring_assignments[i];
    comp_size = component.first;
    bind_1 = component.second;
    aromatic = aromaticity[i];
    WLNSymbol *path = ring->locants[bind_1];
    
    unsigned int predefined = 1;
    std::deque<unsigned char> ring_path; 
    
    // --- PSD BRIDGE ONLY ---
    if(!pseudo_lookup[i].empty()){

      indexed_pair psd_pair = pseudo_lookup[i].front(); // should only be 1
      
      bind_1 = psd_pair.bind_1;
      bind_2 = psd_pair.bind_2;

      // the binding is easy, its just figuring out what the ring path would be without
      // bidirectional travel - dfs is unreliable with this struct design

      // sometimes the notation lists the psd bridges that can be implied from the rings

      if(!search_edge(ring->locants[bind_2],ring->locants[bind_1],false)){
        WLNEdge *edge = AllocateWLNEdge(ring->locants[bind_2],ring->locants[bind_1],graph);
        if(!edge)
          return false;

        ring_path.push_back(bind_1);
        ring_path.push_back(bind_2);

        if(opt_debug){
          fprintf(stderr,"  %d  fusing: %c <-- %c   [",fuses,bind_2,bind_1);
          for (unsigned char ch : ring_path){
            fprintf(stderr," %c(%d)",ch,ch);
          }
          fprintf(stderr," ]\n");
        }
      }
      
      fuses++;
      continue;      
    }


    while( (shared_rings[bind_1] >= 2 && broken_lookup[bind_1].empty())
          || (bridge_locants[bind_1] && ring->locants[bind_1]->num_edges >= 2) ){

      predefined++;
      bind_1++;
      ring_path.push_back(bind_1);
    }

    // whilst the cpp object is large, the alternative is many more lines
    while(!broken_lookup[bind_1].empty()){
      predefined++;
      unsigned char bloc = broken_lookup[bind_1].front();
      broken_lookup[bind_1].pop_front();
      bind_1 = bloc;
      ring_path.push_back(bind_1);
    }

    // --- MULTI ALGORITHM --- 
    
    ring_path.push_back(ring->locants_ch[path]);

    unsigned char highest_loc = '\0'; 
    for (unsigned int i=0;i<comp_size - predefined; i++){  // 2 points are already defined
      
      highest_loc = '\0'; // highest of each child iteration 

      WLNEdge *lc = 0;
      for(lc = path->bonds;lc;lc = lc->nxt){
        WLNSymbol *child = lc->child;
        unsigned char child_loc = ring->locants_ch[child];

        // skip the broken child if not yet included in a ring
        if(child_loc > 128 && !spawned_broken[child_loc])
          continue;

        if(child_loc >= highest_loc)
          highest_loc = child_loc;
      }

      if(!highest_loc){
        if(locant_to_int(ring->locants_ch[path]) == local_size)
          highest_loc = ring->locants_ch[path];
        else{
          fprintf(stderr,"Error: locant path formation is broken in ring definition - '%c(%d)'\n",ring->locants_ch[path],ring->locants_ch[path]);
          return false;
        }
      }

      path = ring->locants[highest_loc];
      ring_path.push_back(ring->locants_ch[path]);
    }

    bind_2 = highest_loc; 

    // annoying catch needed for bridge notation that is 'implied' 
    if(i == ring_assignments.size() - 1 && bind_2 != int_to_locant(local_size)){
      
      unsigned char back = ring_path.back();
      while(back < int_to_locant(local_size) && !ring_path.empty()){
        back++;
        ring_path.push_back(back);
        ring_path.pop_front();
      }
      bind_2 = back;
      
      if(!ring_path.empty())
        bind_1 = ring_path.front();
    }
    if(opt_debug){
      fprintf(stderr,"  %d  fusing: %c <-- %c   [",fuses,bind_2,bind_1);
      for (unsigned char ch : ring_path){
        fprintf(stderr," %c(%d)",ch,ch);
      }
      fprintf(stderr," ]\n");
    }

    // updates the shared rings approach -> paths must MUST be exact for this to work
    for (unsigned char ch : ring_path){
      shared_rings[ch]++; 
      if(ch > 128)
        spawned_broken[ch] = true;
    }
        
    WLNEdge *edge = AllocateWLNEdge(ring->locants[bind_2],ring->locants[bind_1],graph);
    if(!edge)
      return false;


    if(aromatic){
      if(!assign_aromatics(ring_path,ring))
        return false;
    }

    fuses++;
  }

  return local_size; 
}


unsigned char create_relative_position(unsigned char parent){
  // A = 129
  unsigned int relative = 128 + locant_to_int(parent);
  if(relative > 252){
    fprintf(stderr,"Error: relative position is exceeding 252 allowed space - is this is suitable molecule for WLN notation?\n");
    return '\0';
  }
  else
    return relative;
}


// try to handle if errors are occuring due to path changes
bool handle_post_orders(std::vector<std::pair<unsigned char, unsigned char>> &bonds, 
                        unsigned int mode,
                        unsigned int final_size,
                        WLNRing *ring){

  // post unsaturate bonds
  for (std::pair<unsigned char, unsigned char> bond_pair : bonds){
    unsigned char loc_1 = bond_pair.first;
    unsigned char loc_2 = bond_pair.second;
    if(loc_2 > int_to_locant(final_size)){
      loc_1 = 'A';
      loc_2--;
    }

    WLNEdge *edge = search_edge(ring->locants[loc_2],ring->locants[loc_1],false);
    if(!edge)
      edge = search_edge(ring->locants[loc_1],ring->locants[loc_2],false);
    
    if(!edge)
      edge = search_edge(ring->locants[loc_1]->bonds->child,ring->locants[loc_1],false);

    if(!edge)
      edge = search_edge(ring->locants[loc_1]->previous->bonds->child,ring->locants[loc_1],false);

    if(mode)
      edge = unsaturate_edge(edge,1);
    else
      edge = saturate_edge(edge,1);

    if(!edge)
      return false;
  }

  return true;
}

/* parse the WLN ring block, use ignore for already predefined spiro atoms */
void FormWLNRing(WLNRing *ring,std::string &block, unsigned int start, WLNGraph &graph,unsigned char spiro_atom='\0'){


  enum RingType{ POLY=1, PERI=2, BRIDGED=3, PSDBRIDGED = 4}; 
  const char* ring_strings[] = {"MONO","POLY","PERI","BRIDGED","PSDBRIDGED"};
  unsigned int ring_type = POLY;   // start in mono and climb up

  bool warned             = false;  // limit warning messages to console
  bool heterocyclic       = false;  // L|T designator can throw warnings

  // -- paths -- // 
  // int allows way more description in states

  unsigned int state_multi          = 0; // 0 - closed, 1 - open multi notation, 2 - expect size denotation
  unsigned int state_pseudo         = 0; 
  unsigned int state_aromatics      = 0;
  unsigned int state_chelate        = 0;
  bool implied_assignment_used      = false; // allows a shorthand if wanted, but not mixing
  
  unsigned int expected_locants       = 0;
  unsigned int  evaluating_break      = 0;
  unsigned char ring_size_specifier   = '\0';
  unsigned char positional_locant     = '\0';
  unsigned int last_locant_position   = 0;

  std::string special;  

  std::vector<bool> aromaticity; 
  std::vector<std::pair<unsigned char, unsigned char>>  unsaturations;
  std::vector<std::pair<unsigned char, unsigned char>>  saturations; 

  std::vector<unsigned char>    pseudo_locants;
  std::vector<unsigned int>     pseudo_positions; 
  std::vector<unsigned char>    multicyclic_locants;
  std::set<unsigned char>       broken_locants;
  std::map<unsigned char,bool>  bridge_locants;
  
  // broken locants start at A = 129 for extended ascii 
  // first is the standard 'X-' second is 'X-&', third is 'X--', fourth is 'X--&' and so on

  std::vector<std::pair<unsigned int, unsigned char>>  ring_components;
  std::vector<indexed_pair>                            indexed_bindings;  
  

  unsigned int i = 0;    
  unsigned int len = block.size();

  // somehow need to give 

  const char *block_str = block.c_str(); // this should be now globally alive
  unsigned char ch = *block_str++;

  while(ch){

    switch(ch){

      // specials

      case ' ':

        if(state_multi == 3)
          state_multi = 0;

        if(evaluating_break){
          broken_locants.insert(positional_locant);

          if(state_multi == 1 && expected_locants){
            multicyclic_locants.back() = positional_locant;
            state_multi = 2;
            expected_locants--;
          }
          else if(state_pseudo == 1 && expected_locants){
            pseudo_locants.back() = positional_locant;
            expected_locants--;
          }

          evaluating_break = 0;
        }
        if(expected_locants){
          fprintf(stderr,"Error: %d locants expected before space character\n",expected_locants);
          Fatal(i+start);
        }
        else if(state_multi == 1){
          state_multi = 2;
        }
        state_pseudo = 0;
        positional_locant = '\0'; // hard resets on spaces
        break;

      case '&':
        if (state_aromatics){
          aromaticity.push_back(1);
          break;
        }
        else if (state_multi == 3){
          ring_size_specifier += 23;
        }
        else if (state_pseudo){
          pseudo_locants.back() += 23;
        }
        else if(positional_locant){
          
          // can only be an extension of a positional multiplier for a branch
          if (last_locant_position && last_locant_position == i-1)
            positional_locant += 23;
          else{
            state_aromatics = 1;
            aromaticity.push_back(1);
          }
        }
        else{
          // if unhandled, then it must be an aromatic start
          state_aromatics = 1;
          aromaticity.push_back(1);
        }
        break;

      case '/':
        if(state_aromatics){
          fprintf(stderr,"Error: character '%c' cannot be in the aromaticity assignment block\n",ch);
          Fatal(i+start);
        }

        if(!pseudo_positions.empty() && pseudo_positions.back() == ring_components.size() -1){
          for (unsigned int p=0;p<pseudo_positions.size();p++){
            pseudo_positions[p] += - 1; // back shift;
          }
        }
        pseudo_positions.push_back(ring_components.size() -1); // which ring has the pseudo?
        expected_locants = 2; 
        state_pseudo = 1;
        break; 
      

      // turn this into a look ahead type bias in order to significantly tidy this up
      case '-':{

        // gives us a local working copy
        char *local_arr = new char [strlen(block_str)+1]; 
        memset(local_arr,'\0',strlen(block_str)+1);
        strcpy(local_arr,block_str);
        const char *local = local_arr;

        unsigned char local_ch = *(local)++; // moved it over
        unsigned int gap = 0; 
        bool found_next = false;
        
        while(local_ch != '\0'){
          if(local_ch == ' ')
            break;
          if(local_ch == '-'){
            // this calculates the gap to the next '-'
            found_next = true;
            break;
          }
          special.push_back(local_ch);
          gap++;
          local_ch = *local++;
        }

        if(local_arr){
          delete [] local_arr;
          local = 0;
          local_arr = 0;
        }

        // this will change on metallocenes defintions
        if( (state_multi || state_pseudo) && expected_locants){
          gap = 0;
        }

        // could ignore gap zeros as they will come back round again, therefore need a size check

        if(found_next){
          
          // pointer is moved by 1 at the bottom, so its positions -1 
          switch(gap){
            case 0:

              // we resolve only the first one
              evaluating_break = 1; // on a space, number or character, we push to broken_locants

              if(positional_locant){
                if(positional_locant < 128){
                  positional_locant = create_relative_position(positional_locant); // i believe breaking modifier will then get removed
                  last_locant_position = i;
                  if(!positional_locant)
                    Fatal(i+start);
                }
                else{
                  // this means its already been moved, so we move the locant 23+23 across
                  if(positional_locant + 46 > 252){
                    fprintf(stderr,"Error: branching locants are exceeding the 252 space restriction on WLN notation, is this a reasonable molecule?\n");
                    Fatal(start+i);
                  }
                  positional_locant += 46;
                  last_locant_position = i;
                  // no need to move the global pointer here
                }
              }
              else{
                fprintf(stderr,"Error: trying to branch out character without starting point\n");
                Fatal(start+i);  
              }

              break;
            case 1:{
              if(!implied_assignment_used && !positional_locant){
                implied_assignment_used = true;
                positional_locant = 'A';
              }
              // this can only be hypervalent element
  
              if(positional_locant){
                if(spiro_atom){
                  if(positional_locant == spiro_atom){
                    positional_locant++;
                    block_str+=2; 
                    i+=2;
                    break;
                  }
                  else if(ring->locants[positional_locant]){
                    positional_locant++;
                    if(positional_locant == spiro_atom){
                      positional_locant++;
                      block_str+=2; 
                      i+=2;
                      break;
                    }
                  }
                }
                else if(ring->locants[positional_locant])
                  positional_locant++;


                WLNSymbol* new_locant = assign_locant(positional_locant,define_hypervalent_element(special[0],graph),ring);  // elemental definition 
                if(!new_locant)
                  Fatal(i+start);

                graph.string_positions[start+i + 1] = new_locant; // attaches directly

                if(opt_debug)
                  fprintf(stderr,"  assigning hypervalent %c to position %c\n",special[0],positional_locant);
              }
              else{
                fprintf(stderr,"Error: trying to assign element without starting point\n");
                Fatal(start+i);  
              }
              block_str+=2; 
              i+=2; 
              break;
            }
            case 2:{
              if(!implied_assignment_used && !positional_locant){
                implied_assignment_used = true;
                positional_locant = 'A';
              }

              if(std::isdigit(special[0])){
                for(unsigned char dig_check : special){
                  if(!std::isdigit(dig_check)){
                    fprintf(stderr,"Error: mixing numerical and alphabetical special defintions is not allowed\n");
                    Fatal(start+i);
                  }
                }
                if(positional_locant)
                  ring_components.push_back({std::stoi(special),positional_locant}); //big ring
                else
                  ring_components.push_back({std::stoi(special),'A'});
              }
              else{

                if(positional_locant){

                  if(spiro_atom){
                    if(positional_locant == spiro_atom){
                      positional_locant++;
                      block_str+=3; 
                      i+=3;
                      break;
                    }
                    else if(ring->locants[positional_locant]){
                      positional_locant++;
                      if(positional_locant == spiro_atom){
                        positional_locant++;
                        block_str+=3; 
                        i+=3;
                        break;
                      }
                    }
                  }
                  else if(ring->locants[positional_locant])
                    positional_locant++;
                  
                  WLNSymbol* new_locant = assign_locant(positional_locant,define_element(special,graph),ring);  // elemental definition
                  if(!new_locant)
                    Fatal(i+start);

                  graph.string_positions[start+i + 1] = new_locant; // attaches directly to the starting letter

                  if(opt_debug)
                    fprintf(stderr,"  assigning element %s to position %c\n",special.c_str(),positional_locant);
                }
                else{
                  fprintf(stderr,"Error: trying to assign element without starting point\n");
                  Fatal(start+i);  
                }
              }

              block_str+=3; 
              i+=3;              
              break;
            }
            default:
              fprintf(stderr,"Error: %d numerals incased in '-' brackets is unreasonable for WLN to create\n",gap);
              Fatal(start+i);
          }


        }
        else if(i > 0 && block[i-1] == '&')
          state_aromatics = 1;
        else{

          // if there wasnt any other symbol, it must be a notation extender
          evaluating_break = 1; // on a space, number or character, we push to broken_locants

          if(positional_locant){
            if(positional_locant < 128){
              positional_locant = create_relative_position(positional_locant); // i believe breaking modifier will then get removed
              last_locant_position = i;
              if(!positional_locant)
                Fatal(i+start);
            }
            else{
              // this means its already been moved, so we move the locant 23+23 across
              if(positional_locant + 46 > 252){
                fprintf(stderr,"Error: branching locants are exceeding the 252 space restriction on WLN notation, is this a reasonable molecule?\n");
                Fatal(start+i);
              }
              positional_locant += 46;
              last_locant_position = i;
            }
          }
          else{
            fprintf(stderr,"Error: trying to branch out character without starting point\n");
            Fatal(start+i);  
          }

        }


        special.clear();
        break;
      }

      // numerals - easy access

      case '0':
        // place the minus charges on the last ring seen
        if(ring_components.size() == 1){
          ring->post_charges.push_back({'B',-1});
        }else{
          unsigned int track = 0;
          for (unsigned int rn = 0; rn<ring_components.size()-1;rn++)
            track += ring_components[rn].first;
          ring->post_charges.push_back({int_to_locant(track + 1),-1});  // post is needed as this is before pointer assignment
        }
        
        break;

      case '1':
      case '2':
      case '3':
      case '4':
      case '5':
      case '6':
      case '7':
      case '8':
      case '9':
        if(state_aromatics){
          fprintf(stderr,"Error: character '%c' cannot be in the aromaticity assignment block\n",ch);
          Fatal(i+start);
        }

        if(evaluating_break){
          broken_locants.insert(positional_locant);

          if(state_multi == 1 && expected_locants){
            multicyclic_locants.back() = positional_locant;
            expected_locants--;
          }
          else if(state_pseudo == 1 && expected_locants){
            pseudo_locants.back() = positional_locant;
            expected_locants--;
          }

          evaluating_break = 0;
        }
        if (i > 1 && block[i-1] == ' '){
          state_multi   = 1; // enter multi state
          expected_locants = ch - '0';
        }
        else{

          if(positional_locant)
            ring_components.push_back({ch-'0',positional_locant});
          else
            ring_components.push_back({ch-'0','A'});

          positional_locant = '\0';
        }
        break;

      case 'A':
      case 'B':
      case 'C':
      case 'D':
      case 'E':
      case 'F':
      case 'G':
      case 'H':
      case 'I':
      case 'K':
      case 'M':
      case 'N':
      case 'O':
      case 'P':
      case 'Q':
      case 'R':
      case 'S':
      case 'U':
      case 'V':
      case 'W':
      case 'X':
      case 'Y':
      case 'Z':
        if(i == 0 && ch == 'D'){
          state_chelate = 1;
          heterocyclic = true;
          break;
        }

        if(state_aromatics){
          fprintf(stderr,"Error: character '%c' cannot be in the aromaticity assignment block\n",ch);
          Fatal(i+start);
        }

        if(evaluating_break){
          broken_locants.insert(positional_locant);

          if(state_multi == 1 && expected_locants)
            multicyclic_locants.back() = positional_locant;
          else if(state_pseudo == 1 && expected_locants)
            pseudo_locants.back() = positional_locant;
          
          evaluating_break = 0;
        }

        if(expected_locants){

          positional_locant = ch; // use for look back
          expected_locants--;

          if(state_multi)
            multicyclic_locants.push_back(ch);
          else if (state_pseudo)
            pseudo_locants.push_back(ch);
          else{
            fprintf(stderr,"Error: unhandled locant rule\n");
            Fatal(start+i);
          }
          break;
        }
        else if(state_multi == 2){
          ring_size_specifier = ch;
          state_multi = 3;
        }
        else if (positional_locant){

          if(spiro_atom && positional_locant == spiro_atom){
            positional_locant++;
            break;
          }
            
          if (opt_debug)
            fprintf(stderr,"  assigning WLNSymbol %c to position %c\n",ch,positional_locant);

          WLNSymbol *new_locant = 0; 

          switch(ch){
            
            case 'D':
              if(!state_chelate){
                fprintf(stderr,"Error: %c is not allowed as a atom assignment within ring notation\n",ch);
                Fatal(start+i);
              }
              // means open chelating bond
              break;

            case 'S':
            case 'P':
              if(!heterocyclic)
                warned = true;

              if(ring->locants[positional_locant])
                positional_locant++; 
              
              if(spiro_atom && positional_locant == spiro_atom){
                positional_locant++;
                break;
              }  

              new_locant = AllocateWLNSymbol(ch,graph);
              new_locant = assign_locant(positional_locant,new_locant,ring);
              new_locant->set_edge_and_type(5,RING);
              break;

            case 'X':
            case 'Y':
              if(ring->locants[positional_locant])
                positional_locant++;

              if(spiro_atom && positional_locant == spiro_atom){
                positional_locant++;
                break;
              }  

              new_locant = AllocateWLNSymbol(ch,graph);
              new_locant = assign_locant(positional_locant,new_locant,ring);
              new_locant->set_edge_and_type(4,RING);
              break;

            case 'Z': // treat as NH2
            case 'N':
              if(!heterocyclic)
                warned = true;

              if(ring->locants[positional_locant])
                positional_locant++;

              if(spiro_atom && positional_locant == spiro_atom){
                positional_locant++;
                break;
              }  

              new_locant = AllocateWLNSymbol(ch,graph);
              new_locant = assign_locant(positional_locant,new_locant,ring);
              new_locant->set_edge_and_type(3,RING);
              break;

            case 'V':
              if(ring->locants[positional_locant])
                positional_locant++;

              //  if this is a chelting oxygen, treat as nOU without any restrictions
              new_locant = AllocateWLNSymbol(ch,graph);
              new_locant = assign_locant(positional_locant,new_locant,ring);
              new_locant->set_edge_and_type(2,RING);
              break;

            case 'M':
            case 'O':
              if(!heterocyclic)
                warned = true;

              if(ring->locants[positional_locant])
                positional_locant++; 

              new_locant = AllocateWLNSymbol(ch,graph);
              new_locant = assign_locant(positional_locant,new_locant,ring);
              new_locant->set_edge_and_type(2,RING);
              break;

            case 'K':
              if(!heterocyclic)
                warned = true;

              if(ring->locants[positional_locant])
                positional_locant++; 
              
              if(spiro_atom && positional_locant == spiro_atom){
                positional_locant++;
                break;
              }  

              new_locant = AllocateWLNSymbol(ch,graph);
              new_locant = assign_locant(positional_locant,new_locant,ring);
              new_locant->set_edge_and_type(4,RING);
              break;

            case 'U':
              // no need to put this in implied, it has to be specified
              if(i < len - 3 && block[i+1] == '-' && block[i+2] == ' '){
                unsaturations.push_back({positional_locant,block[i+3]});
                block_str += 3;
                i += 3;
              }
              else
                unsaturations.push_back({positional_locant,positional_locant+1});
              break;

            case 'W':
              switch(ring->locants[positional_locant]->ch){
                case 'K':
                  ring->locants[positional_locant]->allowed_edges++;
                  break;
                default:
                  break;
              }
              if(!add_diazo(ring->locants[positional_locant],graph))
                Fatal(i+start);
              break;

            // has the effect of unsaturating a bond
            case 'H':
              // no need to put this in implied, it has to be specified
              saturations.push_back({positional_locant,positional_locant+1});
              break;


            default:
              fprintf(stderr,"Error: %c is not allowed as a atom assignment within ring notation\n",ch);
              Fatal(start+i);
          }

          graph.string_positions[start+i] = new_locant;
        }
        else{
          if( (i>0 && i<len-1 && block[i-1] == ' ') && (block[i+1] == ' ' || block[i+1] == 'T' || block[i+1] == 'J') ){
            if(ring_components.empty()){
              fprintf(stderr,"Error: assigning bridge locants without a ring\n");
              Fatal(start+i);
            }
            else
              bridge_locants[ch] = true;
          }
          else if(i>0 && block[i-1] == ' '){
            positional_locant = ch;
            last_locant_position = i;
          }
          else{
            implied_assignment_used = true;
            positional_locant = 'A';

            if(spiro_atom && positional_locant == spiro_atom){
              positional_locant++;
              break;
            }

            if (opt_debug)
              fprintf(stderr,"  assigning WLNSymbol %c to position %c\n",ch,positional_locant);

            WLNSymbol *new_locant = 0; 

            switch(ch){
              
              case 'D':
              if(!state_chelate){
                fprintf(stderr,"Error: %c is not allowed as a atom assignment within ring notation\n",ch);
                Fatal(start+i);
              }
              break;

              case 'S':
              case 'P':
                if(!heterocyclic)
                  warned = true;

                if(ring->locants[positional_locant])
                  positional_locant++; 

                if(spiro_atom && positional_locant == spiro_atom){
                  positional_locant++;
                  break;
                }  

                new_locant = AllocateWLNSymbol(ch,graph);
                new_locant = assign_locant(positional_locant,new_locant,ring);
                new_locant->set_edge_and_type(5,RING);
                break;

              case 'X':
              case 'Y':
                if(ring->locants[positional_locant])
                  positional_locant++;

                if(spiro_atom && positional_locant == spiro_atom){
                  positional_locant++;
                  break;
                }  

                new_locant = AllocateWLNSymbol(ch,graph);
                new_locant = assign_locant(positional_locant,new_locant,ring);
                new_locant->set_edge_and_type(4,RING);
                break;


              case 'Z':
              case 'N':
                if(!heterocyclic)
                  warned = true;
                
                if(ring->locants[positional_locant])
                  positional_locant++; 

                if(spiro_atom && positional_locant == spiro_atom){
                  positional_locant++;
                  break;
                }  

                new_locant = AllocateWLNSymbol(ch,graph);
                new_locant = assign_locant(positional_locant,new_locant,ring);
                new_locant->set_edge_and_type(3,RING);
                break;

              case 'V':
                if(ring->locants[positional_locant])
                  positional_locant++; 

                new_locant = AllocateWLNSymbol(ch,graph);
                new_locant = assign_locant(positional_locant,new_locant,ring);
                new_locant->set_edge_and_type(2,RING);
                break;

              case 'M':
              case 'O':
                if(!heterocyclic)
                  warned = true;
                
                if(ring->locants[positional_locant])
                  positional_locant++; 

                new_locant = AllocateWLNSymbol(ch,graph);
                new_locant = assign_locant(positional_locant,new_locant,ring);
                new_locant->set_edge_and_type(2,RING);
                break;

              case 'K':
                if(!heterocyclic)
                  warned = true;

                if(ring->locants[positional_locant])
                  positional_locant++; 
                
                if(spiro_atom && positional_locant == spiro_atom){
                  positional_locant++;
                  break;
                }  

                new_locant = AllocateWLNSymbol(ch,graph);
                new_locant = assign_locant(positional_locant,new_locant,ring);
                new_locant->set_edge_and_type(4,RING);
                break;

              case 'U':
                unsaturations.push_back({positional_locant,positional_locant+1});
              break;

                // has the effect of unsaturating a bond
              case 'H':
                // no need to put this in implied, it has to be specified
                saturations.push_back({positional_locant,positional_locant+1});
                break;

              default:
                fprintf(stderr,"Error: %c is not allowed as a atom assignment within ring notation\n",ch);
                Fatal(start+i);
            }

            graph.string_positions[start+i] = new_locant;
          }
        }

        break;

      case 'L':
        if(state_aromatics){
          fprintf(stderr,"Error: character '%c' cannot be in the aromaticity assignment block\n",ch);
          Fatal(i+start);
        }

        if(evaluating_break){
          broken_locants.insert(positional_locant);

          if(state_multi == 1 && expected_locants)
            multicyclic_locants.back() = positional_locant;
          else if(state_pseudo == 1 && expected_locants)
            pseudo_locants.back() = positional_locant;

          evaluating_break = 0;
        }

        if(i==0){
          heterocyclic = false; 
          break;
        }
        if(expected_locants){

          positional_locant = ch; // use for look back
          expected_locants--;

          if(state_multi)
            multicyclic_locants.push_back(ch);
          else if (state_pseudo)
            pseudo_locants.push_back(ch);
          else{
            fprintf(stderr,"Error: unhandled locant rule\n");
            Fatal(start+i);
          }

          break;
        }
        else if(state_multi == 2){
          ring_size_specifier = ch;
          state_multi = 3;
        }
        else{
          if(i>0 && block[i-1] == ' '){
            positional_locant = ch;
            last_locant_position = i;
          }
          else{
            fprintf(stderr,"Error: symbol '%c' is in an unhandled state, please raise issue if this notation is 100%% correct\n",ch);
            Fatal(i+start);
          }
        }
      
        break;


      case 'T':
        if(state_aromatics){
          aromaticity.push_back(0);
          break;
        }
      
        if(evaluating_break){
          broken_locants.insert(positional_locant);

          if(state_multi == 1 && expected_locants)
            multicyclic_locants.back() = positional_locant;
          else if(state_pseudo == 1 && expected_locants)
            pseudo_locants.back() = positional_locant;

          evaluating_break = 0;
        }

        if(i==0){
          heterocyclic = true; 

          break;
        }

        if(expected_locants){

          positional_locant = ch; // use for look back
          expected_locants--;

          if(state_multi)
            multicyclic_locants.push_back(ch);
          else if (state_pseudo)
            pseudo_locants.push_back(ch);
          else{
            fprintf(stderr,"Error: unhandled locant rule\n");
            Fatal(start+i);
          }
          break;
        }
        else if(state_multi == 2){
          ring_size_specifier = ch;
          state_multi = 3;
        }
        else{
          if(i>0 && block[i-1] == ' ' && block[i+1] != 'J'){
            positional_locant = ch;
            last_locant_position = i;
          }
          else{
            // this must be an aromatic state right?
            state_aromatics = 1;
            aromaticity.push_back(0);
          }
        }
          

        break;

      
      // CLOSE

      case 'J':
        if(state_aromatics)
          state_aromatics = 0;
        
        if(evaluating_break){
          broken_locants.insert(positional_locant);

          if(state_multi == 1 && expected_locants)
            multicyclic_locants.back() = positional_locant;
          else if(state_pseudo == 1 && expected_locants)
            pseudo_locants.back() = positional_locant;

          evaluating_break = 0;
        }

        if (i == block.size()-1){
          
          if(ring_components.empty()){
            fprintf(stderr,"Error: error in reading ring components, check numerals in ring notation\n");
            Fatal(start+i);
          }

          if (pseudo_locants.size() > 0)
            ring_type = PSDBRIDGED;

          if (multicyclic_locants.size() > 0 && ring_type < PSDBRIDGED)
            ring_type = PERI;

          if (aromaticity.size() == 1 && aromaticity[0] == false){
            while(aromaticity.size() < ring_components.size())
              aromaticity.push_back(false);
          }
          else if (aromaticity.empty()){
            while(aromaticity.size() < ring_components.size())
              aromaticity.push_back(true);
          }

          // perform the aromatic denotion check
          if (ring_components.size() != aromaticity.size()){
            fprintf(stderr,"Error: mismatch between number of rings and aromatic assignments - %ld vs expected %ld\n",aromaticity.size(),ring_components.size());
            Fatal(i+start);
          }

          // create the bindings needed for pseudo bridges
          for (unsigned int i=0; i< pseudo_positions.size();i++){
            indexed_pair pseudo; 
            pseudo.set(pseudo_locants[i+i],pseudo_locants[i+i+1],pseudo_positions[i]);
            indexed_bindings.push_back(pseudo);
          }

          break;
        }
        if(expected_locants){

          positional_locant = ch; // use for look back
          expected_locants--;

          if(state_multi)
            multicyclic_locants.push_back(ch);
          else if (state_pseudo)
            pseudo_locants.push_back(ch);
          else{
            fprintf(stderr,"Error: unhandled locant rule\n");
            Fatal(start+i);
          }
          break;
        }
        else if(state_multi == 2){
          ring_size_specifier = ch;
          state_multi = 3;
        }
        else{
          if(i>0 && block[i-1] == ' '){
            positional_locant = ch;
            last_locant_position = i;
          }
          else{
            fprintf(stderr,"Error: symbol '%c' is in an unhandled state, please raise issue if this notation is 100%% correct\n",ch);
            Fatal(i+start);
          }
        }
        
        break;

      default:
        break;
    }
    
    i++;
    ch = *(block_str++);
  }

  if(warned)
    fprintf(stderr,"Warning: heterocyclic ring notation required for inter atom assignment, change starting 'L' to 'T'\n");
  

  // debug here
  if (opt_debug){
    
    fprintf(stderr,"  ring type: %s\n",ring_strings[ring_type]);

    fprintf(stderr,"  ring components: ");
    for (std::pair<unsigned int, unsigned char> comp : ring_components){
      
      if(comp.second > 'Z')
        fprintf(stderr,"%d(%d) ",comp.first,comp.second);
      else
        fprintf(stderr,"%d(%c) ",comp.first,comp.second);
    } 
      
    fprintf(stderr,"\n");

    fprintf(stderr,"  aromaticity: ");
    for (bool aromatic : aromaticity)
      fprintf(stderr,"%d ",aromatic);
    fprintf(stderr,"\n");

    fprintf(stderr,"  multicyclic points: ");
    for (unsigned char loc : multicyclic_locants){
      if(loc > 'Z')
        fprintf(stderr,"%d ",loc);
      else
        fprintf(stderr,"%c ",loc);
    }
    fprintf(stderr,"\n");

    fprintf(stderr,"  broken path points: ");
    for (unsigned char loc : broken_locants){
      fprintf(stderr,"%d ",loc);
    }
    fprintf(stderr,"\n");

    fprintf(stderr,"  bridge points: ");
    for (unsigned int i=0;i<252;i++){
      if(bridge_locants[i])
        fprintf(stderr,"%c ",i);
    }
    fprintf(stderr,"\n");

    fprintf(stderr,"  pseudo bridge points: ");
    for (unsigned int i=0; i< pseudo_positions.size();i++){
      fprintf(stderr,"(%d)[%c <-- %c] ",pseudo_positions[i],pseudo_locants[i+i],pseudo_locants[i+i+1]);
    }
    fprintf(stderr,"\n");

    fprintf(stderr,"  size denotion: %d\n",ring_size_specifier ? locant_to_int(ring_size_specifier) : 0);
    fprintf(stderr,"  heterocyclic: %s\n", heterocyclic ? "yes":"no");
  }
  
  unsigned int final_size = 0;
  switch(ring_type){
    case POLY:
      final_size = CreatePolyCyclic(ring_components,aromaticity,bridge_locants,ring,graph);
      break;
    case PERI:
    case PSDBRIDGED:
      final_size = CreateMultiCyclic(ring_components,aromaticity,
                                multicyclic_locants,indexed_bindings,
                                broken_locants,
                                bridge_locants,
                                ring_size_specifier,
                                ring,
                                graph);
    break;
  }

  for (std::pair<unsigned char,int> &post : ring->post_charges)
    graph.charge_additions[ring->locants[post.first]] += post.second;
  

  if (!final_size)
    Fatal(start+i);

  
  if( !handle_post_orders(unsaturations,1,final_size,ring) 
      || !handle_post_orders(saturations,0,final_size,ring))
    Fatal(start+i);

}












/* must be performed before sending to obabel graph*/
bool ExpandWLNSymbols(WLNGraph &graph){

  unsigned int stop = graph.symbol_count;
  for (unsigned int i=1;i<=stop;i++){
    WLNSymbol *sym = graph.SYMBOLS[i];

    switch(sym->ch){

      case 'Y':
      case 'X':
      case 'K':
        resolve_methyls(sym,graph);
        break;

      case 'V':{
        sym->ch = 'C';
        sym->set_edge_and_type(4);
        WLNSymbol *oxygen = AllocateWLNSymbol('O',graph);
        oxygen->set_edge_and_type(2);
        WLNEdge *e = AllocateWLNEdge(oxygen,sym,graph);
        e = unsaturate_edge(e,1);
        if(!e)
          return false;
        break;
      }
      
      case 'W':
        sym->ch = 'C';
        sym->set_edge_and_type(4);
        if(!add_diazo(sym,graph))
          return false;
        break;


      default:
        break; // ignore
    }
  }
  
  return true; 
}


/* backwards search for tentative ionic rule procedures */
unsigned int search_ionic(const char *wln_ptr, unsigned int len,
                          std::vector<std::pair<unsigned int, int>> &charges)
{
  unsigned int first_instance = 0;

  for (unsigned int i=0;i<len;i++){

    // these are required in blocks of 5
    if(wln_ptr[i] == ' ' && wln_ptr[i+1] == '&')
    {
      
      std::string position_1;
      std::string position_2;

      unsigned int local_search = i+2;

      if(std::isdigit(wln_ptr[i+2])){

        while(std::isdigit(wln_ptr[local_search])){
          position_1.push_back(wln_ptr[local_search]);
          local_search++;

          if(local_search > len)
            return first_instance;
        }
      }
      else 
        continue;

      // local search should now be pointing at the '\'
      if(wln_ptr[local_search] == '/')
        local_search++;
      else
        continue;

      if(std::isdigit(wln_ptr[local_search])){
        while(std::isdigit(wln_ptr[local_search])){
          position_2.push_back(wln_ptr[local_search]);
          local_search++;

          if(local_search > len)
            return first_instance;
        }
      }
      else 
        continue;

      
      if(std::stoi(position_1) != 0)
        charges.push_back({std::stoi(position_1),1});
      
      if(std::stoi(position_2) != 0)
        charges.push_back({std::stoi(position_2),-1});

      if(!first_instance)
        first_instance = i;
    }
  }

  return first_instance;
}


/* uses the global position map */
bool AssignCharges(std::vector<std::pair<unsigned int, int>> &charges,WLNGraph &graph){
  if(charges.empty())
    return true;

  for (std::pair<unsigned int, int> pos_charge : charges){
    WLNSymbol *assignment = graph.string_positions[pos_charge.first - 1]; // reindex as wln 1 is string 0
    if(!assignment){
      fprintf(stderr,"Error: trying to assign ionic charge to unavaliable element, check that character %d is avaliable for assignment\n",pos_charge.first);
      return false;
    }
    else{
      graph.charge_additions[assignment] += pos_charge.second;

      if(opt_debug){
        fprintf(stderr, "  character at position [%d] has the following charge addition - %d\n",pos_charge.first,pos_charge.second);
      }
    }
  }
  return true;
}


/**********************************************************************
                         High Level Parser Functions
**********************************************************************/


/* returns the head of the graph, parse all normal notation */
bool ParseWLNString(const char *wln_ptr, WLNGraph &graph) 
{
  
  // keep the memory alive

  if (opt_debug)
    fprintf(stderr, "Parsing WLN notation: %s\n",wln_ptr);

  ObjectStack branch_stack;   // access to both rings and symbols
  branch_stack.reserve(100);  // reasonable size given

  std::vector<std::pair<unsigned int, int>> ionic_charges;
  
  WLNSymbol *curr       = 0;
  WLNSymbol *prev       = 0;
  WLNEdge   *edge       = 0;
  WLNRing   *ring       = 0;
  WLNRing   *wrap_ring  = 0;

  bool pending_locant           = false;
  bool pending_J_closure        = false;
  bool pending_inline_ring      = false;
  bool pending_spiro            = false;
  bool pending_diazo            = false;
  bool pending_ring_in_ring     = false; // rings in rings


  unsigned char on_locant = '\0';         // locant tracking
  unsigned int pending_unsaturate = 0;    // 'U' style bonding
  bool j_skips = false;                   // handle skipping of 'J' if in cyclic notation legitimately 

  std::string special;  // special elemental definitions
  
  // allows consumption of notation after block parses
  unsigned int block_start = 0;
  unsigned int block_end = 0;

  unsigned int len = strlen(wln_ptr);
  unsigned int zero_position = search_ionic(wln_ptr,len,ionic_charges);

  unsigned int i=0;
  unsigned char ch = *wln_ptr;
  
  while(ch)
  {  
    
    // dont read any ionic notation
    if(zero_position && zero_position == i)
      break;

    switch (ch)
    {

    case '0': // cannot be lone, must be an addition to another num
      if(pending_J_closure)
        break;

      else if (pending_locant){
        
        if(prev && prev->type != RING)
          graph.charge_additions[prev]++;

        prev = 0;
        on_locant = '0';
        pending_locant = false;
      }
      else{
        fprintf(stderr,"Error: a lone zero mark is not allowed without positive numerals either side\n");
        Fatal(i);
      }
      break;

    case '1':
    case '2':
    case '3':
    case '4':
    case '5':
    case '6':
    case '7':
    case '8':
    case '9':
      if (pending_J_closure){
        // small addition to allow J handling in points
        if(i > 0 && wln_string[i-1] == ' ')
          j_skips = true;
        
        break;
      }
        
      else if(pending_locant){  // handle all multiplier contractions
        
        // block to hunt the multiplier maximum --> lookahead
        // look for spaces followed by ints. 

        std::string int_sequence;
        int_sequence.push_back(ch);

        while(i < len - 2){
          if(wln_string[i+1] == ' ' && std::isdigit(wln_string[i+2])){
            int_sequence.push_back(wln_string[i+2]);
            i+=2;
            wln_ptr += 2;
          }
          else
            break;
        }

        // pointer is moved to the last number, multiplier value is calculated
        unsigned int multiplier_value = std::stoi(int_sequence);

        fprintf(stderr,"Error: multipliers are not currently supported\n");
        Fatal(i);

        pending_locant = false;
        on_locant = ch;
      }
      else if(pending_ring_in_ring && pending_inline_ring){
          // onlocant holds the char needed to wrap the ring back, 
        
        if(on_locant != '0'){
          curr = wrap_ring->locants[on_locant];
          if(!curr){
            fprintf(stderr,"Error: cannot access looping ring structure\n");
            Fatal(i);
          }

          if(prev){
            edge = AllocateWLNEdge(curr,prev,graph);
            if(pending_unsaturate){
              edge = unsaturate_edge(edge,pending_unsaturate);
              pending_unsaturate = 0;
            }
            if(!edge)
              Fatal(i);
          }
          else
            Fatal(i);
        }
        
        // last notation is not neccessary
        while(wln_ptr){
          if(*wln_ptr == 'J')
            break;
          wln_ptr++;
          i++;
        }

        pending_ring_in_ring = false;
        pending_inline_ring = false;
      } 
      else{
        on_locant = '\0';

        curr = AllocateWLNSymbol('1',graph);
        curr->set_edge_and_type(4);

        if(pending_diazo){
          if(!add_diazo(curr,graph))
            Fatal(i-1);
          pending_diazo = false;
        }

        if(prev){
          edge = AllocateWLNEdge(curr,prev,graph);
          if(!edge)
            Fatal(i);

          if(pending_unsaturate){
            edge = unsaturate_edge(edge,pending_unsaturate);
            pending_unsaturate = 0;
          }
          
        }
        
        std::string int_sequence;
        int_sequence.push_back(ch);

        while(*(wln_ptr+1)){
          if(!std::isdigit(*(wln_ptr+1)))
            break;

          int_sequence.push_back(*(wln_ptr+1));
          wln_ptr++;
          i++;
        }

        curr = create_carbon_chain(curr,std::stoi(int_sequence),graph);
        if(!curr){
          fprintf(stderr,"Error: error in creating carbon chain, raise algorithm issue\n");
          Fatal(i);
        }


        prev = curr;
        break;
      }
      break;
    

    case 'Y':
      if (pending_J_closure)
        break;
      else if (pending_locant)
      {
        fprintf(stderr,"Error: '%c' cannot be a locant assignment, please expand [A-W] with &\n",ch);
        Fatal(i);
      }
      else
      {
        on_locant = '\0';
        
        curr = AllocateWLNSymbol(ch,graph);
        curr->set_edge_and_type(3);

        if(pending_diazo){
          if(!add_diazo(curr,graph))
            Fatal(i-1);
          pending_diazo = false;
        }

        if(prev){
          edge = AllocateWLNEdge(curr,prev,graph);
          if(!edge)
            Fatal(i);

          if(pending_unsaturate){
            edge = unsaturate_edge(edge,pending_unsaturate);
            pending_unsaturate = 0;
          }
          
        }
        
        branch_stack.push({0,curr});
        graph.string_positions[i] = curr;
        pending_unsaturate = 0;
        prev = curr;
      }
      break;

    case 'X':
      if (pending_J_closure)
        break;
      else if (pending_locant)
      {
        fprintf(stderr, "Wiswesser Uncertainities will produce multiple smiles per X entry\n"
                        "since the number of these is at least the size of the ring system\n"
                        "its likely to blow memory allocations, as such they are not supported\n");
        Fatal(i);
      }
      else
      {
        on_locant = '\0';
        
        curr = AllocateWLNSymbol(ch,graph);
        curr->set_edge_and_type(4);

        if(pending_diazo){
          if(!add_diazo(curr,graph))
            Fatal(i-1);
          pending_diazo = false;
        }

        if(prev){
          edge = AllocateWLNEdge(curr,prev,graph);
          if(!edge)
            Fatal(i);

          if(pending_unsaturate){
            edge = unsaturate_edge(edge,pending_unsaturate);
            pending_unsaturate = 0;
          }
          
        }
        
        branch_stack.push({0,curr});
        graph.string_positions[i] = curr;
        prev = curr;
      }
      break;

      // oxygens

    case 'O':
      if (pending_J_closure)
        break;
      else if (pending_locant)
      {
        if(!pending_inline_ring){
          ring = branch_stack.ring;
          if(!ring)
            Fatal(i);
        
          
          curr = ring->locants[ch];
          if(!curr){
            fprintf(stderr,"Error: accessing locants out of range\n");
            Fatal(i);
          }
          
          prev = curr;
        }

        pending_locant = false;
        on_locant = ch;
      }
      else
      {
        on_locant = '\0';
        if(pending_diazo){
          fprintf(stderr,"Error: diazo assignment to an oxygen is a disallowed bond type\n");
          Fatal(i);
        }

        curr = AllocateWLNSymbol(ch,graph);
        curr->set_edge_and_type(2);

        if(prev){
          edge = AllocateWLNEdge(curr,prev,graph);
          if(!edge)
            Fatal(i);

          if(pending_unsaturate){
            edge = unsaturate_edge(edge,pending_unsaturate);
            pending_unsaturate = 0;
          }
          
        }

        graph.string_positions[i] = curr;
        prev = curr;
      }
      break;

    case 'Q':
      if (pending_J_closure)
        break;
      else if (pending_locant)
      {
        if(!pending_inline_ring){
          ring = branch_stack.ring;
          if(!ring)
            Fatal(i);

          curr = ring->locants[ch];
          if(!curr){
            fprintf(stderr,"Error: accessing locants out of range\n");
            Fatal(i);
          }
          prev = curr;
        }

        pending_locant = false;
        on_locant = ch;
      }
      else
      {
        on_locant = '\0';
        if(pending_diazo){
          fprintf(stderr,"Error: diazo assignment to an oxygen is a disallowed bond type\n");
          Fatal(i);
        }

        curr = AllocateWLNSymbol(ch,graph);
        curr->set_edge_and_type(1);

        if(prev){
          edge = AllocateWLNEdge(curr,prev,graph);
          
          if(pending_unsaturate){
            edge = unsaturate_edge(edge,pending_unsaturate);
            pending_unsaturate = 0;
          }
          if(!edge)
            Fatal(i);
        }

        graph.string_positions[i] = curr;
        pending_unsaturate = 0;
        prev = return_object_symbol(branch_stack);
        if(!prev)
          prev = curr;
      }
      break;

    case 'V':
      if (pending_J_closure)
        break;
      else if (pending_locant)
      {

        if(!pending_inline_ring){
          ring = branch_stack.ring;
          if(!ring)
            Fatal(i);

          curr = ring->locants[ch];
          if(!curr){
            fprintf(stderr,"Error: accessing locants out of range\n");
            Fatal(i);
          }
          
          prev = curr;
        }

        pending_locant = false;
        on_locant = ch;
      }
      else
      {
        on_locant = '\0';
        if(pending_diazo){
          fprintf(stderr,"Error: diazo assignment to an carbonyl is a disallowed bond type\n");
          Fatal(i);
        }
        curr = AllocateWLNSymbol(ch,graph);
        curr->set_edge_and_type(2,STANDARD);
        if(prev){
          edge = AllocateWLNEdge(curr,prev,graph);
          if(pending_unsaturate){
            edge = unsaturate_edge(edge,pending_unsaturate);
            pending_unsaturate = 0;
          }
          
          if(!edge)
            Fatal(i);
        }

        graph.string_positions[i] = curr;
        prev = curr;
      }
      break;

    case 'W':
      if (pending_J_closure)
        break;
      else if (pending_locant)
      {

        if(!pending_inline_ring){
          ring = branch_stack.ring;
          if(!ring)
            Fatal(i);

          curr = ring->locants[ch];
          if(!curr){
            fprintf(stderr,"Error: accessing locants out of range\n");
            Fatal(i);
          }
          
          prev = curr;
        }
        pending_locant = false;
        on_locant = ch;
      }
      else
      {
        on_locant = '\0';
        if(pending_diazo){
          fprintf(stderr,"Error: double diazo assignment is a disallowed bond type\n");
          Fatal(i);
        }

        if(prev){
          prev->allowed_edges++;
          if(!add_diazo(prev,graph))
            Fatal(i);
        }
        else
          pending_diazo = true;
        
        graph.string_positions[i] = curr;
        pending_unsaturate = 0;
        prev = curr;
      }
      break;

      // nitrogens

    case 'N':
      if (pending_J_closure)
        break;
      else if (pending_locant)
      {

        if(!pending_inline_ring){
          ring = branch_stack.ring;
          if(!ring)
            Fatal(i);

          curr = ring->locants[ch];
          if(!curr){
            fprintf(stderr,"Error: accessing locants out of range\n");
            Fatal(i);
          }
          prev = curr;
        }

        pending_locant = false;
        on_locant = ch;
      }
      else
      {
        on_locant = '\0';
        
        curr = AllocateWLNSymbol(ch,graph);
        curr->set_edge_and_type(3);

        if(pending_diazo){
          curr->allowed_edges++; // special allowance for Nitro
          if(!add_diazo(curr,graph))
            Fatal(i-1);
          pending_diazo = false;
        }

        if(prev){
          edge = AllocateWLNEdge(curr,prev,graph);
          if(pending_unsaturate){
            edge = unsaturate_edge(edge,pending_unsaturate);
            pending_unsaturate = 0;
          }
          if(!edge)
            Fatal(i);
        }
        
        branch_stack.push({0,curr});
        graph.string_positions[i] = curr;
        pending_unsaturate = 0;
        prev = curr;
      }
      break;

    case 'M':
      if (pending_J_closure)
        break;
      else if (pending_locant)
      {
        if(!pending_inline_ring){
          
          ring = branch_stack.ring;
          if(!ring)
            Fatal(i);
          
          curr = ring->locants[ch];
          if(!curr){
            fprintf(stderr,"Error: accessing locants out of range\n");
            Fatal(i);
          }
          
          prev = curr;
        }
        pending_locant = false;
        on_locant = ch;
      }
      else
      {
        on_locant = '\0';
        if(pending_diazo){
          fprintf(stderr,"Error: diazo assignment to NH is a disallowed bond type\n");
          Fatal(i);
        }

        curr = AllocateWLNSymbol(ch,graph);
        curr->set_edge_and_type(2);

        if(prev){
          edge = AllocateWLNEdge(curr,prev,graph);
          if(pending_unsaturate){
            edge = unsaturate_edge(edge,pending_unsaturate);
            pending_unsaturate = 0;
          }
          if(!edge)
            Fatal(i);
          
        }

        graph.string_positions[i] = curr;
        pending_unsaturate = 0;
        prev = curr;
      }
      break;

    case 'K':
      if (pending_J_closure)
        break;
      else if (pending_locant)
      {
        if(!pending_inline_ring){

          ring = branch_stack.ring;
          if(!ring)
            Fatal(i);

          curr = ring->locants[ch];
          if(!curr){
            fprintf(stderr,"Error: accessing locants out of range\n");
            Fatal(i);
          }
          
          prev = curr;
        }

        pending_locant = false;
        on_locant = ch;
      }
      else
      {
        on_locant = '\0';
        
    
        curr = AllocateWLNSymbol(ch,graph);
        curr->set_edge_and_type(4);

        if(pending_diazo){
          if(!add_diazo(curr,graph))
            Fatal(i-1);
          pending_diazo = false;
        }

        if(prev){
          edge = AllocateWLNEdge(curr,prev,graph);
          if(pending_unsaturate){
            edge = unsaturate_edge(edge,pending_unsaturate);
            pending_unsaturate = 0;
          }
          if(!edge)
            Fatal(i);
          
        }
        
        branch_stack.push({0,curr});
        graph.string_positions[i] = curr;
        prev = curr;
      }
      break;

    case 'Z':
      if (pending_J_closure)
        break;
      else if (pending_locant)
      {
        if(!pending_inline_ring){
          ring = branch_stack.ring;
          if(!ring)
            Fatal(i);

          curr = ring->locants[ch];
          if(!curr){
            fprintf(stderr,"Error: accessing locants out of range\n");
            Fatal(i);
          }
          
          prev = curr;
        }
      
        pending_locant = false;
        on_locant = ch;
      }
      else
      { 
        on_locant = '\0';
        if(pending_diazo){
          fprintf(stderr,"Error: diazo assignment to NH2 is a disallowed bond type\n");
          Fatal(i);
        }

        curr = AllocateWLNSymbol(ch,graph);
        curr->set_edge_and_type(1);

        if(prev){
          edge = AllocateWLNEdge(curr,prev,graph);
          if(pending_unsaturate){
            edge = unsaturate_edge(edge,pending_unsaturate);
            pending_unsaturate = 0;
          }
          if(!edge)
            Fatal(i);
          
        }

        graph.string_positions[i] = curr;
        pending_unsaturate = 0;
        prev = return_object_symbol(branch_stack);
        if(!prev)
          prev = curr;
      }
      break;

      // halogens - need to add rules for semi allowed hyper valence in ionions

    case 'E':
    case 'G':
    case 'F':
    case 'I':
      if (pending_J_closure)
        break;
      else if (pending_locant)
      {
        if(!pending_inline_ring){
          ring = branch_stack.ring;
          if(!ring)
            Fatal(i);

          curr = ring->locants[ch];
          if(!curr){
            fprintf(stderr,"Error: accessing locants out of range\n");
            Fatal(i);
          }
          
          prev = curr;
        }
        pending_locant = false;
        on_locant = ch;
      }
      else
      {
        on_locant = '\0';
        if(pending_diazo){
          fprintf(stderr,"Error: diazo assignment to a non expanded valence halogen is a disallowed bond type\n");
          Fatal(i);
        }

        curr = AllocateWLNSymbol(ch,graph);
        curr->set_edge_and_type(1);

        if(prev){
          edge = AllocateWLNEdge(curr,prev,graph);
          if(pending_unsaturate){
            edge = unsaturate_edge(edge,pending_unsaturate);
            pending_unsaturate = 0;
          }
          
          if(!edge)
            Fatal(i);
          
        }

        graph.string_positions[i] = curr;
        pending_unsaturate = 0;
        prev = return_object_symbol(branch_stack);
        if(!prev)
          prev = curr;
      }
      break;

      // inorganics

    case 'B':
      if (pending_J_closure)  
        break;
      else if (pending_locant)
      {
        if(!pending_inline_ring){
          ring = branch_stack.ring;
          if(!ring)
            Fatal(i);

          curr = ring->locants[ch];
          if(!curr){
            fprintf(stderr,"Error: accessing locants out of range\n");
            Fatal(i);
          }
          
          prev = curr;
        }
        pending_locant = false;
        on_locant = ch;
      }
      else
      {
        on_locant = '\0';
  
        curr = AllocateWLNSymbol(ch,graph);
        curr->set_edge_and_type(3);

        if(pending_diazo){
          if(!add_diazo(curr,graph))
            Fatal(i-1);
          pending_diazo = false;
        }

        if(prev){
          edge = AllocateWLNEdge(curr,prev,graph);
          if(pending_unsaturate){
            edge = unsaturate_edge(edge,pending_unsaturate);
            pending_unsaturate = 0;
          }
        
          if(!edge)
            Fatal(i);
          
        }         
        
        branch_stack.push({0,curr});
        graph.string_positions[i] = curr;
        prev = curr;
      }
      break;

    case 'P':
    case 'S':
      if (pending_J_closure)
        break;
      else if (pending_locant)
      {
        if(!pending_inline_ring){
          ring = branch_stack.ring;
          if(!ring)
            Fatal(i);

          curr = ring->locants[ch];
          if(!curr){
            fprintf(stderr,"Error: accessing locants out of range\n");
            Fatal(i);
          }
          
          prev = curr;
        }
        pending_locant = false;
        on_locant = ch;
      }
      else
      {
        on_locant = '\0';


        curr = AllocateWLNSymbol(ch,graph);
        curr->set_edge_and_type(6);

        if(pending_diazo){
          if(!add_diazo(curr,graph))
            Fatal(i-1);
          pending_diazo = false;
        }

        if(prev){
          edge = AllocateWLNEdge(curr,prev,graph);
          if(pending_unsaturate){
            edge = unsaturate_edge(edge,pending_unsaturate);
            pending_unsaturate = 0;
          }
      
          if(!edge)
            Fatal(i);
          
        }
        
        branch_stack.push({0,curr});
        graph.string_positions[i] = curr;
        prev = curr;
      }
      break;

    // multiply bonded carbon, therefore must be at least a double bond
    case 'C':
      if (pending_J_closure)
        break;
      else if (pending_locant)
      {
        if(!pending_inline_ring){
          ring = branch_stack.ring;
          if(!ring)
            Fatal(i);

          curr = ring->locants[ch];
          if(!curr){
            fprintf(stderr,"Error: accessing locants out of range\n");
            Fatal(i);
          }
          
          prev = curr;
        }
        pending_locant = false;
        on_locant = ch;
      }
      else
      {
        on_locant = '\0';

        curr = AllocateWLNSymbol(ch,graph);
        curr->set_edge_and_type(4);

        if(pending_diazo){
          if(!add_diazo(curr,graph))
            Fatal(i-1);
          pending_diazo = false;
        }

        if(prev && i < len - 1){

          edge = AllocateWLNEdge(curr,prev,graph);
          if(pending_unsaturate){
            edge = unsaturate_edge(edge,pending_unsaturate);
            pending_unsaturate = 0;
          }
          if(!edge)
            Fatal(i);
          
          
        }
        else{
          fprintf(stderr,"Error: 'C' carbon designation must be able to reach max valence without hydrogens\n");
          Fatal(i);
        }

        graph.string_positions[i] = curr;
        prev = curr;
      }
      break;

    case 'A':
      if (pending_J_closure)
        break;
      else if (pending_locant)
      {

        if(!pending_inline_ring){
          ring = branch_stack.ring;
          if(!ring)
            Fatal(i);

          curr = ring->locants[ch];
          if(!curr){
            fprintf(stderr,"Error: accessing locants out of range\n");
            Fatal(i);
          }

          prev = curr;
        }

        pending_locant = false;
        on_locant = ch;
      }
      else{
        fprintf(stderr,"Error: locant only symbol used in atomic definition\n");
        Fatal(i);
      }
      break;
        
    // this can start a chelating ring compound, so has the same block as 'L\T'
    case 'D':
      if (pending_J_closure)
        break;
      else if (pending_locant)
      {
        if(!pending_inline_ring){
          ring = branch_stack.ring;
          if(!ring)
            Fatal(i);

          curr = ring->locants[ch];
          if(!curr){
            fprintf(stderr,"Error: accessing locants out of range\n");
            Fatal(i);
          }
          
          prev = curr;
        }
        pending_locant = false;
        on_locant = ch;
      }
      else
      {
        if(i < len - 2 && wln_string[i+1] == '-' && (wln_string[i+2] == 'T' || wln_string[i+2] == 'L')){
          pending_ring_in_ring = true;

          i++;
          wln_ptr++;
          pending_inline_ring = true;
          break;
        }
          
        if (i == 0)
          pending_inline_ring = true;
        
      
        if (!pending_inline_ring)
        {
          fprintf(stderr, "Error: chelating ring notation started without '-' denotion\n");
          Fatal(i);
        }
        
        pending_inline_ring = false;
        block_start = i;
        pending_J_closure = true;
      }
      break;
        
        
    // hydrogens explicit

    case 'H':
      if (pending_J_closure)
        break;
      else if (pending_locant)
      {
        if(!pending_inline_ring){

          ring = branch_stack.ring;
          if(!ring)
            Fatal(i);

          curr = ring->locants[ch];
          if(!curr){
            fprintf(stderr,"Error: accessing locants out of range\n");
            Fatal(i);
          }
          
          prev = curr;
        }
        pending_locant = false;
        on_locant = ch;
      }
      else{
        on_locant = '\0';

        // explicit hydrogens
        curr = AllocateWLNSymbol(ch,graph);
        curr->set_edge_and_type(1);

        
        if(prev){
          // will add with more examples
          edge = AllocateWLNEdge(curr,prev,graph);
          if(pending_unsaturate){
            edge = unsaturate_edge(edge,pending_unsaturate);
            pending_unsaturate = 0;
          }
          if(!edge)
            Fatal(i);

          switch(prev->ch){
            case 'Z':
              graph.charge_additions[prev]++;
              prev->allowed_edges++;
              break;
            
            default:
              break;
          }

          
        }

        graph.string_positions[i] = curr;
        curr = prev;
        // dont update for H
      }
      break;

      // ring notation

    case 'J':
      if(pending_J_closure && j_skips)
        break;
      if (pending_locant)
      {
        if(!pending_inline_ring){

          ring = branch_stack.ring;
          if(!ring)
            Fatal(i);

          curr = ring->locants[ch];
          if(!curr){
            fprintf(stderr,"Error: accessing locants out of range\n");
            Fatal(i);
          }
          
          prev = curr;
        }
        pending_locant = false;
        on_locant = ch;
      }

      else if (pending_J_closure 
              && ( (i<len-1 && (wln_string[i+1] == ' ' || wln_string[i+1] == '&') && wln_string[i-1] != ' ') 
              || i == len -1)
              )     
      {
        block_end = i;
        
        ring = AllocateWLNRing(graph);
        std::string r_notation = get_notation(block_start,block_end);

        if(pending_spiro){

          ring->locants[on_locant] = prev;

          // check for an aromaticity bond move?
          if(prev->allowed_edges - prev->num_edges < 2){

            // spiro would not be possible here, check if a double bond can be shifted
            WLNEdge *e = 0;
            WLNSymbol *shift = 0;
            for (e = prev->bonds;e;e = e->nxt){
              if (e->order == 2){
                e = saturate_edge(e,1);
                if(!e)
                  Fatal(i);

                shift = e->child;
                break;
              }
            }
            unsigned char next_loc = branch_stack.ring->locants_ch[shift]+1;
            if(!next_loc)
              next_loc = 'A'; // must of done the full loop

            e = search_edge(branch_stack.ring->locants[next_loc],shift);
            e = unsaturate_edge(e,1);
            if(!e)
              Fatal(i);
          }
          
          FormWLNRing(ring,r_notation,block_start,graph,on_locant);
        }
        else
          FormWLNRing(ring,r_notation,block_start,graph);
        

        if(pending_ring_in_ring && !wrap_ring)
          wrap_ring = ring; // instant back access


        branch_stack.push({ring,0});

        block_start = 0;
        block_end = 0;

        // does the incoming locant check

        if(pending_spiro)
          pending_spiro = false;
        else if (prev && on_locant && on_locant != '0')
        {
          if (ring->locants[on_locant]){
            edge = AllocateWLNEdge(ring->locants[on_locant],prev,graph);
            if(pending_unsaturate){
              edge = unsaturate_edge(edge,pending_unsaturate);
              pending_unsaturate = 0;
            }
            if(!edge)
              Fatal(i);
          }   
          else
          {
            fprintf(stderr, "Error: attaching inline ring with out of bounds locant assignment\n");
            Fatal(i);
          }
        }

        on_locant = '\0';
        pending_J_closure = false;
      }
      
      break;

    case 'L':
    case 'T':
      if (pending_J_closure)
        break;
      else if (pending_locant)
      {
        if(!pending_inline_ring){
          ring = branch_stack.ring;
          if(!ring)
            Fatal(i);

          curr = ring->locants[ch];
          if(!curr){
            fprintf(stderr,"Error: accessing locants out of range\n");
            Fatal(i);
          }
          
          prev = curr;
        }
        pending_locant = false;
        on_locant = ch;
      }
      else
      {
        if(i < len - 2 && wln_string[i+1] == '-' && (wln_string[i+2] == 'T' || wln_string[i+2] == 'L')){
          pending_ring_in_ring = true;

          // if pending inline ring doesnt get a L or T, we know it
          // the ring wrap symbol. 

          i++;
          wln_ptr++;
          pending_inline_ring = true;
          break;
        }
          
        if (i == 0)
          pending_inline_ring = true;
        
      
        if (!pending_inline_ring)
        {
          fprintf(stderr, "Error: ring notation started without '-' denotion\n");
          Fatal(i);
        }
        
        pending_inline_ring = false;
        block_start = i;
        pending_J_closure = true;
      }
      break;

    case 'R':
      if (pending_J_closure)
        break;
      
      else if (pending_locant)
      {
        if(!pending_inline_ring){

          ring = branch_stack.ring;
          if(!ring)
            Fatal(i);

          curr = ring->locants[ch];
          if(!curr){
            fprintf(stderr,"Error: accessing locants out of range\n");
            Fatal(i);
          }
          
          prev = curr;
        }
        pending_locant = false;
        on_locant = ch;
      }
      else
      {
        on_locant = '\0';
        ring = AllocateWLNRing(graph);

        std::string r_notation = "L6J";
        FormWLNRing(ring,r_notation,i,graph);
        branch_stack.push({ring,0});

        curr = ring->locants['A'];
        if(prev){
          edge = AllocateWLNEdge(curr,prev,graph);
          
          if(pending_unsaturate){
            edge = unsaturate_edge(edge,pending_unsaturate);
            pending_unsaturate = 0;
          }
          if(!edge)
            Fatal(i);
        }

        graph.string_positions[i] = curr;
        prev = curr;;
      }
      break;

      // bonding

    case 'U':
      if (pending_J_closure)
        break;
      else if (pending_locant)
      {
        if(!pending_inline_ring){

          ring = branch_stack.ring;
          if(!ring)
            Fatal(i);

          curr = ring->locants[ch];
          if(!curr){
            fprintf(stderr,"Error: accessing locants out of range\n");
            Fatal(i);
          }
          
          prev = curr;
        }
        pending_locant = false;
        on_locant = ch;
      }
      else if (pending_diazo){
        fprintf(stderr,"Error: diazo assignment followed by a bond increase is a disallowed bond type\n");
        Fatal(i);
      }
      else{
        on_locant = '\0';
        pending_unsaturate++;
      }
      break;

      // specials

    case ' ':
      if (pending_J_closure){
        j_skips = false;
        break;
      }
      else if (pending_diazo){
        fprintf(stderr,"Error: diazo assignment followed by a space seperator is a disallowed bond type\n");
        Fatal(i);
      }

      // only burn the stacks now on ionic clearance

      pending_locant = true;
      break;



    case '&':
      if (pending_diazo){
        fprintf(stderr,"Error: diazo assignment followed by a branch terminator is a disallowed bond type\n");
        Fatal(i);
      }

      if (pending_J_closure)
        break;
      
      if (pending_inline_ring)
      {
        // spiro notation open
        pending_spiro = true;
      }
      else if (pending_locant)
      {
        // ionic species or spiro, reset the linkings
        prev = 0;
        curr = 0;
        pending_locant = false;
        
        branch_stack.clear_all(); // burn stack
      }
      else if(on_locant){
        curr->ch += 23;
      }

      // this is always a ring pop?
      else if(i < len - 1 && wln_string[i+1] == ' '){

        branch_stack.pop(); // forced closure
        ring = branch_stack.ring;
        if(!ring){
          fprintf(stderr,"Error: popping too many rings, check '&' count\n");
          Fatal(i);
        }
        break;
      }
      else
      {
        if(!branch_stack.empty() && branch_stack.top().second){
          WLNSymbol *top = 0;
          top = branch_stack.top().second;

            // this means a <Y|X|..>'&' so handle methyl
          if(prev && prev == top){
            
            switch(prev->ch){
              // methyl contractions
              case 'X':
              case 'Y':
              case 'K':
                if(prev->num_edges < prev->allowed_edges){
                  if(!add_methyl(prev,graph))
                    Fatal(i);

                  prev = return_object_symbol(branch_stack);
                }
                else{ 
                  // we pop, 
                  branch_stack.pop();
                  prev = branch_stack.branch; 
                }
                break;

              // no contractions possible, we pop the stack
              default:

                // default pop
                branch_stack.pop();
                prev = branch_stack.branch; // if prev is nulled, then a ring is active
                break;
            }
          }
          else{
            // means a closure is done, we return to the first avaliable symbol on the branch stack
            prev = return_object_symbol(branch_stack);
          }
        }
        else if(!branch_stack.empty() && branch_stack.top().first){
          branch_stack.pop();
          ring = branch_stack.ring;
          prev = branch_stack.branch; 
        }
        else{
          fprintf(stderr,"Error: '&' punctuation outside of branching chains is disallowed notation\n");
          Fatal(i);
        }
      }
      break;


    case '-':{
      if (pending_J_closure)
        break;

      else if (pending_inline_ring)
      { 

        if(pending_ring_in_ring){
          // onlocant holds the char needed to wrap the ring back, 
          curr = wrap_ring->locants[on_locant];
          if(!curr){
            fprintf(stderr,"Error: cannot access looping ring structure\n");
            Fatal(i);
          }

          if(prev){
            edge = AllocateWLNEdge(curr,prev,graph);
            if(pending_unsaturate){
              edge = unsaturate_edge(edge,pending_unsaturate);
              pending_unsaturate = 0;
            }
            if(!edge)
              Fatal(i);
          }
          else
            Fatal(i);
    
          // last notation is not neccessary
          while(wln_ptr){
            if(*wln_ptr == 'J')
              break;
            wln_ptr++;
            i++;
          }

          pending_ring_in_ring = false;
          pending_inline_ring = false;

        }
        else{
          fprintf(stderr, "Error: only one pending ring can be active, check closures\n");
          Fatal(i);
        }
      }
#ifdef HIGH_LEVEL_ASSIGNMENTS
      else if(on_locant){
        if(prev->ch < 128){
          prev->ch += 128;
        }
        else
          prev->ch += 46; 
      }
#endif
      else{

        // ahh variable length array.. gotcha
        char *local_arr  = new char [strlen(wln_ptr)+1]; 
        memset(local_arr,'\0',strlen(wln_ptr)+1);
        strcpy(local_arr,wln_ptr);
        const char *local = local_arr;
        
        unsigned char local_ch = *(++local); // moved it over
        unsigned int gap = 0; 
        bool found_next = false;

        while(local_ch != '\0'){
          if(local_ch == ' ')
            break;
          if(local_ch == '-'){
            // this calculates the gap to the next '-'
            found_next = true;
            break;
          }
          special.push_back(local_ch);
          gap++;
          local_ch = *(++local);
        }

        if(local_arr){
          delete [] local_arr;
          local = 0;
          local_arr = 0;
        }
        
        if(!found_next){

          pending_inline_ring = true;

          if(!prev && branch_stack.branch)
            prev = branch_stack.branch;
        }
        else{

          if(gap == 1){
            curr = define_hypervalent_element(special[0],graph);
            if(!curr)
              Fatal(i);
          }
          else if(gap == 2){
            curr = define_element(special,graph);
            if(!curr)
              Fatal(i); 
            
            if(on_locant == '0'){
              graph.charge_additions[curr]++;
              //on_locant = '\0';
            }
                            
          }
          else{
            fprintf(stderr,"Error: special '-' must be either 1 or 2 symbols - %d seen\n",gap);
            Fatal(i);
          }

          if(pending_diazo){
            if(!add_diazo(curr,graph))
              Fatal(i-1);
            pending_diazo = false;
          }
          
          if(prev){
            if(!gap && ring)
              edge = AllocateWLNEdge(ring->locants[prev->ch],prev,graph);
            else
              edge = AllocateWLNEdge(curr,prev,graph);

            if(!edge)
              Fatal(i);
            
          }

          branch_stack.push({0,curr});

          i+= gap+1;
          wln_ptr+= gap+1;

          graph.string_positions[i-gap] = curr;
          pending_unsaturate = 0;
          prev = curr;
        }

      }
      break;
    }
    
    case '/':
      if (pending_J_closure){
        j_skips = true;
        break;
      }
      else if (pending_diazo){
        fprintf(stderr,"Error: diazo assignment followed by a multiplier is a disallowed bond type\n");
        Fatal(i);
      }
      else{
        fprintf(stderr,"Error: multipliers are not currently supported\n");
        Fatal(i);
      }
      
      
      break;

    default:
      fprintf(stderr, "Error: unallowed character! - alphabet: [A-Z][0-1][&-/' ']\n");
      Fatal(i);
    }

    i++;
    ch = *(++wln_ptr);
  }
  


  if (pending_J_closure)
  {
    fprintf(stderr, "Error: expected 'J' to close ring\n");
    Fatal(len);
  }

  if (pending_locant)
  {
    fprintf(stderr, "Error: expected locant to attach to ring\n");
    Fatal(len);
  }

  if (pending_inline_ring)
  {
    fprintf(stderr, "Error: expected inline ring to be defined\n");
    Fatal(len);
  }

  if (pending_spiro)
  {
    fprintf(stderr, "Error: expected sprio ring to be defined\n");
    Fatal(len);
  }


  if(!AssignCharges(ionic_charges,graph))
    Fatal(len);

    // use this for recursion on multipliers
  return true;
}


/* dump wln tree to a dotvis file */
void WLNDumpToDot(FILE *fp, WLNGraph &graph)
{  
  fprintf(fp, "digraph WLNdigraph {\n");
  fprintf(fp, "  rankdir = LR;\n");
  for (unsigned int i=0; i<=graph.symbol_count;i++)
  {
    WLNSymbol *node = graph.SYMBOLS[i];
    if(!node)
      continue;

    fprintf(fp, "  %d", graph.index_lookup[node]);
    if (node->ch == '*')
      fprintf(fp, "[shape=circle,label=\"%s\"];\n", node->special.c_str());
    else if (node->type == RING)
      fprintf(fp, "[shape=circle,label=\"%c\",color=green];\n", node->ch);
    else{
      if(std::isdigit(node->ch)){
        if (!node->special.empty())
          fprintf(fp, "[shape=circle,label=\"%s\"];\n", node->special.c_str());
        else
          fprintf(fp, "[shape=circle,label=\"%c\"];\n", node->ch);
      } 
      else
        fprintf(fp, "[shape=circle,label=\"%c\"];\n", node->ch);
    }
  
      
    WLNEdge *edge = 0;
    for (edge = node->bonds;edge;edge = edge->nxt){

      WLNSymbol *child = edge->child;
      unsigned int bond_order = edge->order;

      // aromatic
      if (bond_order > 1){
        for (unsigned int k=0;k<bond_order;k++){
          fprintf(fp, "  %d", graph.index_lookup[node]);
          fprintf(fp, " -> ");
          fprintf(fp, "%d [arrowhead=none]\n", graph.index_lookup[child]);
        }
      }
      else{
        fprintf(fp, "  %d", graph.index_lookup[node]);
        fprintf(fp, " -> ");
        fprintf(fp, "%d [arrowhead=none]\n", graph.index_lookup[child]);
      }
    }
  }

  fprintf(fp, "}\n");
}

bool WriteGraph(WLNGraph &graph){
  fprintf(stderr,"Dumping wln graph to wln-graph.dot:\n");
  FILE *fp = 0;
  fp = fopen("wln-graph.dot", "w");
  if (!fp)
  {
    fprintf(stderr, "Error: could not create dump .dot file\n");
    fclose(fp);
    return false;
  }
  else
    WLNDumpToDot(fp,graph);
  
  fclose(fp);
  fp = 0;
  fprintf(stderr,"  dumped\n");
  return true;
}



/**********************************************************************
                         BABEL Mol Functions
**********************************************************************/


// holds all the functions for WLN graph conversion, mol object is assumed ALIVE AT ALL TIMES
// uses old NM functions from previous methods: Copyright (C) NextMove Software 2019-present
struct BabelGraph{

  std::map<unsigned int,OpenBabel::OBAtom*> babel_atom_lookup;

  BabelGraph(){};
  ~BabelGraph(){};


  OpenBabel::OBAtom* NMOBMolNewAtom(OpenBabel::OBMol* mol, unsigned int elem,int charge=0,unsigned int hcount=0)
  {

    OpenBabel::OBAtom* result = mol->NewAtom();
    
    result->SetAtomicNum(elem);
    result->SetImplicitHCount(hcount);
    result->SetFormalCharge(charge);

    return result;
  }

  void NMOBAtomSetAromatic(OpenBabel::OBAtom* atm, bool arom)
  {
    OpenBabel::OBMol* mol = (OpenBabel::OBMol*)atm->GetParent();
    if (mol && !mol->HasAromaticPerceived())
        mol->SetAromaticPerceived();

    atm->SetAromatic(arom);
  }


  bool NMOBMolNewBond(OpenBabel::OBMol* mol,
                      OpenBabel::OBAtom* s,
                      OpenBabel::OBAtom* e,
                      unsigned int order, bool arom)
  {
    
    if(!s || !e){
      fprintf(stderr,"Error: could not find atoms in bond, bond creation impossible\n");
      return false;
    }

    if(opt_debug)
      fprintf(stderr,"  bonding: atoms %3d --> %3d [%d]\n",s->GetIdx(),e->GetIdx(),order);
    

    if (!mol->AddBond(s->GetIdx(), e->GetIdx(), order)){
      fprintf(stderr, "Error: failed to make bond betweens atoms %d --> %d\n",s->GetIdx(),e->GetIdx());
      return false;
    }
        
    OpenBabel::OBBond* bptr = mol->GetBond(mol->NumBonds() - 1);
    if(!bptr){
      fprintf(stderr,"Error: could not re-return bond for checking\n");
      return false;
    }

    if (arom){
      bptr->SetAromatic();
      NMOBAtomSetAromatic(s,true);
      NMOBAtomSetAromatic(e,true);
    }
    return true;
  }


  bool NMOBSanitizeMol(OpenBabel::OBMol* mol)
  {
    
    mol->SetAromaticPerceived(false);

    if(!OBKekulize(mol)){
      fprintf(stderr,"Error: failed on kekulize mol\n");
      return false;
    }
      
    
    // WLN has no inherent stereochemistry, this can be a flag but should be off by default
    mol->SetChiralityPerceived(true);
    
    return true;
  }


  bool ConvertFromWLN(OpenBabel::OBMol* mol,WLNGraph &graph){

    // aromaticity is handled by reducing the edges

    if(opt_debug)
      fprintf(stderr,"Converting wln to obabel mol object: \n");

    // set up atoms
    for (unsigned int i=1; i<=graph.symbol_count;i++){
      WLNSymbol *sym = graph.SYMBOLS[i];

      OpenBabel::OBAtom *atom = 0;

      int charge = 0; 
      unsigned int atomic_num = 0;
      unsigned int hcount = 0;

      switch(sym->ch){

        case 'H':
          atomic_num = 1;
          hcount = 0;
          break; 

        case 'B':
          atomic_num = 5;
          break;

        case '1': // gives me back the 'C' character for unsaturations
        case 'C':
          atomic_num = 6;
          while(sym->num_edges < sym->allowed_edges){
            hcount++;
            sym->num_edges++;
          }
          break;

        case 'X':
          atomic_num = 6;
          hcount = 0; // this must have 4 + methyl
          break;

        case 'Y':
          atomic_num = 6;
          if(sym->type != RING)
            hcount = 1;
          break;

        case 'N':
          atomic_num = 7;
          if(sym->type == RING && sym->allowed_edges > 3)
            sym->allowed_edges = 3;

          while(sym->num_edges < sym->allowed_edges){
            hcount++;
            sym->num_edges++;
          }
          break;

        case 'M':
          atomic_num = 7;
          hcount = 1;
          break;

        case 'Z':
          atomic_num = 7; 
          hcount = 2;
          break;

        case 'K':
          atomic_num = 7;
          charge = 1; 
          hcount = 0;
          break;

        case 'O':
          atomic_num = 8;
          if(!sym->num_edges)
            charge = -1;
          break;

        case 'Q':
          atomic_num = 8;
          hcount = 1;
          break;

        case 'F':
          atomic_num = 9;
          if(!sym->num_edges)
            charge = -1;
          break;
        
        case 'P':
          atomic_num = 15;
          while(sym->num_edges < 3){
            hcount++;
            sym->num_edges++;
          }
          break;
        
        case 'S':
          atomic_num = 16;
          while(sym->num_edges < 2){
            hcount++;
            sym->num_edges++;
          }
          break;

        case 'G':
          atomic_num = 17;
          if(!sym->num_edges)
            charge = -1;
          break;

        case 'E':
          atomic_num = 35;
          if(!sym->num_edges)
            charge = -1;
          break;

        case 'I':
          atomic_num = 53;
          if(!sym->num_edges)
            charge = -1;
          break;
      
        case '*':
          atomic_num = special_element_atm(sym->special);
          break;

        default:
          fprintf(stderr,"Error: unrecognised WLNSymbol* char in obabel mol build - %c\n",sym->ch);
          return false;
      }

      // ionic notation - overrides any given formal charge
      if(graph.charge_additions[sym]){
        charge = graph.charge_additions[sym];
      }

      atom = NMOBMolNewAtom(mol,atomic_num,charge,hcount);
      if(!atom){
        fprintf(stderr,"Error: formation of obabel atom object\n");
        return false;
      }

      if(sym->type == RING)
        atom->SetInRing();

      babel_atom_lookup[graph.index_lookup[sym]] = atom;
      if(opt_debug)
        fprintf(stderr,"  created: atom[%d] - atomic num(%d), charge(%d)\n",atom->GetIdx(),atomic_num,charge);
    
    }

    // set bonds
    for(unsigned int i=1;i<=graph.symbol_count;i++){
      WLNSymbol *parent = graph.SYMBOLS[i];

      unsigned int parent_id = graph.index_lookup[parent];
      OpenBabel::OBAtom *par_atom = babel_atom_lookup[parent_id];

      WLNEdge *e = 0;
      
      if(parent->bonds){
        
        for (e = parent->bonds;e;e = e->nxt){
          
          WLNSymbol *child = e->child;

          unsigned int bond_order = e->order;  
    
          unsigned int child_id = graph.index_lookup[child];
          OpenBabel::OBAtom *chi_atom = babel_atom_lookup[child_id];
          if(!NMOBMolNewBond(mol,par_atom,chi_atom,bond_order,false))
            return false;
        }

      }
    }

    return true;
  }


};



/**********************************************************************
                         API FUNCTION
**********************************************************************/


bool ReadWLN(const char *ptr, OpenBabel::OBMol* mol)
{   
  if(!ptr){
    fprintf(stderr,"Error: could not read wln string pointer\n");
    return false;
  }
  else 
    wln_string = ptr; 

  WLNGraph wln_graph;
  BabelGraph obabel; 

  bool state = true;

  if(state)
    state = ParseWLNString(ptr,wln_graph);
  
  // create an optional wln dotfile
  if (opt_wln2dot)
    WriteGraph(wln_graph);
  
  if(state)
    state = ExpandWLNSymbols(wln_graph);

  if(state)
    state = obabel.ConvertFromWLN(mol,wln_graph);

  if(state)
    state = obabel.NMOBSanitizeMol(mol);

  return state;
}
