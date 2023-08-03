/**********************************************************************
Copyright (C) 2019 by NextMove Software

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

#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/kekulize.h>

#include <vector>

#define PENDING_NONE   0
#define PENDING_DEPROT 1
#define PENDING_METHYL 2

#define STACK_NORMAL  0
#define STACK_METHYL  1
#define STACK_POLY    2
#define STACK_RING    3

#define STACK_SHIFT   2
#define STACK_MASK    3


unsigned int debug_wln_read = 0;



OpenBabel::OBAtom* NMOBMolNewAtom(OpenBabel::OBMol* mol, unsigned int elem)
{
  OpenBabel::OBAtom* result = mol->NewAtom();
  result->SetAtomicNum(elem);
  return result;
}


OpenBabel::OBBond* NMOBMolNewBond(OpenBabel::OBMol* mol,
  OpenBabel::OBAtom* beg,
  OpenBabel::OBAtom* end,
  unsigned int order, bool arom)
{
  if (!mol->AddBond(beg->GetIdx(), end->GetIdx(), order))
    return nullptr;
  OpenBabel::OBBond* bptr = mol->GetBond(mol->NumBonds() - 1);
  if (arom)
    bptr->SetAromatic();
  return bptr;
}


void NMOBAtomSetAromatic(OpenBabel::OBAtom* atm, bool arom)
{
  OpenBabel::OBMol* mol = (OpenBabel::OBMol*)atm->GetParent();
  if (mol && !mol->HasAromaticPerceived())
    mol->SetAromaticPerceived();

  atm->SetAromatic(arom);
}


bool NMOBSanitizeMol(OpenBabel::OBMol* mol)
{
  if (!OBKekulize(mol))
    return false;
  mol->SetAromaticPerceived(false);
  return true;
}


struct WLNParser {
  OpenBabel::OBMol* mol;
  const char *orig;
  const char *ptr;

  std::vector<unsigned int> stack;
  std::vector<std::vector<OpenBabel::OBAtom*> > rings;
  std::vector<OpenBabel::OBAtom*> atoms;
  unsigned int pending;
  unsigned int state;
  unsigned int order;
  unsigned int slash;
  OpenBabel::OBAtom* prev;

  WLNParser(const char *wln, OpenBabel::OBMol* m) {
    mol = m;
    orig = wln;
    ptr = wln;

    pending = PENDING_NONE;
    prev = nullptr;
    slash = 0;
    order = 0;
    state = 0;
  }

  bool error() {
    fprintf(stderr,"Error: Character %c in %s\n",*ptr,orig);
    unsigned int len = (unsigned int)(ptr-orig)+22;
    for (unsigned int i=0; i<len; i++)
      fputc(' ',stderr);
    fprintf(stderr,"^\n");
    return false;
  }

  OpenBabel::OBAtom* atom(unsigned int elem, unsigned int hcount) {
    OpenBabel::OBAtom* aptr = NMOBMolNewAtom(mol,elem);
    aptr->SetImplicitHCount(hcount);
    atoms.push_back(aptr);
    return aptr;
  }

  void dec_h(OpenBabel::OBAtom* aptr, unsigned int count) {
    unsigned int hcount = aptr->GetImplicitHCount();
    if (hcount > count)
      aptr->SetImplicitHCount(hcount-count);
    else if (hcount)
      aptr->SetImplicitHCount(0);
  }

  void dec_q(OpenBabel::OBAtom* aptr, unsigned int count) {
    int charge = aptr->GetFormalCharge();
    aptr->SetFormalCharge(charge-count);
  }

  void fuse(OpenBabel::OBAtom* src, OpenBabel::OBAtom* dst, unsigned int order)
  {
    dec_h(src,order);
    dec_h(dst,order);
    NMOBMolNewBond(mol,src,dst,order,false);
  }

  void next(OpenBabel::OBAtom* aptr, unsigned int bo) {
    pending = PENDING_DEPROT;
    prev = aptr;
    order = bo;
    if (bo) {
      if (state == 0)
        state = 1;
    } else state = 2;
  }

  void term() {
    if (pending == PENDING_DEPROT) {
      if (order) {
        dec_q(prev,order);
        dec_h(prev,order);
      }
    } else if (pending == PENDING_METHYL) {
      if (order == 1) {
        OpenBabel::OBAtom* temp = atom(6,4);
        fuse(prev,temp,1);
      }
    }
  }

  void push() {
    unsigned int idx = (unsigned int)(atoms.size())-1;
    stack.push_back((idx<<STACK_SHIFT)+STACK_NORMAL);
  }

  void push_methyl() {
    unsigned int idx = (unsigned int)(atoms.size())-1;
    stack.push_back((idx<<STACK_SHIFT)+STACK_METHYL);
  }

  void push_poly() {
    unsigned int idx = (unsigned int)(atoms.size())-1;
    stack.push_back((idx<<STACK_SHIFT)+STACK_POLY);
  }

  void pop_common() {
    unsigned int top = stack.back();
    switch (top & STACK_MASK) {
    case STACK_RING:
      stack.pop_back();
      rings.pop_back();
      state = 2;
      order = 0;
      if (!stack.empty() && stack.back()!=STACK_RING)
        pop_common();
      return;
    case STACK_NORMAL:
      pending = PENDING_DEPROT;
      stack.pop_back();
      break;
    case STACK_POLY:
      pending = PENDING_NONE;
      break;
    case STACK_METHYL:
      pending = PENDING_METHYL;
      stack.pop_back();
      break;
    }
    prev = atoms[top>>STACK_SHIFT];
    order = 1;
    state = 1;
  }

  bool pop() {
    if (stack.empty())
      return error();
    pop_common();
    return true;
  }

  void may_pop() {
    if (stack.empty()|| stack.back()==STACK_RING) {
      state = 2;
      order = 0;
    } else pop_common();
  }

  bool double_pop() {
    if (stack.empty())
      return false;
    unsigned int top = stack.back();
    if ((top&STACK_MASK) == STACK_POLY) {
      stack.pop_back();
    } else if (!pop())
      return false;
    return pop();
  }

  void drain() {
    term();
    while (!stack.empty()) {
      unsigned int top = stack.back();
      if ((top&STACK_MASK) != STACK_POLY) {
        pop_common();
        term();
      } else stack.pop_back();
    }
  }

  /* Two atoms */
  OpenBabel::OBAtom* cyano() {
    OpenBabel::OBAtom* aptr = atom(6,1);
    OpenBabel::OBAtom* temp = atom(7,0);
    NMOBMolNewBond(mol,aptr,temp,3,false);
    return aptr;
  }

  OpenBabel::OBAtom* carbonyl() {
    OpenBabel::OBAtom* aptr = atom(6,2);
    OpenBabel::OBAtom* temp = atom(8,0);
    NMOBMolNewBond(mol,aptr,temp,2,false);
    return aptr;
  }

  void carbon_monoxide() {
    OpenBabel::OBAtom* aptr = atom(6,0);
    OpenBabel::OBAtom* temp = atom(8,0);
    NMOBMolNewBond(mol,aptr,temp,3,false);
    aptr->SetFormalCharge(-1);
    temp->SetFormalCharge(1);
  }

  OpenBabel::OBAtom* thiocarbonyl() {
    OpenBabel::OBAtom* aptr = atom(6,2);
    OpenBabel::OBAtom* temp = atom(16,0);
    NMOBMolNewBond(mol,aptr,temp,2,false);
    return aptr;
  }

  OpenBabel::OBAtom* diazo() {
    OpenBabel::OBAtom* aptr = atom(7,2);
    OpenBabel::OBAtom* temp = atom(7,0);
    NMOBMolNewBond(mol,aptr,temp,2,false);
    temp->SetFormalCharge(-1);
    aptr->SetFormalCharge(1);
    return aptr;
  }

  OpenBabel::OBAtom* diazonio() {
    OpenBabel::OBAtom* aptr = atom(7,1);
    OpenBabel::OBAtom* temp = atom(7,0);
    NMOBMolNewBond(mol,aptr,temp,3,false);
    aptr->SetFormalCharge(1);
    return aptr;
  }

  OpenBabel::OBAtom* nitroso() {
    OpenBabel::OBAtom* aptr = atom(7,1);
    OpenBabel::OBAtom* temp = atom(8,0);
    NMOBMolNewBond(mol,aptr,temp,2,false);
    return aptr;
  }

  OpenBabel::OBAtom* sulfinyl() {
    OpenBabel::OBAtom* aptr = atom(16,2);
    OpenBabel::OBAtom* temp = atom(8,0);
    NMOBMolNewBond(mol,aptr,temp,2,false);
    return aptr;
  }

  /* Three atoms */
  OpenBabel::OBAtom* nitro() {
    OpenBabel::OBAtom* aptr = atom(7,1);
    aptr->SetFormalCharge(1);
    OpenBabel::OBAtom* temp = atom(8,0);
    NMOBMolNewBond(mol,aptr,temp,2,false);
    temp = atom(8,0);
    NMOBMolNewBond(mol,aptr,temp,1,false);
    temp->SetFormalCharge(-1);
    return aptr;
  }

  OpenBabel::OBAtom* sulfonyl() {
    OpenBabel::OBAtom* aptr = atom(16,2);
    OpenBabel::OBAtom* temp = atom(8,0);
    NMOBMolNewBond(mol,aptr,temp,2,false);
    temp = atom(8,0);
    NMOBMolNewBond(mol,aptr,temp,2,false);
    return aptr;
  }

  /* Four atoms */
  OpenBabel::OBAtom* sulfamoyl() {
    OpenBabel::OBAtom* aptr = atom(16,1);
    OpenBabel::OBAtom* temp = atom(7,2);
    NMOBMolNewBond(mol,aptr,temp,1,false);
    temp = atom(8,0);
    NMOBMolNewBond(mol,aptr,temp,2,false);
    temp = atom(8,0);
    NMOBMolNewBond(mol,aptr,temp,2,false);
    return aptr;
  }

  OpenBabel::OBAtom* sulfonato() {
    OpenBabel::OBAtom* aptr = atom(16,1);
    OpenBabel::OBAtom* temp = atom(8,0);
    NMOBMolNewBond(mol,aptr,temp,2,false);
    temp = atom(8,0);
    NMOBMolNewBond(mol,aptr,temp,2,false);
    temp = atom(8,0);
    NMOBMolNewBond(mol,aptr,temp,1,false);
    temp->SetFormalCharge(-1);
    return aptr;
  }

  /* Five atoms */
  OpenBabel::OBAtom* perchlorate() {
    OpenBabel::OBAtom* aptr = atom(17,0);
    OpenBabel::OBAtom* temp = atom(8,0);
    NMOBMolNewBond(mol,aptr,temp,2,false);
    temp = atom(8,0);
    NMOBMolNewBond(mol,aptr,temp,2,false);
    temp = atom(8,0);
    NMOBMolNewBond(mol,aptr,temp,2,false);
    temp = atom(8,0);
    NMOBMolNewBond(mol,aptr,temp,1,false);
    temp->SetFormalCharge(-1);
    return aptr;
  }

  void new_ring(std::vector<OpenBabel::OBAtom*> &ring, unsigned int size) {
    for (unsigned int i=0; i<size; i++) {
      OpenBabel::OBAtom* aptr = atom(6,1);
      NMOBAtomSetAromatic(aptr,true);
      ring.push_back(aptr);
    }
    // push_ring(ring);
    stack.push_back(STACK_RING);
    rings.push_back(ring);
  }

  void new_cycle(std::vector<OpenBabel::OBAtom*> &ring, unsigned int size) {
    new_ring(ring,size);
    for (unsigned int i=0; i<size-1; i++)
      // NMOBMolNewBond(mol,ring[i],ring[i+1],1+(i&1),true);
      NMOBMolNewBond(mol,ring[i],ring[i+1],1,true);
    NMOBMolNewBond(mol,ring[size-1],ring[0],1,true);
  }

  void bridge(OpenBabel::OBAtom* src, OpenBabel::OBAtom* dst) {
    src->SetImplicitHCount(0);
    dst->SetImplicitHCount(0);
    NMOBMolNewBond(mol,src,dst,1,true);
  }

  void new_ring36(std::vector<OpenBabel::OBAtom*> &ring) {
    new_cycle(ring,7);
    bridge(ring[0],ring[2]);
  }

  void new_ring45(std::vector<OpenBabel::OBAtom*> &ring) {
    new_cycle(ring,7);
    bridge(ring[0],ring[3]);
  }

  void new_ring55(std::vector<OpenBabel::OBAtom*> &ring) {
    new_cycle(ring,8);
    bridge(ring[0],ring[4]);
  }

  void new_ring56(std::vector<OpenBabel::OBAtom*> &ring) {
    new_cycle(ring,9);
    bridge(ring[0],ring[4]);
  }

  void new_ring57(std::vector<OpenBabel::OBAtom*> &ring) {
    new_cycle(ring,10);
    bridge(ring[0],ring[4]);
  }

  void new_ring66(std::vector<OpenBabel::OBAtom*> &ring) {
    new_cycle(ring,10);
    bridge(ring[0],ring[5]);
  }

  void new_ring67(std::vector<OpenBabel::OBAtom*> &ring) {
    new_cycle(ring,11);
    bridge(ring[0],ring[5]);
  }

  OpenBabel::OBAtom* benzene() {
    std::vector<OpenBabel::OBAtom*> ring;
    new_cycle(ring,6);
    return ring[0];
  }

  bool term1(OpenBabel::OBAtom* aptr) {
    if (state == 0) {
      next(aptr,1);
    } else if (order == 1) {
      fuse(prev,aptr,1);
      may_pop();
    } else return error();
    return true;
  }

  bool degree1(OpenBabel::OBAtom* aptr) {
    if (state == 0) {
      next(aptr,1);
    } else if (order == 1) {
      fuse(prev,aptr,1);
      next(aptr,0);
    } else return error();
    return true;
  }

  bool degree2(OpenBabel::OBAtom* aptr) {
    if (order == 1) {
      fuse(prev,aptr,1);
      next(aptr,1);
    } else return error();
    return true;
  }

  bool alkane() {
    unsigned int len = (*ptr++ - '0');
    while(*ptr>='0' && *ptr<='9')
      len = 10*len + (*ptr++ - '0');
    for (unsigned int i=0; i<len; i++) {
      OpenBabel::OBAtom* temp = atom(6,4);
      if (order)
        fuse(prev,temp,order);
      prev = temp;
      order = 1;
    }
    next(prev,order);
    pending = PENDING_NONE;
    return true;
  }

  bool poly(unsigned int elem) {
    if (state == 0) {
      prev = atom(elem,0);
      push_poly();
      state = 1;
    } else if (state == 1) {
      OpenBabel::OBAtom* temp = atom(elem,0);
      fuse(prev,temp,order);
      pending = PENDING_NONE;
      push_poly();
      prev = temp;
    } else return error();
    pending = PENDING_NONE;
    order = 1;
    return true;
  }

  void unsaturate(OpenBabel::OBAtom* src, OpenBabel::OBAtom* dst) {
    OpenBabel::OBBond* bptr = src->GetBond(dst);
    if (bptr)
      bptr->SetBondOrder(2);
    NMOBAtomSetAromatic(src,false);
    NMOBAtomSetAromatic(dst,false);
  }

  void unsaturate2(OpenBabel::OBAtom* src, OpenBabel::OBAtom* dst) {
    OpenBabel::OBBond* bptr = src->GetBond(dst);
    if (bptr)
      bptr->SetBondOrder(3);
    NMOBAtomSetAromatic(src,false);
    NMOBAtomSetAromatic(dst,false);
    src->SetImplicitHCount(0);
    dst->SetImplicitHCount(0);
  }

  bool atend(const char *ptr) {
    switch (*ptr) {
    case '\0':
    case '\t':
    case ' ':
    case '&':
      return true;
    }
    return false;
  }

  bool atend() {
    return atend(ptr);
  }

  // ptr points to 'L' or 'T'
  bool parse_ring() {
    unsigned int size = 0;
    if (ptr[1]>='3' && ptr[1]<='9') {
      size = ptr[1]-'0';
      ptr += 2;
    } else if (ptr[1]=='-' &&
               ptr[2]>='1' && ptr[2]<='9' &&
               ptr[3]>='0' && ptr[3]<='9' &&
               ptr[4]=='-') {
      size = 10*(ptr[2]-'0') + (ptr[3]-'0');
      ptr += 5;
    } else return error();

    bool done = false;
    std::vector<OpenBabel::OBAtom*> ring;
    switch (size) {
    case 3:
      if (ptr[0]=='6') {
        new_ring36(ring);
        done = true;
        size = 7;
        ptr++;
      }
      break;
    case 4:
      if (ptr[0]=='5') {
        new_ring45(ring);
        done = true;
        size = 7;
        ptr++;
      }
      break;
    case 5:
      if (ptr[0]=='5') {
        new_ring55(ring);
        done  = true;
        size = 8;
        ptr++;
      } else if (ptr[0]=='6') {
        new_ring56(ring);
        done = true;
        size = 9;
        ptr++;
      } else if (ptr[0]=='7') {
        new_ring57(ring);
        done = true;
        size = 10;
        ptr++;
      }
      break;
    case 6:
      if (ptr[0]=='6') {
        new_ring66(ring);
        done = true;
        size = 10;
        ptr++;
      } else if (ptr[0]=='7') {
        new_ring67(ring);
        done = true;
        size = 11;
        ptr++;
      }
      break;
    }

    if (!done)
      new_cycle(ring,size);

    if (debug_wln_read)
      printf("DEBUG: ring size=%u ptr=%s\n",size,ptr);

    unsigned int loc = 0;
    unsigned int elem = 0;
    unsigned int hcount = 0;
    OpenBabel::OBAtom* aptr;

    for (;;) {
    if (debug_wln_read)
      printf("DEBUG: loc=%u ptr=%s\n",loc,ptr);
    switch(*ptr) {
    case 'H':
      if (loc >= size)
        return error();
      aptr = ring[loc];
      aptr->SetImplicitHCount(2);
      NMOBAtomSetAromatic(aptr,false);
      ptr++;
      if (*ptr != 'J')
        return error();
      ptr++;
      return true;

    case 'J':
      ptr++;
      return true;

    case 'K':
      if (loc >= size)
        return error();
      aptr = ring[loc];
      aptr->SetAtomicNum(7);
      aptr->SetFormalCharge(1);
      aptr->SetImplicitHCount(0);
      loc++;
      ptr++;
      continue;

    case 'M':
      if (loc >= size)
        return error();
      aptr = ring[loc];
      aptr->SetAtomicNum(7);
      NMOBAtomSetAromatic(aptr,false);
      aptr->SetImplicitHCount(1);
      loc++;
      ptr++;
      continue;

    case 'N':
      if (loc >= size)
        return error();
      aptr = ring[loc];
      aptr->SetAtomicNum(7);
      aptr->SetImplicitHCount(0);
      loc++;
      ptr++;
      continue;

    case 'O':
      if (loc >= size)
        return error();
      aptr = ring[loc];
      aptr->SetAtomicNum(8);
      NMOBAtomSetAromatic(aptr,false);
      aptr->SetImplicitHCount(0);
      loc++;
      ptr++;
      continue;

    case 'S':
      if (loc >= size)
        return error();
      aptr = ring[loc];
      aptr->SetAtomicNum(16);
      aptr->SetImplicitHCount(0);
      NMOBAtomSetAromatic(aptr,false);
      loc++;
      ptr++;
      if (*ptr == 'W') {
        NMOBMolNewBond(mol,aptr,atom(8,0),2,false);
        NMOBMolNewBond(mol,aptr,atom(8,0),2,false);
        ptr++;
      }
      continue;

    case 'T':
      for (unsigned int i=0; i<size; i++) {
        aptr = ring[i];
        if (aptr->IsAromatic() &&
            aptr->GetAtomicNum()==6) {
          unsigned int hcount = aptr->GetImplicitHCount();
          aptr->SetImplicitHCount(hcount+1);
          NMOBAtomSetAromatic(aptr,false);
        }
      }
      ptr++;
      if (*ptr != 'J')
        return error();
      ptr++;
      return true;

    case 'U':
      if (loc+1 >= size)
        return error();
      if (ptr[1]=='U') {
        unsaturate2(ring[loc],ring[loc+1]);
        ptr += 2;
      } else {
        unsaturate(ring[loc],ring[loc+1]);
        ptr++;
      }
      continue;

    case 'V':
      if (loc >= size)
        return error();
      aptr = ring[loc];
      fuse(aptr,atom(8,0),2);
      NMOBAtomSetAromatic(aptr,false);
      loc++;
      ptr++;
      continue;

    case 'X':
    case 'Y':
      if (loc >= size)
        return error();
      aptr = ring[loc];
      aptr->SetImplicitHCount(2);
      NMOBAtomSetAromatic(aptr,false);
      loc++;
      ptr++;
      continue;

    case ' ':
      if (ptr[1]>='A' && ptr[1]<='Z') {
        loc = ptr[1]-'A';
        if (loc >= size)
          return error();
        ptr += 2;
        // Check the following character
        switch (*ptr) {
        case 'H':
        case 'K':
        case 'M':
        case 'N':
        case 'O':
        case 'S':
        case 'U':
        case 'V':
        case 'X':
        case 'Y':
        case '-':
          break;
        default:
          return error();
        }
        continue;
      }
      break;
    case '-':
      elem = 0;
      switch (ptr[1]) {
      case 'A':
        if (ptr[2]=='S' && ptr[3]=='-') {  // -AS-
          elem = 33;
          ptr += 4;
        }
        break;
      case 'G':
        if (ptr[2]=='E' && ptr[3]=='-') {  // -GE-
          elem = 32;
          ptr += 4;
        }
        break;
      case 'H':
        if (ptr[2]=='G' && ptr[3]=='-') {  // -HG-
          elem = 80;
          ptr += 4;
        }
        break;
      case 'S':
        if (ptr[2]=='B' && ptr[3]=='-') {  // -SB-
          elem = 51;
          ptr += 4;
        } else
        if (ptr[2]=='E' && ptr[3]=='-') {  // -SE-
          elem = 34;
          ptr += 4;
        } else
        if (ptr[2]=='I' && ptr[3]=='-') {  // -SI-
          elem = 14;
          ptr += 4;
        } else
        if (ptr[2]=='N' && ptr[3]=='-') {  // -SN-
          elem = 50;
          ptr += 4;
        }
        break;
      }
      if (!elem)
        return error();
      if (loc >= size)
        return error();
      hcount = 0;
      if (*ptr == 'H') {
        if (ptr[1] == 'H') {
          hcount = 2;
          ptr += 2;
        }  else {
          hcount = 1;
          ptr++;
        }
      }
      aptr = ring[loc];
      aptr->SetAtomicNum(elem);
      aptr->SetImplicitHCount(hcount);
      NMOBAtomSetAromatic(aptr,false);
      loc++;
      continue;
    }
    return error();
    }
  }

  bool parse() {
    ptr = orig;
    if (ptr[0]=='W' && ptr[1]=='L' && ptr[2]=='N' &&
        ptr[3]==':' && ptr[4]==' ')
      ptr += 5;
    for(;;) {
    if (debug_wln_read)
      printf("DEBUG: state=%u order=%u ptr=%s\n",state,order,ptr);
    switch (*ptr) {
    case '\t':
      mol->SetTitle(ptr+1);
      /* fall through */
    case '\0':
      if (state == 0)
        return error();
      drain();
      return true;

    case 'B':
      if (order==1) {
        OpenBabel::OBAtom* temp = atom(5,3);
        fuse(prev,temp,order);
        pending = PENDING_DEPROT;
        push();
        prev = temp;
        order = 1;
        ptr++;
        continue;
      }
      break;

    case 'C':
      if (order == 3) {
        OpenBabel::OBAtom* temp = atom(6,4);
        fuse(prev,temp,3);
        next(temp,1);
        ptr++;
        continue;
      }
      if (ptr[1]=='N') {  // CN
        if (order == 1) {
          if (!degree1(cyano()))
            return false;
          ptr += 2;
          continue;
        }
      }
      if (ptr[1]=='U') {  // CU
        if (order == 2) {
          OpenBabel::OBAtom* temp = atom(6,0);
          fuse(prev,temp,2);
          next(temp,2);
          ptr += 2;
          continue;
        }
      }
      break;

    case 'E':
      if (!term1(atom(35,1)))
        return false;
      ptr++;
      continue;

    case 'F':
      if (!term1(atom(9,1)))
        return false;
      ptr++;
      continue;

    case 'G':
      if (!term1(atom(17,1)))
        return false;
      ptr++;
      continue;

    case 'H':
      if (order == 1) {
        may_pop();
        ptr++;
        continue;
      }
      break;

    case 'I':
      if (!term1(atom(53,1)))
        return false;
      ptr++;
      continue;

    case 'K':
      if (order == 1) {
        OpenBabel::OBAtom* temp = atom(7,4);
        temp->SetFormalCharge(1);
        fuse(prev,temp,order);
        pending = PENDING_METHYL;
        push_methyl();
        push_methyl();
        prev = temp;
        order = 1;
        ptr++;
        continue;
      }
      break;

    case 'L':
      if (state == 0) {
        if (!parse_ring())
          return false;
        state = 2;
        order = 0;
        continue;
      }
      break;

    case 'M':
      if (order == 2) {
        OpenBabel::OBAtom* temp = atom(7,3);
        fuse(prev,temp,2);
        next(temp,0);
        ptr++;
        continue;
      }
      if (ptr[1]=='U') {  // MU
        if (state == 0) {
          next(atom(7,3),2);
          ptr += 2;
          continue;
        }
      }
      if (!degree2(atom(7,3)))
        return false;
      ptr++;
      continue;

    case 'N':
      if (order == 3) {
        OpenBabel::OBAtom* temp = atom(7,0);
        fuse(prev,temp,3);
        next(temp,0);
        ptr++;
        continue;
      }
      if (order == 2) {
        OpenBabel::OBAtom* temp = atom(7,3);
        fuse(prev,temp,2);
        next(temp,1);
        ptr++;
        continue;
      }
      if (ptr[1]=='C') {  // NC
        if (state == 0) {
          if (!degree1(cyano()))
            return false;
          ptr += 2;
          continue;
        }
      }
      if (ptr[1]=='N') {  // NN
        if (ptr[2]=='U') {  // NNU
          if (state == 0) {
            next(diazo(),2);
            ptr += 3;
            continue;
          }
        }
        if (!degree1(diazonio()))
          return false;
        ptr += 2;
        continue;
      }
      if (ptr[1]=='O') {  // NO
        if (order == 1) {
          if (!degree1(nitroso()))
            return false;
          ptr += 2;
          continue;
        }
      }
      if (ptr[1]=='U') {  // NU
        if (order == 1) {
          OpenBabel::OBAtom* temp = atom(7,0);
          fuse(prev,temp,1);
          next(temp,2);
          ptr += 2;
          continue;
        }
      }
      if (ptr[1]=='W') {  // NW
        if (order == 1) {
          if (!term1(nitro()))
            return false;
          ptr += 2;
          continue;
        }
      }
      if (order == 1) {
        OpenBabel::OBAtom* temp = atom(7,3);
        fuse(prev,temp,order);
        pending = PENDING_DEPROT;
        prev = temp;
        order = 1;
        push();
        ptr++;
        continue;
      }
      break;

    case 'O':
      if (order == 2) {
        OpenBabel::OBAtom* temp = atom(8,0);
        fuse(prev,temp,2);
        next(temp,0);
        ptr++;
        continue;
      }
      if (ptr[1]=='C') {  // OC
        if (state == 0) {
          if (atend(ptr+2)) {
            carbon_monoxide();
            state = 2;
          } else
            next(carbonyl(),2);
          ptr += 2;
          continue;
        }
      }
      if (ptr[1]=='K') {  // OK
        if (state == 0) {
          prev = atom(8,0);
          prev->SetFormalCharge(-1);
          next(prev,1);
          ptr++;
          continue;
        }
      }
      if (ptr[1]=='N') {  // ON
        if (state == 0) {
          next(nitroso(),1);
          ptr += 2;
          continue;
        }
      }
      if (ptr[1]=='S' && ptr[2]=='W') {  // OSW
        if (state == 0) {
          next(sulfonato(),1);
          ptr += 3;
          continue;
        }
      }
      if (ptr[1]=='S') {  // OS
        if (state == 0) {
          prev = atom(16,2);
          push();
          NMOBMolNewBond(mol,prev,atom(8,0),2,false);
          pending = PENDING_DEPROT;
          state = 1;
          order = 1;
          ptr += 2;
          continue;
        }
      }
      if (ptr[1]=='V') {  // OV
        if (state == 0) {
          OpenBabel::OBAtom* temp = atom(8,0);
          temp->SetFormalCharge(-1);
          next(temp,1);
          ptr++;
          continue;
        }
      }
      // Handle O-AS-R and friends
      if (state == 0) {
        next(atom(8,2),2);
        ptr++;
        continue;
      }
      if (!degree2(atom(8,2)))
        return false;
      ptr++;
      continue;

    case 'P':
      if (state == 0) {
        if (!term1(atom(15,3)))
          return false;
      } else if (state == 1) {
        OpenBabel::OBAtom* temp = atom(15,3);
        fuse(prev,temp,order);
        next(temp,1);
      } else break;
      pending = PENDING_NONE;
      push_poly();
      ptr++;
      if (ptr[0]=='O' && ptr[1]=='&') {  // PO&
        OpenBabel::OBAtom* temp = atom(8,0);
        NMOBMolNewBond(mol,prev,temp,2,false);
        ptr += 2;
      } else if (ptr[0]=='S' && ptr[1]=='&') {  // PS&
        OpenBabel::OBAtom* temp = atom(16,0);
        NMOBMolNewBond(mol,prev,temp,2,false);
        ptr += 2;
      }
      continue;

    case 'Q':
      if (!term1(atom(8,2)))
        return false;
      ptr++;
      continue;

    case 'R':
      if (!degree1(benzene()))
        return false;
      ptr++;
      continue;

    case 'S':
      if (order == 2) {
        OpenBabel::OBAtom* temp = atom(16,0);
        fuse(prev,temp,2);
        next(temp,0);
        ptr++;
        continue;
      }
      if (ptr[1]=='C') {
        if (state == 0) {
          next(thiocarbonyl(),2);
          ptr += 2;
          continue;
        }
      }
      if (ptr[1]=='H') {
        if (!degree1(atom(16,2)))
          return false;
        ptr += 2;
        continue;
      }
      if (ptr[1]=='O' && ptr[2]=='&') {  // SO&
        if (!degree2(sulfinyl()))
          return false;
        ptr += 3;
        continue;
      }
      if (ptr[1]=='U') {
        if (state == 0) {
          next(atom(16,2),2);
          ptr += 2;
          continue;
        }
      }
      if (ptr[1]=='W') {
        if (!degree2(sulfonyl()))
          return false;
        ptr += 2;
        continue;
      }
      if (ptr[1]=='Z' && ptr[2]=='W') {  // SZW
        if (order == 1) {
          if (!term1(sulfamoyl()))
            return false;
          ptr += 3;
          continue;
        }
      }
      if (state == 0) {
        next(atom(16,2),2);
        ptr++;
        continue;
      }
      if (!degree2(atom(16,2)))
        return false;
      ptr++;
      continue;

    case 'T':
      if (state == 0) {
        if (!parse_ring())
          return false;
        state = 2;
        order = 0;
        continue;
      }
      break;

    case 'U':
      if (order != 1)
        return error();
      if (ptr[1]=='U') {
        order = 3;
        ptr += 2;
      } else {
        order = 2;
        ptr++;
      }
      continue;

    case 'V':
      if (ptr[1]=='H') {
        if (state == 0) {
          /* degree1 */
          next(carbonyl(),1);
          ptr += 2;
          continue;
        }
      }
      if (!degree2(carbonyl()))
        return false;
      ptr++;
      continue;

    case 'W':
      if (ptr[1]=='G' && ptr[2]=='W') {  // WGW
        if (state == 0) {
          next(perchlorate(),0);
          ptr += 3;
          continue;
        }
      }
      if (ptr[1]=='N') {  // WN
        // degree1(nitro(1))
        if (state == 0) {
          next(nitro(),1);
          ptr += 2;
          continue;
        }
      }
      if (ptr[1]=='P' && ptr[2]=='H') {  // WPH
        if (state == 0) {
          prev = atom(15,1);
          OpenBabel::OBAtom* temp = atom(8,1);
          NMOBMolNewBond(mol,prev,temp,1,false);
          temp = atom(8,0);
          temp->SetFormalCharge(-1);
          NMOBMolNewBond(mol,prev,temp,1,false);
          next(prev,1);
          ptr += 3;
          continue;
        }
      }
      if (ptr[1]=='S' && ptr[2]=='Q') {  // WSQ
        if (state == 0) {
          prev = atom(16,1);
          NMOBMolNewBond(mol,prev,atom(8,0),2,false);
          NMOBMolNewBond(mol,prev,atom(8,0),2,false);
          NMOBMolNewBond(mol,prev,atom(8,1),1,false);
          pending = PENDING_DEPROT;
          state = 1;
          order = 1;
          ptr += 3;
          continue;
        }
      }
      if (ptr[1]=='S') {  // WS
        if (state == 0) {
          prev = atom(16,2);
          push();
          NMOBMolNewBond(mol,prev,atom(8,0),2,false);
          NMOBMolNewBond(mol,prev,atom(8,0),2,false);
          pending = PENDING_DEPROT;
          state = 1;
          order = 1;
          ptr += 2;
          continue;
        }
      }
      break;

    case 'X':
      if (order==1 || order==2) {
        OpenBabel::OBAtom* temp = atom(6,4);
        fuse(prev,temp,order);
        pending = PENDING_METHYL;
        push_methyl();
        push_methyl();
        prev = temp;
        order = 1;
        ptr++;
        continue;
      }
      break;

    case 'Y':
      if (order==1 || order==2) {
        OpenBabel::OBAtom* temp = atom(6,4);
        fuse(prev,temp,order);
        pending = PENDING_METHYL;
        push_methyl();
        prev = temp;
        order = 1;
        ptr++;
        continue;
      }
      break;

    case 'Z':
      if (!term1(atom(7,3)))
        return false;
      ptr++;
      continue;

    case '-':
      switch (ptr[1]) {
      case 'A':
        if (ptr[2]=='G' && ptr[3]=='-') {  // -AG-
          if (state == 0 && atend(ptr+4)) {
            prev = atom(47,0);
            prev->SetFormalCharge(1);
            state = 2;
            ptr += 4;
            continue;
          }
          if (!poly(47))
            return false;
          ptr += 4;
          continue;
        }
        if (ptr[2]=='S' && ptr[3]=='-') {  // -AS-
          if (!poly(33))
            return false;
          ptr += 4;
          continue;
        }
        break;
      case 'B':
        if (ptr[2]=='A' && ptr[3]=='-') {  // -BA-
          if (!poly(56))
            return false;
          ptr += 4;
          continue;
        }
        if (ptr[2]=='I' && ptr[3]=='-') {  // -BI-
          if (!poly(83))
            return false;
          ptr += 4;
          continue;
        }
        break;
      case 'C':
        if (ptr[2]=='A' && ptr[3]=='-') {  // -CA-
          if (state == 0 && atend(ptr+4)) {
            prev = atom(20,0);
            prev->SetFormalCharge(2);
            state = 2;
            ptr += 4;
            continue;
          }
          if (!poly(20))
            return false;
          ptr += 4;
          continue;
        }
        if (ptr[2]=='D' && ptr[3]=='-') {  // -CD-
          if (!poly(48))
            return false;
          ptr += 4;
          continue;
        }
        if (ptr[2]=='U' && ptr[3]=='-') {  // -CU-
          if (!poly(29))
            return false;
          ptr += 4;
          continue;
        }
        break;
      case 'E':
        if (ptr[2]=='-') {  // -E-
          if (!poly(35))
            return false;
          ptr += 3;
          continue;
        }
        break;
      case 'F':
        if (ptr[2]=='E' && ptr[3]=='-') {  // -FE-
          if (!poly(26))
            return false;
          ptr += 4;
          continue;
        }
        break;
      case 'G':
        if (ptr[2]=='-') {  // -G-
          if (!poly(17))
            return false;
          ptr += 3;
          continue;
        }
        if (ptr[2]=='A' && ptr[3]=='-') {  // -GA-
          if (!poly(31))
            return false;
          ptr += 4;
          continue;
        }
        if (ptr[2]=='E' && ptr[3]=='-') {  // -GE-
          if (!poly(32))
            return false;
          ptr += 4;
          continue;
        }
        break;
      case 'H':
        if (ptr[2]=='G' && ptr[3]=='-') {  // -HG-
          if (!poly(80))
            return false;
          ptr += 4;
          continue;
        }
        break;
      case 'I':
        if (ptr[2]=='-') {  // -I-
          if (!poly(53))
            return false;
          ptr += 3;
          continue;
        }
        if (ptr[2]=='N' && ptr[3]=='-') {  // -IN-
          if (!poly(49))
            return false;
          ptr += 4;
          continue;
        }
        break;
      case 'K':
        if (ptr[2]=='A' && ptr[3]=='-') {  // -KA-
          if (state == 0 && atend(ptr+4)) {
            prev = atom(19,0);
            prev->SetFormalCharge(1);
            state = 2;
            ptr += 4;
            continue;
          }
        }
        break;
      case 'L':
        if (ptr[2]=='I' && ptr[3]=='-') {  // -LI-
          if (state == 0 && atend(ptr+4)) {
            prev = atom(3,0);
            prev->SetFormalCharge(1);
            state = 2;
            ptr += 4;
            continue;
          }
        }
        break;
      case 'M':
        if (ptr[2]=='G' && ptr[3]=='-') {  // -MG-
          if (state == 0 && atend(ptr+4)) {
            prev = atom(12,0);
            prev->SetFormalCharge(2);
            state = 2;
            ptr += 4;
            continue;
          }
          if (!poly(12))
            return false;
          ptr += 4;
          continue;
        }
        if (ptr[2]=='N' && ptr[3]=='-') {  // -MN-
          if (!poly(25))
            return false;
          ptr += 4;
          continue;
        }
        break;
      case 'N':
        if (ptr[2]=='A' && ptr[3]=='-') {  // -NA-
          if (state == 0 && atend(ptr+4)) {
            prev = atom(11,0);
            prev->SetFormalCharge(1);
            state = 2;
            ptr += 4;
            continue;
          }
          if (!poly(11))
            return false;
          ptr += 4;
          continue;
        }
        if (ptr[2]=='I' && ptr[3]=='-') {  // -NI-
          if (!poly(28))
            return false;
          ptr += 4;
          continue;
        }
        break;
      case 'P':
        if (ptr[2]=='B' && ptr[3]=='-') {  // -PB-
          if (!poly(82))
            return false;
          ptr += 4;
          continue;
        }
        if (ptr[2]=='T' && ptr[3]=='-') {  // -PT-
          if (!poly(78))
            return false;
          ptr += 4;
          continue;
        }
        break;
      case 'S':
        if (ptr[2]=='B' && ptr[3]=='-') {  // -SB-
          if (!poly(51))
            return false;
          ptr += 4;
          continue;
        }
        if (ptr[2]=='E' && ptr[3]=='-') {  // -SE-
          if (!poly(34))
            return false;
          ptr += 4;
          continue;
        }
        if (ptr[2]=='I' && ptr[3]=='-') {  // -SI-
          if (!poly(14))
            return false;
          ptr += 4;
          continue;
        }
        if (ptr[2]=='N' && ptr[3]=='-') {  // -SN-
          if (!poly(50))
            return false;
          ptr += 4;
          continue;
        }
        if (ptr[2]=='R' && ptr[3]=='-') {  // -SR-
          if (!poly(38))
            return false;
          ptr += 4;
          continue;
        }
        break;
      case 'T':
        if (ptr[2]=='L' && ptr[3]=='-') {  // -TL-
          if (!poly(81))
            return false;
          ptr += 4;
          continue;
        }
        break;
      case 'U':
        if (ptr[2]=='R' && ptr[3]=='-') {  // -UR-
          if (!poly(92))
            return false;
          ptr += 4;
          continue;
        }
        break;
      case 'V':
        if (ptr[2]=='A' && ptr[3]=='-') {  // -VA-
          if (!poly(23))
            return false;
          ptr += 4;
          continue;
        }
        break;
      case 'W':
        if (ptr[2]=='O' && ptr[3]=='-') {  // -WO-
          if (!poly(74))
            return false;
          ptr += 4;
          continue;
        }
        break;
      case 'Z':
        if (ptr[2]=='N' && ptr[3]=='-') {  // -ZN-
          if (!poly(30))
            return false;
          ptr += 4;
          continue;
        }
        break;
      case ' ':
        if (ptr[2]>='A' && ptr[2]<='Z' &&
            (ptr[3]=='L' || ptr[3]=='T')) {  // "- [A-Z][LT]" e.g. "- AL"
          if (order==1 || order==2) {
            const char *loc_ptr = ptr+2;
            unsigned int loc = *loc_ptr-'A';
            ptr += 3;
            if (!parse_ring())
              return false;
            if (loc >= rings.back().size()) {
              ptr = loc_ptr;
              return error();
            }
            fuse(prev,rings.back()[loc],order);
            state = 2;
            order = 0;
            continue;
          }
        }
        break;
      }
      break;

    case '&':
      term();
      if (ptr[1]=='&') {
        if (!double_pop())
          return false;
        ptr += 2;
        continue;
      }
      if (!pop())
        return false;
      ptr++;
      continue;

    case ' ':
      if (ptr[1]=='&') {
        if (state != 0) {
          drain();
          pending = PENDING_NONE;
          prev = nullptr;
          rings.clear();
          state = 0;
          order = 0;
          ptr += 2;
          continue;
        }
      }
      if (ptr[1]>='A' && ptr[1]<='Z') {
        if (rings.empty())
          return error();
        term();
        unsigned int loc = ptr[1]-'A';
        if (loc < rings.back().size()) {
          prev = rings.back()[loc];
          pending = PENDING_METHYL;
          order = 1;
          state = 1;
          ptr += 2;
          continue;
        }
      }
      break;

    default:
      if (*ptr>='1' && *ptr<='9') {
        if (!alkane())
          return false;
        continue;
      }
      break;
    }
    return error();
    }
  }

  // ptr is E, F, G or I
  int parse_inorganic_halide(unsigned int cation, unsigned int count,
                             unsigned int anion) {
    if (count != 1)
      return 0;
    if (ptr[1]>='2' && ptr[1]<='9' && !ptr[2]) {
      count = ptr[1]-'0';
    } else if (!ptr[1]) {
      count = 1;
    } else return 0;

    prev = atom(cation,0);
    for (unsigned int i=0; i<count; i++) {
      OpenBabel::OBAtom* temp = atom(anion,0);
      NMOBMolNewBond(mol,prev,temp,1,false);
    }
    return 1;
  }

  // ptr is O or S, or ptr-1 is SE or TE
  int parse_inorganic_oxide(unsigned int cation, unsigned int count,
                            unsigned int anion) {
    if (count == 1) {
      if (ptr[1]>='2' && ptr[1]<='9' && !ptr[2]) {
        count = ptr[1]-'0';
      } else if (!ptr[1]) {
        count = 1;
      } else return 0;
      prev = atom(cation,0);
      for (unsigned int i=0; i<count; i++) {
        OpenBabel::OBAtom* temp = atom(anion,0);
        NMOBMolNewBond(mol,prev,temp,2,false);
      }
      return 1;
    } else if (count == 2) {
      if (!ptr[1]) {
        prev = atom(anion,0);
        NMOBMolNewBond(mol,prev,atom(cation,0),1,false);
        NMOBMolNewBond(mol,prev,atom(cation,0),1,false);
        return 1;
      } else if (ptr[1]=='3' && !ptr[2]) {
        prev = atom(anion,0);
        OpenBabel::OBAtom* temp = atom(cation,0);
        NMOBMolNewBond(mol,temp,atom(anion,0),2,false);
        NMOBMolNewBond(mol,prev,temp,1,false);
        temp = atom(cation,0);
        NMOBMolNewBond(mol,temp,atom(anion,0),2,false);
        NMOBMolNewBond(mol,prev,temp,1,false);
        return 1;
      } else if (ptr[1]=='5' && !ptr[2]) {
        prev = atom(anion,0);
        OpenBabel::OBAtom* temp = atom(cation,0);
        NMOBMolNewBond(mol,temp,atom(anion,0),2,false);
        NMOBMolNewBond(mol,temp,atom(anion,0),2,false);
        NMOBMolNewBond(mol,prev,temp,1,false);
        temp = atom(cation,0);
        NMOBMolNewBond(mol,temp,atom(anion,0),2,false);
        NMOBMolNewBond(mol,temp,atom(anion,0),2,false);
        NMOBMolNewBond(mol,prev,temp,1,false);
        return 1;
      }
    }
    return 0;
  }

#define WLN_BORATE     1
#define WLN_CARBONATE  2
#define WLN_CARBONYL   3
#define WLN_CYANIDE    4
#define WLN_NITRATE    5
#define WLN_NITRITE    6
#define WLN_SULFATE    7
#define WLN_SULFITE    8

  // ptr is after anion, i.e. multiplier or done.
  int parse_inorganic_salt(unsigned int cation, unsigned int ccount,
                           unsigned int anion, unsigned int acharge)
  {
    unsigned int acount;
    if (ptr[0]=='*' && ptr[1]>='2' && ptr[1]<='9' && !ptr[2]) {
      acount = ptr[1]-'0';
    } else if (!ptr[0]) {
      acount = 1;
    } else return 0;

    if (ccount != acount*acharge) {
      unsigned int ccharge = (acount*acharge)/ccount;
      if (ccount*ccharge != acount*acharge)
        return 0;
      for (unsigned int i=0; i<ccount; i++) {
        prev = atom(cation,0);
        prev->SetFormalCharge(ccharge);
      }
      cation = 0;
    }
    
    OpenBabel::OBAtom* temp;

    for (unsigned int i=0; i<acount; i++) {
      switch (anion) {
      case WLN_BORATE:
        prev = atom(5,0);
        temp = atom(8,0);
        NMOBMolNewBond(mol,prev,temp,1,false);
        if (cation)
          NMOBMolNewBond(mol,temp,atom(cation,0),1,false);
        else temp->SetFormalCharge(-1);
        temp = atom(8,0);
        NMOBMolNewBond(mol,prev,temp,1,false);
        if (cation)
          NMOBMolNewBond(mol,temp,atom(cation,0),1,false);
        else temp->SetFormalCharge(-1);
        temp = atom(8,0);
        NMOBMolNewBond(mol,prev,temp,1,false);
        if (cation)
          NMOBMolNewBond(mol,temp,atom(cation,0),1,false);
        else temp->SetFormalCharge(-1);
        break;
      case WLN_CARBONATE:
        prev = atom(6,0);
        NMOBMolNewBond(mol,prev,atom(8,0),2,false);
        temp = atom(8,0);
        NMOBMolNewBond(mol,prev,temp,1,false);
        if (cation)
          NMOBMolNewBond(mol,temp,atom(cation,0),1,false);
        else temp->SetFormalCharge(-1);
        temp = atom(8,0);
        NMOBMolNewBond(mol,prev,temp,1,false);
        if (cation)
          NMOBMolNewBond(mol,temp,atom(cation,0),1,false);
        else temp->SetFormalCharge(-1);
        break;
      case WLN_SULFATE:
        prev = atom(16,0);
        NMOBMolNewBond(mol,prev,atom(8,0),2,false);
        NMOBMolNewBond(mol,prev,atom(8,0),2,false);
        temp = atom(8,0);
        NMOBMolNewBond(mol,prev,temp,1,false);
        if (cation)
          NMOBMolNewBond(mol,temp,atom(cation,0),1,false);
        else temp->SetFormalCharge(-1);
        temp = atom(8,0);
        NMOBMolNewBond(mol,prev,temp,1,false);
        if (cation)
          NMOBMolNewBond(mol,temp,atom(cation,0),1,false);
        else temp->SetFormalCharge(-1);
        break;
      case WLN_SULFITE:
        prev = atom(16,0);
        NMOBMolNewBond(mol,prev,atom(8,0),2,false);
        temp = atom(8,0);
        NMOBMolNewBond(mol,prev,temp,1,false);
        if (cation)
          NMOBMolNewBond(mol,temp,atom(cation,0),1,false);
        else temp->SetFormalCharge(-1);
        temp = atom(8,0);
        NMOBMolNewBond(mol,prev,temp,1,false);
        if (cation)
          NMOBMolNewBond(mol,temp,atom(cation,0),1,false);
        else temp->SetFormalCharge(-1);
        break;
      }
    }
    return 1;
  }

  // ptr is after anion, i.e. multiplier or done.
  int parse_inorganic_salt1(unsigned int cation, unsigned int ccount,
                            unsigned int anion)
  {
    unsigned int acount;
    if (ptr[0]=='*' && ptr[1]>='2' && ptr[1]<='9' && !ptr[2]) {
      acount = ptr[1]-'0';
    } else if (!ptr[0]) {
      acount = 1;
    } else return 0;

    if (ccount == 1) {
      prev = atom(cation,0);
    } else return 0;

    OpenBabel::OBAtom* temp;
    OpenBabel::OBAtom* temp2;

    for (unsigned int i=0; i<acount; i++) {
      switch (anion) {
      case WLN_CARBONYL:
        temp = atom(6,0);
        NMOBMolNewBond(mol,temp,atom(8,0),2,false);
        break;
      case WLN_CYANIDE:
        temp = atom(6,0);
        NMOBMolNewBond(mol,temp,atom(7,0),3,false);
        break;
      case WLN_NITRATE:
        temp2 = atom(7,0);
        temp2->SetFormalCharge(1);
        NMOBMolNewBond(mol,atom(8,0),temp2,2,false);
        temp = atom(8,0);
        temp->SetFormalCharge(-1);
        NMOBMolNewBond(mol,temp2,temp,1,false);
        temp = atom(8,0);
        NMOBMolNewBond(mol,temp2,temp,1,false);
        break;
      case WLN_NITRITE:
        temp = atom(8,0);
        temp2 = atom(7,0);
        NMOBMolNewBond(mol,temp,temp2,2,false);
        temp = atom(8,0);
        NMOBMolNewBond(mol,temp,temp2,1,false);
        break;
      default:
        /* Internal error */
        return 0;
      }
      if (ccount == 1)
        NMOBMolNewBond(mol,prev,temp,1,false);
      else temp->SetFormalCharge(-1);
    }
    return 1;
  }

  // 1 success, -1 failure, 0 unknown
  int parse_inorganic() {
    ptr = orig;
    if (ptr[0]=='W' && ptr[1]=='L' && ptr[2]=='N' &&
        ptr[3]==':' && ptr[4]==' ')
      ptr += 5;

    unsigned int cation = 0;
    switch (*ptr) {
    case 'A':
      if (ptr[1]=='L')  // AL
        cation = 13;
      else if (ptr[1]=='G')  // AG
        cation = 47;
      else if (ptr[1]=='U')  // AU
        cation = 79;
      break;
    case 'B':
      if (ptr[1]=='A')  // BA
        cation = 56;
      else if (ptr[1]=='E')  // BE
        cation = 4;
      break;
    case 'C':
      if (ptr[1]=='A')  // CA
        cation = 20;
      else if (ptr[1]=='D')  // CD
        cation = 48;
      else if (ptr[1]=='E')  // CE
        cation = 58;
      else if (ptr[1]=='N')  // CN
        cation = 112;
      else if (ptr[1]=='O')  // CO
        cation = 27;
      else if (ptr[1]=='R')  // CR
        cation = 24;
      else if (ptr[1]=='S')  // CS
        cation = 55;
      else if (ptr[1]=='U')  // CU
        cation = 29;
      break;
    case 'D':
      if (ptr[1]=='Y')  // DY
        cation = 66;
      break;
    case 'E':
      if (ptr[1]=='R')  // ER
        cation = 68;
      else if (ptr[1]=='S')  // ES
        cation = 99;
      else if (ptr[1]=='U')  // EU
        cation = 63;
      break;
    case 'F':
      if (ptr[1]=='E')  // FE
        cation = 26;
      break;
    case 'G':
      if (ptr[1]=='A')  // GA
        cation = 31;
      else if (ptr[1]=='D')  // GD
        cation = 64;
      else if (ptr[1]=='E')  // GE
        cation = 32;
      break;
    case 'H':
      if (ptr[1]=='G')  // HG
        cation = 80;
      else if (ptr[1]=='O')  // HO
        cation = 67;
      break;
    case 'I':
      if (ptr[1]=='N')  // IN
        cation = 49;
      break;
    case 'K':
      if (ptr[1]=='A')  // KA
        cation = 19;
      break;
    case 'L':
      if (ptr[1]=='A')  // LA
        cation = 57;
      else if (ptr[1]=='I')  // LI
        cation = 3;
      else if (ptr[1]=='U')  // LU
        cation = 71;
      break;
    case 'M':
      if (ptr[1]=='G')  // MG
        cation = 12;
      else if (ptr[1]=='N')  // MN
        cation = 25;
      else if (ptr[1]=='O')  // MO
        cation = 42;
      break;
    case 'N':
      if (ptr[1]=='A')  // NA
        cation = 11;
      else if (ptr[1]=='D')  // ND
        cation = 60;
      else if (ptr[1]=='I')  // NI
        cation = 28;
      break;
    case 'P':
      if (ptr[1]=='A')  // PA
        cation = 91;
      else if (ptr[1]=='B')  // PB
        cation = 82;
      else if (ptr[1]=='D')  // PD
        cation = 46;
      else if (ptr[1]=='M')  // PM
        cation = 61;
      else if (ptr[1]=='O')  // PO
        cation = 84;
      else if (ptr[1]=='R')  // PR
        cation = 59;
      else if (ptr[1]=='T')  // PT
        cation = 78;
      else if (ptr[1]=='U')  // PU
        cation = 94;
      break;
    case 'R':
      if (ptr[1]=='A')  // RA
        cation = 88;
      else if (ptr[1]=='B')  // RB
        cation = 37;
      else if (ptr[1]=='E')  // RE
        cation = 75;
      else if (ptr[1]=='F')  // RF
        cation = 104;
      else if (ptr[1]=='H')  // RH
        cation = 45;
      else if (ptr[1]=='N')  // RN
        cation = 86;
      else if (ptr[1]=='U')  // RU
        cation = 44;
      break;
    case 'S':
      if (ptr[1]=='B')  // SB
        cation = 51;
      else if (ptr[1]=='C')  // SC
        cation = 21;
      else if (ptr[1]=='E')  // SE
        cation = 34;
      else if (ptr[1]=='G')  // SG
        cation = 106;
      else if (ptr[1]=='I')  // SI
        cation = 14;
      else if (ptr[1]=='M')  // SM
        cation = 62;
      else if (ptr[1]=='N')  // SN
        cation = 50;
      else if (ptr[1]=='R')  // SR
        cation = 38;
      break;
    case 'T':
      if (ptr[1]=='A')  // TA
        cation = 73;
      else if (ptr[1]=='B')  // TB
        cation = 65;
      else if (ptr[1]=='C')  // TC
        cation = 43;
      else if (ptr[1]=='H')  // TH
        cation = 90;
      else if (ptr[1]=='I')  // TI
        cation = 22;
      else if (ptr[1]=='L')  // TL
        cation = 81;
      else if (ptr[1]=='M')  // TM
        cation = 69;
      break;
    case 'U':
      if (ptr[1]=='R')  // UR
        cation = 92;
      break;
    case 'V':
      if (ptr[1]=='A')  // VA
        cation = 23;
      break;
    case 'W':
      if (ptr[1]=='O')  // WO
        cation = 74;
      break;
    case 'X':
      if (ptr[1]=='E')  // XE
        cation = 54;
      break;
    case 'Y':
      if (ptr[1]=='T')  // YT
        cation = 39;
      break;
    case 'Z':
      if (ptr[1]=='N')  // ZN
        cation = 30;
      else if (ptr[1]=='R')  // ZR
        cation = 40;
      break;
    }

    if (!cation)
      return 0;
    unsigned int count;
    if (ptr[2]>='2' && ptr[2]<='9' && ptr[3]==' ') {
      count = ptr[2]-'0';
      ptr += 4;
    } else if (ptr[2]==' ') {
      count = 1;
      ptr += 3;
    } else return 0;

    switch (*ptr) {
    case 'B':
      if (ptr[1]=='-' && ptr[2]=='O' && ptr[3]=='3') {  // B-O3
        ptr += 4;
        return parse_inorganic_salt(cation,count,WLN_BORATE,3);
      }
      break;
    case 'C':
      if (ptr[1]=='-' && ptr[2]=='N') {  // C-N
        ptr += 3;
        return parse_inorganic_salt1(cation,count,WLN_CYANIDE);
      }
      if (ptr[1]=='-' && ptr[2]=='O' && ptr[3]=='3') {  // C-O3
        ptr += 4;
        return parse_inorganic_salt(cation,count,WLN_CARBONATE,2);
      }
      if (ptr[1]=='-' && ptr[2]=='O') {  // C-O
        ptr += 3;
        return parse_inorganic_salt1(cation,count,WLN_CARBONYL);
      }
      if (ptr[1]=='N') {
        ptr += 2;
        return parse_inorganic_salt1(cation,count,WLN_CYANIDE);
      }
      break;
    case 'E':
      return parse_inorganic_halide(cation,count,35);
    case 'F':
      return parse_inorganic_halide(cation,count,9);
    case 'G':
      return parse_inorganic_halide(cation,count,17);
    case 'I':
      return parse_inorganic_halide(cation,count,53);
    case 'N':
      if (ptr[1]=='-' && ptr[2]=='O' && ptr[3]=='2') {
        ptr += 4;
        return parse_inorganic_salt1(cation,count,WLN_NITRITE);
      }
      if (ptr[1]=='-' && ptr[2]=='O' && ptr[3]=='3') {
        ptr += 4;
        return parse_inorganic_salt1(cation,count,WLN_NITRATE);
      }
      break;
    case 'O':
      return parse_inorganic_oxide(cation,count,8);
    case 'S':
      if (ptr[1]=='-' && ptr[2]=='O' && ptr[3]=='3') {  // S-O3
        ptr += 4;
        return parse_inorganic_salt(cation,count,WLN_SULFITE,2);
      }
      if (ptr[1]=='-' && ptr[2]=='O' && ptr[3]=='4') {  // S-O4
        ptr += 4;
        return parse_inorganic_salt(cation,count,WLN_SULFATE,2);
      }
      if (ptr[1]=='E') {
        ptr++;
        return parse_inorganic_oxide(cation,count,34);
      }
      return parse_inorganic_oxide(cation,count,16);
    case 'T':
      if (ptr[1]=='E') {
        ptr++;
        return parse_inorganic_oxide(cation,count,52);
      }
      break;
    }
    return 0;
  }
};


bool NMReadWLN(const char *ptr, OpenBabel::OBMol* mol)
{
  WLNParser wp(ptr,mol);
  int result = wp.parse_inorganic();
  if (result == 0) {
    if (!wp.parse())
      return false;
  } else if (result < 0)
    return false;

  mol->SetDimension(0);
  return NMOBSanitizeMol(mol);
}

