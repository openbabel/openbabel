/**********************************************************************
depict.cpp - 2D Depiction of molecules using OBPainter.

Copyright (C) 2009-2010 by Tim Vandermeersch
Some portions Copyright (C) 2009 by Chris Morley

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
#include <openbabel/bond.h>
#include <openbabel/ring.h>
#include <openbabel/alias.h>
#include <openbabel/generic.h>
#include <openbabel/depict/depict.h>
#include <openbabel/depict/painter.h>
#include <openbabel/elements.h>
#include <algorithm> // std::reverse
#include <iterator> // std::istream_iterator
#include <openbabel/stereo/stereo.h>
#include <openbabel/stereo/cistrans.h>
#include <openbabel/obiter.h>
#include <openbabel/obfunctions.h>

#include <cmath>

#include <iostream>
using namespace std;

namespace OpenBabel
{

  enum {
    Left,
    Right,
    Up,
    Down
  };

  // defined in mol.cpp
  extern bool SortAtomZ(const pair<OBAtom*,double> &a, const pair<OBAtom*,double> &b);

  class OBDepictPrivate
  {
    public:
      OBDepictPrivate() : mol(0), painter(0), bondLength(40.0), penWidth(2.0),
          bondSpacing(6.0), bondWidth(8.0), fontSize(16), subscriptSize(13),
          aliasMode(false), bondColor("black"), options(0){}
      virtual ~OBDepictPrivate(){};

      virtual void DrawSimpleBond(OBAtom *beginAtom, OBAtom *endAtom, int order, bool crossed_bond=false);
      virtual void DrawWedge(OBAtom *beginAtom, OBAtom *endAtom);
      virtual void DrawHash(OBAtom *beginAtom, OBAtom *endAtom);
      virtual void DrawWobblyBond(OBAtom *beginAtom, OBAtom *endAtom);
      virtual void DrawRingBond(OBAtom *beginAtom, OBAtom *endAtom, const vector3 &center, int order);
      virtual void DrawAtom(OBAtom* atom);
      virtual void DrawAtomLabel(const std::string &label, int alignment, const vector3 &pos);

      virtual void DrawRing(OBRing *ring, OBBitVec &drawnBonds);
      virtual void DrawAromaticRing(OBRing *ring, OBBitVec &drawnBonds);

      bool HasLabel(OBAtom *atom);
      void SetWedgeAndHash(OBMol* mol);

      OBMol     *mol;
      OBPainter *painter;
      double     bondLength;
      double     penWidth;
      double     bondSpacing;
      double     bondWidth;
    // for z-scaling the opacity
      double     zScale;
      double     zMin;
      //bool       drawTerminalC;
      int        fontSize, subscriptSize;
      bool       aliasMode;
      std::string fontFamily;
      OBColor    bondColor;
      unsigned   options;
  };

  class OBDepictPrivateBallAndStick : public OBDepictPrivate
  {
  public:
        OBDepictPrivateBallAndStick(bool symbolOnBall) : m_symbolOnBall(symbolOnBall){}

        bool m_symbolOnBall;

        void DrawSimpleBond(OBAtom *beginAtom, OBAtom *endAtom, int order, bool crossed_bond=false);
        void DrawWedge(OBAtom *beginAtom, OBAtom *endAtom);
        void DrawHash(OBAtom *beginAtom, OBAtom *endAtom);
        void DrawWobblyBond(OBAtom *beginAtom, OBAtom *endAtom);
        void DrawRingBond(OBAtom *beginAtom, OBAtom *endAtom, const vector3 &center, int order);
        void DrawAtom(OBAtom *atom);
        void DrawAtomLabel(const std::string &label, int alignment, const vector3 &pos);

        void DrawRing(OBRing *ring, OBBitVec &drawnBonds);
        void DrawAromaticRing(OBRing *ring, OBBitVec &drawnBonds);

  private:
        void DrawAromaticRingBond(OBAtom *prevAtom, OBAtom *beginAtom, OBAtom *endAtom, OBAtom *nextAtom, const vector3 &center, double dist);
        double GetAtomRadius(OBAtom* atom);

  };

  OBDepict::OBDepict(OBPainter *painter) : d(new OBDepictPrivate)
  {
    d->painter = painter;
  }

  OBDepict::OBDepict(OBPainter *painter, bool withBall, bool symbolOnBall)
  {
    if(withBall)
      d = new OBDepictPrivateBallAndStick(symbolOnBall);
    else
      d = new OBDepictPrivate();
    d->painter = painter;
  }

  OBDepict::~OBDepict()
  {
    delete d->mol;
    d->mol = NULL;
    delete d;
  }

  void OBDepict::SetBondLength(double length)
  {
    d->bondLength = length;
  }

  double OBDepict::GetBondLength() const
  {
    return d->bondLength;
  }

  void OBDepict::SetPenWidth(double width)
  {
    d->penWidth = width;
    d->painter->SetPenWidth(width);
  }

  double OBDepict::GetPenWidth() const
  {
    return d->penWidth;
  }

  void OBDepict::SetBondSpacing(double spacing)
  {
    d->bondSpacing = spacing;
  }

  double OBDepict::GetBondSpacing() const
  {
    return d->bondSpacing;
  }

  void OBDepict::SetBondWidth(double width)
  {
    d->bondWidth = width;
  }

  double OBDepict::GetBondWidth() const
  {
    return d->bondWidth;
  }

/*  void OBDepict::SetDrawingTerminalCarbon(bool enabled)
  {
    d->drawTerminalC = enabled;
  }

  bool OBDepict::GetDrawingTerminalCarbon() const
  {
    return d->drawTerminalC;
  }
*/
  void OBDepict::SetOption(unsigned opts)
  {
    d->options |= opts;
  }

  unsigned OBDepict::GetOptions() const
  {
    return d->options;
  }
  void OBDepict::ClearOptions()
  {
    d->options = 0;
  }

  void OBDepict::SetFontFamily(const std::string &family)
  {
    d->fontFamily = family;
    d->painter->SetFontFamily(family);
  }

  const std::string& OBDepict::GetFontFamily() const
  {
    return d->fontFamily;
  }

  void OBDepict::SetFontSize(int pointSize, bool subscript)
  {
    if (subscript) {
      d->subscriptSize = pointSize;
      return;
    }

    d->fontSize = pointSize;
    d->subscriptSize = (int)(0.85 * pointSize);
  }

  int OBDepict::GetFontSize(bool subscript) const
  {
    if (subscript)
      return d->subscriptSize;
    return d->fontSize;
  }

  void OBDepict::SetAliasMode(bool b)
  {
    d->aliasMode = b;
  }

  //Color is not quite properly integrated into OBDepict, but is needed if
  //element-dependent coloring is to be used.
  void OBDepict::SetBondColor(const std::string& scolor)
  {
    d->bondColor = scolor;
  }

  int GetLabelAlignment(OBAtom *atom)
  {
    // compute the sum of the bond vectors, this gives
    vector3 direction(VZero);
    OBBondIterator i;
    for (OBAtom *nbr = atom->BeginNbrAtom(i); nbr; nbr = atom->NextNbrAtom(i))
      direction += atom->GetVector() - nbr->GetVector();

    const double bias = -0.1; //towards left-alignment, which is more natural
    int alignment = 0;
    if ((atom->GetExplicitDegree() == 2) && (abs(direction.y()) > abs(direction.x()))) {
      if (direction.y() <= 0.0)
        alignment = Up;
      else
        alignment = Down;
    } else {
      if (direction.x() < bias)
        alignment = Right;
      else
        alignment = Left;
    }

    return alignment;
  }

  unsigned int GetAtomSymClass(OBAtom *atom)
  {
    OBPairData *pd = dynamic_cast<OBPairData*>(atom->GetParent()->GetData("OpenBabel Symmetry Classes"));
    if (pd) {

      cout << "same? = " << pd->GetValue() << endl;

      istringstream iss(pd->GetValue());
      std::vector<unsigned int> symmetry_classes;
      copy(istream_iterator<unsigned int>(iss),
           istream_iterator<unsigned int>(),
           back_inserter<vector<unsigned int> >(symmetry_classes));
      // Now find the number of unique elements
      vector<unsigned int> copy_sym = symmetry_classes;
      sort(copy_sym.begin(), copy_sym.end());
      vector<unsigned int>::iterator end_pos = unique(copy_sym.begin(), copy_sym.end()); // Requires sorted elements
      int nclasses = end_pos - copy_sym.begin();

      cout << "sym_class[" << atom->GetIndex() << "] = " << symmetry_classes.at(atom->GetIndex()) << endl;
      return symmetry_classes.at(atom->GetIndex());
    }

    return 99;
  }

  bool OBDepict::AddAtomLabels(AtomLabelType type)
  {
    d->painter->SetPenColor(OBColor("red"));
    d->painter->SetFillColor(OBColor("red"));
    d->painter->SetFontSize((int)(GetFontSize() * 0.8));// smaller text
    OBAtomIterator i;
    for (OBAtom *atom = d->mol->BeginAtom(i); atom; atom = d->mol->NextAtom(i)) {
      vector3 pos(atom->GetVector());
      std::stringstream ss;
      switch (type) {
        case AtomId:
          ss << atom->GetId();
          d->painter->DrawText(pos.x(), pos.y(), ss.str());
          break;
        case AtomSymmetryClass:
          ss << GetAtomSymClass(atom);
          d->painter->DrawText(pos.x(), pos.y(), ss.str());
          break;
        case AtomIndex:
          ss << atom->GetIdx();
          d->painter->DrawText(pos.x(), pos.y(), ss.str());
          break;

        default:
          break;
      }
    }

    return true;
  }

  void OBDepictPrivate::DrawRing(OBRing *ring, OBBitVec &drawnBonds)
  {
    std::vector<int> indexes = ring->_path;
    vector3 center(VZero);
    for (std::vector<int>::iterator l = indexes.begin(); l != indexes.end(); ++l) {
      center += mol->GetAtom(*l)->GetVector();
    }
    center /= indexes.size();

    for (unsigned int l = 0; l < indexes.size(); ++l) {
      OBAtom *begin = mol->GetAtom(indexes[l]);
      OBAtom *end;
      if (l+1 < indexes.size())
        end = mol->GetAtom(indexes[l+1]);
      else
        end = mol->GetAtom(indexes[0]);

      OBBond *ringBond = mol->GetBond(begin, end);
      if (drawnBonds.BitIsSet(ringBond->GetId()))
        continue;

      if((options & OBDepict::internalColor) && ringBond->HasData("color"))
        painter->SetPenColor(OBColor(ringBond->GetData("color")->GetValue()));
      else
        painter->SetPenColor(bondColor);

      DrawRingBond(begin, end, center, ringBond->GetBondOrder());
      drawnBonds.SetBitOn(ringBond->GetId());
    }
  }

  void OBDepictPrivate::DrawAromaticRing(OBRing *ring, OBBitVec &drawnBonds)
  {
    DrawRing(ring, drawnBonds);
  }

  enum Radical { NOT_RADICAL, ONE_DOT, TWO_DOT };

  // Assign 0, 1, or 2 radical dots
  // - if spin is specified, then this determines the number of dots
  // - otherwise, the degree of undervalence determines it
  static Radical AssignRadicalDots(OBAtom* atom)
  {
    unsigned int spin = atom->GetSpinMultiplicity();
    if (spin)
      return spin == 2 ? TWO_DOT : ONE_DOT;

    unsigned int actualvalence = atom->GetTotalValence();
    unsigned int typicalvalence = GetTypicalValence(atom->GetAtomicNum(), actualvalence, atom->GetFormalCharge());
    int diff = typicalvalence - actualvalence;
    if (diff <= 0)
      return NOT_RADICAL;
    return diff == 2 ? TWO_DOT : ONE_DOT;
  }

  bool OBDepict::DrawMolecule(OBMol *mol)
  {
    if (!d->painter)
      return false;

    if (d->mol != NULL)
      delete d->mol;
    d->mol = new OBMol(*mol); // Copy it

    double width=0.0, height=0.0;

    OBAtom *atom;
    OBBondIterator j;
    OBAtomIterator i;


    // Determine which should be wedge and hash bonds...
    // Note: we need to do this before we invert the y-coordinate for depiction
    std::map<OBBond*, enum OBStereo::BondDirection> updown;
    std::map<OBBond*, OBStereo::Ref> from;
    TetStereoToWedgeHash(*d->mol, updown, from);

    if(mol->NumAtoms()>0) {
      // scale bond lengths and invert the y coordinate (both SVG and Cairo use top left as the origin)
      double bondLengthSum = 0.0;
      for (OBBond *bond = mol->BeginBond(j); bond; bond = mol->NextBond(j))
        bondLengthSum += bond->GetLength();
      const double averageBondLength = bondLengthSum / mol->NumBonds();
      double f;
      if(mol->NumBonds()>0) {
        f = d->bondLength / averageBondLength;
      } else if(mol->NumAtoms()>1){
        /* In molecules without bonds but more than one atom,
           use the minimum atom distance as a substitute for average bond length. */
        f = DBL_MAX;
        OBAtomIterator i2;
        OBAtom *atom2;
        vector3 a1pos;
        double currdist;
        for (atom = d->mol->BeginAtom(i); atom; atom = d->mol->NextAtom(i)) {
          a1pos = atom->GetVector();
          for (atom2 = d->mol->BeginAtom(i2); atom2; atom2 = d->mol->NextAtom(i2)) {
            if(atom != atom2) {
              currdist = a1pos.distSq(atom2->GetVector());
              if(currdist < f && currdist != 0)
                f = currdist;
            }
          }
        }
        f = d->bondLength / sqrt(f);
      } else
        f = 1.0;
      for (atom = d->mol->BeginAtom(i); atom; atom = d->mol->NextAtom(i))
        atom->SetVector(atom->GetX() * f, - atom->GetY() * f, atom->GetZ());

      // find min/max values
      double min_x, max_x;
      double min_y, max_y;
      double min_z, max_z;
      atom = d->mol->BeginAtom(i);
      if (atom != NULL) {
        min_x = max_x = atom->GetX();
        min_y = max_y = atom->GetY();
        min_z = max_z = atom->GetZ();
        for (atom = d->mol->NextAtom(i); atom; atom = d->mol->NextAtom(i)) {
          min_x = std::min(min_x, atom->GetX());
          max_x = std::max(max_x, atom->GetX());
          min_y = std::min(min_y, atom->GetY());
          max_y = std::max(max_y, atom->GetY());
          min_z = std::min(min_z, atom->GetZ());
          max_z = std::max(max_z, atom->GetZ());
        }
      }

      double margin;
      if (d->options & noMargin)
        margin = 5.0;
      else
        margin = 40.0;
      // translate all atoms so the bottom-left atom is at margin,margin
      for (atom = d->mol->BeginAtom(i); atom; atom = d->mol->NextAtom(i))
        atom->SetVector(atom->GetX() - min_x + margin, atom->GetY() - min_y + margin, atom->GetZ());

      width  = max_x - min_x + 2*margin;
      height = max_y - min_y + 2*margin;

      d->zScale = max_z - min_z;
      if (fabs(d->zScale) < 1.0e-1)
        d->zScale = 0.0;
      d->zMin = min_z;

      //d->painter->SetPenWidth(d->penWidth);
      //d->painter->SetPenColor(d->pen));
      //d->painter->SetFillColor(OBColor("black"));
    }

    d->painter->NewCanvas(width, height);

    // Identify and remember the ring bonds according to the SSSR
    // - note that OBBond->IsInRing() includes bonds not included in the SSSR as the SSSR excludes very large rings
    std::vector<OBRing*> rings(mol->GetSSSR());
    OBBitVec ringBonds;
    for (std::vector<OBRing*>::iterator k = rings.begin(); k != rings.end(); ++k) {
      OBRing *ring = *k;
      std::vector<int> indexes = ring->_path;
      for (unsigned int l = 0; l < indexes.size(); ++l) {
        OBAtom *begin = d->mol->GetAtom(indexes[l]);
        OBAtom *end;
        if (l+1 < indexes.size())
          end = d->mol->GetAtom(indexes[l+1]);
        else
          end = d->mol->GetAtom(indexes[0]);

        OBBond *ringBond = d->mol->GetBond(begin, end);
        ringBonds.SetBitOn(ringBond->GetId());
      }
    }

    // draw bonds
    for (OBBond *bond = d->mol->BeginBond(j); bond; bond = d->mol->NextBond(j)) {
      OBAtom *begin = bond->GetBeginAtom();
      OBAtom *end = bond->GetEndAtom();
      if((d->options & internalColor) && bond->HasData("color"))
        d->painter->SetPenColor(OBColor(bond->GetData("color")->GetValue()));
      else
        d->painter->SetPenColor(d->bondColor);

      if(from.find(bond)!=from.end()) {
        //is a wedge or hash bond
        if(from[bond]==bond->GetEndAtom()->GetId())
          swap(begin, end);
        if(updown[bond]==OBStereo::UpBond)
          d->DrawWedge(begin, end);
        else if(updown[bond]==OBStereo::DownBond)
          d->DrawHash(begin, end);
        else {
          //This is a bond to a chiral center specified as unknown
          d->DrawWobblyBond(begin, end);
        }
      }
      else if (!ringBonds.BitIsSet(bond->GetId())) { // Ring bonds are handled below
        bool crossed_dbl_bond = false;
        OBStereoFacade sf(d->mol);
        if (sf.HasCisTransStereo(bond->GetId())) {
          OBCisTransStereo *ct = sf.GetCisTransStereo(bond->GetId());
          if (!ct->GetConfig().specified)
            crossed_dbl_bond = true;
        }
        d->DrawSimpleBond(begin, end, bond->GetBondOrder(), crossed_dbl_bond);
      }
    }

    // draw ring bonds
    OBBitVec drawnBonds;
    // draw aromatic rings first, looks better since all double bonds will
    // be inside aromatic rings
    for (std::vector<OBRing*>::iterator k = rings.begin(); k != rings.end(); ++k) {
      OBRing *ring = *k;
      if (ring->IsAromatic())
        d->DrawAromaticRing(ring, drawnBonds);
    }
    // draw aliphatic rings
    for (std::vector<OBRing*>::iterator k = rings.begin(); k != rings.end(); ++k) {
      OBRing *ring = *k;
      if (!ring->IsAromatic())
        d->DrawRing(ring, drawnBonds);
    }

    vector<pair<OBAtom*,double> > zsortedAtoms;
    vector<int> zsorted;
    unsigned int a;
    for (a = 0, atom = d->mol->BeginAtom(i) ; atom ; atom = d->mol->NextAtom(i), ++a)
      {
        pair<OBAtom*,double> entry(atom, atom->GetVector().z());
        zsortedAtoms.push_back(entry);
      }
    sort(zsortedAtoms.begin(), zsortedAtoms.end(), SortAtomZ);
    unsigned int max = zsortedAtoms.size();
    for (a = 0 ; a < max ; a++ ) {
      atom   = zsortedAtoms[a].first;
      double x = atom->GetX();
      double y = atom->GetY();

      d->DrawAtom(atom);

      // draw atom labels
      int alignment = GetLabelAlignment(atom);
      bool rightAligned = false;
      switch (alignment) {
        case Right:
          rightAligned = true;
          /* no break */
        default:
          break;
      }

      if((d->options & internalColor) && atom->HasData("color"))
        d->painter->SetPenColor(OBColor(atom->GetData("color")->GetValue()));
      else if(d->options & bwAtoms)
        d->painter->SetPenColor(d->bondColor);
      else {
        double r, g, b;
        OBElements::GetRGB(atom->GetAtomicNum(), &r, &g, &b);
        d->painter->SetPenColor(OBColor(r, g, b));
      }

      //charge and radical
      int charge = atom->GetFormalCharge();
      Radical radical = AssignRadicalDots(atom); // none, one or two
      if(charge || radical != NOT_RADICAL) {
        OBFontMetrics metrics = d->painter->GetFontMetrics("N");
        double yoffset = d->HasLabel(atom) ? -0.2 * metrics.height : -0.2 * metrics.height;
        /*switch (GetLabelAlignment(atom)) {
          case Up:
          case Left:
          case Right:
            yoffset = - 1.2 * metrics.height;
        }*/
        stringstream ss;
        if(charge) {
          if(abs(charge)!=1)
            ss << abs(charge);
          if(charge>0)
            ss << '+';
          else if (charge<-1) //use underscore for single negative charge and minus if multiple
            ss << '-';
          else
          {
            ss << '_';
            yoffset -= 0.5 * metrics.height;
          }
        }
        d->painter->DrawText(x + 0.4*metrics.width, y+yoffset, ss.str());
        if (radical != NOT_RADICAL) {
          string radchars;
          radchars = radical == ONE_DOT ? "." : "..";
          d->painter->SetFontSize(2 * metrics.fontSize);
          d->painter->DrawText(x + (0.4 + ss.str().size())*metrics.width,
            y + yoffset, radchars);
        }
        d->painter->SetFontSize(metrics.fontSize);//restore
      }

      if (atom->GetAtomicNum() == OBElements::Carbon) {
        if(!(d->options & drawAllC))
        {
          if (atom->GetExplicitDegree() > 1)
            continue;
          if ((atom->GetExplicitDegree() == 1) && !(d->options & drawTermC))//!d->drawTerminalC)
            continue;
        }
      }

      bool written = false;
      stringstream ss;

      //For unexpanded aliases use appropriate form of alias instead of element symbol, Hs, etc
      AliasData* ad = NULL;
      if (d->aliasMode && atom->HasData(AliasDataType))
        ad = static_cast<AliasData*>(atom->GetData(AliasDataType));
      if(ad && !ad->IsExpanded())
      {
        ss <<ad->GetAlias(rightAligned);
        OBColor aliasColor = !ad->GetColor().empty() ? ad->GetColor() : d->bondColor;
          d->painter->SetPenColor(aliasColor);
        written = true;
      }

      if (!written) {
        //Atoms with no AliasData, but 0 atomic num and atomclass==n are output as Rn
        if (atom->GetAtomicNum()==0) {
          OBGenericData *data = atom->GetData("Atom Class");
          if (data) {
            OBPairInteger* acdata = dynamic_cast<OBPairInteger*>(data); // Could replace with C-style cast if willing to live dangerously
            if (acdata) {
              int ac = acdata->GetGenericValue();
              if (ac >= 0) {
                ss << 'R' << ac;
                d->painter->SetPenColor(OBColor("black"));
                written = true;
              }
            }
          }
        }
      }

      if (!written) {
        const char* atomSymbol;
        if(atom->GetAtomicNum() == OBElements::Hydrogen && atom->GetIsotope()>1)
          atomSymbol = atom->GetIsotope()==2 ? "D" : "T";
        else
          atomSymbol = OBElements::GetSymbol(atom->GetAtomicNum());

        unsigned int hCount = atom->GetImplicitHCount();
        // LPW: The allExplicit option will omit the drawing of extra hydrogens
        // to fill the valence.
        if((d->options & allExplicit))
            hCount = 0;
        // rightAligned:
        //   false  CH3
        //   true   H3C
        if (hCount && rightAligned)
          ss << "H";
        if ((hCount > 1) && rightAligned)
          ss << hCount;
        ss << atomSymbol;
        if (hCount && !rightAligned)
          ss << "H";
        if ((hCount > 1) && !rightAligned)
          ss << hCount;
      }
      d->DrawAtomLabel(ss.str(), alignment, vector3(x, y, 0.0));
    }

    return true;
  }

  void OBDepictPrivate::DrawWobblyBond(OBAtom *beginAtom, OBAtom *endAtom)
  {
    vector3 begin = beginAtom->GetVector();
    vector3 end = endAtom->GetVector();
    vector3 vb = end - begin;

    if (HasLabel(beginAtom))
      begin += 0.33 * vb;
    if (HasLabel(endAtom))
      end -= 0.33 * vb;

    vb = end - begin; // Resize the extents of the vb vector

    vector3 orthogonalLine = cross(vb, VZ);
    orthogonalLine.normalize();
    orthogonalLine *= 0.5 * bondWidth;

    double lines[6] = { 0.20, 0.36, 0.52, 0.68, 0.84, 1.0 };

    // This code is adapted from DrawWedge():
    // What we do is just join up the opposite ends of each of the wedge strokes
    // to create a zig-zag bond

    double oldx, oldy, newx, newy;
    oldx = begin.x();
    oldy = begin.y();
    int sign = 1;
    for (int k = 0; k < 6; ++k) {
      double w = lines[k];
      newx = begin.x() + vb.x() * w + sign * orthogonalLine.x() * w;
      newy = begin.y() + vb.y() * w + sign * orthogonalLine.y() * w;
      painter->DrawLine(oldx, oldy, newx, newy);
      oldx = newx;
      oldy = newy;
      sign = -sign;
    }
  }

  void OBDepictPrivate::DrawWedge(OBAtom *beginAtom, OBAtom *endAtom)
  {
    vector3 begin = beginAtom->GetVector();
    vector3 end = endAtom->GetVector();
    vector3 vb = end - begin;

    if (HasLabel(beginAtom))
      begin += 0.33 * vb;
    if (HasLabel(endAtom))
      end -= 0.33 * vb;

    vector3 orthogonalLine = cross(vb, VZ);
    orthogonalLine.normalize();
    orthogonalLine *= 0.5 * bondWidth;
    std::vector<std::pair<double,double> > points;

    points.push_back(std::pair<double,double>(begin.x(), begin.y()));
    points.push_back(std::pair<double,double>(end.x() + orthogonalLine.x(),
                                              end.y() + orthogonalLine.y()));
    points.push_back(std::pair<double,double>(end.x() - orthogonalLine.x(),
                                              end.y() - orthogonalLine.y()));
    painter->DrawPolygon(points);
  }

  void OBDepictPrivate::DrawHash(OBAtom *beginAtom, OBAtom *endAtom)
  {
    vector3 begin = beginAtom->GetVector();
    vector3 end = endAtom->GetVector();
    vector3 vb = end - begin;

    if (HasLabel(beginAtom))
      begin += 0.33 * vb;
    if (HasLabel(endAtom))
      end -= 0.33 * vb;

    vb = end - begin; // Resize the extents of the vb vector

    vector3 orthogonalLine = cross(vb, VZ);
    orthogonalLine.normalize();
    orthogonalLine *= 0.5 * bondWidth;

    double lines[6] = { 0.20, 0.36, 0.52, 0.68, 0.84, 1.0 };
    double oldwidth = painter->GetPenWidth();
    painter->SetPenWidth(1);
    for (int k = 0; k < 6; ++k) {
      double w = lines[k];
      painter->DrawLine(begin.x() + vb.x() * w + orthogonalLine.x() * w,
                        begin.y() + vb.y() * w + orthogonalLine.y() * w,
                        begin.x() + vb.x() * w - orthogonalLine.x() * w,
                        begin.y() + vb.y() * w - orthogonalLine.y() * w);
    }
    painter->SetPenWidth(oldwidth);
  }

  void OBDepictPrivate::DrawSimpleBond(OBAtom *beginAtom, OBAtom *endAtom, int order, bool crossed_dbl_bond)
  {
    vector3 begin = beginAtom->GetVector();
    vector3 end = endAtom->GetVector();
    vector3 vb = end - begin;

    vb.normalize();

    if (HasLabel(beginAtom))
      begin += 13. * vb; // Length is normally 40
    if (HasLabel(endAtom))
      end -= 13. * vb;

    if (order == 1) {
      painter->DrawLine(begin.x(), begin.y(), end.x(), end.y());
    } else if (order == 2) {
      vector3 orthogonalLine = cross(end - begin, VZ).normalize();

      bool useAsymmetricDouble = options & OBDepict::asymmetricDoubleBond;
      if (HasLabel(beginAtom) && HasLabel(endAtom))
        useAsymmetricDouble = false;
      if (HasLabel(beginAtom) && endAtom->GetExplicitDegree() == 3)
        useAsymmetricDouble = false;
      if (HasLabel(endAtom) && beginAtom->GetExplicitDegree() == 3)
        useAsymmetricDouble = false;
      if (crossed_dbl_bond)
        useAsymmetricDouble = false; // The bond looks very strange otherwise in the case of cis

      if (!useAsymmetricDouble) {
        // style1
        //
        // -----------
        // -----------
        vector3 offset = orthogonalLine * 0.5 * bondSpacing;
        if (!crossed_dbl_bond) {
          painter->DrawLine(begin.x() + offset.x(), begin.y() + offset.y(),
                            end.x() + offset.x(), end.y() + offset.y());
          painter->DrawLine(begin.x() - offset.x(), begin.y() - offset.y(),
                            end.x() - offset.x(), end.y() - offset.y());
        }
        else {
          painter->DrawLine(begin.x() + offset.x(), begin.y() + offset.y(),
                            end.x() - offset.x(), end.y() - offset.y());
          painter->DrawLine(begin.x() - offset.x(), begin.y() - offset.y(),
                            end.x() + offset.x(), end.y() + offset.y());
        }
      } else {
        // style2
        //
        //   -------
        // -----------
        vector3 offset1 = orthogonalLine * /*0.5 * */ bondSpacing;
        vector3 offset2 = vb * /*0.5 * */ bondSpacing;
        vector3 offset3 = - vb * /*0.5 * */ bondSpacing;

        if (HasLabel(beginAtom))
          offset2 = VZero;
        if (HasLabel(endAtom))
          offset3 = VZero;

        painter->DrawLine(begin.x(), begin.y(), end.x(), end.y());
        painter->DrawLine(begin.x() + offset1.x() + offset2.x(),
                          begin.y() + offset1.y() + offset2.y(),
                          end.x() + offset1.x() + offset3.x(),
                          end.y() + offset1.y() + offset3.y());
      }
    } else if (order == 3) {
      vector3 orthogonalLine = cross(end - begin, VZ).normalize();
      vector3 offset = orthogonalLine * 0.7 * bondSpacing;
      painter->DrawLine(begin.x(), begin.y(), end.x(), end.y());
      painter->DrawLine(begin.x() + offset.x(), begin.y() + offset.y(),
                        end.x() + offset.x(), end.y() + offset.y());
      painter->DrawLine(begin.x() - offset.x(), begin.y() - offset.y(),
                        end.x() - offset.x(), end.y() - offset.y());
    }
  }

  void OBDepictPrivate::DrawRingBond(OBAtom *beginAtom, OBAtom *endAtom, const vector3 &center, int order)
  {
    if (order != 2) {
      DrawSimpleBond(beginAtom, endAtom, order);
      return;
    }

    vector3 begin = beginAtom->GetVector();
    vector3 end = endAtom->GetVector();

    vector3 vb = (end - begin).normalize();
    vector3 orthogonalLine = cross(vb, VZ)/*.normalize()*/;
    vector3 spacing = orthogonalLine * bondSpacing * 1.2;
    vector3 offset = vb * bondSpacing;
    if ((begin + spacing - center).length() > (begin - spacing - center).length())
      spacing *= -1.0;

    vector3 vbb = end - begin;
    if (HasLabel(beginAtom))
      begin += 0.33 * vbb;
    if (HasLabel(endAtom))
      end -= 0.33 * vbb;
    painter->DrawLine(begin.x(), begin.y(), end.x(), end.y());

    if (HasLabel(beginAtom))
      begin -= 0.10 * vbb;
    if (HasLabel(endAtom))
      end += 0.10 * vbb;
    painter->DrawLine(begin.x() + spacing.x() + offset.x(), begin.y() + spacing.y() + offset.y(),
                      end.x() + spacing.x() - offset.x(), end.y() + spacing.y() - offset.y());
  }

  void OBDepictPrivate::DrawAtom(OBAtom *atom)
  {
  }

  void OBDepictPrivate::DrawAtomLabel(const std::string &label, int alignment, const vector3 &pos)
  {
   /*
    cout << "FontMetrics(" << label << "):" << endl;
    cout << "  ascent = " << metrics.ascent << endl;
    cout << "  descent = " << metrics.descent << endl;
    cout << "  width = " << metrics.width << endl;
    cout << "  height = " << metrics.height << endl;

    painter->SetFillColor(OBColor("white"));
    painter->SetPenColor(OBColor("white"));
    painter->DrawCircle(pos.x(), pos.y(), metrics.ascent / 2);
    painter->SetPenColor(OBColor("black"));
    */

    // compute the total width
    double totalWidth = 0.0;
    if ((alignment == Right) || (alignment == Left) || (label.find("H") == std::string::npos)) {
      for (unsigned int i = 0; i < label.size(); ++i) {
        if (!isalpha(label[i])) {
          painter->SetFontSize(subscriptSize);
          totalWidth += painter->GetFontMetrics(label.substr(i, 1)).width;
        } else {
          painter->SetFontSize(fontSize);
          totalWidth += painter->GetFontMetrics(label.substr(i, 1)).width;
        }
      }
    } else {
      painter->SetFontSize(fontSize);
      totalWidth = painter->GetFontMetrics(label.substr(0, label.find("H"))).width;
      double width = 0.0;
      for (unsigned int i = label.find("H"); i < label.size(); ++i) {
        if (!isalpha(label[i])) {
          painter->SetFontSize(subscriptSize);
          width += painter->GetFontMetrics(label.substr(i, 1)).width;
        } else {
          painter->SetFontSize(fontSize);
          width += painter->GetFontMetrics(label.substr(i, 1)).width;
        }
      }

      if (width > totalWidth)
        totalWidth = width;
    }

    painter->SetFontSize(fontSize);
    OBFontMetrics metrics = painter->GetFontMetrics(label);


    std::string str, subscript;
    // compute the horizontal starting position
    double xOffset, yOffset, yOffsetSubscript;
    switch (alignment) {
      case Right:
        xOffset = 0.5 * painter->GetFontMetrics(label.substr(0, 1)).width -
                  painter->GetFontMetrics(label).width;
        break;
      case Left:
        xOffset = - 0.5 * painter->GetFontMetrics(label.substr(label.size()-1, 1)).width;
        break;
      case Up:
      case Down:
        if (label.find("H") != std::string::npos)
          xOffset = - 0.5 * painter->GetFontMetrics(label.substr(0, label.find("H"))).width;
        else
          xOffset = - 0.5 * totalWidth;
        break;
      default:
        xOffset = - 0.5 * totalWidth;
        break;
    }

    // compute the vertical starting position
    yOffset = 0.5 * (metrics.ascent /*- metrics.descent*/);
    yOffsetSubscript = yOffset - metrics.descent;
    double xInitial = xOffset;

    for (unsigned int i = 0; i < label.size(); ++i) {
      if (label[i] == 'H') {
        if ((alignment == Up) || (alignment == Down))
          if (!str.empty()) {
            // write the current string
            painter->SetFontSize(fontSize);
            painter->DrawText(pos.x() + xOffset, pos.y() + yOffset, str);
            if (alignment == Down) {
              yOffset += metrics.fontSize;
              yOffsetSubscript += metrics.fontSize;
            } else {
              yOffset -= metrics.fontSize;
              yOffsetSubscript -= metrics.fontSize;
            }
            xOffset = xInitial;
            str.clear();
          }
      }


      if (!isalpha(label[i])) {
        if (!str.empty()) {
          // write the current string
          painter->SetFontSize(fontSize);
          OBFontMetrics metrics = painter->GetFontMetrics(str);
          painter->DrawText(pos.x() + xOffset, pos.y() + yOffset, str);
          xOffset += metrics.width;
          str.clear();
        }

        subscript += label.substr(i, 1);
      } else {
        if (!subscript.empty()) {
          // write the current subscript
          painter->SetFontSize(subscriptSize);
          OBFontMetrics metrics = painter->GetFontMetrics(subscript);
          painter->DrawText(pos.x() + xOffset, pos.y() + yOffsetSubscript, subscript);
          xOffset += metrics.width;
          subscript.clear();
        }

        str += label.substr(i, 1);
      }
    }
    if (!str.empty()) {
      painter->SetFontSize(fontSize);
      OBFontMetrics metrics = painter->GetFontMetrics(str);
      painter->DrawText(pos.x() + xOffset, pos.y() + yOffset, str);
    }
    if (!subscript.empty()) {
      painter->SetFontSize(subscriptSize);
      OBFontMetrics metrics = painter->GetFontMetrics(subscript);
      double yOffset = ispunct(subscript[subscript.size()-1]) || ispunct(subscript[0])
        || (subscript.size()>1 && ispunct(subscript[1]))
        ? -yOffsetSubscript : yOffsetSubscript;
      painter->DrawText(pos.x() + xOffset, pos.y() + yOffset, subscript);
    }

  }

  bool OBDepictPrivate::HasLabel(OBAtom *atom)
  {
    if (atom->GetAtomicNum() != OBElements::Carbon)
      return true;
    if ((options & OBDepict::drawAllC) || ((options & OBDepict::drawTermC) && (atom->GetExplicitDegree() == 1)))
      return true;
    return false;
  }

  void OBDepictPrivate::SetWedgeAndHash(OBMol* mol)
  {
    // Remove any existing wedge and hash bonds
    FOR_BONDS_OF_MOL(b,mol)  {
      b->SetWedge(false);
      b->SetHash(false);
    }

    std::map<OBBond*, enum OBStereo::BondDirection> updown;
    std::map<OBBond*, OBStereo::Ref> from;
    std::map<OBBond*, OBStereo::Ref>::const_iterator from_cit;
    TetStereoToWedgeHash(*mol, updown, from);

    for(from_cit=from.begin();from_cit!=from.end();++from_cit) {
      OBBond* pbond = from_cit->first;
      if(updown[pbond]==OBStereo::UpBond)
        pbond->SetHash();
      else if(updown[pbond]==OBStereo::DownBond)
        pbond->SetWedge();
    }
  }
  void OBDepictPrivateBallAndStick::DrawSimpleBond(OBAtom* beginAtom,
  OBAtom* endAtom, int order, bool crossed_bond)
  {
    const vector3 begin = beginAtom->GetVector();
    const vector3 end = endAtom->GetVector();
    const vector3 vb = (end - begin).normalize();

    // calculate the expected opacity and width
    // to simulate perspective
    penWidth = 3.0;
    bondColor.alpha = 1.0;

    if (fabs(zScale) > 1.0e-1) {
      double beginAtomScale = (beginAtom->GetZ() - zMin) / zScale;
      double endAtomScale = (endAtom->GetZ() - zMin) / zScale;
      double averageScale = (beginAtomScale + endAtomScale)/2.0;
      if (averageScale < 0.15)
        averageScale = 0.15;

      penWidth = 3.0 * averageScale;
      bondColor.alpha = averageScale;
    }
    painter->SetPenWidth(penWidth);

    if (order == 1) {
      painter->DrawLine(begin.x(), begin.y(), end.x(), end.y());
    } else if (order == 2) {
      const vector3 orthogonalLine = cross(end - begin, VZ).normalize();

        // style1
        //
        // -----------
        // -----------
        const vector3 offset = orthogonalLine * 0.5 * bondSpacing;
        painter->DrawLine(begin.x() + offset.x(), begin.y() + offset.y(),
                            end.x() + offset.x(), end.y() + offset.y());
        painter->DrawLine(begin.x() - offset.x(), begin.y() - offset.y(),
                            end.x() - offset.x(), end.y() - offset.y());
    } else if (order == 3) {
      const vector3 orthogonalLine = cross(end - begin, VZ).normalize();
      const vector3 offset = orthogonalLine * 0.7 * bondSpacing;
      painter->DrawLine(begin.x(), begin.y(), end.x(), end.y());
      painter->DrawLine(begin.x() + offset.x(), begin.y() + offset.y(),
                        end.x() + offset.x(), end.y() + offset.y());
      painter->DrawLine(begin.x() - offset.x(), begin.y() - offset.y(),
                        end.x() - offset.x(), end.y() - offset.y());
    }

  }

  inline void OBDepictPrivateBallAndStick::DrawWedge(OBAtom* beginAtom,
  OBAtom* endAtom)
  {
    const vector3 begin = beginAtom->GetVector();
    const vector3 end = endAtom->GetVector();
    const vector3 vb = end - begin;
    const vector3 orthogonalLine = cross(vb, VZ).normalize() * 0.5 * bondWidth;
    std::vector<std::pair<double,double> > points;

    points.push_back(std::pair<double,double>(begin.x(), begin.y()));
    points.push_back(std::pair<double,double>(end.x() + orthogonalLine.x(),
                                              end.y() + orthogonalLine.y()));
    points.push_back(std::pair<double,double>(end.x() - orthogonalLine.x(),
                                              end.y() - orthogonalLine.y()));
    painter->DrawPolygon(points);

  }

  inline void OBDepictPrivateBallAndStick::DrawHash(OBAtom* beginAtom,
  OBAtom* endAtom)
  {
    const vector3 begin = beginAtom->GetVector();
    const vector3 end = endAtom->GetVector();
    const vector3 vb = end - begin;
    const vector3 orthogonalLine = cross(vb, VZ).normalize() * 0.5 * bondWidth;
    const double lines[6] = { 0.20, 0.36, 0.52, 0.68, 0.84, 1.0 };
    const double oldwidth = painter->GetPenWidth();
    painter->SetPenWidth(1);
    for (int k = 0; k < 6; ++k) {
      const double w = lines[k];
      painter->DrawLine(begin.x() + vb.x() * w + orthogonalLine.x() * w,
                        begin.y() + vb.y() * w + orthogonalLine.y() * w,
                        begin.x() + vb.x() * w - orthogonalLine.x() * w,
                        begin.y() + vb.y() * w - orthogonalLine.y() * w);
    }
    painter->SetPenWidth(oldwidth);

  }

  inline void OBDepictPrivateBallAndStick::DrawWobblyBond(OBAtom* beginAtom,
  OBAtom* endAtom)
  {
    DrawSimpleBond(beginAtom, endAtom, 1);
  }

  inline void OBDepictPrivateBallAndStick::DrawRingBond(OBAtom* beginAtom,
  OBAtom* endAtom, const vector3& center, int order)
  {
    OBDepictPrivate::DrawRingBond(beginAtom, endAtom, center, order);
  }


  inline void OBDepictPrivateBallAndStick::DrawRing(OBRing* ring,
  OBBitVec& drawnBonds)
  {
    OBDepictPrivate::DrawRing(ring, drawnBonds);
  }

inline void OBDepictPrivateBallAndStick::DrawAromaticRing(OBRing* ring,
OBBitVec& drawnBonds)
  {

    const std::vector<int> indexes = ring->_path;
    const size_t ringSize = indexes.size();
    vector3 center(VZero);
    double maxdist =0.;
    for (std::vector<int>::const_iterator l = indexes.begin(); l != indexes.end(); ++l) {
      center += mol->GetAtom(*l)->GetVector();
      maxdist = max(maxdist, GetAtomRadius(mol->GetAtom(*l)));
    }
    center /= ringSize;

    for (unsigned int l = 0; l < indexes.size(); ++l) {
      OBAtom *prev  = mol->GetAtom(indexes[l]);
      OBAtom *begin = mol->GetAtom(indexes[(l+1) % ringSize]);
      OBAtom *end   = mol->GetAtom(indexes[(l+2) % ringSize]);
      OBAtom *next  = mol->GetAtom(indexes[(l+3) % ringSize]);

      OBBond *ringBond = mol->GetBond(begin, end);

      // calculate the expected opacity and width
      // to simulate perspective
      penWidth = 3.0;
      bondColor.alpha = 1.0;

      if (fabs(zScale) > 1.0e-1) {
        double beginAtomScale = (begin->GetZ() - zMin) / zScale;
        double endAtomScale = (end->GetZ() - zMin) / zScale;
        double averageScale = (beginAtomScale + endAtomScale)/2.0;
        if (averageScale < 0.15)
          averageScale = 0.15;

        penWidth = 3.0 * averageScale;
        bondColor.alpha = averageScale;
      }
      painter->SetPenWidth(penWidth);

      if((options & OBDepict::internalColor) && ringBond->HasData("color"))
        painter->SetPenColor(OBColor(ringBond->GetData("color")->GetValue()));
      else
        painter->SetPenColor(bondColor);

      DrawAromaticRingBond(prev,begin, end, next, center, maxdist);
      drawnBonds.SetBitOn(ringBond->GetId());
    }
  }

  inline void OBDepictPrivateBallAndStick::DrawAromaticRingBond(OBAtom *prevAtom, OBAtom *beginAtom, OBAtom *endAtom, OBAtom *nextAtom, const vector3 &center, double dist)
  {
    const vector3 prev  = prevAtom->GetVector();
    const vector3 begin = beginAtom->GetVector();
    const vector3 end   = endAtom->GetVector();
    const vector3 next  = nextAtom->GetVector();

    const vector3 orthogonalLine = cross(end - begin, VZ).normalize();
    const vector3 offset = orthogonalLine * 0.5 * bondSpacing;
    painter->DrawLine(begin.x() - offset.x(), begin.y() - offset.y(),
                      end.x() - offset.x(), end.y() - offset.y());

    static const float dashpattern[] = {5., 5.};
    static const vector<double> pat = vector<double>(dashpattern,dashpattern + sizeof(dashpattern)/sizeof(double));
    painter->DrawLine(begin.x() + offset.x(), begin.y() + offset.y(),
                      end.x() + offset.x(), end.y() + offset.y(), pat);
  }

  void OBDepictPrivateBallAndStick::DrawAtom(OBAtom *atom)
  {
    double r, g, b;
    OBElements::GetRGB(atom->GetAtomicNum(), &r, &g, &b);
    OBColor atomColor = OBColor(r, g, b);
    double opacity = 1.0;
    if (fabs(zScale) > 1.0e-1)
      opacity = sqrt((atom->GetZ() - zMin) / zScale);
    if (opacity < 0.2)
      opacity = 0.2;

    painter->SetFillRadial(OBColor("white"),atomColor);
    painter->DrawBall(atom->GetVector().x(), atom->GetVector().y(),GetAtomRadius(atom), opacity);
  }

  double OBDepictPrivateBallAndStick::GetAtomRadius(OBAtom *atom)
  {
    double radius = OBElements::GetCovalentRad(atom->GetAtomicNum());
    double perspective = 1.0;
    if (fabs(zScale) > 1.0e-1)
      perspective = (atom->GetZ() - zMin) / zScale;
    if (perspective < 0.5)
      perspective = 0.5;

    return perspective * radius * bondLength / 1.1;
  }

  void OBDepictPrivateBallAndStick::DrawAtomLabel(const std::string &label, int alignment, const vector3 &pos)
    {
      if (m_symbolOnBall)
        OBDepictPrivate::DrawAtomLabel(label,alignment,pos);
    }

}

/// @file depict.cpp
/// @brief 2D depiction of molecules using OBPainter.
