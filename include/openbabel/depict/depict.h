/**********************************************************************
depict.h - 2D Depiction of molecules using OBPainter.

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
#ifndef OB_DEPICT_H
#define OB_DEPICT_H

#include <openbabel/babelconfig.h>
#include <openbabel/math/vector3.h>

namespace OpenBabel
{

  class OBMol;
  class OBPainter;
  class OBDepictPrivate;

  #ifndef OBDEPICT
    #define OBDEPICT
  #endif

  /**
   * @class OBDepict depict.h <openbabel/depict/depict.h>
   * @brief 2D depiction of molecules using OBPainter.
   * @since version 2.3
   */
  class OBDEPICT OBDepict
  {
    public:
      enum AtomLabelType {
        AtomId = 1,
        AtomIndex,
        AtomSymmetryClass,
        AtomValence,
        AtomTetrahedralStereo
      };

      enum OptionType{
        bwAtoms              = 0x0001,
        internalColor        = 0x0002,
        noMargin             = 0x0004,
        drawTermC            = 0x0010,
        drawAllC             = 0x0020,
        noWedgeHashGen       = 0x0100,
        asymmetricDoubleBond = 0x0200,
        allExplicit          = 0x0400
      };

      /**
       * Constructors.
       */
      OBDepict(OBPainter *painter);
      OBDepict(OBPainter *painter, bool withBall, bool symbolOnBall=false);
      /**
       * Destructor.
       */
      ~OBDepict();
      /**
       * Draw @p mol using the painter previously stored in the constructor.
       *
       * @return True if successful.
       */
      bool DrawMolecule(OBMol *mol);
      /**
       * Draw atom labels of a specified @p type.
       *
       * @return True if successful.
       */
      bool AddAtomLabels(AtomLabelType type);

      void SetBondLength(double length);
      double GetBondLength() const;

      void SetPenWidth(double length);
      double GetPenWidth() const;

      void SetBondSpacing(double spacing);
      double GetBondSpacing() const;

      void SetBondWidth(double width);
      double GetBondWidth() const;

      //void SetDrawingTerminalCarbon(bool enabled);
      //bool GetDrawingTerminalCarbon() const;

      void SetOption(unsigned opts); //extendable with binary compatibility
      unsigned GetOptions() const;
      void ClearOptions();

      void SetFontFamily(const std::string &family);
      const std::string& GetFontFamily() const;

      void SetFontSize(int pointSize, bool subscript = false);
      int GetFontSize(bool subscript = false) const;

      void SetAliasMode(bool b=true);

      void SetBondColor(const std::string& scolor);

    private:
      OBDepictPrivate * d;
  };

}

#endif

/// @file depict.h
/// @brief 2D depiction of molecules using OBPainter.
