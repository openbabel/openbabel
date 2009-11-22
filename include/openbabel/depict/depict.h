/**********************************************************************
depict.h - 2D Depiction of molecules

Copyright (C) 2009 by Tim Vandermeersch
Some portions Copyright (C) 2009 by Chris Morley

This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>

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

      /**
       * Constructor.
       */
      OBDepict(OBPainter *painter);
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

      void SetDrawingTerminalCarbon(bool enabled);
      bool GetDrawingTerminalCarbon() const;

      void SetFontFamily(const std::string &family); 
      const std::string& GetFontFamily() const;

      void SetFontSize(int pointSize, bool subscript = false);
      int GetFontSize(bool subscript = false) const;

      void SetAliasMode(bool b=true);

      void SetBondColor(const std::string& scolor);

    private:
      OBDepictPrivate * const d;
  };

}

#endif
