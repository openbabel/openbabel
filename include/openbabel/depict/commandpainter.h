/**********************************************************************
commandpainter.h  - Render as depiction commands

Copyright (C) 2012 by Noel O'Boyle

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
#ifndef OB_COMMANDPAINTER_H
#define OB_COMMANDPAINTER_H

#include <openbabel/depict/painter.h>

namespace OpenBabel
{
  class CommandPainter : public OBPainter
  {
    public:
      CommandPainter(std::ostream& ofs);
      ~CommandPainter();
      //! @name OBPainter methods
      //@{
      void NewCanvas(double width, double height);
      bool IsGood() const;
      void SetFontFamily(const std::string &fontFamily) {} // FIXME
      void SetFontSize(int pointSize);
      void SetFillColor(const OBColor &color);
      void SetFillRadial(const OBColor &start, const OBColor &end);
      void SetPenColor(const OBColor &color);
      void SetPenWidth(double width);
      double GetPenWidth();
      void DrawLine(double x1, double y1, double x2, double y2, const std::vector<double> & dashes=std::vector<double>());
      void DrawPolygon(const std::vector<std::pair<double,double> > &points);
      void DrawCircle(double x, double y, double r);
      void DrawBall(double x, double y, double r, double opacity = 1.0);
      void DrawText(double x, double y, const std::string &text);
      OBFontMetrics GetFontMetrics(const std::string &text);
      //@}

      //! @name CommandPainter specific
      //@{
      //@}

    private:
      std::ostream& m_ofs;
      std::streamsize m_oldprec;
      double m_pen_width;
  };

}

#endif
