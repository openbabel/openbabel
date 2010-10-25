/**********************************************************************
svgpainter.h - Rendering in SVG

Copyright (C) 2009 by Chris Morley

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
#ifndef OB_SVGPAINTER_H
#define OB_SVGPAINTER_H

#include <openbabel/babelconfig.h>
#include <iostream>
#include <openbabel/depict/painter.h>

namespace OpenBabel
{
  class OBDEPICT SVGPainter : public OBPainter
  {
    public:
      SVGPainter();
      SVGPainter(std::ostream& ofs, bool withViewBox=false,
        double width=0.0, double height=0.0, double x=0.0, double y=0.0);
      ~SVGPainter();
      //! @name OBPainter methods
      //@{
      void NewCanvas(double width, double height);
      bool IsGood() const;
      void SetFontFamily(const std::string &fontFamily);
      void SetFontSize(int pointSize);
      void SetFillColor(const OBColor &color);
      void SetPenColor(const OBColor &color);
      void SetPenWidth(double width);
      void DrawLine(double x1, double y1, double x2, double y2);
      void DrawPolygon(const std::vector<std::pair<double,double> > &points);
      void DrawCircle(double x, double y, double r);
      void DrawText(double x, double y, const std::string &text);
      OBFontMetrics GetFontMetrics(const std::string &text);
      //@}

      //! @name CairoPainter specific
      //@{
      void WriteImage(const std::string &filename);
      //@}
    private:
      std::string MakeRGB(OBColor color);

    private:
      std::ostream& m_ofs;
      bool m_withViewBox;
      double m_width, m_height, m_x, m_y;
      OBColor m_Pencolor;
      OBColor m_OrigBondcolor;
      OBColor m_Fillcolor;
      double m_PenWidth;
      int m_fontPointSize;
      std::string m_fontFamily;
  };

}

#endif

//! \file svgpainter.h
//! \brief Generate 2D depictions in the SVG vector graphics format.
