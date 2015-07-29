/**********************************************************************
painter.h - Abstract base class for rendering

Copyright (C) 2009 by Tim Vandermeersch
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
#ifndef OB_PAINTER_H
#define OB_PAINTER_H

#include <openbabel/babelconfig.h>
#include <string>
#include <vector>
#include <sstream>
#ifndef OBDEPICT
  #define OBDEPICT
#endif

namespace OpenBabel
{
  /**
   * @class OBColor painter.h <openbabel/depict/painter.h>
   * @brief Color class used by OBDepict.
   * @since version 2.3
   */
  struct OBDEPICT OBColor
  {
  public:
    OBColor()
    {
      *this = OBColor(0.0, 0.0, 0.0);
    }
    OBColor(double _red, double _green, double _blue, double _alpha = 1.0) :
        red(_red), green(_green), blue(_blue), alpha(_alpha)
    {
    }
    OBColor(const std::string &color)
    {
      if (color[0]=='#')
      {
        std::stringstream ss(color.substr(1));
        unsigned c;
        ss >> std::hex >> c;
        *this = OBColor((c/0x10000)/256.0, ((c%0x10000)/0x100/256.0), (c%0x100)/256.0);
        return;
      }
      if (color == "black")
        *this = OBColor(0.0, 0.0, 0.0);
      else if (color == "white")
        *this = OBColor(1.0, 1.0, 1.0);
      else if (color == "red")
        *this = OBColor(1.0, 0.0, 0.0);
      else if (color == "green")
        *this = OBColor(0.0, 1.0, 0.0);
      else if (color == "blue")
        *this = OBColor(0.0, 0.0, 1.0);
      else if (color == "yellow")
        *this = OBColor(1.0, 1.0, 0.0);
      else if (color == "gray")
        *this = OBColor(0.3, 0.3, 0.3);
      else if (color == "cyan")
        *this = OBColor(1.0, 0.0, 1.0);
      else if (color == "purple")
        *this = OBColor(0.5, 0.0, 0.5);
      else if (color == "teal")
        *this = OBColor(0.0, 0.5, 0.5);
      else if (color == "olive")
        *this = OBColor(0.5, 0.5, 0.0);
      else if (color == "none")
        *this = OBColor(0.0, 0.0, 0.0, 0.0);
      else
        *this = OBColor(0.5, 0.5, 0.5);
    }

    OBColor(std::vector<double> vec) : red(vec[0]), green(vec[1]), blue(vec[2]), alpha(1.0){}

    bool operator !=(const OBColor& other)
    {
      return red!=other.red || green!=other.green || blue!=other.blue;
    }

    bool operator <(const OBColor& other) const
    {
      if (red < other.red)
        return true;
      else if (red > other.red)
        return false;
      else if (green < other.green)
        return true;
      else if (green > other.green)
        return false;
      else if (blue<other.blue)
        return true;
      return false;
    }

    double red, green, blue, alpha;
  };

  /**
   * @class OBFontMetrics painter.h <openbabel/depict/painter.h>
   * @brief Font metrics class used by OBDepict.
   * @since version 2.3
   */
  struct OBDEPICT OBFontMetrics
  {
    int    fontSize;
    double ascent, descent;
    double width, height;
  };

  /**
   * @class OBPainter painter.h <openbabel/depict/painter.h>
   * @brief Abstract painter base class used by OBDepict.
   * @since version 2.3
   */
  class OBDEPICT OBPainter
  {
    public:
      /**
       * Virtual destructor to be inhereted by subclasses
       */
      virtual ~OBPainter() {}

      /**
       * Create a new canvas to paint on with size @p width x @p height.
       * OBDepict will always call NewCanvas before performing any drawing
       * operations. Painters that are capable of drawing on a previously
       * unspecified area don't need to implement this.
       */
      virtual void NewCanvas(double width, double height) = 0;
      /**
       * Before OBDepict performes any drawing operation, this method is called
       * to check if the painter is ready to start drawing. If this method
       * returns false, drawing is aborted.
       */
      virtual bool IsGood() const = 0;
      /**
       * Set the painter's font family.
       */
      virtual void SetFontFamily(const std::string &fontFamily) = 0;
      /**
       * Set the painter's font point size.
       */
      virtual void SetFontSize(int pointSize) = 0;
      /**
       * Set the painter's fill color.
       */
      virtual void SetFillColor(const OBColor &color) = 0;
      /**
       * Set the painter's fill as a radial gradient.
       */
      virtual void SetFillRadial(const OBColor &start, const OBColor &end) = 0;
      /**
       * Set the painter's pen color.
       */
      virtual void SetPenColor(const OBColor &color) = 0;
      /**
       * Set the painter's pen width.
       */
      virtual void SetPenWidth(double width) = 0;
      /**
       * Get the painter's pen width.
       */
      virtual double GetPenWidth() = 0;
      /**
       * Draw a line from @p x1, @p y1 to @p x2, @p y2. The line is drawn using
       * the current pen color and width.
       */
      virtual void DrawLine(double x1, double y1, double x2, double y2, const std::vector<double> & dashes=std::vector<double>()) = 0;
      virtual void DrawCircle(double x, double y, double r) = 0;
      /**
       * Draw a polygon by connecting consecutive points. The last point will be
       * connected to the first one. The lines are drawn using the current pen
       * color and width. The area inside the polygon is filled with the current
       * fill color.
       */
      virtual void DrawPolygon(const std::vector<std::pair<double,double> > &points) = 0;

      /**
       * Draw a pseudo 3D ball, centered on (@p x,@p y) of radius @p r
       * filled with the current Fill color or gradient
       */
      virtual void DrawBall(double x, double y, double r, double opacity = 1.0) = 0;

      virtual void DrawText(double x, double y, const std::string &text) = 0;
      virtual OBFontMetrics GetFontMetrics(const std::string &text) = 0;
  };

}

#endif

//! \file painter.h
//! \brief Base class for graphical 2D depiction classes (e.g., "paint" to SVG)
