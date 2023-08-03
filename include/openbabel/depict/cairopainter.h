#ifndef OB_CAIROPAINTER_H
#define OB_CAIROPAINTER_H

#include <openbabel/depict/painter.h>
#include "cairo.h"

namespace OpenBabel
{
  class CairoPainter : public OBPainter
  {
    public:
      CairoPainter();
      ~CairoPainter();
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
      void DrawBall(double x, double y, double r, double opacity);
      void DrawText(double x, double y, const std::string &text);
      OBFontMetrics GetFontMetrics(const std::string &text);
      //@}

      //! @name CairoPainter specific
      //@{
      void WriteImage(const std::string &filename);
      void WriteImage(std::ostream& ofs);
      void SetWidth(int width) {m_width=width;}
      void SetHeight(int height) {m_height=height;}
      void SetTitle(std::string title) {m_title=title;}
      void SetIndex(int index) {m_index=index;}
      void SetTableSize(int nrows, int ncols) {m_nrows=nrows; m_ncols=ncols;}
      void SetBackground(std::string color) {m_fillcolor=color;}
      void SetBondColor(std::string color) {m_bondcolor=color;}
      void SetTransparent(bool tr) {m_transparent=tr;}
      void SetCropping(bool cr) {m_cropping=cr;}
     //@}

    private:
      cairo_surface_t *m_surface;
      cairo_t *m_cairo;
      int m_fontPointSize;
      int m_width;
      int m_height;
      double m_pen_width;
      std::string m_title;
      int m_index;        // Index of current molecule in a table
      int m_ncols, m_nrows; // number of rows and columns
      std::string m_fillcolor, m_bondcolor;
      bool m_transparent, m_cropping;

  };

}

#endif
