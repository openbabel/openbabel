#ifndef OB_CAIROPAINTER_H
#define OB_CAIROPAINTER_H

#include <openbabel/depict/painter.h>
#include <cairo.h>

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
      void WriteImage(std::ostream& ofs, int newWidth=300, int newHeight=300);
      //@}
 
    private:
      cairo_surface_t *m_surface;
      cairo_t *m_cairo;
      int m_fontPointSize;
  };

}

#endif
