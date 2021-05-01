#include <openbabel/obutil.h>
#include <openbabel/depict/cairopainter.h>

#include <iostream>
using namespace std;

namespace OpenBabel
{

  // Class definition of CairoPainter
  CairoPainter::CairoPainter() : m_surface(nullptr), m_cairo(nullptr),
    m_fontPointSize(12), m_width(0), m_height(0), m_pen_width(1), m_title(""), m_index(1),
    m_fillcolor("white"), m_bondcolor("black"), m_transparent(false)
  {
  }

  CairoPainter::~CairoPainter()
  {
    if (m_cairo)
      cairo_destroy(m_cairo);
    if (m_surface)
      cairo_surface_destroy(m_surface);
  }

  void CairoPainter::NewCanvas(double width, double height)
  {
    double titleheight = m_title.empty() ? 0.0 : 16.0;
    if (m_index == 1) {
      // create new surface to paint on
      if(m_cropping) {
        double ratio = width / height;
        if(ratio > 1.0)
          m_height = m_height / ratio;
        else
          m_width = m_width * ratio;
      }
      m_surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, static_cast<int> (m_width), static_cast<int> (m_height));
      m_cairo = cairo_create(m_surface);
      if(m_transparent)
        cairo_set_source_rgba (m_cairo, 0.0, 0.0, 0.0, 0.0);
      else {
        OBColor bg = OBColor(m_fillcolor);
        cairo_set_source_rgb (m_cairo, bg.red, bg.green, bg.blue);
      }

      cairo_paint (m_cairo);
      cairo_set_line_width(m_cairo, m_pen_width);
    }
    else {
      // reset transformation matrix
      cairo_identity_matrix(m_cairo);
    }

    // Work out some things!
    double cellwidth = m_width/m_ncols;
    double cellheight = m_height/m_nrows;
    int row = (m_index - 1)/m_ncols + 1;
    int col = m_index - ((row-1)*m_ncols);

    // Work out the scaling factor
    double scale_x = cellwidth / (double) width;
    double scale_y = (cellheight-titleheight) / (double) height; // Leave some extra space for the title if present
    double scale = std::min(scale_x, scale_y);

    // Add the title
    if (!m_title.empty()) {
      this->SetPenColor(OBColor(m_bondcolor));
      this->SetFontSize(static_cast<int>(16.0));
      OBFontMetrics fm = this->GetFontMetrics(m_title);
      this->DrawText(cellwidth/2.0 - fm.width/2.0 + cellwidth*(col-1),
                     cellheight - fm.height * 0.25 + cellheight*(row-1), m_title);
    }

    // Translate the over-scaled dimension into the centre
    if (scale < scale_y)
      cairo_translate(m_cairo, 0 + cellwidth*(col-1), cellheight/2.0 - scale*height/2.0 + cellheight*(row-1));
    else
      cairo_translate(m_cairo, cellwidth/2.0 - scale*width/2.0 + cellwidth*(col-1), 0 + cellheight*(row-1));
    cairo_scale(m_cairo, scale, scale); // Set a scaling transformation
  }

  bool CairoPainter::IsGood() const
  {
    if (!m_cairo)
      return false;
    if (!m_surface)
      return false;
    return true;
  }

  void CairoPainter::SetFontSize(int pointSize)
  {
    m_fontPointSize = pointSize;
    cairo_set_font_size(m_cairo, pointSize);
  }

  void CairoPainter::SetFillColor(const OBColor &color)
  {
    cairo_set_source_rgb(m_cairo, color.red, color.green, color.blue);
  }

  void CairoPainter::SetFillRadial(const OBColor &start, const OBColor &end)
  {
    cairo_set_source_rgb(m_cairo, end.red, end.green, end.blue);
  }

  void CairoPainter::SetPenColor(const OBColor &color)
  {
    cairo_set_source_rgb(m_cairo, color.red, color.green, color.blue);
  }

  void CairoPainter::SetPenWidth(double width)
  {
    m_pen_width = width;
  }

  double CairoPainter::GetPenWidth()
  {
    return m_pen_width;
  }

  void CairoPainter::DrawLine(double x1, double y1, double x2, double y2, const std::vector<double>& dashes)
  {
    cairo_set_line_width(m_cairo, m_pen_width);
    cairo_set_line_cap(m_cairo, CAIRO_LINE_CAP_ROUND);
    cairo_set_dash(m_cairo, (dashes.size() ? &dashes[0] : nullptr), dashes.size(), 0.0);
    cairo_move_to(m_cairo, x1, y1);
    cairo_line_to(m_cairo, x2, y2);
    cairo_stroke(m_cairo);
  }

  void CairoPainter::DrawPolygon(const std::vector<std::pair<double,double> > &points)
  {
    std::vector<std::pair<double,double> >::const_iterator i;
    for (i = points.begin(); i != points.end(); ++i)
      cairo_line_to(m_cairo, i->first, i->second); // note: when called without previous point,
                                                   //       this function behaves like cairo_move_to
    cairo_line_to(m_cairo, points.begin()->first, points.begin()->second);
    cairo_fill(m_cairo);
  }

  void CairoPainter::DrawCircle(double x, double y, double r)
  {
    cairo_arc(m_cairo, x, y, r, 0, 2 * M_PI);
    cairo_stroke(m_cairo);
  }

  void CairoPainter::DrawText(double x, double y, const std::string &text)
  {
    cairo_move_to(m_cairo, x, y);
    cairo_show_text(m_cairo, text.c_str());
  }

  OBFontMetrics CairoPainter::GetFontMetrics(const std::string &text)
  {
    cairo_font_extents_t fe;
    cairo_font_extents(m_cairo, &fe);
    cairo_text_extents_t te;
    cairo_text_extents(m_cairo, text.c_str(), &te);

    OBFontMetrics metrics;
    metrics.fontSize = m_fontPointSize;
    metrics.ascent = fe.ascent;
    metrics.descent = -fe.descent;
    metrics.width = te.x_advance;//te.width;
    metrics.height = te.height;
    return metrics;
  }

  void CairoPainter::WriteImage(const std::string &filename)
  {
    if (!m_cairo || !m_surface)
      return;

    cairo_surface_write_to_png(m_surface, filename.c_str());
  }

  static cairo_status_t writeFunction(void* closure, const unsigned char* data, unsigned int length)
  {
    vector<char>* in = reinterpret_cast<vector<char>*>(closure);
    for (unsigned int i = 0; i < length; ++i)
      in->push_back(data[i]);
    return CAIRO_STATUS_SUCCESS;
  }

  void CairoPainter::DrawBall(double x, double y, double r, double opacity)
  {

  }

  void CairoPainter::WriteImage(std::ostream& ofs)
  {
    if (!m_cairo || !m_surface)
      return;
    vector<char> in;
    cairo_surface_write_to_png_stream(m_surface, writeFunction, &in);
    for (unsigned int i = 0; i < in.size(); ++i)
      ofs << in.at(i);
  }

}
