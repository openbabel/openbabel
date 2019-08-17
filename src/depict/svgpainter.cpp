/**********************************************************************
svgpainter.cpp  - Implementation of rendering in SVG

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
#include <openbabel/depict/svgpainter.h>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <math.h>
using namespace std;

#if defined(_MSC_VER) && _MSC_VER >= 1200 && _MSC_VER < 1800 // Between VC++ 6.0 and VC++ 11.0
#include <float.h>
#define isfinite _finite
#endif

namespace OpenBabel
{

  SVGPainter::SVGPainter(ostream& ofs, std::set<ColorGradient>* gradients, bool withViewBox,
    double width, double height)
    :  m_ofs(ofs), m_withViewBox(withViewBox), m_width(width), m_height(height),
       m_Pencolor("black"), m_Fillcolor("white"), m_Gradientcolor(make_pair(OBColor("white"),OBColor("white"))), m_PenWidth(1),
       m_fontPointSize(16), m_isFillcolor(true), m_Gradients(gradients)  {}

  SVGPainter::~SVGPainter()
  {

  }

  void SVGPainter::NewCanvas(double width, double height)
  {
    if(m_withViewBox)
      m_ofs << "<svg width=\"" << m_width << "\" height=\"" << m_height << "\" "
            << "x=\"0\" y=\"0\" "
            << "viewBox=\"0 0 " << width << ' ' << height << "\"\n";
    else
      m_ofs << "<svg width=\"" << width << "\" height=\"" << height << "\" "
            << "x=\"0\" y=\"0\" ";

    //Bond color and width are the initial m_Pencolor and m_PenWidth
    m_ofs << "font-family=\"" << m_fontFamily << "\" stroke=" << MakeRGB(m_Pencolor)
          << "stroke-width=\"" << m_PenWidth << "\"  stroke-linecap=\"round\"" << ">\n";

    if(!m_withViewBox && m_Fillcolor.alpha!=0.0)//Background color for single molecule. Handled by outer svg when table.
      m_ofs << "<rect x=\"0%\" y=\"0%\" width=\"100%\" height=\"100%\" stroke-width=\"0\" fill="
            << MakeRGB(m_Fillcolor) << " />\n";
    m_OrigBondcolor = m_Pencolor;
  }

  void SVGPainter::EndCanvas()
  {
    m_ofs << "</svg>\n";
  }
  bool SVGPainter::IsGood() const
  {
    return true;
  }

  void SVGPainter::SetFontSize(int pointSize)
  {
    m_fontPointSize = pointSize;
  }

  void SVGPainter::SetFontFamily(const std::string &fontFamily)
  {
    m_fontFamily = fontFamily;
  }

  void SVGPainter::SetFillColor(const OBColor &color)
  {
    m_Fillcolor = color; //value when NewCanvas called used for background
    m_isFillcolor = true;
  }

  void SVGPainter::SetFillRadial(const OBColor &start, const OBColor &end)
  {
    m_Gradientcolor = make_pair(start,end);
    m_Gradients->insert(m_Gradientcolor);
    m_isFillcolor = false;
  }

  void SVGPainter::SetPenColor(const OBColor &color)
  {
    m_Pencolor = color; //value when NewCanvas called used for bonds
  }

  void SVGPainter::SetPenWidth(double width)
  {
    m_PenWidth = width; //value when NewCanvas called used for bonds
  }

  double SVGPainter::GetPenWidth()
  {
    return m_PenWidth;
  }

  void SVGPainter::DrawLine(double x1, double y1, double x2, double y2, const std::vector<double>& dashes)
  {
    streamsize oldprec = m_ofs.precision(1);
    m_ofs << fixed << "<line x1=\"" << x1 << "\" y1=\"" << y1 << "\" x2=\""
      << x2 << "\" y2=\"" << y2 << "\"";
    m_ofs << " opacity=\"" << m_Pencolor.alpha << "\"";
    // if(m_Pencolor!=m_OrigBondcolor) // TODO: Bring this line back once Pybel is fine with this
      m_ofs << " stroke=" << MakeRGB(m_Pencolor);
    m_ofs << " stroke-width=\"" << m_PenWidth << "\"";
    if (!dashes.empty()) {
      std::vector<double>::const_iterator it = dashes.begin();
      m_ofs << " stroke-dasharray=\"" << *it;
      for (; it!=dashes.end() ; ++it)
        m_ofs << "," << *it;
      m_ofs << "\"";

    }
    m_ofs << "/>\n";
    m_ofs.precision(oldprec);
  }

  void SVGPainter::DrawPolygon(const std::vector<std::pair<double,double> > &points)
  {
    m_ofs << "<polygon points=\"";
      std::vector<std::pair<double,double> >::const_iterator i;
    for (i = points.begin(); i != points.end(); ++i)
      m_ofs << i->first << ' ' << i->second << ' ';
    m_ofs << "\"";
    m_ofs << " stroke-width=\"" << m_PenWidth << "\"";
    m_ofs << " fill=" << MakeRGB(m_Pencolor);
    m_ofs << " stroke=" << MakeRGB(m_Pencolor);
    m_ofs << "/>\n";
  }

  void SVGPainter::DrawCircle(double x, double y, double r)
  {
    m_ofs << "<circle cx=\"" << x << "\" y=\"" << y << "\" r=\"" << r << "\" />\n";
  }

  void SVGPainter::DrawText(double x, double y, const std::string &text)
  {
    m_ofs << "<text x=\"" << x << "\" y=\"" << y << "\""
      << " fill=" << MakeRGB(m_Pencolor) << "stroke-width='0' font-weight='bold' "
      << "font-size=\"" << m_fontPointSize << "\" >"
      << text << "</text>\n";
  }

  OBFontMetrics SVGPainter::GetFontMetrics(const std::string &text)
  {
    OBFontMetrics metrics;
    metrics.fontSize = m_fontPointSize;
    metrics.ascent   = m_fontPointSize;
    metrics.descent  = m_fontPointSize * -0.23; // Offset from baseline of bottom of text
    metrics.height   = m_fontPointSize *  1.26; // Distance between successive lines of text
    metrics.width    = 0.0;
    for(string::size_type i=0;i<text.size();++i)
      metrics.width += m_fontPointSize * (isalpha(text[i]) ? 0.75 : 0.5);

    return metrics;
  }

  void SVGPainter::WriteImage(const std::string &filename)
  {
  }

  void SVGPainter::DrawBall(double x, double y, double r, double opacity)
  {
    if (!isfinite(opacity))
      opacity = 1.0;
    if (opacity < 0.2)
      opacity = 0.2;

    m_ofs << "<circle cx=\"" << x << "\" cy=\"" << y << "\" r=\"" << r << "\" ";
    m_ofs << "opacity=\"" << opacity << "\" ";
    if (m_isFillcolor) {
      m_ofs << "style=\"stroke:black;fill:" << MakeRGB(m_Fillcolor) << "\"/>\n";
    } else {
      m_ofs << "style=\"stroke:black;stroke-width:0.5;fill:url(#radial";
      m_ofs << RGBcode(m_Gradientcolor.first)<< RGBcode(m_Gradientcolor.second) << ")\"/>\n";
    }
  }

  void   SVGPainter::WriteDefs()
  {
    if (!m_Gradients->empty()) {
      m_ofs << "<defs>\n";
      for (std::set<ColorGradient>::iterator it=m_Gradients->begin(); it!=m_Gradients->end(); ++it) {
        m_ofs << "<radialGradient id='radial";
        m_ofs << RGBcode(it->first)<< RGBcode(it->second) << "' ";
        m_ofs << "cx='50%' cy='50%' r='50%' fx='30%' fy='30%'>\n";
        m_ofs << "  <stop offset=' 0%' stop-color=" << MakeRGB(it->first) << " stop-opacity='1.0'/>\n";
        m_ofs << "  <stop offset='100%' stop-color=" << MakeRGB(it->second) << " stop-opacity ='1.0'/>\n";
        m_ofs << "</radialGradient>\n";
      }
      m_ofs << "</defs>\n";
    }
  }

  string SVGPainter::MakeRGB(OBColor color)
  {
    stringstream ss;
    ss << "\"rgb(" << (int)(255*color.red) << ',' << (int)(255*color.green)
       << ',' << (int)(255*color.blue) << ")\" ";
    return ss.str();
  }

  string SVGPainter::RGBcode(OBColor color)
  {
    stringstream ss;
    ss << std::hex << std::setfill('0') << std::setw(2) << (int)(255*color.red) << (int)(255*color.green)
       <<  (int)(255*color.blue);
    return ss.str();

  }
}
