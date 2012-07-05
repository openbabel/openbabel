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
#include <sstream>
using namespace std;

namespace OpenBabel
{

  SVGPainter::SVGPainter(ostream& ofs, bool withViewBox,
    double width, double height, double x, double y)
    :  m_ofs(ofs), m_withViewBox(withViewBox), m_width(width), m_height(height),
       m_x(x), m_y(y), m_Pencolor("black"), m_Fillcolor("white"), m_PenWidth(1),
       m_fontPointSize(16)  {}

  SVGPainter::~SVGPainter()
  {
    m_ofs << "</svg>\n";
    if(m_withViewBox)
      m_ofs << "</g>\n";
  }

  void SVGPainter::NewCanvas(double width, double height)
  {
    //Using withViewBox to supress xml header and xmlns attributes. May need another way.
    if(!m_withViewBox)
      m_ofs << "<?xml version=\"1.0\"?>\n";

    if(m_withViewBox)
      m_ofs << "<g transform=\"translate(" << m_x << "," << m_y << ")\">\n";

    m_ofs << "<svg ";
    if(!m_withViewBox)
      m_ofs << "xmlns=\"http://www.w3.org/2000/svg\"\n"
               "xmlns:xlink=\"http://www.w3.org/1999/xlink\" "
               "xmlns:cml=\"http://www.xml-cml.org/schema\" ";
    if(m_withViewBox)
      m_ofs << "width=\"" << m_width << "\" height=\"" << m_height << "\" "
            << "x=\"0\" y=\"0\" "
            << "viewBox=\"0 0 " << width << ' ' << height << "\"\n";
    else
      m_ofs << "width=\"" << width << "\" height=\"" << height << "\" "
            << "x=\"0\" y=\"0\" ";

    //Bond color and width are the initial m_Pencolor and m_PenWidth
    m_ofs << "font-family=\"" << m_fontFamily << "\" stroke=" << MakeRGB(m_Pencolor)
          << "stroke-width=\"" << m_PenWidth << "\"  stroke-linecap=\"round\"" << ">\n";

    if(!m_withViewBox && m_Fillcolor.alpha!=0.0)//Background color for single molecule. Handled by outer svg when table.
      m_ofs << "<rect x=\"0%\" y=\"0%\" width=\"100%\" height=\"100%\" stroke-width=\"0\" fill="
            << MakeRGB(m_Fillcolor) << " />\n";
    m_OrigBondcolor = m_Pencolor;
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
  }

  void SVGPainter::SetPenColor(const OBColor &color)
  {
    m_Pencolor = color; //value when NewCanvas called used for bonds
  }

  void SVGPainter::SetPenWidth(double width)
  {
    m_PenWidth = width; //value when NewCanvas called used for bonds
  }

  void SVGPainter::DrawLine(double x1, double y1, double x2, double y2)
  {
    streamsize oldprec = m_ofs.precision(1);
    m_ofs << fixed << "<line x1=\"" << x1 << "\" y1=\"" << y1 << "\" x2=\""
      << x2 << "\" y2=\"" << y2 << "\"";
    // if(m_Pencolor!=m_OrigBondcolor) // TODO: Bring this line back once Pybel is fine with this
      m_ofs << " stroke=" << MakeRGB(m_Pencolor);
    m_ofs << " stroke-width=\"" << m_PenWidth << "\"";
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
      << " fill=" << MakeRGB(m_Pencolor) << " stroke=" << MakeRGB(m_Pencolor) << "stroke-width=\"1\" "
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
    metrics.width = 0.0;
    for(string::size_type i=0;i<text.size();++i)
      metrics.width += m_fontPointSize * (isalpha(text[i]) ? 0.75 : 0.5);

    return metrics;
  }

  void SVGPainter::WriteImage(const std::string &filename)
  {
  }

  string SVGPainter::MakeRGB(OBColor color)
  {
    stringstream ss;
    ss << "\"rgb(" << (int)(255*color.red) << ',' << (int)(255*color.green)
       << ',' << (int)(255*color.blue) << ")\" ";
    return ss.str();
  }
}

