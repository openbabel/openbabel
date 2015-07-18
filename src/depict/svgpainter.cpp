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
using namespace std;

namespace OpenBabel
{

  SVGPainter::SVGPainter(ostream& ofs, bool withViewBox,
    double width, double height, double x, double y)
    :  m_ofs(ofs), m_withViewBox(withViewBox), m_width(width), m_height(height),
       m_x(x), m_y(y), m_Pencolor("black"), m_Fillcolor("white"), m_Gradientcolor(make_pair(OBColor("white"),OBColor("white"))), m_PenWidth(1),
       m_fontPointSize(16), m_isFillcolor(true)  {}

  SVGPainter::~SVGPainter()
  {
    WriteDefs();
    m_ofs << m_mems.str();
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
      m_mems << "<rect x=\"0%\" y=\"0%\" width=\"100%\" height=\"100%\" stroke-width=\"0\" fill="
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
    m_isFillcolor = true;
  }

  void SVGPainter::SetFillRadial(const OBColor &start, const OBColor &end)
  {
    m_Gradientcolor = make_pair(start,end);
    m_Gradients.insert(m_Gradientcolor);
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

  void SVGPainter::DrawLine(double x1, double y1, double x2, double y2)
  {
    streamsize oldprec = m_mems.precision(1);
    m_mems << fixed << "<line x1=\"" << x1 << "\" y1=\"" << y1 << "\" x2=\""
      << x2 << "\" y2=\"" << y2 << "\"";
    // if(m_Pencolor!=m_OrigBondcolor) // TODO: Bring this line back once Pybel is fine with this
      m_mems << " stroke=" << MakeRGB(m_Pencolor);
    m_mems << " stroke-width=\"" << m_PenWidth << "\"";
    m_mems << "/>\n";
    m_mems.precision(oldprec);
  }

  void SVGPainter::DrawPolygon(const std::vector<std::pair<double,double> > &points)
  {
    m_mems << "<polygon points=\"";
      std::vector<std::pair<double,double> >::const_iterator i;
    for (i = points.begin(); i != points.end(); ++i)
      m_mems << i->first << ' ' << i->second << ' ';
    m_mems << "\"";
    m_mems << " stroke-width=\"" << m_PenWidth << "\"";
    m_mems << " fill=" << MakeRGB(m_Pencolor);
    m_mems << " stroke=" << MakeRGB(m_Pencolor);
    m_mems << "/>\n";
  }

  void SVGPainter::DrawCircle(double x, double y, double r)
  {
    m_mems << "<circle cx=\"" << x << "\" y=\"" << y << "\" r=\"" << r << "\" />\n";
  }

  void SVGPainter::DrawText(double x, double y, const std::string &text)
  {
    m_mems << "<text x=\"" << x << "\" y=\"" << y << "\""
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

  void SVGPainter::DrawBall(double x, double y, double r)
  {
    m_mems << "<circle cx=\"" << x << "\" cy=\"" << y << "\" r=\"" << r << "\" ";
    if (m_isFillcolor) {
      m_mems << "style=\"stroke:black;fill:" << MakeRGB(m_Fillcolor) << "\"/>\n";
    } else
    {
      m_mems << "style=\"stroke:black;fill:url(#radial";
      m_mems << RGBcode(m_Gradientcolor.first)<< RGBcode(m_Gradientcolor.second) << ")\"/>\n";
    }
  }

  void   SVGPainter::WriteDefs()
  {
    if (!m_Gradients.empty()) {
      m_ofs << "<defs>\n";
      for (std::set<ColorGradient>::iterator it=m_Gradients.begin(); it!=m_Gradients.end(); ++it) {
        m_ofs << "<radialGradient id='radial";
        m_ofs << RGBcode(it->first)<< RGBcode(it->second) << "' ";
        m_ofs << "cx='50%' cy='50%' r='50%' fx='30%' fy='30%'>\n";
        m_ofs << "  <stop offset=' 0%' start-color=" << MakeRGB(it->first) << " stop-opacity='0.8'/>\n";
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

