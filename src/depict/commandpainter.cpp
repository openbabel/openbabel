/**********************************************************************
commandpainter.cpp  - Render as depiction commands

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
#include <openbabel/obutil.h>
#include <openbabel/depict/commandpainter.h>

#include <iostream>
using namespace std;

namespace OpenBabel
{

  // Class definition of CommandPainter
  CommandPainter::CommandPainter(ostream& ofs): m_ofs(ofs)
  {
    m_oldprec = m_ofs.precision(1);
    m_ofs << fixed;
  }

  CommandPainter::~CommandPainter()
  {
    m_ofs.precision(m_oldprec);
    m_ofs.unsetf(ios_base::fixed);
  }

  void CommandPainter::NewCanvas(double width, double height)
  {
    m_ofs << "NewCanvas " << width << " " << height << endl;
  }

  bool CommandPainter::IsGood() const
  {
    return true;
  }

  void CommandPainter::SetFontSize(int pointSize)
  {
    m_ofs << "SetFontSize " << pointSize << endl;
  }

  void CommandPainter::SetFillColor(const OBColor &color)
  {
    m_ofs << "SetFillColor " << color.red << " " << color.green << " " << color.blue << " " << color.alpha << " (rgba)" << endl;
  }

  void CommandPainter::SetFillRadial(const OBColor &start, const OBColor &end)
  {
    m_ofs << "SetFillRadial" << start.red << " " << start.green << " " << start.blue << " " << start.alpha << " (rgba) to ";
    m_ofs << end.red << " " << end.green << " " << end.blue << " " << end.alpha << " (rgba)" << endl;
  }


  void CommandPainter::SetPenColor(const OBColor &color)
  {
    m_ofs << "SetPenColor " << color.red << " " << color.green << " " << color.blue << " " << color.alpha << " (rgba)" << endl;
  }

  void CommandPainter::SetPenWidth(double width)
  {
    m_pen_width = width;
    m_ofs << "SetPenWidth " << width << endl;
  }

  double CommandPainter::GetPenWidth()
  {
    return m_pen_width;
  }

  void CommandPainter::DrawLine(double x1, double y1, double x2, double y2,  const std::vector<double> & dashes)
  {

    m_ofs << fixed << "DrawLine " << x1 << " " << y1 << " to "
                                  << x2 << " " << y2;
    if (!dashes.empty()) {
      std::vector<double>::const_iterator it;
      m_ofs << " dashes";
      for (it=dashes.begin(); it!=dashes.end() ; ++it)
        m_ofs << " " << *it;

    }
    m_ofs << endl;
  }

  void CommandPainter::DrawPolygon(const std::vector<std::pair<double,double> > &points)
  {
    m_ofs << "DrawPolygon ";
    std::vector<std::pair<double,double> >::const_iterator i;
    for (i = points.begin(); i != points.end(); ++i) {
      if (i != points.begin())
        m_ofs << " to ";
      m_ofs << i->first << " " << i->second;
    }
    m_ofs << endl;
  }

  void CommandPainter::DrawCircle(double x, double y, double r)
  {
    m_ofs << "DrawCircle " << x << " " << y << " radius " << r << endl;
  }

  void CommandPainter::DrawBall(double x, double y, double r, double opacity)
  {
    m_ofs << "DrawBall " << x << " " << y << " radius " << r << endl;
  }

  void CommandPainter::DrawText(double x, double y, const std::string &text)
  {
    m_ofs << "DrawText " << x << " " << y << " \"" << text << "\"" << endl;
  }

  OBFontMetrics CommandPainter::GetFontMetrics(const std::string &text)
  {
    OBFontMetrics metrics;
    metrics.fontSize = 0;
    metrics.ascent = 0;
    metrics.descent = 0;
    metrics.width = 0;
    metrics.height = 0;
    return metrics;
  }
}
