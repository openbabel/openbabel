/**********************************************************************
asciipainter.cpp  - Render molecule as ASCII

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
#include <openbabel/depict/asciipainter.h>

#include <cstdlib> // for abs and fabs
#include <iostream>
using namespace std;

namespace OpenBabel
{

  ASCIIPainter::ASCIIPainter(int width, int height, double aspect): m_width(width), m_height(height), m_aspect(aspect), m_scale(1)
  {
    m_buf.clear();
    m_buf.reserve(m_height);
    for (int i=0; i<m_height; ++i) {
      vector<char> tmp(m_width, ' ');
      m_buf.push_back(tmp);
    }
  }

  ASCIIPainter::~ASCIIPainter()
  {
  }

  void ASCIIPainter::NewCanvas(double width, double height)
  {
    double sx = m_width/width;
    double sy = m_height*m_aspect/height;
    m_scale = std::min(sx ,sy);
  }

  bool ASCIIPainter::IsGood() const
  {
    return true;
  }

  void ASCIIPainter::SetFontSize(int pointSize)
  {
  }

  void ASCIIPainter::SetFillColor(const OBColor &color)
  {
  }

  void ASCIIPainter::SetFillRadial(const OBColor &start, const OBColor &end)
  {
  }


  void ASCIIPainter::SetPenColor(const OBColor &color)
  {
  }

  void ASCIIPainter::SetPenWidth(double width)
  {
  }

  double ASCIIPainter::GetPenWidth()
  {
    return 0.0;
  }

  void ASCIIPainter::DrawLine(double x1, double y1, double x2, double y2, const std::vector<double> & dashes)
  {
    vector<pair<int, int> > coords;
    vector<pair<int, int> >::iterator vp_it;
    string symbols = Bresenham(round(x1*m_scale), round(y1*m_scale/m_aspect),
                               round(x2*m_scale), round(y2*m_scale/m_aspect), coords);
    string::iterator s_it = symbols.begin();
    for(vp_it=coords.begin();vp_it!=coords.end();++vp_it,++s_it) {
      int x = vp_it->first;
      int y = vp_it->second;
      if (x>=0 && x<m_width && y>=0 && y<m_height)
        m_buf.at(y).at(x) = *s_it;
    }
  }

  void ASCIIPainter::DrawPolygon(const std::vector<std::pair<double,double> > &points)
  {
    if (points.size() < 2) return;

    vector<pair<double,double> >::const_iterator vp_it;
    for(vp_it=points.begin();vp_it!=(points.end() - 1);++vp_it) {
      DrawLine(vp_it->first, vp_it->second, (vp_it+1)->first, (vp_it+1)->second);
    }
    DrawLine(vp_it->first, vp_it->second, points.begin()->first, points.begin()->second);
  }

  void ASCIIPainter::DrawCircle(double x, double y, double r)
  {
  }

  void ASCIIPainter::DrawBall(double x, double y, double r, double opacity)
  {
  }

  void ASCIIPainter::DrawText(double x, double y, const std::string &text)
  {
    int sx = round(x*m_scale);
    int sy = round(y*m_scale/m_aspect);
    for (int i=0; i<text.size(); ++i) {
      if (sy>=0 && sy<m_height && sx+i>=0 && sx+i<m_width)
        m_buf.at(sy).at(sx+i) = text.at(i);
    }
  }

  OBFontMetrics ASCIIPainter::GetFontMetrics(const std::string &text)
  {
    OBFontMetrics metrics;
    // The following line tries to workaround the fact that fontSize is an int
    // -  the correct value is 1.0 / m_scale * m_aspect, and we round to the nearest integer.
    metrics.fontSize = round(1.0 / m_scale * m_aspect + 0.5);
    metrics.ascent = 0;
    metrics.descent = 0;
    metrics.width = 1.0 / m_scale;
    metrics.height = 1;
    return metrics;
  }
  void ASCIIPainter::Write(std::ostream &ofs)
  {
    vector<vector<char> >::iterator vvc_it;
    vector<char>::iterator vc_it;
    for(vvc_it=m_buf.begin();vvc_it!=m_buf.end();++vvc_it) {
      for(vc_it=vvc_it->begin();vc_it!=vvc_it->end();++vc_it) {
        ofs << *vc_it;
      }
      ofs << endl;
    }
  }

  int ASCIIPainter::round(double r)
  {
    return static_cast<int>( (r > 0.0) ? r + 0.5 : r - 0.5 );
  }

  // Helper function for Bresenham
  int getdelta(int x, int y, int x2, int y2)
  {
    int ans = 0;
    if (x2 == x) return ans; // Not applicable to vertical lines
    if (y2 > y) {
      ans = 1;
      double slope = (y2-y) / double(x2-x);
      if (fabs(slope) > 1.0)
        ans = 0;
    }
    return ans;
  }

  // Helper function for Bresenham
  string getsymbols(int x, int y, int x2, int y2)
  {
    if (x2 < x) {
      swap(x2, x);
      swap(y2, y);
    }
    double slope;
    if (x == x2)
      slope = 1e99;
    else
      slope = (y2-y) / double(x2-x);

    string ans;
    if (slope > 0)
      ans = (slope > 1.0) ? "|\\": "_\\";
    else
      ans = (slope >-1.0) ? "_/": "|/";
    return ans;
  }

  // (Translated from the Python prototype at https://gist.github.com/2331153)
  // Implementation of Bresenham line generation, a standard way of rasterising lines.
  // Returns a string with symbols for the line along with the coords.
  string ASCIIPainter::Bresenham(int x, int y, int x2, int y2, vector<pair<int, int> > &coords)
  {
    string symbols = getsymbols(x, y, x2, y2);
    int delta = getdelta(x, y, x2, y2);

    int dx = abs(x2 - x);
    int sx = (x2 - x > 0) ? 1: -1;
    int dy = abs(y2 - y);
    int sy = (y2 - y > 0) ? 1: -1;
    bool steep = false;
    if (dy > dx) {
      steep = true;
      swap(x, y);
      swap(dx, dy);
      swap(sx, sy);
    }

    string ans;
    int d = 2*dy - dx;
    for(int i=0; i<dx; ++i) {
      ans.append(d>=0 ? symbols.substr(1,1): symbols.substr(0,1));
      if (steep)
        coords.push_back(pair<int,int>(y, x));
      else {
        int tmp_delta = d>=0 ? delta : 0; // Correction for / and \ in some cases
        coords.push_back(pair<int,int>(x, y + tmp_delta));
      }
      while (d >= 0) {
        y += sy;
        d -= 2*dx;
      }
      x += sx;
      d += 2*dy;
    }
    ans.append(d>=0 ? symbols.substr(1,1): symbols.substr(0,1));
    int tmp_delta = d>=0 ? delta : 0; // Correction for / and \ in some cases
    coords.push_back(pair<int,int>(x2, y2 + tmp_delta));

    // Correction for top end of vertical bars
    if (coords.size() > 1) {
      if (coords.at(0).second < coords.at(1).second) {
        if (ans.at(0) == '|') {
          ans = ans.substr(1);
          coords.erase(coords.begin());
        }
      }
      else if (coords.at(coords.size()-1).second < coords.at(coords.size()-2).second) {
        if (ans.at(ans.size()-1) == '|') {
          ans = ans.substr(0, ans.size()-1);
          coords.erase(coords.end() - 1);
        }
      }
    }
    return ans;
  }
}
