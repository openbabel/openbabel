/**********************************************************************
Copyright (C) 2011 by Noel O'Boyle

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/op.h>
#include <openbabel/depict/depict.h>
#include <openbabel/depict/cairopainter.h>

using namespace std;



namespace OpenBabel
{
  /*
  // Class definition of CairoPainter
  CairoPainter::CairoPainter() : m_surface(0), m_cairo(0), m_fontPointSize(12)
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
    // clean up
    if (m_cairo)
      cairo_destroy(m_cairo);
    if (m_surface)
      cairo_surface_destroy(m_surface);
    // create new surface to paint on
    m_surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, static_cast<int> (width), static_cast<int> (height));
    m_cairo = cairo_create(m_surface);
    cairo_set_source_rgb (m_cairo, 255, 255, 255);
    cairo_paint (m_cairo);
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

  void CairoPainter::SetPenColor(const OBColor &color)
  {
    cairo_set_source_rgb(m_cairo, color.red, color.green, color.blue);
  }
      
  void CairoPainter::SetPenWidth(double width)
  {
    cairo_set_line_width(m_cairo, width);
  }

  void CairoPainter::DrawLine(double x1, double y1, double x2, double y2)
  {
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
    cairo_stroke(m_cairo);
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
    metrics.width = te.width;
    metrics.height = fe.height;
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
    for(int i=0;i<length;++i)
      in->push_back(data[i]);
    return CAIRO_STATUS_SUCCESS;
  }

  static cairo_surface_t *
  scale_surface (cairo_surface_t *old_surface, 
                int old_width, int old_height,
                int new_width, int new_height)
  {
    cairo_surface_t *new_surface = cairo_surface_create_similar(old_surface, CAIRO_CONTENT_COLOR_ALPHA, new_width, new_height);
    cairo_t *cr = cairo_create (new_surface);

    // Draw white background
    cairo_set_source_rgb (cr, 255, 255, 255);
    cairo_rectangle (cr, 0, 0, new_width, new_height);
    cairo_fill (cr);

    // Work out the scaling factor
    double scale_x = new_width / (double) old_width;
    double scale_y = new_height / (double) old_height;
    double scale = std::min(scale_x, scale_y);

    // Translate the over-scaled dimension into the centre
    if (scale < scale_y)
      cairo_translate(cr, 0, new_height/2.0 - scale*old_height/2.0);
    else
      cairo_translate(cr, new_width/2.0 - scale*new_width/2.0, 0);
    
    cairo_scale(cr, scale, scale); // Scale the drawing

    cairo_set_source_surface(cr, old_surface, 0, 0); // Redraw the old surface onto the new

    cairo_paint(cr); // Do the actual drawing

    cairo_destroy(cr);

    return new_surface;
  }

  void CairoPainter::WriteImage(std::ostream& ofs)
  {
    if (!m_cairo || !m_surface)
      return;
    vector<char> in;
    int width = cairo_image_surface_get_width(m_surface);
    int height = cairo_image_surface_get_height(m_surface);
    cairo_surface_t *new_surface = scale_surface (m_surface, width, height, 300, 300);
    cairo_surface_write_to_png_stream(new_surface, writeFunction, &in);
    for(int i=0; i<in.size(); ++i)
      ofs << in.at(i);
  }*/
  // End of class definition of CairoPainter



class PNG2Format : public OBMoleculeFormat
// Derive directly from OBFormat for objects which are not molecules.
{
public:
  //Register this format type ID in the constructor
  PNG2Format()
  {
    OBConversion::RegisterFormat("png2",this);

  }

  /* The first line of the description should be a brief identifier, <40 chars, because
     it is used in dropdown lists, etc. in some user interfaces. The rest is optional.

     Describe any format specific options here. This text is parsed to provide
     checkboxes, etc for the GUI (for details click the control menu),
     so please try to keep to a similar form.

     Write options are the most common, and the "Write" is optional.
     The option f takes a text parameter, so that it is essential that the option
     is registered in the constructor of the class.
     Finish the options with a blank line as shown, if there are more than one
     group of options, or if there are further comments after them.
  */
  virtual const char* Description() //required
  {
    return
    "XXX format\n"
    "Some comments here, on as many lines as necessay\n"
    "Write Options e.g. -xf3 \n"
    "  f# Number of (fictional) levels\n"
    "  n  Omit (virtual) title\n\n"
    
    "Read Options e.g. -as\n"
    "  s  Consider single bonds only\n"
    ;
  };

  //Optional URL where the file format is specified
  virtual const char* SpecificationURL(){return
     "http://www.mdl.com/downloads/public/ctfile/ctfile.jsp";};

  //Optional
  virtual const char* GetMIMEType()
  { return "chemical/x-xxx"; };


  /* Flags() can return be any of the following combined by |
     or be omitted if none apply
     NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY  DEFAULTFORMAT
     READBINARY  WRITEBINARY  READXML  ZEROATOMSOK*/
  virtual unsigned int Flags()
  {
      return NOTREADABLE | WRITEBINARY | WRITEONEONLY;
  };

   /* This optional function is for formats which can contain more than one
     molecule. It is used to quickly position the input stream after the nth
     molecule without have to convert and discard all the n molecules.
     See obconversion.cpp for details and mdlformat.cpp for an example.*/
  virtual int SkipObjects(int n, OBConversion* pConv)
  {
    return 0;
  };

  ////////////////////////////////////////////////////
  /// Declarations for the "API" interface functions. Definitions are below
  virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
  virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

private:
  /* Add declarations for any local function or member variables used.
     Generally only a single instance of a format class is used. Keep this in
     mind if you employ member variables. */
};
  ////////////////////////////////////////////////////

//Make an instance of the format class
PNG2Format thePNG2Format;

/////////////////////////////////////////////////////////////////

bool PNG2Format::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{
  OBMol* pmol = pOb->CastAndClear<OBMol>();
  if(pmol==NULL)
      return false;

  istream& ifs = *pConv->GetInStream();

  pmol->BeginModify();

  /** Parse the input stream and use the OpenBabel API to populate the OBMol **/

  // To use an input option
  if(pConv->IsOption("s",OBConversion::INOPTIONS))
  {
    //Code for when -as is specified
  }

  /* If the molecule has other than 3D coordinates for its atoms, it
  is necessary to set the dimension to 0, or 2 */
  int dim;
  dim = 3;
  pmol->SetDimension(dim);

  pmol->EndModify();

  /* For multi-molecule formats, leave the input stream at the start of the
     next molecule, ready for this routine to be called again.

  /* Return true if ok. Returning false means discard the OBMol and stop
     converting, unless the -e option is set. With a multi-molecule inputstream
     this will skip the current molecule and continue with the next, if SkipObjects()
     has been defined. If it has not, and continuation after errors is still required,
     it is necessary to leave the input stream at the beginning of next object when
     returning false;*/
  return true;
}

////////////////////////////////////////////////////////////////

bool PNG2Format::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
  OBMol* pmol = dynamic_cast<OBMol*>(pOb);
  if(pmol==NULL)
      return false;

  ostream& ofs = *pConv->GetOutStream();

  //*** Coordinate generation ***
  //Generate coordinates only if no existing 2D coordinates
  if(!pmol->Has2D(true))
  {
    OBOp* pOp = OBOp::FindType("gen2D");
    if(!pOp)
    {
      obErrorLog.ThrowError("PNG2Format", "gen2D not found", obError, onceOnly);
      return false;
    }
    if(!pOp->Do(pmol))
    {
      obErrorLog.ThrowError("PNG2Format", string(pmol->GetTitle()) + "- Coordinate generation unsuccessful", obError);
      return false;
    }
  }
  if(!pmol->Has2D() && pmol->NumAtoms()>1)
  {
    string mes("Molecule ");
    mes += pmol->GetTitle();
    mes += " needs 2D coordinates to display in PNG2format";
    obErrorLog.ThrowError("PNG2Format", mes, obError);
    return false;
  }

  CairoPainter painter;
  OBDepict depictor(&painter);
  depictor.DrawMolecule(pmol);
  painter.WriteImage(ofs);

  return true; //or false to stop converting
}

} //namespace OpenBabel

