/**********************************************************************
Copyright (C) 2004 by Chris Morley

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

/* This is a heavily commented template for a OpenBabel format class.

Format classes are plugins: no modifications are needed to existing code
to indroduce a new format. The code just needs to be compiled and linked
with the rest of the OpenBabel code.
Alternatively, they can be built (either singly or in groups) as DLLs
or shared libraries. [**Extra build info**]

Each file may contain more than one format.

This compilable, but non-functional example is for a format which
converts a molecule to and from OpenBabel's internal format OBMol.
The conversion framework can handle other types of object, provided
they are derived from OBBase, such as OBReaction, OBText.

For XML formats, extra support for the parsing is provided, see pubchem.cpp
as an example.
*/

#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>

using namespace std;
namespace OpenBabel
{

class XXXFormat : public OBMoleculeFormat
// Derive directly from OBFormat for objects which are not molecules.
{
public:
	//Register this format type ID in the constructor
  XXXFormat()
	{
		/* XXX is the file extension and is case insensitive. A MIME type can be
		   added as an optional third parameter.
		   Multiple file extensions can be registered by adding extra statements.*/
		OBConversion::RegisterFormat("XXX",this);

		/* If there are any format specific options they should be registered here
		   so that the commandline interface works properly.
		   The first parameter is the option name. If it is a single letter it can be
		   concatinated with other	single letter options. For output options it can be
		   multicharcter and is then	written as --optionname on the command line.
		   The third parameter is the number of parameters the option takes. Currently
		   this is either 1 or 0 and if it is 0 can be omitted for output options.
			 The parameter is always text and needs to be parsed to extract a number.

		   Options can apply when writing - 4th parameter is OBConversion::OUTOPTIONS
		   or can be omitted as shown. A single letter output option is preceded
		   by -x on the command line.
		   Or options can apply to the input format - the 4th parameter is
		   OBConversion::INOPTIONS. They are then  preceded by -a on the command line.

		   Each option letter may be reused in other formats, but within the same group,
		   INOPTIONS or OUTOPTIONS, must take the same number of parameters (0 or 1).
		   There will be an error message when OpenBabel	runs if there are conflicts
		   between formats. A list of formats currently used (which may not be
		   comprehensive) is in docs/options.html.
		*/
		OBConversion::RegisterOptionParam("f", this, 1);
		OBConversion::RegisterOptionParam("n", this);
		OBConversion::RegisterOptionParam("s", this, 0, OBConversion::INOPTIONS);

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
		"	f# Number of (fictional) levels\n"
		"	n  Omit (virtual) title\n\n"
		
		"Read Options e.g. -as\n"
		"	s  Consider single bonds only\n"
		;
  };

  //Optional URL where the file format is specified
	virtual const char* SpecificationURL(){ return ""; }

  //Optional
	virtual const char* GetMIMEType()
  { return "chemical/x-xxx"; };


  /* Flags() can return be any of the following combined by |
	   or be omitted if none apply
     NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY  DEFAULTFORMAT
     READBINARY  WRITEBINARY  READXML  ZEROATOMSOK*/
  virtual unsigned int Flags()
  {
      return READONEONLY;
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
XXXFormat theXXXFormat;

/////////////////////////////////////////////////////////////////

bool XXXFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
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

bool XXXFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
  OBMol* pmol = dynamic_cast<OBMol*>(pOb);
  if(pmol==NULL)
      return false;

  ostream& ofs = *pConv->GetOutStream();

	/** Write the representation of the OBMol molecule to the output stream **/

	//To use an output option
	if(!pConv->IsOption("n")) //OBConversion::OUTOPTIONS is the default
		ofs << "Title = " << pmol->GetTitle() << endl;
	//or if the option has a parameter
	int levels=0;
	const char* p = pConv->IsOption("f");
	if(p) //p==NULL if f option absent
		levels = atoi(p);

	// To find out whether this is the first molecule to be output...
	if(pConv->GetOutputIndex()==1)
		ofs << "The contents of this file were derived from " << pConv->GetInFilename() << endl;
	// ... or the last
	if(!pConv->IsLast())
		ofs << "$$$$" << endl;

	return true; //or false to stop converting
}

} //namespace OpenBabel

