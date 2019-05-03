/**********************************************************************
obfit = Fit molecules according to a SMART pattern
Copyright (C) 2003 Fabien Fontaine
Some portions Copyright (C) 2004-2005 Geoffrey R. Hutchison

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

/*
  Require a fixed molecule, a set of molecules to move, 
  a SMART pattern to match the fixed and the moving molecules
  If the SMART is not found in the fixed molecule the program exits
  If the SMART is not found in a moving molecule, the molecule is not moved 
  example of command line:
  obfit "[nh]1c2c(=O)n(C)c(=O)n(C)c2cc1" testref.sdf testmv.sdf 
*/


// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/parsmart.h>
#include <openbabel/obconversion.h>
#include <openbabel/math/vector3.h>
#include <openbabel/generic.h>
#include <openbabel/atom.h>
#ifndef _MSC_VER
  #include <unistd.h>
#endif
#include <cstdlib>

using namespace std;
using namespace OpenBabel;


//find the center of mass of a list of atoms
vector3 mass_c( vector<int> &aindex, OBMol &mol);

///////////////////////////////////////////////////////////////////////////////
//! \brief superimpose a set of molecules on the atoms of a reference molecule
//! The atoms used for the overlay are defined by the SMART pattern
int main(int argc,char **argv)
{

  int errflg=0;
  char *FileRef=NULL, *FileMove=NULL, *Pattern=NULL;
  string err;
  char *program_name=argv[0];

  // parse the command line
  if (argc!=4)
    {
      errflg++;
    }
  else
    {
      FileRef = argv[2];
      FileMove = argv[3];
      Pattern = argv[1];
    }

  if (errflg)
    {
      err = "Usage: ";
      err += program_name;
      err += " \"PATTERN\" <fixed_structure> <moving_structures>\n";
      ThrowError(err);
      exit(-1);
    }


  // create the pattern
  OBSmartsPattern sp;
  if (!sp.Init(Pattern))
    {
      err = program_name;
      err += ": Unable to read the SMART: ";
      err += Pattern;
      ThrowError(err);
      exit(-1);
    }

  // Find Input filetypes
  OBConversion conv;
  OBFormat *refFormat = conv.FormatFromExt(FileRef);
  if (!refFormat || !conv.SetInAndOutFormats(refFormat,refFormat))
    {
      cerr << program_name << ": cannot read fixed molecule format!" << endl;
      exit (-1);
    }

  OBFormat *moveFormat = conv.FormatFromExt(FileMove);
  if (!moveFormat || !conv.SetInAndOutFormats(moveFormat,moveFormat))
    {
      cerr << program_name << ": cannot read moving molecule(s) format!" << endl;
      exit (-1);
    }

  conv.SetInAndOutFormats(refFormat,moveFormat);

  ifstream ifsref;
  OBMol molref;
  vector< vector <int> > maplist;      // list of matched atoms
  vector< vector <int> >::iterator i;  // and its iterators
  vector< int >::iterator j;
  vector <int> refatoms;

  //Read the reference structure
  ifsref.open(FileRef);
  if (!ifsref)
    {
      cerr << program_name << ": cannot read fixed molecule file: "
           << FileRef << endl;
      exit (-1);
    }

  molref.Clear();
  conv.Read(&molref,&ifsref);

  // and check if the SMART match
  sp.Match(molref);
  maplist = sp.GetUMapList(); // get unique matches
  if (maplist.empty())
    {
      err = program_name;
      err += ": Unable to map SMART: ";
      err += Pattern;
      err += " in reference molecule: ";
      err += FileRef;
      ThrowError(err);
      exit(-1);
    }

  // Find the matching atoms
  for (i = maplist.begin(); i != maplist.end(); ++i)
    {
      refatoms.clear(); // Save only the last set of atoms
      for(j= (*i).begin(); j != (*i).end(); ++j)
        {
          refatoms.push_back(*j);
        }
    }
  if(refatoms.size() > molref.NumAtoms())
    {  
      cerr << program_name 
           << ": The SMARTS pattern produces more matching atoms than are in the reference molecule"
           << endl;
      exit(-1);
    }

  // set the translation vector
  vector3 tvref(0,0,0);
  OBAtom *atom;
  unsigned int c;

  tvref = mass_c(refatoms, molref);
  // center the molecule
  molref.Translate(-tvref);

  //get the coordinates of the SMART atoms
  double *refcoor = new double[refatoms.size()*3];

  for(c=0; c<refatoms.size(); ++c)
    {
      atom = molref.GetAtom(refatoms[c]);
      refcoor[c*3] = atom->x();
      refcoor[c*3+1] = atom->y();
      refcoor[c*3+2] = atom->z();
    }

  conv.SetInAndOutFormats(moveFormat,moveFormat);

  ifstream ifsmv;
  OBMol molmv;
  vector <int> mvatoms;
  vector3 tvmv;
  unsigned int size=0;
  double rmatrix[3][3];

  //Read the moving structures
  ifsmv.open(FileMove);
  if (!ifsmv)
    {
      cerr << program_name << ": cannot read file: "
           << FileMove << endl;
      exit (-1);
    }

  for (;;)
    {
      molmv.Clear();
      conv.Read(&molmv,&ifsmv);                   // Read molecule
      if (molmv.Empty())
        break;

      if (sp.Match(molmv))          // if match perform rotation
        {

          maplist = sp.GetMapList(); // get all matches

          // Find the matching atoms
	    
          // Looping over all matches to find best match
          double rmsd;
          double best_rmsd = 999.999;
          vector <int> best_mvatoms;

          for (i = maplist.begin(); i != maplist.end(); ++i)
            {
              mvatoms.clear(); // Save only the last set of atoms
              for(j= (*i).begin(); j != (*i).end(); ++j)
                {
                  mvatoms.push_back(*j);
                }

              tvmv = mass_c(mvatoms, molmv);
              // center the molecule
              molmv.Translate(-tvmv);

              //Find the rotation matrix
              size = mvatoms.size();
              if (size != refatoms.size())
                {
                  err = program_name;
                  err += ": Error: not the same number of SMART atoms";
                  ThrowError(err);
                  exit(-1);
                }

              double *mvcoor = new double[size*3];

              for(c=0; c < size; ++c)
                {
                  atom = molmv.GetAtom(mvatoms[c]);
                  mvcoor[c*3] = atom->x();
                  mvcoor[c*3+1] = atom->y();
                  mvcoor[c*3+2] = atom->z();
                }

              // quaternion fit
              qtrfit(refcoor, mvcoor, size, rmatrix);

              //rotate all the atoms
              molmv.Rotate(rmatrix);

              // update mvcoor after rotation
              for(c=0; c < size; ++c)
                {
                  atom = molmv.GetAtom(mvatoms[c]);
                  mvcoor[c*3] = atom->x();
                  mvcoor[c*3+1] = atom->y();
                  mvcoor[c*3+2] = atom->z();
                }

              rmsd = calc_rms(refcoor,mvcoor,size);
              if ( rmsd < best_rmsd )
                {
                  best_rmsd = rmsd;
                  best_mvatoms.clear();
                  best_mvatoms.resize(mvatoms.size());
                  best_mvatoms = mvatoms;
                }
                 
              delete[] mvcoor;
            } // loop through matches
 
          // Refit molecule using best match
          mvatoms.clear();
          mvatoms.resize(best_mvatoms.size());
          mvatoms = best_mvatoms;

          tvmv = mass_c(mvatoms, molmv);
          // center the molecule
          molmv.Translate(-tvmv);

          //Find the rotation matrix
          size = mvatoms.size();
          if (size != refatoms.size())
            {
              err = program_name;
              err += ": Error: not the same number of SMART atoms";
              ThrowError(err);
              exit(-1);
            }

          double *mvcoor = new double[size*3];
          for(c=0; c<size; ++c)
            {
              atom = molmv.GetAtom(mvatoms[c]);
              mvcoor[c*3] = atom->x();
              mvcoor[c*3+1] = atom->y();
              mvcoor[c*3+2] = atom->z();
            }

          // quaternion fit
          qtrfit(refcoor, mvcoor, size, rmatrix);

          //rotate all the atoms
          molmv.Rotate(rmatrix);

          for(c=0; c<size; ++c)
            {
              atom = molmv.GetAtom(mvatoms[c]);
              mvcoor[c*3] = atom->x();
              mvcoor[c*3+1] = atom->y();
              mvcoor[c*3+2] = atom->z();
            }
	    
          rmsd = calc_rms(refcoor,mvcoor,size);

          char rmsd_string[80];
          sprintf(rmsd_string,"%f", best_rmsd);
  
          OBPairData *dp = new OBPairData;
          string field_name = "RMSD";
	    
          dp->SetAttribute(field_name);
          dp->SetValue(rmsd_string);
          dp->SetOrigin(external);
          molmv.SetData(dp);

          cerr << "RMSD: " << rmsd_string << endl;

          //translate the rotated molecule
          molmv.Translate(tvref);

          delete[] mvcoor;
        }
      conv.Write(&molmv,&cout);
    }

  delete[] refcoor;
  return(0);
}

///////////////////////////////////////////////////////////////////////////////
//find the center of mass of a list of atoms
vector3 mass_c( vector<int> &aindex, OBMol &mol)
{
    vector3 center(0,0,0);
    vector< int >::iterator j;
    OBAtom *atom;

    for(j= aindex.begin(); j != aindex.end(); ++j)
    {
        atom = mol.GetAtom(*j);
        center += atom->GetVector();
    }

    center /= (float) aindex.size();
    return (center);
}

/* obfit man page*/
/** \page obfit superimpose two molecules based on a pattern
*
* \n
* \par SYNOPSIS
*
* \b obfit <SMARTS-pattern> \<fixed-file\> \<outfile\>
*
* \par DESCRIPTION
*
* Superimpose two molecules using a quaternion fit. The atoms used to fit the
* two molecules are defined by the SMARTS pattern given by the user. It is
* useful to align congeneric series of molecules on a common structural
* scaffold for 3D-QSAR studies. It can also be useful for displaying the
* results of conformational generation.
* \n\n
* Any molecules matching the supplied SMARTS pattern will be rotated and
* translated to provide the smallest possible RMSD between the matching
* regions. If a molecule does not match the SMARTS pattern, it will be output
* with no transformation.
*
* \par EXAMPLES
*  - Align all the molecules in 'moving.sdf' on a single molecule of 'fixed.sdf'
*    by superimposing them on its N-methylpiperidyl portion \n
*       obfit "[nh]1c2c(=O)n(C)c(=O)n(C)c2cc1" testref.sdf testmv.sdf
*
* \par AUTHORS
*
* The obfit program was contributed by \b Fabien \b Fontaine.
*
* Open Babel is currently maintained by \b Geoff \b Hutchison, \b Chris \b Morley and \b Michael \b Banck.
*
* For more contributors to Open Babel, see http://openbabel.org/THANKS.shtml
*
* \par COPYRIGHT
*  Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
*  Some portions Copyright (C) 2001-2005 by Geoffrey R. Hutchison \n \n
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation version 2 of the License.\n \n
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
* \par SEE ALSO
*   The web pages for Open Babel can be found at http://openbabel.org/ \n
*   A guide for constructing SMARTS patterns can be found at http://www.daylight.com/dayhtml/doc/theory/theory.smarts.html
**/
