/**********************************************************************
  (C) 2008-2010 by Jens Thomas

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
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/elements.h>
#include <openbabel/generic.h>
#include <openbabel/internalcoord.h>


#include <algorithm>
#include <cmath>

#ifdef _MSC_VER
#include <regex>
#else
#include <regex.h>
#endif

using namespace std;

namespace OpenBabel
{
#define BOHR_TO_ANGSTROM 0.529177249
#define ANGSTROM_TO_BOHR 1.889725989

  class GAMESSUKFormat
  {
    /* Base class for GAMESS-UK readers with various utility functions that are used by
     * both input and output readers.
     *
     * The most important is ReadGeometry, which takes a list of strings defining the
     * geometry (in the style of a GAMESS-UK input deck) and creates the OBMol from it.
     * This routine supports both Zmatrix and Cartesian formats, although currently not
     * mixed decks.
     */

  public:
    bool ReadGeometry(OBMol &mol, vector<string> &geomList);
    bool ReadVariables(istream &ifs, double factor, string stopstr);
    bool ReadLineCartesian(OBAtom *atom, vector<string> &tokens, double factor);
    bool ReadLineZmatrix(OBMol &mol, OBAtom *atom, vector<string> &tokens, double factor, int *zmatLineCount);
    double Rescale(string text);
    bool IsUnits(string text);
    /**
     * Converts a string to a numerical type
     * This purloined from: http://www.codeguru.com/forum/showthread.php?t=231054
     */
    template <class T>
    bool from_string(T& t, const std::string& s,
                     std::ios_base& (*f)(std::ios_base&))
    {
      std::istringstream iss(s);
      return !(iss >> f >> t).fail();
    }

    // Variables
    enum ReadMode_t {CARTESIAN, ZMATRIX, VARIABLES, CONSTANTS, SKIP};
    ReadMode_t ReadMode;
    char buffer[BUFF_SIZE];
    stringstream errorMsg;

  private:
    map<string, double> variables; // map from variable name to value
    vector<OBInternalCoord*> vic; // Holds lists of internal coordinates
    int LabelToAtomicNumber(string label);
  };


  bool GAMESSUKFormat::ReadGeometry(OBMol &mol, vector<string> &geomList)
  {

    /* Read a geometry from a list. Any variables that appear in the geometry need
     * to be in the variables map that should have been populated before this is called.
     */

    if (geomList.size()==0){
      obErrorLog.ThrowError(__FUNCTION__,
                            "Problems reading a GAMESS-UK Input file: ReadGeometry got empty list",
                            obWarning);
      return false;
    }

    vector<string> tokens; // list of lines and list of tokens on a line
    string line; // For convenience so we can refer to lines from the iterator as 'line'
    double factor=BOHR_TO_ANGSTROM; // The coordinate conversion factor for handling bohr/angstrom issues

    mol.BeginModify();
    // Clear out any existing information
    mol.Clear();
    vic.clear();

    ReadMode=SKIP;
    bool ContainsZmatrix=false;
    int zmatLineCount=0;

    /*
      cerr << "ReadGeometry got geometry list: \n";
      for (vector<string>::iterator i=geomList.begin(); i !=geomList.end(); i++) {

      // Alias the line
      line = *i;
      cerr << "line: " << line << endl;
      }
    */

    for (vector<string>::iterator i=geomList.begin(); i !=geomList.end(); ++i) {

      // Alias the line
      line = *i;

      //cerr << "ReadGeometry line is: " << line << endl;

      // Check for commas & split with that as the separator if necessary
      if (line.find(',')!=string::npos) {
        tokenize(tokens, line, ",");
      } else {
        tokenize(tokens, line, " \t\n");
      }


      // Set the mode
      if (line.compare(0, 4, "zmat")==0 || line.compare(0, 4, "inte")==0) {
        ReadMode=ZMATRIX;
        //cout << "ZMATRIX mode " << ReadMode << endl;
        //cout << "tokens.size()" << tokens.size() << endl;
        if (tokens.size()>1) if (IsUnits(tokens[1])) factor=Rescale(tokens[1]);
        ContainsZmatrix=true;
        vic.push_back((OBInternalCoord*)NULL); // OBMol indexed from 1 -- potential atom index problem
      } else if (line.compare(0, 4, "coor")==0 || line.compare(0, 4, "cart")==0 ||line.compare(0, 4, "geom")==0) {
        ReadMode=CARTESIAN;
        //cout << "CARTESIAN mode " << ReadMode << endl;
        if (tokens.size()>1) if (IsUnits(tokens[1])) factor=Rescale(tokens[1]);

        /*
          We need to have read the variables first
          } else if (line.compare(0, 4, "vari")==0) {
          ReadMode=VARIABLES;
          //cout << "VARIABLES mode "<< ReadMode << endl;
          if (tokens.size() == 2) factor=Rescale(tokens[1]);
          //cout << "Factor now " << factor << endl;
          } else if (line.compare(0, 4, "cons")==0) {
          ReadMode=CONSTANTS;
          //cout << "CONSTANTS mode\n";
          if (tokens.size() == 2)
          factor=Rescale(tokens[1]);
          //cout << "Factor now " << factor << endl;
          */

      } else if (line.compare(0, 3, "end")==0) {
        ReadMode=SKIP;
        //cout << "SKIP mode " << ReadMode << endl;
      } else {
        if (ReadMode==SKIP) continue;
        if (ReadMode==ZMATRIX) {
          // Create an atom
          OBAtom *atom = mol.NewAtom();
          // Read the ZMatrix definition line
          if (! ReadLineZmatrix(mol,atom,tokens,factor,&zmatLineCount) )
            {
              errorMsg << "Problems reading a GAMESS-UK Input file: "
                       << "Could not read zmat line: " << line;
              obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() ,
                                    obWarning);
              return (false);
            }

        } // End ReadMode ZMATRIX

        if (ReadMode==CARTESIAN) {
          OBAtom *atom = mol.NewAtom();
          if (! ReadLineCartesian(atom,tokens,factor) )
            {
              errorMsg << "Problems reading a GAMESS-UK Input file: "
                       << "Could not read xyz line: " << line;
              obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() ,
                                    obWarning);
              return (false);
            }

        } // End ReadMode CARTESIAN


      } // End Test for first chars on line
    } // End loop over lines


    if (ContainsZmatrix)InternalToCartesian(vic,mol);
    mol.EndModify();

    return true;
  } // End Read Geometry

  bool GAMESSUKFormat::IsUnits(string text)
  {
    /* See if the supplied string specifies a unit */

    if ( text.compare(0, 4, "angs")==0 ||
         text.compare(0, 4, "bohr")==0 ||
         text.compare(0, 4, "a.u.")==0 ||
         text.compare(0, 2, "au")==0) {
      return true;
    } else {
      return false;
    }
  }

  double GAMESSUKFormat::Rescale(string text)
  {
    /* Return the correct scale factor given a string identifying the units */

    if (! IsUnits(text) ){
      errorMsg << "Problems reading GUK input - bad scale factor: " << text;
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      return -1.0;
    }

    if (text.compare(0, 4, "angs")==0) {
      return 1.0;
    } else if (text.compare(0, 4, "bohr")==0||text.compare(0, 4, "a.u.")==0
               ||text.compare(0, 2, "au")==0) {
      return BOHR_TO_ANGSTROM;
    } else {
      return -1.0;
    }
  }

  int GAMESSUKFormat::LabelToAtomicNumber(string label)
  {
    /*
     * Given a string with the label for an atom return the atomic number
     * As we are using the GetAtomicNum function case is not important
     */

    // See if the first 2 characters give us a valid atomic #
    int Z=OBElements::GetAtomicNum(label.substr(0,2).c_str());

    // If not try the first one
    if (Z==0) Z=OBElements::GetAtomicNum(label.substr(0,1).c_str());

    if (Z==0){
      // Check if it's an x (dummy) atom
      if(  label.substr(0,1) != "x" && label.substr(0,1) != "X" )
        {
          // Houston...
          errorMsg << "LabelToAtomicNumber got bad Label: " << label << std::endl;
          obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
        }
    }
    return Z;
  }

  bool GAMESSUKFormat::ReadLineCartesian(OBAtom *atom, vector<string> &tokens, double factor)
  {

    /*  Read a line defining the Cartesian coordinates for an atom
     * This assumes the line is formatted in GAMESS-UK input style as:
     * x y z AtomicNumber Label
     */
    
    bool ok=false;
    int Z;
    double x,y,z;

    // 4th field is the atomic number
    ok = from_string<int>(Z, tokens.at(3), std::dec);
    atom->SetAtomicNum(Z);

    // Read the atom coordinates
    ok = from_string<double>(x, tokens.at(0), std::dec);
    if ( ! ok)
      {
        // Can't convert to double so see if it's in the variables
        if (variables.find(tokens[0])==variables.end()) return false;
        x = variables[tokens[0]];
      }

    ok = from_string<double>(y, tokens.at(1), std::dec);
    if ( ! ok)
      {
        // Can't convert to double so see if it's in the variables
        if (variables.find(tokens[1])==variables.end()) return false;
        y = variables[tokens[1]];
      }

    ok = from_string<double>(z, tokens.at(2), std::dec);
    if ( ! ok)
      {
        // Can't convert to double so see if it's in the variables
        if (variables.find(tokens[2])==variables.end()) return false;
        z = variables[tokens[2]];
      }

    // Convert to Angstroms
    x=x*factor;
    y=y*factor;
    z=z*factor;
    atom->SetVector(x, y, z); //set coordinates
    return true;
  }

  bool GAMESSUKFormat::ReadLineZmatrix(OBMol &mol, OBAtom *atom, vector<string> &tokens, double factor, int *zmatLineCount)
  {
    /*
     * Read a line from a GAMESS-UK input defining an atom in Inernal Coordinates.
     * We create a list of OBInternalCoords that match the list of atoms in the molecule
     * that are defined using internal coordinates
     */

    double var;
    bool ok=false;
    int n;

    vic.push_back(new OBInternalCoord);
    atom->SetAtomicNum(LabelToAtomicNumber(tokens[0]));

    switch (*zmatLineCount) {
    case 0:
      break;

    case 1:
      if (tokens.size() < 3) {return false;}

      // Specify the atom that defines the distance to this one
      ok = from_string<int>(n, tokens.at(1), std::dec);
      vic[*zmatLineCount]->_a = mol.GetAtom(n);

      // Get the distance
      ok = from_string<double>(var, tokens.at(2), std::dec);
      if ( !ok )
        {
          // Can't convert to double so see if it's in the variables
          if (variables.find(tokens[2])==variables.end()) return false;
          var = variables[tokens[2]];
        }
      vic[*zmatLineCount]->_dst = var;
      break;

    case 2:
      if (tokens.size() < 5) {return false;}

      // Specify the atom that defines the distance to this one
      ok = from_string<int>(n, tokens.at(1), std::dec);
      vic[*zmatLineCount]->_a = mol.GetAtom(n);

      // Get the distance
      ok = from_string<double>(var, tokens.at(2), std::dec);
      if ( !ok )
        {
          // Can't convert to double so see if it's in the variables
          if (variables.find(tokens[2])==variables.end()) return false;
          var = variables[tokens[2]];
        }
      vic[*zmatLineCount]->_dst = var;

      // Specify atom defining angle
      ok = from_string<int>(n, tokens.at(3), std::dec);
      vic[*zmatLineCount]->_b = mol.GetAtom(n);
      // Get the angle
      ok = from_string<double>(var, tokens.at(4), std::dec);
      if ( !ok )
        {
          // Can't convert to double so see if it's in the variables
          if (variables.find(tokens[4])==variables.end()) return false;
          var = variables[tokens[4]];
        }
      vic[*zmatLineCount]->_ang = var;
      break;

    default:
      if (tokens.size() < 7) {return false;}

      ok = from_string<int>(n, tokens.at(1), std::dec);
      vic[*zmatLineCount]->_a = mol.GetAtom(n);
      // Get the distance
      ok = from_string<double>(var, tokens.at(2), std::dec);
      if ( !ok )
        {
          // Can't convert to double so see if it's in the variables
          if (variables.find(tokens[2])==variables.end()) return false;
          var = variables[tokens[2]];
        }
      vic[*zmatLineCount]->_dst = var;

      ok = from_string<int>(n, tokens.at(3), std::dec);
      vic[*zmatLineCount]->_b = mol.GetAtom(n);
      // Get the angle
      ok = from_string<double>(var, tokens.at(4), std::dec);
      if ( !ok )
        {
          // Can't convert to double so see if it's in the variables
          if (variables.find(tokens[4])==variables.end()) return false;
          var = variables[tokens[4]];
        }
      vic[*zmatLineCount]->_ang = var;
      
      ok = from_string<int>(n, tokens.at(5), std::dec);
      vic[*zmatLineCount]->_c = mol.GetAtom(n);
      // Get the torsion angle
      ok = from_string<double>(var, tokens.at(6), std::dec);
      if ( !ok )
        {
          // Can't convert to double so see if it's in the variables
          if (variables.find(tokens[6])==variables.end()) return false;
          var = variables[tokens[6]];
        }
      vic[*zmatLineCount]->_tor = var;
    }

    (*zmatLineCount)++;
    return true;
  }

  bool GAMESSUKFormat::ReadVariables(istream &ifs, double factor, string stopstr)
  {
    /*
     * This takes an input stream that is positioned where the list of variables
     * starts and the reads the variables into the supplied map
     *
     * This is different to ReadGeometry (which takes a vector of strings as input) because
     * currently the variables always need to be read after the geometry, so we need to save the
     * geometry and then read the variables. However this means that we can parse the variables
     * directly into a map and don't need to keep a copy of the specifcation as strings.
     *
     * stopstr is a string that defines when we stop reading
     */

    string line;
    vector<string> tokens;
    bool ok=false;
    double var;

    // Now read in all the varibles
    while (ifs.good() && ifs.getline(buffer, BUFF_SIZE)) {

      // Skip commnents
      if (EQn(buffer, "#", 1) || EQn(buffer, "?", 1))
        continue;

      // Copy line to a C++ string and convert to lower case
      // & remove leading and trailing spaces
      line = buffer;
      // transform(method.begin(), method.end(), method.begin(), ::tolower);
      ToLower(line);
      Trim(line);

      // Check for end of variables
      if (line.length()==0 && stopstr.length()==0) break;
      if (stopstr.length()>0 && line.compare(0, stopstr.length(), stopstr)==0) break;

      // Check for commas & split with that as the separator if necessary
      if (line.find(',')!=string::npos) {
        tokenize(tokens, line, ",");
      } else {
        tokenize(tokens, line, " \t\n");
      }

      ok = from_string<double>(var, tokens.at(3), std::dec);
      if ( !ok )
        {
          errorMsg << "Problems reading a GAMESS-UK  file: "
                   << "Could not read variable line: " << line;
          obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
          return false;
        }
      // Add to list of variables
      variables[tokens[0]]=var*factor;
    }

    /*
      cerr << "Got list of variables: " << endl;
      for (map<string,double>::iterator i=variables.begin(); i
      != variables.end(); i++) {
      cerr << "Name: " << i->first << " Value: " << i->second << endl;
      }
    */

    return true;

  } // end Read Variables


  class GAMESSUKInputFormat : public OBMoleculeFormat, public GAMESSUKFormat
  {
  public:
    //Register this format type ID
    GAMESSUKInputFormat()
    {
      OBConversion::RegisterFormat("gukin",this, "chemical/x-gamess-input");
      // Command-line keywords
      //OBConversion::RegisterOptionParam("k", NULL, 1, OBConversion::OUTOPTIONS);
      // Command-line keyword file
      //OBConversion::RegisterOptionParam("f", NULL, 1, OBConversion::OUTOPTIONS);
    }


    virtual const char* Description() //required
    {
      return
        "GAMESS-UK Input\n";
    };

    virtual const char* SpecificationURL()
    {return "http://www.cfs.dl.ac.uk";}; //optional

    virtual const char* GetMIMEType()
    { return "chemical/x-gamessuk-input"; };

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
      return READONEONLY; // | NOTREADABLE;
    };

    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
  };

  //Make an instance of the format class
  GAMESSUKInputFormat theGAMESSUKInputFormat;


  bool GAMESSUKInputFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)

  {
    /*
     * Stuff to think about:
     * - At outset check whether we are in zmatrix, cartesian or nw-chem format
     *   (we currently only suppot homogeneous formats - not mixed).
     *
     * For each line need to check:
     * - Is this a comment (?,#)
     * - Are the tokens separated by commas, if so use tokenize to split at commas
     * - Is there an 'end' token on the line
     *
     * For each line we want to check that we haven't hit a change from c to zm or vv
     *
     */

    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if (pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream& ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;

    // Get a default title as the filename
    const char* title = pConv->GetTitle();
    mol.BeginModify();
    mol.SetTitle(title);
    mol.EndModify();

    vector<string> geomList, tokens; // list of lines and list of tokens on a line
    string line; // For convenience so we can refer to lines from the iterator as 'line'
    ReadMode_t ReadMode=SKIP;
    double factor=BOHR_TO_ANGSTROM;

    // Read File and copy geometry specification into geomList
    while (ifs.good() && ifs.getline(buffer, BUFF_SIZE)) {

      // Skip commnents
      if (EQn(buffer, "#", 1) || EQn(buffer, "?", 1)) continue;

      // Copy line to a C++ string and convert to lower case
      // & remove leading and trailing spaces
      line = buffer;
      // transform(method.begin(), method.end(), method.begin(), ::tolower);
      ToLower(line);
      Trim(line);

      // Start of coordinate specifiation
      if (line.compare(0, 4, "zmat")==0)
	{
	  ReadMode=ZMATRIX;
	  geomList.push_back(line);
	  continue;
	}
      else if (line.compare(0, 4, "geom")==0)
	{
	  ReadMode=CARTESIAN;
	  geomList.push_back(line);
	  continue;
	}

      // Reading the coordinate specification into the list
      if (ReadMode==ZMATRIX || ReadMode==CARTESIAN)
	{

	  // Variables specification - process directly from filestream
	  // and then remove from the geometry specification
	  if (line.compare(0, 4, "vari")==0 || line.compare(0, 4, "const")==0)
	    {

	      // Check for commas & split with that as the separator if necessary
	      if (line.find(',')!=string::npos)
		tokenize(tokens, line, ",");
	      else
		tokenize(tokens, line, " \t\n");

	      // See if we need to rescale
	      if (IsUnits(tokens[1])) factor=Rescale(tokens[1]);

	      if (! ReadVariables(ifs, factor, "end")) return false;

	      ReadMode=SKIP;
	      geomList.push_back("end\n");
	      continue;
	    }

	  if (line.compare(0, 3, "end")==0) ReadMode=SKIP;

	  geomList.push_back(line);
	}

    }// End while reading loop

    // Now go and process the coordinate specification if we got any
    bool ok = ReadGeometry(mol, geomList);

    if (mol.NumAtoms() == 0) { // e.g., if we're at the end of a file PR#1737209
      mol.EndModify();
      return false;
    } else {
      if (!pConv->IsOption("b",OBConversion::INOPTIONS))
        mol.ConnectTheDots();
      if (!pConv->IsOption("s",OBConversion::INOPTIONS) && !pConv->IsOption("b",OBConversion::INOPTIONS))
        mol.PerceiveBondOrders();
      return ok;
    }

  } // End ReadMolecule

  ////////////////////////////////////////////////////////////////

  bool GAMESSUKInputFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    char buffer[BUFF_SIZE];

    ofs << "title" << endl;
    ofs << mol.GetTitle() << endl << endl;

    ofs << "#" << endl;
    ofs << "# NB: Class I directives (e.g. memory, multiplicity, charge etc) go here" << endl;
    ofs << "#" << endl;
    ofs << "# For more information see: http://www.cfs.dl.ac.uk/docs/index.shtml" << endl;
    ofs << "#" << endl;
    ofs << endl;

    ofs << "geometry angstrom" << endl;
    FOR_ATOMS_OF_MOL(atom, mol)
      {
        snprintf(buffer, BUFF_SIZE, "%15.8f %15.8f %15.8f %3d %3s\n",
		 atom->GetX(),
		 atom->GetY(),
		 atom->GetZ(),
		 atom->GetAtomicNum(),
		 OBElements::GetSymbol(atom->GetAtomicNum())
		 );
	ofs << buffer;
      }
    ofs << "end" << endl << endl;

    ofs << endl;
    ofs << "basis 6-31G" << endl;
    ofs << endl;

    ofs << "#" << endl;
    ofs << "# NB: Class II directives go here" << endl;
    ofs << "#" << endl;
    ofs << "# To perform a dft calculation with b3lyp and medium quadrature uncomment the below" << endl;
    ofs << "# dft b3lyp" << endl;
    ofs << "# dft quadrature medium" << endl;
    ofs << "#" << endl;
    ofs << endl;

    ofs << "runtype scf" << endl;
    ofs << endl;
    ofs << "enter" << endl;

    return(true);
  } //End WriteMolecule


  class GAMESSUKOutputFormat : public OBMoleculeFormat, public GAMESSUKFormat
  {
  public:
    //Register this format type ID
    GAMESSUKOutputFormat()
    { OBConversion::RegisterFormat("gukout",this, "chemical/x-gamess-output"); }

    virtual const char* Description() //required
    { return "GAMESS-UK Output\n"; };

    virtual const char* SpecificationURL()
    {return "http://www.cfs.dl.ac.uk";}; //optional

    virtual const char* GetMIMEType()
    { return "chemical/x-gamessuk-output"; };

    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);

  private:
    enum RunType_t { UNKNOWN, SINGLEPOINT, OPTXYZ, OPTZMAT, SADDLE, FREQUENCIES };
    vector<string> tokens, geomList; // list of lines and list of tokens on a line
    string line; // For convenience so we can refer to lines from the iterator as 'line'
    bool ReadInputZmatrix( OBMol &mol, std::istream &ifs );
    bool ReadInitialCartesian( OBMol &mol, std::istream &ifs );
    bool ReadOptGeomXyz1( OBMol &mol, std::istream &ifs );
    bool ReadOptGeomXyz2( OBMol &mol, std::istream &ifs );
    bool ReadNormalModesHessian( OBMol &mol, std::istream &ifs);
    bool ReadNormalModesForce( OBMol &mol, std::istream &ifs);
  };

  //Make an instance of the format class
  GAMESSUKOutputFormat theGAMESSUKOutputFormat;

  bool GAMESSUKOutputFormat::ReadInputZmatrix( OBMol &mol, std::istream &ifs )
  {
    /* The zmatrix entered by the user
     * REM:  need to add stuff for "automatic z-matrix generation" as we currently
     * ignore the zmatrix & just read the cartesian coordinates
     */
    geomList.clear();
    
    // skip 2 lines
    ifs.getline(buffer, BUFF_SIZE) && ifs.getline(buffer, BUFF_SIZE);
    
    // Stick a header line first
    geomList.push_back("zmatrix bohr");
    
    // Read zmatrix into list until blank line
    while (ifs.good() && ifs.getline(buffer, BUFF_SIZE) && strlen(buffer) != 0)
      {
        line = buffer;
        // transform(method.begin(), method.end(), method.begin(), ::tolower);
        ToLower(line);
        Trim(line);
        geomList.push_back(line);
      }
      
    // Skip 2 lines
    ifs.getline(buffer, BUFF_SIZE);
    ifs.getline(buffer, BUFF_SIZE);
      
    // Check if line is variables line
    if (strstr(buffer,"name            input  type     hessian         minima") != NULL)
      {
        // Skip additional line to be where variables are printed
        ifs.getline(buffer, BUFF_SIZE);
        // Read in the variables till we hit blank line
        if (! ReadVariables(ifs, 1.0, "")) return false;
      }
      
    // Now go and process the geometry
    return ReadGeometry(mol, geomList);
  } // ReadInputZmatrix

  bool GAMESSUKOutputFormat::ReadInitialCartesian( OBMol &mol, std::istream &ifs )
  {
    bool ok=false;
    double x,y,z;
    int n;

    // Skip 3 lines
    ifs.getline(buffer, BUFF_SIZE) &&
      ifs.getline(buffer, BUFF_SIZE) &&
      ifs.getline(buffer, BUFF_SIZE);

    // Create regex for the coords
    //                     ------label--------   -------charge-------- < seems enough for a match
    string pattern(" *\\* *[a-zA-Z]{1,2}[0-9]* *[0-9]{1,3}\\.[0-9]{1}");
    bool iok;
#ifdef _MSC_VER
    std::tr1::regex myregex;
    try {
      myregex.assign(pattern,
                     std::tr1::regex_constants::extended |
                     std::tr1::regex_constants::nosubs);
      iok = true;
    } catch (std::tr1::regex_error ex) {
      iok = false;
    }
#else
    regex_t *myregex = new regex_t;
    iok = regcomp(myregex, pattern.c_str(), REG_EXTENDED | REG_NOSUB)==0;
#endif
    if (!iok) cerr << "Error compiling regex in GUK OUTPUT!\n";

    // Read in the coordinates - we process them directly rather
    // then use ReadGeometry as we probably should do...
    mol.BeginModify();
    while (ifs.good() && ifs.getline(buffer, BUFF_SIZE)){

      // End of geometry block
      if (strstr(buffer,"*************************")!=NULL)break;
#ifdef _MSC_VER
      if (std::tr1::regex_search(buffer, myregex)) {
#else
        if (regexec(myregex, buffer, 0, 0, 0)==0) {
#endif
          //cerr << "Got Coord line: " << buffer << endl;
          OBAtom *atom = mol.NewAtom();
          tokenize(tokens,buffer," ");

          ok = from_string<int>(n, tokens.at(2), std::dec);
          atom->SetAtomicNum(n);
          ok = from_string<double>(x, tokens.at(3), std::dec);
          x=x*BOHR_TO_ANGSTROM;
          ok = from_string<double>(y, tokens.at(4), std::dec);
          y=y*BOHR_TO_ANGSTROM;
          ok = from_string<double>(z, tokens.at(5), std::dec);
          z=z*BOHR_TO_ANGSTROM;
          atom->SetVector(x, y, z);
        }
      }
      mol.EndModify();
#ifndef _MSC_VER
      regfree(myregex);
#endif
      return true;
    } // End ReadInitalCartesian


    bool GAMESSUKOutputFormat::ReadOptGeomXyz1( OBMol &mol, std::istream &ifs )
    {
      bool ok=false;
      double x,y,z;
      int n;
      
      // Clear the Molecule as we're going to start from scratch again.
      mol.BeginModify();
      mol.Clear();

      // FF to start of coordinate specification
      while (ifs.good() && ifs.getline(buffer, BUFF_SIZE)) {
        if (strstr(buffer,
                   "atom     znuc       x             y             z") != NULL) break;
      }

      // Skip 2 lines - should then be at the coordinates
      ifs.getline(buffer, BUFF_SIZE) &&
        ifs.getline(buffer, BUFF_SIZE);

      // Read in the coordinates - we process them directly rather
      // then use ReadGeometry as we probably should do...
      while (ifs.good() && ifs.getline(buffer, BUFF_SIZE)){

        // End of geometry block
        if (strstr(buffer,"*************************")!=NULL)break;

        //cerr << "Got Coord line: " << buffer << endl;
        OBAtom *atom = mol.NewAtom();
        tokenize(tokens,buffer," ");

        ok = from_string<int>(n, tokens.at(2), std::dec);
        atom->SetAtomicNum(n);
        ok = from_string<double>(x, tokens.at(3), std::dec);
        x=x*BOHR_TO_ANGSTROM;
        ok = from_string<double>(y, tokens.at(4), std::dec);
        y=y*BOHR_TO_ANGSTROM;
        ok = from_string<double>(z, tokens.at(5), std::dec);
        z=z*BOHR_TO_ANGSTROM;
        atom->SetVector(x, y, z);

      }

      mol.EndModify();
      return true;
    } // End ReadOptGeomXyz

    bool GAMESSUKOutputFormat::ReadOptGeomXyz2( OBMol &mol, std::istream &ifs )
    {
      bool ok=false;
      double x,y,z;
      int n;

      // Clear the Molecule as we're going to start from scratch again.
      mol.BeginModify();
      mol.Clear();

      while (ifs.good() && ifs.getline(buffer, BUFF_SIZE)) {
        if (strstr(buffer,
                   "       x              y              z            chg  tag") != NULL) break;
      }

      // Skip 1 line - should then be at the coordinates
      ifs.getline(buffer, BUFF_SIZE);

      // Read in the coordinates - we process them directly rather
      // then use ReadGeometry as we probably should do...
      while (ifs.good() && ifs.getline(buffer, BUFF_SIZE)){

        // End of geometry block
        if (strstr(buffer,"============================================================")!=NULL)break;

        //cerr << "Got Coord line: " << buffer << endl;
        OBAtom *atom = mol.NewAtom();
        tokenize(tokens,buffer," ");

        ok = from_string<int>(n, tokens.at(3), std::dec);
        atom->SetAtomicNum(n);
        ok = from_string<double>(x, tokens.at(0), std::dec);
        x=x*BOHR_TO_ANGSTROM;
        ok = from_string<double>(y, tokens.at(1), std::dec);
        y=y*BOHR_TO_ANGSTROM;
        ok = from_string<double>(z, tokens.at(2), std::dec);
        z=z*BOHR_TO_ANGSTROM;
        atom->SetVector(x, y, z);
      }

      mol.EndModify();
      return true;

    } // End ReadOptGeomZmat

    bool GAMESSUKOutputFormat::ReadNormalModesHessian( OBMol &mol, std::istream &ifs)
    {

      bool ok=false;
      double dtmp,dtmp2;

      int ncols = 8; // Think this is always the case
      int natoms = mol.NumAtoms();
      int maxroot = natoms*3;
    
      // Create data structures
      std::vector< double > frequencies, intensities;
      std::vector< std::vector< vector3 > > Lx;
      
      // Set up data structures with null data
      for( int i=0; i<maxroot; i++ )
        {
          std::vector< vector3 > atoml;
          for( int j=0; j < natoms; j++ )
            {
              atoml.push_back( vector3(0.0,0.0,0.0) );
            }
          Lx.push_back( atoml );
        }
      
      ifs.getline(buffer, BUFF_SIZE); // skip "===============" line

      int root7;
      for ( int root1=0; root1 < maxroot; root1+=ncols )
        {
          root7 = root1 + ncols;
          root7 = min(maxroot,root7);

          //Skip 6 lines to col header with frequencies
          for( int j=0; j < 6; j++ )
            ifs.getline(buffer, BUFF_SIZE);

          tokenize(tokens,buffer," \t\n");
          for( std::size_t si=0; si < tokens.size(); si++ )
            {
              ok = from_string<double>(dtmp, tokens.at(si), std::dec);
              frequencies.push_back(dtmp);
              intensities.push_back( 0.0 ); // Add placeholder data
            }

          // Skip 2 lines to where data matrix starts
          ifs.getline(buffer, BUFF_SIZE);
          ifs.getline(buffer, BUFF_SIZE);

          int mycols=root7-root1;
          // Loop over atoms & the x,y,z
          int atomcount=0;
          for ( int i=0; i<maxroot; i+=3 )
            {
              for ( int j=0; j<3; j++ )
                {
                  ifs.getline(buffer, BUFF_SIZE);
                  //std::cout << "GOT LINE:" << buffer <<std::endl;
                  tokenize(tokens,buffer," \t\n");
                  int start=1;
                  if ( j == 0 )
                    start=3;
                  for ( int k=0; k<mycols; k++ )
                    {
                      //std::cout << "Lx[ " << root1+k << " ]" <<
                      //  "][ " << atomcount << " ] [ " << j << " ] = " << tokens.at(start+k) << std::endl;
                      ok = from_string<double>(dtmp, tokens.at(start+k), std::dec);
                      if ( j==0)
                        Lx[ root1+k ][ atomcount ].SetX( dtmp );
                      else if ( j==1 )
                        Lx[ root1+k ][ atomcount ].SetY( dtmp );
                      else if ( j==2 )
                        Lx[ root1+k ][ atomcount ].SetZ( dtmp );
                    }
                } // End j loop
              atomcount+=1;
            } // End loop over atoms
        } // loop over root1

      // Now skip down to read in intensities
      for( int i=0; i<7; i++ )
        {
          ifs.getline(buffer, BUFF_SIZE);
        }
      
      // loop until we've read them all in
      for( std::size_t si=0; si<frequencies.size(); si++ )
        {
          ifs.getline(buffer, BUFF_SIZE);
          // End of info
          if (strstr(buffer,"============")!=NULL)break;
          tokenize(tokens,buffer," \t\n");
          
          ok = from_string<double>(dtmp, tokens.at(1), std::dec);
          ok = from_string<double>(dtmp2, tokens.at(6), std::dec);
          // Now match them up
          for( std::size_t sj=0; sj<frequencies.size(); sj++ )
            {
              if ( std::abs( frequencies.at(sj) - dtmp ) < 0.01 )
                {
                  intensities[sj]= dtmp2;
                  continue;
                }
            }
        }

      //for (int i=0; i< frequencies.size(); i++ )
      //  std::cout << "Frequency: " << frequencies.at(i) << " : " << intensities.at(i) << std::endl;
      
      if(frequencies.size()>0)
        {
          OBVibrationData* vd = new OBVibrationData;
          vd->SetData(Lx, frequencies, intensities);
          vd->SetOrigin(fileformatInput);
          mol.SetData(vd);
        }
      return ok;
    } // End ReadNormalModesHessian

    bool GAMESSUKOutputFormat::ReadNormalModesForce( OBMol &mol, std::istream &ifs)
    {

      bool ok=false;
      double dtmp;

      int ncols = 9; // Think this is always the case
      int natoms = mol.NumAtoms();
      int maxroot = natoms*3;
      int start,mycols;

      // Create data structures
      std::vector< double > frequencies, intensities;
      std::vector< std::vector< vector3 > > Lx;
      
      // Set up data structures with null data
      for( int i=0; i<maxroot; i++ )
        {
          std::vector< vector3 > atoml;
          for( int j=0; j < natoms; j++ )
            {
              atoml.push_back( vector3(0.0,0.0,0.0) );
            }
          Lx.push_back( atoml );
        }
      
      ifs.getline(buffer, BUFF_SIZE); // skip "===============" line

      int root7;
      for ( int root1=0; root1 < maxroot; root1+=ncols )
        {
          root7 = root1 + ncols;
          root7 = min(maxroot,root7);
          mycols=root7-root1;

          //Skip 6 lines to col header with frequencies
          for( int j=0; j < 6; j++ )
            ifs.getline(buffer, BUFF_SIZE);

          line=buffer;
          // Need  to manually chop up line
          start=20; // Numbers start at col 20 & are 11 characters long
          for( int i=0; i<mycols; i++)
            {
              ok = from_string<double>(dtmp, line.substr(start,12), std::dec);
              frequencies.push_back(dtmp);
              intensities.push_back( 10.0 ); // Intensities aren't printed so just use 10
              start+=12;
            }

          // Skip 2 lines to where data matrix starts
          ifs.getline(buffer, BUFF_SIZE);
          ifs.getline(buffer, BUFF_SIZE);

          // Loop over atoms & the x,y,z
          int atomcount=0;
          for ( int i=0; i<maxroot; i+=3 )
            {
              //for j in range(3):
              for ( int j=0; j<3; j++ )
                {
                  ifs.getline(buffer, BUFF_SIZE);
                  //std::cout << "GOT LINE:" << buffer <<std::endl;
                  tokenize(tokens,buffer," \t\n");
                  start=1;
                  if ( j == 0 )
                    start=3;
                  for ( int k=0; k<mycols; k++ )
                    {
                      // std::cout << "Lx[ " << root1+k << " ]" <<
                      //  "][ " << atomcount << " ] [ " << j << " ] = " << tokens.at(start+k) << std::endl;
                      ok = from_string<double>(dtmp, tokens.at(start+k), std::dec);
                      if ( j==0)
                        Lx[ root1+k ][ atomcount ].SetX( dtmp );
                      else if ( j==1 )
                        Lx[ root1+k ][ atomcount ].SetY( dtmp );
                      else if ( j==2 )
                        Lx[ root1+k ][ atomcount ].SetZ( dtmp );
                    }
                } // End j loop
              atomcount+=1;
            } // End loop over atoms
        } // loop over root1


      //for (int i=0; i< frequencies.size(); i++ )
      //  std::cout << "Frequency: " << frequencies.at(i) << " : " << intensities.at(i) << std::endl;

      if(frequencies.size()>0)
        {
          OBVibrationData* vd = new OBVibrationData;
          vd->SetData(Lx, frequencies, intensities);
          vd->SetOrigin(fileformatInput);
          mol.SetData(vd);
        }
      return ok;
    } // End ReadNormalModesForce

    /*
      bool GAMESSUKOutputFormat::ReadOptGeomZmat( OBMol &mol, std::istream &ifs )
      {
  
      //Below was for reading in zmatricies - ignore for the time being

      // Original geometry specification should still be in geomList
      // So just update the variables
      //cerr << "Got converged for OPTZMAT\n";

      // FF to variable specification
      while (ifs.good() && ifs.getline(buffer, BUFF_SIZE)) {
      if (strstr(buffer,
      " variable           value                hessian") != NULL) break;
      }
      // Skip a line - should then be at variable specification
      ifs.getline(buffer, BUFF_SIZE);

      // Process them
      if (! ReadVariables(ifs, BOHR_TO_ANGSTROM,
      "===============================================")) return false;

      // Now go and process with the geometry we read before
      return ReadGeometry(mol, geomList);

      } //ReadOptGeomZmat
    */

    bool GAMESSUKOutputFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv) {

      /*
        Read a GAMESS-UK output file. The reader is currently set up to only return one molecule, i.e
        if the run is some sort of optimisation run, then only the optimised geometry is returned.
        
        Previously we parsed in the z-matrix - and the code to do that is still here and should work.
        However, the zmatrix is not actually used anywhere in OpenBabel currently - it's just converted
        to Cartesians, and this code appears to be buggy, so for the time being we just stick to reading
        Cartesians.
      */

      OBMol *pmol = dynamic_cast<OBMol*>(pOb);
      if (pmol==NULL)
        return false;

      //Define some references so we can use the old parameter names
      istream& ifs = *pConv->GetInStream();
      OBMol &mol = *pmol;

      // Get a default title as the filename
      const char* title = pConv->GetTitle();
      mol.BeginModify();
      mol.SetTitle(title);
      mol.EndModify();

      RunType_t RunType=UNKNOWN;
      bool ok;
      std::string runt;

      while (ifs.good() && ifs.getline(buffer, BUFF_SIZE))
        {

          if (strstr(buffer,"                              input z-matrix") != NULL)
            {
              /* OpenBabel's handling of zmatricies is currently too buggy and the zmatrix
               * read in isn't currently used - it's just converted to cartesians, so we
               * can skip this for the time being
               */
              continue;
              /*
                ok = ReadInputZmatrix( mol, ifs );
                // Set Runtype to SINGLEPOINT so we don't read in the cartesians
                RunType=SINGLEPOINT;
              */
            } // End Reading user z-matrix
          
          // Read the cartesian coordinates if we've not read in the ZMATRIX
          if (strstr(buffer,"*            charge       x             y              z       shells") != NULL &&
              RunType==UNKNOWN)
            ok = ReadInitialCartesian( mol, ifs );
          
          // Determine the RunType - affects how we move on from here.
          if (strstr(buffer," * RUN TYPE") != NULL)
            {
              tokenize(tokens,buffer," \t\n");
              runt=tokens[3].substr(0,5);
              if(runt=="optxy") RunType=OPTXYZ;
              else if (runt=="optim") RunType=OPTZMAT;
              else if (runt=="saddl") RunType=SADDLE;
              continue;
            } // End RUNTYPE
          
          // Read the optimised geometry
          if (strstr(buffer,"optimization converged") != NULL)
            {
              if (RunType==OPTXYZ)
                ok = ReadOptGeomXyz1( mol, ifs );
              else if (RunType==OPTZMAT || RunType==SADDLE)
                ok = ReadOptGeomXyz2( mol, ifs );
            } // End read optimised geometry

          // Frequencies for runtype hessian
          if (strstr(buffer,"cartesians to normal") != NULL)
            ok = ReadNormalModesHessian( mol, ifs);

          // Frequencies for runtype force
          if (strstr(buffer,"eigenvectors of cartesian") != NULL)
            ok = ReadNormalModesForce( mol, ifs);
          
        } // End Reading loop
    
      if (mol.NumAtoms() == 0) { // Something went wrong
        mol.EndModify();
        return false;
      } else {
        mol.BeginModify();
        if (!pConv->IsOption("b",OBConversion::INOPTIONS))
          mol.ConnectTheDots();
        if (!pConv->IsOption("s",OBConversion::INOPTIONS) && !pConv->IsOption("b",OBConversion::INOPTIONS))
          mol.PerceiveBondOrders();
        mol.EndModify();
        return true;
      }
    
    } // End GAMESSUKOutputFormat::ReadMolecule
  
  
  } //namespace OpenBabel
