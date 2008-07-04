/**********************************************************************
  Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
  Some portions Copyright (C) 2004 by Chris Morley
  Some portions Copyright (C) 2006 by Donald E. Curtis

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

#include <algorithm>

#include <regex.h>

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
	  const double Rescale(string text);
	  bool IsUnits(string text);
	  enum ReadMode_t {CARTESIAN, ZMATRIX, VARIABLES, CONSTANTS, SKIP};
	  ReadMode_t ReadMode;	
	  char buffer[BUFF_SIZE];
	  stringstream errorMsg;
		
private:
	  map<string, double> variables; // map from variable name to value
	  vector<OBInternalCoord*> vic; // Holds lists of internal coordinates
	  const int LabelToAtomicNumber(string label);
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
	
	for (vector<string>::iterator i=geomList.begin(); i !=geomList.end(); i++) {
		
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
}

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


const double GAMESSUKFormat::Rescale(string text)
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

const int GAMESSUKFormat::LabelToAtomicNumber(string label)
{
	/* 
	 * Given a string with the label for an atom return the atomic number 
	 * As we are using the GetAtomicNum function case is not important
	 */

	// See if the first 2 characters give us a valid atomic #
	int Z=etab.GetAtomicNum(label.substr(0,2).c_str());
	
	// If not try the first one
	if (Z==0) Z=etab.GetAtomicNum(label.substr(0,1).c_str());

	// Houston...
	if (Z==0){
		errorMsg << "LabelToAtomicNumber got bad Label: " << label;
		obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
	}
	return Z;
}

bool GAMESSUKFormat::ReadLineCartesian(OBAtom *atom, vector<string> &tokens, double factor)
{
	
	/*  Read a line defining the Cartesian coordinates for an atom
	 * This assumes the line is formatted in GAMESS-UK input style as:
	 * x y z AtomicNumber Label
	 */

	// 4th field is the atomic number
	atom->SetAtomicNum(atoi(tokens[3].c_str()));

	// Read the atom coordinates
	char *endptr;
	double x = strtod((char*)tokens[0].c_str(), &endptr);
	if (endptr == (char*)tokens[0].c_str()) {
		// Can't convert to double so see if it's in the variables
		if (variables.find(tokens[0])==variables.end()) return false;
		x = variables[tokens[0]];
	}

	double y = strtod((char*)tokens[1].c_str(), &endptr);
	if (endptr == (char*)tokens[1].c_str()) {
		// Can't convert to double so see if it's in the variables
		if (variables.find(tokens[1])==variables.end()) return false;
		y = variables[tokens[1]];
	}

	double z = strtod((char*)tokens[2].c_str(), &endptr);
	if (endptr == (char*)tokens[2].c_str()) {
		// Can't convert to double so see if it's in the variables
		if (variables.find(tokens[2])==variables.end())  return false;
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
	
	char *endptr;
	double var;
	
	vic.push_back(new OBInternalCoord);
	atom->SetAtomicNum(LabelToAtomicNumber(tokens[0]));
	
	switch (*zmatLineCount) {
       case 0:
    	   break;

       case 1:
    	   if (tokens.size() < 3) {return false;}
    	   
    	   // Specify the atom that defines the distance to this one
    	   vic[*zmatLineCount]->_a = mol.GetAtom(atoi(tokens[1].c_str()));
    	   
    	   // Get the distance
    	   	var = strtod((char*)tokens[2].c_str(), &endptr);
    		if (endptr == (char*)tokens[2].c_str()) {
    			// Can't convert to double so see if it's in the variables
    			if (variables.find(tokens[2])==variables.end()) return false;
    			var = variables[tokens[2]];
    		}
    		vic[*zmatLineCount]->_dst = var;
    	   break;

       case 2:
    	   if (tokens.size() < 5) {return false;}
    	   
    	   // Specify the atom that defines the distance to this one
    	   vic[*zmatLineCount]->_a = mol.GetAtom(atoi(tokens[1].c_str()));
    	   // Get the distance
    	   	var = strtod((char*)tokens[2].c_str(), &endptr);
    		if (endptr == (char*)tokens[2].c_str()) {
    			// Can't convert to double so see if it's in the variables
    			if (variables.find(tokens[2])==variables.end()) return false;
    			var = variables[tokens[2]];
    		}
    		vic[*zmatLineCount]->_dst = var;
    		
    		// Specify atom defining angle
    	   vic[*zmatLineCount]->_b = mol.GetAtom(atoi(tokens[3].c_str()));
    	   // Get the angle
    	   	var = strtod((char*)tokens[4].c_str(), &endptr);
    		if (endptr == (char*)tokens[4].c_str()) {
    			// Can't convert to double so see if it's in the variables
    			if (variables.find(tokens[4])==variables.end()) return false;
    			var = variables[tokens[4]];
    		}
    		vic[*zmatLineCount]->_ang = var;   	   
    	   break;

       default:
    	   if (tokens.size() < 7) {return false;}
    	   
    	   vic[*zmatLineCount]->_a = mol.GetAtom(atoi(tokens[1].c_str()));
    	   // Get the distance
    	   	var = strtod((char*)tokens[2].c_str(), &endptr);
    		if (endptr == (char*)tokens[2].c_str()) {
    			// Can't convert to double so see if it's in the variables
    			if (variables.find(tokens[2])==variables.end()) return false;
    			var = variables[tokens[2]];
    		}
    		vic[*zmatLineCount]->_dst = var;    	   
    	   
    	   vic[*zmatLineCount]->_b = mol.GetAtom(atoi(tokens[3].c_str()));
    	   // Get the angle
    	   	var = strtod((char*)tokens[4].c_str(), &endptr);
    		if (endptr == (char*)tokens[4].c_str()) {
    			// Can't convert to double so see if it's in the variables
    			if (variables.find(tokens[4])==variables.end()) return false;
    			var = variables[tokens[4]];
    		}
    		vic[*zmatLineCount]->_ang = var;      	   
    	   
    	   vic[*zmatLineCount]->_c = mol.GetAtom(atoi(tokens[5].c_str()));	
    	   // Get the torsion angle
    	   	var = strtod((char*)tokens[6].c_str(), &endptr);
    		if (endptr == (char*)tokens[6].c_str()) {
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
		if (line.length()==0 and stopstr.length()==0) break;
		if (stopstr.length()>0 && line.compare(0, stopstr.length(), stopstr)==0) break;
		
		// Check for commas & split with that as the separator if necessary
		if (line.find(',')!=string::npos) {
			tokenize(tokens, line, ",");
		} else {
			tokenize(tokens, line, " \t\n");
		}

		char *endptr;
		double var = strtod((char*)tokens[1].c_str(), &endptr);
		if (endptr == (char*)tokens[1].c_str()) {
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
	
}


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
  //virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
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
		if (line.compare(0, 4, "zmat")==0) {
			ReadMode=ZMATRIX;
			geomList.push_back(line);
			continue;
		} else if (line.compare(0, 4, "geom")==0) {
			ReadMode=CARTESIAN;
			geomList.push_back(line);
			continue;
		}

		
		// Reading the coordinate specification into the list
		if (ReadMode==ZMATRIX || ReadMode==CARTESIAN) {
		
			// Variables specification - process directly from filestream
			// and then remove from the geometry specification
			if (line.compare(0, 4, "vari")==0 || line.compare(0, 4, "const")==0) {
				
				// Check for commas & split with that as the separator if necessary
				if (line.find(',')!=string::npos) {
					tokenize(tokens, line, ",");
				} else {
					tokenize(tokens, line, " \t\n");
				}	

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
    	return ok;
    }

} // End ReadMolecule	
	
	class GAMESSUKOutputFormat : public OBMoleculeFormat, public GAMESSUKFormat
	{
	public:
	  //Register this format type ID
	  GAMESSUKOutputFormat()
	  {
	    OBConversion::RegisterFormat("gukout",this, "chemical/x-gamess-output");
	    // Command-line keywords
	    //OBConversion::RegisterOptionParam("k", NULL, 1, OBConversion::OUTOPTIONS);
	    // Command-line keyword file
	    //OBConversion::RegisterOptionParam("f", NULL, 1, OBConversion::OUTOPTIONS);
	  }


	  virtual const char* Description() //required
	  {
	    return
	      "GAMESS-UK Output\n";
	  };

	  virtual const char* SpecificationURL()
	  {return "http://www.cfs.dl.ac.uk";}; //optional

	  virtual const char* GetMIMEType() 
	  { return "chemical/x-gamessuk-output"; };


	  ////////////////////////////////////////////////////
	  /// The "API" interface functions
	  //virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
	  virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
	};

	//Make an instance of the format class
	GAMESSUKOutputFormat theGAMESSUKOutputFormat;

	
	
  bool GAMESSUKOutputFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv) {

	/*
	 
	Read in coordinates after any reorientation due to the symmetry code
	- if there is a field called "input z-matrix" then we read the initial zmatrix in here
	  This geometry is not needed if we are optimising as we can use the optimised geometry
	  However, we need to keep this geometry as it's the only one we have if it's not an opt.
	  
	- if there is no zmat, we need to read the "molecular geometry" field. This geometry
	  only needs to be used if we are not in an optimisation.
	  
	 Read the RUN TYPE field to work out whether we need to use the first geometry.
	 
	 If it's a single point caculation, we can return the molecule at this point
	 
	 If it's some form of structure search, we need to go and find the final structure

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

	vector<string> tokens, geomList; // list of lines and list of tokens on a line
	string line; // For convenience so we can refer to lines from the iterator as 'line'
	//ReadMode_t ReadMode=SKIP;
	
	enum RunType_t { UNKNOWN, SINGLEPOINT, OPTXYZ, OPTZMAT, SADDLE };
	RunType_t RunType=UNKNOWN;
	bool ok;
	
	while (ifs.good() && ifs.getline(buffer, BUFF_SIZE)) {
		
		/* The zmatrix entered by the user
		 * REM:  need to add stuff for "automatic z-matrix generation" as we currently
		 * ignore the zmatrix & just read the cartesian coordinates 
		 */
		
		if (strstr(buffer,"                              input z-matrix") != NULL){
			
			// Set Runtype to SINGLEPOINT so we don't read in the cartesians
			RunType=SINGLEPOINT;
			
			geomList.clear();
			
			// skip 2 lines
			ifs.getline(buffer, BUFF_SIZE) && ifs.getline(buffer, BUFF_SIZE);
			
			// Stick a header line first
			geomList.push_back("zmatrix bohr");
			
			// Read zmatrix into list until blank line 
			while (ifs.good() && ifs.getline(buffer, BUFF_SIZE) && strlen(buffer) != 0){
				line = buffer;
				// transform(method.begin(), method.end(), method.begin(), ::tolower);
				ToLower(line);
				Trim(line);
				geomList.push_back(line);
			}
			
			// Skip 3 lines
			ifs.getline(buffer, BUFF_SIZE) && 
			ifs.getline(buffer, BUFF_SIZE) && 
			ifs.getline(buffer, BUFF_SIZE);
			
			// Read in the variables till we hit blank line
			if (! ReadVariables(ifs, BOHR_TO_ANGSTROM, "")) return false;
			
			// Now go and process the geometry
			ok = ReadGeometry(mol, geomList);
					
		} // End Reading user z-matrix
		
		// Read the cartesian coordinates if we've not read in the ZMATRIX
		if (strstr(buffer,"*            charge       x             y              z       shells") != NULL &&
				RunType==UNKNOWN){
			
			// Skip 3 lines
			ifs.getline(buffer, BUFF_SIZE) && 
			ifs.getline(buffer, BUFF_SIZE) && 
			ifs.getline(buffer, BUFF_SIZE);
			
			// Create regex for the coords
			regex_t *myregex = new regex_t;
			int iok;
			iok = regcomp( myregex,
					//     ------label--------   -------charge-------- < seems enough for a match
					" *\\* *[a-zA-Z]{1,2}[0-9]* *[0-9]{1,3}\\.[0-9]{1}",
					REG_EXTENDED | REG_NOSUB);
			if (iok !=0) cerr << "Error compiling regex in GUK OUTPUT!\n";
			
			// Read in the coordinates - we process them directly rather 
			// then use ReadGeometry as we probably should do...
			mol.BeginModify();
			while (ifs.good() && ifs.getline(buffer, BUFF_SIZE)){
				
				// End of geometry block
				if (strstr(buffer,"*************************")!=NULL)break;
				
				if (regexec( myregex, buffer, 0, 0, 0)==0) {
					//cerr << "Got Coord line: " << buffer << endl;
					OBAtom *atom = mol.NewAtom();
					tokenize(tokens,buffer," ");
					atom->SetAtomicNum(atoi(tokens[2].c_str()));
					double x=atof(tokens[3].c_str())*BOHR_TO_ANGSTROM;
					double y=atof(tokens[4].c_str())*BOHR_TO_ANGSTROM;
					double z=atof(tokens[5].c_str())*BOHR_TO_ANGSTROM;
					atom->SetVector(x, y, z);
				}
			}
			mol.EndModify();			
			regfree(myregex);
			
		} // End Read Cartesian Coords
		
		
		// Determine the RunType - affects how we move on from here.
		if (strstr(buffer," * RUN TYPE") != NULL){
			tokenize(tokens,buffer," \t\n");
			
			if(tokens[3].compare(0,6,"optxyz")==0){
				//cerr << "runtype is optxyz\n";
				RunType=OPTXYZ;
				break;
			} else if (tokens[3].compare(0,8,"optimize")==0){
				//cerr << "runtype is optimise\n";
				RunType=OPTZMAT;
				break;
			} else if (tokens[3].compare(0,6,"saddle")==0){
				//cerr << "runtype is optimise\n";
				RunType=SADDLE;
				break;
			} else {
				RunType=SINGLEPOINT;
				break;
			}
		}
	} // End First Reading loop
	
	
	if(RunType==SINGLEPOINT){
		// We can return the molecule that we've read in
		if (mol.NumAtoms() == 0) { // e.g., if we're at the end of a file PR#1737209
			mol.EndModify();
			return false;
		} else {
			return true;
		}
	}
	
	
	// Clear the Molecule as we're going to start from scratch again.
	mol.BeginModify();
	mol.Clear();
	mol.EndModify();
	
	// Start trundling through the file again - just get the last geometry
	while (ifs.good() && ifs.getline(buffer, BUFF_SIZE)) {
		if (strstr(buffer,"optimization converged") != NULL)
		{	
			if (RunType==OPTXYZ){
				//cerr << "Got converged for OPTXYZ\n";
				
				// FF to start of coordinate specification
				while (ifs.good() && ifs.getline(buffer, BUFF_SIZE)) {
					if (strstr(buffer,
							"atom     znuc       x             y             z") != NULL) break;
				}
				
				// Skip 3 lines - should then be at the coordinates
				ifs.getline(buffer, BUFF_SIZE) && 
				ifs.getline(buffer, BUFF_SIZE) && 
				ifs.getline(buffer, BUFF_SIZE);
				
				// Read in the coordinates - we process them directly rather 
				// then use ReadGeometry as we probably should do...
				mol.BeginModify();
				while (ifs.good() && ifs.getline(buffer, BUFF_SIZE)){
					
					// End of geometry block
					if (strstr(buffer,"*************************")!=NULL)break;
					
					//cerr << "Got Coord line: " << buffer << endl;
					OBAtom *atom = mol.NewAtom();
					tokenize(tokens,buffer," ");
					atom->SetAtomicNum(atoi(tokens[2].c_str()));
					double x=atof(tokens[3].c_str())*BOHR_TO_ANGSTROM;
					double y=atof(tokens[4].c_str())*BOHR_TO_ANGSTROM;
					double z=atof(tokens[5].c_str())*BOHR_TO_ANGSTROM;
					atom->SetVector(x, y, z);
				}
				
				mol.EndModify();
				return true;
				
				
			} else if (RunType==OPTZMAT || RunType==SADDLE) {
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
				
			}
		}
		
	} // End Second Reading loop		

	return true;

} // End GAMESSUKOutputFormat::ReadMolecule


} //namespace OpenBabel
