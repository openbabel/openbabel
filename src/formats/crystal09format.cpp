/**********************************************************************
Copyright (C) 2012 Barry Moore <moore0557@gmail.com>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/elements.h>
#include <openbabel/generic.h>

#include <openbabel/math/vector3.h>

using namespace std;
bool isParentheses (char c)
{
  switch(c)
  {
    case '(':
    case ')':
      return true;
    default:
      return false;
  }
}

namespace OpenBabel
{

class Crystal09Format : public OBMoleculeFormat
{
public:

  //Register this format type ID in the constructor
  Crystal09Format()
  {
    OBConversion::RegisterFormat("c09out",this);
    OBConversion::RegisterOptionParam("b", this, 0, OBConversion::INOPTIONS);
    OBConversion::RegisterOptionParam("s", this, 0, OBConversion::INOPTIONS);
  }

  virtual const char* Description() //required
  {
    return
        "Crystal 09 output format\n"

        "Read Options e.g. -as \n"
        "  s  Consider single bonds only\n"
        "  b  Disable bonding entirely\n" ;
  }

  //Optional URL where the file format is specified
  virtual const char* SpecificationURL()
  {
    return "http://www.crystal.unito.it/";
  }

  /* Flags() can return be any of the following combined by |
     or be omitted if none apply
     NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY  DEFAULTFORMAT
     READBINARY  WRITEBINARY  READXML  ZEROATOMSOK*/
  virtual unsigned int Flags()
  {
    return READONEONLY | NOTWRITABLE;
  }

  /// Declarations for the "API" interface functions
  virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
};

//Make an instance of the format class
Crystal09Format theCrystal09Format;

////////////////-> READ FUNCTIONALITY <-///////////////////////////

bool Crystal09Format::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{
  OBMol* pmol = pOb->CastAndClear<OBMol>();
  istream& ifs = *pConv->GetInStream();
  pmol->BeginModify();
  pmol->SetDimension(3);

  int numAtoms = 0;
  int checkAtoms = 0;
  string line;
  vector <string> vs;
  vector <vector3> one;
  vector <vector3> two;
  vector <vector3> three;
  vector <vector3> four;
  vector <vector3> five;
  vector <vector3> six;
  vector <vector <vector3> > displacements;
  vector <double> xtmp;
  vector <double> ytmp;
  vector <double> ztmp;
  vector <double> freq;
  vector <double> intensity;
  int Iter = 0;
  int extraIter = 0;



  while(!getline(ifs,line).eof()){
    // IF Statement to find number of atoms in asymmetric unit
    if ( line.find("ATOMS IN THE ASYMMETRIC UNIT") != string::npos &&
         checkAtoms == 0 ){
      vector<string> vs;
      tokenize(vs,line);
      int numTokens = vs.size();
      checkAtoms = atoi(vs[(numTokens-1)].c_str());
    }
    // IF statement to input cartesian coordinates of primitive cell
    if ( line.find("CARTESIAN COORDINATES - PRIMITIVE CELL") != string::npos &&
         numAtoms == 0){
      double x,y,z;
      // Skip Three Lines after match
      getline(ifs,line);
      getline(ifs,line);
      getline(ifs,line);

      while (getline(ifs,line)){
        tokenize(vs,line);

        if ( vs.size() < 6 ){
          if ( vs.size() > 0 ){
            //Implies input is missing for some reason
            cerr << "Error with line: " << line << endl;
            cerr << "Structure should be: AtomNumber AtomicNumber Element "
                    "XCoord YCoord ZCoord" << endl;
            break; //Missing Input
          }
          else{
            //Implies end of block i.e. a blank line, checked against
            //checkAtoms for equal numbers
            break;
          }
        }
        else {
          OBAtom * atom = pmol->NewAtom();
          atom->SetAtomicNum(atoi(vs[1].c_str()));
          x = strtod ((char*)vs[3].c_str(), NULL);
          y = strtod ((char*)vs[4].c_str(), NULL);
          z = strtod ((char*)vs[5].c_str(), NULL);
          atom->SetVector(x,y,z);
          numAtoms++;
        }
      }
    }
    //////-> Unit Cell Vector Parser-<///////////
    if (line.find("DIRECT LATTICE VECTORS CARTESIAN COMPONENTS") != string::npos){
      vector3 xvec,yvec,zvec; //vector3 classes to handle lattice vectors
      getline(ifs,line); //skip one line
      
      //First Line to Parse use SetX
      getline(ifs,line);
      tokenize(vs,line);
      xvec.SetX( strtod ((char*)vs[0].c_str(), NULL) );
      yvec.SetX( strtod ((char*)vs[1].c_str(), NULL) );
      zvec.SetX( strtod ((char*)vs[2].c_str(), NULL) );

      //Second Line to Parse use SetY
      getline(ifs,line);
      tokenize(vs,line);
      xvec.SetY( strtod ((char*)vs[0].c_str(), NULL) );
      yvec.SetY( strtod ((char*)vs[1].c_str(), NULL) );
      zvec.SetY( strtod ((char*)vs[2].c_str(), NULL) );

      //Third Line to Parse use SetZ
      getline(ifs,line);
      tokenize(vs,line);
      xvec.SetZ( strtod ((char*)vs[0].c_str(), NULL) );
      yvec.SetZ( strtod ((char*)vs[1].c_str(), NULL) );
      zvec.SetZ( strtod ((char*)vs[2].c_str(), NULL) );

      //Declare a pointer for the UnitCell Data. Set the Unit Cell for
      //OBUnitCell and OBMol
      OBUnitCell *cell = new OBUnitCell;
      cell->SetData(xvec,yvec,zvec);
      pmol->SetData(cell);

    }
    ///////////////-> VIBRATIONAL FREQUENCIES HERE <-/////////////////
    if ( line.find("MODES         EIGV          FREQUENCIES") != string::npos) {

      double tmp1,tmp2;

      getline(ifs,line);//Skip First line

      while (getline(ifs,line)) {
        // Strip out parenthesis, as these may or may not have whitespace
        // around them and can make the tokenize output unpredictable:
        line.erase(remove_if(line.begin(), line.end(), isParentheses), line.end());
        tokenize(vs,line);

        if ( vs.size() < 11 ){
          if ( vs.size() > 0 ){
            //Implies input is missing for some reason
            cerr << "Error with line: " << line << endl;
            break; //Missing Input
          }
          else{
            break; //Implies end of block i.e. a blank line
          }
        }
        else {
          tmp1 = strtod((char*)vs[3].c_str(), NULL);
          tmp2 = strtod((char*)vs[7].c_str(), NULL);
          freq.push_back(tmp1);
          intensity.push_back(tmp2);
        }
      }
      //while Still in this loop we need to skip to the displacement vectors
      int numFreq = freq.size();
      if(numFreq > 0){
        Iter = numFreq / 6;
        extraIter = numFreq % 6;
      }
      else{
        cerr << "Couldn't Parse Frequencies, Check Input" << endl;
        break;
      }
      //Skip two lines to get to displacements
      getline(ifs,line);
      getline(ifs,line);

      for(int i=0; i<Iter; i++){
        getline(ifs,line);
        getline(ifs,line);
        for(int j=0; j<numAtoms ;j++){
          //First we need to construct the vector3 for each frequency and
          //store vector of vector3 in displacements
          getline(ifs,line);
          tokenize(vs,line);
          for(int l=0; l<6; l++){
            xtmp.push_back(strtod((char*)vs[l+4].c_str(), NULL));
          }

          getline(ifs,line);
          tokenize(vs,line);
          for(int l=0; l<6; l++){
            ytmp.push_back(strtod((char*)vs[l+1].c_str(), NULL));
          }

          getline(ifs,line);
          tokenize(vs,line);
          for(int l=0; l<6; l++){
            ztmp.push_back(strtod((char*)vs[l+1].c_str(), NULL));
          }

          one.push_back(vector3(xtmp[0],ytmp[0],ztmp[0]));
          two.push_back(vector3(xtmp[1],ytmp[1],ztmp[1]));
          three.push_back(vector3(xtmp[2],ytmp[2],ztmp[2]));
          four.push_back(vector3(xtmp[3],ytmp[3],ztmp[3]));
          five.push_back(vector3(xtmp[4],ytmp[4],ztmp[4]));
          six.push_back(vector3(xtmp[5],ytmp[5],ztmp[5]));

          xtmp.clear();
          ytmp.clear();
          ztmp.clear();
        }

        displacements.push_back(one);
        displacements.push_back(two);
        displacements.push_back(three);
        displacements.push_back(four);
        displacements.push_back(five);
        displacements.push_back(six);

        one.clear();
        two.clear();
        three.clear();
        four.clear();
        five.clear();
        six.clear();

        getline(ifs,line);
      }

      if(extraIter > 0){
        for(int l=0; l < extraIter; l++){
          displacements.pop_back();
        }
      }

      OBVibrationData* vd = new OBVibrationData;
      vd->SetData(displacements, freq, intensity);
      pmol->SetData(vd);
    }
  }

  //Now We have the frequencies, intensities, and a vector of vector3 (Linear
  //Blocks of numAtoms) and we need to set the vibration data for pmol
  if ( checkAtoms != numAtoms ){
    cerr << "Number of Atoms Specified in Input Does Not Equal Number of "
            "Atoms Read From File!" << endl;
    pmol->EndModify();
    return false;
  }
  else{
    //Input Options
    if( !pConv->IsOption( "b", OBConversion::INOPTIONS ) )
      pmol->ConnectTheDots();

    if( !pConv->IsOption( "s", OBConversion::INOPTIONS) &&
        !pConv->IsOption( "b", OBConversion::INOPTIONS ) ) {
      pmol->PerceiveBondOrders();
    }

    pmol->EndModify();
    return true;
  }
}
}


