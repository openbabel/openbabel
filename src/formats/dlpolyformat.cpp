/**********************************************************************
Copyright  (C) 2010 by Jens Thomas
 
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

#include <iomanip>

namespace OpenBabel
{

class DlpolyConfigFormat : public OBMoleculeFormat
{
public:
  //Register this format type ID
  DlpolyConfigFormat()
  {
    OBConversion::RegisterFormat("CONFIG",this);
  }
  
  virtual const char* Description() //required
  {
    return "DL-POLY CONFIG\n";
  };
  
  virtual const char* SpecificationURL()
  { 
    return "http://www.cse.scitech.ac.uk/ccg/software/DL_POLY";
  };
  
  //Flags() can return be any the following combined by | or be omitted if none apply
  // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
  virtual unsigned int Flags()
  {
    return WRITEONEONLY;
  };
  
  ////////////////////////////////////////////////////
  /// The "API" interface functions
  virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
  virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);


private:

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

  int LabelToAtomicNumber(std::string label);

  std::stringstream errorMsg;
  
};
  
  
  //Make an instance of the format class
  DlpolyConfigFormat theDlpolyConfigFormat;
  
  /////////////////////////////////////////////////////////////////
  int DlpolyConfigFormat::LabelToAtomicNumber(std::string label)
  {
    /*
     * Given a string with the label for an atom return the atomic number
     * As we are using the GetAtomicNum function case is not important
     */
    
    // See if the first 2 characters give us a valid atomic number
    int Z=etab.GetAtomicNum(label.substr(0,2).c_str());

    // If not try the first one
    if (Z==0) Z=etab.GetAtomicNum(label.substr(0,1).c_str());

    if (Z==0){
      // Houston...
      errorMsg << "LabelToAtomicNumber got bad Label: " << label << std::endl;
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
    }
    return Z;

  } //End LabelToAtomicNumber


  bool DlpolyConfigFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    
    std::vector<std::string> tokens; // list of lines and list of tokens on a line
    std::string line; // For convenience so we can refer to lines from the iterator as 'line'
    
    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;
    
    //Define some references so we can use the old parameter names
    std::istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    //const char* title = pConv->GetTitle();
    
    bool ok;
    std::string title;
    int levcfg,imcon;
    std::string AtomLabel;
    int AtomIndex,AtomicNumber;
    double x,y,z;
  
    char buffer[BUFF_SIZE];
    OBAtom *atom;
  
    // Title line
    if ( ! ifs.getline(buffer,BUFF_SIZE) )
      {
        obErrorLog.ThrowError(__FUNCTION__, "Problem reading title line", obWarning);
        return false;
      }  
    title = buffer;
  
    // levcfg & imcon
    if ( ! ifs.getline(buffer,BUFF_SIZE) )
      {
        line=buffer;
        line="Problem reading levcfg line: " + line;
        obErrorLog.ThrowError(__FUNCTION__, line, obWarning);
        return false;
      }
  
    tokenize(tokens, buffer, " \t\n");
    if ( ! tokens.size() >= 2 || ! ( from_string<int>(levcfg, tokens.at(0), std::dec) && from_string<int>(imcon, tokens.at(1), std::dec) ) )
      {
        line=buffer;
        line="Problem reading levcfg line: " + line;
        obErrorLog.ThrowError(__FUNCTION__, line, obWarning);
        return false;
      }

    // If imcon > 0 then there are 3 lines with the cell vectors - for now we just skip them
    if ( imcon > 0 )
      ifs.getline(buffer,BUFF_SIZE) && ifs.getline(buffer,BUFF_SIZE) && ifs.getline(buffer,BUFF_SIZE);

    // Start of coordinates - just loop through reading in data
    int readindex=0;
    AtomicNumber=-1; // Set to -1 so we can check if we've read one in or need to get it from the label.

    mol.BeginModify();
    while (ifs.getline(buffer,BUFF_SIZE))
      {

        //std::cout << "Got line " << buffer << std::endl;

        // Line with AtomLabel, AtomIndex & AtomicNumber - only AtomName required
        if ( readindex == 0 )
          {

            tokenize(tokens, buffer, " \t\n");
            AtomLabel = tokens.at(0);
            if ( tokens.size() >= 2 )
              {
                ok = from_string<int>(AtomIndex, tokens.at(1), std::dec);
              }
            if ( tokens.size() == 3 )
              {
                ok = from_string<int>(AtomicNumber, tokens.at(2), std::dec);
              }
          
            // Got data so read in  coordinates on next line
            readindex=1;
            continue;
          
          } // readindex == 0

        // Coordinate line
        if ( readindex == 1 )
          {
            tokenize(tokens, buffer, " \t\n");
            if ( tokens.size() != 3 || 
                 ! ( from_string<double>(x, tokens.at(0), std::dec) 
                     && from_string<double>(y, tokens.at(1), std::dec) 
                     && from_string<double>(z, tokens.at(2), std::dec) ) )
              {
                line=buffer;
                line="Problem reading coordinate line: " + line;
                obErrorLog.ThrowError(__FUNCTION__, line, obWarning);
                return false;
              }
          
            // Check if we've read in the atomic charge or need to try and extract it from the label
            if ( AtomicNumber == -1 ) AtomicNumber = LabelToAtomicNumber( AtomLabel );

            //std::cout << "Got charge " << AtomicNumber << std::endl;

            atom = mol.NewAtom();
            atom->SetAtomicNum(AtomicNumber);
            atom->SetVector(x, y, z); //set coordinates

            // Reset Atomic Number
            AtomicNumber=-1;
          
            // levcfg == 0, next line will be atom info
            if (levcfg == 0 )
              readindex=0;
            else
              readindex=2;

            continue;
          
          } // readindex == 1

        // Skip velocities line for now
        if ( readindex == 2 )
          {
            if (levcfg == 1 )
              readindex=0;
            else
              readindex=3;

            continue;
          }

        // Skip forces line for now
        if ( readindex == 3 )
          {
            readindex=0;
            continue;
          }

        // Should never get here  
        std::cerr << "MUMMY!\n";

      } // End while reading loop

    mol.EndModify();
    mol.SetTitle(title);
  
    return(true);

  } // End ReadMolecule

  
  bool DlpolyConfigFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {

    /**
     * Write a DL-POLY CONFIG file. Ints are 10 chars wide, floats are formatted as 20.15f
     */

    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;
    
    //Define some references so we can use the old parameter names
    std::ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;
    
    char buffer[BUFF_SIZE];

    // We only print the coordinates and labels
    int levcfg=0, imcon=0;
    
    int idx=0; // For counting molecule index
    std::string title = mol.GetTitle();

    // 80 char title
    ofs << title.substr(0,80) << std::endl;
    
    // Set levcfg & imcon
    ofs << std::setw(10) << levcfg << std::setw(10) << imcon << std::endl;
    
    // Now loop through molecule
    FOR_ATOMS_OF_MOL(atom, mol)
      {
        
        ofs << std::setw(8) << etab.GetSymbol(atom->GetAtomicNum()) << std::setw(10) << ++idx << std::setw(10) << atom->GetAtomicNum() << std::endl;
        snprintf(buffer, BUFF_SIZE, "%20.15f %20.15f %20.15f\n",
                 atom->GetX(),
                 atom->GetY(),
                 atom->GetZ()
                 );
        ofs << buffer;
      }

    return(true);
    
  } //End WriteMolecule
  
  
} //namespace OpenBabel
