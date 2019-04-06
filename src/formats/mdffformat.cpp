/**********************************************************************
Copyright (C) 2004 by Chris Morley for template
Copyright (C) 2009 by David C. Lonie for VASP template
Copyright (C) 2014 by Kirill Okhotnikov for MDFF format 

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
#include <openbabel/elements.h>
#include <openbabel/generic.h>
#include <openbabel/obiter.h>


#include <limits.h>
#include <locale> // For isalpha(int)
#include <map>
#include <stdexcept>
#include <cstdlib>
#include <algorithm>

#ifdef _MSC_VER
#define INFINITY (DBL_MAX+DBL_MAX)
#define NAN (INFINITY-INFINITY)
#endif

using namespace std;
namespace OpenBabel {
  class MDFFFormat : public OBMoleculeFormat
  {
  public:

    MDFFFormat()
    {
      OBConversion::RegisterFormat("POSFF",this);
      OBConversion::RegisterFormat("CONTFF",this);
      OBConversion::RegisterFormat("MDFF",this);      
    }

    virtual const char* Description()
    {
      return
        "MDFF format\n"
        "The format used in the POSFF and CONTFF files used by MDFF\n\n"

        "POSFF and CONTFF are read to obtain information from MDFF calculations.\n"
        "The program will try to read the IONS.POT file if the name of the\n"
        "input file is POSFF or CONTFF.\n"

        "Write Options e.g. -xw\n"
        "  w Sort atoms by atomic number\n"
        "  u <elementlist> Sort atoms by list of element symbols provided in comma-separated string w/o spaces\n"
        "  i Write IONS.POT file\n"              
        ;

    };

    virtual const char* SpecificationURL(){return "https://code.google.com/p/mdff/";};

    /* Flags() can return be any of the following combined by |
       or be omitted if none apply
       NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY  DEFAULTFORMAT
       READBINARY  WRITEBINARY  READXML  ZEROATOMSOK */
    virtual unsigned int Flags()
    {
      return READONEONLY;
    };

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
  MDFFFormat theMDFFFormat;

  /////////////////////////////////////////////////////////////////
  
  struct atm_t_prop
  {
    int num_of_atoms;
    int atom_etab_num;
    string atom_symbol;
    double atom_charge;
    atm_t_prop(): num_of_atoms(0), atom_symbol(""), atom_charge(0.0) {};
  };
  

  bool MDFFFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    // Move stream to EOF, some apps check ifs position to check for multimolecule files.
    // MDFF does not support this, and this parser makes its own streams.
    istream &ifs = *pConv->GetInStream();
    ifs.seekg (0, ios::end);

    char buffer[BUFF_SIZE];
    double x,y,z;
    unsigned int totalAtoms_fl = 0, totalTypes = 0;
    unsigned int totalAtoms = 0, atomCount = 0;
    OBAtom *atom;
    bool cartesian;
    vector<string> vs;
    vector<atm_t_prop> atom_t_prop;

    // Get path of input file:
    //    ifs_path.getline(buffer,BUFF_SIZE);
    //    path = buffer;
    string full_path = pConv->GetInFilename();
    size_t found = full_path.rfind("/");
    string path = (found == string::npos) ? "" : path.substr(0, found);
    string short_fn = full_path.substr(path.length(), string::npos);

    // Open files
    string posff_filename = pConv->GetInFilename(); 
    ifstream ifs_posff (posff_filename.c_str());
    if (!ifs_posff) {
      return false; // No geometry file?
    }

    bool process_ions = (short_fn == "CONTFF") || (short_fn == "POSFF");
    string ionspot_filename = (path == "") ? "IONS.POT" : path + "/IONS.POT";    
    ifstream ifs_ions;
    if( process_ions )
    {  
      ifs_ions.open(ionspot_filename.c_str());
      process_ions = ifs_ions.is_open();
    }  
    
    pmol->BeginModify();

    // Start working on POSFF:
    ifs_posff.getline(buffer,BUFF_SIZE); // Total Number of atoms
    totalAtoms_fl = atoi(buffer);
    ifs_posff.getline(buffer,BUFF_SIZE); // Comment line
    pmol->SetTitle(buffer);

    ifs_posff.getline(buffer,BUFF_SIZE); // X_Vec vector
    tokenize(vs, buffer);
    x = atof(vs.at(0).c_str());
    y = atof(vs.at(1).c_str());
    z = atof(vs.at(2).c_str());
    vector3 x_vec (x,y,z);

    ifs_posff.getline(buffer,BUFF_SIZE); // Y_Vec vector
    tokenize(vs, buffer);
    x = atof(vs.at(0).c_str());
    y = atof(vs.at(1).c_str());
    z = atof(vs.at(2).c_str());
    vector3 y_vec (x,y,z);

    ifs_posff.getline(buffer,BUFF_SIZE); // Z_Vec vector
    tokenize(vs, buffer);
    x = atof(vs.at(0).c_str());
    y = atof(vs.at(1).c_str());
    z = atof(vs.at(2).c_str());
    vector3 z_vec (x,y,z);

    // Build unit cell
    OBUnitCell *cell = new OBUnitCell;
    cell->SetData(x_vec, y_vec, z_vec);
    cell->SetSpaceGroup(1);
    pmol->SetData(cell);
    
    // get number of different atoms types.
    ifs_posff.getline(buffer,BUFF_SIZE);
    totalTypes = atoi(buffer);
    atom_t_prop.resize(totalTypes);

    ifs_posff.getline(buffer, BUFF_SIZE);
    tokenize(vs, buffer);
    
    if(vs.size() != atom_t_prop.size() )
    {
      obErrorLog.ThrowError(__FUNCTION__, "Number of types is wrong. Format is not read.", obError);
      pmol->EndModify();
      return false;
    }  
    for (size_t i = 0; i < atom_t_prop.size(); ++i) 
    {  
      atom_t_prop[i].atom_symbol = vs[i];
      atom_t_prop[i].atom_etab_num = OpenBabel::OBElements::GetAtomicNum(atom_t_prop[i].atom_symbol.c_str());
    }  
 
    // Fetch next line to get stoichiometry
    ifs_posff.getline(buffer,BUFF_SIZE);
    tokenize(vs, buffer);
    
    if(vs.size() != atom_t_prop.size() )
    {
      obErrorLog.ThrowError(__FUNCTION__, "Number of types is wrong. Format is not readed.", obError);
      pmol->EndModify();
      return false;
    }  
    
    // Extract and sum the atom counts. The sum is used to parse the atomic
    // coordinates
    totalAtoms = 0;
    for (size_t i = 0; i < atom_t_prop.size(); ++i) 
    {  
      int currentCount = atoi(vs.at(i).c_str());
      atom_t_prop[i].num_of_atoms =  currentCount;
      totalAtoms += currentCount;
    }
    
    //Read Ions properties
    if( process_ions )
    {
      while(!ifs_ions.eof())
      {
        vector<string> vs;
        ifs_ions.getline(buffer, BUFF_SIZE); 
        tokenize(vs, buffer);
        vs.erase(find(vs.begin(), vs.end(), "!") ,vs.end());
        vs.erase(remove(vs.begin(), vs.end(), "="));
                
        if(vs.size() == 0)
          continue;
        if(vs[0] == "qch")
        {
          if(vs.size() - 1 != atom_t_prop.size())
            obErrorLog.ThrowError(__FUNCTION__, "Number of charges for atom is wrong. Skipped", obWarning);
          else
            for(int i = 0; i < atom_t_prop.size(); i++)
            {
              string tk = vs[i + 1];
              string str = tk.substr(tk.length() - 2) == "d0" ? tk.substr(0, tk.length() - 2) : tk;
              atom_t_prop[i].atom_charge = atof(tk.c_str());
            }  
        } 
      }
    }  

    // Cartesian or fractional?
    ifs_posff.getline(buffer,BUFF_SIZE);
    // [C|c|K|k] indicates cartesian coordinates, anything else (including
    // an empty string, buffer[0] == 0) indicates fractional coordinates
    cartesian = ( buffer[0] == 'C' || buffer[0] == 'c' || 
                  buffer[0] == 'K' || buffer[0] == 'k' );

    atomCount = 0;
    for(unsigned int i = 0; i < atom_t_prop.size(); i++)
    {
      bool err_break = false;
      for(unsigned int j = 0; j < atom_t_prop[i].num_of_atoms; j++)
      {
        ifs_posff.getline(buffer, BUFF_SIZE); // atom location
        // Parse the buffer now.
        tokenize(vs, buffer);
        atom = pmol->NewAtom();
        atom->SetAtomicNum(atom_t_prop[i].atom_etab_num);
        
        err_break = ( vs.size() < 4 ) || ( vs[0] != atom_t_prop[i].atom_symbol );
        if (err_break)
          break;

        x = atof(vs[1].c_str());
        y = atof(vs[2].c_str());
        z = atof(vs[3].c_str());
        vector3 coords (x,y,z);
        if (!cartesian)
          coords = cell->FractionalToCartesian( coords );
        atom->SetVector(coords);
        atom->SetFormalCharge(atom_t_prop[i].atom_charge);
        atomCount++;
      }
      if (err_break)
        break;
    }
    
    if ( atomCount != totalAtoms_fl )
    {
      stringstream errorMsg;      
      errorMsg << "Problems reading a MDFF file: "
               << "The number of atoms read is not right."
               <<  "atomCount = " << atomCount << " totalAtoms_fl = " << totalAtoms_fl << endl;
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
      
      pmol->EndModify();
      return false;
    }

    // There is some trailing garbage, but AFAIK it's not useful for anything.
    ifs_posff.close();
    pmol->EndModify();
    
    return true;
  }
  
  class aindx
  {
  public:
    int index_param;
    //int etab_number;
    int atom_index;
    bool operator<(const aindx &left) const
    {
      vector<int> vs;
      vs.push_back(index_param - left.index_param);
      //vs.push_back(etab_number - left.etab_number);
      vs.push_back(atom_index  - left.atom_index);
      
      bool less = false;
      for(int i = 0; i < vs.size(); i++)      
      {
        if( vs[i] != 0)
        {  
          less = vs[i] < 0;
          break;
        }  
      }
      return less;  
    }
  };

  bool MDFFFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    //No surprises in this routine, cartesian coordinates are written out
    //and if at least a single atom has information about constraints,
    //then selective dynamics is used and the info is written out.
    //The atoms are ordered according to their atomic number so that the
    //output looks nice, this can be reversed by using command line flag "-xw".
    //
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if (pmol == NULL) {
      return false;
    }

    ostream& ofs = *pConv->GetOutStream();
    OBMol mol(*pmol);
    
    if(mol.HasData(OBGenericDataType::UnitCell))
    {
      OBUnitCell *uc = static_cast<OBUnitCell*>(mol.GetData(OBGenericDataType::UnitCell));
      uc->FillUnitCell(&mol);
    }            

    char buffer[BUFF_SIZE];
    OBUnitCell *uc = NULL;
    vector<vector3> cell;

    const char * sortAtoms     = pConv->IsOption("w", OBConversion::OUTOPTIONS);
    const char * sortAtomsList = pConv->IsOption("u", OBConversion::OUTOPTIONS);
    const char * writeIONS     = pConv->IsOption("i", OBConversion::OUTOPTIONS);    

    // Create a list of ids. These may be sorted by atomic number depending
    // on the value of keepOrder.
    
    map<int, int> indl; 
    
    if (sortAtoms != NULL) 
    {
      indl.clear();
      for(int i = 0; i < 200; i++)
        indl[i] = i;
    }
    
    if (sortAtomsList != NULL) 
    {
      indl.clear();
      vector<string> vs;
      tokenize(vs, sortAtomsList);
      for(int i = 0; i < vs.size(); i++)
        indl[OBElements::GetAtomicNum(vs[i].c_str())] = i;
    }
    
    map<aindx, OBAtom *> amap;
    
    FOR_ATOMS_OF_MOL(atom, mol) 
    {
      aindx ndx;
      ndx.index_param = indl[atom->GetAtomicNum()];
      ndx.atom_index = atom->GetIndex();
      amap[ndx] = &(*atom);
    }
    
    //Set elements array
    vector< pair<string, unsigned int> > atypes_def;
    string last_atom_smb = "";
    for(map<aindx, OBAtom *>::const_iterator it = amap.begin(); it != amap.end(); ++it)
    {
      string curr_atom_smb = OpenBabel::OBElements::GetSymbol(it->second->GetAtomicNum());
      if( last_atom_smb != curr_atom_smb )
      {  
        last_atom_smb = curr_atom_smb;
        atypes_def.push_back( pair<string, unsigned int>(curr_atom_smb, 0) );
      }
      atypes_def[atypes_def.size() - 1].second++;
    }

    // write number of atoms
    ofs << mol.NumAtoms() << endl;
    // write title
    ofs << mol.GetTitle() << endl;
    if (!mol.HasData(OBGenericDataType::UnitCell)) {
      // the unit cell has not been defined. Leave as all zeros so the user
      // can fill it in themselves
      for (int ii = 0; ii < 3; ii++) {
        snprintf(buffer, BUFF_SIZE, "0.0  0.0  0.0");
        ofs << buffer << endl;
      }
    }
    else
    {
      // there is a unit cell, write it out
      uc = static_cast<OBUnitCell*>(mol.GetData(OBGenericDataType::UnitCell));
      cell = uc->GetCellVectors();
      for (vector<vector3>::const_iterator i = cell.begin();
           i != cell.end(); ++i) {
        snprintf(buffer, BUFF_SIZE, "%20.15f%20.15f%20.15f",
                 i->x(), i->y(), i->z());
        ofs << buffer << endl;
      }
    }
    //Print the number of atoms types
    ofs <<  atypes_def.size() << endl;
    
    for (int i = 0; i < atypes_def.size(); i++)
    {
      snprintf(buffer, BUFF_SIZE, "%-3s ", atypes_def[i].first.c_str());
      ofs << buffer ;
    }
    ofs << endl;
    
    for (int i = 0; i < atypes_def.size(); i++)
    {
      snprintf(buffer, BUFF_SIZE, "%-3u ", atypes_def[i].second);
      ofs << buffer ;
    }
    ofs << endl;

    // and test if any of the atoms has constraints
    // print the atomic coordinates in \AA
    ofs << "Cartesian" << endl;
    
    map<string, double> charge_smb;

    for (map<aindx, OBAtom *>::const_iterator it  = amap.begin(); 
                                              it != amap.end(); ++it)
    {  
      // Print coordinates
      string smb = OpenBabel::OBElements::GetSymbol(it->second->GetAtomicNum());
      snprintf(buffer, BUFF_SIZE, "%-3s %26.19f %26.19f %26.19f", smb.c_str(),
               it->second->GetX(), it->second->GetY(), it->second->GetZ());
      
      if(charge_smb.find(smb) == charge_smb.end() )
        charge_smb[smb] = it->second->GetFormalCharge();
      else
        if(charge_smb[smb] != it->second->GetFormalCharge())
          charge_smb[smb] = NAN;
      
      ofs << buffer << endl;
    }
    
    if( writeIONS == NULL)
      return true;
    
    //Write IONS.POT
    string path = pConv->GetOutFilename();
    size_t found = path.rfind("/");
    path = (found == string::npos) ? "" : path.substr(0, found);
    string ionspot_filename = (path == "") ? "IONS.POT" : path + "/IONS.POT";
    ofstream ofs_ions;
    ofs_ions.open(ionspot_filename.c_str(), fstream::out);
    
    for(int i = 0; i < 4; i++)
    {  
      if( (i == 0) || (i == 2) )
        ofs_ions << " ! " ;
      else if (i == 1)
        ofs_ions << " mass = " ;
      else if (i == 3)
        ofs_ions << " qch = " ;
        
      for(int j = 0; j < atypes_def.size(); j++)
      {
        if( (i == 0) || (i == 2) )
          ofs_ions << atypes_def[j].first << "  ";
        else if (i == 1)
          ofs_ions << OBElements::GetMass(OBElements::GetAtomicNum(atypes_def[j].first.c_str())) << "d0 ";
        else if (i == 3)
          ofs_ions << charge_smb[atypes_def[j].first] << "d0 ";
      }
      ofs_ions << endl;
    }

    return true;
  }

} //namespace OpenBabel


