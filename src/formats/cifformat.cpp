/**********************************************************************
cifformat.cpp - Implementation of subclass of OBFormat for conversion of OBMol.

Copyright (C) 2006 Vincent Favre-Nicolin

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

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/elements.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/generic.h>
#include <cstdlib>

#include <openbabel/math/spacegroup.h>

#include <sstream>
#include <vector>
#include <list>
#include <map>
#include <set>

#define NOCHARGE FLT_MAX

#ifdef _MSC_VER
 #pragma warning( disable : 4503 )
 // The decorated name was longer than the compiler limit (4096), and was truncated.
 // This is due to the use of templates specialized on templates repeatedly.
 // The correctness of the program, however, is unaffected by the truncated name,
 // but if you get link time errors on a truncated symbol, it will be more difficult
 // to determine the type of the symbol in the error. Debugging will also be more difficult;
 // the debugger will also have difficultly mapping symbol name to type name.
#endif

using namespace std;

namespace OpenBabel
{
  class CIFFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    CIFFormat()
    {
      RegisterFormat("cif", "chemical/x-cif");
    }

    virtual const char* Description() //required
    {
      return
        "Crystallographic Information File\n"
        "The CIF file format is the standard interchange format for small-molecule crystal structures\n\n"
        "Fractional coordinates are converted to cartesian ones using the following convention:\n\n"

        "- The x axis is parallel to a\n"
        "- The y axis is in the (a,b) plane\n"
        "- The z axis is along c*\n\n"

        "Ref: Int. Tables for Crystallography (2006), vol. B, sec 3.3.1.1.1\n"
        "  (the matrix used is the 2nd form listed)\n\n"

        "Read Options e.g. -ab:\n"
        "  s  Output single bonds only\n"
        "  b  Disable bonding entirely\n"
        "  B  Use bonds listed in CIF file from _geom_bond_etc records (overrides option b) \n\n";
    };

    virtual const char* SpecificationURL()
    {return "http://www.iucr.org/iucr-top/cif/spec/";}; //optional

    //*** This section identical for most OBMol conversions ***
    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);
  };
  //############################## Case-insensituve string####################################################
  // :@todo: This duplicates normal case-insensitive string comparison in OpenBabel
  /// Case-insensitive string class
  /// From: Guru of the Week #29
  /// e.g.: http://gcc.gnu.org/onlinedocs/libstdc++/21_strings/gotw29a.txt
  ///
  /// Public domain
  struct ci_char_traits : public std::char_traits<char>
  {
    static bool eq( char c1, char c2 );

    static bool ne( char c1, char c2 );

    static bool lt( char c1, char c2 );

    static int compare(const char* s1,const char* s2,size_t n );

    static const char* find( const char* s, int n, char a );
  };

  typedef std::basic_string<char, ci_char_traits> ci_string;
  int strnicmp(const char *s1, const char *s2, int len)
  {
    unsigned char c1, c2;
    while (len)
      {
        c1 = *s1; c2 = *s2;
        s1++; s2++;
        if (!c1) return c2 ? -1 : 0;
        if (!c2) return 1;
        if (c1 != c2)
          {
            c1 = tolower(c1);
            c2 = tolower(c2);
            if (c1 != c2) return c1 < c2 ? -1 : 1;
          }
        len--;
      }
    return 0;
  }

  bool ci_char_traits::eq( char c1, char c2 )
  {return tolower(c1) == tolower(c2);}

  bool ci_char_traits::ne( char c1, char c2 )
  {return tolower(c1) != tolower(c2);}

  bool ci_char_traits::lt( char c1, char c2 )
  {return tolower(c1) < tolower(c2);}

  int ci_char_traits::compare(const char* s1,const char* s2,size_t n )
  {return strnicmp( s1, s2, n );}

  const char* ci_char_traits::find( const char* s, int n, char a )
  {
    while( n-- > 0 && tolower(*s) != tolower(a) ) ++s;
    return s;
  }
  //############################## CIF CLASSES headers####################################################
  /** The CIFData class holds all the information from a \e single data_ block from a cif file.
   *
   * It is a placeholder for all comments, item and loop data, as raw strings copied from
   * a cif file.
   *
   * It is also used to interpret this data to extract parts of the cif data, i.e.
   * only part of the core cif dictionnary are recognized. CIF tags currently recognized
   * include ("tag1 > tag2" means tag1 is preferred to tag2 when extracting the info, only one is reported):
   *  - crystal name: _chemical_name_systematic > _chemical_name_mineral > _chemical_name_structure_type > _chemical_name_common
   *  - crystal formula: _chemical_formula_analytical > _chemical_formula_structural > _chemical_formula_iupac > _chemical_formula_moiety
   *  - unit cell:  _cell_length_{a,b,c} ; _cell_angle_{alpha,beta,gamma}
   *  - spacegroup number: _space_group_IT_number > _symmetry_Int_Tables_number
   *  - spacegroup Hall symbol: _space_group_name_Hall > _symmetry_space_group_name_Hall
   *  - spacegroup Hermann-Mauguin symbol:_space_group_name_H-M_alt > _symmetry_space_group_name_H-M
   *  - atom coordinates: _atom_site_fract_{x} ; _atom_site_Cartn_{x,y,z}
   *  - atom occupancy: _atom_site_occupancy
   *  - atom label & symbol: _atom_site_type_symbol ; _atom_site_label
   *
   * Cartesian coordinates are stored in Angstroems, angles in radians.
   *
   * If another data field is needed, it is possible to directly access the string data
   * (CIFData::mvComment , CIFData::mvItem and CIFData::mvLoop) to search for the correct tags.
   */
  class CIFData
  {
  public:
    CIFData();

    /// Extract lattice parameters, spacegroup (symbol or number), atomic positions,
    /// chemical name and formula if available.
    /// All other data is ignored
    void ExtractAll();
    /// Extract name & formula for the crystal
    void ExtractName();
    /// Extract unit cell
    void ExtractUnitCell();
    /// Extract spacegroup number or symbol
    void ExtractSpacegroup();
    /// Extract all atomic positions. Will generate cartesian from fractional
    /// coordinates or vice-versa if only cartesian coordinates are available.
    void ExtractAtomicPositions();
    /// Extract listed bond distances, from _geom_bond_* loops
    void ExtractBonds();
    //// Extract Charges information from cif file and assign it to atoms
    void ExtractCharges();
    /// Generate fractional coordinates from cartesian ones for all atoms
    /// CIFData::CalcMatrices() must be called first
    void Cartesian2FractionalCoord();
    /// Generate cartesian coordinates from fractional ones for all atoms
    /// CIFData::CalcMatrices() must be called first
    void Fractional2CartesianCoord();
    /// Convert from fractional to cartesian coordinates
    /// CIFData::CalcMatrices() must be called first
    void f2c(float &x,float &y, float &z);
    /// Convert from cartesia to fractional coordinates
    /// CIFData::CalcMatrices() must be called first
    void c2f(float &x,float &y, float &z);
    /// Calculate real space transformation matrices
    /// requires unit cell parameters
    void CalcMatrices();
    /// Comments from CIF file, in the order they were read
    std::list<std::string> mvComment;
    /// Individual CIF items
    std::map<ci_string,std::string> mvItem;
    /// CIF Loop data
    std::map<std::set<ci_string>,std::map<ci_string,std::vector<std::string> > > mvLoop;
    /// Lattice parameters, in ansgtroem and degrees - vector size is 0 if no
    /// parameters have been obtained yet.
    std::vector<float> mvLatticePar;
    /// Spacegroup number from International Tables (_space_group_IT_number), or -1.
    unsigned int mSpacegroupNumberIT;
    /// Spacegroup Hall symbol (or empty string) (_space_group_name_Hall)
    std::string mSpacegroupSymbolHall;
    /// Spacegroup Hermann-Mauguin symbol (or empty string) (_space_group_name_H-M_alt)
    std::string mSpacegroupHermannMauguin;
    /// Crystal name. Or empty string if none is available.
    std::string mName;
    /// Formula. Or empty string if none is available.
    std::string mFormula;
    /// Atom record
    struct CIFAtom
    {
      CIFAtom();
      /// Label of the atom, or empty string (_atom_site_label).
      std::string mLabel;
      /// Symbol of the atom, or empty string (_atom_type_symbol or _atom_site_type_symbol).
      std::string mSymbol;
      /// Fractionnal coordinates (_atom_site_fract_{x,y,z}) or empty vector.
      std::vector<float> mCoordFrac;
      /// Cartesian coordinates in Angstroem (_atom_site_Cartn_{x,y,z}) or empty vector.
      /// Transformation to fractionnal coordinates currently assumes
      /// "a parallel to x; b in the plane of y and z" (see _atom_sites_Cartn_transform_axes)
      std::vector<float> mCoordCart;
      /// Site occupancy, or -1
      float mOccupancy;
      //charge from oxydation
      float mCharge;
    };
    /// Atoms, if any are found
    std::vector<CIFAtom> mvAtom;
    /// Bond distance record
    struct CIFBond
    {
      /// Label of the two bonded atoms
      std::string mLabel1,mLabel2;
      /// distance
      float mDistance;
    };
    /// Atoms, if any are found
    std::vector<CIFBond> mvBond;
    /// Fractionnal2Cartesian matrix
    float mOrthMatrix[3][3];
    /// Cartesian2Fractionnal matrix
    float mOrthMatrixInvert[3][3];
    const SpaceGroup *mSpaceGroup;
    /// Name for the CIF data block
    std::string mDataBlockName;
  };
  /** Main CIF class - parses the stream and separates data blocks, comments, items, loops.
   * All values are stored as string, and Each CIF block is stored in a separate CIFData object.
   * No interpretaion is made here - this must be done from all CIFData objects.
   */
  class CIF
  {
  public:
    /// Creates the CIF object from a stream
    ///
    /// \param interpret: if true, interpret all data blocks. See CIFData::ExtractAll()
    CIF(std::istream &in, const bool interpret=true);
    //private:
    /// Separate the file in data blocks and parse them to sort tags, loops and comments.
    /// All is stored in the original strings.
    ///
    /// Returns the name of the next data block
    void Parse(std::istream &in);
    /// The data blocks, after parsing. The key is the name of the data block
    std::map<std::string,CIFData> mvData;
    /// Global comments, outside and data block
    std::list<std::string> mvComment;
  };
  /// Convert one CIF value to a floating-point value
  /// Return 0 if no value can be converted (e.g. if '.' or '?' is encountered)
  float CIFNumeric2Float(const std::string &s);
  /// Convert one CIF value to a floating-point value
  /// Return 0 if no value can be converted (e.g. if '.' or '?' is encountered)
  int CIFNumeric2Int(const std::string &s);

  template <typename T> string to_string(T pNumber)
  {
    ostringstream oOStrStream;
    oOStrStream << pNumber;
    return oOStrStream.str();
  }

  bool is_double(const std::string& s, double& r_double);

  //############################## CIF CLASSES CODE ####################################################
  CIFData::CIFAtom::CIFAtom():
    mLabel(""),mSymbol(""),mOccupancy(1.0f)
  {}

  CIFData::CIFData()
  {}

  void CIFData::ExtractAll()
  {
    {
      stringstream ss;
      ss<<"CIF: interpreting data block: "<<mDataBlockName;
      obErrorLog.ThrowError(__FUNCTION__, ss.str(), obInfo);
    }
    if(mDataBlockName=="data_global")
      { // :KLUDGE: this data block name is used for journal &author information
        // for IUCr journals, so do not generate an error if the block contains
        // no structural information
        bool empty_iucrjournal_block=true;
        if(mvItem.find("_cell_length_a"   )!=mvItem.end()) empty_iucrjournal_block=false;
        if(mvItem.find("_cell_length_b"   )!=mvItem.end()) empty_iucrjournal_block=false;
        if(mvItem.find("_cell_length_c"   )!=mvItem.end()) empty_iucrjournal_block=false;
        for(map<set<ci_string>,map<ci_string,vector<string> > >::const_iterator loop=mvLoop.begin();
            loop!=mvLoop.end();++loop)
          {
            if(loop->second.find("_atom_site_fract_x")!=loop->second.end()) empty_iucrjournal_block=false;
            if(loop->second.find("_atom_site_fract_y")!=loop->second.end()) empty_iucrjournal_block=false;
            if(loop->second.find("_atom_site_fract_z")!=loop->second.end()) empty_iucrjournal_block=false;
            if(loop->second.find("_atom_site_Cartn_x")!=loop->second.end()) empty_iucrjournal_block=false;
            if(loop->second.find("_atom_site_Cartn_y")!=loop->second.end()) empty_iucrjournal_block=false;
            if(loop->second.find("_atom_site_Cartn_z")!=loop->second.end()) empty_iucrjournal_block=false;
          }
        if(empty_iucrjournal_block)
          {
            stringstream ss;
            ss << "CIF WARNING: found en empty 'data_global' block - SKIPPING\n"
               << "  (you can safely ignore this if reading a CIF file from an IUCr journal)";
            obErrorLog.ThrowError(__FUNCTION__, ss.str(), obWarning);
            return;
          }
    }
    // :@todo: Take care of values listed as "." and "?" instead of a real value.
    this->ExtractName();
    // Spacegroup must be extracted before unit cell
    this->ExtractSpacegroup();
    this->ExtractUnitCell();
    this->ExtractAtomicPositions();
    if(mvAtom.size()==0)
      {
        stringstream ss;
        ss << "CIF Error: no atom found ! (in data block:"<<mDataBlockName<<")";
        obErrorLog.ThrowError(__FUNCTION__, ss.str(), obError);
      }
    this->ExtractBonds();
    this->ExtractCharges();
  }

  void CIFData::ExtractUnitCell()
  {
    // Use spacegroup to determine missing angle or length
    const int spgid= mSpaceGroup->GetId();
    if(  (mvItem.find("_cell_length_a")!=mvItem.end())
       ||(mvItem.find("_cell_length_b")!=mvItem.end())
       ||(mvItem.find("_cell_length_c")!=mvItem.end()) )
      {
        mvLatticePar.resize(6);
        for(unsigned int i=0;i<6;i++) mvLatticePar[i]=float(0);
        map<ci_string,string>::const_iterator positem;
        positem=mvItem.find("_cell_length_a");
        if(positem!=mvItem.end())
          mvLatticePar[0]=CIFNumeric2Float(positem->second);
        positem=mvItem.find("_cell_length_b");
        if(positem!=mvItem.end())
          mvLatticePar[1]=CIFNumeric2Float(positem->second);
        positem=mvItem.find("_cell_length_c");
        if(positem!=mvItem.end())
          mvLatticePar[2]=CIFNumeric2Float(positem->second);
        positem=mvItem.find("_cell_angle_alpha");
        if(positem!=mvItem.end())
          mvLatticePar[3]=CIFNumeric2Float(positem->second);
        positem=mvItem.find("_cell_angle_beta");
        if(positem!=mvItem.end())
          mvLatticePar[4]=CIFNumeric2Float(positem->second);
        positem=mvItem.find("_cell_angle_gamma");
        if(positem!=mvItem.end())
          mvLatticePar[5]=CIFNumeric2Float(positem->second);
        stringstream ss;
        ss << "Found Lattice parameters:" << mvLatticePar[0] << " , " << mvLatticePar[1] << " , " << mvLatticePar[2]
           << " , " << mvLatticePar[3] << " , " << mvLatticePar[4] << " , " << mvLatticePar[5];
        obErrorLog.ThrowError(__FUNCTION__, ss.str(), obDebug);
        mvLatticePar[3] = static_cast<float> (mvLatticePar[3] * DEG_TO_RAD);// pi/180
        mvLatticePar[4] = static_cast<float> (mvLatticePar[4] * DEG_TO_RAD);
        mvLatticePar[5] = static_cast<float> (mvLatticePar[5] * DEG_TO_RAD);

        // Fill values depending on spacegroup, *only* when missing
        if((spgid>2)&&(spgid<=15))
        {// :TODO: monoclinic spg, depending on unique axis....
        }
        if((spgid>15)&&(spgid<=142))
        {// orthorombic & tetragonal
          if(mvLatticePar[3]==0) mvLatticePar[3]=M_PI/2;
          if(mvLatticePar[4]==0) mvLatticePar[4]=M_PI/2;
          if(mvLatticePar[5]==0) mvLatticePar[5]=M_PI/2;
        }
        if((spgid>74)&&(spgid<=142))
        {// Tetragonal, make sure a=b if one is missing
          if(mvLatticePar[1]==0) mvLatticePar[1]=mvLatticePar[0];
          if(mvLatticePar[0]==0) mvLatticePar[0]=mvLatticePar[1];
        }
        if((spgid>142)&&(spgid<=194))
        {// trigonal/ rhomboedric...
          const string::size_type pos = mSpaceGroup->GetHallName().find('R');
          if(pos==std::string::npos)
          {//rhomboedric cell, a=b=c, alpha=beta=gamma
            float a=0;
            if(mvLatticePar[0]>a) a=mvLatticePar[0];
            if(mvLatticePar[1]>a) a=mvLatticePar[1];
            if(mvLatticePar[2]>a) a=mvLatticePar[2];
            if(mvLatticePar[0]==0) mvLatticePar[0]=a;
            if(mvLatticePar[1]==0) mvLatticePar[1]=a;
            if(mvLatticePar[2]==0) mvLatticePar[2]=a;

            float alpha=0;
            if(mvLatticePar[3]>alpha) alpha=mvLatticePar[3];
            if(mvLatticePar[4]>alpha) alpha=mvLatticePar[4];
            if(mvLatticePar[5]>alpha) alpha=mvLatticePar[5];
            if(mvLatticePar[3]==0) mvLatticePar[3]=alpha;
            if(mvLatticePar[4]==0) mvLatticePar[4]=alpha;
            if(mvLatticePar[5]==0) mvLatticePar[5]=alpha;
          }
          else
          {//hexagonal cell, a=b & alpha=beta=pi/2, gamma= 2*pi/3
            if(mvLatticePar[1]==0) mvLatticePar[1]=mvLatticePar[0];
            if(mvLatticePar[0]==0) mvLatticePar[0]=mvLatticePar[1];
            if(mvLatticePar[3]==0) mvLatticePar[3]=M_PI/2;
            if(mvLatticePar[4]==0) mvLatticePar[4]=M_PI/2;
            if(mvLatticePar[5]==0) mvLatticePar[5]=2*M_PI/3;
          }
        }
        if(spgid>194)
        {
          if(mvLatticePar[3]==0) mvLatticePar[3]=M_PI/2;
          if(mvLatticePar[4]==0) mvLatticePar[4]=M_PI/2;
          if(mvLatticePar[5]==0) mvLatticePar[5]=M_PI/2;
          // In case some idiot cif only supplies one value, make sure a=b=c
          float a=0;
          if(mvLatticePar[0]>a) a=mvLatticePar[0];
          if(mvLatticePar[1]>a) a=mvLatticePar[1];
          if(mvLatticePar[2]>a) a=mvLatticePar[2];
          if(mvLatticePar[0]==0) mvLatticePar[0]=a;
          if(mvLatticePar[1]==0) mvLatticePar[1]=a;
          if(mvLatticePar[2]==0) mvLatticePar[2]=a;
        }
        // Handle missing values
        if(mvLatticePar[3]<1e-6)
        {
            stringstream ss;
            ss << "CIF WARNING: missing alpha value, defaulting to 90 degrees (in data block:"<<mDataBlockName<<")";
            obErrorLog.ThrowError(__FUNCTION__, ss.str(), obWarning);
            mvLatticePar[3]=90*DEG_TO_RAD;
        }
        if(mvLatticePar[4]<1e-6)
        {
            stringstream ss;
            ss << "CIF WARNING: missing beta value, defaulting to 90 degrees (in data block:"<<mDataBlockName<<")";
            obErrorLog.ThrowError(__FUNCTION__, ss.str(), obWarning);
            mvLatticePar[4]=90*DEG_TO_RAD;
        }
        if(mvLatticePar[5]<1e-6)
        {
            stringstream ss;
            ss << "CIF WARNING: missing gamma value, defaulting to 90 degrees (in data block:"<<mDataBlockName<<")";
            obErrorLog.ThrowError(__FUNCTION__, ss.str(), obWarning);
            mvLatticePar[5]=90*DEG_TO_RAD;
        }
        if(mvLatticePar[1]<1e-6)
        {
            stringstream ss;
            ss << "CIF Error: missing b lattice parameter - cannot interpret structure ! (in data block:"<<mDataBlockName<<")";
            obErrorLog.ThrowError(__FUNCTION__, ss.str(), obError);
        }
        if(mvLatticePar[2]<1e-6)
        {
            stringstream ss;
            ss << "CIF Error: missing c lattice parameter - cannot interpret structure ! (in data block:"<<mDataBlockName<<")";
            obErrorLog.ThrowError(__FUNCTION__, ss.str(), obError);
        }

        this->CalcMatrices();
      }
      else
      {
         stringstream ss;
         ss << "CIF Error: missing a,b and c value - cannot interpret structure ! (in data block:"<<mDataBlockName<<")";
         obErrorLog.ThrowError(__FUNCTION__, ss.str(), obError);
      }
  }

  void CIFData::ExtractSpacegroup()
  {
    map<ci_string,string>::const_iterator positem;
    bool found = false;
    positem=mvItem.find("_space_group_IT_number");
    if(positem!=mvItem.end())
      {
        mSpacegroupNumberIT=CIFNumeric2Int(positem->second);
        found = true;
        stringstream ss;
        ss << "Found spacegroup IT number:" << mSpacegroupNumberIT;
        obErrorLog.ThrowError(__FUNCTION__, ss.str(), obDebug);
      }
    else
      {
        positem=mvItem.find("_symmetry_Int_Tables_number");
        if(positem!=mvItem.end())
          {
            mSpacegroupNumberIT=CIFNumeric2Int(positem->second);
            found = true;
            stringstream ss;
            ss << "Found spacegroup IT number (with OBSOLETE CIF #1.0 TAG):" << mSpacegroupNumberIT;
            obErrorLog.ThrowError(__FUNCTION__, ss.str(), obDebug);
          }
        else {
          positem=mvItem.find("_symmetry_group_IT_number");
          if(positem!=mvItem.end())
          {
            mSpacegroupNumberIT=CIFNumeric2Int(positem->second);
            found = true;
            stringstream ss;
            ss << "Found spacegroup IT number (with NON-STANDARD CIF TAG):" << mSpacegroupNumberIT;
            obErrorLog.ThrowError(__FUNCTION__, ss.str(), obDebug);
          }
          else
            mSpacegroupNumberIT=0;
        }
      }

    positem=mvItem.find("_space_group_name_Hall");
    if(positem!=mvItem.end())
      {
        mSpacegroupSymbolHall=positem->second;
        found = true;
        obErrorLog.ThrowError(__FUNCTION__, "Found spacegroup Hall symbol:"+mSpacegroupSymbolHall, obDebug);
      }
    else
      {
        positem=mvItem.find("_symmetry_space_group_name_Hall");
        if(positem!=mvItem.end())
          {
            mSpacegroupSymbolHall=positem->second;
            found = true;
            obErrorLog.ThrowError(__FUNCTION__, "Found spacegroup Hall symbol (with OBSOLETE CIF #1.0 TAG):"+mSpacegroupSymbolHall, obDebug);
          }
      }

    positem=mvItem.find("_space_group_name_H-M_alt");
    if(positem!=mvItem.end())
      {
        mSpacegroupHermannMauguin=positem->second;
        found = true;
        obErrorLog.ThrowError(__FUNCTION__, "Found spacegroup Hermann-Mauguin symbol:"+mSpacegroupHermannMauguin, obDebug);
      }
    else
      {
        positem=mvItem.find("_symmetry_space_group_name_H-M");
        if(positem!=mvItem.end())
          {
            mSpacegroupHermannMauguin=positem->second;
            found = true;
            obErrorLog.ThrowError(__FUNCTION__, "Found spacegroup Hermann-Mauguin symbol (with OBSOLETE CIF #1.0 TAG):"+mSpacegroupHermannMauguin, obDebug);
          }
      }
    // DDL2 tag is "_space_group.IT_coordinate_system_code", converted by the cif reader to "_space_group_IT_coordinate_system_code"
    positem=mvItem.find("_space_group_IT_coordinate_system_code");
    if(positem!=mvItem.end())
      {
        obErrorLog.ThrowError(__FUNCTION__, "Found spacegroup IT_coordinate_system_code:"+positem->second, obDebug);
        if((mSpacegroupHermannMauguin.length()>0) && (positem->second=="1" || positem->second=="2"))
        {
          // this is a HACK which will work as long as the HM symbols in spacegroups.txt have the ":1" or ":2" extension listed, when needed
          mSpacegroupHermannMauguin=mSpacegroupHermannMauguin+string(":")+positem->second;
        }
        else
        {
          stringstream ss;
          ss << "CIF Error: found DDL2 tag _space_group.IT_coordinate_system_code ("<<positem->second<<")"<<endl
             <<"            but could not interpret it ! Origin choice or axis may be incorrect.";
          obErrorLog.ThrowError(__FUNCTION__, ss.str(), obWarning);
        }
      }

    mSpaceGroup=NULL;
    // be forgiving - if spg not found, try again
    // Prefer Hall > HM == number, as Hall symbol is truly unique
    if (mSpacegroupSymbolHall.length() > 0) {
      //Make sure there are no leading spaces before Hall symbol (kludge)
      for(std::string::iterator pos=mSpacegroupSymbolHall.begin();pos!=mSpacegroupSymbolHall.end();)
      {
        if((char)(*pos)==' ')  pos=mSpacegroupSymbolHall.erase(pos);
        else ++pos;
      }
      mSpaceGroup = SpaceGroup::GetSpaceGroup(mSpacegroupSymbolHall);
    }
    if((mSpaceGroup == NULL)&& (mSpacegroupHermannMauguin.length() > 0)) {
      mSpaceGroup = SpaceGroup::GetSpaceGroup(mSpacegroupHermannMauguin);
    }
    if((mSpaceGroup == NULL)&&(mSpacegroupNumberIT != 0)) {
      mSpaceGroup = SpaceGroup::GetSpaceGroup(mSpacegroupNumberIT);
    }
    if(mSpaceGroup == NULL) {
      SpaceGroup *sg = new SpaceGroup();
      positem=mvItem.find("_space_group_symop_operation_xyz");
      if(positem==mvItem.end())
        positem=mvItem.find("_symmetry_equiv_pos_as_xyz");
      if(positem!=mvItem.end())
        {
          sg->AddTransform (positem->second);
          found = true;
        }
      else {
        for(map<set<ci_string>,map<ci_string,vector<string> > >::const_iterator loop=mvLoop.begin();
            loop!=mvLoop.end();++loop)
          {
            map<ci_string,vector<string> >::const_iterator pos;
            unsigned i, nb;
            pos=loop->second.find("_space_group_symop_operation_xyz");
            if (pos==loop->second.end())
              pos=loop->second.find("_symmetry_equiv_pos_as_xyz");
            if (pos!=loop->second.end())
              {
                nb=pos->second.size();
                found = true;
                for (i = 0; i < nb; i++)
                  sg->AddTransform(pos->second[i]);
                break; // found the transforms, so we have done with them
              }
          }
        if (found)
          mSpaceGroup = SpaceGroup::Find(sg);
        if (mSpaceGroup == NULL && sg->IsValid())
          mSpaceGroup = sg;
        else
          delete sg;
      }
    }
    if(mSpaceGroup == NULL)
    {
        stringstream ss;
        ss << "CIF Error: missing spacegroup description: defaulting to P1... (in data block:"<<mDataBlockName<<")";
        obErrorLog.ThrowError(__FUNCTION__, ss.str(), obWarning);
        mSpaceGroup = SpaceGroup::GetSpaceGroup(1);
    }
    // set the space group name to Hall symbol
    mSpacegroupSymbolHall = mSpaceGroup->GetHallName();
  }

  void CIFData::ExtractName()
  {
    map<ci_string,string>::const_iterator positem;
    positem=mvItem.find("_chemical_name_systematic");
    if(positem!=mvItem.end())
      {
        mName=positem->second;
        obErrorLog.ThrowError(__FUNCTION__, "Found chemical name:"+mName, obDebug);
      }
    else
      {
        positem=mvItem.find("_chemical_name_mineral");
        if(positem!=mvItem.end())
          {
            mName=positem->second;
            obErrorLog.ThrowError(__FUNCTION__, "Found chemical name:"+mName, obDebug);
          }
        else
          {
            positem=mvItem.find("_chemical_name_structure_type");
            if(positem!=mvItem.end())
              {
                mName=positem->second;
                obErrorLog.ThrowError(__FUNCTION__, "Found chemical name:"+mName, obDebug);
              }
            else
              {
                positem=mvItem.find("_chemical_name_common");
                if(positem!=mvItem.end())
                  {
                    mName=positem->second;
                    obErrorLog.ThrowError(__FUNCTION__, "Found chemical name:"+mName, obDebug);
                  }
              }
          }
      }
    /// Crystal formula
    positem=mvItem.find("_chemical_formula_analytical");
    if(positem!=mvItem.end())
      {
        mFormula=positem->second;
        obErrorLog.ThrowError(__FUNCTION__, "Found chemical formula:"+mFormula, obDebug);
      }
    else
      {
        positem=mvItem.find("_chemical_formula_structural");
        if(positem!=mvItem.end())
          {
            mFormula=positem->second;
            obErrorLog.ThrowError(__FUNCTION__, "Found chemical formula:"+mFormula, obDebug);
          }
        else
          {
            positem=mvItem.find("_chemical_formula_iupac");
            if(positem!=mvItem.end())
              {
                mFormula=positem->second;
                obErrorLog.ThrowError(__FUNCTION__, "Found chemical formula:"+mFormula, obDebug);
              }
            else
              {
                positem=mvItem.find("_chemical_formula_moiety");
                if(positem!=mvItem.end())
                  {
                    mFormula=positem->second;
                    obErrorLog.ThrowError(__FUNCTION__, "Found chemical formula:"+mFormula, obDebug);
                  }
              }
          }
      }
  }

  void CIFData::ExtractAtomicPositions()
  {
    map<ci_string,string>::const_iterator positem;
    for(map<set<ci_string>,map<ci_string,vector<string> > >::const_iterator loop=mvLoop.begin();
        loop!=mvLoop.end();++loop)
      {
        if(mvAtom.size()>0) break;// only extract ONE list of atoms, preferably fractional coordinates
        map<ci_string,vector<string> >::const_iterator posx,posy,posz,poslabel,possymbol,posoccup;
        posx=loop->second.find("_atom_site_fract_x");
        posy=loop->second.find("_atom_site_fract_y");
        posz=loop->second.find("_atom_site_fract_z");
        unsigned int nb = 0;
        if( (posx!=loop->second.end()) && (posy!=loop->second.end()) && (posz!=loop->second.end()))
          {
            nb=posx->second.size();
            mvAtom.resize(nb);
            for(unsigned int i=0;i<nb;++i)
              {
                mvAtom[i].mCoordFrac.resize(3);
                mvAtom[i].mCoordFrac[0]=CIFNumeric2Float(posx->second[i]);
                mvAtom[i].mCoordFrac[1]=CIFNumeric2Float(posy->second[i]);
                mvAtom[i].mCoordFrac[2]=CIFNumeric2Float(posz->second[i]);
              }
            this->Fractional2CartesianCoord();
          }
        else
          {
            posx=loop->second.find("_atom_site_Cartn_x");
            posy=loop->second.find("_atom_site_Cartn_y");
            posz=loop->second.find("_atom_site_Cartn_z");
            if( (posx!=loop->second.end()) && (posy!=loop->second.end()) && (posz!=loop->second.end()))
              {
                nb=posx->second.size();
                mvAtom.resize(nb);
                for(unsigned int i=0;i<nb;++i)
                  {
                    mvAtom[i].mCoordCart.resize(3);
                    mvAtom[i].mCoordCart[0]=CIFNumeric2Float(posx->second[i]);
                    mvAtom[i].mCoordCart[1]=CIFNumeric2Float(posy->second[i]);
                    mvAtom[i].mCoordCart[2]=CIFNumeric2Float(posz->second[i]);
                  }
                this->Cartesian2FractionalCoord();
              }
          }
        if(mvAtom.size()>0)
          {// Got the atoms, get names and symbols
            possymbol=loop->second.find("_atom_site_type_symbol");
            if(possymbol!=loop->second.end())
              for(unsigned int i=0;i<nb;++i)
                mvAtom[i].mSymbol=possymbol->second[i];
            poslabel=loop->second.find("_atom_site_label");
            if(poslabel!=loop->second.end())
              for(unsigned int i=0;i<nb;++i)
                {
                  mvAtom[i].mLabel=poslabel->second[i];
                  if(possymbol==loop->second.end())
                    {// There was no symbol, use the labels to guess it
                      int nbc=0;
                      if(mvAtom[i].mLabel.size()==1)
                        if(isalpha(mvAtom[i].mLabel[0])) nbc=1;
                      if(mvAtom[i].mLabel.size()>=2)
                        {
                          if(isalpha(mvAtom[i].mLabel[0]) && isalpha(mvAtom[i].mLabel[1])) nbc=2;
                          else if(isalpha(mvAtom[i].mLabel[0])) nbc=1;
                        }
                      if(nbc>0) mvAtom[i].mSymbol=mvAtom[i].mLabel.substr(0,nbc);
                      else mvAtom[i].mSymbol="H";//Something wen wrong, no symbol !
                    }
                }
            // Occupancy ?
            posoccup=loop->second.find("_atom_site_occupancy");
            if(posoccup!=loop->second.end())
              for(unsigned int i=0;i<nb;++i)
                {
                  mvAtom[i].mOccupancy=CIFNumeric2Float(posoccup->second[i]);
                  if( (mvAtom[i].mOccupancy <= 0.0) || (mvAtom[i].mOccupancy > 1.0) )
                    mvAtom[i].mOccupancy = 1.0;
                }
            // Now be somewhat verbose
            stringstream ss;
            ss << "Found "<<nb<<" atoms."<<endl;
            for(unsigned int i=0;i<nb;++i)
              {
                ss<<mvAtom[i].mLabel<<" "<<mvAtom[i].mSymbol;
                if(mvAtom[i].mCoordFrac.size()>0)
                  {
                    ss<<" , Fractional: ";
                    for(unsigned int j=0;j<mvAtom[i].mCoordFrac.size();++j)
                      ss<<mvAtom[i].mCoordFrac[j]<<" ";
                  }
                if(mvAtom[i].mCoordCart.size()>0)
                  {
                    ss<<" , Cartesian: ";
                    for(unsigned int j=0;j<mvAtom[i].mCoordCart.size();++j)
                      ss<<mvAtom[i].mCoordCart[j]<<" ";
                  }
                ss<<" , Occupancy= "<<mvAtom[i].mOccupancy<<endl;
              }
            obErrorLog.ThrowError(__FUNCTION__, ss.str(), obDebug);
          }
      }
  }

  void CIFData::ExtractBonds()
  {
    map<ci_string,string>::const_iterator positem;
    for(map<set<ci_string>,map<ci_string,vector<string> > >::const_iterator loop=mvLoop.begin(); loop!=mvLoop.end();++loop)
      {
        //if(mvBond.size()>0) break;// Only allow one bond list
        map<ci_string,vector<string> >::const_iterator poslabel1,poslabel2,posdist;
        poslabel1=loop->second.find("_geom_bond_atom_site_label_1");
        poslabel2=loop->second.find("_geom_bond_atom_site_label_2");
        posdist=loop->second.find("_geom_bond_distance");
        if( (poslabel1!=loop->second.end()) && (poslabel2!=loop->second.end()) && (posdist!=loop->second.end()))
          {
            obErrorLog.ThrowError(__FUNCTION__, "Found _geom_bond* record...", obDebug);
            const unsigned long nb=poslabel1->second.size();
            mvBond.resize(nb);
            for(unsigned int i=0;i<nb;++i)
              {
                mvBond[i].mLabel1=poslabel1->second[i];
                mvBond[i].mLabel2=poslabel2->second[i];
                mvBond[i].mDistance=CIFNumeric2Float(posdist->second[i]);
                stringstream ss;
                ss << "  d(" << mvBond[i].mLabel1 << "-" << mvBond[i].mLabel2 << ")=" << mvBond[i].mDistance;
                obErrorLog.ThrowError(__FUNCTION__, ss.str(), obDebug);
              }
          }
      }
  }

  void CIFData::ExtractCharges()
  {
    map<ci_string,string>::const_iterator positem;

    map<std::string, double> lbl2ox;
    for(map<set<ci_string>, map<ci_string, vector<string> > >::const_iterator loop=mvLoop.begin(); loop!=mvLoop.end(); ++loop)
    {
      //if(mvBond.size()>0) break;// Only allow one bond list
      map<ci_string,vector<string> >::const_iterator pos_symbol, pos_ox_number, posdist;
      pos_symbol    =loop->second.find("_atom_type_symbol");
      pos_ox_number =loop->second.find("_atom_type_oxidation_number");
      if( (pos_symbol != loop->second.end()) && (pos_ox_number != loop->second.end()) )
      {
        obErrorLog.ThrowError(__FUNCTION__, " Found _atom_type* record with oxydation number...", obDebug);
        const unsigned long nl = pos_symbol->second.size();

        for(unsigned int i = 0; i < nl; i++)
        {
          lbl2ox[pos_symbol->second[i]] = CIFNumeric2Float(pos_ox_number->second[i]);
          obErrorLog.ThrowError(__FUNCTION__, " has oxydation "+pos_ox_number->second[i], obDebug);
        }
      }
    }

    for (std::vector<CIFAtom>::iterator it = mvAtom.begin() ; it != mvAtom.end(); ++it)
    {
      string label = (*it).mLabel;

      if( lbl2ox.count(label) > 0 )
        (*it).mCharge = lbl2ox[label];
      else
      {
        (*it).mCharge = NOCHARGE;
        obErrorLog.ThrowError(__FUNCTION__, "Charge for label: "+label+" cannot be found.", obDebug);
      }
    }
  }

  void CIFData::CalcMatrices()
  {
    if(mvLatticePar.size()==0) return;//:@todo: throw error
    float a,b,c,alpha,beta,gamma;//direct space parameters
    float aa,bb,cc,alphaa,betaa,gammaa;//reciprocal space parameters
    float v;//volume of the unit cell
    a=mvLatticePar[0];
    b=mvLatticePar[1];
    c=mvLatticePar[2];
    alpha=mvLatticePar[3];
    beta=mvLatticePar[4];
    gamma=mvLatticePar[5];

    v=sqrt(1-cos(alpha)*cos(alpha)-cos(beta)*cos(beta)-cos(gamma)*cos(gamma)
           +2*cos(alpha)*cos(beta)*cos(gamma));

    aa=sin(alpha)/a/v;
    bb=sin(beta )/b/v;
    cc=sin(gamma)/c/v;

    alphaa=acos( (cos(beta )*cos(gamma)-cos(alpha))/sin(beta )/sin(gamma) );
    betaa =acos( (cos(alpha)*cos(gamma)-cos(beta ))/sin(alpha)/sin(gamma) );
    gammaa=acos( (cos(alpha)*cos(beta )-cos(gamma))/sin(alpha)/sin(beta ) );

    mOrthMatrix[0][0]=a;
    mOrthMatrix[0][1]=b*cos(gamma);
    mOrthMatrix[0][2]=c*cos(beta);

    mOrthMatrix[1][0]=0;
    mOrthMatrix[1][1]=b*sin(gamma);
    mOrthMatrix[1][2]=-c*sin(beta)*cos(alphaa);

    mOrthMatrix[2][0]=0;
    mOrthMatrix[2][1]=0;
    mOrthMatrix[2][2]=1/cc;

    // Invert upper triangular matrix
    float cm[3][3];
    cm[0][0]=mOrthMatrix[0][0];
    cm[0][1]=mOrthMatrix[0][1];
    cm[0][2]=mOrthMatrix[0][2];

    cm[1][0]=mOrthMatrix[1][0];
    cm[1][1]=mOrthMatrix[1][1];
    cm[1][2]=mOrthMatrix[1][2];

    cm[2][0]=mOrthMatrix[2][0];
    cm[2][1]=mOrthMatrix[2][1];
    cm[2][2]=mOrthMatrix[2][2];
    for(long i=0;i<3;i++)
      for(long j=0;j<3;j++)
        if(i==j) mOrthMatrixInvert[i][j]=1;
        else mOrthMatrixInvert[i][j]=0;
    for(long i=0;i<3;i++)
      {
        float a;
        for(long j=i-1;j>=0;j--)
          {
            a=cm[j][i]/cm[i][i];
            for(long k=0;k<3;k++) mOrthMatrixInvert[j][k] -= mOrthMatrixInvert[i][k]*a;
            for(long k=0;k<3;k++) cm[j][k] -= cm[i][k]*a;
          }
        a=cm[i][i];
        for(long k=0;k<3;k++) mOrthMatrixInvert[i][k] /= a;
        for(long k=0;k<3;k++) cm[i][k] /= a;
      }
      stringstream ss;
      ss <<"Fractional2Cartesian matrix:"<<endl
           <<mOrthMatrix[0][0]<<" "<<mOrthMatrix[0][1]<<" "<<mOrthMatrix[0][2]<<endl
           <<mOrthMatrix[1][0]<<" "<<mOrthMatrix[1][1]<<" "<<mOrthMatrix[1][2]<<endl
           <<mOrthMatrix[2][0]<<" "<<mOrthMatrix[2][1]<<" "<<mOrthMatrix[2][2]<<endl<<endl;
      ss <<"Cartesian2Fractional matrix:"<<endl
           <<mOrthMatrixInvert[0][0]<<" "<<mOrthMatrixInvert[0][1]<<" "<<mOrthMatrixInvert[0][2]<<endl
           <<mOrthMatrixInvert[1][0]<<" "<<mOrthMatrixInvert[1][1]<<" "<<mOrthMatrixInvert[1][2]<<endl
           <<mOrthMatrixInvert[2][0]<<" "<<mOrthMatrixInvert[2][1]<<" "<<mOrthMatrixInvert[2][2];
      obErrorLog.ThrowError(__FUNCTION__, ss.str(), obDebug);
  }

  void CIFData::f2c(float &x,float &y, float &z)
  {
    const float x0=x,y0=y,z0=z;
    x=mOrthMatrix[0][0]*x0+mOrthMatrix[0][1]*y0+mOrthMatrix[0][2]*z0;
    y=mOrthMatrix[1][0]*x0+mOrthMatrix[1][1]*y0+mOrthMatrix[1][2]*z0;
    z=mOrthMatrix[2][0]*x0+mOrthMatrix[2][1]*y0+mOrthMatrix[2][2]*z0;
  }

  void CIFData::c2f(float &x,float &y, float &z)
  {
    const float x0=x,y0=y,z0=z;
    x=mOrthMatrixInvert[0][0]*x0+mOrthMatrixInvert[0][1]*y0+mOrthMatrixInvert[0][2]*z0;
    y=mOrthMatrixInvert[1][0]*x0+mOrthMatrixInvert[1][1]*y0+mOrthMatrixInvert[1][2]*z0;
    z=mOrthMatrixInvert[2][0]*x0+mOrthMatrixInvert[2][1]*y0+mOrthMatrixInvert[2][2]*z0;
  }

  void CIFData::Cartesian2FractionalCoord()
  {
    if(mvLatticePar.size()==0) return;//:@todo: report error
    for(vector<CIFAtom>::iterator pos=mvAtom.begin();pos!=mvAtom.end();++pos)
      {
        pos->mCoordFrac.resize(3);
        pos->mCoordFrac[0]=pos->mCoordCart.at(0);
        pos->mCoordFrac[1]=pos->mCoordCart.at(1);
        pos->mCoordFrac[2]=pos->mCoordCart.at(2);
        c2f(pos->mCoordFrac[0],pos->mCoordFrac[1],pos->mCoordFrac[2]);
      }
  }

  void CIFData::Fractional2CartesianCoord()
  {
    if(mvLatticePar.size()==0) return;//:@todo: report error
    for(vector<CIFAtom>::iterator pos=mvAtom.begin();pos!=mvAtom.end();++pos)
      {
        pos->mCoordCart.resize(3);
        pos->mCoordCart[0]=pos->mCoordFrac.at(0);
        pos->mCoordCart[1]=pos->mCoordFrac.at(1);
        pos->mCoordCart[2]=pos->mCoordFrac.at(2);
        f2c(pos->mCoordCart[0],pos->mCoordCart[1],pos->mCoordCart[2]);
      }
  }

  /////


  CIF::CIF(istream &is, const bool interpret)
  {
    bool found_atoms=false;
    while(!found_atoms)
    {
      // :TODO: we don't need a vector of CIFData, since only one block is read at a time
      mvData.clear();
      this->Parse(is);
      // Extract structure from 1 block
      if(interpret)
        for(map<string,CIFData>::iterator posd=mvData.begin();posd!=mvData.end();++posd)
        {
          posd->second.ExtractAll();
          if(posd->second.mvAtom.size()>0) found_atoms=true;
        }
    }
  }

  bool iseol(const char c) { return ((c=='\n')||(c=='\r'));}

  /// Read one value, whether it is numeric, string or text
  string CIFReadValue(istream &in,char &lastc)
  {
    bool vv=false;//very verbose ?
    string value("");
    while(!isgraph(in.peek())) in.get(lastc);
    while(in.peek()=='#')
      {//discard these comments for now
        string tmp;
        getline(in,tmp);
        lastc='\r';
        while(!isgraph(in.peek())) in.get(lastc);
      }
    if(in.peek()=='_') {
      stringstream errorMsg;
      errorMsg << "Warning: Trying to read a value but found a new CIF tag !";
      obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obError);
      return value;
    }
    if(in.peek()==';')
      {//SemiColonTextField
        bool warning=!iseol(lastc);
        if(warning){
          stringstream errorMsg;
          errorMsg << "Warning: Trying to read a SemiColonTextField but last char is not an end-of-line char !";
          obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obError);
        }
        value="";
        in.get(lastc);
        while(in.peek()!=';')
          {
            if (in.peek() == '_') {
              stringstream errorMsg;
              errorMsg << "Warning: Trying to read a value but found a new CIF tag !";
              obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obError);
              warning = true;
              break;
            }
            string tmp;
            getline(in,tmp);
            value+=tmp+" ";
          }
        if (!warning)
          in.get(lastc);
        if(vv) obErrorLog.ThrowError(__FUNCTION__, "SemiColonTextField:"+value, obDebug);
        if(warning && !vv) obErrorLog.ThrowError(__FUNCTION__, "SemiColonTextField:"+value, obDebug);
        return value;
      }
    if((in.peek()=='\'') || (in.peek()=='\"'))
      {//QuotedString
        char delim;
        in.get(delim);
        value="";
        while(!((lastc==delim)&&(!isgraph(in.peek()))) )
          {
            in.get(lastc);
            value+=lastc;
          }
        if(vv) obErrorLog.ThrowError(__FUNCTION__, "QuotedString:"+value, obDebug);
        return value.substr(0,value.size()-1);
      }
    // If we got here, we have an ordinary value, numeric or unquoted string
    in>>value;
    if(vv) obErrorLog.ThrowError(__FUNCTION__, "NormalValue:"+value, obDebug);
    return value;
  }

  void CIF::Parse(istream &in)
  {
    bool vv=false;//very verbose ?
    char lastc=' ';
    string block="";// Current block data
    while(!in.eof())
      {
        while(!isgraph(in.peek()) && !in.eof()) in.get(lastc);
        if(in.peek()=='#')
          {//Comment
            string tmp;
            getline(in,tmp);
            if(block=="") mvComment.push_back(tmp);
            else mvData[block].mvComment.push_back(tmp);
            lastc='\r';
            continue;
          }
        if(in.peek()=='_')
          {//Tag
            string tag,value;
            in>>tag;
            // Convert all dots to underscores to cover much of DDL2 with this DDL1 parser.
            for (string::size_type pos = tag.find('.'); pos != string::npos; pos = tag.find('.', ++ pos))
              tag.replace(pos, 1, 1, '_');
            value=CIFReadValue(in,lastc);
            mvData[block].mvItem[ci_string(tag.c_str())]=value;
            if(vv)
            {
              stringstream ss;
              ss<<"New Tag:"<<tag<<" ("<<value.size()<<"):"<<value;
              obErrorLog.ThrowError(__FUNCTION__, ss.str(), obDebug);
            }
            continue;
          }
        if((in.peek()=='d') || (in.peek()=='D'))
          {// Data
            if(!mvData.empty()) return; // We want just a single data block

            string tmp;
            in>>tmp;
            block=tmp.substr(5);
            if(vv)
            {
              stringstream ss;
              ss<<endl<<endl<<"NEW BLOCK DATA: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ->"<<block<<endl<<endl;
              obErrorLog.ThrowError(__FUNCTION__, ss.str(), obDebug);
            }
            mvData[block]=CIFData();
            mvData[block].mDataBlockName=tmp;
            continue;
          }
        if((in.peek()=='l') || (in.peek()=='L'))
          {// loop_
            vector<ci_string> tit;
            string tmp;
            in>>tmp; //should be loop_
            if(vv) obErrorLog.ThrowError(__FUNCTION__, "LOOP : "+tmp, obDebug);
            while(true)
              {//read titles
                while(!isgraph(in.peek()) && !in.eof()) in.get(lastc);
                if(in.peek()=='#')
                  {
                    getline(in,tmp);
                    if(block=="") mvComment.push_back(tmp);
                    else mvData[block].mvComment.push_back(tmp);
                    continue;
                  }
                if(in.peek()!='_')
                  {
                    stringstream ss;
                    ss << "End of loop titles:"<<(char)in.peek();
                    if(vv) obErrorLog.ThrowError(__FUNCTION__, ss.str(), obDebug);
                    break;
                  }
                in>>tmp;
                // Convert all dots to underscores to cover much of DDL2 with this DDL1 parser.
                for (string::size_type pos = tmp.find('.'); pos != string::npos; pos = tmp.find('.', ++ pos))
                  tmp.replace(pos, 1, 1, '_');
                tit.push_back(ci_string(tmp.c_str()));
                if(vv) obErrorLog.ThrowError(__FUNCTION__, " , "+tmp, obDebug);
              }
            map<ci_string,vector<string> > lp;
            while(true)
              {
                std::ios::pos_type pos0=in.tellg();
                while(!isgraph(in.peek()) && !in.eof()) in.get(lastc);
                if(in.eof()) break;
                if(in.peek()=='_') break;
                if(in.peek()=='#')
                  {// Comment (in a loop ??)
                    //const std::ios::pos_type pos=in.tellg();
                    string tmp;
                    getline(in,tmp);
                    pos0=in.tellg();
                    if(block=="") mvComment.push_back(tmp);
                    else mvData[block].mvComment.push_back(tmp);
                    lastc='\r';
                    if(vv) obErrorLog.ThrowError(__FUNCTION__, "Comment in a loop (?):"+tmp, obDebug);
                    //in.seekg(pos);
                    break;
                  }
                //in>>tmp;
                tmp=CIFReadValue(in,lastc);
                if(ci_string(tmp.c_str())=="loop_")
                  {//go back and continue
                    in.clear();
                    in.seekg(pos0,std::ios::beg);
                    stringstream ss;
                    ss <<"END OF LOOP :"<<tmp<<","<<(char)in.peek()<<","<<in.tellg();
                    if(vv) obErrorLog.ThrowError(__FUNCTION__, ss.str(), obDebug);
                    break;
                  }
                if(tmp.size()>=5)
                  if(ci_string(tmp.substr(0,5).c_str())=="data_")
                    {//go back and continue
                      in.clear();
                      in.seekg(pos0,std::ios::beg);
                      stringstream ss;
                      ss <<"END OF LOOP :"<<tmp<<","<<(char)in.peek()<<","<<in.tellg();
                      if(vv) obErrorLog.ThrowError(__FUNCTION__, ss.str(), obDebug);
                      break;
                    }
                for(unsigned int i=0;i<tit.size();++i)
                  {//Read all values
                    if(i>0) tmp=CIFReadValue(in,lastc);
                    lp[tit[i]].push_back(tmp);
                    stringstream ss;
                    ss <<" LOOP VALUE    #"<<lp[tit[i]].size()<<","<<i<<" :  "<<tmp;
                    if(vv) obErrorLog.ThrowError(__FUNCTION__, ss.str(), obDebug);
                  }
              }
            // The key to the mvLoop map is the set of column titles
            set<ci_string> stit;
            for(unsigned int i=0;i<tit.size();++i) stit.insert(tit[i]);
            mvData[block].mvLoop[stit]=lp;
            continue;
          }
        // If we get here, something went wrong ! Discard till end of line...
        // It is OK if this is just a blank line though
        string junk;
        getline(in,junk);

        if(junk.size()>0)
        {
          stringstream errorMsg;
          errorMsg << "Warning: one line could not be interpreted while reading a CIF file:"<<endl
                   << " -> line contents:" << junk;
          obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
        }
      }
  }

  float CIFNumeric2Float(const string &s)
  {
    if((s==".") || (s=="?")) return 0.0;
    float v;
    const int n=sscanf(s.c_str(),"%f",&v);
    if(n!=1) return 0.0;
    return v;
  }

  int CIFNumeric2Int(const string &s)
  {
    if((s==".") || (s=="?")) return 0;
    int v;
    const int n=sscanf(s.c_str(),"%d",&v);
    if(n!=1) return 0;
    return v;
  }

  bool is_double(const std::string& s, double& r_double)
  {
    std::istringstream i(s);

    if (i >> r_double)
      return true;

    r_double = 0.0;
    return false;
  }


  //################ END CIF CLASSES######################################

  //Make an instance of the format class
  CIFFormat theCIFFormat;

  // Helper function for CorrectFormatCharges
  // Is this atom an oxygen in a water molecule
  // We know the oxygen is connected to one ion, but check for non-hydrogens
  // Returns: true if the atom is an oxygen and connected to two hydrogens and up to one other atom
  bool CIFisWaterOxygen(OBAtom *atom)
  {
    if (atom->GetAtomicNum() != OBElements::Oxygen)
      return false;

    int nonHydrogenCount = 0;
    int hydrogenCount = 0;
    FOR_NBORS_OF_ATOM(neighbor, *atom) {
      if (neighbor->GetAtomicNum() != OBElements::Hydrogen)
        nonHydrogenCount++;
      else
        hydrogenCount++;
    }

    return (hydrogenCount == 2 && nonHydrogenCount <= 1);
  }

  // Look for lone ions, and correct their formal charges
  void CorrectFormalCharges(OBMol *mol)
  {
    if (!mol)
      return;

    // First look for NR4, PR4 ions,
    // or bare halides, alkali and alkaline earth metal ions
    FOR_ATOMS_OF_MOL(atom, *mol) {

      if ((atom->GetAtomicNum() == 7 || atom->GetAtomicNum() == 15)
          && atom->GetExplicitValence() == 4) {
        // check if we should make a positive charge?
        // i.e., 4 non-metal neighbors
        bool nonMetalNeighbors = true;
        FOR_NBORS_OF_ATOM(neighbor, &*atom) {
          switch (neighbor->GetAtomicNum()) {
          case 1:
          case 5: case 6: case 7: case 8: case 9:
          case 14: case 15: case 16: case 17:
          case 33: case 34: case 35:
          case 53:
            continue; // good non-metals
          default:
            nonMetalNeighbors = false;
            break; // stop looking
          }
        }
        if (nonMetalNeighbors) // 4 non-metals, e.g. NH4+
          atom->SetFormalCharge(+1);
      }

      // Now look for simple atomic ions like Na, Li, F, Cl, Br...
      // If we have an existing formal charge, keep going
      if (atom->GetFormalCharge() != 0)
        continue;

      // If we're connected to anything besides H2O, keep going
      if (atom->GetExplicitDegree() != 0) {
        int nonWaterBonds = 0;
        FOR_NBORS_OF_ATOM(neighbor, &*atom) {
          if (!CIFisWaterOxygen(&*neighbor)) {
            nonWaterBonds = 1;
            break;
          }
        }
        if (nonWaterBonds)
          continue; // look at another atom
      }

      switch(atom->GetAtomicNum()) {
      case 3: case 11: case 19: case 37: case 55: case 87:
        // Alkali ions
        atom->SetFormalCharge(+1);
        break;
      case 4: case 12: case 20: case 38: case 56: case 88:
        // Alkaline earth ions
        atom->SetFormalCharge(+2);
        break;
      case 9: case 17: case 35: case 53: case 85:
        // Halides
        atom->SetFormalCharge(-1);
        break;
      }
    }
  }

  /////////////////////////////////////////////////////////////////
  bool CIFFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    // If installed, use the mmCIF parser to read CIF
    OBFormat *obformat = OBFormat::FindType("mmcif");
    if (obformat) { return obformat->ReadMolecule(pOb, pConv); }
    obErrorLog.ThrowError(__FUNCTION__, "mmCIF parser not found. Using CIF parser.", obDebug);

    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    CIF cif(*pConv->GetInStream(),true);
    // Loop on all data blocks until we find one structure :@todo: handle multiple structures
    for(map<string,CIFData>::iterator pos=cif.mvData.begin();pos!=cif.mvData.end();++pos)
      if(pos->second.mvAtom.size()>0)
        {
          pmol->BeginModify();
          if(pos->second.mvLatticePar.size()==6)
            {// We have one unit cell
              string spg=pos->second.mSpacegroupSymbolHall;
              if(spg=="") spg=pos->second.mSpacegroupHermannMauguin;
              if(spg=="") spg=pos->second.mSpacegroupNumberIT;
              if(spg=="") spg="P1";
              OBUnitCell *pCell=new OBUnitCell;
              pCell->SetOrigin(fileformatInput);
              pCell->SetData(pos->second.mvLatticePar[0],
                             pos->second.mvLatticePar[1],
                             pos->second.mvLatticePar[2],
                             pos->second.mvLatticePar[3]/DEG_TO_RAD,
                             pos->second.mvLatticePar[4]/DEG_TO_RAD,
                             pos->second.mvLatticePar[5]/DEG_TO_RAD);
              pCell->SetSpaceGroup(spg);
              pCell->SetSpaceGroup(pos->second.mSpaceGroup);
              pmol->SetData(pCell);
            }
          if(pos->second.mName!="") pmol->SetTitle(pos->second.mName);
          else
            if(pos->second.mFormula!="") pmol->SetTitle(pos->second.mFormula);
            else pmol->SetTitle(pConv->GetTitle());

          if(pos->second.mFormula!="") pmol->SetFormula(pos->second.mFormula);

          // Keep a map linking the cif atom label to the obatom*, for bond interpretation later
          std::map<std::string,OBAtom *> vLabelOBatom;

          const unsigned int nbatoms=pos->second.mvAtom.size();
          pmol->ReserveAtoms(nbatoms);
          for(vector<CIFData::CIFAtom>::const_iterator posat=pos->second.mvAtom.begin();posat!=pos->second.mvAtom.end();++posat)
            {
              // Problem: posat->mSymbol is not guaranteed to actually be a symbol
              // see http://www.iucr.org/iucr-top/cif/cifdic_html/1/cif_core.dic/Iatom_type_symbol.html
              // Try to strip the string to have a better chance to have a valid symbol
              // This is not guaranteed to work still, as the CIF standard allows about any string...
              string tmpSymbol=posat->mSymbol;
              unsigned int nbc=0;
              if((tmpSymbol.size()==1) && isalpha(tmpSymbol[0])) nbc=1;
              else if(tmpSymbol.size()>=2)
                {
                  if(isalpha(tmpSymbol[0]) && isalpha(tmpSymbol[1])) nbc=2;
                  else if(isalpha(tmpSymbol[0])) nbc=1;
                }

              OBAtom *atom  = pmol->NewAtom();

              vLabelOBatom.insert(make_pair(posat->mLabel,atom));

              if(tmpSymbol.size()>nbc)
                {// Try to find a formal charge in the symbol
                  int charge=0;
                  int sign=0;
                  for(unsigned int i=nbc;i<tmpSymbol.size();++i)
                    {// Use first number found as formal charge
                      if(isdigit(tmpSymbol[i]) && (charge==0)) charge=atoi(tmpSymbol.substr(i,1).c_str());
                      if('-'==tmpSymbol[i]) sign-=1;
                      if('+'==tmpSymbol[i]) sign+=1;
                    }
                  if(0!=sign) // no sign, no charge
                    {
                      if(charge==0) charge=1;
                      stringstream ss;
                      ss << tmpSymbol<<" / symbol="<<tmpSymbol.substr(0,nbc)<<" charge= "<<sign*charge;
                      obErrorLog.ThrowError(__FUNCTION__, ss.str(), obDebug);
                      atom->SetFormalCharge(sign*charge);
                    }
                }

              if(nbc>0) tmpSymbol=tmpSymbol.substr(0,nbc);
              else tmpSymbol="C";//Something went wrong, no symbol ! Default to C ??

              int atomicNum = OBElements::GetAtomicNum(tmpSymbol.c_str());
              // Test for some oxygens with subscripts
              if (atomicNum == 0 && tmpSymbol[0] == 'O') {
                atomicNum = 8; // e.g. Ob, OH, etc.
              }

              atom->SetAtomicNum(atomicNum); //set atomic number, or '0' if the atom type is not recognized
              atom->SetType(tmpSymbol); //set atomic number, or '0' if the atom type is not recognized
              atom->SetVector(posat->mCoordCart[0],posat->mCoordCart[1],posat->mCoordCart[2]);
              if(posat->mLabel.size()>0)
              {
                OBPairData *label = new OBPairData;
                label->SetAttribute("_atom_site_label");
                label->SetValue(posat->mLabel);
                label->SetOrigin(fileformatInput);
                atom->SetData(label);
              }

              OBPairFloatingPoint *occup_data = new OBPairFloatingPoint;
              occup_data->SetAttribute("_atom_site_occupancy");
              occup_data->SetValue(posat->mOccupancy);
              occup_data->SetOrigin(fileformatInput);
              atom->SetData(occup_data);

              if( posat->mCharge != NOCHARGE )
              {
                OBPairFloatingPoint *charge_data = new OBPairFloatingPoint;
                charge_data->SetAttribute("input_charge");
                charge_data->SetValue(posat->mCharge);
                charge_data->SetOrigin(fileformatInput);
                atom->SetData(charge_data);
              }
            }
          if (!pConv->IsOption("b",OBConversion::INOPTIONS))
            pmol->ConnectTheDots();
          if (pConv->IsOption("B",OBConversion::INOPTIONS))
            {
              for(vector<CIFData::CIFBond>::const_iterator posbond=pos->second.mvBond.begin();posbond!=pos->second.mvBond.end();++posbond)
                {// Add bonds present in the cif and not detected by ConnectTheDots()
                  std::map<std::string,OBAtom *>::iterator posat1,posat2;
                  posat1=vLabelOBatom.find(posbond->mLabel1);
                  posat2=vLabelOBatom.find(posbond->mLabel2);
                  if(posat1!=vLabelOBatom.end() && posat2!=vLabelOBatom.end())
                    {
                      stringstream ss;
                      ss << "  Adding cif bond ? "<<posat1->first<<"-"<<posat2->first;
                      obErrorLog.ThrowError(__FUNCTION__, ss.str(), obDebug);
                      if(pmol->GetBond(posat1->second,posat2->second)==NULL)
                        {
                           obErrorLog.ThrowError(__FUNCTION__, "  :Bond added !", obDebug);
                           OBBond * bond=pmol->NewBond();
                           bond->SetBegin(posat1->second);
                           bond->SetEnd(posat2->second);
                           bond->SetBondOrder(1);
                           bond->SetLength(double(posbond->mDistance));
                        }
                       else obErrorLog.ThrowError(__FUNCTION__, "  :Bond already present.. ", obDebug);
                    }
                }
            }
          if (!pConv->IsOption("s",OBConversion::INOPTIONS) && !pConv->IsOption("b",OBConversion::INOPTIONS))
            pmol->PerceiveBondOrders();
          pmol->EndModify();
          pmol->SetAutomaticFormalCharge(false); // we should have set formal charges
          CorrectFormalCharges(pmol); // Look for lone Na -> Na+, etc.
          return true;
        }

    // If we got here, no structure was found
    obErrorLog.ThrowError(__FUNCTION__, "Problems reading a CIF file: no structure found !" , obWarning);
    return(false);
  }

  ////////////////////////////////////////////////////////////////

  bool CIFFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;
    ostream &ofs = *pConv->GetOutStream();

    char buffer[BUFF_SIZE];

    ofs <<"# CIF file generated by openbabel "<<BABEL_VERSION<<", see http://openbabel.sf.net"<<endl;

    ofs << "data_I"<<endl;
    // Use pmol->GetTitle() as chemical name, though it will probably be the file name
    ofs <<"_chemical_name_common '"<<pmol->GetTitle()<<"'"<<endl;
    // Print Unit cell if we have it
    OBUnitCell *pUC=NULL;
    if (pmol->HasData(OBGenericDataType::UnitCell))
      {
        pUC = (OBUnitCell*)pmol->GetData(OBGenericDataType::UnitCell);
        ofs << "_cell_length_a " << pUC->GetA() << endl
            << "_cell_length_b " << pUC->GetB() << endl
            << "_cell_length_c " << pUC->GetC() << endl
            << "_cell_angle_alpha " << pUC->GetAlpha() << endl
            << "_cell_angle_beta "  << pUC->GetBeta() << endl
            << "_cell_angle_gamma " << pUC->GetGamma() << endl;
        // Save the space group if known
        const SpaceGroup* pSG = pUC->GetSpaceGroup();
        if (pSG != NULL)
          {
            // Do we have an extended HM symbol, with origin choice as ":1" or ":2" ? If so, remove it.
            size_t n=pSG->GetHMName().find(":");
            if(n==string::npos)
              ofs << "_space_group_name_H-M_alt '" << pSG->GetHMName() << "'" << endl;
            else
              ofs << "_space_group_name_H-M_alt '" << pSG->GetHMName().substr(0,n) << "'" << endl;
            ofs << "_space_group_name_Hall '" << pSG->GetHallName() << "'" << endl;
            ofs << "loop_" <<endl
                << "    _symmetry_equiv_pos_as_xyz" << endl;
            transform3dIterator ti;
            const transform3d *t = pSG->BeginTransform(ti);
            while(t)
              {
                ofs << "    " << t->DescribeAsString() << endl;
                t = pSG->NextTransform(ti);
              }
          }
      }

    ofs << "loop_"                      << endl
        << "    _atom_site_label"       << endl
        << "    _atom_site_type_symbol" << endl
        << "    _atom_site_fract_x"     << endl
        << "    _atom_site_fract_y"     << endl
        << "    _atom_site_fract_z"     << endl
        << "    _atom_site_occupancy"   << endl;
    unsigned int i = 0;
    FOR_ATOMS_OF_MOL(atom, *pmol)
      {
         double X, Y, Z; //atom coordinates
         vector3 v = atom->GetVector();
         if (pUC != NULL) {
           v = pUC->CartesianToFractional(v);
           v = pUC->WrapFractionalCoordinate(v);
         }
         X = v.x();
         Y = v.y();
         Z = v.z();
         string label_str;
         double occup;

         if (atom->HasData("_atom_site_occupancy"))
           {
             occup = (dynamic_cast<OBPairFloatingPoint *> (atom->GetData("_atom_site_occupancy")))->GetGenericValue();
           }
         else occup = 1.0;

         if (atom->HasData("_atom_site_label"))
           {
             OBPairData *label = dynamic_cast<OBPairData *> (atom->GetData("_atom_site_label"));
             label_str = label->GetValue().c_str();
           }
         else
           {
             label_str = OBElements::GetSymbol(atom->GetAtomicNum()) + to_string(i);
             i++;
           }

         snprintf(buffer, BUFF_SIZE, "    %-8s%-5s%.5f%10.5f%10.5f%8.3f\n",
                  label_str.c_str(), OBElements::GetSymbol(atom->GetAtomicNum()),
                  X, Y, Z, occup);

         ofs << buffer;
      }
    return true;
  }//WriteMolecule
} //namespace OpenBabel
