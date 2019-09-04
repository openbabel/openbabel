//
// Copyright (C) 2010 David C. Lonie
// Copyright (C) 2018 Patrick Avery
//
// Molekel - Molecular Visualization Program
// Copyright (C) 2006, 2007 Swiss National Supercomputing Centre (CSCS)
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
// MA  02110-1301, USA.
//
// $Author$
// $Date$
// $Revision$
//

// STD
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cstring>

#include <openbabel/obconversion.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/elements.h>
#include <openbabel/generic.h>
#include <cstdlib>


#include <openbabel/griddata.h>

#define EV_TO_KCAL_PER_MOL 23.060538

using namespace std;
using namespace OpenBabel;


static const double BOHR_TO_ANGSTROM = 0.529177249;

namespace OpenBabel {

  class ADFOutputFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    ADFOutputFormat()
    {
      OBConversion::RegisterFormat("adfout",this);
    }

    virtual const char* Description() //required
    {
      return
        "ADF output format\n"
        "Read Options e.g. -as\n"
        "  s  Output single bonds only\n"
        "  b  Disable bonding entirely\n\n";
    };

    virtual const char* SpecificationURL()
    {return "http://www.scm.com/";}; //optional

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
      return READONEONLY | NOTWRITABLE;
    };

    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
  };
  //***

  //Make an instance of the format class
  ADFOutputFormat theADFOutputFormat;

  /////////////////////////////////////////////////////////////////
  bool ADFOutputFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title = pConv->GetTitle();

    char buffer[BUFF_SIZE];
    string str,str1;
    double x,y,z;
    OBAtom *atom;
    vector<string> vs;

    int charge = 0;
    unsigned int spin = 1;
    bool hasPartialCharges = false;

    mol.BeginModify();

    while	(ifs.getline(buffer,BUFF_SIZE))
      {
        if(strstr(buffer,"Coordinates (Cartesian)") != NULL)
          {
            mol.Clear();
            mol.BeginModify();
            ifs.getline(buffer,BUFF_SIZE);	// =========
            ifs.getline(buffer,BUFF_SIZE);	// blank line
            ifs.getline(buffer,BUFF_SIZE);  // Column headings (Atom, etc.)
            ifs.getline(buffer,BUFF_SIZE);  // Column headings 2nd line (X Y Z)
            ifs.getline(buffer,BUFF_SIZE);  // ---------

            ifs.getline(buffer,BUFF_SIZE);  // actual data
            tokenize(vs,buffer);
            while (strstr(buffer, "----") == NULL && vs.size() >= 8)
              {
                atom = mol.NewAtom();
                atom->SetAtomicNum(OBElements::GetAtomicNum(vs[1].c_str())); // atom number, then symbol
                // columns 2, 3, 4 = coordinates in bohr
                x = atof((char*)vs[5].c_str());
                y = atof((char*)vs[6].c_str());
                z = atof((char*)vs[7].c_str());
                atom->SetVector(x,y,z);

                if (!ifs.getline(buffer,BUFF_SIZE))
                  break;
                tokenize(vs,buffer);
              }
          }
        else if(strstr(buffer,"Dipole Moment  ***") != NULL)
          {
            ifs.getline(buffer,BUFF_SIZE);	// =========
            ifs.getline(buffer,BUFF_SIZE);	// blank line
            ifs.getline(buffer,BUFF_SIZE); // actual components  Vector: ###  #### ###
            tokenize(vs,buffer);
            if (vs.size() >= 5)
              {
                OBVectorData *dipoleMoment = new OBVectorData;
                dipoleMoment->SetAttribute("Dipole Moment");
                double x, y, z;
                x = atof(vs[2].c_str());
                y = atof(vs[3].c_str());
                z = atof(vs[4].c_str());
                dipoleMoment->SetData(x, y, z);
                dipoleMoment->SetOrigin(fileformatInput);
                mol.SetData(dipoleMoment);
              }
            if (!ifs.getline(buffer,BUFF_SIZE)) break;
          }
          else if(strstr(buffer,"M U L L I K E N") != NULL)
            {
              ifs.getline(buffer,BUFF_SIZE);	// ========
              ifs.getline(buffer,BUFF_SIZE);	// (blank)
              ifs.getline(buffer,BUFF_SIZE);	// The survey below
              ifs.getline(buffer,BUFF_SIZE);	// a)
              ifs.getline(buffer,BUFF_SIZE);	// b)
              ifs.getline(buffer,BUFF_SIZE);	// c)
              ifs.getline(buffer,BUFF_SIZE);	//
              ifs.getline(buffer,BUFF_SIZE);	// Atom (column headings)
              ifs.getline(buffer,BUFF_SIZE);  // ---

              ifs.getline(buffer,BUFF_SIZE); // Actual data!
              tokenize(vs,buffer);
              while (vs.size() >= 3)
                {
                  atom = mol.GetAtom(atoi(vs[0].c_str()));
                  if (atom) {
                    atom->SetPartialCharge(atof(vs[2].c_str()));
                    hasPartialCharges = true;
                  }

                  if (!ifs.getline(buffer,BUFF_SIZE))
                    break;
                  tokenize(vs,buffer);
                }
            }
        else if(strstr(buffer,"Net Charge") != NULL)
          {
            tokenize(vs, buffer);
            if (vs.size() > 3) // Net Charge: ##
              {
                charge = atoi(vs[2].c_str());
              }
          }
        else if (strstr(buffer,"Bond Energy") != NULL)
          {
            double energy = 0;
            for (;;) {
              if (strstr(buffer, "a.u.")) {
                ifs.getline(buffer,BUFF_SIZE);
                continue;
              }
              else if (strstr(buffer, "eV")) {
                tokenize(vs, buffer);
                // The energies have a variable number of
                // prefixes...rather than check each possible case,
                // just try to convert tokens 1 thru vs.size() to
                // double and keep the last one that works.
                for (unsigned int i = 1; i < vs.size(); i++) {
                  if (double d = atof(vs.at(i).c_str())) {
                    energy = d;
                  }
                }
                ifs.getline(buffer,BUFF_SIZE);
              }
              else break;
            }
            mol.SetEnergy(energy * EV_TO_KCAL_PER_MOL);
          }

      } // end while

    if (mol.NumAtoms() == 0) { // e.g., if we're at the end of a file PR#1737209
      mol.EndModify();
      return false;
    }

    if (!pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.ConnectTheDots();
    if (!pConv->IsOption("s",OBConversion::INOPTIONS) && !pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();

    mol.EndModify();

    if (hasPartialCharges) {
      mol.SetPartialChargesPerceived();
      // Annotate that partial charges come from Q-Chem Mulliken
      OBPairData *dp = new OBPairData;
      dp->SetAttribute("PartialCharges");
      dp->SetValue("Mulliken");
      dp->SetOrigin(perceived);
      mol.SetData(dp);
    }
    mol.SetTotalCharge(charge);
    mol.SetTotalSpinMultiplicity(spin);

    mol.SetTitle(title);
    return(true);
  }

  class ADFInputFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    ADFInputFormat()
    {
      OBConversion::RegisterFormat("adf", this);
    }

    virtual const char* Description() //required
    {
      return
        "ADF cartesian input format\n"
        "Read Options e.g. -as\n"
        "  s  Output single bonds only\n"
        "  b  Disable bonding entirely\n\n";
    };

    virtual const char* SpecificationURL()
    {return "http://www.scm.com/Doc/Doc2007.01/ADF/ADFUsersGuide/page32.html";}; //optional

    //*** This section identical for most OBMol conversions ***
    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv)
      { return false; }
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

    virtual unsigned int Flags()
    {
      return NOTREADABLE | WRITEONEONLY;
    };
  };
  //***

  //Make an instance of the format class
  ADFInputFormat theADFInputFormat;

  bool ADFInputFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    char buffer[BUFF_SIZE];

    snprintf(buffer, BUFF_SIZE, "TITLE %s\n\n", mol.GetTitle());
    ofs << buffer;

    // Output CHARGE and spin
    // Note that ADF expects the spin parameter to be the # of unpaired spins
    // (i.e., singlet = 0, doublet = 1, etc.)
    snprintf(buffer, BUFF_SIZE, "CHARGE %d  %d\n\n",
             mol.GetTotalCharge(),
             mol.GetTotalSpinMultiplicity() - 1);
    ofs << buffer;

    // Cartesian input -- change this if you want a z-matrix format
    snprintf(buffer, BUFF_SIZE, "Number of atoms\n %d\n\n", mol.NumAtoms());
    ofs << buffer;

    ofs << "ATOMS Cartesian\n";
    FOR_ATOMS_OF_MOL(atom, mol)
      {
        snprintf(buffer, BUFF_SIZE, "%-3s%15.5f%15.5f%15.5f\n",
                 OBElements::GetSymbol(atom->GetAtomicNum()),
                 atom->GetX(),
                 atom->GetY(),
                 atom->GetZ());
        ofs << buffer;
      }
    ofs << "End\n\n";

    // command-line keywords (-xk "blah")
    const char *keywords = pConv->IsOption("k",OBConversion::OUTOPTIONS);
    const char *keywordFile = pConv->IsOption("f",OBConversion::OUTOPTIONS);

    // If the user specified a full file, pick that over anything else
    if (keywordFile) {
        ifstream kfstream(keywordFile);
        string keyBuffer;
        if (kfstream)
          {
            while (getline(kfstream, keyBuffer))
              ofs << keyBuffer << endl;
          }
      }
    else if (keywords) {
      ofs << keywords << endl;
    }
    else {
      ofs << "Basis\n";
      ofs << "End\n\n";

      ofs << "Geometry\n";
      ofs << "End\n\n";
    }

    ofs << endl; // one final blank line

    return true;
  }

  class ADFBandFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    ADFBandFormat()
    {
      OBConversion::RegisterFormat("adfband",this);
    }

    virtual const char* Description() //required
    {
      return "ADF Band output format\n";
    };

    virtual const char* SpecificationURL()
    {return "https://www.scm.com/product/band_periodicdft/";}; //optional

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
      return READONEONLY | NOTWRITABLE;
    };

    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
  };
  //***

  //Make an instance of the format class
  ADFBandFormat theADFBandFormat;

  /////////////////////////////////////////////////////////////////
  bool ADFBandFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if (!pmol)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title = pConv->GetTitle();

    char buffer[BUFF_SIZE];
    vector<string> vs;

    // This will remain as 1.0 if we have units of Angstroms
    double lengthConversion = 1.0;

    mol.BeginModify();

    while (ifs.getline(buffer, BUFF_SIZE)) {
      if (strstr(buffer, "length Bohr") || strstr(buffer, "length BOHR") ||
          strstr(buffer, "length bohr")) {
        // We have units of Bohr!
        lengthConversion = BOHR_TO_ANGSTROM;
      }
      else if (strstr(buffer,
                      "G E O M E T R Y    I N    X - Y - Z    F O R M A T")) {
        // We need to clear the atoms before proceeding. Since this comes
        // before all the other data, we can just clear the whole molecule
        mol.Clear();
        mol.BeginModify();

        ifs.getline(buffer, BUFF_SIZE); // ===========
        ifs.getline(buffer, BUFF_SIZE); // *blank line*
        while (ifs.getline(buffer, BUFF_SIZE)) {
          tokenize(vs, buffer);
          if (vs.size() < 4 || vs[0] == "VEC1")
            break;

          OBAtom* atom = mol.NewAtom();
          atom->SetAtomicNum(OBElements::GetAtomicNum(vs[0].c_str()));
          double x = atof(vs[1].c_str()) * lengthConversion;
          double y = atof(vs[2].c_str()) * lengthConversion;
          double z = atof(vs[3].c_str()) * lengthConversion;
          atom->SetVector(x, y, z);
        }
      }
      else if (strstr(buffer, "REAL SPACE LATTICE VECTORS")) {
        ifs.getline(buffer, BUFF_SIZE); // ------------

        std::vector<vector3> vectors;
        for (int i = 0; i < 3; ++i) {
          ifs.getline(buffer, BUFF_SIZE);
          tokenize(vs, buffer);
          if (vs.size() < 5)
            break;

          // These are in Bohrs
          double x = atof(vs[1].c_str()) * BOHR_TO_ANGSTROM;
          double y = atof(vs[2].c_str()) * BOHR_TO_ANGSTROM;
          double z = atof(vs[3].c_str()) * BOHR_TO_ANGSTROM;
          vectors.push_back(vector3(x, y, z));
        }

        while (vectors.size() < 3)
          vectors.push_back(vector3(0.0, 0.0, 0.0));

        // Build unit cell
        OBUnitCell* cell = new OBUnitCell;
        cell->SetData(vectors[0], vectors[1], vectors[2]);
        cell->SetSpaceGroup(1);
        pmol->SetData(cell);
      }
      else if (strstr(buffer, "E N E R G Y   A N A L Y S I S")) {
        // Final bond energy line looks like this:
        //  Final bond energy (LDA)                                      -1.09942745     -29.9169     -689.90
        while (ifs.getline(buffer, BUFF_SIZE)) {
          if (strstr(buffer, "Final bond energy")) {
            tokenize(vs, buffer);

            // Line should be of size 7
            if (vs.size() != 7)
              break;

            // Units of the final column should be in kcal/mol
            mol.SetEnergy(atof(vs[6].c_str()));
            break;
          }
        }
      }
    }

    if (mol.NumAtoms() == 0) { // e.g., if we're at the end of a file
      mol.EndModify();
      return false;
    }

    mol.EndModify();

    mol.SetTitle(title);
    return true;
  }

  class ADFDftbFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    ADFDftbFormat()
    {
      OBConversion::RegisterFormat("adfdftb",this);
    }

    virtual const char* Description() //required
    {
      return "ADF DFTB output format\n";
    };

    virtual const char* SpecificationURL()
    {return "https://www.scm.com/product/dftb/";}; //optional

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
      return READONEONLY | NOTWRITABLE;
    };

    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
  };
  //***

  //Make an instance of the format class
  ADFDftbFormat theADFDftbFormat;

  /////////////////////////////////////////////////////////////////
  bool ADFDftbFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if (!pmol)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title = pConv->GetTitle();

    char buffer[BUFF_SIZE];
    vector<string> vs;

    mol.BeginModify();

    while (ifs.getline(buffer, BUFF_SIZE)) {
      if (strcmp(buffer, "Geometry") == 0) {
        // We need to clear the atoms before proceeding. Since this comes
        // before all the other data, we can just clear the whole molecule
        mol.Clear();
        mol.BeginModify();

        ifs.getline(buffer, BUFF_SIZE); // ------------
        ifs.getline(buffer, BUFF_SIZE); // Atoms
        ifs.getline(buffer, BUFF_SIZE);

        // Make sure it is the correct line
        if (strstr(buffer, "Index") && strstr(buffer, "Symbol")) {
          double lengthConversion = 1.0;
          // Check the units
          if (strstr(buffer, "bohr"))
            lengthConversion = BOHR_TO_ANGSTROM;

          while (ifs.getline(buffer, BUFF_SIZE)) {
            tokenize(vs, buffer);
            // Should be of size 5
            if (vs.size() < 5)
              break;

            OBAtom* atom = mol.NewAtom();
            atom->SetAtomicNum(OBElements::GetAtomicNum(vs[1].c_str()));
            double x = atof(vs[2].c_str()) * lengthConversion;
            double y = atof(vs[3].c_str()) * lengthConversion;
            double z = atof(vs[4].c_str()) * lengthConversion;
            atom->SetVector(x, y, z);
          }
        }

        // Now read the lattice vectors
        ifs.getline(buffer, BUFF_SIZE);
        if (strstr(buffer, "Lattice vectors")) {
          double lengthConversion = 1.0;
          // Check the units
          if (strstr(buffer, "bohr"))
            lengthConversion = BOHR_TO_ANGSTROM;

          std::vector<vector3> vectors;
          for (short i = 0; i < 3; ++i) {
            ifs.getline(buffer, BUFF_SIZE);
            tokenize(vs, buffer);
            if (vs.size() != 4)
              break;

            double x = atof(vs[1].c_str()) * lengthConversion;
            double y = atof(vs[2].c_str()) * lengthConversion;
            double z = atof(vs[3].c_str()) * lengthConversion;
            vectors.push_back(vector3(x, y, z));
          }

          while (vectors.size() < 3)
            vectors.push_back(vector3(0.0, 0.0, 0.0));

          // Build unit cell
          OBUnitCell* cell = new OBUnitCell;
          cell->SetData(vectors[0], vectors[1], vectors[2]);
          cell->SetSpaceGroup(1);
          pmol->SetData(cell);
        }
      }
      else if (strcmp(buffer, "Energies") == 0 ||
               strcmp(buffer, "Energy Decomposition") == 0) {
        // Final energy line looks like this:
        // Total Energy (eV)                   -220.34976964
        while (ifs.getline(buffer, BUFF_SIZE)) {
          if (strstr(buffer, "Total Energy (eV)")) {
            tokenize(vs, buffer);

            // Line should be of size 4
            if (vs.size() != 4)
              break;

            mol.SetEnergy(atof(vs[3].c_str()) * EV_TO_KCAL_PER_MOL);
            break;
          }
        }
      }
    }

    if (mol.NumAtoms() == 0) { // e.g., if we're at the end of a file
      mol.EndModify();
      return false;
    }

    mol.EndModify();

    mol.SetTitle(title);
    return true;
  }


class OBT41Format : public OBMoleculeFormat
{
public:
    /// Constructor: register 't41' and "T41" format.
    OBT41Format()
    {
        OBConversion::RegisterFormat( "t41", this );
        OBConversion::RegisterFormat( "T41", this );
    }

    /// Return description.
    virtual const char* Description() //required
    {
        return
        "ADF TAPE41 format\n\n"

        "Currently the ADF Tape41 support reads grids from\n"
        "TAPE41 text files. To generate an ASCII version from\n"
        "the default binary, use the dmpkf program.\n\n"

        "Read Options e.g. -as\n"
        "  s  Output single bonds only\n"
        "  b  Disable bonding entirely\n\n";
    }

    /// Return a specification url, not really a specification since
    /// I couldn't find it but close enough.
    virtual const char* SpecificationURL()
    {
        return "http://www.scm.com/Doc/Doc2006.01/ADF/Analysis/page8.html";
    }

    /// Return MIME type, NULL in this case.
    virtual const char* GetMIMEType() { return 0; };

      /// Return read/write flag: read only.
    virtual unsigned int Flags()
    {
        return READONEONLY | READBINARY | NOTWRITABLE;
    };

    /// Skip to object: used for multi-object file formats.
    virtual int SkipObjects( int n, OBConversion* pConv ) { return 0; }

    /// Read.
    virtual bool ReadMolecule( OBBase* pOb, OBConversion* pConv );

    bool ReadASCII( OBBase* pOb, OBConversion* pConv);
    bool ReadBinary( OBBase* pOb, OBConversion* pConv);

    /// Write: always returns false.
    virtual bool WriteMolecule( OBBase* , OBConversion* )
    { return false; }

private:
    ///Utility function that eats all the remaining characters on the current and next line.
    void eol( istream& is ) const { string s; getline( is, s ); getline( is, s ); }

    ///Advance to next tag.
    bool NextTag( istream& is, const std::string& tag ) const
    {
        string buf = "";
        while( is >> buf ) if( buf == tag ) return true;
        return false;
    }

    /// 3x3 Matrix - 3d vector inplace multiply.
    void MatVecMul( const double xColumn[ 3 ],
                    const double yColumn[ 3 ],
                    const double zColumn[ 3 ],
                    double v[ 3 ] )
    {
        double t[ 3 ];
        t[ 0 ] = v[ 0 ];
        t[ 1 ] = v[ 1 ];
        t[ 2 ] = v[ 2 ];
        v[ 0 ] = xColumn[ 0 ] * t[ 0 ] + yColumn[ 0 ] * t[ 1 ] + zColumn[ 0 ] * t[ 2 ];
        v[ 1 ] = xColumn[ 1 ] * t[ 0 ] + yColumn[ 1 ] * t[ 1 ] + zColumn[ 1 ] * t[ 2 ];
        v[ 2 ] = xColumn[ 2 ] * t[ 0 ] + yColumn[ 2 ] * t[ 1 ] + zColumn[ 2 ] * t[ 2 ];
    }


    /// Add grids from SCF
    void AddSCFGrids( istream& is, OBGridData& t41 ) {}

    ///Inner class used to hold atomic number, coordinate, charge data
    struct AtomData
    {
        int atomicNum;
        double coord[ 3 ];
        double charge;
        AtomData() : atomicNum( 0 ), charge( 0. ) {}
        AtomData( int an ) : atomicNum( an ) {}
        AtomData( int an, const double ac[ 3 ], double c ) : atomicNum( an ), charge( c )
        {
            coord[ 0 ] = ac[ 0 ];
            coord[ 1 ] = ac[ 1 ];
            coord[ 2 ] = ac[ 2 ];
        }
    };

    ///Inner class used to store grid data info before adding it
    ///into an OBMol instance.
    struct T41GridData
    {
        bool valid;
        double startPoint[ 3 ];
        int numPoints[ 3 ];
        double xAxis[ 3 ];
        double yAxis[ 3 ];
        double zAxis[ 3 ];
        int numSymmetries;
        std::vector< std::string > labels;
        bool unrestricted;
        T41GridData() : valid( false ) {}
        operator bool() const { return valid; }
    };

    typedef T41GridData GridData;

    ///Read grid data.
    GridData ReadGridData( istream& is ) const;

    ///Read SCF grids.
    bool ReadSCFGrid( istream& is, OBGridData& t41Data ) const;

    ///Read SCF orbital grids.
    bool ReadSCFOrbitalGrid( istream& is, OBGridData& t41Data ) const;

    ///Read SumFrag grids.
    bool ReadSumFragGrid( istream& is, OBGridData& t41Data ) const;

                OBGridData *NewData(const GridData &gd);
};

//------------------------------------------------------------------------------

// Global variable used to register Tape41 format.
OBT41Format t41Format__;

//------------------------------------------------------------------------------


//==============================================================================

OBGridData *OBT41Format::NewData(const T41GridData &gd)
{
        OBGridData *t41Data = new OBGridData;
  t41Data->SetNumberOfPoints( gd.numPoints[ 0 ], gd.numPoints[ 1 ], gd.numPoints[ 2 ] );
        t41Data->SetLimits( gd.startPoint, gd.xAxis, gd.yAxis, gd.zAxis );
  t41Data->SetUnrestricted( gd.unrestricted );
  t41Data->SetNumSymmetries( gd.numSymmetries );

        return t41Data;
}

//------------------------------------------------------------------------------
bool OBT41Format::ReadMolecule( OBBase* pOb, OBConversion* pConv )
{
  istream& ifs = *pConv->GetInStream();
  // Check if the file is ASCII or Binary
  // Binary TAPE41 files start with "SUPERINDEX"
  if (ifs.peek() == 'S')
    return ReadBinary(pOb, pConv);
  else
    return ReadASCII(pOb, pConv);
}

bool OBT41Format::ReadBinary( OBBase* pOb, OBConversion* pConv )
{
  obErrorLog.ThrowError( __FUNCTION__, "OpenBabel does not currently support the TAPE41 binary format. Please use dmpkf to convert to ASCII.", obError );
  return false;
}

bool OBT41Format::ReadASCII( OBBase* pOb, OBConversion* pConv )
{
      OBMol* pmol = dynamic_cast< OBMol* >(pOb);
      if( pmol == 0 ) return false;

      istream& ifs = *pConv->GetInStream();

      GridData gd;
      gd = ReadGridData( ifs );

      OBGridData* t41Data = 0;
      if( gd )
      {
         streampos current = ifs.tellg();

         // We create new data for each grid
         // If we don't find any legitimate data, we'll end and still have "OBGridData"
         // rather than a legitimate label
                         t41Data = NewData(gd);
         while( ReadSCFOrbitalGrid( ifs, *t41Data ) );
         if (t41Data->GetAttribute() == "GridData") {
           delete t41Data;
         } else
           pmol->SetData( t41Data );

         ifs.clear();
         ifs.seekg( current, ios::beg );

                         t41Data = NewData(gd);
         while( ReadSCFGrid( ifs, *t41Data ) );
         if (t41Data->GetAttribute() == "GridData") {
           delete t41Data;
         } else
           pmol->SetData( t41Data );

         ifs.clear();
         ifs.seekg( current, ios::beg );

                         t41Data = NewData(gd);
         while( ReadSumFragGrid( ifs, *t41Data ) );
         if (t41Data->GetAttribute() == "GridData") {
           delete t41Data;
         } else
           pmol->SetData( t41Data );

         ifs.clear();
         ifs.seekg( current, ios::beg );
      }

      string buf;
      // nuuc
      while( buf != "Geometry" ) ifs >> buf; cout << buf << endl;
      ifs >> buf; cout << buf << endl;
      if( buf != "nnuc" )
      {
          obErrorLog.ThrowError( __FUNCTION__, "no 'nuuc' after first Geometry tag" );
          return false;
      }
      eol( ifs );
      unsigned int numAtoms = 0;
      ifs >> numAtoms; cout << numAtoms << endl;
      buf  = "";

      // labels
      while( buf != "Geometry" ) ifs >> buf; cout << buf << endl;
      ifs >> buf; cout << buf << endl;
      if( buf != "labels" )
      {
          obErrorLog.ThrowError( __FUNCTION__, "no 'labels' after second Geometry tag" );
          return false;
      }
      eol( ifs );
      std::vector< AtomData > atoms;
      atoms.reserve( numAtoms );
      for (unsigned int i = 0; i != numAtoms; ++i)
      {
          ifs >> buf; cout << buf << endl;
          atoms.push_back( OBElements::GetAtomicNum( buf.c_str() ) );
      }
      if( atoms.size() != numAtoms )
      {
          obErrorLog.ThrowError( __FUNCTION__, "wrong number of atoms" );
          return false;
      }
      //coordinates
      buf = "";
      while( buf != "Geometry" ) ifs >> buf; cout << buf << endl;
      ifs >> buf; cout << buf << endl;
      if( buf != "xyznuc" )
      {
          obErrorLog.ThrowError( __FUNCTION__, "no 'xyznuc' after third Geometry tag" );
          return false;
      }
      eol( ifs );
      for (unsigned int i = 0; i != numAtoms; ++i)
      {
          ifs >> atoms[ i ].coord[ 0 ] >> atoms[ i ].coord[ 1 ] >> atoms[ i ].coord[ 2 ];
          cout << atoms[ i ].coord[ 0 ] << ' ' << atoms[ i ].coord[ 1 ] << ' ' << atoms[ i ].coord[ 2 ] << endl;
      }
      //charge
      buf = "";
      while( buf != "Geometry" ) ifs >> buf; cout << buf << endl;
      ifs >> buf; cout << buf << endl;
      if( buf != "qtch" )
      {
          obErrorLog.ThrowError( __FUNCTION__, "no 'qtch' after fourth Geometry tag" );
          return false;
      }
      eol( ifs );
      for (unsigned int i = 0; i != numAtoms; ++i)
      {
          ifs >> atoms[ i ].charge;
      }

      // unit of length
      buf = "";
      while( buf != "Geometry" ) ifs >> buf; cout << buf << endl;
      ifs >> buf >> buf >> buf; cout << buf << endl;
      if( buf != "length" )
      {
          obErrorLog.ThrowError( __FUNCTION__, "no 'unit of length' after fifth Geometry tag" );
          return false;
      }
      eol( ifs );
      double scale = 1.0;
      ifs >> scale;
      /// @todo multply coordinates by axis length;
      for (unsigned int i = 0; i != numAtoms; ++i)
      {

          atoms[ i ].coord[ 0 ] *= BOHR_TO_ANGSTROM;
          atoms[ i ].coord[ 1 ] *= BOHR_TO_ANGSTROM;
          atoms[ i ].coord[ 2 ] *= BOHR_TO_ANGSTROM;
  //        atoms[ i ].coord[ 0 ] *= scale;
  //        atoms[ i ].coord[ 1 ] *= scale;
  //        atoms[ i ].coord[ 2 ] *= scale;
  //        MatVecMul( gd.xAxis, gd.yAxis, gd.zAxis, atoms[ i ].coord );
  //        atoms[ i ].coord[ 0 ] += gd.startPoint[ 0 ];
  //        atoms[ i ].coord[ 1 ] += gd.startPoint[ 1 ];
  //        atoms[ i ].coord[ 2 ] += gd.startPoint[ 2 ];
      }

      // build OB molecule

      pmol->BeginModify();

      pmol->SetDimension( 3 );

      pmol->ReserveAtoms( numAtoms );

      for (unsigned int i = 0; i < numAtoms; ++i)
      {
          OBAtom *atom = pmol->NewAtom();
          atom->SetAtomicNum( atoms[ i ].atomicNum );
          atom->SetVector( atoms[ i ].coord[ 0 ],
                           atoms[ i ].coord[ 1 ],
                           atoms[ i ].coord[ 2 ] );
          atom->SetPartialCharge( atoms[ i ].charge );
      }

      if( !pConv->IsOption( "b", OBConversion::INOPTIONS ) ) pmol->ConnectTheDots();
      if (!pConv->IsOption( "s", OBConversion::INOPTIONS )
          && !pConv->IsOption( "b", OBConversion::INOPTIONS ) )
      {
          pmol->PerceiveBondOrders();
      }
      pmol->EndModify();

      return true;
}

//------------------------------------------------------------------------------
OBT41Format::GridData OBT41Format::ReadGridData( istream& is ) const
{
    GridData gd;
    string buf;
    // Start_point
    if( !NextTag( is, "Grid" ) ) return gd;
    is >> buf;
    if( buf != "Start_point" ) return gd;
    eol( is );
    is >> gd.startPoint[ 0 ] >> gd.startPoint[ 1 ] >> gd.startPoint[ 2 ];

    gd.startPoint[ 0 ] *= BOHR_TO_ANGSTROM;
    gd.startPoint[ 1 ] *= BOHR_TO_ANGSTROM;
    gd.startPoint[ 2 ] *= BOHR_TO_ANGSTROM;

    // nr of points x
    if( !NextTag( is, "Grid" ) ) return gd;
    is >> buf >> buf >> buf >> buf;
    if( buf != "x" ) return gd;
    eol( is );
    is >> gd.numPoints[ 0 ];
    // nr of points y
    if( !NextTag( is, "Grid" ) ) return gd;
    is >> buf >> buf >> buf >> buf;
    if( buf != "y" ) return gd;
    eol( is );
    is >> gd.numPoints[ 1 ];
    // nr of points z
    if( !NextTag( is, "Grid" ) ) return gd;
    is >> buf >> buf >> buf >> buf;
    if( buf != "z" ) return gd;
    eol( is );
    is >> gd.numPoints[ 2 ];
    // total nr of points
    if( !NextTag( is, "Grid" ) ) return gd;
    is >> buf >> buf >> buf >> buf;
    if( buf != "points" ) return gd;
    eol( is );
    int n = 0;
    is >> n;
    if( gd.numPoints[ 0 ] * gd.numPoints[ 1 ] * gd.numPoints[ 2 ] != n ) return gd;
    //x-vector
    if( !NextTag( is, "Grid" ) ) return gd;
    is >> buf;
    if( buf != "x-vector" ) return gd;
    eol( is );
    is >> gd.xAxis[ 0 ] >> gd.xAxis[ 1 ] >> gd.xAxis[ 2 ];

    gd.xAxis[ 0 ] *= BOHR_TO_ANGSTROM;
    gd.xAxis[ 1 ] *= BOHR_TO_ANGSTROM;
    gd.xAxis[ 2 ] *= BOHR_TO_ANGSTROM;

    //y-vector
    if( !NextTag( is, "Grid" ) ) return gd;
    is >> buf;
    if( buf != "y-vector" ) return gd;
    eol( is );
    is >> gd.yAxis[ 0 ] >> gd.yAxis[ 1 ] >> gd.yAxis[ 2 ];

    gd.yAxis[ 0 ] *= BOHR_TO_ANGSTROM;
    gd.yAxis[ 1 ] *= BOHR_TO_ANGSTROM;
    gd.yAxis[ 2 ] *= BOHR_TO_ANGSTROM;

    //z-vector
    if( !NextTag( is, "Grid" ) ) return gd;
    is >> buf;
    if( buf != "z-vector" ) return gd;
    eol( is );
    is >> gd.zAxis[ 0 ] >> gd.zAxis[ 1 ] >> gd.zAxis[ 2 ];

    gd.zAxis[ 0 ] *= BOHR_TO_ANGSTROM;
    gd.zAxis[ 1 ] *= BOHR_TO_ANGSTROM;
    gd.zAxis[ 2 ] *= BOHR_TO_ANGSTROM;

    //nr of symmetries
    if( !NextTag( is, "Grid" ) ) return gd;
    is >> buf >> buf >> buf;
    if( buf != "symmetries" ) return gd;
    eol( is );
    is >> gd.numSymmetries;
    //labels ///@warning only one label supported
    if( !NextTag( is, "Grid" ) ) return gd;
    is >> buf;
    if( buf != "labels" ) return gd;
    eol( is );
    is >> buf; gd.labels.push_back( buf );
    //unrestricted
    if( !NextTag( is, "Grid" ) ) return gd;
    is >> buf;
    if( buf != "unrestricted" ) return gd;
    eol( is );
    char c;
    is >> c;
    gd.unrestricted = ( c == 'T' );

    gd.valid = true;
    return gd;
}

//------------------------------------------------------------------------------
inline bool IsNum( const string& s )
{
    bool isnum = true;
    for (unsigned int i = 0; i != s.size(); ++i)
    {
        if( !isdigit( s[ i ] ) )
        {
            isnum = false;
            break;
        }
    }
    return isnum;
}

bool OBT41Format::ReadSCFOrbitalGrid( istream& is, OBGridData& t41Data ) const
{
    //find next tag starting with 'SCF'
    //if tag starts with SCF_ check next line
    //  advance to next SCF_ tag until tag on next line is a number
    //  read number = orbital id
    //  skip line
    //  read grid data
    if( !is ) return false;
    string buf;
    while( is >> buf ) if( buf.find( "SCF", 0 ) == 0 && buf.size() > 3 ) break;
    if( !is ) return false;
    // SCF_
    const string scf = buf;
    buf = "";
    is >> buf;
    if( !IsNum( buf ) ) // number ?
    {
        while( is >> buf ) // not a number keep on reading
        {
            if( buf == scf ) // SCF_X_X
            {
                is >> buf; // read tag on next line
                if( IsNum( buf ) ) break; // break if (orbital) number
            }
        }
    }
    if( !is ) return false; // eof -> return
    // if we get here it means we are past the orbital number so read
    // read grid values
    const string label = scf + ' ' + buf; cout << label << endl;
    const int numPoints = t41Data.GetNumberOfPoints();
    vector< double > grid( numPoints );
    eol( is );
    if( !is ) return false;
    for( int i = 0; i != numPoints; ++i )
    {
        is >> grid[ i ];
    }

    // Now to translate from a vector stored in z, y, x to one stored in xyz
    // ADF format order to OpenBabel GridData order!
    int voxels[3];
    t41Data.GetNumberOfPoints(voxels[0], voxels[1], voxels[2]);
    for (int k = 0; k < voxels[2]; ++k)
      for (int j = 0; j < voxels[1]; ++j)
        for (int i = 0; i < voxels[0]; ++i)
          {
            t41Data.SetValue(i, j, k,
                             grid[k*voxels[0]*voxels[1] + j*voxels[0] + i]);
          }

                t41Data.SetAttribute( label );
    return true;
}

//------------------------------------------------------------------------------
bool OBT41Format::ReadSCFGrid( istream& is, OBGridData& t41Data ) const
{
        if( !is ) return false;
    string buf;
    while( is >> buf ) if( buf.find( "SCF", 0 ) == 0 && buf.size() == 3 ) break;
    if( !is ) return false;
    // if tag = SCF read next line then skip line and read grid data
    const string scf = buf; // SCF
    is >> buf; // tag on line after SCF
    const string label = scf + ' ' + buf; cout << label << endl;
    eol( is );
    if( !is ) return false;
    // read grid data
    const int numPoints = t41Data.GetNumberOfPoints();
    vector< double > grid( numPoints );
    int i = 0;
    for( ; i != numPoints; ++i ) is >> grid[ i ];

    // Now to translate from a vector stored in z, y, x to one stored in xyz
    // ADF format order to OpenBabel GridData order!
    // Now to translate from a vector stored in z, y, x to one stored in xyz
    // ADF format order to OpenBabel GridData order!
    int voxels[3];
    t41Data.GetNumberOfPoints(voxels[0], voxels[1], voxels[2]);
    for (int k = 0; k < voxels[2]; ++k)
      for (int j = 0; j < voxels[1]; ++j)
        for (int i = 0; i < voxels[0]; ++i)
          {
            t41Data.SetValue(i, j, k,
                             grid[k*voxels[0]*voxels[1] + j*voxels[0] + i]);
          }

                t41Data.SetAttribute( label );
    return true;
}


//------------------------------------------------------------------------------
bool OBT41Format::ReadSumFragGrid( istream& is, OBGridData& t41Data ) const
{
    if( !is ) return false;
    string buf;
    while( is >> buf ) if( buf == "SumFrag" ) break; // look for SumFrag string
    if( !is ) return false; // not found -> return
    const string sumfrag = buf; // found read next line then skip one line and read data
    is >> buf;
    const string label = sumfrag + ' ' + buf; cout << label << endl;
    eol( is );
    if( !is ) return false;
    const int numPoints = t41Data.GetNumberOfPoints();
    vector< double > grid( numPoints );
    int i = 0;
    for( ; i != numPoints; ++i ) is >> grid[ i ];

    // Now to translate from a vector stored in z, y, x to one stored in xyz
    // ADF format order to OpenBabel GridData order!
    int voxels[3];
    t41Data.GetNumberOfPoints(voxels[0], voxels[1], voxels[2]);
    for (int k = 0; k < voxels[2]; ++k)
      for (int j = 0; j < voxels[1]; ++j)
        for (int i = 0; i < voxels[0]; ++i)
          {
            t41Data.SetValue(i, j, k,
                             grid[k*voxels[0]*voxels[1] + j*voxels[0] + i]);
          }

                t41Data.SetAttribute( label );
    return true;
}

} // end namespace OpenBabel
