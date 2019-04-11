/**********************************************************************
Copyright (C) 2007 by Maxim Fedorovsky, University of Fribourg (Switzerland).

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include <numeric>
#include <typeinfo>
#include <functional>
#include <cstdlib>
#include <algorithm>

// No diagnoalization yet. Perhaps for 2.2 -GRH
// #include <eigen/matrix.h>

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/elements.h>
#include <openbabel/generic.h>

#include <openbabel/math/matrix3x3.h>

#define BOHR2ANGSTROM 0.5291772083
#define HARTREE2INVCM 219474.631371
#define AMU2AU 1822.88848
#define REDDIPSTR2INT 974.864

using namespace std;

namespace OpenBabel
{
  //! \brief Class for parsing Gaussian formatted checkpoint files
  class FCHKFormat : public OBMoleculeFormat
  {
  public :
    /* Register this format type ID */
    FCHKFormat()
    {
      OBConversion::RegisterFormat("fchk", this,
                                   "chemical/x-gaussian-checkpoint");
      OBConversion::RegisterFormat("fch", this,
                                   "chemical/x-gaussian-checkpoint");
      OBConversion::RegisterFormat("fck", this,
                                   "chemical/x-gaussian-checkpoint");
    }

    virtual const char * Description()
    {
      return "Gaussian formatted checkpoint file format\n"
             "A formatted text file containing the results of a Gaussian calculation\n"
             "Currently supports reading molecular geometries from fchk files. More to come.\n\n"

             "Read options e.g. -as\n"
             " s  Single bonds only\n"
             " b  No bond perception\n\n";
      // Vibrational analysis not yet supported in OB-2.1.
      //                    v  Do not perform the vibrational analysis\n\n";
    };

    virtual const char * GetMIMEType()
    {
      return "chemical/x-gaussian-checkpoint";
    };

    virtual unsigned int Flags()
    {
      return NOTWRITABLE;
    };

    virtual bool ReadMolecule(OBBase *, OBConversion *);

  private :
    bool static read_int(const char * const, int * const);
    template<class T>  static bool read_numbers(const char * const,
                                                vector<T> &,
                                                const unsigned int width = 0);
    template <class T> static bool read_section(const char * const,
                                                vector<T> &,
                                                const unsigned int,
                                                bool * const,
                                                const char * const,
                                                const unsigned int,
                                                const unsigned int width = 0);
    bool static validate_section(const char * const,
                                 const int,
                                 const char * const,
                                 const unsigned int);
    bool static validate_number(const int,
                                const char * const,
                                const unsigned int);
    void static construct_mol(OBMol * const,
                              OBConversion * const,
                              const unsigned int,
                              const vector<int> &,
                              const vector<double> &,
                              const int,
                              const vector<int> &,
                              const vector<int> &);
    //       void static vibana(OBMol * const,
    //                          const vector<double> &,
    //                          const vector<double> &);
    //       bool static generate_tr_rot(OBMol * const,
    //                                   double * const,
    //                                   unsigned int * const);
    //       void static rotmatrix_axis(const double * const, const double,
    //                                  matrix3x3 &);
    //       double static cosine(const double * const, const double * const);
    //       void static normalize_set(double * const,
    //                                 const unsigned int, const unsigned int);
    //       void static orthogonalize_set(double * const,
    //                                     const unsigned int, const unsigned int,
    //                                     const bool);
    //       void static rests(const unsigned int, unsigned int * const,
    //                         unsigned int * const, const unsigned int);

    //       void static retrieve_hessian_indices(const unsigned int,
    //                                            unsigned int * const,
    //                                            unsigned int * const,
    //                                            unsigned int * const,
    //                                            unsigned int * const);
  };

  /* Make an instance of the format class */
  FCHKFormat theFCHKFormat;


  /*!
  **\brief Extract the data.
  **
  **The data are :
  **  comment
  **  charge
  **  multiplicity
  **  geometry
  **  connection table
  */
  bool FCHKFormat::ReadMolecule(OBBase * pOb, OBConversion * pConv)
  {
    OBMol * const pmol = dynamic_cast<OBMol*>(pOb);

    if (!pmol)
      return false;

    /* input stream */
    istream * const pifs = pConv->GetInStream();

    /* for logging errors */
    stringstream  error_msg;

    /* help variables */
    char buff[BUFF_SIZE];
    unsigned int lineno = 0;
    int  intvar;
    bool finished;
    bool atomnos_found = false;
    bool coords_found = false;
    bool nbond_found = false;
    bool ibond_found = false;
    bool hessian_found = false;
    bool dipder_found = false;
    bool alphaorb_found = false;
    bool betaorb_found = false;

    /* variables for storing the read data */
    int  Natoms = -1, MxBond = -1;
    int  numAOrb = -1, numBOrb = -1;
    int  numAElec = -1, numBElec = -1;
    vector<int>    atomnos, NBond, IBond;
    vector<double> coords, hessian, dipder, alphaorb, betaorb;

    /*** start reading ***/
    pmol->BeginModify();

    /* first line : comment */
    if (!pifs->getline(buff, BUFF_SIZE))
      {
        obErrorLog.ThrowError("FCHKFormat::ReadMolecule()",
                              "Failed to read the comment.",
                              obError);
        return false;
      }

    ++lineno;
    pmol->SetTitle(buff);

    /* extract all the information here */
    while (pifs->getline(buff, BUFF_SIZE))
      {
        ++lineno;

        /* Number of atoms */
        if (buff == strstr(buff, "Number of atoms"))
          {
            if (!FCHKFormat::read_int(buff, &Natoms))
              {
                error_msg << "Could not read the number of atoms from line #"
                          << lineno << ".";

                obErrorLog.ThrowError("FCHKFormat::ReadMolecule()",
                                      error_msg.str(),
                                      obError);
                return false;
              }
            if (0 >= Natoms)
              {
                error_msg << "Invalid number of atoms : " << Natoms << ".";
                obErrorLog.ThrowError("FCHKFormat::ReadMolecule()",
                                      error_msg.str(),
                                      obError);
                return false;
              }
            continue;
          }

        /* Dipole Moment */
        if (buff == strstr(buff, "Dipole Moment"))
          {
            if (!pifs->getline(buff, BUFF_SIZE)) {
              error_msg << "Could not read the dipole moment after line #"
                        << lineno << ".";
              obErrorLog.ThrowError("FCHKFormat::ReadMolecule()",
                                    error_msg.str(),
                                    obError);
              return false;
            }
            ++lineno;
            vector<double> moments;
            if (!FCHKFormat::read_numbers(buff, moments) || moments.size() != 3) {
              error_msg << "Could not read the dipole moment from line #"
                        << lineno << ".";
              obErrorLog.ThrowError("FCHKFormat::ReadMolecule()",
                                    error_msg.str(),
                                    obError);
              return false;
            }

            OBVectorData *dipoleMoment = new OBVectorData;
            dipoleMoment->SetAttribute("Dipole Moment");
            dipoleMoment->SetData(moments[0], moments[1], moments[2]);
            dipoleMoment->SetOrigin(fileformatInput);
            pmol->SetData(dipoleMoment);
            continue;
          }

        /* Charge */
        if (buff == strstr(buff, "Charge"))
          {
            if (!FCHKFormat::read_int(buff, &intvar))
              {
                error_msg << "Could not read the charge from line #"
                          << lineno << ".";
                obErrorLog.ThrowError("FCHKFormat::ReadMolecule()",
                                      error_msg.str(),
                                      obError);
                return false;
              }

            pmol->SetTotalCharge(intvar);
            continue;
          }

        /* Multiplicity */
        if (buff == strstr(buff, "Multiplicity"))
          {
            if (!FCHKFormat::read_int(buff, &intvar))
              {
                error_msg << "Could not read the multiplicity from line #"
                          << lineno << ".";
                obErrorLog.ThrowError("FCHKFormat::ReadMolecule()",
                                      error_msg.str(),
                                      obError);
                return false;
              }

            pmol->SetTotalSpinMultiplicity(intvar);
            continue;
          }

        /* start of the "Atomic numbers" section */
        if (buff == strstr(buff, "Atomic numbers"))
          {
            if (!FCHKFormat::validate_number(Natoms, "number of atoms", lineno))
              return false;

            if (!FCHKFormat::validate_section(buff,
                                              Natoms,
                                              "number of atomic numbers",
                                              lineno))
              return false;

            atomnos_found = true;
            continue;
          }

        /* start of the "Current cartesian coordinates" section */
        if (buff == strstr(buff, "Current cartesian coordinates"))
          {
            if (!FCHKFormat::validate_number(Natoms, "number of atoms", lineno))
              return false;

            if (!FCHKFormat::validate_section(buff,
                                              3 * Natoms,
                                              "number of coordinates",
                                              lineno))
              return false;

            coords_found = true;
            continue;
          }

        /* MxBond - largest number of bonds to an atom */
        if (buff == strstr(buff, "MxBond"))
          {
            if (!FCHKFormat::read_int(buff, &MxBond))
              {
                error_msg << "Could not read the number of bonds to each atom"
                          << " (MxBond) from line #" << lineno << ".";

                obErrorLog.ThrowError("FCHKFormat::ReadMolecule()",
                                      error_msg.str(),
                                      obError);
                return false;
              }
            if (0 >= MxBond)
              {
                error_msg << "Invalid number of bonds to each atom : "
                          << MxBond << ".";
                obErrorLog.ThrowError("FCHKFormat::ReadMolecule()",
                                      error_msg.str(),
                                      obError);
                return false;
              }
            continue;
          }

        /* start of the "NBond" section */
        if (buff == strstr(buff, "NBond"))
          {
            if (!FCHKFormat::validate_number(Natoms, "number of atoms", lineno))
              return false;

            if (!FCHKFormat::validate_number(MxBond,
                                             "number of bonds to each atom",
                                             lineno))
              return false;

            if (!FCHKFormat::validate_section(buff,
                                              Natoms,
                                              "number of bonds to each atom",
                                              lineno))
              return false;

            nbond_found = true;
            continue;
          }

        /* start of the "IBond" section */
        if (buff == strstr(buff, "IBond"))
          {
            if (!FCHKFormat::validate_number(Natoms, "number of atoms", lineno))
              return false;

            if (!FCHKFormat::validate_number(MxBond,
                                             "number of bonds to each atom",
                                             lineno))
              return false;

            if (!FCHKFormat::validate_section(buff,
                                              MxBond * Natoms,
                                              "number of bonds",
                                              lineno))
              return false;

            ibond_found = true;
            continue;
          }

        /* start of the "Cartesian Force Constants" section */
        if (buff == strstr(buff, "Cartesian Force Constants"))
          {
            if (!FCHKFormat::validate_number(Natoms, "number of atoms", lineno))
              return false;

            if (!FCHKFormat::validate_section(buff,
                                              (3 * Natoms) * (3 * Natoms + 1) / 2,
                                              "number of force constants",
                                              lineno))
              return false;

            hessian_found = true;
            continue;
          }

        /* start of the "Dipole Derivatives" section */
        if (buff == strstr(buff, "Dipole Derivatives"))
          {
            if (!FCHKFormat::validate_number(Natoms, "number of atoms", lineno))
              return false;

            if (!FCHKFormat::validate_section(buff,
                                              9 * Natoms,
                                              "number of dipole derivatives",
                                              lineno))
              return false;

            dipder_found = true;
            continue;
          }

        if (buff == strstr(buff, "Number of alpha electrons"))
          {
            if (!FCHKFormat::read_int(buff, &numAElec))
              {
                error_msg << "Could not read the number of alpha electrons"
                          << " from line #" << lineno << ".";

                obErrorLog.ThrowError("FCHKFormat::ReadMolecule()",
                                      error_msg.str(),
                                      obError);
                return false;
              }
            if (0 >= numAElec)
              {
                error_msg << "Invalid number of alpha electrons: "
                          << MxBond << ".";
                obErrorLog.ThrowError("FCHKFormat::ReadMolecule()",
                                      error_msg.str(),
                                      obError);
                return false;
              }
            continue;
          }

        if (buff == strstr(buff, "Number of beta electrons"))
          {
            if (!FCHKFormat::read_int(buff, &numBElec))
              {
                error_msg << "Could not read the number of beta electrons"
                          << " from line #" << lineno << ".";

                obErrorLog.ThrowError("FCHKFormat::ReadMolecule()",
                                      error_msg.str(),
                                      obError);
                return false;
              }
            if (0 >= numBElec)
              {
                error_msg << "Invalid number of beta electrons: "
                          << MxBond << ".";
                obErrorLog.ThrowError("FCHKFormat::ReadMolecule()",
                                      error_msg.str(),
                                      obError);
                return false;
              }
            continue;
          }

        /* start of the alpha orbitals section */
        if (buff == strstr(buff, "Alpha Orbital Energies"))
          {
            if (!FCHKFormat::read_int(buff, &numAOrb))
              {
                error_msg << "Could not read the number of alpha orbitals"
                          << " from line #" << lineno << ".";

                obErrorLog.ThrowError("FCHKFormat::ReadMolecule()",
                                      error_msg.str(),
                                      obError);
                return false;
              }
            if (0 >= numAOrb)
              {
                error_msg << "Invalid number of alpha orbitals: "
                          << MxBond << ".";
                obErrorLog.ThrowError("FCHKFormat::ReadMolecule()",
                                      error_msg.str(),
                                      obError);
                return false;
              }

            alphaorb_found = true;
            continue;
          }

        /* start of the beta orbitals section */
        if (buff == strstr(buff, "Beta Orbital Energies"))
          {
            if (!FCHKFormat::read_int(buff, &numBOrb))
              {
                error_msg << "Could not read the number of beta orbitals"
                          << " from line #" << lineno << ".";

                obErrorLog.ThrowError("FCHKFormat::ReadMolecule()",
                                      error_msg.str(),
                                      obError);
                return false;
              }
            if (0 >= numBOrb)
              {
                error_msg << "Invalid number of beta orbitals: "
                          << MxBond << ".";
                obErrorLog.ThrowError("FCHKFormat::ReadMolecule()",
                                      error_msg.str(),
                                      obError);
                return false;
              }

            betaorb_found = true;
            continue;
          }

        /* reading the atomic numbers */
        if (atomnos_found && Natoms > (int)atomnos.size())
          {
            if (!FCHKFormat::read_section(buff, atomnos, Natoms, &finished,
                                          "atomic numbers", lineno))
              return false;

            atomnos_found = !finished;
            continue;
          }

        /* reading the coordinates */
        if (coords_found && 3 * Natoms > (int)coords.size())
          {
            if (!FCHKFormat::read_section(buff, coords, 3 * Natoms, &finished,
                                          "coordinates", lineno, 16))
              return false;

            coords_found = !finished;
            continue;
          }

        /* reading the numbers of bonds to each atom */
        if (nbond_found && Natoms > (int)NBond.size())
          {
            if (!FCHKFormat::read_section(buff, NBond, Natoms, &finished,
                                          "number of bonds to each atom", lineno))
              return false;

            nbond_found = !finished;
            continue;
          }

        /* reading the list of atoms bound to each atom */
        if (ibond_found && MxBond * Natoms > (int)IBond.size())
          {
            if (!FCHKFormat::read_section(buff, IBond, MxBond * Natoms, &finished,
                                          "atom bonds", lineno))
              return false;

            ibond_found = !finished;
            continue;
          }

        /* reading the hessian */
        if (hessian_found &&
            (3 * Natoms) * (3 * Natoms + 1) / 2 > (int)hessian.size())
          {
            if (!FCHKFormat::read_section(buff, hessian,
                                          (3 * Natoms) * (3 * Natoms + 1) / 2,
                                          &finished,
                                          "hessian", lineno))
              return false;

            hessian_found = !finished;
            continue;
          }

        /* reading the dipole derivatives */
        if (dipder_found && 9 * Natoms > (int)dipder.size())
          {
            if (!FCHKFormat::read_section(buff, dipder,
                                          9 * Natoms,
                                          &finished,
                                          "dipole derivatives", lineno))
              return false;

            dipder_found = !finished;
            continue;
          }

        /* reading the alpha orbital energies */
        if (alphaorb_found && (numAOrb > alphaorb.size()))
          {
            if (!FCHKFormat::read_section(buff, alphaorb,
                                          numAOrb,
                                          &finished,
                                          "alpha orbital energies", lineno))
              return false;

            alphaorb_found = !finished; // if we've read once, don't re-read
            continue;
          }

        /* reading the beta orbital energies */
        if (betaorb_found &&  (numBOrb > betaorb.size()))
          {
            if (!FCHKFormat::read_section(buff, betaorb,
                                          numBOrb,
                                          &finished,
                                          "beta orbital energies", lineno))
              return false;

            betaorb_found = !finished;
            continue;
          }

      } /* while */

    /*** validating ***/
    if (0 >= Natoms)
      {
        error_msg << "Number of atoms could not be read.";
        obErrorLog.ThrowError("FCHKFormat::ReadMolecule()",
                              error_msg.str(),
                              obError);
        return false;
      }

    if (atomnos.empty())
      {
        error_msg << "\"Atomic numbers\" section was not found.";
        obErrorLog.ThrowError("FCHKFormat::ReadMolecule()",
                              error_msg.str(),
                              obError);
        return false;
      }

    if (coords.empty())
      {
        error_msg << "\"Current cartesian coordinates\" section was not found.";
        obErrorLog.ThrowError("FCHKFormat::ReadMolecule()",
                              error_msg.str(),
                              obError);
        return false;
      }

    /* connectivity : if MxBond is set, all the sections must be supplied */
    if (-1 != MxBond)
      {
        if (NBond.empty() || IBond.empty())
          {
            error_msg << "If MxBond is set, then the \"NBond\" and \"IBond\""
                      << " sections *must* be supplied.";
            obErrorLog.ThrowError("FCHKFormat::ReadMolecule()",
                                  error_msg.str(),
                                  obError);
            return false;
          }

        /* not nbonds > MxBond or <= 0
           no atom numbers < 0 or > Natoms */
        if (NBond.end() != find_if(NBond.begin(),
                                   NBond.end(),
                                   bind2nd(less_equal<int>(), 0)) ||
            NBond.end() != find_if(NBond.begin(),
                                   NBond.end(),
                                   bind2nd(greater<int>(), MxBond)) ||
            IBond.end() != find_if(IBond.begin(),
                                   IBond.end(),
                                   bind2nd(less<int>(), 0)) ||
            IBond.end() != find_if(IBond.begin(),
                                   IBond.end(),
                                   bind2nd(greater<int>(), Natoms)))
          {
            error_msg << "Invalid connectivity : check the \"NBond\" and/or"
                      << " \"IBond\" section(s).";
            obErrorLog.ThrowError("FCHKFormat::ReadMolecule()",
                                  error_msg.str(),
                                  obWarning);
            MxBond = -1;
          }
      }

    /*** finalizing ***/
    FCHKFormat::construct_mol(pmol, pConv, Natoms, atomnos, coords,
                              MxBond, NBond, IBond);

    //    if (!hessian.empty() && !pConv->IsOption("v", OBConversion::INOPTIONS))
    //    {
    //      FCHKFormat::vibana(pmol, hessian, dipder);
    //    }

    /*** finished ***/
    pmol->EndModify();

    if (numAOrb > 0 && alphaorb.size() == numAOrb) {
      OBOrbitalData *od = new OBOrbitalData; // create new store
      vector<string> symmetries;

      if (numBOrb > 0 && betaorb.size() == numBOrb) {  // open shell calculation
        od->LoadAlphaOrbitals(alphaorb, symmetries, numAElec);
        od->LoadBetaOrbitals(betaorb, symmetries, numBElec);
      } else {
        od->LoadClosedShellOrbitals(alphaorb, symmetries, numAElec);
      }
      od->SetOrigin(fileformatInput);
      pmol->SetData(od);
    }

    return true;
  }

  /*!
  **\brief Convert the last token in a string to integer
  **
  **Return false if the conversion failed.
  **The read integer is stored in the num argument.
  */
  bool FCHKFormat::read_int(const char * const line, int * const num)
  {
    vector<string>  vs;
    char * endptr;

    tokenize(vs, line);
    *num = strtol(vs.back().c_str(), &endptr, 10);

    return endptr != vs.back().c_str();
  }

  /*!
  **\brief Read integers or doubles into a vector from a string
  **\param line string
  **\param v    vector to read in
  */
  template<class T>
  bool FCHKFormat::read_numbers(const char * const line, vector<T> & v, const unsigned int width)
  {
    if (width == 0) {
      vector<string>  vs;
      tokenize(vs, line);

      if (0 < vs.size())
        {
          char * endptr;
          T num;

          for (vector<string>::const_iterator it=vs.begin(); vs.end() != it; ++it)
            {
              if (typeid(double) == typeid(T))
                num = static_cast<T>(strtod((*it).c_str(), &endptr));
              else
                num = static_cast<T>(strtol((*it).c_str(), &endptr, 10));

              /* quit if the value cannot be read */
              if (endptr == (*it).c_str())
                return false;

              v.push_back(num);
            }
        }
    } else {
      // fixed width records (e.g., coordinates = 16 characters)
      string lineCopy(line), substring;
      T num;
      char * endptr;
      int maxColumns = 80 / width;

      for (int column = 0; column < maxColumns; ++column) {
        substring = lineCopy.substr(column * width, width);

        if (typeid(double) == typeid(T))
          num = static_cast<T>(strtod(substring.c_str(), &endptr));
        else
          num = static_cast<T>(strtol(substring.c_str(), &endptr, 10));

        /* quit if the value cannot be read */
        if (endptr == substring.c_str())
          break;

        v.push_back(num);
      }
    }

    return true;
  }

  /*!
  **\brief Validate a section
  **
  **Try to convert the last token of a string to integer and compare it to a
  **given number. If they are not equal report an error and return false.
  **Otherwise return true.
  */
  bool FCHKFormat::validate_section(const char * const line,
                                    const int nreq,
                                    const char * const desc,
                                    const unsigned int lineno)
  {
    int intvar;
    stringstream error_msg;

    if (!FCHKFormat::read_int(line, &intvar))
      {
        error_msg << "Could not read the " << desc
                  << " from line #" << lineno << ".";
        obErrorLog.ThrowError("FCHKFormat::ReadMolecule()",
                              error_msg.str(),
                              obError);
        return false;
      }
    if (intvar != nreq)
      {
        error_msg << desc << " must be exactly " << nreq
                  << ", found " << intvar << ".";
        obErrorLog.ThrowError("FCHKFormat::ReadMolecule()",
                              error_msg.str(),
                              obError);
        return false;
      }

    return true;
  }

  /*!
  **\brief Validate a number
  **
  **Compare it with a given value and report an error unless they are equal.
  */
  bool FCHKFormat::validate_number(const int num,
                                   const char * const desc,
                                   const unsigned int lineno)
  {
    stringstream error_msg;

    if (-1 == num)
      {
        error_msg << desc   << " must be already read before line #"
                  << lineno << ".";
        obErrorLog.ThrowError("FCHKFormat::ReadMolecule()",
                              error_msg.str(),
                              obError);
        return false;
      }
    return true;
  }

  /*!
  **\brief Read numbers in a section
  **\param line     line
  **\param v        vector to read in
  **\param nreq     maximum number of elements in v
  **\param finished whether the reading is finished
  **\param desc     name of the section
  **\param lineno   line number
  */
  template <class T>
  bool FCHKFormat::read_section(const char * const line,
                                vector<T> & v,
                                const unsigned int nreq,
                                bool * const finished,
                                const char * const desc,
                                const unsigned int lineno,
                                const unsigned int width)
  {
    stringstream error_msg;

    *finished = false;

    if(!FCHKFormat::read_numbers(line, v, width))
      {
        error_msg << "Expecting " << desc << " in line #" << lineno << ".";
        obErrorLog.ThrowError("FCHKFormat::ReadMolecule()",
                              error_msg.str(),
                              obError);
        return false;
      }
    if (nreq <= v.size())
      {
        *finished = true;
        if (nreq < v.size())
          {
            error_msg << "Ignoring the superfluous " << desc
                      << "in line #" << lineno << ".";
            obErrorLog.ThrowError("FCHKFormat::ReadMolecule()",
                                  error_msg.str(),
                                  obWarning);
          }
      }

    return true;
  }

  /*!
  **\brief Construct the molecule
  **\param pmol    molecule
  **\param Natoms  number of atoms
  **\param atomnos atomic numbers
  **\param coords  coordinates
  **\param MxBond  largest number of bonds to an atom (usually 4)
  **\param NBond   number of bonds to each atom
  **\param IBond   list of atoms bounded to each atom
  */
  void FCHKFormat::construct_mol(OBMol * const pmol,
                                 OBConversion * const pConv,
                                 const unsigned int Natoms,
                                 const vector<int> & atomnos,
                                 const vector<double> & coords,
                                 const int MxBond,
                                 const vector<int> & NBond,
                                 const vector<int> & IBond)
  {
    pmol->ReserveAtoms(Natoms);

    OBAtom * atom;
    for (unsigned int a = 0; Natoms > a; ++a)
      {
        atom = pmol->NewAtom();

        atom->SetAtomicNum(atomnos[a]);
        atom->SetVector(coords[0 + 3 * a] * BOHR2ANGSTROM,
                        coords[1 + 3 * a] * BOHR2ANGSTROM,
                        coords[2 + 3 * a] * BOHR2ANGSTROM);
      }

    /* unless suppressed */
    if (!pConv->IsOption("b", OBConversion::INOPTIONS))
      {
        /* if the connectivity is provided */
        if (-1 != MxBond)
          {
            /* atoms */
            for (unsigned int a = 0; Natoms > a; ++a)
              {
                /* bonds for an atom */
                for (unsigned int i = 0; (unsigned int)NBond[a] > i; ++i)
                  {
                    pmol->AddBond(1 + a, IBond[MxBond*a + i], 1);
                  }
              }
          }
        else
          {
            pmol->ConnectTheDots();
          }
      }

    if (!pConv->IsOption("s", OBConversion::INOPTIONS) &&
        !pConv->IsOption("b", OBConversion::INOPTIONS))
      pmol->PerceiveBondOrders();
  }

} /* namespace OpenBabel */
