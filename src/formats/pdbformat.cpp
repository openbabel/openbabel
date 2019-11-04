/**********************************************************************
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2003-2006 Geoffrey R. Hutchison
Some portions Copyright (C) 2004 by Chris Morley

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
#include <openbabel/obfunctions.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/elements.h>
#include <openbabel/generic.h>
#include <openbabel/data.h>

#include <vector>
#include <map>
#include <cstdlib>
#include <algorithm>

#include <sstream>

using namespace std;
namespace OpenBabel
{
  class PDBFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    PDBFormat()
    { 
      OBConversion::RegisterFormat("pdb",this, "chemical/x-pdb");
      OBConversion::RegisterFormat("ent",this, "chemical/x-pdb");

      OBConversion::RegisterOptionParam("s", this, 0, OBConversion::INOPTIONS);
      OBConversion::RegisterOptionParam("b", this, 0, OBConversion::INOPTIONS);
      OBConversion::RegisterOptionParam("c", this, 0, OBConversion::INOPTIONS);

      OBConversion::RegisterOptionParam("o", this, 0, OBConversion::OUTOPTIONS);
      OBConversion::RegisterOptionParam("n", this, 0, OBConversion::OUTOPTIONS);
    }

    virtual const char* Description() //required
    {
      return
        "Protein Data Bank format\n"
        "Read Options e.g. -as\n"
        "  s  Output single bonds only\n"
        "  b  Disable bonding entirely\n"
        "  c  Ignore CONECT records\n\n"

        "Write Options, e.g. -xo\n"
        "  n  Do not write duplicate CONECT records to indicate bond order\n"
        "  o  Write origin in space group label (CRYST1 section)\n\n";
    };

    virtual const char* SpecificationURL()
    { return "http://www.wwpdb.org/docs.html";};

    virtual const char* GetMIMEType()
    { return "chemical/x-pdb"; };

    //*** This section identical for most OBMol conversions ***
    ////////////////////////////////////////////////////
    /// The "API" interface functions
  	virtual int SkipObjects(int n, OBConversion* pConv);
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

  };
  //***

  //Make an instance of the format class
  PDBFormat thePDBFormat;

  ////////////////////////////////////////////////////
  /// Utility functions
  static void fixRhombohedralSpaceGroupWriter(string &strHM);
  static void fixRhombohedralSpaceGroupReader(string &strHM);
  static bool parseAtomRecord(char *buffer, OBMol & mol, int chainNum);
  static bool parseConectRecord(char *buffer, OBMol & mol);
  static bool readIntegerFromRecord(char *buffer, unsigned int columnAsSpecifiedInPDB, long int *target);

  //extern OBResidueData    resdat; now in mol.h

  /////////////////////////////////////////////////////////////////
 	int PDBFormat::SkipObjects(int n, OBConversion* pConv)
  {
    if (n == 0)
      ++ n;
    istream &ifs = *pConv->GetInStream();
    char buffer[BUFF_SIZE];
    while (n && ifs.getline(buffer,BUFF_SIZE))
      {
        if (EQn(buffer,"ENDMDL",6))
          -- n;
      }

    return ifs.good() ? 1 : -1;
  }
  /////////////////////////////////////////////////////////////////
   template <typename T> string to_string(T pNumber)
  {
    ostringstream oOStrStream;
    oOStrStream << pNumber;
    return oOStrStream.str();
  }

  /////////////////////////////////////////////////////////////////
  bool PDBFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title = pConv->GetTitle();

    int chainNum = 1;
    char buffer[BUFF_SIZE] = {0,};
    string line, key, value;
    OBPairData *dp;

    mol.SetTitle(title);
    // We need to prevent chains perception routines from running while
    // we are adding residues from the PDB file
    mol.SetChainsPerceived();

    mol.BeginModify();
    bool ateend = false;
    while (ifs.good() && ifs.getline(buffer,BUFF_SIZE))
      {
        if (EQn(buffer,"ENDMDL",6)) {
          ateend = true;
          break;
        }
        if (EQn(buffer,"END",3)) {
          // eat anything until the next ENDMDL
          while (ifs.getline(buffer,BUFF_SIZE) && !EQn(buffer,"ENDMDL",6));
          ateend = true;
          break;
        }
        if (EQn(buffer,"TER",3)) {
          chainNum++;
          continue;
        }
        if (EQn(buffer,"ATOM",4) || EQn(buffer,"HETATM",6))
          {
            if( ! parseAtomRecord(buffer,mol,chainNum))
              {
                stringstream errorMsg;
                errorMsg << "WARNING: Problems reading a PDB file\n"
                         << "  Problems reading a ATOM/HETATM record.\n";
                obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obError);
              }
            continue;
          }

        if (EQn(buffer,"CONECT",6)) {
          // Don't parse a CONECT record if the user tells us to ignore them
          if (!pConv->IsOption("c",OBConversion::INOPTIONS)) {
            parseConectRecord(buffer,mol);
            continue;
          }
        }

        // crystal cells
        if (EQn(buffer,"CRYST1",6)) {
          float a, b, c, alpha, beta, gamma;
          string group = "";

          sscanf (&(buffer[6]), "%9f%9f%9f%7f%7f%7f", &a, &b, &c,
                  &alpha, &beta, &gamma);
          buffer[66] = '\0';
          group += &(buffer[55]);
          Trim (group);
          fixRhombohedralSpaceGroupReader(group);

          OBUnitCell *pCell=new OBUnitCell;
          pCell->SetOrigin(fileformatInput);
          pCell->SetData(a,b,c,alpha,beta,gamma);
          pCell->SetSpaceGroup(group);
          pmol->SetData(pCell);
          continue;
        }

        // another record type, add it as an OBPairData entry
        line = buffer;
        // if the file is valid, all lines should have more than 6 characters
        if (line.length() < 6)
          {
            stringstream errorMsg;
            errorMsg << "ERROR: not a valid PDB file" << endl;
            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obError);
            return false;
          }
        key = line.substr(0,6); // the first 6 characters are the record name
        Trim(key);
        value = line.substr(6);

        // We haven't found this record yet
        if (!mol.HasData(key)) {
          dp = new OBPairData;
          dp->SetAttribute(key);
          dp->SetValue(value);
          dp->SetOrigin(fileformatInput);
          mol.SetData(dp);
        }
        // Add on additional lines
        else {
          dp = static_cast<OBPairData*>(mol.GetData(key));
          line = dp->GetValue();
          line += '\n';
          line += value;
          dp->SetValue(line);
        }
      }

    if (!mol.NumAtoms()) { // skip the rest of this processing
      mol.EndModify();
      return ateend; //explictly empty molecules are not invalid
    }

    resdat.AssignBonds(mol);
    /*assign hetatm bonds based on distance*/

    mol.EndModify();
    // Clear all virtual bond data
    vector<OBGenericData*> vbonds = mol.GetAllData(OBGenericDataType::VirtualBondData);
    mol.DeleteData(vbonds);

    if (!pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.ConnectTheDots();

    if (!pConv->IsOption("s",OBConversion::INOPTIONS) && !pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();

    // EndModify() blows away the chains perception flag so we set it again here
    mol.SetChainsPerceived();

    // Guess how many hydrogens are present on each atom based on typical valencies
    FOR_ATOMS_OF_MOL(matom, mol)
      OBAtomAssignTypicalImplicitHydrogens(&*matom);

    // clean out remaining blank lines
    std::streampos ipos;
    do
    {
      ipos = ifs.tellg();
      ifs.getline(buffer,BUFF_SIZE);
    }
    while(strlen(buffer) == 0 && !ifs.eof() );
    ifs.seekg(ipos);

    return(true);
  }

  /////////////////////////////////////////////////////////////////////////
  //! Utility function to read a 5-digit integer starting from a specified column
  /*! This function reads a 5-digit integer, starting from column
    columnAsSpecifiedInPDB from the buffer, converts it to a long
    integer, and returns either false or true, if the conversion was
    successful or not. If the conversion was not successful, the target
    is set to a random value.

    For instance, the PDB Format Description for a CONECT record specifies

    COLUMNS        DATA TYPE        FIELD           DEFINITION
    ---------------------------------------------------------------------------------
    1 -  6         Record name      "CONECT"
    7 - 11         Integer          serial          Atom serial number
    ...

    To read the Atom serial number, you would call

    long int target;
    if ( readIntegerFromRecord(buffer, 7, &target) == false ) {
    cerr << "Could not parse" << endl;
    }

    This function does not check the length of the buffer, or
    strlen(buffer). If the buffer is not long enough => SEGFAULT.
  */
  static bool readIntegerFromRecord(char *buffer, unsigned int columnAsSpecifiedInPDB, long int *target)
  {
    char integerBuffer[6];
    integerBuffer[5] = '\0';

    strncpy(integerBuffer, buffer+columnAsSpecifiedInPDB-1, 5);

    char *errorCheckingEndPtr;
    *target = strtol(integerBuffer, &errorCheckingEndPtr, 10);
    if (integerBuffer == errorCheckingEndPtr)
      return(false);
    return(true);
  }

  //! Read a CONECT record
  /*! This function reads a CONECT record, as specified
    http://www.rcsb.org/pdb/docs/format/pdbguide2.2/guide2.2_frame.html,
    in short:

    COLUMNS         DATA TYPE        FIELD           DEFINITION
    ---------------------------------------------------------------------------------
    1 -  6         Record name      "CONECT"
    7 - 11         Integer          serial          Atom serial number
    12 - 16         Integer          serial          Serial number of bonded atom
    17 - 21         Integer          serial          Serial number of bonded atom
    22 - 26         Integer          serial          Serial number of bonded atom
    27 - 31         Integer          serial          Serial number of bonded atom
    32 - 36         Integer          serial          Serial number of hydrogen bonded atom
    37 - 41         Integer          serial          Serial number of hydrogen bonded atom
    42 - 46         Integer          serial          Serial number of salt bridged atom
    47 - 51         Integer          serial          Serial number of hydrogen bonded atom
    52 - 56         Integer          serial          Serial number of hydrogen bonded atom
    57 - 61         Integer          serial          Serial number of salt bridged atom

    Hydrogen bonds and salt bridges are ignored. --Stefan Kebekus.
  */

  bool parseConectRecord(char *buffer,OBMol &mol)
  {
    stringstream errorMsg;
    string clearError;

    // Setup strings and string buffers
    vector<string> vs;
    buffer[70] = '\0';
    if (strlen(buffer) < 70)
      {
        errorMsg << "WARNING: Problems reading a PDB file\n"
                 << "  Problems reading a CONECT record.\n"
                 << "  According to the PDB specification,\n"
                 << "  the record should have 70 columns, but OpenBabel found "
                 << strlen(buffer) << " columns." << endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obInfo);
        errorMsg.str(clearError);
      }

    // Serial number of the first atom, read from column 7-11 of the
    // connect record, to which the other atoms connect to.
    long int startAtomSerialNumber;
    // A pointer to the first atom.
    OBAtom *firstAtom = NULL;
    // Serial numbers of the atoms which bind to firstAtom, read from
    // columns 12-16, 17-21, 22-27 and 27-31 of the connect record. Note
    // that we reserve space for 5 integers, but read only four of
    // them. This is to simplify the determination of the bond order;
    // see below.
    long int boundedAtomsSerialNumbers[5]  = {0,0,0,0,0};
    // Bools which tell us which of the serial numbers in
    // boundedAtomsSerialNumbers are read from the file, and which are
    // invalid
    bool boundedAtomsSerialNumbersValid[5] = {false, false, false, false, false};

    // Pragmatic approach -- too many non-standard PDB files out there
    // (including some old ones from us)
    // So if we have a small number of atoms, then try to break by spaces
    // Otherwise (i.e., NumAtoms() > 9,999 we need to go by position)
    // We'll switch back and forth a few times to save duplicating common code

    if (mol.NumAtoms() <= 9999)
      {
        // make sure we don't look at salt bridges or whatever, so cut the buffer short
        buffer[32] = '\0';
        tokenize(vs,buffer);
        if( vs.empty() || vs.size() < 2)
          return false;
        vs.erase(vs.begin()); // remove "CONECT"

        startAtomSerialNumber = atoi(vs[0].c_str());
      }
    else
      {
        if (readIntegerFromRecord(buffer, 7, &startAtomSerialNumber) == false)
          {
            errorMsg << "WARNING: Problems reading a PDB file\n"
                     << "  Problems reading a CONECT record.\n"
                     << "  According to the PDB specification,\n"
                     << "  columns 7-11 should contain the serial number of an atom.\n"
                     << "  THIS CONECT RECORD WILL BE IGNORED." << endl;
            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
            return(false);
          }
      }

    vector<OBAtom*>::iterator i;
    for (OBAtom *a1 = mol.BeginAtom(i);a1;a1 = mol.NextAtom(i)) {
      // atoms may not have residue information, but if they do,
      // check serial numbers
      if (a1->GetResidue() != NULL &&
          static_cast<long int>(a1->GetResidue()->
                                GetSerialNum(a1)) == startAtomSerialNumber)
        {
          firstAtom = a1;
          break;
        }
    }

    if (firstAtom == NULL)
      {
        errorMsg << "WARNING: Problems reading a PDB file:\n"
                 << "  Problems reading a CONECT record.\n"
                 << "  According to the PDB specification,\n"
                 << "  columns 7-11 should contain the serial number of an atom.\n"
                 << "  No atom was found with this serial number.\n"
                 << "  THIS CONECT RECORD WILL BE IGNORED." << endl;
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
        return(false);
      }

    if (mol.NumAtoms() < 9999)
      {
        if (vs.size() > 1) boundedAtomsSerialNumbers[0] = atoi(vs[1].c_str());
        if (vs.size() > 2) boundedAtomsSerialNumbers[1] = atoi(vs[2].c_str());
        if (vs.size() > 3) boundedAtomsSerialNumbers[2] = atoi(vs[3].c_str());
        if (vs.size() > 4) boundedAtomsSerialNumbers[3] = atoi(vs[4].c_str());

        unsigned int limit = 4;
        if (vs.size() <= 4)
          limit = vs.size() - 1;

        for (unsigned int s = 0; s < limit; ++s)
          boundedAtomsSerialNumbersValid[s] = true;
      }
    else
      {
        // Now read the serial numbers. If the first serial number is not
        // present, this connect record probably contains only hydrogen
        // bonds and salt bridges, which we ignore. In that case, we just
        // exit gracefully.
        boundedAtomsSerialNumbersValid[0] = readIntegerFromRecord(buffer, 12, boundedAtomsSerialNumbers+0);
        if (boundedAtomsSerialNumbersValid[0] == false)
          return(true);
        boundedAtomsSerialNumbersValid[1] = readIntegerFromRecord(buffer, 17, boundedAtomsSerialNumbers+1);
        boundedAtomsSerialNumbersValid[2] = readIntegerFromRecord(buffer, 22, boundedAtomsSerialNumbers+2);
        boundedAtomsSerialNumbersValid[3] = readIntegerFromRecord(buffer, 27, boundedAtomsSerialNumbers+3);
      }

    // Now iterate over the VALID boundedAtomsSerialNumbers and connect
    // the atoms.
    for(unsigned int k=0; boundedAtomsSerialNumbersValid[k]; k++)
      {
        // Find atom that is connected to, write an error message
        OBAtom *connectedAtom = 0L;
        for (OBAtom *a1 = mol.BeginAtom(i);a1;a1 = mol.NextAtom(i)) {
          // again, atoms may not have residues, but if they do, check serials
          if (a1->GetResidue() != NULL &&
              static_cast<long int>(a1->GetResidue()->
                                    GetSerialNum(a1)) == boundedAtomsSerialNumbers[k])
            {
              connectedAtom = a1;
              break;
            }
        }
        if (connectedAtom == 0L)
          {
            errorMsg << "WARNING: Problems reading a PDB file:\n"
                     << "  Problems reading a CONECT record.\n"
                     << "  According to the PDB specification,\n"
                     << "  Atoms with serial #" << startAtomSerialNumber
                     << " and #" << boundedAtomsSerialNumbers[k]
                     << " should be connected\n"
                     << "  However, an atom with serial #" << boundedAtomsSerialNumbers[k] << " was not found.\n"
                     << "  THIS CONECT RECORD WILL BE IGNORED." << endl;
            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str() , obWarning);
            return(false);
          }

        // Figure the bond order
        unsigned char order = 0;
        while(boundedAtomsSerialNumbersValid[k+order+1] && (boundedAtomsSerialNumbers[k+order]
                                                            == boundedAtomsSerialNumbers[k+order+1]))
          order++;
        k += order;

        // Generate the bond
        if (firstAtom->GetIdx() < connectedAtom->GetIdx()) { // record the bond 'in one direction' only
          OBBond *bond = mol.GetBond(firstAtom, connectedAtom);
          if (!bond)
            mol.AddBond(firstAtom->GetIdx(), connectedAtom->GetIdx(), order+1);
          else // An additional CONECT record with the same firstAtom that references
               // a bond created in the previous CONECT record.
               // For example, the 1136->1138 double bond in the following:
               //   CONECT 1136 1128 1137 1137 1138
               //   CONECT 1136 1138 1139
            bond->SetBondOrder(bond->GetBondOrder() + order+1);
        }

      }
    return(true);
  }

  //////////////////////////////////////////////////////////////////////////////
  bool PDBFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    unsigned int i;
    char buffer[BUFF_SIZE];
    char type_name[10], padded_name[10];
    char the_res[10];
    char the_chain = ' ';
    const char *element_name;
    int res_num;
    char the_insertioncode = ' ';
    bool het=true;
    int model_num = 0;
    const int MAX_HM_NAME_LEN = 11;

    if (!pConv->IsLast() || pConv->GetOutputIndex() > 1)
      { // More than one molecule record
        model_num = pConv->GetOutputIndex(); // MODEL 1-based index
        snprintf(buffer, BUFF_SIZE, "MODEL %8d", model_num);
        ofs << buffer << endl;
      }

    // write back all fields (REMARKS, HELIX, SHEET, SITE, ...)
    bool compndWritten = false;
    bool authorWritten = false;
    std::vector<OBGenericData*> pairData = mol.GetAllData(OBGenericDataType::PairData);
    for (std::vector<OBGenericData*>::iterator data = pairData.begin(); data != pairData.end(); ++data) {
      OBPairData *pd = static_cast<OBPairData*>(*data);
      string attr = pd->GetAttribute();

      // filter to make sure we are writing pdb fields only
      if (attr != "HEADER" && attr != "OBSLTE" && attr != "TITLE" && attr != "SPLIT" &&
          attr != "CAVEAT" && attr != "COMPND" && attr != "SOURCE" && attr != "KEYWDS" &&
          attr != "EXPDTA" && attr != "NUMMDL" && attr != "MDLTYP" && attr != "AUTHOR" &&
          attr != "REVDAT" && attr != "SPRSDE" && attr != "JRNL" && attr != "REMARK" &&
          attr != "DBREF" && attr != "DBREF1" && attr != "DBREF2" && attr != "SEQADV" &&
          attr != "SEQRES" && attr != "MODRES" && attr != "HET" && attr != "HETNAM" &&
          attr != "HETSYN" && attr != "FORMUL" && attr != "HELIX" && attr != "SHEET" &&
          attr != "SSBOND" && attr != "LINK" && attr != "CISPEP" && attr != "SITE" &&
          attr != "ORIGX1" && attr != "ORIGX2" && attr != "ORIGX3" && attr != "SCALE1" &&
          attr != "SCALE2" && attr != "SCALE3" && attr != "MATRIX1" && attr != "MATRIX2" &&
          attr != "MATRIX3" && attr != "MODEL")
        continue;

      if (attr == "COMPND")
        compndWritten = true;
      if (attr == "AUTHOR")
        authorWritten = true;

      // compute spacing needed. HELIX, SITE, HET, ... are trimmed when reading
      int nSpacing = 6 - attr.size();
      for (int i = 0; i < nSpacing; ++i)
        attr += " ";


      std::string lines = pd->GetValue();
      string::size_type last = 0;
      string::size_type pos = lines.find('\n');
      while (last != string::npos) {
        string line = lines.substr(last, pos - last);
        if (pos == string::npos)
          last = string::npos;
        else
          last = pos + 1;
        pos = lines.find('\n', last);

        ofs << attr << line << endl;
      }
    }

    if (!compndWritten) {
      if (strlen(mol.GetTitle()) > 0)
        snprintf(buffer, BUFF_SIZE, "COMPND    %s ",mol.GetTitle());
      else
        snprintf(buffer, BUFF_SIZE, "COMPND    UNNAMED");
      ofs << buffer << endl;
    }

    if (!authorWritten) {
      snprintf(buffer, BUFF_SIZE, "AUTHOR    GENERATED BY OPEN BABEL %s",BABEL_VERSION);
      ofs << buffer << endl;
    }

    // Write CRYST1 record, containing unit cell parameters, space group
    // and Z value (supposed to be 1)
    if (pmol->HasData(OBGenericDataType::UnitCell))
      {
        OBUnitCell *pUC = (OBUnitCell*)pmol->GetData(OBGenericDataType::UnitCell);
        if(pUC->GetSpaceGroup()){
          string tmpHM=pUC->GetSpaceGroup()->GetHMName();
          fixRhombohedralSpaceGroupWriter(tmpHM);

          // Do we have an extended HM symbol, with origin choice as ":1" or ":2" ? If so, remove it.
          size_t n=tmpHM.find(":");
          if(n!=string::npos) tmpHM=tmpHM.substr(0, n);

          if (pConv->IsOption("o", OBConversion::OUTOPTIONS))
            {
              unsigned int origin = pUC->GetSpaceGroup()->GetOriginAlternative();
              if (origin == pUC->GetSpaceGroup()->HEXAGONAL_ORIGIN)
                tmpHM[0] = 'H';
              else if (origin > 0)
                tmpHM += ":" + to_string(origin);

              if (tmpHM.length() > MAX_HM_NAME_LEN)
              {
                tmpHM.erase(std::remove(tmpHM.begin(), tmpHM.end(), ' '),
                            tmpHM.end());
              }
            }

          snprintf(buffer, BUFF_SIZE,
                   "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s 1",
                   pUC->GetA(), pUC->GetB(), pUC->GetC(),
                   pUC->GetAlpha(), pUC->GetBeta(), pUC->GetGamma(),
                   tmpHM.c_str());
        }
        else
          snprintf(buffer, BUFF_SIZE,
                   "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s 1",
                   pUC->GetA(), pUC->GetB(), pUC->GetC(),
                   pUC->GetAlpha(), pUC->GetBeta(), pUC->GetGamma(),
                   "P1");

        ofs << buffer << endl;
      }

    // before we write any records, we should check to see if any coord < -1000
    // which will cause errors in the formatting

    double minX, minY, minZ;
    minX = minY = minZ = -999.0f;
    FOR_ATOMS_OF_MOL(a, mol)
      {
        if (a->GetX() < minX)
          minX = a->GetX();
        if (a->GetY() < minY)
          minY = a->GetY();
        if (a->GetZ() < minZ)
          minZ = a->GetZ();
      }
    vector3 transV = VZero;
    if (minX < -999.0)
      transV.SetX(-1.0*minX - 900.0);
    if (minY < -999.0)
      transV.SetY(-1.0*minY - 900.0);
    if (minZ < -999.0)
      transV.SetZ(-1.0*minZ - 900.0);

    // if minX, minY, or minZ was never changed, shift will be 0.0f
    // otherwise, move enough so that smallest coord is > -999.0f
    mol.Translate(transV);

    OBAtom *atom;
    OBResidue *res;
    for (i = 1; i <= mol.NumAtoms(); i++)
      {
        atom = mol.GetAtom(i);
        strncpy(type_name, OBElements::GetSymbol(atom->GetAtomicNum()), sizeof(type_name));
        type_name[sizeof(type_name) - 1] = '\0';

        //two char. elements are on position 13 and 14 one char. start at 14
        if (strlen(type_name) > 1)
          type_name[1] = toupper(type_name[1]);
        else
          {
            char tmp[10];
            strncpy(tmp, type_name, 9); // make sure to null-terminate tmp
            snprintf(type_name, sizeof(type_name), " %-3s", tmp);
          }

        if ( (res = atom->GetResidue()) != 0 )
          {
            het = res->IsHetAtom(atom);
            snprintf(the_res,4,"%s",(char*)res->GetName().c_str());
            the_res[4] = '\0';
            snprintf(type_name,5,"%s",(char*)res->GetAtomID(atom).c_str());
            the_chain = res->GetChain();

            //two char. elements are on position 13 and 14 one char. start at 14
            if (strlen(OBElements::GetSymbol(atom->GetAtomicNum())) == 1)
              {
                if (strlen(type_name) < 4)
                  {
                    char tmp[10];
                    strncpy(tmp, type_name, 9); // make sure to null-terminate tmp
                    snprintf(padded_name, sizeof(padded_name), " %-3s", tmp);
                    strncpy(type_name,padded_name,4);
                    type_name[4] = '\0';
                  }
                else
                  {
                    /*
                      type_name[4] = type_name[3];
                      type_name[3] = type_name[2];
                      type_name[2] = type_name[1];
                      type_name[1] = type_name[0];
                      type_name[0] = type_name[4];
                    */
                    type_name[4] = '\0';
                  }
              }
            res_num = res->GetNum();
            the_insertioncode = res->GetInsertionCode();
            if (0 == the_insertioncode) the_insertioncode=' ';
          }
        else
          {
            strcpy(the_res,"UNK");
            the_res[3] = '\0';
            snprintf(padded_name,sizeof(padded_name), "%s",type_name);
            strncpy(type_name,padded_name,4);
            type_name[4] = '\0';
            res_num = 1;
            the_insertioncode=' ';
          }

        element_name = OBElements::GetSymbol(atom->GetAtomicNum());

        int charge = atom->GetFormalCharge();
        char scharge[3] = { ' ', ' ', '\0' };
        if(0 != charge)
          {
            snprintf(scharge, 3, "%+d", charge);
            char tmp = scharge[1];
            scharge[1] = scharge[0];
            scharge[0] = tmp;
          }

        double occup = 1.0;
        if (atom->HasData("_atom_site_occupancy"))
        {
         OBPairFloatingPoint *occup_fp = dynamic_cast<OBPairFloatingPoint*> (atom->GetData("_atom_site_occupancy"));
         occup = occup_fp->GetGenericValue();
        }

        snprintf(buffer, BUFF_SIZE, "%s%5d %-4s %-3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f  0.00          %2s%2s\n",
                 het?"HETATM":"ATOM  ",
                 i,
                 type_name,
                 the_res,
                 the_chain,
                 res_num,
                 the_insertioncode,
                 atom->GetX(),
                 atom->GetY(),
                 atom->GetZ(),
                 occup,
                 element_name,
                 scharge);
        ofs << buffer;
      }

    OBAtom *nbr;
    vector<OBBond*>::iterator k;
    for (i = 1; i <= mol.NumAtoms(); i ++)
      {
        atom = mol.GetAtom(i);
        if (atom->GetExplicitDegree() == 0)
          continue; // no need to write a CONECT record -- no bonds

        // Write out up to 4 real bonds per line PR#1711154
        int currentValence = 0;
        for (nbr = atom->BeginNbrAtom(k);nbr;nbr = atom->NextNbrAtom(k))
          {
            OBBond *bond = mol.GetBond(atom, nbr);
            if(!bond) continue;
            unsigned bondorder = bond->GetBondOrder();
            if(bondorder == 0 || pConv->IsOption("n", OBConversion::OUTOPTIONS)) 
              bondorder = 1;
            //a non-standard convention is to store bond orders by
            //replicating conect records
            for(unsigned bo = 0; bo < bondorder; bo++) {
              if ((currentValence % 4) == 0) {
                if (currentValence > 0) {
                  // Add the trailing space to finish the previous record
                  ofs << "                                       \n";
                }
                // write the start of a new CONECT record
                snprintf(buffer, BUFF_SIZE, "CONECT%5d", i);
                ofs << buffer;
              }
              currentValence++;
              snprintf(buffer, BUFF_SIZE, "%5d", nbr->GetIdx());
              ofs << buffer;
            }
          }

        // Add trailing spaces
        while ((currentValence % 4) != 0) {
          ofs << "     ";
          currentValence++;
        }
        ofs << "                                       \n";
      }

    snprintf(buffer, BUFF_SIZE, "MASTER        0    0    0    0    0    0    0    0 ");
    ofs << buffer;
    snprintf(buffer, BUFF_SIZE, "%4d    0 %4d    0\n",mol.NumAtoms(),mol.NumAtoms());
    ofs << buffer;

    if (model_num) {
      ofs << "ENDMDL" << endl;
	  if (pConv->IsLast()) {
	    ofs << "END\n";
	  }
    }
	else {
	  ofs << "END\n";
	}

    return(true);
  }

  ////////////////////////////////////////////////////////////////
  static void fixRhombohedralSpaceGroupWriter(string &strHM)
  {
    /* This is due to the requirment of PDB to name rhombohedral groups
       with H (http://deposit.rcsb.org/adit/docs/pdb_atom_format.html) */
    const int SIZE = 7;
    const char* groups[SIZE]  =   {"R 3:H",
                                   "R -3:H",
                                   "R 3 2:H",
                                   "R 3 m:H",
                                   "R 3 c:H",
                                   "R -3 m:H",
                                   "R -3 c:H"};

    std::vector<string> vec(groups, groups + SIZE);
    if(std::find(vec.begin(), vec.end(), strHM) != vec.end())
    {
      strHM[0] = 'H';
    }
  }

  static void fixRhombohedralSpaceGroupReader(string &strHM)
  {
    /* This is due to the requirment of PDB to name rhombohedral groups
       with H (http://deposit.rcsb.org/adit/docs/pdb_atom_format.html) */
    const int SIZE = 7;
    const char* groups[SIZE]  =   {"H 3",
                                   "H -3",
                                   "H 3 2",
                                   "H 3 m",
                                   "H 3 c",
                                   "H -3 m",
                                   "H -3 c"};

    std::vector<string> vec(groups, groups + SIZE);

    if(std::find(vec.begin(), vec.end(), strHM) != vec.end())
    {
      strHM[0] = 'R';
      strHM += ":H";
    }
  }

  /*

     From http://deposit.rcsb.org/adit/docs/pdb_atom_format.html

	COLUMNS        DATA TYPE       CONTENTS
	--------------------------------------------------------------------------------
	 1 -  6        Record name     "ATOM  "
	 7 - 11        Integer         Atom serial number.
	13 - 16        Atom            Atom name.
	17             Character       Alternate location indicator.
	18 - 20        Residue name    Residue name.
	22             Character       Chain identifier.
	23 - 26        Integer         Residue sequence number.
	27             AChar           Code for insertion of residues.
	31 - 38        Real(8.3)       Orthogonal coordinates for X in Angstroms.
	39 - 46        Real(8.3)       Orthogonal coordinates for Y in Angstroms.
	47 - 54        Real(8.3)       Orthogonal coordinates for Z in Angstroms.
	55 - 60        Real(6.2)       Occupancy.
	61 - 66        Real(6.2)       Temperature factor (Default = 0.0).
	73 - 76        LString(4)      Segment identifier, left-justified.
	77 - 78        LString(2)      Element symbol, right-justified.
	79 - 80        LString(2)      Charge on the atom.
  */
  static bool parseAtomRecord(char *buffer, OBMol &mol,int /*chainNum*/)
  /* ATOMFORMAT "(i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3,2f6.2,a2,a2)" */
  {
    string sbuf = &buffer[6];
    if (sbuf.size() < 48)
      return(false);

    bool hetatm = (EQn(buffer,"HETATM",6)) ? true : false;
    bool elementFound = false; // true if correct element found in col 77-78

    /* serial number */
    string serno = sbuf.substr(0,5);

    /* atom name */
    string atmid = sbuf.substr(6,4);

    /* chain */
    char chain = sbuf.substr(15,1)[0];

    /* insertion code */
    char insertioncode = sbuf.substr(27-6-1,1)[0];
    if (' '==insertioncode) insertioncode=0;
    /* element */
    string element = "  ";
    if (sbuf.size() > 71)
      {
        element = sbuf.substr(70,2);
        if (isalpha(element[1]))
          {
            if (element[0] == ' ')
              {
                element.erase(0, 1);
                elementFound = true;
              }
            else if (isalpha(element[0]))
              {
                elementFound = true;
                element[1] = tolower(element[1]);
              }
          }
      }

    if (!elementFound)
      {
        stringstream errorMsg;
        errorMsg << "WARNING: Problems reading a PDB file\n"
                 << "  Problems reading a HETATM or ATOM record.\n"
                 << "  According to the PDB specification,\n"
                 << "  columns 77-78 should contain the element symbol of an atom.\n"
                 << "  but OpenBabel found '" << element << "' (atom " << mol.NumAtoms()+1 << ")";
        obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
      }

    // charge - optional
    string scharge;
    if (sbuf.size() > 73)
      {
        scharge = sbuf.substr(72,2);
      }

    //trim spaces on the right and left sides
    while (!atmid.empty() && atmid[0] == ' ')
      atmid = atmid.erase(0, 1);

    while (!atmid.empty() && atmid[atmid.size()-1] == ' ')
      atmid = atmid.substr(0,atmid.size()-1);

    /* residue name */
    string resname = sbuf.substr(11,3);
    if (resname == "   ")
      resname = "UNK";
    else
      {
        while (!resname.empty() && resname[0] == ' ')
          resname = resname.substr(1,resname.size()-1);

        while (!resname.empty() && resname[resname.size()-1] == ' ')
          resname = resname.substr(0,resname.size()-1);
      }

    string type;
    if (!elementFound) {
      // OK, we have to fall back to determining the element from the atom type
      // This is unreliable, but there's no other choice
      if (EQn(buffer,"ATOM",4)) {
        type = atmid.substr(0,2);
        if (isdigit(type[0])) {
          // sometimes non-standard files have, e.g 11HH
          if (!isdigit(type[1])) type = atmid.substr(1,1);
          else type = atmid.substr(2,1);
        } else if ((sbuf[6] == ' ' &&
                   strncasecmp(type.c_str(), "Zn", 2) != 0 &&
                   strncasecmp(type.c_str(), "Fe", 2) != 0) ||
                   isdigit(type[1]))	//type[1] is digit in Platon
          type = atmid.substr(0,1);     // one-character element


        if (resname.substr(0,2) == "AS" || resname[0] == 'N') {
          if (atmid == "AD1")
            type = "O";
          if (atmid == "AD2")
            type = "N";
        }
        if (resname.substr(0,3) == "HIS" || resname[0] == 'H') {
          if (atmid == "AD1" || atmid == "AE2")
            type = "N";
          if (atmid == "AE1" || atmid == "AD2")
            type = "C";
        }
        if (resname.substr(0,2) == "GL" || resname[0] == 'Q') {
          if (atmid == "AE1")
            type = "O";
          if (atmid == "AE2")
            type = "N";
        }
        // fix: #2002557
        if (atmid[0] == 'H' &&
            (atmid[1] == 'D' || atmid[1] == 'E' ||
             atmid[1] == 'G' || atmid[1] == 'H' ||
             atmid[1] == 'N')) // HD, HE, HG, HH, HN...
          type = "H";

        if (type.size() == 2)
          type[1] = tolower(type[1]);

      } else { //must be hetatm record
        if (isalpha(element[1]) && (isalpha(element[0]) || (element[0] == ' '))) {
          if (isalpha(element[0]))
            type = element.substr(0,2);
          else
            type = element.substr(1,1);

          if (type.size() == 2)
            type[1] = tolower(type[1]);
        } else { // no element column to use
          if (isalpha(atmid[0])) {
            if (atmid.size() > 2)
              type = atmid.substr(0,2);
            else if (atmid[0] == 'A') // alpha prefix
              type = atmid.substr(1, atmid.size() - 1);
            else
              type = atmid.substr(0,1);
          } else if (atmid[0] == ' ')
            type = atmid.substr(1,1); // one char element
          else
            type = atmid.substr(1,2);

          // Some cleanup steps
          if (atmid == resname) {
            type = atmid;
            if (type.size() == 2)
              type[1] = tolower(type[1]);
          } else
            if (resname == "ADR" || resname == "COA" || resname == "FAD" ||
                resname == "GPG" || resname == "NAD" || resname == "NAL" ||
                resname == "NDP" || resname == "ABA") {
              if (type.size() > 1)
                type = type.substr(0,1);
              //type.erase(1,type.size()-1);
            } else // other residues
              if (isdigit(type[0])){
                type = type.substr(1,1);
              }
              else
                if (type.size() > 1 && isdigit(type[1]))
                  type = type.substr(0,1);
                else
                  if (type.size() > 1 && isalpha(type[1])) {
                    if (type[0] == 'O' && type[1] == 'H')
                      type = type.substr(0,1); // no "Oh" element (e.g. 1MBN)
                    else if(isupper(type[1])) {
                      type[1] = tolower(type[1]);
                    }
                  }
        }

      } // HETATM records
    } // no element column to use

    OBAtom atom;
    /* X, Y, Z */
    string xstr = sbuf.substr(24,8);
    string ystr = sbuf.substr(32,8);
    string zstr = sbuf.substr(40,8);
    vector3 v(atof(xstr.c_str()),atof(ystr.c_str()),atof(zstr.c_str()));
    atom.SetVector(v);

    double occupancy = atof(sbuf.substr(48, 6).c_str());
    OBPairFloatingPoint* occup = new OBPairFloatingPoint;
    occup->SetAttribute("_atom_site_occupancy");
    if (occupancy <= 0.0 || occupancy > 1.0){
      occupancy = 1.0;
    }
    occup->SetValue(occupancy);
    occup->SetOrigin(fileformatInput);
    atom.SetData(occup);

    // useful for debugging unknown atom types (e.g., PR#1577238)
    //    cout << mol.NumAtoms() + 1  << " : '" << element << "'" << " " << OBElements::GetAtomicNum(element.c_str()) << endl;
    if (elementFound)
      atom.SetAtomicNum(OBElements::GetAtomicNum(element.c_str()));
    else { // use our old-style guess from athe atom type
      unsigned int atomic_num = OBElements::GetAtomicNum(type.c_str());
      if (atomic_num ==  0) { //try one character if two character element not found
        type = type.substr(0,1);
        atomic_num = OBElements::GetAtomicNum(type.c_str());
      }
      atom.SetAtomicNum(atomic_num);
    }

    if ( (! scharge.empty()) && "  " != scharge )
      {
        if ( isdigit(scharge[0]) && ('+' == scharge[1] || '-' == scharge[1]) )
          {
            const char reorderCharge[3] = { scharge[1], scharge[0], '\0' };
            const int charge = atoi(reorderCharge);
            atom.SetFormalCharge(charge);
          }
        else
          {
            stringstream errorMsg;
            errorMsg << "WARNING: Problems reading a PDB file\n"
                     << "  Problems reading a HETATM or ATOM record.\n"
                     << "  According to the PDB specification,\n"
                     << "  columns 79-80 should contain charge of the atom\n"
                     << "  but OpenBabel found '" << scharge << "' (atom " << mol.NumAtoms()+1 << ").";
            obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obWarning);
          }
      }
    else {
      atom.SetFormalCharge(0);
    }

    /* residue sequence number */
    string resnum = sbuf.substr(16,4);
    OBResidue *res  = (mol.NumResidues() > 0) ? mol.GetResidue(mol.NumResidues()-1) : NULL;
    if (res == NULL
        || res->GetName() != resname
        || res->GetNumString() != resnum
        || res->GetChain() != chain
        || res->GetInsertionCode() != insertioncode)
      {
        vector<OBResidue*>::iterator ri;
        for (res = mol.BeginResidue(ri) ; res ; res = mol.NextResidue(ri))
          if (res->GetName() == resname
              && res->GetNumString() == resnum
              && static_cast<int>(res->GetChain()) == chain
              && static_cast<int>(res->GetInsertionCode()) == insertioncode) {
            if (insertioncode) fprintf(stderr,"I: identified residue wrt insertion code: '%c'\n",insertioncode);
            break;
          }

        if (res == NULL) {
          res = mol.NewResidue();
          res->SetChain(chain);
          res->SetName(resname);
          res->SetNum(resnum);
          res->SetInsertionCode(insertioncode);
        }
      }

    if (!mol.AddAtom(atom))
      return(false);
    else {
      OBAtom *atom = mol.GetAtom(mol.NumAtoms());

      res->AddAtom(atom);
      res->SetSerialNum(atom, atoi(serno.c_str()));
      res->SetAtomID(atom, sbuf.substr(6,4));
      res->SetHetAtom(atom, hetatm);

      return(true);
    }
  } // end reading atom records

} //namespace OpenBabel
