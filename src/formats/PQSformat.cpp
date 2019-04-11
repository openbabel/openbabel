/**********************************************************************
Copyright (C) 2003 by Pawel Wolinski
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
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/elements.h>
#include <cstdlib>

#include <ctype.h>

#if HAVE_STRINGS_H
#include <strings.h>
#endif

using namespace std;
namespace OpenBabel
{

  class PQSFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    PQSFormat()
    {
      OBConversion::RegisterFormat("pqs",this);
    }

    virtual const char* Description() //required
    {
      return
        "Parallel Quantum Solutions format\n"
        "No comments yet\n";
    };

    virtual const char* SpecificationURL()
    {return "http://www.pqs-chem.com/";};

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
      return READONEONLY | WRITEONEONLY;
    };

    //*** This section identical for most OBMol conversions ***
    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

  };
  //***

  //Make an instance of the format class
  PQSFormat thePQSFormat;

  /////////////////////////////////////////////////////////////////
  /* Lower first 5 characters of each word
   * (words separated by ' ' or '=')
   * Omit the filename after 'file=' card */
  void lowerit(char *s)
  {
    char tmp[6];
    unsigned int i, do_lower=5;
    for (i=0; i<strlen(s); i++)
      {
        if (s[i]==' ')
          do_lower=5;
        if (s[i]=='=')
          {
            strncpy(tmp,&s[i-4],5);
            tmp[5]='\0';
            if (strcmp(tmp,"file=")!=0)
              do_lower=5;
          }
        else
          {
            if (do_lower)
              {
                s[i]=tolower(s[i]);
                do_lower--;
              }
          }
      }
  }

  bool card_found(char *s)
  {
    int i;
    const char *input_cards[] ={"titl","file","cpu ","text","calc", // 0-4
                          "geom","basi","inte","gues","scf ", // 5-9
                          "forc","intc","freq","nbo ","pop ", //10-14
                          "pop=","semi","opti","mass","nmr ", //15-19
                          "lmp2","numh","rest","nucl","mp2 ", //20-24
                          "mem=","%mem","jump","clea","stop", //25-29
                          "mtst","dyna","anfc","corr","ffld", //30-34
                          "hess","path","scan","chk=","save", //35-39
                          "scr=","thre","iter","diis","lvsh", //40-44
                          "pseu","sthr","nodd","virt","fact", //45-49
                          "gran","anne","prin","loca","dft=", //50-54
                          "cuto","preo"};         //55-56
    lowerit(s);
    for (i=0; i<56; i++)
      if (strstr(s,input_cards[i])!=NULL)
        return true;
    return false;
  }


  int ReadPQS_geom(istream &ifs, OBMol &mol, const char *title,
                   int input_style, double bohr_to_angstrom)
  {
    int atom_count=0;
    double x, y, z;
    char buffer[BUFF_SIZE];
    string str;
    OBAtom *atom;
    vector<string> vs;

    mol.Clear();
    mol.BeginModify();

    while (ifs.getline(buffer,BUFF_SIZE) && !card_found(buffer))
      {
        if (buffer[0]!='$')
          {
            tokenize(vs, buffer);
            if (vs.size() < 1) return false; // timvdm 18/06/2008
            atom=mol.NewAtom();
            str=vs[0];
            if (input_style==0)
              {
                if (vs.size() < 4) return false; // timvdm 18/06/2008
                atom->SetAtomicNum(OBElements::GetAtomicNum(str.c_str()));
                x=atof((char*) vs[1].c_str())*bohr_to_angstrom;
                y=atof((char*) vs[2].c_str())*bohr_to_angstrom;
                z=atof((char*) vs[3].c_str())*bohr_to_angstrom;
              }
            else
              {
                if (vs.size() < 5) return false; // timvdm 18/06/2008
                str.replace (0,2,"");
                atom->SetAtomicNum(OBElements::GetAtomicNum(str.c_str()));
                x=atof((char*) vs[2].c_str())*bohr_to_angstrom;
                y=atof((char*) vs[3].c_str())*bohr_to_angstrom;
                z=atof((char*) vs[4].c_str())*bohr_to_angstrom;
              }
            atom->SetVector(x, y, z);
            atom_count++;
          }
      }

    mol.ConnectTheDots();
    mol.PerceiveBondOrders();

    mol.EndModify();
    mol.SetTitle(title);

    return atom_count;
  }


  /////////////////////////////////////////////////////////////////
  bool PQSFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title = pConv->GetTitle();

    char buffer[BUFF_SIZE];
    char coord_file[256];
    char full_coord_path[256]="\0";
    ifstream coordFileStream;
    double bohr_to_angstrom=1.0;
    unsigned int input_style, atom_count=0; //CM i removed
    bool geom_found;

    geom_found=false;
    while (!geom_found && ifs.getline(buffer,BUFF_SIZE))
      {
        lowerit(buffer);      //look for geom except in title or text
        if (strstr(buffer,"geom")!=NULL &&
            (strncmp(buffer,"text",4)!=0 && strncmp(buffer,"titl",4)!=0))
          {
            geom_found=true;
            lowerit(buffer);

            if (strstr(buffer,"bohr")!=NULL)
              bohr_to_angstrom=0.529177249;
            else
              bohr_to_angstrom=1.0;
            input_style=0;
            if (strstr(buffer,"=tx90")!=NULL)
              input_style=1;
            if (strstr(buffer,"=tx92")!=NULL)
              input_style=0;
            if (strstr(buffer,"=pqs" )!=NULL)
              input_style=0;

            if (strstr(buffer,"file=")!=NULL)
              {  //external geometry file
                strncpy(coord_file,strstr(buffer,"file=")+5, sizeof(coord_file));
                coord_file[sizeof(coord_file) - 1] = '\0';
                if (strrchr(coord_file,' ')!=NULL)
                  *strrchr(coord_file,' ')='\0';
                if (coord_file[0]!='/')
                  {
                    strncpy(full_coord_path,title, sizeof(full_coord_path));
                    full_coord_path[sizeof(full_coord_path)-1] = '\0';
                    if (strrchr(full_coord_path,'/')!=NULL)
                      *(strrchr(full_coord_path,'/')+1)='\0';
                    else
                      full_coord_path[0] = '\0';
                  }
                strcat(full_coord_path,coord_file);
                full_coord_path[sizeof(full_coord_path) - 1] = '\0';
                stringstream errorMsg;
                errorMsg <<"External geometry file referenced: "<< \
                  full_coord_path<<endl;
                obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obInfo);

                coordFileStream.open(full_coord_path);
                if (!coordFileStream)
                  {
                    obErrorLog.ThrowError(__FUNCTION__, "Cannot read external geometry file!", obError);
                    return(false);
                    //                    exit (-1);
                  }
                else
                  {
                    ifs.seekg(0, ios::end); //move .inp file pointer to the end of file

                    //New framework mods
                    OBConversion coordconv(&coordFileStream);
                    OBFormat* pFormat;
                    if (strstr(buffer,"=car" )!=NULL)
                      pFormat =OBConversion::FindFormat("BIOSYM");
                    if (strstr(buffer,"=hin" )!=NULL)
                      pFormat = OBConversion::FindFormat("HIN");
                    if (strstr(buffer,"=pdb" )!=NULL)
                      pFormat = OBConversion::FindFormat("PDB");
                    if (strstr(buffer,"=mop" )!=NULL)
                      pFormat = OBConversion::FindFormat("MOPAC");
                    return pFormat->ReadMolecule(&mol,&coordconv);

                    /*         if (strstr(buffer,"=car" )!=NULL)
                               return ReadBiosymCAR(coordFileStream, mol, title);
                               if (strstr(buffer,"=hin" )!=NULL)
                               return ReadHIN(coordFileStream, mol, title);
                               if (strstr(buffer,"=pdb" )!=NULL)
                               return ReadPDB(coordFileStream, mol, title);
                               if (strstr(buffer,"=mop" )!=NULL)
                               return ReadMOPAC(coordFileStream, mol, title);
                    */
                    //end new framework mods

                    //probably pqs's own xyz format
                    atom_count=ReadPQS_geom(coordFileStream,mol,title,
                                            input_style,bohr_to_angstrom);
                  }
              }
          }
      }

    if (geom_found)
      {
        if (atom_count==0)        //read directly form .inp file
          atom_count=ReadPQS_geom(ifs,mol,title,input_style,bohr_to_angstrom);
        if (atom_count==0)
          {   //try .coord file
            strncpy(coord_file,title, sizeof(coord_file));
            coord_file[sizeof(coord_file) - 1] = '\0';
            if (strrchr(coord_file,'.')!=NULL)
              *strrchr(coord_file,'.')='\0';
            strcat(coord_file,".coord");
            coordFileStream.open(coord_file);
            if (!coordFileStream)
              {
                stringstream errorMsg;
                errorMsg <<"ReadPQS: cannot read external "<<coord_file<<" file!"<<endl;
                obErrorLog.ThrowError(__FUNCTION__, errorMsg.str(), obError);
                return(false);
                //                exit (-1);
              }
            else
              atom_count=ReadPQS_geom(coordFileStream,mol,title,0,
                                      bohr_to_angstrom);
          }
      }
    else
      obErrorLog.ThrowError(__FUNCTION__, "Error reading PQS file.  GEOM card not found!", obWarning);

    ifs.seekg(0, ios::end); //move .inp file pointer to the end of file
    if (atom_count>0)
      return true;
    else
      return false;
  }

  ////////////////////////////////////////////////////////////////

  bool PQSFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    unsigned int i;
    char buffer[BUFF_SIZE];
    OBAtom *atom;
    ofs<<"TEXT="<<mol.GetTitle()<<endl;
    ofs<<"GEOM=PQS"<<endl;
    for (i=1; i<=mol.NumAtoms(); i++)
      {
        atom=mol.GetAtom(i);
        snprintf(buffer, BUFF_SIZE, "%s           %10.6lf   %10.6lf   %10.6lf",
                OBElements::GetSymbol(atom->GetAtomicNum()),
                atom->GetX(),
                atom->GetY(),
                atom->GetZ());
        ofs<<buffer<<endl;
      }
    return(true);
  }

} //namespace OpenBabel
