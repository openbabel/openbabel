/**********************************************************
 *                                                        *
 *                    Y  A  S  A  R  A                    *
 *                                                        *
 * Yet Another Scientific Artificial Reality Application  *
 *                                                        *
 **********************************************************
 *         www.YASARA.org - Watching Nature@Work          *
 **********************************************************
 *   Reader and Writer for Yasara Objects (*.YOB files)   *
 **********************************************************
 *            (C) 2002-2006 by Elmar Krieger              *
 **********************************************************

This reader/writer for Yasara objects is free software;
you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

***********************************************************/
#include <openbabel/babelconfig.h>

/* Adapt on some 64bit machines! */
typedef int int32;
typedef unsigned int uint32;
typedef short int16;
typedef unsigned short uint16;

#include <openbabel/obmolecformat.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/elements.h>
#include <openbabel/bond.h>
#include <cstdlib>


#define mem_alloc malloc
#define mem_free free
#define mem_set memset
#define mem_inc(POINTER,INC) (*((char**)&(POINTER))+=(INC))
#define str_copy strcpy
#define STR_____ 0x20202020
#define FMTOA 0.00001
#define ATOFM 100000.

using namespace std;

namespace OpenBabel {

/* DECODE A SIGNED INT32 STORED IN LITTLE ENDIAN FORMAT
   ==================================================== */
int int32lemem(char *data)
{ int value,intmino2;

  value=(int)((unsigned char*)data)[0]+((int)((unsigned char*)data)[1]<<8)+
        ((int)((unsigned char*)data)[2]<<16)+((int)((unsigned char*)data)[3]<<24);
  intmino2=-1073741824;
  if ((unsigned int) INT_MIN!=0x80000000&&(value&0x80000000)) value=(intmino2*2)+(value&0x7fffffff);
  return(value); }

int int32le(int32 value)
{ return(int32lemem((char*)&value)); }

/* DECODE AN UNSIGNED INT32 STORED IN LITTLE ENDIAN FORMAT
   ======================================================= */
unsigned int uint32lemem(char *data)
{ return((int)((unsigned char*)data)[0]+((int)((unsigned char*)data)[1]<<8)+
         ((int)((unsigned char*)data)[2]<<16)+((int)((unsigned char*)data)[3]<<24)); }

unsigned int uint32le(uint32 value)
{ return(uint32lemem((char*)&value)); }

/* STORE AN INT32 IN LITTLE ENDIAN FORMAT
   ====================================== */
void storeint32le(char *data,int32 value)
{ int i;

  for (i=0;i<4;i++) data[i]=(value>>(i*8))&255; }

/* CONVERT MAXIMALLY n CHARACTERS TO INTEGER
   ========================================= */
int str_natoi(char *string,int number)
{ char ch;
  int i;

  for (i=0;i<number;i++) if (!string[i]) return(atoi(string));
  ch=string[number];
  string[number]=0;
  i=atoi(string);
  string[number]=ch;
  return(i); }

/* COPY MAXIMALLY n CHARACTERS AND TERMINATE WITH ZERO
   =================================================== */
void str_ncopy(char *string1,char *string2,int len)
{ int i;
  char ch;

  for (i=0;i<len;i++)
  { ch=string2[i];
    string1[i]=ch;
    if (!ch) break; }
  string1[i]=0; }

/* ELEMENT SYMBOLS */
#define MOB_ELEMENTS 128
const char *mob_elementsym[MOB_ELEMENTS]=
{"?",
 "H","He",
 "Li","Be","B","C","N","O","F","Ne",
 "Na","Mg","Al","Si","P","S","Cl","Ar",
 "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
 "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe",
 "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu",
 "Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W", "Re","Os","Ir","Pt","Au",
 "Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U", "Np","Pu","Am",
 "Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","?","?",
 "?","?","?","?","?","?","?","?","?","?","?","?","?","?","?","?" };

/* HEADER */
#define MOB_ATOMS 0
#define MOB_TERMINI 1
#define MOB_TERLIST 2
#define MOB_HEADERLEN 3
#define MOB_HEADERSIZE 12
/* ATOM ENTRIES */
#define MOB_LINKS 0
#define MOB_LINKSMASK 15
#define MOB_IDLEN 1
#define MOB_ELEMENT 2
#define MOB_ELEMENTMASK 127
#define MOB_HETATOMFLAG 128
#define MOB_COLOR 3
#define MOB_UNDEFCOLOR 64
#define MOB_POSX 4
#define MOB_POSY 8
#define MOB_POSZ 12
#define MOB_LINK1 16
#define MOB_LINKTYPEMASK 0xff000000
#define MOB_LINKTYPESHIFTS 24
#define MOB_LINKATOMMASK 0xffffff
#define MOB_PRINTRESIDSIZE 10
#define MOB_NOATOM -1
/* LINK TYPES */
#define MOB_LINKTYPES 10
#define MOB_LINKNONE 0
#define MOB_LINKSINGLE 1
#define MOB_LINKDOUBLE 2
#define MOB_LINKTRIPLE 3
#define MOB_LINKRESONANCE50 4
#define MOB_LINKRESONANCE25 5
#define MOB_LINKRESONANCE33 6
#define MOB_LINKRESONANCE66 7
#define MOB_LINKRESONANCE75 8
#define MOB_LINKQUADRUPLE 9
/* ATOM ID STRUCTURE */
#define MOB_RESNAMEMASK 0xffffff
#define MOB_CHAINMASK 0xff000000
#define MOB_CHAINSHIFTS 24
/* FLAGS OF ID FIELDS PRESENT IN MOB DATA */
#define MOB_ATOMFLAG 1
#define MOB_RESFLAG 2
#define MOB_INSCODEALTLOCFLAG 4
#define MOB_OCCUPANCYFLAG 8
#define MOB_TEMPFACFLAG 16
#define MOB_SEGMENTFLAG 32
#define MOB_PROPERTYFLAG 0x2000
#define MOB_AROMATICFLAG 0x4000
#define MOB_DOTFLAG 0x40000
#define MOB_SURFENVFLAG 0x80000
#define MOB_FIXEDFLAG 0x800000
#define MOB_ATOMSIZEMAX 44         // Header+3*POS+Flags+Atom+ResChain+Numb+AltIns+Occupany+Bfac
#define MOB_INFOHISTORY 6
#define MOB_INFOHISTORYSIZE 136
#define MOB_INFOEND 0x7fffffff
#define MOB_INFOENDSIZE 8
#define mob_element(ATOM) ((ATOM)->header[MOB_ELEMENT]&MOB_ELEMENTMASK)
#define mob_links(ATOM) ((ATOM)->header[MOB_LINKS]&MOB_LINKSMASK)

typedef int32 mobdata;

struct mobatom                     /* Modeling object atom structure */
{ unsigned char header[4];
  int32 posx;
  int32 posy;
  int32 posz;
  int32 link[4]; };

struct atomid                      /* The order of the fields is not important */
{ uint32 atom;
  uint32 resnamechain;
  uint32 resno;
  uint32 restype;
  uint16 inscode;
  uint16 altloc;
  uint32 flags;                    /* Use the same MOB_ masks as mobatom */
  uint32 secstr;
  uint32 segment;
  float occupancy;
  float bfactor;
  float property; };

/* SET POINTER TO DATA START
   ========================= */
struct mobatom *mob_start(mobdata *mob)
{ int termini;

  termini=uint32le(mob[MOB_TERMINI]);
  return((struct mobatom*)(mob+termini+2)); }

/* GET SIZE OF ATOM
   ================ */
int mob_atomlen(struct mobatom *atom)
{ return(4+mob_links(atom)+atom->header[MOB_IDLEN]); }

int mob_atomsize(struct mobatom *atom)
{ return(mob_atomlen(atom)*sizeof(mobdata)); }

/* MOVE POINTER TO NEXT ATOM OF YASARA OBJECT STRUCTURE
   ==================================================== */
struct mobatom *mob_next(struct mobatom *atom)
{ mem_inc(atom,mob_atomsize(atom));
  return(atom); }

void mob_setnext(struct mobatom **atomadd)
{ *atomadd=mob_next(*atomadd); }

/* GET ATOM ID
   =========== */
void mob_getid(struct atomid *id,struct mobatom *atom)
{ int idpos,links,flags,inscodealtloc;

  links=mob_links(atom);
  flags=int32le(atom->link[links]);
  idpos=links+1;
  id->atom=atom->link[idpos++];
  id->resnamechain=atom->link[idpos++];
  id->resno=atom->link[idpos++];
  if (flags&MOB_INSCODEALTLOCFLAG)
  { inscodealtloc=int32le(atom->link[idpos++]);
    id->inscode=inscodealtloc&0xffff;
    id->altloc=inscodealtloc&0xffff; }
  else id->inscode=id->altloc=0;
  if (flags&MOB_OCCUPANCYFLAG) id->occupancy=((float*)atom->link)[idpos++];
  else id->occupancy=1.0f;
  if (flags&MOB_TEMPFACFLAG) id->bfactor=((float*)atom->link)[idpos++];
  else id->bfactor=0;
  if (flags&MOB_SEGMENTFLAG) id->segment=atom->link[idpos++];
  else id->segment=0;
  if (flags&MOB_PROPERTYFLAG) id->property=((float*)atom->link)[idpos++];
  else id->property=0;
  id->flags=flags&(MOB_SURFENVFLAG|MOB_DOTFLAG); }

/* CLEAR ATOMID
   ============ */
#define SEC_UNDEFINED 4
void mob_clearid(struct atomid *id)
{ id->atom=id->resnamechain=id->resno=STR_____;
  id->inscode=id->altloc=0;
  id->occupancy=1.0;
  id->bfactor=0.0;
  id->secstr=SEC_UNDEFINED; }

/* INVALIDATE ATOMID
   ================= */
void mob_invid(struct atomid *id)
{ id->atom=id->resnamechain=id->resno=0xffffffff;
  id->inscode=id->altloc=0;
  id->secstr=SEC_UNDEFINED; }

/* COMPARE TWO RESIDUE NAMES
   ========================= */
int mob_issameres(struct atomid *id1,struct atomid *id2)
{ if (id1->resnamechain==id2->resnamechain&&id1->resno==id2->resno&&
      id1->inscode==id2->inscode) return(1);
  return(0); }

/* CHECK IF ATOM BELONGS TO GIVEN RESIDUE
   ====================================== */
int mob_hasres(struct mobatom *atom,struct atomid *id2)
{ struct atomid id1;

  mob_getid(&id1,atom);
  return(mob_issameres(&id1,id2)); }

/* GET LENGTH OF RESIDUE
   =====================
   atomsleft SPECIFIES THE NUMBER OF ATOMS LEFT
   TILL END OF STRUCTURE IS REACHED */
int mob_reslen(struct mobatom *atom,int atomsleft)
{ int i;
  struct atomid id;

  mob_getid(&id,atom);
  i=0;
  while (atomsleft--)
  { if (!mob_hasres(atom,&id)) break;
    i++;
    atom=mob_next(atom); }
  return(i); }

/* CHECK IF ATOM HAS THE GIVEN NAME
   ================================ */
int mob_hasname(struct mobatom *atom,struct atomid *id2)
{ struct atomid id1;

  mob_getid(&id1,atom);
  if (id1.atom==id2->atom&&id1.altloc==id2->altloc) return(1); else return(0); }


class YOBFormat : public OBMoleculeFormat
{
public:
    //Register this format type ID
    YOBFormat()
    {
        OBConversion::RegisterFormat("yob",this);
    }

    virtual const char* Description() //required
    {
        return
            "YASARA.org YOB format\n"
            "The native YASARA format.\n";
    };

    virtual const char* SpecificationURL(){return
            "http://www.yasara.org";}; //optional

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
        return READONEONLY|READBINARY|WRITEBINARY;
    };

    //*** This section identical for most OBMol conversions ***
    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);


};
//***

//Make an instance of the format class
YOBFormat theYOBFormat;

/* READ YASARA OBJECT
   ==================
   IF options CONTAINS 'f' (COMMAND LINE -af), ATOM NAMES WILL STAY FIXED (I.E. LEAVING THE SPACE IN THE FIRST COLUMN) */
bool YOBFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{
  OBMol* pmol = pOb->CastAndClear<OBMol>();
  if(pmol==NULL)
      return false;

  //Define some references so we can use the old parameter names
  istream &ifs = *pConv->GetInStream();
  OBMol &mol = *pmol;
  //   const char* title = pConv->GetTitle();

  bool hetatom;
  char buffer[8],resname[4],atomname[5];
  string str;
  unsigned int i,j/*,m,q*/;
  unsigned int /*resno,chainNum,*/link,linked,linktype,atoms,element,links,chain;
  int /*samenames,*/reslen,charged;
  unsigned int infosize,size;
  mobdata *mob;
  struct mobatom *srcatom,/**atom2,*/*resstart;
  struct atomid id;
  OBAtom *dstatom;
  OBResidue *res;
  bool has_residue_info = false;

  /* VERIFY FILE FORMAT */
  ifs.read(buffer,8);
  if (!EQn(buffer,"YMOB",4)) return(false);
  /* IGNORE ENTIRE INFO STRUCTURE */
  infosize=uint32lemem(&buffer[4]);
  for (i=0;i<infosize;i++) ifs.read(buffer,1);
  /* READ ATOM DATA */
  ifs.read(buffer,4);
  size=uint32lemem(buffer);
  mob=(mobdata*)mem_alloc(size);
  if (!mob) return(false);
  ifs.read((char*)mob,size);
  /* PARSE THE ATOMS */
  mol.Clear();
  mol.BeginModify();
  mob_invid(&id);
  atoms=uint32le(mob[MOB_ATOMS]);
  srcatom=mob_start(mob);
  res=NULL;
  charged=0;
  for (i=0;i<atoms;i++)
  { /* GET ELEMENT, TYPE AND POSITION */
    element=mob_element(srcatom);
    dstatom=mol.NewAtom();
    dstatom->SetAtomicNum(element);
    dstatom->SetType(mob_elementsym[element]);
    vector3 pos(int32le(srcatom->posx)*-FMTOA,int32le(srcatom->posy)*FMTOA,
                int32le(srcatom->posz)*FMTOA);
    dstatom->SetVector(pos);
    if (!mob_hasres(srcatom,&id))
    { /* NEW RESIDUE FOUND */
      has_residue_info = true;
      resstart=srcatom;
      reslen=mob_reslen(resstart,atoms-i);
      mob_getid(&id,srcatom);
      res=mol.NewResidue();
      *((int32*)resname)=id.resnamechain;
      chain=resname[3];
      if (chain>='0'&&chain<='9') chain-=48;
      else if (chain>='A'&&chain<='Z') chain-=64;
      else if (chain>='a'&&chain<='z') chain-=96;
      resname[3]=0;
      res->SetChainNum(chain);
      str=resname;
      res->SetName(str);
      res->SetNum(str_natoi((char*)&id.resno,4)); }
    else mob_getid(&id,srcatom);
    /* SET THE CHARGE, WHICH IS PASSED IN THE property FIELD */
    dstatom->SetPartialCharge(id.property);
    if (id.property!=0) charged=1;
    res->AddAtom(dstatom);
    /* WHO KNOWS WHY HALF OF THE ATOM PARAMETERS CAN BE SET DIRECTLY
       ABOVE, AND THE SECOND HALF NEEDS A RESIDUE OBJECT IN BETWEEN? */
    res->SetSerialNum(dstatom,i+1);
    *((int32*)atomname)=id.atom;
    atomname[4]=0;
    /* SHIFT ATOM NAME BY ONE CHARACTER TO REMOVE LEADING SPACE */
    if (atomname[0]==' '&&(!pConv->IsOption("f",OBConversion::INOPTIONS))) memcpy(atomname,atomname+1,4);
    /* RENAME TERMINAL OXYGENS */
    str=atomname;
    if (str=="OT1") str="O";
    if (str=="OT2") str="OXT";
    res->SetAtomID(dstatom,str);
    if (srcatom->header[MOB_ELEMENT]&MOB_HETATOMFLAG) hetatom=true;
    else hetatom=false;
    res->SetHetAtom(dstatom,hetatom);
    /* NOW ADD THE BONDS */
    links=srcatom->header[MOB_LINKS];
    for (j=0;j<links;j++)
    { link=uint32le(srcatom->link[j]);
      linked=link&MOB_LINKATOMMASK;
      if (linked<i)
      { linktype=link>>MOB_LINKTYPESHIFTS;
        if (linktype==MOB_LINKQUADRUPLE) linktype=4;
        else if (linktype>MOB_LINKTRIPLE) linktype=5;
        mol.AddBond(i+1,linked+1,linktype); } }
    mob_setnext(&srcatom); }
  mem_free(mob);

  mol.EndModify();
  if (charged) mol.SetPartialChargesPerceived();
  if (has_residue_info)
    mol.SetChainsPerceived();
  if (!mol.NumAtoms()) return(false);
  return(true); }

/* WRITE YASARA OBJECT
   ===================
   IF options CONTAINS 'f' (COMMAND LINE -xf), ATOM NAMES WILL STAY FIXED (I.E. ADDING NO SPACE IN THE FIRST COLUMN)
   IF options CONTAINS 'a' (COMMAND LINE -xa), AROMATIC ATOMS WILL BE FLAGGED */
bool YOBFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
  OBMol* pmol = dynamic_cast<OBMol*>(pOb);
  if(pmol==NULL)
      return false;

  //Define some references so we can use the old parameter names
  ostream &ofs = *pConv->GetOutStream();
  OBMol &mol = *pmol;

  //  bool hetatom;
  char buffer[32],/*resname[4],*/atomname[5];
  unsigned char double1[8]={0,0,0,0,0,0,0xf0,0x3f};
  //   char *str;
  int i,j,/*m,q,*/pos;
  int /*resno,chainNum,link,linktype,*/atoms,element,links/*,chain*/;
  int /*samenames,reslen,*/bondorder,flags;
  unsigned int /*infosize,*/size;
  //   mobdata *mob;
  //   struct mobatom *atom2,*resstart;
  //   struct atomid id;
  OBAtom *srcatom,*linkedatom;
  OBResidue *res;
  OBBond *bond;
  vector<OBAtom*>::iterator e;
  vector<OBBond*>::iterator iter;

  atoms=mol.NumAtoms();
  if (!atoms) return(false);
  /* WRITE FILE FORMAT */
  ofs<<"YMOB";
  /* WRITE AN INFO STRUCTURE WITH EMPTY TRANSFORMATION HISTORY (UNIT MATRIX) */
  storeint32le(buffer,MOB_INFOHISTORYSIZE+MOB_INFOENDSIZE);
  ofs.write(buffer,4);
  storeint32le(buffer,MOB_INFOHISTORY);
  storeint32le(&buffer[4],MOB_INFOHISTORYSIZE);
  ofs.write(buffer,8);
  mem_set(buffer,0,8);
  for (i=0;i<4;i++)
  { for (j=0;j<4;j++)
    { if (i==j) ofs.write((char*)double1,8);
      else ofs.write(buffer,8); } }
  storeint32le(buffer,MOB_INFOEND);
  storeint32le(&buffer[4],MOB_INFOENDSIZE);
  ofs.write(buffer,8);
  /* DETERMINE SIZE OF ATOM DATA */
  size=MOB_HEADERSIZE;
  for (i=1;i<=atoms;i++)
  { srcatom=mol.GetAtom(i);
    /* COUNT BONDS */
    links=0;
    for (bond=srcatom->BeginBond(iter);bond;bond=srcatom->NextBond(iter)) links++;
    size+=links*4+32; }
  /* WRITE DATA HEADER */
  storeint32le(buffer,size);
  storeint32le(&buffer[4],atoms);
  storeint32le(&buffer[8],1);
  storeint32le(&buffer[12],atoms-1);
  ofs.write(buffer,16);
  /* WRITE ATOMS */
  for (i=1;i<=atoms;i++)
  { srcatom=mol.GetAtom(i);
    /* COUNT BONDS AGAIN (EASY ACCESS TO THE NUMBER OF BONDS IS TOO WELL HIDDEN ;-) */
    links=0;
    for (bond=srcatom->BeginBond(iter);bond;bond=srcatom->NextBond(iter)) links++;
    /* WRITE ATOM HEADER+POSITION */
    buffer[MOB_LINKS]=links;
    buffer[MOB_IDLEN]=4;
    element=srcatom->GetAtomicNum();
    buffer[MOB_ELEMENT]=element;
    buffer[MOB_COLOR]=MOB_UNDEFCOLOR;
    storeint32le(&buffer[4],(int)(srcatom->GetX()*-ATOFM));
    storeint32le(&buffer[8],(int)(srcatom->GetY()*ATOFM));
    storeint32le(&buffer[12],(int)(srcatom->GetZ()*ATOFM));
    ofs.write(buffer,16);
    /* WRITE BONDS */
    //printf("Babel atom %d with %d links:\n",i,links);
    for (linkedatom=srcatom->BeginNbrAtom(iter);linkedatom;linkedatom=srcatom->NextNbrAtom(iter))
    { storeint32le(buffer,linkedatom->GetIdx()-1);
      bondorder=(*iter)->GetBondOrder();
      //printf("  Order %d\n",bondorder);
      if (bondorder==4) bondorder=MOB_LINKQUADRUPLE;
      else if (bondorder==5) bondorder=MOB_LINKRESONANCE50;
      buffer[3]=bondorder;
      ofs.write(buffer,4); }
    /* WRITE ATOM ID */
    mem_set(buffer,0,sizeof(buffer));
    flags=MOB_ATOMFLAG|MOB_RESFLAG;
    if (pConv->IsOption("a",OBConversion::OUTOPTIONS)&&srcatom->IsAromatic()) flags|=MOB_AROMATICFLAG;
    storeint32le(buffer,flags);
    if (srcatom->HasResidue())
    { res=srcatom->GetResidue();
      str_ncopy(atomname,(char*)res->GetAtomID(srcatom).c_str(),4);
      pos=4;
      if (!pConv->IsOption("f",OBConversion::OUTOPTIONS))
      { /* GUESS IF THE ATOM NAME NEEDS TO BE SHIFTED */
        if (strlen(mob_elementsym[element])==1||strncasecmp(mob_elementsym[element],atomname,2)) pos=5; }
      str_copy(&buffer[pos],atomname);
      str_copy(&buffer[8],(char*)res->GetName().c_str());
      snprintf(&buffer[12], 4, "%4d",res->GetNum()); }
    else
    { str_copy(&buffer[4],OBElements::GetSymbol(element));
      str_copy(&buffer[8],"UNK    1"); }
    for (j=4;j<16;j++) if (!buffer[j]) buffer[j]=' ';
    ofs.write(buffer,16); }
  return(true); }

}
