/**********************************************************************
Copyright (C) 2007 by Mike N. Burnett, burnettmn@ornl.gov
Copyright (C) 2007 by Sergei V. Trepalin sergey_trepalin@chemical-block.com
Copyright (C) 2007 by Andrei Gakh andrei.gakh@nnsa.doe.gov

Code translation from Java

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
#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/elements.h>
#include <openbabel/bond.h>

#include <openbabel/mcdlutil.h>
#include <cstdlib>

using namespace std;
namespace OpenBabel
{

class MCDLFormat : public OBMoleculeFormat
{
public:
  MCDLFormat()
  {
    OBConversion::RegisterFormat("mcdl",this);
    init();
  }

  virtual const char* Description() //required
  {
    return
    "MCDL format\n"
    "Modular Chemical Descriptor Language\n\n"
    "As described in [gb2001]_.\n\n"

".. [gb2001] A.A. Gakh and M.N. Burnett. **Modular Chemical Descriptor\n"
"            Language (MCDL): Composition, Connectivity and\n"
"            Supplementary Modules.**\n"
"            *J. Chem. Inf. Comput. Sci.*, **2004**, *41*, 1491-1499.\n"
"            [`Link <https://doi.org/10.1021/ci000108y>`_]\n\n"

"Here's an example conversion from SMILES to MCDL::\n\n"
"  obabel -:\"CC(=O)Cl\" -omcdl\n"
"  CHHH;COCl[2]\n";
  }

  virtual const char* SpecificationURL(){return
     "http://pubs.acs.org/cgi-bin/abstract.cgi/jcisd8/2001/41/i06/abs/ci000108y.html";}

  virtual const char* GetMIMEType()
  { return "chemical/x-MCDL"; }

  /* Flags() can return be any of the following combined by |
     or be omitted if none apply
     NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY  DEFAULTFORMAT
     READBINARY  WRITEBINARY  READXML  ZEROATOMSOK*/
  virtual unsigned int Flags()
  {
      return 0;
  }

  virtual int SkipObjects(int n, OBConversion* pConv)
  {
      if(n==0) n++;
      string temp;
      istream& ifs = *pConv->GetInStream();
      do {
          if(ifs.good())
            getline(ifs, temp);
        }while(ifs.good() && --n);
      return ifs.good() ? 1 : -1;
  }

  ////////////////////////////////////////////////////
  /// Declarations for the "API" interface functions. Definitions are below
  virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
  virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

private:

  string fsastart;
  string fsbstart;
  string fchstart;
  string fradstart;
    //coordinates string
  string fnastart;
  string fnbstart;
  string fzcoorstart;
  string fablockstart;
  string fbblockstart;
  string fchargeblockstart;
  string fstereobondstart;
  string ftitlestart;

//  int fragNo;
  int maxdepth;
  int kflag;

  int ntatoms;
  int nbonds;
  string  finalstr;

  int qx [MAXFRAGS];
  int qa [MAXBONDS][4];
//  int HVal[NELEMMAX];

private:
  void init();
  void initGlobals();
  void solve(int ntypes, int z[MAXBONDS][4], int depth);
  string constring(int conntab [MAXBONDS][4], char * tstr);
  string intToStr(int k);
  string getMCDL(OBMol* pmol, bool includeCoordinates);
  void restoreFullMCDL(string value, OBMol* pmol);
  void setMCDL(const string lineToParse, OBMol* pmol, string & sout);
  void assignCharges(const std::vector <int> aPosition, const std::vector <int> iA1,
    const std::vector <int> iA2, std::vector <int>& aCharges, std::vector <int>& charges,
    std::vector <int>& bondOrder, int aPos, int nPrev, int nt, int acount, int bcount);
  int indexOf(const string instring, const string substring, int fromPos=0);
  int lastIndexOf(const string instring, const string substring);
  bool parseFormula(const string formulaString, std::vector <int>& enumber);
  string getMolTitle(string & line);

};



  void MCDLFormat::init(){
    fsastart="{SA:";
    fsbstart="{SB:";
    fchstart="{CZ:";
    fradstart="{RA:";
    //coordinates string
    fnastart="{NA:";
    fnbstart="{NB:";
    fzcoorstart="{ZV:";
    fablockstart="{CC:";
    fbblockstart="{BB:";
    fchargeblockstart="{MM:CHG,";
    fstereobondstart="{BS:";
    ftitlestart="{CN:}";

   };

  void MCDLFormat::initGlobals(){
    int i,j;

    maxdepth=0;
    kflag=0;

    ntatoms=0;
    nbonds=0;
    finalstr="";

  for (i=0; i<MAXFRAGS; i++) {
    qx[i]=0;
    for (j=0; j<4; j++) qa[i][j]=0;
  };

  };



/* Create descriptor connectivity string from connectivity table of fragments */
  string MCDLFormat::constring(int conntab [MAXBONDS][4], char * tstr)
{
    int  i,j,k,n,nn,icons[6],comma;
    char line[82],semis[100];
    string result="";

    result="[";
    semis[0] = '\0';

    for (i=0; i<ntatoms; i++)
    {
        if (i > 0)
            strcat(semis,";");

        n = 0;
        for (j=0; j<nbonds; j++)
            if (conntab[j][2] == (i+1))
                icons[n++] = conntab[j][3];
        if (n > 1)
        {
            for (j=0; j<n-1; j++)
                for (k=j+1; k<n; k++)
                    if (icons[k] < icons[j])
                    {
                        nn = icons[j];
                        icons[j] = icons[k];
                        icons[k] = nn;
                    }

        }
        comma = 0;
        for (j=0; j<n; j++)
        {
            if (icons[j] > (i+1) && comma == 0)
            {
                sprintf(line,"%s%d",semis,icons[j]);
                result=result+line;//strcat(sp,line);
                comma = 1;
                semis[0] = '\0';
            }
            else if (icons[j] > (i+1) && comma == 1)
            {
                sprintf(line,",%d",icons[j]);
                result=result+line;//strcat(sp,line);
            }

        }
    }

    result=result+"]";
    return result;
}


  void MCDLFormat::solve(int ntypes, int z[MAXBONDS][4], int depth)
{
    int  i,dupnum,j,k,ktype,nn,nt;
    bool iflag,newone;
    int *  nsum[MAXBONDS];
    bool lflag[MAXFRAGS];
    char strg[MAXFRAGS+1];
	char * strngs[MAXFRAGS+1];
    char tstr[MAXFRAGS+1];
    int  numdups, dupfrag, jump;
    bool jflag;
    int  ix[MAXFRAGS],conntab[MAXBONDS][4],cx[MAXFRAGS];
    int  mx[MAXFRAGS];

	//stack overflow message-move data from stack to heap
	for (i=0; i<=MAXFRAGS; i++) strngs[i]=(char *)malloc(MAXFRAGS);
    for (i=0; i<MAXBONDS; i++)  nsum[i]=(int *) malloc(MAXFRAGS);

    // depth = recursion level
    if (depth > 10)
    {
        printf("Ten recursion levels exceeded.\n");
        return;
        //exitprog(freeflag);
    }
    if (depth > maxdepth)
        maxdepth = depth;

    for (i=0; i<nbonds; i++)
    {
        for (j=0; j<4; j++)
            conntab[i][j] = z[i][j];
        mx[conntab[i][0]-1] = conntab[i][2];
    }

    jflag = true;
    if (depth == 1)
    {
        jflag = false;
        for (i=0; i<ntatoms; i++)
        {
            ix[i] = mx[i];
            if (ix[i] != mx[0])
                jflag = true;
        }
    }

    iflag = true;
    while (iflag)
    {
        if (jflag)
            for (i=0; i<ntatoms; i++)
                ix[i] = 0;

            nt = 1;
            ktype = 1;
            iflag = false;

           while (nt <= ntypes && jflag)
           {
            for (j=0; j<ntatoms; j++)
            {
                   lflag[j] = false;
                for (i=0; i<nbonds; i++)
                       nsum[i][j] = 0;
            }

            for (i=0; i<nbonds; i++)
            if (conntab[i][2] == nt)
            {
                j = conntab[i][0]-1;
                k = conntab[i][3]-1;
                nsum[j][k]++;
                lflag[j] = true;
            }

                nn = 0;

            for (i=0; i<ntatoms; i++)
            if (lflag[i])
            {
                for (j=0; j<ntypes; j++)
                    strg[j] = nsum[i][j]+48;
                strg[ntypes] = '\0';
                newone = true;
                for (j=0; j<nn; j++)
                    if (strcmp(strngs[j],strg) == 0)
                    {    newone = false;
                        break;
                    }
                if (newone)
                    strcpy(strngs[nn++],strg);
            }

            if (nn > 1)
            {
            iflag = true;

            // sort strings
            for (i=0; i<nn-1; i++)
                for (j=i+1; j<nn; j++)
                    if (strcmp(strngs[j],strngs[i]) > 0)
                    {
                    strcpy(tstr,strngs[i]);
                    strcpy(strngs[i],strngs[j]);
                    strcpy(strngs[j],tstr);
                    }
                }

            for (i=0; i<ntatoms; i++)
            if (lflag[i])
            {
                for (j=0; j<ntypes; j++)
                    strg[j] = nsum[i][j]+48;
                strg[ntypes] = '\0';
                for (j=0; j<nn; j++)
                    if (strcmp(strg,strngs[j]) == 0)
                        ix[i] = ktype+j;
            }

                nt++;
                ktype += nn;
           }

        ntypes = 0;
        for (i=0; i<nbonds; i++)
        {
            conntab[i][2] = ix[conntab[i][0]-1];
            conntab[i][3] = ix[conntab[i][1]-1];
            if (conntab[i][2] > ntypes)
                ntypes = conntab[i][2];
        }

            if (iflag)
            continue;

        if (ntypes < ntatoms)
        {
            for (j=0; j<ntypes; j++)
            {
            numdups = 0;
            for (i=0; i<ntatoms; i++)
                if (ix[i] == j+1)
                {
                    numdups++;
                    jump = ix[i] - mx[i];
                }
            if (numdups > 1)
            {
                dupfrag = j+1;
                break;
            }
            }

            for (i=0; i<ntatoms; i++)
                cx[i] = ix[i];

            for (j=0; j<numdups; j++)
            {
                for (i=0; i<ntatoms; i++)
                    ix[i] = cx[i];

            dupnum = 0;
                for (i=0; i<ntatoms; i++)
            {
                if (ix[i] > dupfrag)
                    ix[i] = mx[i] + jump + 1;
                else if (ix[i] == dupfrag && j != dupnum)
                {
                    dupnum++;
                    ix[i] = mx[i] + jump + 1;
                }
                else if (ix[i] == dupfrag && j == dupnum)
                    dupnum++;
            }

            ntypes = 0;
                for (i=0; i<nbonds; i++)
                {
                conntab[i][2] = ix[conntab[i][0]-1];
                conntab[i][3] = ix[conntab[i][1]-1];
                if (conntab[i][2] > ntypes)
                    ntypes = conntab[i][2];
                }
                solve(ntypes,conntab,depth+1);
            }
        }
    }

    if (ntypes == ntatoms)
    {
            //if (true) fprintf(o_verbose,"%s\n",constring(conntab,tstr));

        if (!kflag)
        {
            kflag++;
      finalstr=constring(conntab,tstr);
            for (i=0; i<ntatoms; i++)
                qx[i] = ix[i];
            for (i=0; i<nbonds; i++)
                for (j=0; j<4; j++)
                    qa[i][j] = conntab[i][j];
        }
        else if (finalstr != constring(conntab,tstr))
        {
            for (i=0; i<ntatoms; i++)
            {
                tstr[0] = '\0';
                strg[0] = '\0';
                for (j=0; j<ntatoms; j++)
                {
                    strcat(tstr,"0");
                    strcat(strg,"0");
                }
                for (j=0; j<nbonds; j++)
                {
                    if (conntab[j][2] == (i+1))
                           tstr[conntab[j][3]-1] = '1';
                    if (qa[j][2] == (i+1))
                           strg[qa[j][3]-1] = '1';
                }
                if (strcmp(tstr,strg) < 0)
                       break;
                if (strcmp(tstr,strg) > 0)
                {
                    finalstr=constring(conntab,tstr);
                    for (j=0; j<ntatoms; j++)
                        qx[j] = ix[j];
                    for (j=0; j<nbonds; j++)
                        for (k=0; k<4; k++)
                            qa[j][k] = conntab[j][k];
                    break;
                }
            }
        }
    }
	//freeing resources
	for (i=0; i<=MAXFRAGS; i++) free(strngs[i]);
	for (i=0; i< MAXBONDS; i++) free(nsum[i]);
}


  string MCDLFormat::intToStr(int k) {
    char temp[16];
    sprintf(temp,"%d",k);
  string line=temp;
  return line;
  };



  string MCDLFormat::getMCDL(OBMol* pmol, bool includeCoordinates) {
    int  i=0;
    int  j=0;
    int  k=0;
    int  n=0;
    int  nt=0;
    int  nn=0;
    int  ntypes=0;
    int  charge=0;
    int  netcharge=0;
    int  netradical=0;
    int  chgarray[60][2];
    std::vector<int> ifragnum(MAXFRAGS);
    int  ltypes=0;
    int  ia0=0;
    int  ia1=0;
    int  chgflag=0;
    int  radicalflag=0;
    int  icons[96];
    int  newone=0;
    std::vector<int> nchg(MAXFRAGS);
    std::vector<int> nrad(MAXFRAGS);
    int naStore,nbStore;
    string constr="";
    string frag="";
    string fragment [MAXFRAGS];
    string line="";
    string semis="";
    int  comma;
    std::vector<int> ix(MAXFRAGS);
    int ia[MAXBONDS][4];
    string data="";
    int z[MAXBONDS][4];
    //FRAGCON
    int nHydr[MAXFRAGS];
    string aSymb [MAXFRAGS];
    int aCharge[MAXFRAGS];
    int aRadical[MAXFRAGS];
    std::vector<int> aNumber(MAXFRAGS);
    int bonds[MAXBONDS][4];
    int nConn[MAXFRAGS];
    //string s;
    int l;//,m;
    string astereo="";
    string bstereo="";
    string s1,s2;
    std::vector <int> eqList;
    std::vector <int>  bondStereoList;
    std::vector <int>  atomStereoList;
    std::vector <int>  anumStereo;
    string as1,as2,as3,as4;
    string linestereo;
    OBAtom *atom;
    OBBond *bond;

/* read in fragments and find different types */
/* ntypes: number of fragment types */
/* ntatoms: total number of fragments */
/* fragment[i]: atomic symbol string of fragment type i */
/* ifragnum[i]: number of fragments of type i */

  initGlobals();
    ntypes = 0;
    ntatoms = 0;
    for (i=0; i<MAXFRAGS; i++)
        ifragnum[i] = 0;

    chgflag = 0;
    netcharge = 0;
    radicalflag=0;
    netradical=0;

  netcharge=pmol->GetTotalCharge();
  netradical=pmol->GetTotalSpinMultiplicity()-1;

//Here FRAGCON conversion...
  pmol->DeleteHydrogens();
  createStereoLists(pmol,bondStereoList,atomStereoList,eqList);
//    return "started"; //passed
  naStore=pmol->NumAtoms();
  nbStore=pmol->NumBonds();
  for (i=1; i<=naStore; i++) {
    atom=pmol->GetAtom(i);
    nHydr[i-1]=atom->GetImplicitHCount()+atom->ExplicitHydrogenCount();
    aCharge[i-1]=atom->GetFormalCharge();
    aRadical[i-1]=atom->GetSpinMultiplicity();
    aSymb[i-1]=OBElements::GetSymbol(atom->GetAtomicNum());
    nConn[i-1]=atom->GetHvyDegree();
    aNumber[i-1]=i-1;
  };
  for (i=1; i<=nbStore; i++) {
    bond=pmol->GetBond(i-1);
    bonds[i-1][0]=bond->GetBeginAtomIdx();
    bonds[i-1][1]=bond->GetEndAtomIdx();
      bonds[i-1][2]=bond->GetBondOrder();  //type of bond store
      bonds[i-1][3]=i;                     //bond number in sm store
    };
    //FRAGCON
    i=0;
    while (i < naStore) {
      if ((nConn[i] == 1) && (nHydr[i] == 0)) {
        k=-1;
        for (j=0; j<nbStore; j++) if ((bonds[j][0] == (i+1)) || (bonds[j][1] == (i+1))) {
          k=j;
          break;
        }
        if (k>=0) {
          if (bonds[k][0] == (i+1)) nn=bonds[k][1]-1; else nn=bonds[k][0]-1;
          aCharge[nn]=aCharge[nn]+aCharge[i];
          aRadical[nn]=aRadical[nn]+aRadical[i];
          aSymb[nn]=aSymb[nn]+aSymb[i];
          for (j=i; j<(naStore-1); j++) {
            nHydr[j]=nHydr[j+1];
            aCharge[j]=aCharge[j+1];
            aRadical[j]=aRadical[j+1];
            aSymb[j]=aSymb[j+1];
            nConn[j]=nConn[j+1];
            aNumber[j]=aNumber[j+1];
          };
          for (j=0; j<nbStore; j++) {
            if (bonds[j][0] > (i+1)) bonds[j][0]--;
            if (bonds[j][1] > (i+1)) bonds[j][1]--;
          };
          for (j=k; j<(nbStore-1); j++) {
            bonds[j][0]=bonds[j+1][0];
            bonds[j][1]=bonds[j+1][1];
            bonds[j][2]=bonds[j+1][2];
            bonds[j][3]=bonds[j+1][3];
          };
          i--;
          naStore--;
          nbStore--;
        };
      };
      i++;
    }
//Ald without fragcon
  for (i=1; i<=naStore; i++) {
      charge=aCharge[i-1];
      if (charge != 0) chgflag = 1;
      if (aRadical[i-1] != 0) radicalflag=1;
      ntatoms++;
      newone = 1;
      frag=aSymb[i-1];
      j=nHydr[i-1];
      if (j>0) for (l=0; l<j; l++) {
        frag=frag+"H";
      };
      for (j=0; j<ntypes; j++) if (fragment[j] == frag) {
        ifragnum[j]++;
        newone = 0;
        break;
      }
      if (newone != 0) {
        fragment[ntypes]=frag;
        ifragnum[ntypes]++;
        ntypes++;
      }
    }

    ltypes=ntypes;
// order fragment types in ASCII order
    for (i=0; i<ntypes-1; i++) for (j=i+1; j<ntypes; j++) if (fragment[j] < fragment[i]) {
      frag=fragment[i];
      fragment[i]=fragment[j];
      fragment[j]=frag;
      k = ifragnum[i];
      ifragnum[i] = ifragnum[j];
      ifragnum[j] = k;
    }
// reread the fragments
// ix[i]: fragment type of fragment i
    for (i=1; i<=ntatoms; i++) {
      frag=aSymb[i-1];
      j=nHydr[i-1];
      if (j>0) for (l=0; l<j; l++) {
        frag=frag+"H";
      };
      for (j=0; j<ntypes; j++) if (frag == fragment[j]) ix[i-1] = j+1;
    }

//FRAGCON format atoms conversion
// read in the connections
// ia[i][0]: from fragment; ia[i][1]: to fragment of connection i
    nbonds = 0;
    for (i=1; i<=nbStore; i++) {
      // Possible, DRAGON format can contains multiply bonds between 2 atoms. Bond order?
      ia0=bonds[i-1][0];
      ia1=bonds[i-1][1];

      newone = 1;
      for (j=0; j<nbonds; j++)
          if ((ia[j][0] == ia0 && ia[j][1] == ia1) ||
              (ia[j][0] == ia1 && ia[j][1] == ia0))
          {
              newone = 0;
              break;
          }
      if (newone == 1) {
        ia[nbonds][0] = ia0;
        ia[nbonds][1] = ia1;
        nbonds=nbonds+1;
        ia[nbonds][0] = ia1;
        ia[nbonds][1] = ia0;
        nbonds=nbonds+1;
      }
    };
    for (i=0; i<nbonds; i++) {
        ia[i][2] = ix[ia[i][0]-1];
        ia[i][3] = ix[ia[i][1]-1];
    };
    for (i=0; i<nbonds; i++) for (j=0; j<4; j++) {
      z[i][j] = ia[i][j];
    };
    maxdepth = 1;

        // data ready for solution
    solve(ntypes,z,1);

    for (i=0; i<ntatoms; i++) ix[i] = qx[i];
    for (i=0; i<nbonds; i++) for (j=0; j<4; j++) ia[i][j] = qa[i][j];

//atoms stereo getting
    astereo=getAtomMCDL(pmol,ntatoms,ix,aNumber,atomStereoList,eqList);
//bond stereo getting
    bstereo=getBondMCDL(pmol,nbStore,ntatoms,ix,aNumber,bonds,bondStereoList,eqList);

    // original and final fragment numbers
    if (chgflag != 0) for (i=0; i<ntatoms; i++) nchg[ix[i]] = aCharge[i];
    if (radicalflag != 0) for (i=0; i<ntatoms; i++) nrad[ix[i]] = aRadical[i];

    // net charge
    data="";
    if (netcharge != 0)
    {
        if (netcharge > 0)
            data=data+intToStr(netcharge)+"+;";
        else data=data+intToStr(abs(netcharge))+"-;";
    }
    if (netradical != 0) data=data+intToStr(netradical)+"*;";
    // fragment portion of linear descriptor
    for (i=0; i<ltypes; i++) {
      if (i > 0) {
        data=data+';';
      }
      if (ifragnum[i] > 1) {
        data=data+intToStr(ifragnum[i]);
      }
      data=data+fragment[i];
    }
    constr="[";
    // connectivity portion of linear descriptor
    semis="";
    nt = 1;
    for (i=0; i<ntatoms; i++) {
      if (i > 0) {
        semis=semis+";";
      }
      n = 0;
      for (j=0; j<nbonds; j++) if (ia[j][2] == (i+1)) icons[n++] = ia[j][3];
      if (n > 1) {
        for (j=0; j<n-1; j++) for (k=j+1; k<n; k++) if (icons[k] < icons[j]) {
          nn = icons[j];
          icons[j] = icons[k];
          icons[k] = nn;
        }
      }
      comma = 0;
      for (j=0; j<n; j++) {
        if ((icons[j] > (i+1)) && (comma == 0)) {
          constr=constr+semis;
          constr=constr+intToStr(icons[j]);
          comma = 1;
          semis="";
        }  else if ((icons[j] > (i+1)) && (comma == 1)) {
          constr=constr+","+intToStr(icons[j]);
        }
      }
      nt++;
    };
    line=line+"]";

    constr=constr+line;
    data=data+constr;
    // fragment charges
    if (chgflag != 0) {
      nt = 0;
      for (n=-3; n<4; n++) for (i=0; i<ntatoms; i++) if ((nchg[ix[i]] == n) && (n != 0)) {
        chgarray[nt][0] = n;
        chgarray[nt++][1] = ix[i];
      }
      for (i=0; i<nt-1; i++) for (j=i; j<nt; j++) if (chgarray[i][0] > chgarray[j][0]) {
        k = chgarray[i][0];
        n = chgarray[i][1];
        chgarray[i][0] = chgarray[j][0];
        chgarray[i][1] = chgarray[j][1];
        chgarray[j][0] = k;
        chgarray[j][1] = n;
      }
      for (i=0; i<nt-1; i++) for (j=i; j<nt; j++) if ((chgarray[i][0] == chgarray[j][0]) && (chgarray[i][1] > chgarray[j][1])) {
        n = chgarray[i][1];
        chgarray[i][1] = chgarray[j][1];
        chgarray[j][1] = n;
      }
      line="";
      for (i=-1; i>-4; i--) {
        n = 0;
        for (j=0; j<nt; j++) if (chgarray[j][0] == i) n++;
        if (n > 0) {
          for (j=0; j<nt; j++) if (chgarray[j][0] == i) {
            if (chgarray[j][1] == 0)
              frag=""+intToStr(chgarray[j][1])+","+intToStr(abs(i))+"-";//Math.abs(i)+"-1";
            else
              frag=""+intToStr(chgarray[j][1])+","+intToStr(abs(i))+"-";//Math.abs(i)+"-"+chgarray[j][1];
            if (line.length() > 0) line=line+";";
            line=line+frag;
          }
        }
      }
      for (i=1; i<4; i++){
        n = 0;
        for (j=0; j<nt; j++) if (chgarray[j][0] == i) n++;
        if (n > 0) {
          for (j=0; j<nt; j++) if (chgarray[j][0] == i) {
            if (chgarray[j][1] == 0)
              frag=""+intToStr(chgarray[j][1])+","+intToStr(i)+"+";//i+"+1";
            else
              frag=""+intToStr(chgarray[j][1])+","+intToStr(i)+"+";
            if (line.length() > 0) line=line+";";
            line=line+frag;
          }
        }
      }
      if (ntypes>1) {
        constr=fchstart;//"{CZ:";
        constr=constr+line;
        constr=constr+"}";
        data=data+constr;
      };

    };
    //radical processing
    if (radicalflag != 0) {
      for (i=0; i<ntatoms; i++) {
        chgarray[i][0]=0;
        chgarray[i][1]=0;
      };
      nt = 0;
      for (n=1; n<3; n++) for (i=0; i<ntatoms; i++) if (nrad[ix[i]] == n) {
        chgarray[nt][0] = n;
        chgarray[nt++][1] = ix[i];
      }
      for (i=0; i<nt-1; i++) for (j=i; j<nt; j++) if (chgarray[i][0] > chgarray[j][0]) {
        k = chgarray[i][0];
        n = chgarray[i][1];
        chgarray[i][0] = chgarray[j][0];
        chgarray[i][1] = chgarray[j][1];
        chgarray[j][0] = k;
        chgarray[j][1] = n;
      }
      for (i=0; i<nt-1; i++) for (j=i; j<nt; j++) if ((chgarray[i][0] == chgarray[j][0]) && (chgarray[i][1] > chgarray[j][1])) {
        n = chgarray[i][1];
        chgarray[i][1] = chgarray[j][1];
        chgarray[j][1] = n;
      }
      line="";
      for (i=1; i<3; i++){
        n = 0;
        for (j=0; j<nt; j++) if (chgarray[j][0] == i) n++;
        if (n > 0) {
          for (j=0; j<nt; j++) if (chgarray[j][0] == i) {
            if (chgarray[j][1] == 0)
              frag=""+intToStr(chgarray[j][1])+","+intToStr(i)+"*";
            else
              frag=""+intToStr(chgarray[j][1])+","+intToStr(i)+"*";
            if (line.length() > 0) line=line+";";
            line=line+frag;
          }
        }
      }
      if (ntypes>1) {
        constr=fradstart;//"{CZ:";
        constr=constr+line;
        constr=constr+"}";
        data=data+constr;
      };

    };

//*****************************************************************************
    data=data+astereo;
    data=data+bstereo;
//*****************************************************************************
  /*   Block should be saved for future using

    if (includeCoordinates && (sm.fAtom.maxIndex() > 0))  {
      //Atom and bond blocks include
      data=data+fnastart+sm.fAtom.maxIndex()+"}";//};
      data=data+fnbstart+sm.fBond.maxIndex()+"}";//};
      data=data+fzcoorstart+"N}";
      //analizing minimal and maximal values of atomic coordinates
      xMin=sm.fAtom.getRX(1); xMax=sm.fAtom.getRX(1); yMin=sm.fAtom.getRY(1); yMax=sm.fAtom.getRY(1);
      for (i=1; i<=sm.fAtom.maxIndex(); i++) {
        if (sm.fAtom.getRX(i) < xMin) xMin=sm.fAtom.getRX(i);
        if (sm.fAtom.getRX(i) > xMax) xMax=sm.fAtom.getRX(i);
        if (sm.fAtom.getRY(i) < yMin) yMin=sm.fAtom.getRY(i);
        if (sm.fAtom.getRY(i) > yMax) yMax=sm.fAtom.getRY(i);
      };
      scale=0;
      if (xMax != xMin) scale=1/(xMax-xMin);
      if (yMax != yMin) if ((1/(yMax-yMin)) < scale) scale=1/(yMax-yMin);
      if (scale == 0) scale=1;
      line=fablockstart;
      nCharge=0;
      for (i=1; i<=sm.fAtom.maxIndex(); i++) {
        r=sm.fAtom.getRX(i);
        r=9.99*(r-xMin)*scale;
        s=CommonRout.str(r,0,2);
        line=line+s+",";
        r=sm.fAtom.getRY(i);
        r=9.99*(r-yMin)*scale;
        s=CommonRout.str(r,0,2);
        line=line+s+Atom.symbolofAtom(sm.fAtom.getNA(i));
        if (i < sm.fAtom.maxIndex()) line=line+";"; else line=line+"}";//};
        if (sm.fAtom.getNC(i) != 0) nCharge++;
      };
      data=data+line;
      line=fbblockstart;
      //two lines - adition for stereo handling
      test=false;
      linestereo=fstereobondstart;
      for (i=1; i<=sm.fBond.maxIndex(); i++) {
        n=sm.fBond.getTB(i);
        s="s";
        if (n > 3) n=1;
        if (n == 2) s="d"; else if (n == 3) s="t";
        s=""+sm.fBond.getAT(i,1)+s+sm.fBond.getAT(i,2);
        line=line+s;
        if (i < sm.fBond.maxIndex()) line=line+";"; else line=line+"}"; //};
        //stereoanalizing
        n=sm.fBond.getTB(i);
        if ((n >= 9) && (n <= 11)) {
          if (test) linestereo=linestereo+";";
          test=true;
          s="p";
          if (n == 10) s="o"; else if (n == 11) s="e";
          s=""+sm.fBond.getAT(i,1)+s+sm.fBond.getAT(i,2);
          linestereo=linestereo+s;
        };
      };
      data=data+line;
      if (test) {
        linestereo=linestereo+"}"; //};
        data=data+linestereo;
      };
      if (nCharge != 0) {
        line=fchargeblockstart+nCharge;
        for (i=1; i<= sm.fAtom.maxIndex(); i++) if (sm.fAtom.getNC(i) != 0) line=line+","+i+","+sm.fAtom.getNC(i);
        line=line+"}";
        data=data+line;
      };
    };
  */
    return data;
  };

  void MCDLFormat::restoreFullMCDL(string value, OBMol* pmol) {

  };

  void MCDLFormat::setMCDL(const string lineToParse, OBMol* pmol, string & sout) {

  std::vector <int> nH(MAXFRAGS);
  std::vector <int> nF(MAXFRAGS);
//  int nb;//,na;
  unsigned int i, j, nt, n1, n2;//,n3, nfrag;
  bool test;
  OBAtom sa;
  string mf="";
  string s="";
  string temp="";
  string ss="";
  string sstore="";
  int m,n,k;
  std::vector <int> nHydr(MAXFRAGS);
  string astereo="";
  string bstereo="";
  string chargestring="";
  string radicalstring="";
  int nPrev;
  string sa1="";
  string sa2="";
  string sF="";
  std::vector <int> charges(MAXFRAGS);
  std::vector <int> radicals(MAXFRAGS);
  std::vector <int> fragIndex(MAXFRAGS);
  std::vector <int> sbSpecial(MAXBONDS);
  std::vector <int> enumber(NELEMMAX);  //number of elements...
  std::vector <int> bondOrders(MAXBONDS);
  std::vector <double> rx(MAXFRAGS);
  std::vector <double> ry(MAXFRAGS);
  std::vector <string> names(MAXFRAGS);
  int netcharge=0;
  int netradical=0;
  unsigned int nelements=0;
  int kk;
  string value="";
  unsigned int acount, bcount, flags;

  std::vector <int> iA1(MAXBONDS);
  std::vector <int> iA2(MAXBONDS);
  std::vector <int> aPosition(MAXFRAGS);
  std::vector <int> aCharge(MAXFRAGS);
  std::vector <int> aRad(MAXFRAGS);
  std::vector <int> stereoBonds(MAXBONDS);

  value=lineToParse;
  sout="";

  if ((indexOf(value,fablockstart) >= 0) && (indexOf(value,fbblockstart) > 0)) {
    restoreFullMCDL(value,pmol);
    return;
  };
  for (i=0; i<MAXBONDS; i++) {
    stereoBonds[i]=0;
    sbSpecial[i]=0;
    iA1[i]=0;
    iA2[i]=0;
  };

  for (i=0; i<MAXFRAGS; i++) {
    charges[i]=0;
    radicals[i]=0;
    fragIndex[i]=0;
    //special[i]=0;
    nHydr[i]=0;
    names[i]="";
    aPosition[i]=6;
    aCharge[i]=0;
    aRad[i]=0;
  };
  //removing net charges and radiacal
  test=false;
  n=1;
  while ((value.length() > 0) && (value.at(0) >= '0') && (value.at(0) <= '9') && (! test) && (n>0)) {
    n=indexOf(value,";");
    if (n>0) {
      s=value.substr(0,n);
      n1=indexOf(s,"+");
      if (n1 < 0) n1=indexOf(s,"-");
      if (n1 < 0) n1=indexOf(s,"*");
      if (n1 < 0) {
        n=0;
        test=true;
      };
    };
    if (n > 0) {
      s=value.substr(0,n);
      if ((s.at(s.length()-1) == '+') || (s.at(s.length()-1) == '-')) {
        //total charge processing
        if (s.at(s.length()-1) == '+') netcharge=1; else netcharge=-1;
        s=s.substr(0,s.length()-1);
        n1=atoi(s.c_str());
        netcharge=netcharge*n1;
      } else if (s.at(s.length()-1) == '*') {
        //radical processing
        s=s.substr(0,s.length()-1);
        n1=atoi(s.c_str());
        netradical=n1;
      };
      value=value.substr(n+1,value.length());
      if (value.length() > 0) if (value.at(0) == ';'){
        value=value.substr(1);
      };
    };
  };

  //stereo string extraction

  n1=indexOf(value,fsastart);
  if (n1>=0) {
    n2=indexOf(value,"}",n1);
    if (n2>0) {
      astereo=""+value.substr(n1+fsastart.length(),n2);
      value=value.substr(0,n1)+value.substr(n2+1);
    }
  };
  n1=indexOf(value,fsbstart);
  if (n1>=0) {
    n2=indexOf(value,"}",n1);
    if (n2>0) {
      bstereo=""+value.substr(n1+fsbstart.length(),n2);
      value=value.substr(0,n1)+value.substr(n2+1);
    }
  };
  //charges processing
  n1=indexOf(value,fchstart);
  if (n1>=0) {
    n2=indexOf(value,"}",n1);
    if (n2>0) {
      chargestring=""+value.substr(n1+fchstart.length(),n2);
      value=value.substr(0,n1)+value.substr(n2+1);
    };
  };
  if ((netcharge != 0) && (chargestring == "")) {
    chargestring = "1," +intToStr(abs(netcharge));
    if (netcharge < 0) chargestring=chargestring+"-"; else chargestring=chargestring+"+";
  };
  //radical processing
  n1=indexOf(value,fradstart);
  if (n1>=0) {
    n2=indexOf(value,"}",n1);
    if (n2>0) {
      radicalstring=""+value.substr(n1+fradstart.length(),n2);
      value=value.substr(0,n1)+value.substr(n2+1);
    };
  };
  if ((netradical != 0) && (radicalstring == "")) {
    radicalstring= "1," + intToStr(abs(netcharge)) + "*";
  };

  n1=indexOf(value,"]");
  if (n1 > 0) {
    value=value.substr(0,n1);  //by unknown reason it was n1-1
  };

  //Atom array analizing
  test=true;
  nt=0;
  nelements=0;


  while (test) {
    n1=indexOf(value,";");
    n2=indexOf(value,"[");

    if (n1<0) {
      n1=n2;
    } else {
      if ((n2<n1) && (n2>=0)) n1=n2;
    };
    mf=value.substr(0,n1);
    value=value.substr(n1,value.length());
    if (value.at(0)==';') value=value.substr(1,value.length());
    n1=1;
    if ((mf.at(0)>='0') && (mf.at(0)<='9')) {
      test=true;
      temp="";
      while (test) {
        temp=temp+mf.at(0);
        mf=mf.substr(1,mf.length());
        test=(mf.at(0)>='0') && (mf.at(0)<='9');
      };
      n1=atoi(temp.c_str());
    };
    //do not change n1==nFrag!!!
    nt=nt+n1;
  //Conversion of XxHHHH to XxHn
    n2=lastIndexOf(mf,"HH");
    while (n2>0) {
      if (n2<(mf.length()-2)) {
        ss=mf.substr(mf.length()-1,mf.length());
        k=atoi(ss.c_str())+1;//Integer.parseInt(ss)+1;
      } else k=2;
      mf=mf.substr(0,n2+1)+intToStr(k);
      n2=lastIndexOf(mf,"HH");
    };
    //passed to this point...
    //End convsrsion
    n2=indexOf(mf,"H");
    temp="0";
    if (n2>0)  {
      if (n2<(mf.length()-1)) {
        temp=mf.substr(n2+1,mf.length());
      } else temp="1";
      mf=mf.substr(0,n2);
    };
    n2=atoi(temp.c_str());//Integer.parseInt(temp);
    //do not change n2 - number of hydrogens...

    nH[nelements]=n2;
    nF[nelements]=n1;
    names[nelements]=mf;
    nelements++;

    test=value.at(0) != '[';
  };

  sa.Clear();

  acount=nt;

  //parsing and analizing...
  for (i=0; i<MAXFRAGS; i++) charges[i]=0;
  if (chargestring != "") while (chargestring.length() > 0) {
    n1=indexOf(chargestring,";");
    if (n1>0) {
      s=chargestring.substr(0,n1);
      chargestring=chargestring.substr(n1+1);
    } else {
      s=chargestring;
      chargestring="";
    };
    n1=indexOf(s,",");
    if (n1 > 0) {
      sa1=s.substr(0,n1);
      s=s.substr(n1+1);
      n2=1;
      if (indexOf(s,"-") > 0) n2=-1;
      s=s.substr(0,s.length()-1);
      n1=atoi(s.c_str());
      n1=n1*n2;
      n2=atoi(sa1.c_str());
      charges[n2-1]=n1;
    };
  };
  for (i=0; i<MAXFRAGS; i++) radicals[i]=0;
  if (radicalstring != "") while (radicalstring.length() > 0) {
    n1=indexOf(radicalstring,";");
    if (n1>0) {
      s=radicalstring.substr(0,n1);
      radicalstring=radicalstring.substr(n1+1);
    } else {
      s=radicalstring;
      radicalstring="";
    };
    n1=indexOf(s,",");
    if (n1 > 0) {
      sa1=s.substr(0,n1);
      s=s.substr(n1+1);
      n2=1;
      if (indexOf(s,"-") > 0) n2=-1;
      s=s.substr(0,s.length()-1);
      n2=atoi(sa1.c_str());
      radicals[n2-1]=1;
    };
  };
//passed



  nt=0;
  ss="";
  bcount=0;

  for (i=0; i<nelements; i++) {
    s=names[i];
    n1=1;
    if (s.length()>1) if ((s.at(1)>='a') && (s.at(1)<='z')) n1=2;
    temp=s.substr(0,n1);
    if (n1<s.length()) sstore=s.substr(n1,s.length()); else sstore="";
    n1=nF[i];
    n2=OBElements::GetAtomicNum(temp.c_str());//Atom.positionofAtom(temp);
    nPrev=acount;

    for (j=1; j<=n1; j++) {
      fragIndex[nt]=nt+1;//sa.fragIndex=fragNo;
      nHydr[nt]=nH[i];      //sa.special=(byte)nH.getValue(i+1);
      nt++;

      aPosition[nt-1]=n2;

      if (sstore.length()>0) if (parseFormula(sstore,enumber)) { //bf=Formula.parseString(sstore);
        for (k=1; k<NELEMMAX; k++) if (enumber[k] > 0) for (m=1; m<=enumber[k]; m++) {
          fragIndex[acount]=nt+1;//sa.fragIndex=fragNo;  //Do not increment!
          nHydr[acount]=0;//special[acount]=0;
          aPosition[acount]=k;
          acount++;

          kk=hydrogenValency(k);//HVal[k];
          if (kk == 0) kk=1;
          bcount++;
          sbSpecial[bcount-1]=1;
          bondOrders[bcount-1]=kk;
          iA1[bcount-1]=nt-1;
          iA2[bcount-1]=acount-1;
        };
        //checking for nitrogen...
        //length analizing-extracting and addition of other atoms...
      };
      //charges analizing and restore

    if (charges[nt-1] != 0) {      //search for appropriate atom to put charge..
    if (nPrev == acount) {
      aCharge[nt-1]=charges[nt-1];
    } else  {
      if (charges[nt-1] < 0) {
        assignCharges(aPosition,iA1,iA2,aCharge,charges,bondOrders,8,nPrev,nt,acount,bcount);
        assignCharges(aPosition,iA1,iA2,aCharge,charges,bondOrders,7,nPrev,nt,acount,bcount);
        assignCharges(aPosition,iA1,iA2,aCharge,charges,bondOrders,16,nPrev,nt,acount,bcount);
        if (charges[nt-1] != 0) {
          aCharge[nt-1]=charges[nt-1];
        };
      } else {  //positive charge-at central atom
        aCharge[nt-1]=charges[nt-1];
      };
    };
  };
  if (radicals[nt-1] != 0) {      //search for appropriate atom to put charge..
    aRad[nt-1]=radicals[nt-1];
  };
  nPrev=acount;
};

  };

  //Bond analizing....
  value=value.substr(1,value.length());

  nt=0;
  ss="";
  test=true;
  while ((value.length()>0) && test) {
    n1=indexOf(value,";");
    n2=indexOf(value,"]");
    if (n1<0) n1=n2; else if ((n2<n1) && (n2>=0)) n1=n2;
    mf="";
    if ((n1>=0) && (value.length() > n1)) {
      mf=value.substr(0,n1);
      value=value.substr(n1+1,value.length());
    } else {
      test=false;
      mf=value;
      value="";
    };
    nt++;
    //parsing each fragment
    while (mf.length()>0) {
      n1=indexOf(mf,",");
      if (n1>=0) {
        s=mf.substr(0,n1);
        mf=mf.substr(n1+1,mf.length());
      } else {
        s=mf;
        mf="";
      };
      n1=atoi(s.c_str());
      bcount++;
      bondOrders[bcount-1]=0;
      iA1[bcount-1]=nt-1;
      iA2[bcount-1]=n1-1;
    };

    if (value.length()>0) if (value.at(0)==']') value=""; //end parsing
  };

  //Alternation
  alternate(aPosition,aCharge,aRad,nHydr,iA1,iA2,bondOrders,acount,bcount);
  generateDiagram(iA1,iA2,rx,ry,acount,bcount);
  for (i=0; i<bcount; i++) if (bondOrders[i] == 1) stereoBonds[i]=-1; //flags for bonds, which might be stereo
  implementAtomStereo(iA1,iA2,stereoBonds,rx,ry,acount,bcount,astereo);
  implementBondStereo(iA1,iA2,rx,ry,acount,bcount,bstereo);


  if (pmol != NULL) {
    for (i=0; i<acount; i++) {
      sa.Clear();
      sa.SetAtomicNum(aPosition[i]);
      if (aCharge[i] != 0) sa.SetFormalCharge(aCharge[i]);
      if (aRad[i] != 0) sa.SetSpinMultiplicity(1);
      sa.SetVector(rx[i],ry[i],0.0);
      pmol->AddAtom(sa);
    };
    for (i=0; i<bcount; i++) {
      flags=0;
      if (stereoBonds[i] == 1) flags=OB_WEDGE_BOND; else
      if (stereoBonds[i] == 2) flags=OB_HASH_BOND;
      pmol->AddBond(iA1[i]+1,iA2[i]+1,bondOrders[i],flags);
    };
  };
};

  void MCDLFormat::assignCharges(const std::vector <int> aPosition, const std::vector <int> iA1,
  const std::vector <int> iA2, std::vector <int>& aCharges, std::vector <int>& charges,
  std::vector <int>& bondOrder, int aPos, int nPrev, int nt, int acount, int bcount){
  //set negative charges
    int k,j,n;

  //return;

  for (k=nPrev; k<acount; k++) {
    if (aPosition[k] == aPos) {
        aCharges[k]=-1;
        charges[nt-1]=charges[nt-1]+1;
    for (j=0; j<bcount; j++) {
      if ((((iA1[j]+1) == nt) && (iA2[j] == k)) || ((iA1[j] == k) && ((iA2[j]+1) == nt))) {
      n=bondOrder[j];
            if (n > 1) {
              n=n-1;
        bondOrder[j]=n;
          };
          };
        };
    };
      if (charges[nt-1] == 0) break;
    };
  };


  int MCDLFormat::indexOf(const string instring, const string substring, int fromPos) {
    int result=instring.find(substring,fromPos);
    if (result == string::npos) result=-1;
    if (result >= instring.length()) result=-1;
      return result;
  };

  int MCDLFormat::lastIndexOf(const string instring, const string substring) {
    int result,n;
  bool test;

  result=-1;
  test=true;
  n=-1;
  while (test) {
    n=instring.find(substring,n+1);
    if (n == string::npos)
      test=false;
    else {
        result=n;
    };
  };
    //int result=instring.find_last_of(substring);
    //if  (result == string::npos) result=-1;
    //if  (result > (instring.length()-substring.length())) result=-1;
    return result;
  };

bool MCDLFormat::parseFormula(const string formulaString, std::vector <int>& enumber) {
  //vector<string> items;
  unsigned int i, n, k, n1, n2;//,j,nStart;
  string s;
  bool test;
  string asym;
  string value=formulaString;

  for (i = 0; i<NELEMMCDL; i++) enumber[i] = 0;

  for (i = 1; i<NELEMMCDL; i++) if (strlen(OBElements::GetSymbol(i)) == 2) {
      test=true;
    asym=OBElements::GetSymbol(i);
      while (test) {
        test=false;
        n=indexOf(value,asym);
        if (n>=0) {
          test=true;
          value=value.substr(0,n)+value.substr(n+asym.length(),value.length());
          k=1;
          if (n<value.length()) if ((value.at(n)>='0') && (value.at(n)<='9')) {
            n1=n;
            n2=n;
            while ((n2<(value.length()-1)) && (value.at(n2)>='0') && (value.at(n2)<='9')) n2++;
            if (! ((value.at(n2)>='0') && (value.at(n2)<='9'))) n2--;
            s=value.substr(n1,n2+1);
            k=atoi(s.c_str());
            value=value.substr(0,n1)+value.substr(n2+1,value.length());
          };
          enumber[i]=enumber[i]+k;
        };
      };
    };
  for (i = 1; i<NELEMMCDL; i++) if (strlen(OBElements::GetSymbol(i)) == 1) {
      test=true;
    asym=OBElements::GetSymbol(i);
      while (test) {
        test=false;
        n=indexOf(value,asym);
        if (n>=0) {
          test=true;
          value=value.substr(0,n)+value.substr(n+asym.length(),value.length());
          k=1;
          if (n<value.length()) if ((value.at(n)>='0') && (value.at(n)<='9')) {
            n1=n;
            n2=n;
            while ((n2<(value.length()-1)) && (value.at(n2)>='0') && (value.at(n2)<='9')) n2++;
            if (! ((value.at(n2)>='0') && (value.at(n2)<='9'))) n2--;
            s=value.substr(n1,n2+1);
      k=atoi(s.c_str());
            value=value.substr(0,n1)+value.substr(n2+1,value.length());
          };
          enumber[i]=enumber[i]+k;
        };
      };
    };
  return (value.length() == 0);
};

  string MCDLFormat::getMolTitle(string & line) {
    string::size_type n,k;
  string result;
  n=line.find(ftitlestart);
  if (n != string::npos) {
    k=line.find("}",n+ftitlestart.size());
    if (k != string::npos) {
    result=line.substr(n+ftitlestart.length(),k-n-ftitlestart.length());
    line=line.substr(0,n+1)+line.substr(k+1);
    }
  }
  return result;
  }


  ////////////////////////////////////////////////////

//Make an instance of the format class
MCDLFormat theMCDLFormat;

/////////////////////////////////////////////////////////////////

bool MCDLFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
{
  OBMol* pmol = pOb->CastAndClear<OBMol>();
  if(pmol==NULL)
      return false;

  istream& ifs = *pConv->GetInStream();

  pmol->BeginModify();

  int dim=0;
  pmol->SetDimension(dim);

  string line="";
    if(ifs.good()) getline(ifs, line);

  string molTitle=getMolTitle(line);

  if (molTitle.length() > 0) pmol->SetTitle(molTitle.c_str());
  if (line.length() > 0) setMCDL(line,pmol,molTitle);

  pmol->EndModify();

  return true;
}

////////////////////////////////////////////////////////////////

bool MCDLFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
  OBMol* pmol = dynamic_cast<OBMol*>(pOb);
  if(pmol==NULL) return false;

  std::ostream & ofs = *pConv->GetOutStream();

  //To use an output option
  string title=pmol->GetTitle();
  if (title.length() > 0) title=ftitlestart+title+"}";
  ofs << getMCDL(pmol,false) << title << endl;
  // prepareTest(pmol,ofs);
  //generateDiagram(pmol,ofs);

  return true; //or false to stop converting
}

} //namespace OpenBabel

