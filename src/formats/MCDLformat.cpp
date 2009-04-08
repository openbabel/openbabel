/**********************************************************************
Copyright (C) 2007 by Sergei V. Trepalin sergey_trepalin@chemical-block.com
Some portions Copyright (C) 2007 by Mike N. Burnett, burnettmn@ornl.gov
Code translation from Java

This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>
 
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
#include <openbabel/chiral.h>

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
        "Gakh A.A., Burnett M.N.,\n"
        "Modular Chemical Descriptor Language (MCDL):\n"
        "Composition, Connectivity and Supplementary Modules\n"
        "J.Chem.Inf.Comput.Sci, 2001, 41, 1494-1499"; 
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
    static const int MAXBONDS=300;
    static const int MAXFRAGS=200;
    static const int MAXCHARS=1000;
    static const int MAX_DEPTH=10;
    static const int NELEMMAX=120;

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

  //The following are temporarily in this file so that all MCDL stuff is in one file.
  //Normally in mcdlutil.h

  // Return valency by hydrogen for given atomic position in the Periodic Table
  int hydrogenValency(int na);
  int maxValency(int na);
  int alternate(OBMol * pmol, const int nH[], int bondOrders []);

  int alternate(const int aPosition[],const int aCharge[], const int aRad[], const int nHydr[],
                const int iA1[],const int iA2[],int bondOrders[], int nAtoms, int nBonds);

  int alternate(const std::vector<int> aPosition,const std::vector<int> aCharge, 
                const std::vector<int> aRad,const std::vector<int> nHydr, const std::vector<int> iA1,
                const std::vector<int> iA2, std::vector<int> & bondOrders, int nAtoms, int nBonds);


  void MCDLFormat::init(){
    /*
      for (int i=1; i<NELEMMAX; i++) HVal[i]=0;
      HVal[1]=1;
      HVal[5]=3;
      HVal[6]=4;
      HVal[7]=3;
      HVal[8]=2;
      HVal[9]=1;
      HVal[13]=3;
      HVal[14]=4;
      HVal[15]=3;
      HVal[16]=2;
      HVal[17]=1;
      HVal[32]=4;
      HVal[33]=3;
      HVal[32]=2;
      HVal[33]=1;
      HVal[50]=2;
      HVal[51]=3;
      HVal[52]=2;
      HVal[53]=1;
      HVal[82]=2;
      HVal[84]=2;
      HVal[85]=1;
      HVal[87]=1;
      HVal[88]=2;
    */	  
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
  }

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
  }



  /* Create descriptor connectivity string from connectivity table of fragments */
  string MCDLFormat::constring(int conntab [MAXBONDS][4], char * tstr)
  {
    int  i,j,k,n,nn,icons[6],comma;
    char line[82],semis[100];
    //char *sp = tstr;
    string result="";

    result="[";//snprintf(sp, 1, "[");
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
                snprintf(line, 82, "%s%d",semis,icons[j]);
                result=result+line;//strcat(sp,line);
                comma = 1;
                semis[0] = '\0';
              }
            else if (icons[j] > (i+1) && comma == 1)
              {
                snprintf(line, 82, ",%d",icons[j]);
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
    //int  nsum[MAXBONDS][MAXFRAGS];
    int *  nsum[MAXBONDS];
    bool lflag[MAXFRAGS];
    char strg[MAXFRAGS+1];
    //char strngs[MAXFRAGS][MAXFRAGS+1];
    char * strngs[MAXFRAGS+1];
    //std::vector <string> strngs[MAXFRAGS+1];
    char tstr[MAXFRAGS+1];
    //char *constring(int z[][4], char *);
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
    snprintf(temp, 16, "%d", k);
    string line = temp;
    return line;
  }

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
    int  ifragnum [MAXFRAGS];
    int  ltypes=0;
    int  ia0=0;
    int  ia1=0;
    int  chgflag=0;
    int  radicalflag=0;
    int  icons[96];
    int  newone=0;
    int  nchg[MAXFRAGS];
    int  nrad[MAXFRAGS];
    int naStore,nbStore;
    string constr="";//char constr[MAXCHARS];
    string frag="";//char frag[10];
    string fragment [MAXFRAGS];//char fragment[MAXFRAGS][10];
    string line="";//char line[82];
    string semis="";//char semis[100];
    int  comma;
    int ix[MAXFRAGS];
    //   int nH[MAXFRAGS];
    int ia[MAXBONDS][4];
    string data="";
    int z[MAXBONDS][4];
    //FRAGCON
    int nHydr[MAXFRAGS];
    string aSymb [MAXFRAGS];
    int aCharge[MAXFRAGS];
    int aRadical[MAXFRAGS];
    int aNumber[MAXFRAGS];
    int bonds[MAXBONDS][4];
    int nConn[MAXFRAGS];
    //string s;
    int l;//,m;
    //   int neighbours[4];
    //   int bondAccumulator[4];
    //   Listar priority=new Listar(4);
    //   int atn;
    string astereo="";
    string bstereo="";
    //Vector v=new Vector();
    string s1,s2;
    //   int nh,n1,n2,n3,n4,an1,an2,an3,an4,nCharge;
    //   bool testDouble;
    //    Listar eqList=new Listar(0);
    //    Listar bondStereoList=null;
    //    Listar atomStereoList=null;
    //    Listar anumStereo=new Listar(0);
    bool test;//, testParity;
    string as1,as2,as3,as4;
    //   double xMin,xMax,yMin,yMax,scale,r;
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

    //    sm.setMolfile(i_file);
    //    sm.defineConn();
    //    sm.determineFormula();

    //    if (! sm.singleFragment()) return "";  //Cannot detrmine MCDL descriptors for non-connected fragments
    /*
    //removing bad stereobonds-addition from 29 June 2006. This block was not tested in 200K test
    BondConnection bk=new BondConnection(Bond.nBondsMax);
    sm.defineBondConnection(bk);
    for (i=1; i<=pmol->_nbonds; i++) if ((pmol->.fBond.getTB(i) == 9) || (sm.fBond.getTB(i) == 10)) {
    test=true;
    n=sm.fBond.getAT(i,1);
    for (j=1; j<=bk.getNB(n); j++) {
    m=bk.getBN(n,j);
    if ((sm.fBond.getTB(m) == 2) || (sm.fBond.getTB(m) == 3)) test=false;
    };
    if (! test) sm.fBond.setTB(i,(byte)1);
    };
    //end addition from 29 June
    sm.allAboutCycles();
    */  
    pmol->FindChiralCenters();
    test=pmol->IsChiral();  //Instead of above lines...
    netcharge=pmol->GetTotalCharge();
    netradical=pmol->GetTotalSpinMultiplicity()-1;
    //return "almos start";  OKAY, passed

    //    for (i=1; i<=sm.fAtom.maxIndex(); i++) {
    //      netcharge=netcharge+sm.fAtom.getNC(i);
    //      netradical=netradical+sm.fAtom.getRL(i);
    //    };
    /*
    //checking if multiply stereocenters or acyclic double=bonds are present
    n=0;
    m=0;
    for (i=1; i<=sm.fBond.maxIndex(); i++) {
    if ((sm.fBond.getTB(i) == 9) || (sm.fBond.getTB(i) == 10)) n++;
    if ((sm.fBond.getTB(i) == 2) || (sm.fBond.getDB(i) <2)) m++;
    if ((m > 0) || (n > 1)) break;
    };
    if ((m > 0) || (n > 1)) ChainRotate.makeEquivalentList(sm,eqList,false);
    //Formation of bonds, which have stereo notations....
    bondStereoList=new Listar(sm.fBond.maxIndex());
    for (i=1; i<=sm.fBond.maxIndex(); i++) if ((sm.fBond.getTB(i) == 2) && (sm.fBond.getDB(i) < 2)) {
    //possible YES!
    n=sm.fBond.getAT(i,1);
    test=(sm.fAtom.getNB(n) > 1) && (sm.fAtom.getNB(n) <= 3);
    if (test) {
    n=sm.fBond.getAT(i,2);
    test=(sm.fAtom.getNB(n) > 1) && (sm.fAtom.getNB(n) <= 3);
    };
    n=sm.fBond.getAT(i,1);
    m=sm.fBond.getAT(i,2);
    if (test && (sm.fAtom.getNB(n) == 3)) {  //checking if equivalent substitutors are present....
    n1=0;
    n2=0;
    for (j=1; j<=sm.fAtom.getNB(n); j++) if (sm.fAtom.getAC(n,j) != m) {
    if (n1 == 0) n1=sm.fAtom.getAC(n,j); else n2=sm.fAtom.getAC(n,j);
    };
    if ((n1 > 0) && (n2 > 0)) {
    if (eqList.getValue(n1) == eqList.getValue(n2)) test=false;
    };
    };
    n=sm.fBond.getAT(i,2);
    m=sm.fBond.getAT(i,1);
    if (test && (sm.fAtom.getNB(n) == 3)) {  //checking if equivalent substitutors are present....
    n1=0;
    n2=0;
    for (j=1; j<=sm.fAtom.getNB(n); j++) if (sm.fAtom.getAC(n,j) != m) {
    if (n1 == 0) n1=sm.fAtom.getAC(n,j); else n2=sm.fAtom.getAC(n,j);
    };
    if ((n1 > 0) && (n2 > 0)) {
    if (eqList.getValue(n1) == eqList.getValue(n2)) test=false;
    };
    };
    if (test) bondStereoList.setValue(i,(byte)1);
    };  //end of double bond stereo list formation

    //Formation of atoms, which have stereo notation....
    atomStereoList=new Listar(sm.fAtom.maxIndex());
    for (i=1; i<=sm.fBond.maxIndex(); i++) if ((sm.fBond.getTB(i) == 9) || (sm.fBond.getTB(i) == 10)) {
    n=sm.fBond.getAT(i,1);
    atomStereoList.setValue(n,(byte)1);
    };
    */

    //Here FRAGCON conversion...
    pmol->DeleteHydrogens();
    //    return "started"; //passed	
    naStore=pmol->NumAtoms();
    nbStore=pmol->NumBonds();
    for (i=1; i<=naStore; i++) {
      atom=pmol->GetAtom(i);
      nHydr[i-1]=atom->ImplicitHydrogenCount()+atom->ExplicitHydrogenCount();
      aCharge[i-1]=atom->GetFormalCharge();
      aRadical[i-1]=atom->GetSpinMultiplicity();
      aSymb[i-1]=etab.GetSymbol(atom->GetAtomicNum());
      nConn[i-1]=atom->GetHvyValence();
      aNumber[i-1]=i;
    };
    for (i=1; i<=nbStore; i++) {
      bond=pmol->GetBond(i-1);
      bonds[i-1][0]=bond->GetBeginAtomIdx();
      bonds[i-1][1]=bond->GetEndAtomIdx();
      bonds[i-1][2]=bond->GetBondOrder();//sm.fBond.getTB(i);     //type of bond store
      bonds[i-1][3]=i;                     //bond number in sm store
    };
    //FRAGCON
    //	return "Bonds OK";  
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
      frag=aSymb[i-1];//Atom.symbolofAtom(sm.fAtom.getNA(i));
      j=nHydr[i-1];//sm.getNH(i);
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
    //z = new int [nbonds][4];//(RowArray *) malloc(nbonds*4*sizeof(int));
    for (i=0; i<nbonds; i++) for (j=0; j<4; j++) {
        z[i][j] = ia[i][j];
      };
    //qa = new int [nbonds][4];//(RowArray *) malloc(nbonds*4*sizeof(int));
    //qx = new int [ntatoms];//(int *) malloc(ntatoms*sizeof(int));
    maxdepth = 1;

    // data ready for solution
    solve(ntypes,z,1);
    //    z=null;//free(z);

    for (i=0; i<ntatoms; i++) ix[i] = qx[i];
    //    qx=null;//free(qx);
    for (i=0; i<nbonds; i++) for (j=0; j<4; j++) ia[i][j] = qa[i][j];
    //    qa=null;// free(qa);
    //atom symbols settings-will be used for stereodescriptors....
    /*     Stereo-handling code is below
           for (i=1; i<=sm.fAtom.maxIndex(); i++) sm.fAtom.setANUM(i,"");
           for (i=0; i<ntatoms; i++) {
           n=ix[i];        //new numeration
           atn=aNumber[i]; //numeration in sm
           sm.fAtom.setANUM(atn,""+n);
           };

           for (i=0; i<ntatoms; i++) {
           atn=aNumber[i]; //atom number in sm (initial molecule)
           if (atomStereoList.getValue(atn) != 0) {  //Contains stereo information
           n=ix[i];        //central atom, new numeration
           data=""+n;
           for (j=1; j<=priority.maxIndex(); j++) priority.setValue(j,0);
           as1=""; as2="";
           for (j=1; j<=4; j++) {
           s=getAtomSymbol(sm,atn,-1,j,"00");
           //Search through bonds...
           n1=-1;
           s1="";
           for (k=1; k<=sm.fBond.maxIndex(); k++) if ((sm.fBond.getAT(k,1) == atn) || (sm.fBond.getAT(k,2) == atn)) {
           if (sm.fBond.getAT(k,1) == atn) m=sm.fBond.getAT(k,2); else m=sm.fBond.getAT(k,1);
           s1=getAtomSymbol(sm,m);
           if (s.compareTo(s1) == 0) {
           n1=k;
           break;
           };
           };
           if (n1 > 0) priority.setValue(j,(byte)n1); else if (j<4) priority.setValue(j,-1);
           if (j <= 2) data=data+","+s;  else if (j == 3) as1=s; else as2=s;
           };
           try {
           s=""+priority.getValue(1)+" "+priority.getValue(2)+" "+priority.getValue(3)+" "+priority.getValue(4);
           m=analizeRSRotate(sm,atn,priority);
           if (m == 1) data=data+","+as2+","+as1; else data=data+","+as1+","+as2;
           v.addElement(data);
           anumStereo.append(atn);  //accumulation of atomic numbers in sm to check meso-forms...
           } catch (Exception ex) {
           ex.printStackTrace();
           }
           data=new String();
           };
           };
           //meso-isomer handling
           String ss="";
           boolean presentOtherStereo=false;
           if (anumStereo.maxIndex() > 1) for (i=1; i<anumStereo.maxIndex(); i++) if (anumStereo.getValue(i) > 0) {
           n=0;
           for (j=i+1; j<=anumStereo.maxIndex(); j++) if ((anumStereo.getValue(i) > 0) && (anumStereo.getValue(j) > 0)) if (eqList.getValue(anumStereo.getValue(i)) == eqList.getValue(anumStereo.getValue(j))) n++;
           if ((n % 2) != 1) presentOtherStereo=true;
           if (n>=1) {  //meso-forms search...
           s=(String)v.elementAt(i-1);
           k=analizeParity(s);
           testParity=false;
           for (j=i+1; j<=anumStereo.maxIndex(); j++) if ((anumStereo.getValue(i) > 0) && (anumStereo.getValue(j) > 0)) if (eqList.getValue(anumStereo.getValue(i)) == eqList.getValue(anumStereo.getValue(j))) {
           s1=(String)v.elementAt(j-1);
           if (! testParity) {
           m=analizeParity(s1);
           if (k*m == -1) testParity=true;
           };
           if (s.compareTo(s1) > 0) s=s1;
           };
           if (testParity) {
           if ((ss.length() == 0) || (ss.compareTo(s) > 0)) ss=s;
           for (j=i+1; j<=anumStereo.maxIndex(); j++) if ((anumStereo.getValue(i) > 0) && (anumStereo.getValue(j) > 0)) if ((anumStereo.getValue(i) > 0) && (anumStereo.getValue(j) > 0)) if (eqList.getValue(anumStereo.getValue(i)) == eqList.getValue(anumStereo.getValue(j))) anumStereo.setValue(j,-1);
           anumStereo.setValue(i,-1);
           } else {
           for (j=i+1; j<=anumStereo.maxIndex(); j++) if ((anumStereo.getValue(i) > 0) && (anumStereo.getValue(j) > 0)) if ((anumStereo.getValue(i) > 0) && (anumStereo.getValue(j) > 0)) if (eqList.getValue(anumStereo.getValue(i)) == eqList.getValue(anumStereo.getValue(j))) anumStereo.setValue(j,0);
           anumStereo.setValue(i,0);
           presentOtherStereo=true;
           };
           } else anumStereo.setValue(i,(byte)0);
           };
           //Do not analyze meso-forms for complex molecules....
           if (presentOtherStereo) ss="";

           if (ss.length() > 0) {  //some meso-isomers were found-analyze their parity and change if necessary
           k=analizeParity(ss);
           if (k == 1) for (i=1; i<=anumStereo.maxIndex(); i++) if (anumStereo.getValue(i) < 0) {  //all parities should be changed to satisfy..
           s=(String)v.elementAt(i-1);
           s1=changeParity(s);
           v.setElementAt(s1,i-1);
           };
           };
           //end meso-isomers handling
           //sorting vectors...
           if (v.size() > 1) for (i=0; i<(v.size() - 1); i++) for (j=i+1; j<v.size(); j++) {
           s1=(String)v.elementAt(i);
           s2=(String)v.elementAt(j);
           if (compareStringsNumbers(s1,s2) > 0) {
           v.setElementAt(s1,j);
           v.setElementAt(s2,i);
           };
           };
           if (v.size() > 0) {
           astereo=fsastart;
           for (i=0; i < v.size(); i++) {
           s1=(String)v.elementAt(i);
           while (s1.indexOf("zz") > 0) {  //Substitite temporary lowest-priority 00 to zero
           n=s1.indexOf("zz");
           s1=s1.substring(0,n)+"0"+s1.substring(n+2);
           };
           while (s1.indexOf("00") > 0) {  //Substitite temporary lowest-priority 00 to zero
           n=s1.indexOf("00");
           s1=s1.substring(0,n)+"0"+s1.substring(n+2);
           };
           if (i > 0) astereo=astereo+";";
           s1=removeZeroeth(s1);
           astereo=astereo+s1;
           };
           astereo=astereo+"}";
           };
           //    if ((flab != null) && (astereo != null)) flab.setText(astereo);
           //stereobonds analizing
           v.removeAllElements();
           anumStereo.redim(0);
           for (i=0; i<nbStore; i++) if (bonds[i][2] == 2) {
           n=bonds[i][3]; //Old bond number
           testDouble=bondStereoList.getValue(n) == 1;  //bondStereoList is filled on input-stereonotation can be used for acyclic double-bond with non-equivalent substitutors....
           if (testDouble) {
           n1=ix[bonds[i][0]-1];
           n2=ix[bonds[i][1]-1];
           if ( n1>n2) {
           k=bonds[i][0];
           bonds[i][0]=bonds[i][1];
           bonds[i][1]=k;
           };
           //1-st atom has minimal number now. Search for minimal number, attached to first
           n1=ix[bonds[i][0]-1];
           n2=ix[bonds[i][1]-1];

           data=""+n1+"d"+n2;
           n1=-1;  //to 1-st atom
           n2=-1;
           n3=-1;  //to 2-nd atom
           n4=-1;
           an1=0; an2=0; an3=0; an4=0;
           // search for the attached atoms...
           for (j=0; j<nbStore; j++) if (i != j) {
           if ((bonds[j][0] == bonds[i][0]) || (bonds[j][1] == bonds[i][0])) {
           if (n1 == -1) n1=j; else n2=j;
           if (bonds[j][0] == bonds[i][0]) {
           if (an1 == 0) an1=ix[bonds[j][1]-1]; else an2=ix[bonds[j][1]-1];
           } else {
           if (an1 == 0) an1=ix[bonds[j][0]-1]; else an2=ix[bonds[j][0]-1];
           };
           };
           if ((bonds[j][0] == bonds[i][1]) || (bonds[j][1] == bonds[i][1])) {
           if (n3 == -1) n3=j; else n4=j;
           if (bonds[j][0] == bonds[i][1]) {
           if (an3 == 0) an3=ix[bonds[j][1]-1]; else an4=ix[bonds[j][1]-1];
           } else {
           if (an3 == 0) an3=ix[bonds[j][0]-1]; else an4=ix[bonds[j][0]-1];
           };
           };
           };
           if ((an1 > 0) && (an2 > 0)) if (an1 > an2) {
           k=n1; n1=n2; n2=k;
           k=an1; an1=an2; an2=k;
           };
           //Now n1,n2 contains bond numbers in bonds array, connected to bonds[i][0] atom. n3,n4-bond numbers, connected to bonds[i][1] atom.
           //I have to use corresponding elements bonds[n1][3] to determine bond numbers in initial molecule sm.
           if (an1 > 0) as1=""+an1; else {
           //Hydrogens MUST not be...
           as1=getAtomSymbol(sm,aNumber[bonds[i][0]-1],aNumber[bonds[i][1]-1],1,"zz");
           };
           if (an2 > 0) as2=""+an2; else {
           as2=getAtomSymbol(sm,aNumber[bonds[i][0]-1],aNumber[bonds[i][1]-1],2,"zz");
           if (as2.compareTo("zz") == 0) {
           as2="00";
           };
           };
           if (an3 > 0) as3=""+an3; else {
           //Hydrogens MUST not be...
           as3=getAtomSymbol(sm,aNumber[bonds[i][1]-1],aNumber[bonds[i][0]-1],1,"zz");
           };
           if (an4 > 0) as4=""+an4; else {
           as4=getAtomSymbol(sm,aNumber[bonds[i][1]-1],aNumber[bonds[i][0]-1],2,"zz");
           if (as4.compareTo("zz") == 0) {
           as4="00";
           };
           };
           if (n1 < 0) {  //bond is connected with undefined atom - search for it....
           an1=aNumber[bonds[i][0]-1];
           an2=aNumber[bonds[i][1]-1]; //excluded atom...
           for (j=1; j<=sm.fBond.maxIndex(); j++) {
           an3=sm.fBond.getAT(j,1);
           an4=sm.fBond.getAT(j,2);
           if ((an3 == an1) || (an4 == an1)) {
           if (an3 == an1) k=an4; else k=an3;
           if (k != an2) {
           s=getAtomSymbol(sm,k);
           if (s.compareTo(as1) == 0) n1=j;
           };
           };
           if (n1 > 0) break;
           };
           //Searching for anum!
           } else n1=bonds[n1][3];
           if (n3 < 0) {  //bond is connected with undefined atom - search for it....
           an1=aNumber[bonds[i][1]-1];
           an2=aNumber[bonds[i][0]-1]; //excluded atom...
           for (j=1; j<=sm.fBond.maxIndex(); j++) {
           an3=sm.fBond.getAT(j,1);
           an4=sm.fBond.getAT(j,2);
           if ((an3 == an1) || (an4 == an1)) {
           if (an3 == an1) k=an4; else k=an3;
           if (k != an2) {
           s=getAtomSymbol(sm,k);
           if (s.compareTo(as3) == 0) n3=j;
           };
           };
           if (n3 > 0) break;
           };
           //Searching for anum!
           } else n3=bonds[n3][3];
           if ((n1 < 0) || (n3 < 0)) {
           //my error in logic. raise an exception
           //or do nothing...
           data="";
           } else {
           //sign determine...
           //Then I can to determine scalar product
           //and angle between bonds n1 and n3
           k=sproduct(sm,bonds[i][3],n1,n3);
           if (as2.compareTo("00") == 0) {
           data=data+","+as2;
           if (k == 1) {  //cis
           data=data+","+as4;
           data=data+","+as3;//an3;
           data=data+","+as1;
           } else if (k == 2) {  //trans
           data=data+","+as3;
           data=data+","+as4;
           data=data+","+as1;
           } else data="";  //collinear

           } else {
           data=data+","+as1;//an1;
           if (k == 1) {  //cis
           data=data+","+as3;//an3;
           data=data+","+as4;
           data=data+","+as2;
           } else if (k == 2) {  //trans
           data=data+","+as4;
           data=data+","+as3;
           data=data+","+as2;
           } else data="";  //collinear
           };
           };
           if (data != null) if (data.length() > 0) {
           v.addElement(data);
           anumStereo.append(n);
           };
           data=new String();
           }
           };
           //Identical double-bonds search, like cis,trans-hexa-2,4-diene
           if (anumStereo.maxIndex() > 1) for (i=1; i<anumStereo.maxIndex(); i++) if (anumStereo.getValue(i) > 0) {
           n=0;
           for (j=i+1; j<=anumStereo.maxIndex(); j++) if (bondEquivalent(anumStereo.getValue(i),anumStereo.getValue(j),eqList,sm)) n++;
           if (n>=1) {  //bond searching...
           s=(String)v.elementAt(i-1);
           k=analizeParityBond(s);

           testParity=false;
           for (j=i+1; j<=anumStereo.maxIndex(); j++) if (bondEquivalent(anumStereo.getValue(i),anumStereo.getValue(j),eqList,sm)) {
           s1=(String)v.elementAt(j-1);
           if (! testParity) {
           m=analizeParityBond(s1);
           if (k*m == -1) testParity=true;
           };
           if (compareStringsNumbers(s,s1)>0) s=s1;
           };
           if (testParity) {
           k=analizeParityBond(s);
           if (k == 1) {  //all parities should be changed to satisfy..
           s=(String)v.elementAt(i-1);
           s1=changeParityBond(s);
           v.setElementAt(s1,i-1);
           for (j=i+1; j<=anumStereo.maxIndex(); j++) if (bondEquivalent(anumStereo.getValue(i),anumStereo.getValue(j),eqList,sm)) {
           s=(String)v.elementAt(j-1);
           s1=changeParityBond(s);
           v.setElementAt(s1,j-1);
           };
           };
           };
           // set all anumStereo processed to 0 to indicate, that no processing is required
           for (j=i+1; j<=anumStereo.maxIndex(); j++) if (bondEquivalent(anumStereo.getValue(i),anumStereo.getValue(j),eqList,sm)) anumStereo.setValue(j,(byte)0);
           };
           anumStereo.setValue(i,(byte)0);
           }; //end symmetrical bonds

           if (v.size() > 1) for (i=0; i<(v.size() - 1); i++) for (j=i+1; j<v.size(); j++) {
           s1=(String)v.elementAt(i);
           s2=(String)v.elementAt(j);
           if (compareStringsNumbers(s1,s2) > 0) {
           v.setElementAt(s1,j);
           v.setElementAt(s2,i);
           };
           };
           if (v.size() > 0) {   //stereobonds present
           bstereo=fsbstart;
           for (i=0; i < v.size(); i++) {
           s1=(String)v.elementAt(i);
           while (s1.indexOf("zz") > 0) {  //Substitite temporary lowest-priority 00 to zero
           n=s1.indexOf("zz");
           s1=s1.substring(0,n)+"0"+s1.substring(n+2);
           };
           while (s1.indexOf("00") > 0) {  //Substitite temporary lowest-priority 00 to zero
           n=s1.indexOf("00");
           s1=s1.substring(0,n)+"0"+s1.substring(n+2);
           };
           if (i > 0) bstereo=bstereo+";";
           s1=removeZeroeth(s1);
           bstereo=bstereo+s1;
           };
           bstereo=bstereo+"}";
           } else bstereo=null;
    */
    // Output

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

    /*
    //*****************************************************************************
    if (astereo != null) data=data+astereo;
    if (bstereo != null) data=data+bstereo;
    //*****************************************************************************

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
  }

  void MCDLFormat::restoreFullMCDL(string value, OBMol* pmol) {

  }

  void MCDLFormat::setMCDL(const string lineToParse, OBMol* pmol, string & sout) {

    std::vector <int> nH(MAXFRAGS);
    std::vector <int> nF(MAXFRAGS);
    //  int nb;//,na;
    int i,j,nt,n1,n2;//,n3, nfrag;
    bool test;
    OBAtom sa;
    string mf="";
    string s="";
    string temp="";
    string ss="";
    string sstore="";
    int m,n,k;
    //  bool testReplace;
    std::vector <int> nHydr(MAXFRAGS);//  int nHydr[MAXFRAGS];
    //int nHydr[MAXFRAGS];
    string astereo="";
    string bstereo="";
    string chargestring="";
    string radicalstring="";
    int nPrev;//,atn,bn;
    string sa1="";
    string sa2="";
    string sF="";
    std::vector <int> charges(MAXFRAGS);
    std::vector <int> radicals(MAXFRAGS);
    std::vector <int> fragIndex(MAXFRAGS);
    //int special[MAXFRAGS]; - nHydr instead using...
    std::vector <int> sbSpecial(MAXBONDS);
    std::vector <int> enumber(NELEMMAX);  //number of elements...
    std::vector <int> bondOrders(MAXBONDS);
    std::vector <string> names(MAXFRAGS);
    int netcharge=0;
    int netradical=0;
    //  bool test1;
    //  bool threeCoor;
    int nelements=0;
    int kk;
    string value="";
    int acount,bcount;

    std::vector <int> iA1(MAXBONDS);
    std::vector <int> iA2(MAXBONDS);
    std::vector <int> aPosition(MAXFRAGS);
    std::vector <int> aCharge(MAXFRAGS);
    std::vector <int> aRad(MAXFRAGS);

    value=lineToParse;
    sout="";

    if ((indexOf(value,fablockstart) >= 0) && (indexOf(value,fbblockstart) > 0)) {
      restoreFullMCDL(value,pmol);
      return;
    };
    for (i=0; i<MAXBONDS; i++) {
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
      chargestring="1,"+abs(netcharge);
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
      radicalstring="1,"+abs(netcharge)+'*';
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

    //for (i=0; i<nt; i++) pmol->AddAtom(sa);
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
    //acount=pmol->NumAtoms();
    bcount=0;

    for (i=0; i<nelements; i++) {
      s=names[i];
      n1=1;
      if (s.length()>1) if ((s.at(1)>='a') && (s.at(1)<='z')) n1=2;
      temp=s.substr(0,n1);
      if (n1<s.length()) sstore=s.substr(n1,s.length()); else sstore="";
      n1=nF[i];
      int iso=0;
      n2=etab.GetAtomicNum(temp.c_str(),iso);//Atom.positionofAtom(temp);
      nPrev=acount;

      //sout=sout+"NA="+intToStr(n2)+" NF="+intToStr(n1)+" frag="+sstore+" temp="+temp+" ! ";
      for (j=1; j<=n1; j++) {
        fragIndex[nt]=nt+1;//sa.fragIndex=fragNo;
        nHydr[nt]=nH[i];      //sa.special=(byte)nH.getValue(i+1);
        nt++;

        //saptr=pmol->GetAtom(nt);
        //saptr->Clear();
        //saptr->SetAtomicNum(n2);// sa.nA=(byte)n2;
        aPosition[nt-1]=n2;

        if (sstore.length()>0) if (parseFormula(sstore,enumber)) { //bf=Formula.parseString(sstore);
            for (k=1; k<NELEMMAX; k++) if (enumber[k] > 0) for (m=1; m<=enumber[k]; m++) {
                  fragIndex[acount]=nt+1;//sa.fragIndex=fragNo;  //Do not increment!
                  nHydr[acount]=0;//special[acount]=0;
                  aPosition[acount]=k;
                  acount++; 

                  //sa.Clear();
                  //sa.SetAtomicNum(k);// sa.nA=bf.eList[k];
                  //if (pmol != NULL) pmol->AddAtom(sa);
          
                  kk=hydrogenValency(k);//HVal[k];
                  if (kk == 0) kk=1;
                  //if (pmol != NULL) pmol->AddBond(nt,acount,kk);
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
            //saptr=pmol->GetAtom(nt);
            //saptr->SetFormalCharge(charges[nt-1]);
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
              //saptr=pmol->GetAtom(nt);
              //saptr->SetFormalCharge(charges[nt-1]);
              aCharge[nt-1]=charges[nt-1];
            };
          };
        };
        if (radicals[nt-1] != 0) {      //search for appropriate atom to put charge..
          //saptr=pmol->GetAtom(nt);
          //saptr->SetSpinMultiplicity(1);
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
        //if (pmol != NULL) pmol->AddBond(nt,n1,1);
        bcount++;
        bondOrders[bcount-1]=0;
        iA1[bcount-1]=nt-1;
        iA2[bcount-1]=n1-1;
      };

      if (value.length()>0) if (value.at(0)==']') value=""; //end parsing
    };

    //Alternation
    //alternate(aPosition,aCharge,aRad,nHydr,iA1,iA2,bondOrders);
    alternate(aPosition,aCharge,aRad,nHydr,iA1,iA2,bondOrders,acount,bcount);
    if (pmol != NULL) {
      for (i=0; i<acount; i++) {
        sa.Clear();
        sa.SetAtomicNum(aPosition[i]);
        if (aCharge[i] != 0) sa.SetFormalCharge(aCharge[i]);
        if (aRad[i] != 0) sa.SetSpinMultiplicity(1);
        pmol->AddAtom(sa);
      };
      for (i=0; i<bcount; i++) pmol->AddBond(iA1[i]+1,iA2[i]+1,bondOrders[i]);
    };
  
    //pmol using....
    /*
      test=true;
      while (test) {
      test=BondAlternate.isolateSixMembered(result,nHydr);
      test1=BondAlternate.alternateSinglets(result,nHydr,null);
      if (test1) test=true;
      };

      //new StructureDlg("",result);

      SimpleMolecule smInitial=new SimpleMolecule();  //is required for future alternations
      smInitial.moleculeCopy(result);

      BondAlternate.alternate(result,null,nHydr);
      test=true;
      while (test) {
      test=BondAlternate.isolateSixMembered(result,nHydr);
      test1=BondAlternate.alternateSinglets(result,nHydr,null);
      if (test1) test=true;
      };

      test=BondAlternate.checkCorrect(result,nHydr);

      if (! test) {
      result.moleculeCopy(smInitial);
      BondAlternate.specialAlternate(result,nHydr);
      test=true;
      while (test) {
      test=BondAlternate.isolateSixMembered(result,nHydr);
      test1=BondAlternate.alternateSinglets(result,nHydr,null);
      if (test1) test=true;
      };
      test=BondAlternate.checkCorrect(result,nHydr);
      if (! test) {
      result.moleculeCopy(smInitial);
      test=BondAlternate.allAlternate(result,nHydr);
      test=true;
      while (test) {
      test=BondAlternate.isolateSixMembered(result,nHydr);
      test1=BondAlternate.alternateSinglets(result,nHydr,null);
      if (test1) test=true;
      };
      };
      };
      //remove Special flag
      for (i=1; i<=result.fAtom.maxIndex(); i++) result.fAtom.setSpecial(i,(byte)0);
      for (i=1; i<=result.fBond.maxIndex(); i++) result.fBond.setSpecial(i,(byte)0);
      result.defineConn();
      result.determineFormula();
    */


    //Stereo atoms analizing...
    /*
      if (astereo != null) {
      ss=astereo;
      astereo=addZeroeth(astereo,"-1");
      while (astereo.length() >0 ) {
      s="";
      n1=astereo.indexOf(";");
      if (n1 > 0) {
      s=astereo.substring(0,n1);
      if (n1 < (astereo.length()-1)) astereo=astereo.substring(n1+1); else astereo="";
      } else {
      s=astereo;
      astereo="";
      };
      //analize s
      if (s != null) if (s.length() > 0) {
      //save data in Astereo... Bond reconfiguring is possible only after re-drawing and chain rotations
      CommonRout.RemoveSpaces(s);
      n1=s.indexOf(",");
      if (n1 > 0) {
      threeCoor=(s.indexOf("-1") > 0);

      temp=s.substring(0,n1);
      s=s.substring(n1+1);
      atn=Integer.parseInt(temp);
      n1=s.indexOf(",");
      if (n1 > 0) s=s.substring(n1+1);
      n1=s.indexOf(",");
      if (n1 > 0) s=s.substring(n1+1);
      n1=s.indexOf("H");
      if (n1<0) {
      n1=s.indexOf("0");
      if (n1 > 0) n1=s.indexOf(",0"); //occupy last position
      };
      if (n1 == 0) {
      result.fAtom.setAStereo(atn,(byte)1);  //set "R" configuration - priorities are taken from atomic numbers...
      } else if (n1 > 0) {
      result.fAtom.setAStereo(atn,(byte)2);  //set S configuration   - priorities are take from atomic numbers
      } else {
      //Priority analizing
      n1=s.indexOf(",");
      if (n1 > 0) {
      temp=s.substring(0,n1);
      s=s.substring(n1+1);
      if (temp.compareTo("0") == 0) temp="zz";
      if (s.compareTo("0") == 0) s="zz";
      if (threeCoor) {
      if (compareStringsNumbers(temp,s) < 0) result.fAtom.setAStereo(atn,(byte)1); else result.fAtom.setAStereo(atn,(byte)2);
      } else {
      if (compareStringsNumbers(temp,s) > 0) result.fAtom.setAStereo(atn,(byte)1); else result.fAtom.setAStereo(atn,(byte)2);
      };
      };
      };
      };
      };
      };
      };
    */
    //Stereo bonds analizing...
    /*
      String s1;
      if (bstereo != null) {
      ss=bstereo;
      bstereo=addZeroeth(bstereo,"0");
      while (bstereo.length() >0 ) {
      s="";
      n1=bstereo.indexOf(";");
      if (n1 > 0) {
      s=bstereo.substring(0,n1);
      if (n1 < (bstereo.length()-1)) bstereo=bstereo.substring(n1+1); else bstereo="";
      } else {
      s=bstereo;
      bstereo="";
      };
      //analize s
      if (s != null) if (s.length() > 0) {
      //save data in Bstereo... Bond reconfiguring is possible only after re-drawing
      CommonRout.RemoveSpaces(s);
      n1=s.indexOf(",");
      if (n1 > 0) {
      temp=s.substring(0,n1);
      s=s.substring(n1+1);
      n1=temp.indexOf("d");
      if (n1 < 0) n1=temp.indexOf("D");
      if (n1 > 0) {
      n2=Integer.parseInt(temp.substring(0,n1));
      n1=Integer.parseInt(temp.substring(n1+1));
      //search for bond...
      bn=0;
      for (i=1; i<=result.fBond.maxIndex(); i++) if (((result.fBond.getAT(i,1) == n1) && (result.fBond.getAT(i,2) == n2)) || ((result.fBond.getAT(i,1) == n2) && (result.fBond.getAT(i,2) == n1))) {
      bn=i;
      break;
      };
      n1=-1;
      n2=-1;
      k=s.indexOf(",");
      sF="";
      if (k > 0)  {
      sF=s.substring(0,k);
      s=s.substring(k+1);  //removing 1-st fragment, added to an1
      };
      k=s.indexOf(",");
      sa1=""; sa2="";
      if (k > 0) {                    //and analizing of order of two remaining fragments, connected to an2
      sa1=s.substring(0,k);
      s=s.substring(k+1);
      try {
      n1=Integer.parseInt(sa1);
      if (n1 == 0) {
      n1=-1;
      sa1="00";//"zz";
      };
      } catch (Exception ex) {
      n1=-1;
      };
      };
      k=s.indexOf(",");
      if (k > 0) {
      sa2=s.substring(0,k);
      s=s.substring(k+1);
      try {
      n2=Integer.parseInt(sa2);
      if (n2 == 0) {
      n2=-1;
      sa2="00";//"zz";
      };
      } catch (Exception ex) {
      n2=-1;
      };
      };
      if (bn > 0) {
      if ((n1 < 0) && (n2 < 0)) {
      if (compareStringsNumbers(sa1,sa2) > 0) k=1; else k=2;
      } else {
      if (n1 < 0) k=1; else //E
      if (n2 < 0) k=2; else { //Z
      if (n1 < n2) k=2; else k=1;
      };
      };
      if (sF.compareTo("0") == 0) k=3-k;
      result.fBond.setBStereo(bn,(byte)k);
      if (flab != null) {
      s1=flab.getText();
      s1=s1+" "+n1+" "+n2+" "+s;
      flab.setText(s1);
      };
      };
      };
      };
      };
      };
      };
    */
  }

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
  }


  int MCDLFormat::indexOf(const string instring, const string substring, int fromPos) {
    int result=instring.find(substring,fromPos);
    if (result == string::npos) result=-1; 
    if (result >= instring.length()) result=-1;
    return result;
  }

  int MCDLFormat::lastIndexOf(const string instring, const string substring) {
    int result,n;
    bool test;

    result=-1;
    test=true;
    n=-1;
    while (test) {
      n=instring.find(substring,n+1);
      if (n == string::npos) test=false; else {
        result=n;
      };
    };
    //int result=instring.find_last_of(substring);
    //if  (result == string::npos) result=-1;
    //if  (result > (instring.length()-substring.length())) result=-1;
    return result;
  }

  bool MCDLFormat::parseFormula(const string formulaString, std::vector <int>& enumber) {
    //vector<string> items;
    int i,n,k,n1,n2;//,j,nStart;
    string s;
    bool test;
    string asym;
    string value=formulaString;
    
    for (i=0; i<etab.GetNumberOfElements(); i++) enumber[i]=0;

    for (i=1; i<etab.GetNumberOfElements(); i++) if (strlen(etab.GetSymbol(i))==2) {
        test=true;
        asym=etab.GetSymbol(i);
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
    for (i=1; i<etab.GetNumberOfElements(); i++) if (strlen(etab.GetSymbol(i))==1) {
        test=true;
        asym=etab.GetSymbol(i);
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
  }

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
    if(pmol==NULL)
      return false;

    ostream& ofs = *pConv->GetOutStream();

    //To use an output option
    string title=pmol->GetTitle();
    if (title.length() > 0) title=ftitlestart+title+"}";
    ofs << getMCDL(pmol,false) << title << endl;

    return true; //or false to stop converting
  }

  //**************************************************************
  //The following are temporarily in this file so that all MCDL stuff is in one file.
  //Normally in and mcdlutil.cpp

#define NELEMMCDL 121
  //Hydrogen valencies. Zero dummy element is the first 
	const int hVal[NELEMMCDL] = {  
    0,1,0,0,0,3,4,3,2,1,
    0,0,0,3,4,3,2,1,0,0,
    0,0,0,0,0,0,0,0,0,0,
    0,0,4,3,2,1,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,
    2,3,2,1,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,
    0,0,2,0,2,1,0,1,2,0,
    0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,1,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0};
	const int maxVal[NELEMMCDL] = {
    0,1,0,1,2,4,4,5,3,1,
    0,1,2,4,4,6,6,7,0,1,
    2,3,4,5,6,7,6,4,4,2,
    2,3,4,5,6,7,8,1,2,3,
    4,5,6,7,8,6,6,2,2,3,
    4,5,6,7,8,1,2,3,4,4,
    3,3,3,3,3,4,3,3,3,3,
    3,3,4,5,6,7,8,6,6,3,
    2,3,4,5,6,7,8,1,2,3,
    4,5,6,6,6,6,3,4,3,3,
    3,3,1,1,1,0,0,0,0,0,
    0,8,1,8,5,0,0,0,0,0,0};



  string intToStr(int k) {
    char temp[16];
    snprintf(temp, 16, "%d", k);
    string line = temp;
    return line;
  }


  int hydrogenValency(int na) {   //Hydrogen valency
    int result=0;
    if (na < NELEMMCDL) result=hVal[na];
    return result;
  }

  int maxValency(int na) {        //Maximal valency of a specified element
    int result=8;
    if (na < NELEMMCDL) result=maxVal[na];
    return result;
  }

  static int findAlternateSinglets(const std::vector<int>iA1, const std::vector<int>iA2, const std::vector<int> nH,
                                   const std::vector<int> hydrogenValency, std::vector<int>& bondOrder, int nAtoms, int nBonds) {
    //return values: 0 - no singlet found, 1 - found and assigned singlet, 2 - assigned singlet only with special

    std::vector<int>nUnassigned(nAtoms);
    std::vector<int>nAssigned(nAtoms);
    std::vector<int>bNumber(nAtoms);
    int i,n,k;
    int result;

    result=0;
    for (i=0; i<nAtoms; i++) {
      nUnassigned[i]=0;
      nAssigned[i]=0;
    };
    for (i=0; i<nBonds; i++) {
      n=bondOrder[i];
      if (n == 0) {
        nUnassigned[iA1[i]]=nUnassigned[iA1[i]]+1;
        nUnassigned[iA2[i]]=nUnassigned[iA2[i]]+1;
        bNumber[iA1[i]]=i;
        bNumber[iA2[i]]=i;
      } else {
        nAssigned[iA1[i]]=nAssigned[iA1[i]]+n;
        nAssigned[iA2[i]]=nAssigned[iA2[i]]+n;
      };
    };
	
    for (i=0; i<nAtoms; i++) if ((hydrogenValency[i] > 0) && (nUnassigned[i] == 1)) {
        n=bNumber[i];
        k=hydrogenValency[i]-(nH[i]+nAssigned[i]);
        if (k < 1) {
          if (k == 0) {  //increment H valence by 2
            k=2;
          } else {
            k=1;
          }; 
          result=2;
        };
        if (k > 3) {
          k=3;
          result=3;
        };
        bondOrder[n]=k;
        if (result == 0) {
          result=1;
        };
      };
    return result;
  }

  static void makeAssignment(const std::vector<int> iA1, const std::vector<int> iA2, const std::vector<int> nH, 
                             const std::vector<int> hydrogenValency,	const std::vector<int> bondAssignment, const std::vector<int> specialFlag, 
                             std::vector<int>& bondOrder, int nAtoms, int nBonds, int & nAss) {
   
    int i,k;

    nAss=0;
    for (i=0; i<nBonds; i++) if (bondOrder[i] == 0) {
        bondOrder[i]=bondAssignment[nAss]+specialFlag[i];
        nAss++;
        k=1;
        while (k != 0) {
          k=findAlternateSinglets(iA1,iA2,nH,hydrogenValency,bondOrder,nAtoms,nBonds);
        };
      };
  }

  static bool analyzeOK(const std::vector<int> iA1,const std::vector<int> iA2, const std::vector<int> nH, 
                        const std::vector<int> hydrogenValency,	const std::vector<int> maxValency, const std::vector<int> bondOrder, 
                        const std::vector<int> atomCheckFlag, int nAtoms, int nBonds, int & nGtMax, int & nNEH, int & nOddEven, 
                        bool testExceedHydrogen, bool oddEvenCheck) {
  
    std::vector<int>nBondsValency(nAtoms);  //dynamically allocation
    int i,k;
    bool result;
  
    nGtMax=0;
    nNEH=0;
    nOddEven=0;
    for (i=0; i<nAtoms; i++) nBondsValency[i]=0;
    for (i=0; i<nBonds; i++) {
      nBondsValency[iA1[i]]=nBondsValency[iA1[i]]+bondOrder[i];
      nBondsValency[iA2[i]]=nBondsValency[iA2[i]]+bondOrder[i];
    };
    for (i=0; i<nAtoms; i++) if (atomCheckFlag[i] == 1) {
        if (nBondsValency[i] > maxValency[i]) nGtMax++;
        if (testExceedHydrogen) {
          if ((nH[i]+nBondsValency[i]) != hydrogenValency[i]) nNEH++;
        } else {
          if ((hydrogenValency[i] > 0) && (nH[i] > 0)) if ((nH[i]+nBondsValency[i]) != hydrogenValency[i]) nNEH++;
          if ((hydrogenValency[i] > 0) && (nH[i] ==0)) if (nBondsValency[i] < hydrogenValency[i]) nNEH++;
          if (oddEvenCheck) {
            k=nH[i]+nBondsValency[i];
            if ((k % 2) != (maxValency[i] % 2)) nOddEven++;
          };
        };
      };
    result=((nGtMax == 0) && (nNEH == 0) && (nOddEven == 0));
    return result;
  }

  static bool incrementAssignment(std::vector<int>& bondAssignment, int nAss) {
    int i,j;
    bool result;
 
    result=false;
    for (i=nAss-1; i>=0; i--) if (bondAssignment[i] == 1) {
        bondAssignment[i]++;
        if (i < (nAss-1)) for (j=i+1; j<nAss; j++) bondAssignment[j]=1;
        result=true;
        return result;
      };
    return result;
  }

  static int determineBondsOrder(const std::vector<int> iA1, const std::vector<int> iA2, 
                                 const std::vector<int> nH, const std::vector<int> maxValency,std::vector<int>& bondOrder, 
                                 std::vector<int>& hydrogenValency, int nAtoms, int nBonds, bool oddEvenViolate) {
    //On input BortOrder has to be initialized. Real bond orders have to be putted. 0 means this order should be determined
    //MaxValency and HydrogenValency and NH are required only for those atoms, which are included in BondOrder=0
  
    std::vector<int>nNeighbour(nAtoms);
    std::vector<int>bondAssignment(nBonds);  //Should be Max(Atoms, Bonds);
    std::vector<int>bondOrderStore(nBonds);
    std::vector<int>specialFlag(nBonds);

    int i,j,k,n,nAss;
    bool  test;
    int nGtMax;  //number of atoms, for which maximal valency is exausted
    int nNEH;    //bumber of atoms, for which number of calc and desirable hydrogens are disagree
    int nOddEven;//Parity violation - atom with odd maximal valency has even valency and back
    int result;

    k=nAtoms;
    if (nBonds > k) k=nBonds;
    k++;

    result=0;
    for (i=0; i<nAtoms; i++) nNeighbour[i]=0;
    for (i=0; i<nBonds; i++) {
      nNeighbour[iA1[i]]=nNeighbour[iA1[i]]+1;
      nNeighbour[iA2[i]]=nNeighbour[iA2[i]]+1;
      specialFlag[i]=0;
    };
    //Special flag-find allene and cumulene
    for (i=0; i<nBonds; i++) {
      if ((nNeighbour[iA1[i]] == 2) && (nNeighbour[iA2[i]] == 2) && (hydrogenValency[iA1[i]] == 4)
          && (hydrogenValency[iA2[i]] == 4) && (nH[iA1[i]] == 0) && (nH[iA2[i]] == 0)) specialFlag[i]=1;
    };

    //correct hydrogen valency by 2 if required....
    for (i=0; i<nAtoms; i++) if (hydrogenValency[i] > 0) {
        n=nH[i]+nNeighbour[i];
        if (n > hydrogenValency[i]) hydrogenValency[i]=hydrogenValency[i]+2;
        if (n > hydrogenValency[i]) hydrogenValency[i]=hydrogenValency[i]+2;
      };
    //Label all bonds which have order 1 exactly
    for (i=0; i<nAtoms; i++) if ((hydrogenValency[i] > 0) && (nH[i] > 0)) {
        n=nH[i]+nNeighbour[i];
        if (n == hydrogenValency[i]) for (j=0; j<nBonds; j++) if (((iA1[j] == i) || (iA2[j] == i)) && (bondOrder[j] == 0)) {
              bondOrder[j]=1;
            };
      };

    k=1;
    while (k != 0) {
      k=findAlternateSinglets(iA1,iA2,nH,hydrogenValency,bondOrder,nAtoms,nBonds);
    };
    //repeat default no. hydrogens correction
    for (i=0; i<nAtoms; i++) nNeighbour[i]=0;
    for (i=0; i<nBonds; i++) {
      k=bondOrder[i];
      if (k == 0) k=1;
      nNeighbour[iA1[i]]=nNeighbour[iA1[i]]+k;
      nNeighbour[iA2[i]]=nNeighbour[iA2[i]]+k;
    };
    //correct hydrogen valency by 2 if required....


    for (i=0; i<nAtoms; i++) if (hydrogenValency[i] > 0) {
        n=nH[i]+nNeighbour[i];
        if (n > hydrogenValency[i]) hydrogenValency[i]=hydrogenValency[i]+2;
        if (n > hydrogenValency[i]) hydrogenValency[i]=hydrogenValency[i]+2;
      };

    //below array NNeighbour is used for picking atoms, for which valency have to be checked..
    for (i=0; i<nAtoms; i++) nNeighbour[i]=0;
    for (i=0; i<nBonds; i++) if (bondOrder[i] == 0) {
        nNeighbour[iA1[i]]=1;
        nNeighbour[iA2[i]]=1;
      };
    test=false;
    for (i=0; i<nBonds; i++) {
      bondAssignment[i]=1;
      bondOrderStore[i]=bondOrder[i];
      if (bondOrder[i] == 0) test=true;   //bad bonds exists-must to assign...
    };
    if (! test) return result;    //All fine-not necessary especially alternate
    test=false;
    while (! test) {              //Try to assign with fina valencies
      for (i=0; i<nBonds; i++) bondOrder[i]=bondOrderStore[i];
      makeAssignment(iA1,iA2,nH,hydrogenValency,bondAssignment,specialFlag,bondOrder,nAtoms,nBonds,nAss);
      test=analyzeOK(iA1,iA2,nH,hydrogenValency,maxValency,bondOrder,nNeighbour,nAtoms,nBonds,nGtMax,nNEH,nOddEven,true,true);
      if (! test) {
        test=! incrementAssignment(bondAssignment,nAss);
        if (test) result=1;
      };
    };
    if (result == 1) {        //try to assign with hydrogen valency exceeded
      result=0;
      for (i=0; i<nBonds; i++) bondAssignment[i]=1;
      test=false;
      while (! test) {
        for (i=0; i<nBonds; i++) bondOrder[i]=bondOrderStore[i];
        makeAssignment(iA1,iA2,nH,hydrogenValency,bondAssignment,specialFlag,bondOrder,nAtoms,nBonds,nAss);
        test=analyzeOK(iA1,iA2,nH,hydrogenValency,maxValency,bondOrder,nNeighbour,nAtoms,nBonds,nGtMax,nNEH,nOddEven,false,true);
        if (! test) {
          test=! incrementAssignment(bondAssignment,nAss);
          if (test) result=1;
        };
      };
    };
    if (result == 1) {        //hydrogen valency is still bad. Try to analyze withoud Odd/Even checking (metals, etc).
      result=0;
      for (i=0; i<nBonds; i++) bondAssignment[i]=1;
      test=false;
      while (! test) {
        for (i=0; i<nBonds; i++) bondOrder[i]=bondOrderStore[i];
        makeAssignment(iA1,iA2,nH,hydrogenValency,bondAssignment,specialFlag,bondOrder,nAtoms,nBonds,nAss);
        test=analyzeOK(iA1,iA2,nH,hydrogenValency,maxValency,bondOrder,nNeighbour,nAtoms,nBonds,nGtMax,nNEH,nOddEven,false,false);
        if (! test) {
          test=! incrementAssignment(bondAssignment,nAss);
          if (test) result=1;
        };
      };
    };
    return result;
  }

  int alternate(const std::vector<int> aPosition,const std::vector<int> aCharge, 
                const std::vector<int> aRad,const std::vector<int> nHydr, const std::vector<int> iA1,
                const std::vector<int> iA2, std::vector<int> & bondOrders, int nAtoms, int nBonds) {

    std::vector<int> hVal(nAtoms);
    std::vector<int> maxVal(nAtoms);
    int i,result;

    for (i=0; i<nAtoms; i++) {
      hVal[i]=hydrogenValency(aPosition[i]);
      if (hVal[i] > 0) {
        if (aRad[i] != 0) hVal[i]=hVal[i]-1;  //returns 0 for non-radical and 2 for radical
        if (aPosition[i] == 5) hVal[i]=hVal[i]-aCharge[i];           //B 
        else if (aPosition[i] == 6) hVal[i]=hVal[i]-abs(aCharge[i]); //C
        else hVal[i]=hVal[i]+aCharge[i];  //Heteroatoms
        if (hVal[i] < 0) hVal[i]=0;
      };
      maxVal[i]=maxValency(aPosition[i]);
      if (aCharge[i] != 0) maxVal[i]=maxVal[i]+1;
    };

    result=determineBondsOrder(iA1,iA2,nHydr,maxVal,bondOrders,hVal,nAtoms,nBonds,true);

    return result;
  }

  int alternate(OBMol * pmol, const std::vector<int> nH, std::vector<int>& bondOrders) {
    //This overloaded method does not work. By unknown reason I cannot extract 
    //connection matrix from pmol

    //number of hydogens must be filled on input.
    std::vector<int> hVal(pmol->NumAtoms());
    std::vector<int> maxVal(pmol->NumAtoms());
    std::vector<int>iA1(pmol->NumBonds());
    std::vector<int>iA2(pmol->NumBonds());
    int nAtoms,nBonds;
    int i,k,na;
    OBAtom * sa;
    OBBond * sb;
    int result=0;

    pmol->AssignSpinMultiplicity();

    nAtoms=pmol->NumAtoms();
    nBonds=pmol->NumBonds();

    for (i=0; i<nBonds; i++) {
      sb=pmol->GetBond(i);
      iA1[i]=sb->GetBeginAtomIdx()-1;
      iA2[i]=sb->GetEndAtomIdx()-1;
    };


    for (i=1; i<=nAtoms; i++) {
      sa=pmol->GetAtom(i);
      na=sa->GetAtomicNum();
      hVal[i-1]=hydrogenValency(na);
      if (hVal[i-1] > 0) {
        if (sa->GetSpinMultiplicity() != 0) hVal[i-1]=hVal[i-1]-1;  //returns 0 for non-radical and 2 for radical
        k=sa->GetFormalCharge();
        if (sa->IsHeteroatom()) hVal[i-1]=hVal[i-1]+k;
        else if (na == 6) hVal[i-1]=hVal[i-1]-abs(k);
        else hVal[i-1]=hVal[i-1]-k;
        if (hVal[i-1] < 0) hVal[i-1]=0;
      };
      maxVal[i-1]=maxValency(na);
      if (sa->GetFormalCharge() != 0) maxVal[i-1]=maxVal[i-1]+1;
    };

    result=determineBondsOrder(iA1,iA2,nH,maxVal,bondOrders,hVal,nAtoms,nBonds,true);
    for (i=0; i<nBonds; i++) {
      sb=pmol->GetBond(i);
      sb->SetBondOrder(bondOrders[i]);
    }; 
    return result;
  }

} //namespace OpenBabel

