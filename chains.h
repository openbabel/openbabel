#ifndef OBLIB_CHAINS_H
#define OBLIB_CHAINS_H

#define MaxMonoAtom 20
#define MaxMonoBond 20

namespace OpenBabel {

class OBAtom;
class OBMol;

typedef struct 
{
    int flag;
    short elem, count;
    int n1, n2, n3, n4;
} Template;

class OBChainsParser
{
public:

	OBChainsParser(void);
	~OBChainsParser(void);

	bool PerceiveChains(OBMol &);

private: // methods

	bool  DetermineHetAtoms(OBMol &);
	bool  DetermineConnectedChains(OBMol &);
	bool  DeterminePeptideBackbone(OBMol &);
	bool  DeterminePeptideSidechains(OBMol &);
	bool  DetermineNucleicBackbone(OBMol &);
	bool  DetermineNucleicSidechains(OBMol &);
	bool  DetermineHydrogens(OBMol &);

	void  SetupMol(OBMol &);
    void  SetResidueInformation(OBMol &);
    void  ClearResidueInformation(OBMol &);
	void  CleanupMol(void);

	void  AssignResidue(OBMol &, int, int, int);
	int   IdentifyResidue(void *, OBMol &, int, int); // ByteCode *

	void  DefineMonomer(void **, int, char *); // ByteCode **
	int   IdentifyElement(char *);

	bool  MatchConstraint(OBAtom *, int);
	bool  Match2Constraints(Template *, OBAtom *, OBAtom *);
	bool  Match3Constraints(Template *, OBAtom *, OBAtom *, OBAtom *);
	bool  Match4Constraints(Template *, OBAtom *, OBAtom *, OBAtom *, OBAtom *);

	void  ConstrainBackbone(OBMol &, Template *, int);

	int   RecurseChain(OBMol &, int, int);
	void  TraceNucleicChain(OBMol &, int, int); 
	void  TracePeptideChain(OBMol &, int, int);

	char *ParseSmiles(char *, int);

private: // members

	void *PDecisionTree; // ByteCode *
	void *NDecisionTree; // ByteCode *

    int   ResMonoAtom[MaxMonoAtom];
    int   ResMonoBond[MaxMonoBond];

	unsigned short *bitmasks;
	unsigned char  *resids;
	unsigned char  *flags;
	bool           *hetflags;
	short          *atomids;
	short          *resnos;
	short          *sernos;
	char           *hcounts;
	char           *chains;
};

}
#endif // OBLIB_CHAINS_H

