#ifndef OELIB_CHAINS_H
#define OELIB_CHAINS_H

#define MaxMonoAtom 20
#define MaxMonoBond 20

namespace OpenEye {

class OEAtom;
class OEMol;

typedef struct 
{
    int flag;
    short elem, count;
    int n1, n2, n3, n4;
} Template;

class OEChainsParser
{
public:

	OEChainsParser(void);
	~OEChainsParser(void);

	bool PerceiveChains(OEMol &);

private: // methods

	bool  DetermineHetAtoms(OEMol &);
	bool  DetermineConnectedChains(OEMol &);
	bool  DeterminePeptideBackbone(OEMol &);
	bool  DeterminePeptideSidechains(OEMol &);
	bool  DetermineNucleicBackbone(OEMol &);
	bool  DetermineNucleicSidechains(OEMol &);
	bool  DetermineHydrogens(OEMol &);

	void  SetupMol(OEMol &);
    void  SetResidueInformation(OEMol &);
    void  ClearResidueInformation(OEMol &);
	void  CleanupMol(void);

	void  AssignResidue(OEMol &, int, int, int);
	int   IdentifyResidue(void *, OEMol &, int, int); // ByteCode *

	void  DefineMonomer(void **, int, char *); // ByteCode **
	int   IdentifyElement(char *);

	bool  MatchConstraint(OEAtom *, int);
	bool  Match2Constraints(Template *, OEAtom *, OEAtom *);
	bool  Match3Constraints(Template *, OEAtom *, OEAtom *, OEAtom *);
	bool  Match4Constraints(Template *, OEAtom *, OEAtom *, OEAtom *, OEAtom *);

	void  ConstrainBackbone(OEMol &, Template *, int);

	int   RecurseChain(OEMol &, int, int);
	void  TraceNucleicChain(OEMol &, int, int); 
	void  TracePeptideChain(OEMol &, int, int);

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
#endif // OELIB_CHAINS_H

