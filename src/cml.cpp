/**********************************************************************
Copyright (C) 2002-2003 Peter Murray-Rust.
Some portions Copyright (c) 2003 by Geoffrey R. Hutchison

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

/** NOTE.
The intention of this software is to provide a definitive toolkit
for reading and writing CML documents (or documents which contain
CML fragments).

CML itself is defined by the published specification(s) and must not be changed
or extended - this is a key aspect of markup languages. If you do wish to create
documents with functionality beyond CML you will need to develop your
own namespace.

It is my intention to support developments in CML and this code for
enhancement of the functionality and precision of CML. This activity
MUST be coordinated; uncoordinated development of information protocols
leads to chaos. We have already seen documents purporting to be "CML"
which are not even well-formed XML. Some others do not adhere to the
published DTD or Schema. In the spirit of OpenSource I wish to encourage
innovation and collaboration and ask that if you wish to develop and
extend the functionality of CML code, please contact me at pm286@cam.ac.uk.

It is extremely important that "CML" code processes and produces
conformant CML (which also implies conformant XML). There are many subtleties
in XML which must be addressed, such as:ignorable whitespace, standalone and
default values, encodings, entity referencing rules. CDATA sections,
namespace scope, character encodings, and MANY aspects of XML Schema such as
inheritance, schema locations, default qualification, nillability,
datatyping, etc. If you do not understand these you are likely to break
conformance with CML.

The following requests are addressed to anyone or organization making changes to this code
who has not contacted me. The GPL prevents them being mandatory but they should be
respected by anyone who believes in interoperability

The code itself is OpenSource. You may make any changes under
the terms of the GPL above.

HOWEVER: ANY alteration to the code means that your system may deviate from
a compliant CML system. CML has been trademarked precisely to protect its implementation.
(The International Union of Crystallography has trademarked their CIF
specification, so this policy is in keeping with accepted practice in
Open Specifications).

YOU are therefore requested to announce that your code has been modified from that deposited
by myself or supplied from the current OpenSource site. You should be very careful
in any claims that it complies with any of the published
CML specifications. You are requested to make it clear to users that your use of
"CML" the <cml> tag and the CML namespaces posted on http://www.xml-cml.org has
been made without the knowledge of the authors of CML.

This notice must be included with any modified software.

Peter Murray-Rust, 2002, 2003
*/

#include "mol.h"

# include <string>

# include <time.h>
/* ---- Size of time-string ---- */
#define TIME_STR_SIZE 64

using namespace std;
namespace OpenBabel {

// default molecule size
int ATOMSIZE = 100;

// ----------------lexical processing ------------------
// this is a crude approximate to a SAX2-parser which omits all DTD events
// the XML parsing should be contained within this
bool ReadXML(istream &ifs);
// private lexical processing for start tag
string startTag(string s);
// private function. processes all <...>. Does not do comments or DTDs. Should do
void tag(string s);
// gets attribute
string getAttribute(vector <pair<string, string> > &atts, string name);
// checks attributes for validity and namespaces
void processAttributes(vector <pair<string, string> > &atts);
// processes namespaces
void processNamespace(string name, string value);
// checks for valid names
bool isXMLName(string s);
// split concatenated attributes
void splitAttributes(string s, vector <pair <string, string> > &atts);
// interludes
void endElement(string s);
// interlude
void startElement(string name, vector<pair<string,string> > &atts);
// check elements against allowed list
void makeAllowedElementLists();
// check attributes against allowed list
void makeAllowedAttributeLists();
// find any undeclared attributes
vector <string> getUnknownAttributes(vector <string> &allowed, vector <pair <string, string> > &atts);
// convert XML entities
string processXMLEntities(string s);
// output XML entities
string escapeXMLEntities(string s);
// output PCDATA
bool writePCDATA(ostream &ofs, string value);
// manage prefixed namespaces
pair <string, string> getNamespacePair(string name);
// split a list of floats
void processFloatTokens(vector <double> &v, vector<double>::size_type n, string att);
// split a list of ints
void processIntegerTokens(vector <int> &v, vector<int>::size_type n, string att);
// split a list of ints
void processStringTokens(vector <string> &v, vector<string>::size_type n, string att);

// ----------------SAX2-like callbacks ------------------
// SAX2 is an OpenSource API for event-based processing of XML documents.
// See http://www.saxproject.org and SourceForge. David Brownell has written
// an O'Reilly book on SAX2.
// NOTE: There is NO agreed SAX2 API in any language other than Java.
// This binding is created by direct "translation" of the Java API
// which is included.

// void characters(char[] ch, int start, int length)
//          Receive notification of character data.
void characters(string length);
// void endDocument()
//          Receive notification of the end of a document.
void endDocument();
// void endElement(java.lang.String namespaceURI, java.lang.String localName, java.lang.String qName)
//          Receive notification of the end of an element.
void endElement(string namespaceURI, string localName, string qName);
// void endPrefixMapping(java.lang.String prefix)
//          End the scope of a prefix-URI mapping.
void endPrefixMapping(string prefix);
// void ignorableWhitespace(char[] ch, int start, int length)
//          Receive notification of ignorable whitespace in element content.
void ignorableWhitespace(char* ch, int start, int length);
// void processingInstruction(java.lang.String target, java.lang.String data)
//          Receive notification of a processing instruction.
void processingInstruction(string target, string data);
// void setDocumentLocator(Locator locator)
//          Receive an object for locating the origin of SAX document events.
// void setDocumentLocator() NOT implemented
// void skippedEntity(java.lang.String name)
//          Receive notification of a skipped entity.
void skippedEntity(string name);
// void startDocument()
//          Receive notification of the beginning of a document.
void startDocument();
// void startElement(java.lang.String namespaceURI, java.lang.String localName, java.lang.String qName, Attributes atts)
//          Receive notification of the beginning of an element.
void startElement(string namespaceURI, string localName, string qName, vector <pair<string,string> > &atts);
// void startPrefixMapping(java.lang.String prefix, java.lang.String uri)
//          Begin the scope of a prefix-URI Namespace mapping.
void startPrefixMapping(string prefix, string uri);

// ----------------XML utilities ----------------

// writes foo="bar" (preceded by space)
bool writeAttribute(ostream& ofs, string attname, string value);
// writes <foo or <bar:foo
bool writeStartTagStart(ostream& ofs, string name);
// writes >
bool writeStartTagEnd(ostream& ofs);
// writes </foo> or </bar:foo>
bool writeEndTag(ostream& ofs, string name);
// writes />
bool writeCombinedTagEnd(ostream& ofs);
// normalizes string
string getNormalizedString(char* ch);
// normalizes string
string getNormalizedString(string s);

// ----------------utilities ----------------
// trims whitespace from start and end of string
string trim(string s);
// to upper
string toUpperCase(string s);
// to lower
string toLowerCase(string s);
// is a string in a vector of strings
bool isInStringVector(vector <string> v, string s);
// list of unused elementNames
vector <string> unusedElementNameVector;
// gets datetime in s
bool getTimestr(string& s);

// ----------------CML elements----------------
const string CML1_ELEMENT_NAMES = "angle atom atomArray bond bondArray coordinate2 coordinate3 crystal electron feature float floatArray floatMatrix integer integerArray link list molecule reaction sequence string stringArray torsion";

const char* CML1_NAMESPACE = "http://www.xml-cml.org/dtd/cml_1_0_1.dtd";
const char* CML2_NAMESPACE = "http://www.xml-cml.org/schema/cml2/core";
const char* STMML_NAMESPACE = "http://www.xml-cml.org/schema/stmml";
const char* C_PREFIX           = "cml";

const char* C_2D = "2D";
const char* C_3D = "3D";

vector <string> CML_ELEMENT_VECTOR;
string CML_ELEMENT_NAMES = "amount angle atom atomArray atomParity bond bondArray cml crystal electron formula length molecule stereo substance substanceList symmetry torsion ";

vector <string> STMML_ELEMENT_VECTOR;
string STMML_ELEMENT_NAMES = "action actionList alternative annotation array appinfo definition description dictionary dimension documentation entry enumeration link list matrix metadataList metadata object observation relatedEntry scalar stmml table unit unitList unitType ";

const char* C_ATOMID           = "atomID";
const char* C_ATOMREF          = "atomRef";
const char* C_ATOMREF1        = "atomRef1";
const char* C_ATOMREF2        = "atomRef2";
const char* C_ATOMREFS         = "atomRefs";
const char* C_ATOMREFS1        = "atomRefs1";
const char* C_ATOMREFS2        = "atomRefs2";
const char* C_ATOMREFS3        = "atomRefs3";
const char* C_ATOMREFS4        = "atomRefs4";
const char* C_BUILTIN          = "builtin";
const char* C_CONVENTION       = "convention";
const char* C_CONTENT          = "content";
const char* C_DATATYPE         = "dataType";
const char* C_DICTREF          = "dictRef";
const char* C_FORMALCHARGE     = "formalCharge";
const char* C_ELEMENTTYPE      = "elementType";
const char* C_HYDROGENCOUNT    = "hydrogenCount";
const char* C_ID               = "id";
const char* C_NAME             = "name";
const char* C_OCCUPANCY        = "occupancy";
const char* C_ORDER            = "order";
const char* C_POINTGROUP       = "pointgroup";
const char* C_SPACEGROUP       = "spacegroup";
const char* C_TITLE            = "title";
const char* C_UNITS            = "units";
const char* C_X2               = "x2";
const char* C_Y2               = "y2";
const char* C_X3               = "x3";
const char* C_Y3               = "y3";
const char* C_Z3               = "z3";
const char* C_XFRACT           = "xFract";
const char* C_YFRACT           = "yFract";
const char* C_ZFRACT           = "zFract";
const char* C_XY2              = "xy2";
const char* C_XYZ3             = "xyz3";
const char* C_XYZFRACT         = "xyzFract";

// metadata
const char* DC_DESCRIPTION     = "dc:description";
const char* DC_IDENTIFIER      = "dc:identifier";
const char* DC_CONTENT         = "dc:content";
const char* DC_RIGHTS          = "dc:rights";
const char* DC_TYPE            = "dc:type";
const char* DC_CONTRIBUTOR     = "dc:contributor";
const char* DC_CREATOR         = "dc:creator";
const char* DC_DATE            = "dc:date";

const char* CMLM_STRUCTURE     = "cmlm:structure";

// CML and STMML elements                
const char* C_ANGLE            = "angle";
const char* C_ARRAY            = "array";
const char* C_ATOM             = "atom";
const char* C_ATOMARRAY        = "atomArray";
const char* C_ATOMPARITY       = "atomParity";
const char* C_BOND             = "bond";
const char* C_BONDARRAY        = "bondArray";
const char* C_CML              = "cml";
const char* C_COORDINATE2      = "coordinate2";
const char* C_COORDINATE3      = "coordinate3";
const char* C_CRYSTAL          = "crystal";
const char* C_ELECTRON         = "electron";
const char* C_FEATURE          = "feature";
const char* C_FLOAT            = "float";
const char* C_FLOATARRAY       = "floatArray";
const char* C_FLOATMATRIX      = "floatMatrix";
const char* C_FORMULA          = "formula";
const char* C_INTEGER          = "integer";
const char* C_INTEGERARRAY     = "integerArray";
const char* C_LENGTH           = "length";
const char* C_MATRIX           = "matrix";
const char* C_METADATA         = "metadata";
const char* C_METADATALIST     = "metadataList";
const char* C_MOLECULE         = "molecule";
const char* C_REACTION         = "reaction";
const char* C_SCALAR           = "scalar";
const char* C_SEQUENCE         = "sequence";
const char* C_STEREO           = "stereo";
const char* C_STRING           = "string";
const char* C_STRINGARRAY      = "stringArray";
const char* C_SYMMETRY         = "symmetry";
const char* C_TORSION          = "torsion";


const char* _COLON             = ":";
const char* _EMPTY             = "";
const char* _EQUALS            = "=";
const char* _LANGLE            = "<";
const char* _NEWLINE           = "\n";
const char* _QUERY             = "?";
const char* _QUOTE             = "\"";
const char* _RANGLE            = ">";
const char* _SLASH             = "/";
const char* _SPACE             = " ";
const char* _SPACE_NEWLINE     = " \n";

const char* X_QUOT             = "quot";
const char* X_APOS             = "apos";
const char* X_LT               = "lt";
const char* X_GT               = "gt";
const char* X_AMP              = "amp";
const char* E_TAGO             = "</";
const char* S_XMLDECL          = "<?xml";
const char* E_PI               = "?>";
const char* S_PI               = "<?";
const char* X_DOCTYPE          = "<!DOCTYPE";
const char* X_ENCODING         = "encoding";
const char* X_STANDALONE       = "standalone";
const char* X_SYSTEM           = "SYSTEM";
const char* X_VERSION          = "version";
const char* X_XMLNS            = "xmlns";
const char* X_XML              = "xml";
const char* S_COMMENT          = "<!--";
const char* E_COMMENT          = "-->";
const char* S_CDATA            = "<![CDATA[";
const char* E_CDATA            = "]]>";

const string C_CML1              = "CML1";
const string C_CML2              = "CML2";

// <angle>
bool startAngle(vector <pair<string,string> > &atts);
// </angle>
bool endAngle();
bool WriteAngle(ostream &ofs, pair <vector<OBAtom*>, double> angle);
string ANGLE_ATTRIBUTES =  "id title convention atomRefs atomRefs3 units";
vector <string> ANGLE_ATTRIBUTE_VECTOR;

// <array>
bool startArray(vector <pair<string,string> > &atts);
// </array>
bool endArray();
string ARRAY_ATTRIBUTES =  "id title convention size";
vector <string> ARRAY_ATTRIBUTE_VECTOR;

// <atom>
bool startAtom(vector <pair<string,string> > &atts);
// process builtin children of <atom>
bool processAtomBuiltin();
// </atom>
bool endAtom();
//bool WriteAtom(ostream &ofs, OBAtom* atom, int count);
vector <string> ATOMATTRIBUTE_VECTOR;
string ATOMATTRIBUTES = "id title convention dictRef";
string ATOMBUILTINS =
	"x2 y2 x3 y3 z3 xy2 xyz3 xFract yFract zFract elementType formalCharge hydrogenCount";

// <atomArray>
bool startAtomArray(vector <pair<string,string> > &atts);
// </atomArray>
bool endAtomArray();
//bool WriteAtomArray(ostream &ofs);
vector <string> ATOMARRAY_ATTRIBUTE_VECTOR;
string ATOMARRAY_ATTRIBUTES = "id title convention dictRef";

bool startAtomParity(vector <pair<string,string> > &atts);
bool endAtomParity(vector <pair<string,string> > &atts);

// <bond>
bool startBond(vector <pair<string,string> > &atts);
// process builtin children of <bond>
bool processBondBuiltin();
// </bond>
bool endBond();
//bool WriteBond(ostream &ofs, OBBond* bond);
vector <string> BOND_ATTRIBUTE_VECTOR;
string BOND_ATTRIBUTES = "id title convention dictRef";
string BOND_BUILTINS =  "atomRef atomRefs2 order stereo";

// <bondArray>
bool startBondArray(vector <pair<string,string> > &atts);
// </bondArray>
bool endBondArray();
//bool WriteBondArray(ostream &ofs);
vector <string> BONDARRAY_ATTRIBUTE_VECTOR;
string BONDARRAY_ATTRIBUTES = "id title convention dictRef";

// <cml>
bool startCML(vector <pair<string,string> > &atts);
// </cml>
bool endCML();
string CML_ATTRIBUTES =  "id title convention";
vector <string> CML_ATTRIBUTE_VECTOR;

// <coordinate2>
bool startCoordinate2(vector <pair<string,string> > &atts);
// </coordinate2>
bool endCoordinate2();
string COORDINATE2_ATTRIBUTES =  "id title convention";
vector <string> COORDINATE2_ATTRIBUTE_VECTOR;

// <coordinate3>
bool startCoordinate3(vector <pair<string,string> > &atts);
// </coordinate3>
bool endCoordinate3();
string COORDINATE3_ATTRIBUTES =  "id title convention";
vector <string> COORDINATE3_ATTRIBUTE_VECTOR;

// <crystal>
bool startCrystal(vector <pair<string,string> > &atts);
// </crystal>
bool endCrystal();
bool WriteCrystal(ostream &ofs);
string CRYSTAL_ATTRIBUTES =  "id title convention";
vector <string> CRYSTAL_ATTRIBUTE_VECTOR;

// <electron>
bool startElectron(vector <pair<string,string> > &atts);
// </electron>
bool endElectron();
bool WriteElectron(ostream &ofs);
string ELECTRON_ATTRIBUTES =  "id title convention";
vector <string> ELECTRON_ATTRIBUTE_VECTOR;

// <feature>
bool startFeature(vector <pair<string,string> > &atts);
// </feature>
bool endFeature();
bool WriteFeature(ostream &ofs);
string FEATURE_ATTRIBUTES =  "id title convention";
vector <string> FEATURE_ATTRIBUTE_VECTOR;

// <formula>
bool startFormula(vector <pair<string,string> > &atts);
// </formula>
bool endFormula();
bool WriteFormula(ostream &ofs);
string FORMULA_ATTRIBUTES =  "id title convention";
vector <string> FORMULA_ATTRIBUTE_VECTOR;

// <floatMatrix>
bool startFloatMatrix(vector <pair<string,string> > &atts);
// </floatMatrix>
bool endFloatMatrix();
bool WriteFloatMatrix(ostream &ofs);
string FLOATMATRIX_ATTRIBUTES =  "id title convention";
vector <string> FLOATMATRIX_ATTRIBUTE_VECTOR;

// <length>
bool startLength(vector <pair<string,string> > &atts);
// </length>
bool endLength();
//bool WriteLength(ostream &ofs, pair <vector<OBAtom*>, double> length);
string LENGTH_ATTRIBUTES =  "id title convention atomRefs2 units";
vector <string> LENGTH_ATTRIBUTE_VECTOR;

// <molecule>
bool startMolecule(vector <pair<string,string> > &atts);
// </molecule>
bool endMolecule();
//bool WriteMolecule(ostream &ofs);
vector <string> MOLECULE_ATTRIBUTE_VECTOR;
string MOLECULE_ATTRIBUTES =  "id title count convention";

// <metadataList>
bool WriteMetadataList(ostream &ofs);
    
// <reaction>
bool startReaction(vector <pair<string,string> > &atts);
// </reaction>
bool endReaction();
bool WriteReaction(ostream &ofs);
string REACTION_ATTRIBUTES =  "id title convention";
vector <string> REACTION_ATTRIBUTE_VECTOR;

// <scalar>
bool startScalar(vector <pair<string,string> > &atts);
// </scalar>
bool endScalar();
string SCALAR_ATTRIBUTES =  "id title convention dictRef dataType size units";
vector <string> SCALAR_ATTRIBUTE_VECTOR;

// <sequence>
bool startSequence(vector <pair<string,string> > &atts);
// </sequence>
bool endSequence();
bool WriteSequence(ostream &ofs);
string SEQUENCE_ATTRIBUTES =  "id title convention";
vector <string> SEQUENCE_ATTRIBUTE_VECTOR;

bool startStereo(vector <pair<string,string> > &atts);
bool endStereo(vector <pair<string,string> > &atts);

// <string>
bool startString(vector <pair<string,string> > &atts);
// </string>
bool endString();
vector <string> STRINGATTRIBUTE_VECTOR;
string STRINGATTRIBUTES =  "id title count convention builtin";

bool addString();

// <symmetry>
bool startSymmetry(vector <pair<string,string> > &atts);
// </symmetry>
bool endSymmetry();
bool WriteSymmetry(ostream &ofs, pair <vector<OBAtom*>, double> symmetry);
string SYMMETRY_ATTRIBUTES =  "id title convention atomRefs atomRefs4 units";
vector <string> SYMMETRY_ATTRIBUTE_VECTOR;

// <torsion>
bool startTorsion(vector <pair<string,string> > &atts);
// </torsion>
bool endTorsion();
bool WriteTorsion(ostream &ofs, pair <vector<OBAtom*>, double> torsion);
string TORSION_ATTRIBUTES =  "id title convention atomRefs atomRefs4 units";
vector <string> TORSION_ATTRIBUTE_VECTOR;

// normalize PCDATA for builtins
void processBuiltinPCDATA();
// <fooArray> children of atomArray
bool processAtomArrayChild();
// <fooArray> children of bondArray
bool processBondArrayChild();

// ----------------CML utilities---------------
// translates attributes (e.g. atomRefs3) into atomRef vectors
void getAtomRefs(vector<string>::size_type size, vector <OBAtom*> &v, string atomRefString);

// ----------------babel stuff---------------
// clear Molecule workspace every time we need a new one
bool clearMoleculeWorkspace();
// numeric value (babel) from CML order; unknown returns -1
int getBabelBondOrder(string o);
// numeric value (babel) from CML stereo; unknown returns -1
int getBabelBondFlag(string s);
// get serial number of atom (starting at 1); 0 for not found
OBAtom *getAtomPtr(string s);

// debug
void debug(ostream &ofs);

// ----------------variables-----------------
  // \todo Get rid of all the global variables--potential memory leaks.

// variables
// output stream
ostream *ofsPtr;
// tracks level of XML elements
vector <string> elementStack;
// current element name
string parent;
// current parent element name
string currentElem;
// last PCDATA (text) content read
string pcdata;
// current attributes
vector<pair<string,string> > currentAtts;
// namespaces
typedef vector<pair<string,string> > namespaceVector_t;
namespaceVector_t namespaceVector;
// has rootElement been read?
bool readRoot;
// are we in a comment?
bool inComment = false;

// all atoms need to be unique IDs in CML; this links IDs to AtomPtrs
vector <pair <string, OBAtom*> > atomIdVector;
// current Molecule
OBMol *molPtr;
// current Atom
OBAtom *atomPtr;
// current Bond
OBBond *bondPtr;
// dimensionality of coordinates
const char *dimension;
// only CML can separate 2 and 3D
string cmlDimension;

// atom stuff
int natoms;
int atomicNum;
string atomId;
int formalCharge;	// defaults to zero
// coordinates
double currentX, currentY, currentZ;
string elementArray;
string chargeArray;
string idArray;
string x2Array;
string y2Array;
string x3Array;
string y3Array;
string z3Array;

string atomRefs4;

double length;

vector <string> idVector;
vector <string> elementTypeVector;
vector <int> atomicNumVector;
vector <int> formalChargeVector;
vector <int> hydrogenCountVector;
vector <double> x2Vector;
vector <double> y2Vector;
vector <double> x3Vector;
vector <double> y3Vector;
vector <double> z3Vector;

// current reference to lists of atoms
vector <OBAtom*> atomRefs2Vector;
vector <OBAtom*> atomRefs3Vector;
vector <OBAtom*> atomRefs4Vector;

// bond stuff
unsigned int nbonds;	// [ejk] assumed unsigned is OK, where is the global guy set?
// ends of bond
string bondBeginAtom;
string bondEndAtom;
string orderString;
string stereoString;

string atomRef1Array;
string atomRef2Array;
string orderArray;
string stereoArray;

vector <string> atomRef1Vector;
vector <string> atomRef2Vector;
vector <string> orderVector;
vector <string> stereoVector;

// using builtins?
bool useBuiltin;

bool inputNamespace;
char *inputPrefix;
bool inputArray;

string cmlType = "";
bool outputCML1;
bool outputCML2;
bool outputDoctype;
bool outputDeclaration;
bool outputPretty;
bool outputNamespace;
string outputPrefix;
bool outputArray;
bool outputDebug;

string angleUnits;
string lengthUnits;
string torsionUnits;

// scalar
string scalarDataType;
string scalarUnits;

// crystal
bool fractional;
double cellParam[6];
// symmetry
string spacegroup;
string pointgroup;
vector <double*> rotTransVector;
vector <double*> rotVector;

// angle
vector <pair <vector<OBAtom*>, double> > angleVector;
// length
vector <pair <vector<OBAtom*>, double> > lengthVector;
// torsion
vector <pair <vector<OBAtom*>, double> > torsionVector;
// atomParity
vector <pair <vector<OBAtom*>, double> > atomParityVector;
// stereo
vector <pair <vector<OBAtom*>, string> > stereoSVector;

// internal coordinates
typedef vector <OBInternalCoord*> internalVector_t;
internalVector_t internalVector;

/** --------------------- sizes -----------------------*/
// this should be replaced by a better routine for reading
#define XMLBUFF_SIZE 100000

/** --------------------- initialization-----------------------*/
void makeAllowedElementLists() {
	tokenize(CML_ELEMENT_VECTOR, CML_ELEMENT_NAMES, " \n");
	tokenize(STMML_ELEMENT_VECTOR, STMML_ELEMENT_NAMES, " \n");
}

void makeAllowedAttributeLists() {
	tokenize(ANGLE_ATTRIBUTE_VECTOR, ANGLE_ATTRIBUTES, " \n");
	string s = ATOMATTRIBUTES;
	s.append(_SPACE);
	s.append(ATOMBUILTINS);
	tokenize(ATOMATTRIBUTE_VECTOR, s, " \n");
	tokenize(ATOMARRAY_ATTRIBUTE_VECTOR, ATOMARRAY_ATTRIBUTES, " \n");
	s = BOND_ATTRIBUTES;
	s.append(_SPACE);
	s.append(BOND_BUILTINS);
	tokenize(BOND_ATTRIBUTE_VECTOR, s, " \n");
	tokenize(BONDARRAY_ATTRIBUTE_VECTOR, BONDARRAY_ATTRIBUTES, " \n");
	tokenize(CML_ATTRIBUTE_VECTOR, CML_ATTRIBUTES, " \n");
	tokenize(COORDINATE2_ATTRIBUTE_VECTOR, COORDINATE2_ATTRIBUTES, " \n");
	tokenize(COORDINATE3_ATTRIBUTE_VECTOR, COORDINATE3_ATTRIBUTES, " \n");
	tokenize(CRYSTAL_ATTRIBUTE_VECTOR, CRYSTAL_ATTRIBUTES, " \n");
	tokenize(ELECTRON_ATTRIBUTE_VECTOR, ELECTRON_ATTRIBUTES, " \n");
	tokenize(FEATURE_ATTRIBUTE_VECTOR, FEATURE_ATTRIBUTES, " \n");
	tokenize(FORMULA_ATTRIBUTE_VECTOR, FORMULA_ATTRIBUTES, " \n");
	tokenize(LENGTH_ATTRIBUTE_VECTOR, LENGTH_ATTRIBUTES, " \n");
	tokenize(MOLECULE_ATTRIBUTE_VECTOR, MOLECULE_ATTRIBUTES, " \n");
	tokenize(REACTION_ATTRIBUTE_VECTOR, REACTION_ATTRIBUTES, " \n");
	tokenize(SCALAR_ATTRIBUTE_VECTOR, SCALAR_ATTRIBUTES, " \n");
	tokenize(SEQUENCE_ATTRIBUTE_VECTOR, SEQUENCE_ATTRIBUTES, " \n");
	tokenize(SYMMETRY_ATTRIBUTE_VECTOR, SYMMETRY_ATTRIBUTES, " \n");
	tokenize(TORSION_ATTRIBUTE_VECTOR, TORSION_ATTRIBUTES, " \n");
}

/** ------------------------------XML reader--------------------------*/

//void cmlError(int type, string msg, int line) {
//}

void cmlError(string msg) {
     cout << msg << endl;
}

bool ReadXML(istream &ifs) {
//	char buffer[XMLBUFF_SIZE];
	char buffer[BUFF_SIZE];
    string token;
	size_t lt;
	size_t rt;
    int lineCount = 0;

	bool lookForOpenTag = true;

	makeAllowedElementLists();
	makeAllowedAttributeLists();
// reset all variables here (ReadXML is essentially a static routine);    
	startDocument();
//char *fgets(char *s, size_t n, FILE *stream); 
/*--
	while (ifs.getline(buffer,XMLBUFF_SIZE)) {
        if (strlen(buffer) >= XMLBUFF_SIZE - 2) {
            cmlError("******BUG. Probable buffer overflow - sorry - add some newlines to XML input file to reduce linesize. ******");
        }
--*/        
    string buff;
	while (getline(ifs, buff)) {
        lineCount++;
//		string buff(buffer);
// omit whitespace lines
		if (trim(buff) == _EMPTY) continue;
		if (readRoot) {
		  	// cmlError(FATAL, "no nonWhitespace allowed after root element: " + buff;
			return false;
		}
// normalize Newlines to SPACE
		if (token != _EMPTY) token += _SPACE;
		for (;;) {
            if (inComment) {
                lt = buff.find(E_COMMENT);
                if (lt > buff.size()) {
					token += buff;
					buff = _EMPTY;
					break;
                } else {
                    inComment = false;
                    buff = buff.substr(lt+1);
                    continue;
                }
            }
			if (lookForOpenTag) {
				lt = buff.find("<");
// not found, more input...
				if (lt > buff.size()) {
					token += buff;
					buff = _EMPTY;
					break;
				} else {
// found start of tag
					token += buff.substr(0,lt);
// process inter-tag characters
					characters(token);
					buff = buff.substr(lt);
					token = _EMPTY;
					lookForOpenTag = false;
				}
			} else {
				rt = buff.find(">");
// not found, more input...
				if (rt > buff.size()) {
					token += buff;
					buff = _EMPTY;
					break;
				} else {
// found end
					string ss = buff.substr(0,rt+1);
					token += ss;
					tag(token);
					buff = buff.substr(rt+1);
					token = _EMPTY;
					lookForOpenTag = true;
				}
			}
			if (buff == _EMPTY) break;
		}
	}
	endDocument();
	return true;	// [ejk] just guessed at this return value
}

// process anything in balanced <...>
void tag(string s) {
	// attributes
	vector <pair <string, string> > atts;

	string name;
	string::size_type l = s.length();
	string sl = toLowerCase(s);
// XML declaration
	if (sl.substr(0, 5) == S_XMLDECL) {
		if (s.substr(l-2, 2) == E_PI) {
			string ss = s.substr(5, l-7);
			splitAttributes(ss, atts);
			string standalone = getAttribute(atts, X_STANDALONE);
			if (standalone == "no") {
			  cmlError("cannot process standalone='no' yet");
			}
			string version = getAttribute(atts, X_VERSION);
			if (version != "1.0") {
			  cmlError("XML version must be 1.0");
			}
			string encoding = toLowerCase(getAttribute(atts, X_ENCODING));
			if (encoding != "utf-8" && encoding != _EMPTY) {
			  cmlError("Cannot support encoding: " + encoding);
			}
		} else {
		  cmlError("Bad XML declaration: " + s);
		}
// DOCTYPE is not processed
	} else if (s.substr(0,9) == X_DOCTYPE) {
		if (s.find("[") <= s.size()) {
		  			cmlError("cannot process internal subset of DOCTYPE " + s);
		} else {
		  //			cout << "DOCTYPE info ignored" << endl;
		}
// comments are ignored
	} else if (s.substr(0,4) == S_COMMENT) {
		if (s.substr(l-3, 3) == E_COMMENT) {
            inComment = false;
//			cout << "Comment ignored: " << s << endl;
		} else {
            inComment = true;
		  cmlError("Bad comment: " + s);
		}
// Processing instructions
	} else if (s.substr(0,2) == S_PI) {
		if (s.substr(l-2, 2) == E_PI) {
			s = s.substr(2, l-4);
			string::size_type idx = s.find(_SPACE);
			string target = (idx < s.size()) ? s.substr(0, idx) : s;
			string data = (idx < s.size()) ? trim(s.substr(idx)) : string(_EMPTY);
			processingInstruction(target, data);
		} else {
		  cmlError("Bad PI: " + s);
		}
// CDATA sections
	} else if (s.substr(0,9) == S_CDATA) {
		if (s.substr(l-3, 3) == E_CDATA) {
			pcdata += s.substr(9, l-12);
		} else {
		  cmlError("Bad CDATA: " + s);
		}
// end tag
	} else if (s.substr(1,1) == _SLASH) {
		endElement(s.substr(2, l-3));
// empty tag
	} else if (s.substr(l-2, 1) == _SLASH) {
		name = startTag(s.substr(1,l-3));
		endElement(name);
// start tag
	} else {
		startTag(s.substr(1, l-2));
	}
}

string escapeXMLEntities(string s) {
	string ss;
	int ii;
	char *cc = (char*) s.c_str();
	for (string::size_type i = 0; i < s.length(); ++i) {
		ii = (int) cc[i];
		if (cc[i] == '&') {
			ss.append("&");
            ss.append(X_AMP);
            ss.append(";");
		} else if (cc[i] == '"') {
			ss.append("&");
            ss.append(X_QUOT);
            ss.append(";");
		} else if (cc[i] == '\'') {
			ss.append("&");
            ss.append(X_APOS);
            ss.append(";");
		} else if (cc[i] == '<') {
			ss.append("&");
            ss.append(X_LT);
            ss.append(";");
		} else if (cc[i] == '>') {
			ss.append("&");
            ss.append(X_GT);
            ss.append(";");
// characters above 255
		} else if (ii > 255) {
		  cmlError("characters above 255 not supported in CML: " + ii);
// characters 128-255
		} else if (ii > 127) {
			ss.append("&#");
			ss.append(""+ii);
			ss.append(";");
// characters >= 32
		} else if (cc[i] > ' ') {
			ss.append(1, cc[i]);
// white space
		} else if (cc[i] == ' ' || cc[i] == '\t' || cc[i] == '\n' || cc[i] == '\r') {
			ss.append(1, cc[i]);
		} else {
		  cmlError("non-printing characters not suported: " + (int)cc[i]);
		}
	}
	return ss;
}

string processXMLEntities(string s) {
	string s0(s);
	string ss;
	for (;;) {
		string::size_type idx = s.find("&");
		if (idx >= s.length()) {
			ss.append(s);
			break;
		}
		ss.append(s.substr(0, idx));
		s = s.substr(idx+1);
		idx = s.find(";");
		if (idx >= s.length()) {
		  cmlError("entity without closing ; in :" + s0 + _COLON);
		}
		string e = s.substr(0, idx);
		if (e == X_QUOT) {
			ss.append("\"");
		} else if (e == X_APOS) {
			ss.append("'");
		} else if (e == X_LT) {
			ss.append("<");
		} else if (e == X_GT) {
			ss.append(">");
		} else if (e == X_AMP) {
			ss.append("&");
		} else if (e.substr(0, 1) == "#") {
			int i = atoi((char*)e.substr(1).c_str());
			if (i >= 32 && i < 256 || i == 9 || i==10 || i==13) {
				ss.append(1, (char)i);
			} else {
			  cmlError("unsupported character: #" + i);
			}
		} else {
			skippedEntity(e);
		}
		s = s.substr(idx+1);
	}
	return ss;
}

string startTag(string s) {
	vector <pair<string,string> > atts;

	s = trim(s);
	if (s.find("&") <= s.size()) {
	  cmlError("CML reader cannot process entity references (sorry)..." + s);
	}
	string ss = s;
	string name;
	string::size_type idx = s.find(_SPACE);
	if (idx > s.size()) {
		name = s;
		s = _EMPTY;
	} else {
		name = s.substr(0, idx);
		s = trim(s.substr(idx+1));
	}
	splitAttributes(s, atts);
	if (!isXMLName(name)) {
	  cmlError("invalid XML name: " + name);
	}
	startElement(name, atts);
	return name;
}

void splitAttributes(string s, vector <pair <string, string> > &atts) {
	pair<string, string> att;

	while (true) {
		string::size_type idx = s.find("=");
		if (idx > s.size()) {
			if (trim(s) != _EMPTY) {
			  cmlError("Bad attribute at " + s);
			}
			break;
		}
		att.first = trim(s.substr(0, idx));
		s = trim(s.substr(idx+1));
		if (s.length() < 2) {
		  cmlError("Bad attribute value: " + s);
			break;
		}
// quote or X_APOS
		string quoter = s.substr(0, 1);
		if (quoter != "\"" && quoter != "\'") {
		  cmlError("Unquoted attribute value: " + s);
			break;
		}
		s = s.substr(1);
		idx = s.find(quoter);
		if (idx > s.size()) {
		  cmlError("Unbalanced quotes in attribute value: " + s);
			break;
		}
		att.second = processXMLEntities(s.substr(0, idx));
		atts.push_back(att);
		s = trim(s.substr(idx+1));
		if (trim(s) == _EMPTY) break;
	}
}

// check attributes against allowed list; result is unknown attributes
vector <string> getUnknownAttributes(vector <string> &allowed, vector <pair <string, string> > &atts) {
	vector <string> badAtts;
	for (vector<pair <string,string> >::size_type i = 0; i < atts.size(); ++i) {
		string attName = atts[i].first;
		if (attName.substr(0, 5) == X_XMLNS) continue;
		bool ok = false;
		for (vector<string>::size_type j = 0; j < allowed.size(); ++j) {
			if (allowed[j] == attName) {
				ok = true;
				break;
			}
		}
		if (!ok) {
			badAtts.push_back(attName);
		}
	}
	return badAtts;
}

void printVector(vector <string> v, ostream& ofs) {
	for (vector<string>::size_type i = 0; i < v.size(); ++i) {
		ofs << v[i] << _SPACE;
	}
}

void noteUnusedElementName(string name, string msg) {
	if (!isInStringVector(unusedElementNameVector, name)) {
	  //		cout << msg << name << endl;
		unusedElementNameVector.push_back(name);
	}
}


bool writeAttribute(ostream &ofs, string name, int value) {
	ofs << _SPACE << name << _EQUALS << _QUOTE << value << _QUOTE;
	return true; // [ejk] assumed
}

bool writeAttribute(ostream &ofs, string name, double value) {
	ofs << _SPACE << name << _EQUALS << _QUOTE << value << _QUOTE;
	return true; // [ejk] assumed
}

bool writeAttribute(ostream &ofs, string name, string value) {
	value = trim(value);
	if (value != _EMPTY) {
		string value1 = escapeXMLEntities(value);
        ofs << _SPACE << name << _EQUALS << _QUOTE << value1 << _QUOTE;
	}
	return true; // [ejk] assumed
}

bool writeStartTagStart(ostream& ofs, string name) {
    ofs << _LANGLE << outputPrefix << name;
}

bool writeStartTagEnd(ostream& ofs) {
    ofs << _RANGLE;
}

bool writeEndTag(ostream& ofs, string name) {
    ofs << _LANGLE << _SLASH << outputPrefix << name << _RANGLE << endl;
}

bool writeCombinedTagEnd(ostream& ofs) {
    ofs << _SLASH << _RANGLE << endl;
}

bool writeBuiltin(ostream&ofs, string name, int value) {
	ofs << "<" << outputPrefix << "integer builtin=\"" <<  name << "\">" << value << E_TAGO << outputPrefix << "integer>" << endl;
	return true; // [ejk] assumed success
}

bool writeBuiltin(ostream&ofs, string name, double value) {
	ofs << "<" << outputPrefix << "float builtin=\"" <<  name << "\">" << value << E_TAGO << outputPrefix << "float>" << endl;
	return true; // [ejk] assumed success
}

bool writeBuiltin(ostream&ofs, string name, string value) {
	value = trim(value);
	if (value != _EMPTY) {
		value = escapeXMLEntities(value);
		ofs << "<" << outputPrefix << "string builtin=\"" <<  name << "\">" << value << E_TAGO << outputPrefix << "string>" << endl;
	}
	return true; // [ejk] assumed
}

bool appendToArray(string &array, int value) {
	if (array != _EMPTY) array.append(_SPACE);
	char ss[20];
	sprintf(ss, "%i", value);
	string s(ss);
	array.append(trim(ss));
	return true; // [ejk] assumed
}

bool appendToArray(string &array, double value) {
	if (array != _EMPTY) array.append(_SPACE);
	char ss[20];
	sprintf(ss, "%f", value);
	string s(ss);
	array.append(trim(ss));
	return true; // [ejk] assumed
}

bool appendToArray(string &array, string value) {
	value = escapeXMLEntities(value);
	if (array != _EMPTY) array.append(_SPACE);
	array.append(trim(value));
	return true; // [ejk] assumed
}

bool writePCDATA(ostream&ofs, string value) {
	ofs << escapeXMLEntities(value);
	return true; // [ejk] assumed
}

// normalizes string
string getNormalizedString(string s) {
    return getNormalizedString((char*)s.c_str());
}

// normalizes string
string getNormalizedString(char* ch) {
    bool inWhite = true;
    bool start = true;
    string s = "";
    for (int i = 0;; i++) {
        char c = ch[i];
        if (c == 0) break;
        if (c == ' ' || c == '\n' || c == '\t' || c == '\r') {
            inWhite = true;
        } else {
            if (start) {
                start = false;
            } else {
                if (inWhite) s += " ";
            }
            inWhite = false;
            s += c;
        }
    }
    return s;
}

// ------------------------ SAX events -------------------
// SAX-like call back
void startDocument() {
//  cout << "starting CML document; crude XML parser. Assumes well-formed; ignores DTDs and entities" << endl;
    readRoot = false;
// clear all internals
	currentElem = _EMPTY;
	string token = _EMPTY;
	inComment = false;
    cmlDimension = "";

    clearMoleculeWorkspace();

    useBuiltin = false;
    inputNamespace = "";
    inputPrefix = "";
    inputArray = false;
    cmlType = "";
    outputCML1 = false;
    outputCML2= false;
    outputDoctype = "";
    outputDeclaration = false;
    outputPretty= false;
    outputNamespace= false;
    outputPrefix = "";
    outputArray= false;
    outputDebug = false;
    
    angleUnits = "";
    lengthUnits = "";
    torsionUnits = "";
    scalarDataType = "";
    scalarUnits = "";
  
}


// SAX-like call back
void endDocument() {
  //	cout << "read CML document" << endl;
	for (namespaceVector_t::size_type i = 0; i < namespaceVector.size(); ++i) {
	  //		cout << "namespace :" << namespaceVector[i].first << _COLON << namespaceVector[i].second << endl;
	}
}

bool clearMoleculeWorkspace() {    
    natoms = 0;
    atomicNum = 0;
    atomId = "";
    formalCharge = 0;	// defaults to zero
    currentX = 0.0; 
    currentY = 0.0;
    currentZ = 0.0;
    elementArray = "";
    chargeArray = "";
    idArray = "";
    x2Array = "";
    y2Array = "";
    x3Array = "";
    y3Array = "";
    z3Array = "";
    atomRefs4 = "";
    length = 0.0;
    idVector.clear();
    elementTypeVector.clear();
    atomicNumVector.clear();
    formalChargeVector.clear();
    hydrogenCountVector.clear();
    x2Vector.clear();
    y2Vector.clear();
    x3Vector.clear();
    y3Vector.clear();
    z3Vector.clear();
    atomRefs2Vector.clear();
    atomRefs3Vector.clear();
    atomRefs4Vector.clear();
    nbonds = 0;	
    bondBeginAtom = "";
    bondEndAtom = "";
    orderString = "";
    stereoString = "";
    atomRef1Array = "";
    atomRef2Array = "";
    orderArray = "";
    stereoArray = "";
    atomRef1Vector.clear();
    atomRef2Vector.clear();
    orderVector.clear();
    stereoVector.clear();
    fractional = false;
    spacegroup = "";
    pointgroup = "";
    rotTransVector.clear();
    rotVector.clear();
    angleVector.clear();
    lengthVector.clear();
    torsionVector.clear();
    atomParityVector.clear();
    stereoSVector.clear();
}

void startElement(string name, vector<pair<string,string> > &atts) {
	processAttributes(atts);
	pair <string, string> nsPair = getNamespacePair(name);
	name = (nsPair.first == _EMPTY) ? name : name.substr(nsPair.first.length()+1);
	startElement(nsPair.second, name, nsPair.first, atts);
}

pair <string, string> getNamespacePair(string name) {
	pair <string, string> nsPair;
	nsPair.first = _EMPTY;
	nsPair.second = _EMPTY;
	string::size_type idx = name.find(_COLON);
	if (idx < name.length()) {
		nsPair.first = name.substr(0, idx);
		name = name.substr(idx+1);
	}
	for (namespaceVector_t::size_type i = 0; i < namespaceVector.size(); ++i) {
		if (namespaceVector[i].first == nsPair.first) {
			nsPair.second = namespaceVector[i].second;
			break;
		}
	}
	return nsPair;
}

// records whether this is CMLV1.0 or CML2 (Schema)
void setCMLType(string ct) {
	if (cmlType == _EMPTY) {
		cmlType = ct;
	} else if (cmlType != ct) {
  	    cmlError("Cannot mix namespaces, was: " + cmlType);
	}
}

void startElement(string namespaceURI, string localName, string prefix, vector<pair<string,string> > &atts) {

   
	if (currentElem != _EMPTY) elementStack.push_back(currentElem);
	parent = currentElem;
	currentElem = trim(localName);
	currentAtts = atts;
	pcdata = _EMPTY;
	useBuiltin = false;
    if (localName == C_ATOM) {
        startAtom(atts);
	} else if (localName == C_ATOMARRAY) {
        startAtomArray(atts);
	} else if (localName == C_ATOMPARITY) {
		setCMLType(C_CML2);
        startAtomParity(atts);
	} else if (localName == C_BOND) {
        startBond(atts);
	} else if (localName == C_BONDARRAY) {
        startBondArray(atts);
    } else if (localName == C_CML) {
		startCML(atts);
    } else if (localName == C_CRYSTAL) {
		startCrystal(atts);
    } else if (localName == C_ELECTRON) {
		startElectron(atts);
    } else if (localName == C_FEATURE) {
		startFeature(atts);
		setCMLType(C_CML1);
    } else if (localName == C_FORMULA) {
		startFormula(atts);
    } else if (localName == C_MOLECULE) {
		startMolecule(atts);
	} else if (
		localName == C_COORDINATE2 ||
		localName == C_COORDINATE3 ||
		localName == C_FLOAT ||
		localName == C_INTEGER ||
	    localName == C_STRING) {
			setCMLType(C_CML1);
	// delay processing till endElement as we need the pcdata
	} else if (localName == C_FLOATMATRIX) {
		setCMLType(C_CML1);
	} else if (localName == C_FLOATARRAY) {
//        startFloatArray(atts)
		setCMLType(C_CML1);
	} else if (localName == C_INTEGERARRAY) {
//        startIntegerArray(atts)
		setCMLType(C_CML1);
	} else if (localName == C_STRINGARRAY) {
//        startStringArray(atts)
		setCMLType(C_CML1);
	// delay processing till endElement as we need the pcdata
	} else if (localName == C_LENGTH) {
		setCMLType(C_CML2);
		startLength(atts);
	} else if (localName == C_ANGLE) {
		startAngle(atts);
	} else if (localName == C_TORSION) {
		startTorsion(atts);
	} else if (localName == C_SCALAR) {
		startScalar(atts);
		setCMLType(C_CML2);
	} else if (localName == C_STEREO) {
		startStereo(atts);
		setCMLType(C_CML2);
	} else if (localName == C_ARRAY) {
		setCMLType(C_CML2);
	} else if (localName == C_MATRIX) {
		setCMLType(C_CML2);
// other CML2 elements neglected in babel
	} else if (
		localName == "substance" ||
		localName == "substanceList" ||
		localName == "amount") {
		setCMLType(C_CML2);
		noteUnusedElementName(localName, "CML2 element not relevant to babel: ");
	} else if (isInStringVector(CML_ELEMENT_VECTOR, localName)) {
	  //		cout << "[debug] CML element not relevant to babel: " << localName << endl;
// other STMML elements - neglected in Babel
	} else if (
		localName == "actionList" ||
		localName == "action" ||
		localName == "alternative" ||
		localName == "annotation" ||
		localName == "appinfo" ||
		localName == "definition" ||
		localName == "description" ||
		localName == "dictionary" ||
		localName == "dimension" ||
		localName == "documentation" ||
		localName == "entry" ||
		localName == "enumeration" ||
		localName == "list" ||
		localName == "link" ||
		localName == "metadata" ||
		localName == "metadataList" ||
		localName == "object" ||
		localName == "observation" ||
		localName == "relatedEntry" ||
		localName == "table" ||
		localName == "unit" ||
		localName == "unitType" ||
		localName == "unitList") {
		setCMLType(C_CML2);
		noteUnusedElementName(localName, "STMML element not relevant to babel: ");
	} else if (isInStringVector(STMML_ELEMENT_VECTOR, localName)) {
		setCMLType(C_CML2);
		//		cout << "[debug] STMML element not relevant to babel: " << localName << endl;
	} else {
	  //		cout << "start element ignored: " << localName << endl;
	}
}

void processAttributes(vector<pair<string,string> > &atts) {
	for (vector<pair<string,string> >::size_type i = 0; i < atts.size(); ++i) {
		string name = atts[i].first;
		if (!isXMLName(name)) {
		  cmlError("invalid XML name: " + name);
		} else if (name.substr(0,5) == X_XMLNS) {
			processNamespace(name.substr(5), atts[i].second);
		}
	}
}

void processNamespace(string name, string value) {
	pair <string, string> ns;

	string::size_type idx = name.find(_COLON);
	ns.first = (idx < name.size()) ? name.substr(idx) : string(_EMPTY);
	ns.second = value;
	bool nsExists = false;
	for (namespaceVector_t::size_type i = 0; i < namespaceVector.size(); ++i) {
		if (ns.first == namespaceVector[i].first) {
			nsExists = true;
			if (namespaceVector[i].second != value) {
			  cmlError("redefinition of namespace: " +
			      namespaceVector[i].second + " => " + value);
			}
			break;
		}
	}
	if (!nsExists) {
		namespaceVector.push_back(ns);
		if (ns.second == STMML_NAMESPACE) {
			setCMLType(C_CML2);
		} else if (ns.second == CML2_NAMESPACE) {
			setCMLType(C_CML2);
		} else if (ns.second == CML1_NAMESPACE) {
		}
	}
}

void endElement(string name) {
	pair <string, string> nsPair = getNamespacePair(name);
	name = (nsPair.first == _EMPTY) ? name : name.substr(nsPair.first.length()+1);
	endElement(nsPair.second, name, nsPair.first);
}

void endElement(string namespaceURI, string localName, string prefix) {
	vector <string> strings;

	string name = trim(localName);
	if (name != currentElem) {
	  cmlError("unbalanced tags at: " + name);
	}
    if (name == C_MOLECULE) {
		endMolecule();
	} else if (name == C_ATOM) {
		endAtom();
	} else if (name == C_ATOMARRAY) {
		endAtomArray();
	} else if (name == C_ATOMPARITY) {
		//endAtomArray();
	} else if (name == C_BOND) {
		endBond();
	} else if (name == C_BONDARRAY) {
		endBondArray();
	} else if (name == C_CRYSTAL) {
		endCrystal();
	} else if (name == C_ELECTRON) {
		endElectron();
	} else if (name == C_FORMULA) {
		endFormula();
	} else if (name == C_FEATURE) {
//		endFormula();
	} else if (
		name == C_COORDINATE2 ||
		name == C_COORDINATE3) {
		if (C_ATOM == parent) {
			processAtomBuiltin();
		} else {
		  //			cout << "IGNORED <" << name << "> as not child of <atom>" << endl;
		}
	} else if (
		name == C_FLOAT ||
		name == C_INTEGER ||
	    name == C_STRING) {
		if (C_ATOM == parent) {
			processAtomBuiltin();
		} else if (C_BOND == parent) {
			processBondBuiltin();
		} else if (C_MOLECULE == parent) {
			addString();
		} else {
		  //			cout << "IGNORED <" << name << "> as not child of <molecule>, <atom> or <bond>" << endl;
		}
	} else if (localName == C_ARRAY) {
		setCMLType(C_CML2);
	} else if (localName == C_MATRIX) {
		setCMLType(C_CML2);
	} else if (
		localName == C_FLOATMATRIX) {
		setCMLType(C_CML1);
		// no action
	} else if (
        localName == C_FLOATARRAY ||
		localName == C_INTEGERARRAY ||
	    localName == C_STRINGARRAY) {
        setCMLType(C_CML1);
        inputArray = true;
		if (parent == C_ATOMARRAY) {
			processAtomArrayChild();
		} else if (parent == C_BONDARRAY) {
			processBondArrayChild();
		}
	} else if (name == C_LENGTH) {
		endLength();
	} else if (name == C_ARRAY) {
		endAngle();
	} else if (name == C_TORSION) {
		endTorsion();
	} else if (name == C_REACTION) {
		endReaction();
	} else if (name == C_SEQUENCE) {
		endSequence();
	} else {
//		cout << "end element ignored: " << name << endl;
	}
// I assume there is a vector<> function which is neater...
    int ns = elementStack.size();
    if (ns > 0) {
	    currentElem = elementStack[ns-1];
	    parent = (ns <= 1) ? string(_EMPTY) : elementStack[ns-2];
	    elementStack.pop_back();
	} else {
	}
	if (ns == 0) {
		readRoot = true;
	}
	pcdata = _EMPTY;
}

void characters(string s) {
	pcdata = processXMLEntities(s);
}

void processingInstruction(string target, string data) {
  //	cout << "PI: " << target << _SPACE << data << endl;
}

void skippedEntity(string name) {
  //	cout << "skipped entity: " << name << endl;
}

// gets attribute of given name; "" if not found
string getAttribute(vector <pair<string, string> > &atts, string name) {
	string s;
    for (vector<pair<string,string> >::size_type i = 0; i < atts.size(); ++i) {
		if (atts[i].first == name) {
			return atts[i].second;
		}
	}
	return _EMPTY;
}

bool isXMLName(string n) {
	bool ok = true;;

	char *str = (char*) n.c_str();
	char c = *str++;
// first character must be a-zA-Z_
	if (!( (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || c == '_' )) {
		ok = false;
	}
	while ((c = *str++) != 0) {
		if (c >= '0' && c <= '9') {
		} else if (c >= 'a' && c <= 'z') {
		} else if (c >= 'A' && c <= 'Z') {
		} else if (c == '_' || c == ':' || c == '-' || c == '.') {
		} else {
			ok = false;
		}
	}
	if (!ok) {
	  cmlError("invalid XML name: " + n);
	}
	return ok;
}

// -------------------------- utilities --------------------
// trims whitespace from start and end.
string trim(string s) {
	char *ss = (char*) s.c_str();
	int l;

	l = strlen(ss);
	while ((l > 0) && (ss[0] == ' ' || ss[0] == '\n' || ss[0] == '\t')) {
		ss++;
	    --l;
	}

	while ((l > 0) && (ss[l-1] == ' ' || ss[l-1] == '\n' || ss[l-1] == '\t')) {
		ss[l-1] = '\0';
	    --l;
	}
	string sss(ss);
	return sss;
}

string toLowerCase(string s) {
	string ss(s);
	unsigned int i;
	for (i = 0;i < s.size();++i) {
	    ss[i] = tolower(s[i]);
	}
	return ss;
}

string toUpperCase(string s) {
	string ss(s);
	unsigned int i;
	for (i = 0;i < s.size();++i) {
	    ss[i] = toupper(s[i]);
	}
	return ss;
}

bool isInStringVector(vector <string> v, string s) {
	for (vector<string>::size_type i = 0; i < v.size(); ++i) {
		if (v[i] == s) {
			return true;
		}
	}
	return false;
}

// ---------------------- CML utilties -----------------

// utility for all PCATA in <foo builtin="*"> (CML1)
void processBuiltinPCDATA() {
	vector <string> strings;
	useBuiltin = true;
// normalize whitespace
	tokenize(strings, pcdata, " \t\n");
	if (strings.size() != 1) {
	  cmlError("must give value for builtin");
		pcdata = _EMPTY;
	} else {
		pcdata = strings[0];
	}
}

// get pointer to atom
OBAtom *getAtomPtr(string s) {
// this is crude... I expect vector<> has a map utility
	for (string::size_type j = 0; j < atomIdVector.size(); ++j) {
		if (s == atomIdVector[j].first) {
			return atomIdVector[j].second;
		}
	}
	return 0;
}

// process atomRefs2="a1 a2", etc.
// add results to vector
void getAtomRefs(vector<string>::size_type size, vector <OBAtom*> &v, string atomRefString) {
	vector <string> sv;
	atomRefString += _SPACE;
	tokenize(sv, atomRefString, " \n");
	if (sv.size() != size) {
//	  cmlError("unexpected size for atomRefs attribute: " + sv.size() + _SLASH + size);
		return;
	}
	for (vector<string>::size_type i = 0; i < size; ++i) {
		OBAtom* atPtr = getAtomPtr(sv[i]);
		if (atPtr == 0) {
		  cmlError("cannot find atom: " + sv[i]);
			return;
		}
		v.push_back(atPtr);
	}
}

// ---------------------- CML and STMML elements -----------------

// ------------------- <cml> --------------------
bool startCML(vector <pair<string,string> > &atts) {
	vector <string> badAtts = getUnknownAttributes(CML_ATTRIBUTE_VECTOR, atts);
	if (badAtts.size() > 0) {
	  cmlError("unknown attributes on <cml>: ");
		printVector(badAtts, cerr);
	}
	return true; // [ejk] assumed
}

bool endCML() {
	return true; // [ejk] assumed
}
// --------------------<angle>/<angleArray>----------------

bool startAngle(vector <pair<string,string> > &atts) {
	vector <string> badAtts = getUnknownAttributes(ANGLE_ATTRIBUTE_VECTOR, atts);
	if (badAtts.size() > 0) {
	  cmlError("unknown attributes on <angle>: ");
		printVector(badAtts, cerr);
	}
	angleUnits = "degrees";
	atomRefs3Vector.clear();
// check other attributes
	for (vector<pair<string,string> >::size_type i = 0; i < atts.size(); ++i) {
		if (atts[i].first == "id") {
	    } else if (atts[i].first == C_TITLE) {
	    } else if (atts[i].first == C_CONVENTION) {
	    } else if (atts[i].first == C_ATOMREFS) {
			setCMLType(C_CML1);
			getAtomRefs(3, atomRefs3Vector, atts[i].second);
	    } else if (atts[i].first == C_ATOMREFS3) {
			setCMLType(C_CML2);
			getAtomRefs(3, atomRefs3Vector, atts[i].second);
	    } else if (atts[i].first == C_UNITS) {
			angleUnits = atts[i].second;
	    } else {
	      //			cout << "IGNORED angle attribute: " << atts[i].first << endl;
		}
	}
	return true; // [ejk] assumed success
}

bool endAngle() {
	pair <vector<OBAtom *>, double> angle;
	if (atomRefs3Vector.size() != 3) {
	  cmlError("must have defined 3 atoms for angle");
	}
	for (unsigned int i = 0; i < 3; ++i) {
		angle.first.push_back(atomRefs3Vector[i]);
	}
	angle.second = atof((char*)pcdata.c_str());
	angleVector.push_back(angle);
	return true; // [ejk] assumed
}

// not yet finished
bool WriteAngle(ostream &ofs, pair <vector<OBAtom*>, double> angle) {
    writeStartTagStart(ofs, C_ANGLE);
    string atomRefs3 = "a";
    atomRefs3 += angle.first[0]->GetIdx();
    atomRefs3 += " a";
    atomRefs3 += angle.first[1]->GetIdx();
    atomRefs3 += " a";
    atomRefs3 += angle.first[2]->GetIdx();
    writeAttribute(ofs, C_ATOMREFS3, atomRefs3);
    writeStartTagEnd(ofs);
	ofs << angle.second;
    writeEndTag(ofs, C_ANGLE);
	return true; // [ejk] assumed
}

// --------------------<atom>/<atomArray>----------------

bool startAtom(vector <pair<string,string> > &atts) {
	vector <string> badAtts = getUnknownAttributes(ATOMATTRIBUTE_VECTOR, atts);
	if (badAtts.size() > 0) {
	  cmlError("unknown attributes on <atom>: ");
		printVector(badAtts, cerr);
	}
	currentX = currentY = currentZ = 0;
	formalCharge = 0;

	atomicNum = etab.GetAtomicNum((char*)getAttribute(atts, C_ELEMENTTYPE).c_str());
	atomId = getAttribute(atts, C_ID);
	formalCharge = atoi(getAttribute(atts, C_FORMALCHARGE).c_str());
	string x2        = getAttribute(atts, C_X2);
	string y2        = getAttribute(atts, C_Y2);
	string x3        = getAttribute(atts, C_X3);
	string y3        = getAttribute(atts, C_Y3);
	string z3        = getAttribute(atts, C_Z3);
	string xFract    = getAttribute(atts, C_XFRACT);
	string yFract    = getAttribute(atts, C_YFRACT);
	string zFract    = getAttribute(atts, C_ZFRACT);
	string xy2       = getAttribute(atts, C_XY2);
	string xyz3      = getAttribute(atts, C_XYZ3);
	string xyzFract  = getAttribute(atts, C_XYZFRACT);
	if (x3 != _EMPTY) {
		currentX = atof(x3.c_str());
		setCMLType(C_CML2);
	} else if (x2 != _EMPTY) {
		currentX = atof(x2.c_str());
		setCMLType(C_CML2);
	}
	if (y3 != _EMPTY) {
		currentY = atof(y3.c_str());
		setCMLType(C_CML2);
	} else if (y2 != _EMPTY) {
		currentY = atof(y2.c_str());
		setCMLType(C_CML2);
	}
	if (z3 != _EMPTY) {
		currentZ = atof(z3.c_str());
		setCMLType(C_CML2);
	}
	if (xFract != _EMPTY) {
        cmlError("Openbabel does not support fractional coordinates");
        fractional = false;     // to avoid propagation errors
//		currentY = atof(yFract.c_str());
//		setCMLType(C_CML2);
	}
	if (yFract != _EMPTY) {
        cmlError("Openbabel does not support fractional coordinates");
        fractional = false;     // to avoid propagation errors
//		currentY = atof(yFract.c_str());
//		setCMLType(C_CML2);
	}
	if (zFract != _EMPTY) {
        cmlError("Openbabel does not support fractional coordinates");
        fractional = false;     // to avoid propagation errors
//		currentZ = atof(zFract.c_str());
//		setCMLType(C_CML2);
	}
	if (xy2 != _EMPTY) {
        vector <string> sv;
        tokenize(sv, xy2, _SPACE_NEWLINE);
        if (sv.size() != 2) {
            cmlError("xy2 attribute must have 2 floats");
        } else {
            currentX = atof(sv[0].c_str());
            currentY = atof(sv[1].c_str());
            setCMLType(C_CML2);
        }
	}
	if (xyz3 != _EMPTY) {
        vector <string> sv;
        tokenize(sv, xyz3, _SPACE_NEWLINE);
        if (sv.size() != 3) {
            cmlError("xyz3 attribute must have 3 floats");
        } else {
            currentX = atof(sv[0].c_str());
            currentY = atof(sv[1].c_str());
            currentZ = atof(sv[2].c_str());
            setCMLType(C_CML2);
        }
	}
	if (xyzFract != _EMPTY) {
        vector <string> sv;
        tokenize(sv, xyzFract, _SPACE_NEWLINE);
        if (sv.size() != 3) {
            cmlError("xyzFract attribute must have 3 floats");
        } else {
            cmlError("Openbabel does not support fractional coordinates");
            fractional = false;     // to avoid propagation errors
//            currentX = atof(sv[0].c_str());
//            currentY = atof(sv[1].c_str());
//            currentZ = atof(sv[2].c_str());
            setCMLType(C_CML2);
        }
	}
// check other attributes
	for (vector <pair<string,string> >::size_type i = 0; i < atts.size(); ++i) {
		if (atts[i].first == C_ELEMENTTYPE) {
	    } else if (atts[i].first == C_ID) {
	    } else if (atts[i].first == C_FORMALCHARGE) {
	    } else if (atts[i].first == C_X2) {
	    } else if (atts[i].first == C_Y2) {
	    } else if (atts[i].first == C_X3) {
	    } else if (atts[i].first == C_Y3) {
	    } else if (atts[i].first == C_Z3) {
	    } else if (atts[i].first == C_XFRACT) {
	    } else if (atts[i].first == C_YFRACT) {
	    } else if (atts[i].first == C_ZFRACT) {
	    } else if (atts[i].first == C_XY2) {
	    } else if (atts[i].first == C_XYZ3) {
	    } else {
	      //			cout << "IGNORED atom attribute: " << atts[i].first << endl;
		}
	}
	return true; // [ejk] assumed
}

bool processAtomArrayChild() {
	vector <string> strings;

	string builtin = getAttribute(currentAtts, C_BUILTIN);
	if (builtin == _EMPTY) {
	    cmlError("must have builtin attribute on: " + currentElem);
	}
	pcdata += "\n";
	tokenize(strings, pcdata, " \n\t");
	if (natoms == 0) {
		natoms = strings.size();
		if (natoms == 0) {
		    cmlError("no atoms in array: " + pcdata);
		}
	}
	if (static_cast<vector<string>::size_type>(natoms) != strings.size()) {
        cmlError("inconsistent atoms in arrays: " + pcdata);
	}
	for (int i = 0; i < natoms; ++i) {
		if (builtin == C_ELEMENTTYPE) {
			elementTypeVector.push_back(strings[i]);
		} else if (builtin == C_ATOMID) {
			idVector.push_back(strings[i]);
		} else if (builtin == C_FORMALCHARGE) {
			formalChargeVector.push_back(atoi((char*)strings[i].c_str()));
		} else if (builtin == C_X2) {
			x2Vector.push_back(atof((char*)strings[i].c_str()));
		} else if (builtin == C_Y2) {
			y2Vector.push_back(atof((char*)strings[i].c_str()));
		} else if (builtin == C_X3) {
			x3Vector.push_back(atof((char*)strings[i].c_str()));
		} else if (builtin == C_Y3) {
			y3Vector.push_back(atof((char*)strings[i].c_str()));
		} else if (builtin == C_Z3) {
			z3Vector.push_back(atof((char*)strings[i].c_str()));
		} else if (builtin == C_XFRACT) {
            cmlError("Openbabel does not support fractional coordinates");
            fractional = false;     // to avoid propagation errors
//			x3Vector.push_back(atof((char*)strings[i].c_str()));
		} else if (builtin == C_YFRACT) {
            cmlError("Openbabel does not support fractional coordinates");
            fractional = false;     // to avoid propagation errors
//			y3Vector.push_back(atof((char*)strings[i].c_str()));
		} else if (builtin == C_ZFRACT) {
            cmlError("Openbabel does not support fractional coordinates");
            fractional = false;     // to avoid propagation errors
//			z3Vector.push_back(atof((char*)strings[i].c_str()));
		}
	}
	return true; // [ejk] assumed
}

// adds builtin attributes for atom
bool processAtomBuiltin() {
	vector3 v;

	string builtin = getAttribute(currentAtts, C_BUILTIN);
	if (builtin == _EMPTY) {
	  cmlError("No builtin attribute for <atom><" + currentElem + ">");
		return false;
	}
	setCMLType(C_CML1);
	processBuiltinPCDATA();
	if (currentElem == C_COORDINATE2) {
		vector <double> fv;
		processFloatTokens(fv, 2, pcdata);
// 3d takes precedence over 2D
	    if (builtin == C_XY2) {
			if (cmlDimension != C_3D) {
				currentX = fv[0];
				currentY = fv[1];
			}
	    } else {
	      cmlError("IGNORED coordinate2 builtin: " + builtin);
			return false;
	    }
	} else if (currentElem == C_COORDINATE3) {
		vector <double> fv;
		processFloatTokens(fv, 3, pcdata);
// 3d takes precedence over 2D
	    if (builtin == C_XYZ3) {
			currentX = fv[0];
			currentY = fv[1];
			currentZ = fv[2];
	    } else if (builtin == C_XYZFRACT) {
			currentX = fv[0];
			currentY = fv[1];
			currentZ = fv[2];
			fractional = true;
	    } else {
	      cmlError("IGNORED coordinate2 builtin: " + builtin);
			return false;
	    }
	} else if (currentElem == C_FLOAT) {
		double value = atof(pcdata.c_str());
// 3d takes precedence over 2D
	    if (builtin == C_X2) {
			if (cmlDimension != C_3D) currentX = value;
	    } else if (builtin == C_Y2) {
			if (cmlDimension != C_3D) currentY = value;
	    } else if (builtin == C_X3) {
			cmlDimension = C_3D;
			currentX = value;
	    } else if (builtin == C_Y3) {
			cmlDimension = C_3D;
			currentY = value;
	    } else if (builtin == C_Z3) {
			cmlDimension = C_3D;
			currentZ = value;
	    } else if (builtin == C_XFRACT) {
			cmlDimension = C_3D;
			currentX = value;
            fractional = true;
	    } else if (builtin == C_YFRACT) {
			cmlDimension = C_3D;
			currentY = value;
            fractional = true;
	    } else if (builtin == C_ZFRACT) {
			cmlDimension = C_3D;
			currentZ = value;
            fractional = true;
	    } else {
	      cmlError("IGNORED float builtin: " + builtin);
			return false;
	    }
	} else if (currentElem == C_INTEGER) {
		int ival = atoi(pcdata.c_str());
	    if (builtin == C_FORMALCHARGE) {
			formalCharge = ival;
	    } else {
	      cmlError("IGNORED integer builtin: " + builtin);
			return false;
		}
	} else if (currentElem == C_STRING) {
	    if (builtin == C_ELEMENTTYPE) {
  	        atomicNum = etab.GetAtomicNum((char*)pcdata.c_str());
	    } else if (builtin == C_ATOMID) {
			atomId = pcdata;
	    } else {
	      cmlError("IGNORED string builtin: " + builtin);
			return false;
		}
	}
	return true;
}

bool endAtom() {
// seems to be all Openbabel in this routine    
	OBAtom atom;
	pair<string, OBAtom*> at;

	atom.SetAtomicNum(atomicNum);
	atom.SetFormalCharge(formalCharge);
    if (fractional) {
        cmlError("Openbabel does not support fractional coordinates");
        fractional = false;     // to avoid propagation errors
//        atom.SetVector(currentX, currentY, currentZ);
    } else {
        atom.SetVector(currentX, currentY, currentZ);
    }

    molPtr->AddAtom(atom);
    int nat = molPtr->NumAtoms();
    OBAtom* atPtr = molPtr->GetAtom(nat);	// counts from 1?
// store atoms
	at.first = atomId;
	at.second = atPtr;
	atomIdVector.push_back(at);
	return true; // [ejk] assumed
}

bool WriteAtom(ostream &ofs, OBAtom* atom, int count) {
	double x, y, z;
	char* elementType;
	int charge;
	char ids[8];

	charge = atom->GetFormalCharge();
	x = atom->GetX(),
	y = atom->GetY(),
	z = atom->GetZ(),
	elementType = etab.GetSymbol(atom->GetAtomicNum());
//.there must be an easier way to concatenate a number and a string...
	string id = "a";
	sprintf(ids, "%i", count);
	string cs(ids);
	id.append(trim(cs));
	if (!outputArray) {
        writeStartTagStart(ofs, C_ATOM);
		writeAttribute(ofs, C_ID, id);
// CML2
		if (outputCML2) {
			writeAttribute(ofs, C_ELEMENTTYPE, elementType);
			if (charge != 0) writeAttribute(ofs, C_FORMALCHARGE, charge);
			if (molPtr->HasNonZeroCoords()) {
				if (strcmp(dimension, C_2D) == 0) {
					writeAttribute(ofs, C_X2, x);
					writeAttribute(ofs, C_Y2, y);
				} else if (strcmp(dimension, C_3D) == 0) {
// should never get to this until OB is changed
                    if (fractional) {
                        cmlError("Openbabel does not support fractional coordinates");
                        writeAttribute(ofs, C_XFRACT, x);
                        writeAttribute(ofs, C_YFRACT, y);
                        writeAttribute(ofs, C_ZFRACT, z);
                    } else {
                        writeAttribute(ofs, C_X3, x);
                        writeAttribute(ofs, C_Y3, y);
                        writeAttribute(ofs, C_Z3, z);
                    }
				}
			}
            writeCombinedTagEnd(ofs);
// CML1
		} else {
            writeStartTagEnd(ofs);
            ofs << endl;
			writeBuiltin(ofs, C_ELEMENTTYPE, elementType);
			if (charge != 0) writeBuiltin(ofs, C_FORMALCHARGE, charge);
			if (molPtr->HasNonZeroCoords()) {
				if (strcmp(dimension, C_2D) == 0) {
					writeBuiltin(ofs, C_X2, x);
					writeBuiltin(ofs, C_Y2, y);
				} else if (strcmp(dimension, C_3D) == 0) {
                    if (fractional) {
                        writeBuiltin(ofs, C_XFRACT, x);
                        writeBuiltin(ofs, C_YFRACT, y);
                        writeBuiltin(ofs, C_ZFRACT, z);
                    } else {
                        writeBuiltin(ofs, C_X3, x);
                        writeBuiltin(ofs, C_Y3, y);
                        writeBuiltin(ofs, C_Z3, z);
                    }                        
				}
			}
            writeEndTag(ofs, C_ATOM);
		}
	} else {
		appendToArray(idArray, id);
		appendToArray(elementArray, elementType);
		appendToArray(chargeArray, charge);
		if (molPtr->HasNonZeroCoords()) {
			if (strcmp(dimension, C_2D) == 0) {
				appendToArray(x2Array, x);
				appendToArray(y2Array, y);
			} else if (strcmp(dimension, C_3D) == 0) {
                if (fractional) {
// should never get here                    
                    cmlError("Openbabel does not support fractional coordinates");
//                    appendToArray(x3Array, x);
//                    appendToArray(y3Array, y);
//                    appendToArray(z3Array, z);
                } else {
                    appendToArray(x3Array, x);
                    appendToArray(y3Array, y);
                    appendToArray(z3Array, z);
                }
			}
		}
	}
	return true; // [ejk] assumed OK
}

void processStringTokens(vector <string> &v, vector<string>::size_type n, string att) {
	if (att == _EMPTY) return;
	vector <string> sv;
	att += _SPACE;
	tokenize(sv, att, _SPACE_NEWLINE);
	if (sv.size() != n) {
	    cmlError("inconsistent array attribute sizes: ");
        cout << sv.size() << _SLASH << n << endl;
		return;
	}
	for (vector<string>::size_type i = 0; i < n; ++i) {
        //v[i] = sv[i];
	    v.push_back(sv[i]);
    }
}

void processIntTokens(vector <int> &v, vector<int>::size_type n, string att) {
	if (att == _EMPTY) {
        return;
    }
	vector <string> sv;
	att += _SPACE;
	tokenize(sv, att, _SPACE_NEWLINE);
	if (sv.size() != n) {
	  cmlError("inconsistent array attribute sizes: ");
      cerr << sv.size() << _SLASH << n << endl;
		return;
	}
	for (vector<int>::size_type i = 0; i < n; ++i) {
        v.push_back(atoi((char*)sv[i].c_str()));
    }
}

void processFloatTokens(vector <double> &v, vector<double>::size_type n, string att) {
	if (att == _EMPTY) return;
	vector <string> sv;
	att += _SPACE;
	tokenize(sv, att, _SPACE_NEWLINE);
	if (sv.size() != n) {
	  cmlError("inconsistent array attribute sizes: ");
      cerr << sv.size() << _SLASH << n << endl;
		return;
	}
	for (vector<double>::size_type i = 0; i < n; ++i) {
        v.push_back(atof((char*)sv[i].c_str()));
    }
}

bool startAtomArray(vector <pair<string,string> > &atts) {
	vector <string> sv;
	string atomID = getAttribute(atts, C_ATOMID);
    
// atomArray with attributes => CML2+array
// everything else exits here
	if (atomID == _EMPTY) {
        return false;
    }
// only CML2+array gets to here    
	setCMLType(C_CML2);
    inputArray = true;
	atomId += _SPACE;
	tokenize(sv, atomID, _SPACE_NEWLINE);
	int mynatoms = sv.size();
    if (mynatoms == 0) {
        cmlError("startAtomArray: No atoms given");
        return false;
    }
    natoms = mynatoms;
	processStringTokens(idVector, mynatoms, atomID);
	processStringTokens(elementTypeVector, mynatoms, getAttribute(atts, C_ELEMENTTYPE));
	string fCharge = getAttribute(atts, C_FORMALCHARGE);
	processIntTokens(formalChargeVector, mynatoms, fCharge);
	processIntTokens(hydrogenCountVector, mynatoms, getAttribute(atts, C_HYDROGENCOUNT));
	processFloatTokens(x2Vector, mynatoms, getAttribute(atts, C_X2));
	processFloatTokens(y2Vector, mynatoms, getAttribute(atts, C_Y2));
	processFloatTokens(x3Vector, mynatoms, getAttribute(atts, C_X3));
	processFloatTokens(y3Vector, mynatoms, getAttribute(atts, C_Y3));
	processFloatTokens(z3Vector, mynatoms, getAttribute(atts, C_Z3));
	if ("" != getAttribute(atts, C_XYZ3) ||
    "" != getAttribute(atts, C_XY2)) {
        cmlError("attributes xyz3 and xy2 not supported in CML2 array mode");
    }
	return true; // [ejk] assumed
}

bool endAtomArray() {
	if (cmlType == C_CML2 || inputArray) {
		for (int i = 0; i < natoms; ++i) {
			OBAtom atom;
            pair<string, OBAtom*> at;
			if (elementTypeVector.size() > 0) {
                atom.SetAtomicNum(etab.GetAtomicNum((char*)elementTypeVector[i].c_str()));
			}
			if (formalChargeVector.size() > 0) atom.SetFormalCharge(formalChargeVector[i]);
			vector3 v;
			if (x2Vector.size() > 0) v.SetX(x2Vector[i]);
			if (y2Vector.size() > 0) v.SetY(y2Vector[i]);
			if (x3Vector.size() > 0) v.SetX(x3Vector[i]);
			if (y3Vector.size() > 0) v.SetY(y3Vector[i]);
			if (z3Vector.size() > 0) v.SetZ(z3Vector[i]);
			atom.SetVector(v);
		    molPtr->AddAtom(atom);
		    OBAtom* atPtr = molPtr->GetAtom(i+1);	// counts from 1?
// store atoms
			at.first = idVector[i];
			at.second = atPtr;
			atomIdVector.push_back(at);
		}
	}
	return true;
}

bool WriteAtomArray(ostream &ofs) {
	OBAtom *atom;
	vector<OBNodeBase*>::iterator i;

	int count = 0;
    writeStartTagStart(ofs, C_ATOMARRAY);
	if (!outputArray) {
        writeStartTagEnd(ofs);
        ofs << endl;
	}
	for (atom = molPtr->BeginAtom(i);atom;atom = molPtr->NextAtom(i)) {
		WriteAtom(ofs, atom, ++count);
	}
	if (outputArray) {
// CML1 array
		if (outputCML1) {
            writeStartTagEnd(ofs);
            ofs << endl;

            writeStartTagStart(ofs, C_STRINGARRAY);
            writeAttribute(ofs, C_BUILTIN, C_ATOMID);
            writeStartTagEnd(ofs);
            ofs << idArray;
            writeEndTag(ofs, C_STRINGARRAY);

            writeStartTagStart(ofs, C_STRINGARRAY);
            writeAttribute(ofs, C_BUILTIN, C_ELEMENTTYPE);
            writeStartTagEnd(ofs);
            ofs << elementArray;
            writeEndTag(ofs, C_STRINGARRAY);

            writeStartTagStart(ofs, C_INTEGERARRAY);
            writeAttribute(ofs, C_BUILTIN, C_FORMALCHARGE);
            writeStartTagEnd(ofs);
            ofs << chargeArray;
            writeEndTag(ofs, C_INTEGERARRAY);

			if (molPtr->HasNonZeroCoords()) {
				if (strcmp(dimension, C_2D) == 0) {

                    writeStartTagStart(ofs, C_FLOATARRAY);
                    writeAttribute(ofs, C_BUILTIN, C_X2);
                    writeStartTagEnd(ofs);
                    ofs << x2Array;
                    writeEndTag(ofs, C_FLOATARRAY);

                    writeStartTagStart(ofs, C_FLOATARRAY);
                    writeAttribute(ofs, C_BUILTIN, C_Y2);
                    writeStartTagEnd(ofs);
                    ofs << y2Array;
                    writeEndTag(ofs, C_FLOATARRAY);

				} else if (strcmp(dimension, C_3D) == 0) {
                    if (fractional) {
// should never get here in Openbabel
                        cmlError("Openbabel does not support fractional coordinates");
                        writeStartTagStart(ofs, C_FLOATARRAY);
                        writeAttribute(ofs, C_BUILTIN, C_XFRACT);
                        writeStartTagEnd(ofs);
                        ofs << x3Array;
                        writeEndTag(ofs, C_FLOATARRAY);
    
                        writeStartTagStart(ofs, C_FLOATARRAY);
                        writeAttribute(ofs, C_BUILTIN, C_YFRACT);
                        writeStartTagEnd(ofs);
                        ofs << y3Array;
                        writeEndTag(ofs, C_FLOATARRAY);
    
                        writeStartTagStart(ofs, C_FLOATARRAY);
                        writeAttribute(ofs, C_BUILTIN, C_ZFRACT);
                        writeStartTagEnd(ofs);
                        ofs << z3Array;
                        writeEndTag(ofs, C_FLOATARRAY);
                    } else {
                        writeStartTagStart(ofs, C_FLOATARRAY);
                        writeAttribute(ofs, C_BUILTIN, C_X3);
                        writeStartTagEnd(ofs);
                        ofs << x3Array;
                        writeEndTag(ofs, C_FLOATARRAY);
    
                        writeStartTagStart(ofs, C_FLOATARRAY);
                        writeAttribute(ofs, C_BUILTIN, C_Y3);
                        writeStartTagEnd(ofs);
                        ofs << y3Array;
                        writeEndTag(ofs, C_FLOATARRAY);
    
                        writeStartTagStart(ofs, C_FLOATARRAY);
                        writeAttribute(ofs, C_BUILTIN, C_Z3);
                        writeStartTagEnd(ofs);
                        ofs << z3Array;
                        writeEndTag(ofs, C_FLOATARRAY);
                    }
				}
			}
            writeEndTag(ofs, C_ATOMARRAY);
// CML2 array
		} else {
            writeAttribute(ofs, C_ATOMID, idArray);
            writeAttribute(ofs, C_ELEMENTTYPE, elementArray);
            writeAttribute(ofs, C_FORMALCHARGE, chargeArray);
			if (molPtr->HasNonZeroCoords()) {
				if (strcmp(dimension, C_2D) == 0) {
                    writeAttribute(ofs, C_X2, x2Array);
                    writeAttribute(ofs, C_Y2, y2Array);
				} else if (strcmp(dimension, C_3D) == 0) {
                    if (fractional) {
// should never get here in Openbabel
                        cmlError("Openbabel does not support fractional coordinates");
                        writeAttribute(ofs, C_XFRACT, x3Array);
                        writeAttribute(ofs, C_YFRACT, y3Array);
                        writeAttribute(ofs, C_ZFRACT, z3Array);
                    } else {
                        writeAttribute(ofs, C_X3, x3Array);
                        writeAttribute(ofs, C_Y3, y3Array);
                        writeAttribute(ofs, C_Z3, z3Array);
                    }
				}
			}
            writeCombinedTagEnd(ofs);
		}
	} else {
        writeEndTag(ofs, C_ATOMARRAY);
	}
	return true; // [ejk] assumed
}

bool startAtomParity(vector <pair<string,string> > &atts) {
	atomRefs4 = getAttribute(atts, C_ATOMREFS4);
	return true; // [ejk] assumed
}

bool endAtomParity(vector <pair<string,string> > &atts) {
	pair <vector<OBAtom*>, double> ap;
	vector <OBAtom*> atomRef;
	getAtomRefs(4, atomRef, atomRefs4);
	if (atomRef.size() != 4) {
	  cmlError("atomRefs4 must reference 4 atoms");
		return false;
	}
	for (int i = 0; i < 4; ++i) ap.first.push_back(atomRef[i]);
	setCMLType(C_CML2);
	ap.second = atof((char*)pcdata.c_str());
	atomParityVector.push_back(ap);
	return true; // [ejk] assumed
}

bool WriteAtomParity(ostream &ofs) {
  //	cout << "WriteAtomParity NYI" << endl;
	return true; // [ejk] assumed
}

// --------------------<bond><bondArray>----------------
bool startBond(vector <pair<string,string> > &atts) {
	vector <string> badAtts = getUnknownAttributes(BOND_ATTRIBUTE_VECTOR, atts);
	if (badAtts.size() > 0) {
  	    cmlError("unknown attributes on <bond>:");
		printVector(badAtts, cerr);
	}

	vector <string> atomRefs;

	bondBeginAtom = _EMPTY;
	bondEndAtom = _EMPTY;
//    bondIdString = getAttribute(currentAtts, C_ID); // Babel can't store IDs
    orderString = getAttribute(currentAtts, C_ORDER);
    stereoString = getAttribute(currentAtts, C_STEREO);
	tokenize(atomRefs, (char*)getAttribute(currentAtts, C_ATOMREFS2).c_str(), " \n\t,");
	if (atomRefs.size() == 0) {
		return false;
	} else if (atomRefs.size() != 2) {
	  cmlError("must have 2 atom Refs per bond");
		return false;
	} else {
		setCMLType(C_CML2);
	}
	bondBeginAtom = atomRefs[0];
	bondEndAtom = atomRefs[1];
	return true; // [ejk] assumed
}

// adds builtin attributes for bond
bool processBondBuiltin() {

	string builtin = getAttribute(currentAtts, C_BUILTIN);
	if (builtin == _EMPTY) {
	    cmlError("No builtin attribute for <bond><"+currentElem+">");
		return false;
	}
	setCMLType(C_CML1);
	if (currentElem == C_FLOAT) {
		double value = atof(pcdata.c_str());
	    if (builtin == C_LENGTH) {
			length = value;
	    } else {
	      cmlError("IGNORED float builtin for bond: "+builtin);
			return false;
		}
	} else if (currentElem == C_INTEGER) {
		int ival = atoi(pcdata.c_str());
		return false;
	} else if (currentElem == C_STRING) {
	    if (builtin == C_ATOMREF) {
			if (bondBeginAtom == _EMPTY) {
				bondBeginAtom = pcdata;
			} else {
				if (bondEndAtom != _EMPTY) {
				  //					cerr << "too many atomRef builtins" << endl;
					return false;
				}
				bondEndAtom = pcdata;
			}
// maybe id might occur as child, but better as attribute of bond parent
// anyway Babel can't store it
//		} else if (builtin == C_ID) {
//			idString = pcdata;
		} else if (builtin == C_ORDER) {
			orderString = pcdata;
	    } else if (builtin == C_STEREO) {
			stereoString = pcdata;
	    } else {
	      cmlError("IGNORED integer builtin: "+builtin);
			return false;
		}
	}
	return true;
}


bool processBondArrayChild() {
	vector <string> strings;

	string builtin = getAttribute(currentAtts, C_BUILTIN);
	if (builtin == _EMPTY) {
	    cmlError("must have builtin attribute on: "+currentElem);
	}
	pcdata += "\n";
	tokenize(strings, pcdata, " \n\t");
	if (nbonds == 0) {
		nbonds = strings.size();
		if (nbonds == 0) {
		    cmlError("no bonds in array: "+pcdata);
		}
	}
	if (nbonds != strings.size()) {
	    cmlError("inconsistent bonds in arrays: "+pcdata);
	}
	bool atomRef1 = (atomRef1Vector.size() == 0);
	for (unsigned int i = 0; i < nbonds; ++i) {
		if (builtin == "atomRef") {
			if (atomRef1) {
				atomRef1Vector.push_back(strings[i]);
			} else {
				atomRef2Vector.push_back(strings[i]);
			}
// this is the correct way of holding bond Ids in CML1arrays            
		} else if (builtin == C_ID) {
			idVector.push_back(strings[i]);
		} else if (builtin == C_ORDER) {
			orderVector.push_back(strings[i]);
		} else if (builtin == C_STEREO) {
			stereoVector.push_back(strings[i]);
		}
	}
	return true; // [ejk] assumed
}

bool endBond() {
	pair <vector<OBAtom*>, double> len;
	OBBond bond;

    bondPtr = &bond;
	OBAtom* beginAtomPtr = getAtomPtr(bondBeginAtom);
	OBAtom* endAtomPtr = getAtomPtr(bondEndAtom);
	if (beginAtomPtr == 0 || endAtomPtr == 0) {
	    cmlError("could not find atom refs in bond");
		return false;
	}
    bondPtr->SetBegin(beginAtomPtr);
    bondPtr->SetEnd(endAtomPtr);
//    if (idString != _EMPTY) bondPtr->SetID(idString);     // no IDs in Babel
    if (orderString != _EMPTY) bondPtr->SetBO(getBabelBondOrder(orderString));
    if (stereoString == "W") {
		bondPtr->SetUp();
	} else if (stereoString == "H") {
		bondPtr->SetDown();
	}
// length is property of molecule
	if (length >= 0) {
		len.first.push_back(beginAtomPtr);
		len.first.push_back(endAtomPtr);
		len.second = length;
		lengthVector.push_back(len);
	}
    molPtr->AddBond(*bondPtr);
	return true;
}

bool WriteBond(ostream &ofs, OBBond* bond) {
	int a1;
	int a2;
	int bo;
	char bos[8];
	char* boChar;
	a1 = bond->GetBeginAtomIdx();
	a2 = bond->GetEndAtomIdx();
	bo = bond->GetBO();
	switch (bo) {
		case 0: boChar = ""; break;
		case 1: boChar = "1"; break;
		case 2: boChar = "2"; break;
		case 3: boChar = "3"; break;
		case 5: boChar = "A"; break;
		default: boChar = ""; break;
	}
	string atomRef1 = "a";
	sprintf(bos, "%i", a1);
	string bos1(bos);
	atomRef1.append(trim(bos1));
	string atomRef2 = "a";
	sprintf(bos, "%i", a2);
	string bos2(bos);
	atomRef2.append(trim(bos2));
	if (!outputArray) {
        writeStartTagStart(ofs, C_BOND);
		if (outputCML2) {
			string atomRefs2 = atomRef1+_SPACE+atomRef2;
			writeAttribute(ofs, C_ATOMREFS2, atomRefs2);
			writeAttribute(ofs, C_ORDER, boChar);
            writeCombinedTagEnd(ofs);
		} else {
            writeStartTagEnd(ofs);
            ofs << endl;
			writeBuiltin(ofs, "atomRef", atomRef1);
			writeBuiltin(ofs, "atomRef", atomRef2);
			writeBuiltin(ofs, C_ORDER, boChar);
            writeEndTag(ofs, C_BOND);
		}
	} else {
		appendToArray(atomRef1Array, atomRef1);
		appendToArray(atomRef2Array, atomRef2);
		appendToArray(orderArray, boChar);
	}
	return true; // [ejk] assumed
}

bool startBondArray(vector <pair<string,string> > &atts) {
	vector <string> sv;
	string atomRef1 = getAttribute(atts, C_ATOMREFS1);
	if (atomRef1 == _EMPTY) {
        return false;
    }
// only CML2+array gets to here    
	setCMLType(C_CML2);
    inputArray = true;
	atomRef1 += _SPACE;
	tokenize(sv, atomRef1, _SPACE_NEWLINE);
	int mynbonds = sv.size();	// explicitly not the global nbonds
    if (mynbonds == 0) {
        cmlError("startBondArray: No bonds given");
        return false;
    }
	processStringTokens(atomRef1Vector, mynbonds, atomRef1);
	processStringTokens(atomRef2Vector, mynbonds, getAttribute(atts,  C_ATOMREFS2));
	processStringTokens(orderVector, mynbonds, getAttribute(atts, C_ORDER));
	processStringTokens(stereoVector, mynbonds, getAttribute(atts, C_STEREO));
    nbonds = mynbonds;
	return true; // [ejk] assumed
}

bool endBondArray() {
	if (/*cmlType == C_CML2 || */inputArray) {
		if (atomRef1Vector.size() == 0 ||
			atomRef2Vector.size() == 0) {
		  cmlError("atomRef arrays must be given for bonds");
		}
		for (unsigned int i = 0; i < nbonds; ++i) {
			OBBond bond;
			bondPtr = &bond;
			OBAtom* beginAtomPtr = getAtomPtr(atomRef1Vector[i]);
			OBAtom* endAtomPtr = getAtomPtr(atomRef2Vector[i]);
			if (beginAtomPtr == 0 || endAtomPtr == 0) {
			  cmlError("could not find atom refs in bond");
				return false;
			}
		    bondPtr->SetBegin(beginAtomPtr);
		    bondPtr->SetEnd(endAtomPtr);
		    if (orderVector.size() > 0) bondPtr->SetBO(getBabelBondOrder(orderVector[i]));
		    if (stereoVector.size() > 0) {
			    if (stereoVector[i] == "W") {
					bondPtr->SetUp();
				} else if (stereoVector[i] == "H") {
					bondPtr->SetDown();
				}
			}
		    molPtr->AddBond(*bondPtr);
		}
	}
	return true;
}

bool WriteBondArray(ostream &ofs) {
	if (molPtr->NumBonds() == 0) return false;

    writeStartTagStart(ofs, C_BONDARRAY);
	if (!outputArray) {
        writeStartTagEnd(ofs);
        ofs << endl;
	}
	//so the bonds come out sorted

	OBAtom *atom;
	OBAtom *nbr;
	OBBond *bond;
	vector<OBNodeBase*>::iterator i;
	vector<OBEdgeBase*>::iterator j;
	for (atom = molPtr->BeginAtom(i);atom;atom = molPtr->NextAtom(i)) {
		for (nbr = atom->BeginNbrAtom(j);nbr;nbr = atom->NextNbrAtom(j)) {
			if (atom->GetIdx() < nbr->GetIdx()) {
				bond = (OBBond*) *j;
				WriteBond(ofs, bond);
			}
		}
	}
	if (outputArray) {
		if (outputCML1) {
            writeStartTagEnd(ofs);
            ofs << endl;

            writeStartTagStart(ofs, C_STRINGARRAY);
            writeAttribute(ofs, C_BUILTIN, C_ATOMREF1);
            writeStartTagEnd(ofs);
			ofs << atomRef1Array;
            writeEndTag(ofs, C_STRINGARRAY);

            writeStartTagStart(ofs, C_STRINGARRAY);
            writeAttribute(ofs, C_BUILTIN, C_ATOMREF2);
            writeStartTagEnd(ofs);
			ofs << atomRef2Array;
            writeEndTag(ofs, C_STRINGARRAY);

            writeStartTagStart(ofs, C_STRINGARRAY);
            writeAttribute(ofs, C_BUILTIN, C_ORDER);
            writeStartTagEnd(ofs);
			ofs << orderArray;
            writeEndTag(ofs, C_STRINGARRAY);

            writeEndTag(ofs, C_BONDARRAY);
		} else {
            writeAttribute(ofs, C_ATOMREF1, atomRef1Array);
            writeAttribute(ofs, C_ATOMREF2, atomRef2Array);
            writeAttribute(ofs, C_ORDER, orderArray);

            writeCombinedTagEnd(ofs);
		}
	} else {
        writeEndTag(ofs, C_BONDARRAY);
	}
	return true; // [ejk] assumed
}

// numeric value (babel) from CML order; unknown returns -1
int getBabelBondOrder(string o) {
	if (o == "1" || o == "S") return 1;
	if (o == "2" || o == "D") return 2;
	if (o == "3" || o == "T") return 3;
	if (o == "A") return 5;
	return -1;
}

// numeric value (babel) from CML stereo; unknown returns -1
int getBabelBondFlag(string s) {
	if (s == "W") return OB_WEDGE_BOND;
	if (s == "H") return OB_HASH_BOND;
	return -1;
}

// --------------------<crystal>----------------

bool startCrystal(vector <pair<string,string> > &atts) {
	vector <string> badAtts = getUnknownAttributes(CRYSTAL_ATTRIBUTE_VECTOR, atts);
	if (badAtts.size() > 0) {
	  cmlError("unknown attributes on <crystal>: ");
		printVector(badAtts, cerr);
	}
// check other attributes
	for (vector <pair<string,string> >::size_type i = 0; i < atts.size(); ++i) {
		if (atts[i].first == C_ID) {
	    } else if (atts[i].first == C_TITLE) {
	    } else if (atts[i].first == C_CONVENTION) {
	    } else if (atts[i].first == C_SPACEGROUP) {
	    } else if (atts[i].first == C_POINTGROUP) {
	    } else {
	      //			cout << "IGNORED crystal attribute: " << atts[i].first << endl;
		}
	}
	return true; // [ejk] assumed
}

// adds builtin attributes for cryst
bool processCrystalBuiltin() {
	vector3 v;

	string builtin = getAttribute(currentAtts, C_BUILTIN);
	if (builtin == _EMPTY) {
	  cmlError("No builtin attribute for <cryst><"+currentElem+">");
		return false;
	}
	setCMLType(C_CML1);
	processBuiltinPCDATA();
	if (currentElem == C_FLOAT) {
		double f = atof((char*)pcdata.c_str());
		if (currentElem == "acell") {
			cellParam[0] = f;
		} else if (currentElem == "bcell") {
			cellParam[1] = f;
		} else if (currentElem == "ccell") {
			cellParam[2] = f;
		} else if (currentElem == "alpha") {
			cellParam[3] = f;
		} else if (currentElem == "beta") {
			cellParam[4] = f;
		} else if (currentElem == "gamma") {
			cellParam[5] = f;
		} else {
		  cmlError("IGNORED float builtin: "+builtin);
			return false;
		}
	} else {
	  cmlError("IGNORED builtin for "+currentElem+" in crystal; "+builtin);
	}
	return true;
}

bool endCrystal() {
	return true; // [ejk] assumed
}

// not yet finished
bool WriteCrystal(ostream &ofs) {
    writeStartTagStart(ofs, C_CRYSTAL);
    writeStartTagEnd(ofs);
    writeEndTag(ofs, C_CRYSTAL);
	return true; // [ejk] assumed
}

// --------------------<electron>----------------

bool startElectron(vector <pair<string,string> > &atts) {
	vector <string> badAtts = getUnknownAttributes(ELECTRON_ATTRIBUTE_VECTOR, atts);
	if (badAtts.size() > 0) {
	  cmlError("unknown attributes on <electron>: ");
		printVector(badAtts, cerr);
	}
// check other attributes
	for (vector <pair<string,string> >::size_type i = 0; i < atts.size(); ++i) {
		if (atts[i].first == "id") {
	    } else if (atts[i].first == C_TITLE) {
	    } else if (atts[i].first == C_CONVENTION) {
	    } else {
	      //			cout << "IGNORED electron attribute: " << atts[i].first << endl;
		}
	}
	return true; // [ejk] assumed
}

bool endElectron() {
//	pair <vector<OBAtom*>, double> electron;
//	electronVector.push_back(electron);
	return true; // [ejk] assumed
}

// not yet finished
bool WriteElectron(ostream &ofs) {
    writeStartTagStart(ofs, C_ELECTRON);
    writeStartTagEnd(ofs);
    writeEndTag(ofs, C_ELECTRON);
	return true; // [ejk] assumed
}

// --------------------<feature>----------------

bool startFeature(vector <pair<string,string> > &atts) {
	vector <string> badAtts = getUnknownAttributes(FEATURE_ATTRIBUTE_VECTOR, atts);
	if (badAtts.size() > 0) {
	  cmlError("unknown attributes on <feature>: ");
		printVector(badAtts, cerr);
	}
// check other attributes
	for (vector<pair<string,string> >::size_type i = 0; i < atts.size(); ++i) {
		if (atts[i].first == "id") {
	    } else if (atts[i].first == C_TITLE) {
	    } else if (atts[i].first == C_CONVENTION) {
	    } else {
	      //			cout << "IGNORED feature attribute: " << atts[i].first << endl;
		}
	}
	return true; // [ejk] assumed
}

bool endFeature() {
	return true; // [ejk] assumed
}

bool WriteFeature(ostream &ofs) {
    writeStartTagStart(ofs, C_FEATURE);
    writeStartTagEnd(ofs);
    writeEndTag(ofs, C_FEATURE);
	return true; // [ejk] assumed
}

// --------------------<formula>----------------

bool startFormula(vector <pair<string,string> > &atts) {
	vector <string> badAtts = getUnknownAttributes(FORMULA_ATTRIBUTE_VECTOR, atts);
	if (badAtts.size() > 0) {
	  cmlError("unknown attributes on <formula>: ");
		printVector(badAtts, cerr);
	}
// check other attributes
	for (vector <pair<string,string> >::size_type i = 0; i < atts.size(); ++i) {
		if (atts[i].first == "id") {
	    } else if (atts[i].first == C_TITLE) {
	    } else if (atts[i].first == C_CONVENTION) {
	    } else {
	      //			cout << "IGNORED formula attribute: " << atts[i].first << endl;
		}
	}
	return true; // [ejk] assumed
}

bool endFormula() {
	return true; // [ejk] assumed
}

bool WriteFormula(ostream &ofs) {
    writeStartTagStart(ofs, C_FORMULA);
    writeStartTagEnd(ofs);
    writeEndTag(ofs, C_FORMULA);
	return true; // [ejk] assumed
}

// --------------------<length>----------------

bool startLength(vector <pair<string,string> > &atts) {
	vector <string> badAtts = getUnknownAttributes(LENGTH_ATTRIBUTE_VECTOR, atts);
	if (badAtts.size() > 0) {
	  cmlError("unknown attributes on <length>: ");
		printVector(badAtts, cerr);
	}
	lengthUnits = "angstrom";
	atomRefs2Vector.clear();
// check other attributes
	for (vector <pair<string,string> >::size_type i = 0; i < atts.size(); ++i) {
		if (atts[i].first == "id") {
	    } else if (atts[i].first == C_TITLE) {
	    } else if (atts[i].first == C_CONVENTION) {
	    } else if (atts[i].first == C_ATOMREFS2) {
			getAtomRefs(2, atomRefs2Vector, atts[i].second);
	    } else if (atts[i].first == C_UNITS) {
			lengthUnits = atts[i].second;
	    } else {
	      //			cout << "IGNORED length attribute: " << atts[i].first << endl;
		}
	}
	return true; // [ejk] assumed
}

bool endLength() {
	pair <vector<OBAtom*>, double> length;
	if (atomRefs2Vector.size() != 2) {
	  cmlError("must have defined 2 atoms for length");
	}
	for (int i = 0; i < 2; ++i) {
		length.first.push_back(atomRefs2Vector[i]);
	}
	length.second = atof((char*)pcdata.c_str());
	lengthVector.push_back(length);
	return true; // [ejk] assumed
}

// not yet finished
bool WriteLength(ostream &ofs, pair <vector<OBAtom*>, double> length) {
    writeStartTagStart(ofs, C_LENGTH);
    string atomRefs2 = "a";
    atomRefs2 += length.first[0]->GetIdx();
    atomRefs2 += " a";
    atomRefs2 += length.first[1]->GetIdx();
    writeAttribute(ofs, C_ATOMREFS2, atomRefs2);
    writeStartTagEnd(ofs);
    writeEndTag(ofs, C_LENGTH);
	return true; // [ejk] assumed
}

// ------------------ <molecule> ----------------

bool startMolecule(vector <pair<string,string> > &atts) {
	vector <string> badAtts = getUnknownAttributes(MOLECULE_ATTRIBUTE_VECTOR, atts);
	if (badAtts.size() > 0) {
	  cmlError("unknown attributes on <molecule>: ");
	  //		printVector(badAtts, cerr);
	}

	molPtr->BeginModify();
	molPtr->ReserveAtoms(ATOMSIZE);
	molPtr->SetTitle((char*)getAttribute(atts, C_TITLE).c_str());
    
	return true; // [ejk] assumed
}

bool debugMolecule(ostream &ofs) {
    dimension = C_3D;
	OBAtom *atom;
	vector<OBNodeBase*>::iterator i;
	int count = 0;
	for (atom = molPtr->BeginAtom(i);atom;atom = molPtr->NextAtom(i)) {
		WriteAtom(ofs, atom, ++count);
	}
	OBAtom *nbr;
	OBBond *bond;
	vector<OBEdgeBase*>::iterator j;
	for (atom = molPtr->BeginAtom(i);atom;atom = molPtr->NextAtom(i)) {
		for (nbr = atom->BeginNbrAtom(j);nbr;nbr = atom->NextNbrAtom(j)) {
			if (atom->GetIdx() < nbr->GetIdx()) {
				bond = (OBBond*) *j;
				WriteBond(ofs, bond);
			}
		}
	}
	if (lengthVector.size() > 0) {
		ofs << "Lengths: " << endl;
		for (unsigned int i = 0; i < lengthVector.size(); ++i) {
			pair <vector<OBAtom*>, double> length = lengthVector[i];
			WriteLength(ofs, length);
		}
	}
	if (angleVector.size() > 0) {
		ofs << "Angles: " << endl;
		for (unsigned int i = 0; i < angleVector.size(); ++i) {
			pair <vector<OBAtom*>, double> angle = angleVector[i];
			WriteAngle(ofs, angle);
		}
	}
	if (torsionVector.size() > 0) {
		ofs << "Torsions: " << endl;
		for (unsigned int i = 0; i < torsionVector.size(); ++i) {
			pair <vector<OBAtom*>, double> torsion = torsionVector[i];
			WriteTorsion(ofs, torsion);
		}
	}
	return true; // [ejk] assumed
}

// returns index of length between atoms ( starts at 0; -1 = not found)
int getLengthIndex(OBAtom* a0, OBAtom* a1) {
	for (unsigned int i = 0; i < lengthVector.size(); ++i) {
		if (a0 == lengthVector[i].first[0] &&
			a1 == lengthVector[i].first[1]) return i;
		if (a0 == lengthVector[i].first[1] &&
			a1 == lengthVector[i].first[0]) return i;
	}
	return -1;
}

// returns index of angle between atoms ( starts at 0; -1 = not found)
int getAngleIndex(OBAtom* a0, OBAtom* a1, OBAtom* a2) {
	for (unsigned int i = 0; i < angleVector.size(); ++i) {
		if (a0 == angleVector[i].first[0] &&
			a1 == angleVector[i].first[1] &&
			a2 == angleVector[i].first[2]) return i;
		if (a0 == angleVector[i].first[2] &&
			a1 == angleVector[i].first[1] &&
			a2 == angleVector[i].first[0]) return i;
	}
	return -1;
}

// returns index of torsion between atoms ( starts at 1; 0 = not found; negative = "wrong way round")
int getTorsionIndex(OBAtom* a0, OBAtom* a1, OBAtom* a2, OBAtom* a3) {
	for (unsigned int i = 0; i < torsionVector.size(); ++i) {
		if (a0 == torsionVector[i].first[0] &&
			a1 == torsionVector[i].first[1] &&
			a2 == torsionVector[i].first[2] &&
			a3 == torsionVector[i].first[3]) return i+1;
		if (a0 == torsionVector[i].first[3] &&
			a1 == torsionVector[i].first[2] &&
			a2 == torsionVector[i].first[1] &&
			a3 == torsionVector[i].first[0]) return -((int)i+1);
	}
	return 0;
}

//attempts to get torsions in right order for internals
// returns serial of torsion starting at 1. If not found
// returns zero. If torsion is "wrong way round", returns
// negative serial
int getFirstTorsionIndexForAtom(OBAtom* a0) {
	unsigned int k = a0->GetIdx();
	for (unsigned int i = 0; i < torsionVector.size(); ++i) {
		if (a0 == torsionVector[i].first[0]) {
			bool ok = true;
			for (unsigned int j = 1; j <= 3; ++j) {
				OBAtom* atPtr = torsionVector[i].first[j];
				if (atPtr->GetIdx() > k) {
					ok = true;
					break;
				}
			}
			if (ok) return i+1;
		}
 		if (a0 == torsionVector[i].first[3]) {
			bool ok = true;
			for (unsigned int j = 0; j < 3; ++j) {
				OBAtom* atPtr = torsionVector[i].first[j];
				if (atPtr->GetIdx() > k) {
					ok = true;
					break;
				}
			}
			if (ok) return -((int)i+1);
		}
	}
	return 0;
}

// try to generate cartesians from length/angle/torsion
// note that internals have implied order
// this is tacky, but the semantics of internals is messy
void generateInternals() {
	internalVector.clear();
    OBInternalCoord *coord;

	if (molPtr->HasNonZeroCoords()) return;
	unsigned int nTors = torsionVector.size();
	if (nTors == 0) return;
	unsigned int nAng = angleVector.size();
	if (nAng == 0) return;
	unsigned int nLen = lengthVector.size();
	if (nLen == 0) return;
	if (nLen + 1 < molPtr->NumAtoms()) {
	  //		if (nLen > 0) cout << "Not enough lengths to generate all internals" << endl;
	}
	if (nAng + 2 < molPtr->NumAtoms()) {
	  //		if (nAng > 0) cout << "Not enough angles to generate all internals" << endl;
	}
	if (nTors + 3 < molPtr->NumAtoms()) {
	  //		if (nTors > 0) cout << "Not enough torsions to generate all internals" << endl;
	}
// require a torsion with first 4 atoms
	OBAtom* at0 = molPtr->GetAtom(1);
	OBAtom* at1 = molPtr->GetAtom(2);
	OBAtom* at2 = molPtr->GetAtom(3);
// first atom
	coord = new OBInternalCoord();
	internalVector.push_back(coord);
// second (distance only)
	coord = new OBInternalCoord();
	coord->_a = at0;
      int a0 = at0->GetIdx();
      int a1 = at1->GetIdx();
      int a2 = at2->GetIdx();
	int idx = getLengthIndex(at0, at1);
	if (idx == -1) {
//	  cmlError("cannot find length: " + a0 + _SLASH + a1);
		return;
	}
	coord->_dst = lengthVector[idx].second;
	internalVector.push_back(coord);
// third (distance and angle only)
	coord = new OBInternalCoord();
	coord->_a = at1;
	coord->_b = at0;
	idx = getLengthIndex(at1, at2);
	if (idx == -1) {
//	  cmlError("cannot find length: " + a1 + _SLASH + a2 + endl);
		return;
	}
	coord->_dst = lengthVector[idx].second;
	idx = getAngleIndex(at0, at1, at2);
	if (idx == -1) {
//	  cmlError("cannot find angle: " + a0 + _SLASH + a1 + _SLASH + at2);
		return;
	}
	coord->_ang = angleVector[idx].second;
	internalVector.push_back(coord);

	// non-infinite loop: already checked that 3 < molPtr->NumAtoms()
	for (unsigned int i = 3; i < molPtr->NumAtoms(); ++i) {
		OBAtom* at0 = molPtr->GetAtom(i+1);
		idx = getFirstTorsionIndexForAtom(at0);
		if (idx == 0) {
		  cmlError("cannot find torsion... ");
			return;
		}
		int iTor = (idx > 0) ? idx-1 : -idx - 1;
		at0 = torsionVector[iTor].first[0];
		OBAtom* at1 = torsionVector[iTor].first[1];
		OBAtom* at2 = torsionVector[iTor].first[2];
		OBAtom* at3 = torsionVector[iTor].first[3];
      a0 = at0->GetIdx();
      a1 = at1->GetIdx();
      a2 = at2->GetIdx();
      int a3 = at3->GetIdx();
// if torsion is "wrong way round", swap atoms
		if (idx < 0) {
			OBAtom* temp = at0;
			at0 = at3;
			at3 = temp;
			temp = at1;
			at1 = at2;
			at2 = temp;
		}
// distance, and angle
		coord = new OBInternalCoord();
		coord->_a = at1;
		coord->_b = at2;
		coord->_c = at3;
		idx = getLengthIndex(at2, at3);
		if (idx == -1) {
//		  cmlError(("cannot find length: " + a2) + (_SLASH + a3));
			return;
		}
		coord->_dst = lengthVector[idx].second;
		int idx = getAngleIndex(at1, at2, at3);
		if (idx == -1) {
//		  cmlError("cannot find angle: " + a1 + (_SLASH + a2) + (_SLASH + a3));
			return;
		}
		coord->_ang = angleVector[idx].second;
		coord->_tor = torsionVector[iTor].second;
		internalVector.push_back(coord);
	}
	for (internalVector_t::size_type i = 0; i < internalVector.size(); ++i) {
		OBInternalCoord* coord = internalVector[i];
		int aa = (coord->_a != 0) ? coord->_a->GetIdx() : 0;
		int bb = (coord->_b != 0) ? coord->_b->GetIdx() : 0;
		int cc = (coord->_c != 0) ? coord->_c->GetIdx() : 0;
		//		cout << "a" << cc << _SPACE;
		//		cout << "a" << bb << _COLON;
// 		cout << "a" << aa << _COLON;
// 		cout << "a" << (i+1) << _COLON;
// 		cout << coord->_dst << _SPACE;
// 		cout << coord->_ang << _SPACE;
// 		cout << coord->_tor << endl;
	}
}

bool endMolecule() {
  //	debugMolecule(cout);
	generateInternals();
    InternalToCartesian(internalVector, *molPtr);

	molPtr->EndModify();

	molPtr->ConnectTheDots();

	if (outputDebug) {
	  //		debug(cout);
	}
	return true; // [ejk] assumed
}

bool WriteMolecule(ostream &ofs) {
    clearMoleculeWorkspace();
    
	if (outputDeclaration) {
		ofs << _LANGLE << _QUERY << X_XML << _SPACE << X_VERSION << _EQUALS <<  _QUOTE << "1.0" << _QUOTE << _QUERY << _RANGLE << endl;
	}
	if (outputDoctype) {
		ofs << X_DOCTYPE << _SPACE << C_MOLECULE << _SPACE << X_SYSTEM << _QUOTE << CML1_NAMESPACE << _QUOTE << _RANGLE << endl;
	}
	if (outputPretty) {
//		cout << "<!-- imagine the XML is pretty printed -->" << endl;
	}
	outputPrefix = "";
	if (outputNamespace) {
		outputPrefix += C_PREFIX;
		outputPrefix +=_COLON;
	}
    
    writeStartTagStart(ofs, C_MOLECULE);
	if (outputNamespace) {
		ofs << _SPACE << X_XMLNS << _COLON << C_PREFIX << _EQUALS << _QUOTE << CML2_NAMESPACE << _QUOTE << endl;
	}
    string title = getNormalizedString(molPtr->GetTitle());
 	writeAttribute(ofs, C_TITLE, title);
// a mechanism is needed for IDs on elements
	writeStartTagEnd(ofs);
    
    ofs << endl;
    if (outputCML2) {
        WriteMetadataList(ofs);
    }

	if (molPtr->HasData(obCommentData)) {
		OBCommentData *cd = (OBCommentData*)molPtr->GetData(obCommentData);
        string nData = getNormalizedString(cd->GetData());
        if (nData.length() > 0) {
            if (outputCML1) {
                writeStartTagStart(ofs, C_STRING);
                writeAttribute(ofs, C_TITLE, "comment");
                writeStartTagEnd(ofs);
                ofs << nData;
                writeEndTag(ofs, C_STRING);
            } else if (outputCML2) {
                writeStartTagStart(ofs, C_SCALAR);
                writeAttribute(ofs, C_DICTREF, "foo:comment");
                writeStartTagEnd(ofs);
                ofs << nData;
                writeEndTag(ofs, C_SCALAR);
            }
        }
	}

	if (outputDebug) debug(ofs);

	WriteAtomArray(ofs);
	WriteBondArray(ofs);

	vector<OBGenericData*>::iterator k;
	vector<OBGenericData*> vdata = molPtr->GetData();
	for (k = vdata.begin();k != vdata.end();++k) {
		if ((*k)->GetDataType() == obPairData) {
            if (outputCML1) {
                writeStartTagStart(ofs, C_STRING);
                writeAttribute(ofs, C_TITLE, (*k)->GetAttribute());
                writeStartTagEnd(ofs);
                ofs << ((OBPairData*)(*k))->GetValue();
                writeEndTag(ofs, C_STRING);
            } else if (outputCML2) {
                writeStartTagStart(ofs, C_SCALAR);
                writeAttribute(ofs, C_DICTREF, (*k)->GetAttribute());
                writeStartTagEnd(ofs);
                ofs << ((OBPairData*)(*k))->GetValue();
                writeEndTag(ofs, C_SCALAR);
            }
		}
	}

    writeEndTag(ofs, C_MOLECULE);

	return(true);
}

// --------------------<metadata>----------------

bool WriteMetadataList(ostream &ofs) {
    
    writeStartTagStart(ofs, C_METADATALIST);
    writeAttribute(ofs, C_TITLE, "generated automatically from Openbabel");
    writeStartTagEnd(ofs);
    ofs << endl;
    
    writeStartTagStart(ofs, C_METADATA);
    writeAttribute(ofs, C_NAME, DC_CREATOR);
    writeAttribute(ofs, C_CONTENT, "OpenBabel version 1-100.1");
    writeCombinedTagEnd(ofs);
        
    writeStartTagStart(ofs, C_METADATA);
    writeAttribute(ofs, C_NAME, DC_DESCRIPTION);
    writeAttribute(ofs, C_CONTENT, "Conversion of legacy filetype to CML");
    writeCombinedTagEnd(ofs);
        
    writeStartTagStart(ofs, C_METADATA);
    writeAttribute(ofs, C_NAME, DC_IDENTIFIER);
    writeAttribute(ofs, C_CONTENT, "Unknown");
    writeCombinedTagEnd(ofs);
        
    writeStartTagStart(ofs, C_METADATA);
    writeAttribute(ofs, C_NAME, DC_CONTENT);
    string content = cmlType;
    if (inputArray) content += " array";
    writeAttribute(ofs, C_CONTENT, content);
    writeCombinedTagEnd(ofs);
        
    writeStartTagStart(ofs, C_METADATA);
    writeAttribute(ofs, C_NAME, DC_RIGHTS);
    writeAttribute(ofs, C_CONTENT, "unknown");
    writeCombinedTagEnd(ofs);
        
    writeStartTagStart(ofs, C_METADATA);
    writeAttribute(ofs, C_NAME, DC_TYPE);
    writeAttribute(ofs, C_CONTENT, "chemistry");
    writeCombinedTagEnd(ofs);
        
    writeStartTagStart(ofs, C_METADATA);
    writeAttribute(ofs, C_NAME, DC_CONTRIBUTOR);
    writeAttribute(ofs, C_CONTENT, "unknown");
    writeCombinedTagEnd(ofs);
        
    writeStartTagStart(ofs, C_METADATA);
    writeAttribute(ofs, C_NAME, DC_CREATOR);
    writeAttribute(ofs, C_CONTENT, "Openbabel V1-100.1");
    writeCombinedTagEnd(ofs);
        
    writeStartTagStart(ofs, C_METADATA);
    writeAttribute(ofs, C_NAME, DC_DATE);
    string time;
    getTimestr(time);
    writeAttribute(ofs, C_CONTENT, time);
    writeCombinedTagEnd(ofs);
        
    writeStartTagStart(ofs, C_METADATA);
    writeAttribute(ofs, C_NAME, CMLM_STRUCTURE);
    writeAttribute(ofs, C_CONTENT, "yes");
    writeCombinedTagEnd(ofs);
        
    writeEndTag(ofs, C_METADATALIST);
}

bool getTimestr(string& s) {
    time_t akttime;                              /* Systemtime                        */
    char timestr[TIME_STR_SIZE + 1] = "";        /* Timestring                        */
    size_t time_res;                             /* Result of strftime                */
    char *log_name;                              /* Pointer to buffer with login name */
    
    /* ---- Get the system-time ---- */
    akttime = time((time_t *) NULL);
    time_res = strftime(timestr,
                     TIME_STR_SIZE,
                     "%a %b %d %H:%M:%S %Z %Y",
                     localtime((time_t *) &akttime)
                    );
    s = getNormalizedString(timestr);                    
    return true;;                    
}


// --------------------<reaction>----------------

bool startReaction(vector <pair<string,string> > &atts) {
	vector <string> badAtts = getUnknownAttributes(REACTION_ATTRIBUTE_VECTOR, atts);
	if (badAtts.size() > 0) {
	  cmlError("unknown attributes on <reaction>: ");
	  //		printVector(badAtts, cerr);
	}
// check other attributes
	for (vector <pair<string,string> >::size_type i = 0; i < atts.size(); ++i) {
		if (atts[i].first == "id") {
	    } else if (atts[i].first == C_TITLE) {
	    } else if (atts[i].first == C_CONVENTION) {
	    } else {
	      //			cout << "IGNORED reaction attribute: " << atts[i].first << endl;
		}
	}
	return true; // [ejk] assumed
}

bool endReaction() {
	return true; // [ejk] assumed
}

bool WriteReaction(ostream &ofs) {
    writeStartTagStart(ofs, C_REACTION);
    writeStartTagEnd(ofs);
    writeEndTag(ofs, C_REACTION);
	return true; // [ejk] assumed
}

// --------------------<scalar>----------------

bool startScalar(vector <pair<string,string> > &atts) {
	vector <string> badAtts = getUnknownAttributes(SCALAR_ATTRIBUTE_VECTOR, atts);
	if (badAtts.size() > 0) {
	  cmlError("unknown attributes on <scalar>: ");
	  //		printVector(badAtts, cerr);
	}
// check other attributes
	for (vector <pair<string,string> >::size_type i = 0; i < atts.size(); ++i) {
		if (atts[i].first == "id") {
	    } else if (atts[i].first == C_TITLE) {
	    } else if (atts[i].first == C_CONVENTION) {
	    } else if (atts[i].first == C_DATATYPE) {
			scalarDataType = atts[i].second;
	    } else if (atts[i].first == C_UNITS) {
			scalarUnits = atts[i].second;
	    } else {
	      //			cout << "IGNORED scalar attribute: " << atts[i].first << endl;
		}
	}
	return true; // [ejk] assumed
}

bool endScalar() {
	string title = getAttribute(currentAtts, C_TITLE);
// only place scalar is required in CML2
	if (parent == "crystal") {
		double f = atof((char*)pcdata.c_str());
		if (title == "a") cellParam[0] = f;
		if (title == "b") cellParam[1] = f;
		if (title == "c") cellParam[2] = f;
		if (title == "alpha") cellParam[3] = f;
		if (title == "beta") cellParam[4] = f;
		if (title == "gamma") cellParam[5] = f;
	}
	return true; // [ejk] assumed
}

// not yet finished
bool WriteScalar(ostream &ofs) {
    writeStartTagStart(ofs, C_SCALAR);
    writeStartTagEnd(ofs);
    writeEndTag(ofs, C_SCALAR);
	return true; // [ejk] assumed
}

// --------------------<sequence>----------------

bool startSequence(vector <pair<string,string> > &atts) {
	vector <string> badAtts = getUnknownAttributes(SEQUENCE_ATTRIBUTE_VECTOR, atts);
	if (badAtts.size() > 0) {
	  cmlError("unknown attributes on <sequence>: ");
	  //		printVector(badAtts, cerr);
	}
// check other attributes
	for (vector <pair<string,string> >::size_type i = 0; i < atts.size(); ++i) {
		if (atts[i].first == "id") {
	    } else if (atts[i].first == C_TITLE) {
	    } else if (atts[i].first == C_CONVENTION) {
	    } else {
	      //			cout << "IGNORED sequence attribute: " << atts[i].first << endl;
		}
	}
	return true; // [ejk] assumed
}

bool endSequence() {
	return true; // [ejk] assumed
}

// not yet finished
bool WriteSequence(ostream &ofs) {
    writeStartTagStart(ofs, C_SEQUENCE);
    writeStartTagEnd(ofs);
    writeEndTag(ofs, C_SEQUENCE);
	return true; // [ejk] assumed
}

// -------------------------<stereo>-------------------

bool startStereo(vector <pair<string,string> > &atts) {
	atomRefs4 = getAttribute(atts, C_ATOMREFS4);
	return true; // [ejk] assumed
}

bool endStereo(vector <pair<string,string> > &atts) {
	pair <vector<OBAtom*>, string> st;
	vector <OBAtom*> atomRef;
	getAtomRefs(4, atomRef, atomRefs4);
	if (atomRef.size() != 4) {
	  cmlError("atomRefs4 must referemce 4 atoms");
		return false;
	}
	for (unsigned int i = 0; i < 4; ++i) st.first.push_back(atomRef[i]);
	setCMLType(C_CML2);
	st.second = pcdata;
	stereoSVector.push_back(st);
	return true; // [ejk] assumed
}

bool WriteStereo(ostream &ofs) {
    writeStartTagStart(ofs, C_STEREO);
    writeStartTagEnd(ofs);
    writeEndTag(ofs, C_STEREO);
	return true; // [ejk] assumed
}


// -------------------------<symmetry>-------------------

bool startSymmetry(vector <pair<string,string> > &atts) {
	vector <string> badAtts = getUnknownAttributes(SYMMETRY_ATTRIBUTE_VECTOR, atts);
	if (badAtts.size() > 0) {
	  cmlError("unknown attributes on <symmetry>: ");
	  //		printVector(badAtts, cerr);
	}
	spacegroup = getAttribute(atts, C_SPACEGROUP);
	pointgroup = getAttribute(atts, C_POINTGROUP);
// check other attributes
	for (vector <pair<string,string> >::size_type i = 0; i < atts.size(); ++i) {
		if (atts[i].first == "id") {
	    } else if (atts[i].first == C_TITLE) {
	    } else if (atts[i].first == C_CONVENTION) {
	    } else {
	      //			cout << "IGNORED symmetry attribute: " << atts[i].first << endl;
		}
	}
	return true; // [ejk] assumed
}

bool endSymmetry(vector <pair<string,string> > &atts) {
	return true; // [ejk] assumed
}

bool WriteSymmetry(ostream &ofs) {
    writeStartTagStart(ofs, C_SYMMETRY);
    writeAttribute(ofs, C_SPACEGROUP, spacegroup);
    writeAttribute(ofs, C_POINTGROUP, pointgroup);
    writeStartTagEnd(ofs);
    writeEndTag(ofs, C_SYMMETRY);
	return true; // [ejk] assumed
}


// -------------------------<string>-------------------
bool addString() {
	string title = getAttribute(currentAtts, C_TITLE);
	if (title != _EMPTY) {
		OBPairData *dp = new OBPairData;
		dp->SetAttribute(title);
		dp->SetValue(pcdata);
		molPtr->SetData(dp);
	}
	return true; // [ejk] assumed
}

// --------------------<torsion>/<torsionArray>----------------

bool startTorsion(vector <pair<string,string> > &atts) {
	vector <string> badAtts = getUnknownAttributes(TORSION_ATTRIBUTE_VECTOR, atts);
	if (badAtts.size() > 0) {
	  cmlError("unknown attributes on <torsion>: ");
	  //		printVector(badAtts, cerr);
	}
	torsionUnits = "degrees";
	atomRefs4Vector.clear();
// check other attributes
	for (vector <pair<string,string> >::size_type i = 0; i < atts.size(); ++i) {
		if (atts[i].first == "id") {
	    } else if (atts[i].first == C_TITLE) {
	    } else if (atts[i].first == C_CONVENTION) {
	    } else if (atts[i].first == C_ATOMREFS) {
			setCMLType(C_CML1);
			getAtomRefs(4, atomRefs4Vector, atts[i].second);
	    } else if (atts[i].first == C_ATOMREFS4) {
			setCMLType(C_CML2);
			getAtomRefs(4, atomRefs4Vector, atts[i].second);
	    } else if (atts[i].first == C_UNITS) {
			torsionUnits = atts[i].second;
	    } else {
	      //			cout << "IGNORED torsion attribute: " << atts[i].first << endl;
		}
	}
	return true; // [ejk] assumed
}

bool endTorsion() {
	pair <vector<OBAtom*>, double> torsion;
	if (atomRefs4Vector.size() != 4) {
	  cmlError("must have defined 4 atoms for torsion");
	}
	for (int i = 0; i < 4; ++i) {
		torsion.first.push_back(atomRefs4Vector[i]);
	}
	torsion.second = atof((char*)pcdata.c_str());
	torsionVector.push_back(torsion);
	return true; // [ejk] assumed
}

// not yet finished
bool WriteTorsion(ostream &ofs, pair <vector<OBAtom*>, double> torsion) {
    writeStartTagStart(ofs, C_TORSION);
    string atomRefs4 = "a";
    atomRefs4 += torsion.first[0]->GetIdx();
    atomRefs4 += " a";
    atomRefs4 += torsion.first[1]->GetIdx();
    atomRefs4 += " a";
    atomRefs4 += torsion.first[2]->GetIdx();
    atomRefs4 += " a";
    atomRefs4 += torsion.first[3]->GetIdx();
    writeAttribute(ofs, C_ATOMREFS4, atomRefs4);
    writeStartTagEnd(ofs);
	ofs << torsion.second;
    writeEndTag(ofs, C_TORSION);

	return true; // [ejk] assumed
}


/**------------------Babel INPUT OUTPUT --------------------------*/


  // Clear out all the global variable vectors after read/write
void CleanUp()
{
  atomIdVector.clear();
  idVector.clear();
  elementTypeVector.clear();
  atomicNumVector.clear();
  formalChargeVector.clear();
  hydrogenCountVector.clear();
  x2Vector.clear();
  y2Vector.clear();
  x3Vector.clear();
  y3Vector.clear();
  z3Vector.clear();
  atomRefs2Vector.clear();
  atomRefs3Vector.clear();
  atomRefs4Vector.clear();
  atomRef1Vector.clear();
  atomRef2Vector.clear();
  orderVector.clear();
  stereoVector.clear();
  rotTransVector.clear();
  rotVector.clear();
  angleVector.clear();
  lengthVector.clear();
  torsionVector.clear();
  atomParityVector.clear();
  stereoSVector.clear();
  internalVector.clear();
}

bool ReadCML(istream &ifs,OBMol &mol, const char *title) {
	molPtr = &mol;
	ReadXML(ifs);

    outputDebug = false;    
    if (outputDebug) {
        debug(cout);
    }
    
	CleanUp();
	return true;
}

// output routines

bool WriteCML(ostream &ofs,OBMol &mol,const char *dim,const char* xmlOptions)
{
	ofsPtr = &ofs;
	outputCML1 = true;
	outputCML2 = false;
	outputDoctype = false;
	outputPretty = false;
	outputNamespace = false;
	outputPrefix = "";
	outputArray = false;
	outputDebug = false;

	if (xmlOptions != 0) {
		string xo(xmlOptions);
		if (xo.find("2") < xo.length()) {
			outputCML2 = true;
			outputCML1 = false;
		}
		if (xo.find("d") < xo.length()) outputDoctype = true;
		if (xo.find("p") < xo.length()) outputPretty = true;
		if (xo.find("n") < xo.length()) outputNamespace = true;
		if (xo.find("c") < xo.length()) outputPrefix = "cml";
		if (xo.find("a") < xo.length()) outputArray = true;
		if (xo.find("g") < xo.length()) outputDebug = true;
		if (xo.find("v") < xo.length()) outputDeclaration = true;
	}

	molPtr = &mol;
	dimension = dim;
	WriteMolecule(ofs);
	CleanUp();
	return true; // [ejk] assumed
}

void debug(ostream &ofs) {
	ofs << "<table>" << endl;
	ofs << "<tr><td>Flags</td><td>" << molPtr->GetFlags() << "</td></tr>" << endl;
	ofs << "<tr><td>Title</td><td>" << molPtr->GetTitle() << "</td></tr>" << endl;
	ofs << "<tr><td>InputType</td><td>" << molPtr->GetInputType() << "</td></tr>" << endl;
	ofs << "<tr><td>OutputType</td><td>" << molPtr->GetOutputType() << "</td></tr>" << endl;
	ofs << "<tr><td>NumAtoms</td><td>" << molPtr->NumAtoms() << "</td></tr>" << endl;
	ofs << "<tr><td>NumBonds</td><td>" << molPtr->NumBonds() << "</td></tr>" << endl;
	ofs << "<tr><td>NumHvyAtoms</td><td>" << molPtr->NumHvyAtoms() << "</td></tr>" << endl;
	ofs << "<tr><td>NumResidues</td><td>" << molPtr->NumResidues() << "</td></tr>" << endl;
	ofs << "<tr><td>NumRotors</td><td>" << molPtr->NumRotors() << "</td></tr>" << endl;
	ofs << "<tr><td>Energy</td><td>" << molPtr->GetEnergy() << "</td></tr>" << endl;
	ofs << "<tr><td>MolWt</td><td>" << molPtr->GetMolWt() << "</td></tr>" << endl;
	ofs << "<tr><td>IsCompressed</td><td>" << molPtr->IsCompressed() << "</td></tr>" << endl;
	ofs << "<tr><td>2D</td><td>" << molPtr->Has2D() << "</td></tr>" << endl;
	ofs << "<tr><td>3D</td><td>" << molPtr->Has3D() << "</td></tr>" << endl;
	ofs << "<tr><td>NonZeroCoords</td><td>" << molPtr->HasNonZeroCoords() << "</td></tr>" << endl;
	ofs << "<tr><td>AromaticPerceived</td><td>" << molPtr->HasAromaticPerceived() << "</td></tr>" << endl;
	ofs << "<tr><td>SSSRPerceived</td><td>" << molPtr->HasSSSRPerceived() << "</td></tr>" << endl;
	ofs << "<tr><td>RingAtomsAndBondsPerceived</td><td>" << molPtr->HasRingAtomsAndBondsPerceived() << "</td></tr>" << endl;
	ofs << "<tr><td>AtomTypesPerceived</td><td>" << molPtr->HasAtomTypesPerceived() << "</td></tr>" << endl;
	ofs << "<tr><td>ChiralityPerceived</td><td>" << molPtr->HasChiralityPerceived() << "</td></tr>" << endl;
	ofs << "<tr><td>PartialChargesPerceived</td><td>" << molPtr->HasPartialChargesPerceived() << "</td></tr>" << endl;
	ofs << "<tr><td>HybridizationPerceived</td><td>" << molPtr->HasHybridizationPerceived() << "</td></tr>" << endl;
	ofs << "<tr><td>ImplicitValencePerceived</td><td>" << molPtr->HasImplicitValencePerceived() << "</td></tr>" << endl;
	ofs << "<tr><td>KekulePerceived</td><td>" << molPtr->HasKekulePerceived() << "</td></tr>" << endl;
	ofs << "<tr><td>ClosureBondsPerceived</td><td>" << molPtr->HasClosureBondsPerceived() << "</td></tr>" << endl;
	ofs << "<tr><td>ChainsPerceived</td><td>" << molPtr->HasChainsPerceived() << "</td></tr>" << endl;
	ofs << "<tr><td>HydrogensAdded</td><td>" << molPtr->HasHydrogensAdded() << "</td></tr>" << endl;
	ofs << "<tr><td>AromaticCorrected</td><td>" << molPtr->HasAromaticCorrected() << "</td></tr>" << endl;
	ofs << "<tr><td>IsCorrectedForPH</td><td>" << molPtr->IsCorrectedForPH() << "</td></tr>" << endl;
	ofs << "<tr><td>IsChiral</td><td>" << molPtr->IsChiral() << "</td></tr>" << endl;
	ofs << "</table>" << endl;

	ofs << "<h2>Atoms</h2>" << endl;
	ofs << "<table>" << endl;
	ofs << "<tr>" << endl;
    	ofs << "<th>Num</th>" << endl;
    	ofs << "<th>FChg</th>" << endl;
    	ofs << "<th>AtNo</th>" << endl;
    	ofs << "<th>Idx</th>" << endl;
    	ofs << "<th>CrdIdx</th>" << endl;
    	ofs << "<th>CIdx</th>" << endl;
    	ofs << "<th>Val</th>" << endl;
    	ofs << "<th>Hyb</th>" << endl;
    	ofs << "<th>ImpVal</th>" << endl;
    	ofs << "<th>HvVal</th>" << endl;
    	ofs << "<th>HetVal</th>" << endl;
    	ofs << "<th>Typ</th>" << endl;

    	ofs << "<th>X</th>" << endl;
    	ofs << "<th>Y</th>" << endl;
    	ofs << "<th>Z</th>" << endl;

    	ofs << "<th>PChg</th>" << endl;

    	ofs << "<th>FreeOx</th>" << endl;
    	ofs << "<th>ImpHyd</th>" << endl;
    	ofs << "<th>ExpHyd</th>" << endl;
    	ofs << "<th>MembRing</th>" << endl;
    	ofs << "<th>BOSum</th>" << endl;
    	ofs << "<th>KBOSum</th>" << endl;

    	ofs << "<th>Res</th>" << endl;
    	ofs << "<th>H</th>" << endl;
    	ofs << "<th>C</th>" << endl;
    	ofs << "<th>N</th>" << endl;
    	ofs << "<th>O</th>" << endl;
    	ofs << "<th>S</th>" << endl;
    	ofs << "<th>P</th>" << endl;
    	ofs << "<th>Ar</th>" << endl;
    	ofs << "<th>Ring</th>" << endl;
    	ofs << "<th>Het</th>" << endl;
    	ofs << "<th>C=O</th>" << endl;
    	ofs << "<th>P=O</th>" << endl;
    	ofs << "<th>S=O</th>" << endl;
    	ofs << "<th>N=O</th>" << endl;
    	ofs << "<th>AmN</th>" << endl;
    	ofs << "<th>PolH</th>" << endl;
    	ofs << "<th>NonPH</th>" << endl;
    	ofs << "<th>AroNOx</th>" << endl;
    	ofs << "<th>Chir</th>" << endl;
    	ofs << "<th>Ax</th>" << endl;
    	ofs << "<th>Clock</th>" << endl;
    	ofs << "<th>AntiCl</th>" << endl;
    	ofs << "<th>ChirSpec</th>" << endl;
    	ofs << "<th>AlBeta</th>" << endl;
    	ofs << "<th>Sing</th>" << endl;
    	ofs << "<th>Doub</th>" << endl;
    	ofs << "<th>Arom</th>" << endl;
    	ofs << "</tr>" << endl;

	for (unsigned int i = 0; i < molPtr->NumAtoms(); ++i) {
		OBAtom* atPtr = molPtr->GetAtom(i+1);
    	ofs << "<tr>" ;
    	ofs << "<td>" << i+1 << "</td>" ;
    	ofs << "<td>" << atPtr->GetFormalCharge() << "</td>" ;
    	ofs << "<td>" << atPtr->GetAtomicNum() << "</td>" ;
    	ofs << "<td>" << atPtr->GetIdx() << "</td>" ;
    	ofs << "<td>" << atPtr->GetCoordinateIdx() << "</td>" ;
    	ofs << "<td>" << atPtr->GetCIdx() << "</td>" ;
    	ofs << "<td>" << atPtr->GetValence() << "</td>" ;
    	ofs << "<td>" << atPtr->GetHyb() << "</td>" ;
    	ofs << "<td>" << atPtr->GetImplicitValence() << "</td>" ;
    	ofs << "<td>" << atPtr->GetHvyValence() << "</td>" ;
    	ofs << "<td>" << atPtr->GetHeteroValence() << "</td>" ;
    	ofs << "<td>" << atPtr->GetType() << "</td>" ;

    	ofs << "<td>" << atPtr->GetX() << "</td>" ;
    	ofs << "<td>" << atPtr->GetY() << "</td>" ;
    	ofs << "<td>" << atPtr->GetZ() << "</td>" ;

    	ofs << "<td>" << atPtr->GetPartialCharge() << "</td>" ;

    	ofs << "<td>" << atPtr->CountFreeOxygens() << "</td>" ;
    	ofs << "<td>" << atPtr->ImplicitHydrogenCount() << "</td>" ;
    	ofs << "<td>" << atPtr->ExplicitHydrogenCount() << "</td>" ;
    	ofs << "<td>" << atPtr->MemberOfRingCount() << "</td>" ;
    	ofs << "<td>" << atPtr->BOSum() << "</td>" ;
    	ofs << "<td>" << atPtr->KBOSum() << "</td>" ;

    	ofs << "<td>" << atPtr->HasResidue() << "</td>" ;
    	ofs << "<td>" << atPtr->IsHydrogen() << "</td>" ;
    	ofs << "<td>" << atPtr->IsCarbon() << "</td>" ;
    	ofs << "<td>" << atPtr->IsNitrogen() << "</td>" ;
    	ofs << "<td>" << atPtr->IsOxygen() << "</td>" ;
    	ofs << "<td>" << atPtr->IsSulfur() << "</td>" ;
    	ofs << "<td>" << atPtr->IsPhosphorus() << "</td>" ;
    	ofs << "<td>" << atPtr->IsAromatic() << "</td>" ;
    	ofs << "<td>" << atPtr->IsInRing() << "</td>" ;
    	ofs << "<td>" << atPtr->IsHeteroatom() << "</td>" ;
    	ofs << "<td>" << atPtr->IsCarboxylOxygen() << "</td>" ;
    	ofs << "<td>" << atPtr->IsPhosphateOxygen() << "</td>" ;
    	ofs << "<td>" << atPtr->IsSulfateOxygen() << "</td>" ;
    	ofs << "<td>" << atPtr->IsNitroOxygen() << "</td>" ;
    	ofs << "<td>" << atPtr->IsAmideNitrogen() << "</td>" ;
    	ofs << "<td>" << atPtr->IsPolarHydrogen() << "</td>" ;
    	ofs << "<td>" << atPtr->IsNonPolarHydrogen() << "</td>" ;
    	ofs << "<td>" << atPtr->IsAromaticNOxide() << "</td>" ;
    	ofs << "<td>" << atPtr->IsChiral() << "</td>" ;
    	ofs << "<td>" << atPtr->IsAxial() << "</td>" ;
    	ofs << "<td>" << atPtr->IsClockwise() << "</td>" ;
    	ofs << "<td>" << atPtr->IsAntiClockwise() << "</td>" ;
    	ofs << "<td>" << atPtr->HasChiralitySpecified() << "</td>" ;
    	ofs << "<td>" << atPtr->HasAlphaBetaUnsat() << "</td>" ;
    	ofs << "<td>" << atPtr->HasSingleBond() << "</td>" ;
    	ofs << "<td>" << atPtr->HasDoubleBond() << "</td>" ;
    	ofs << "<td>" << atPtr->HasAromaticBond() << "</td>" ;
    	ofs << "</tr>" << endl;
	}
	ofs << "</table>" << endl;

	ofs << "<h2>Bonds</h2>" << endl;
	ofs << "<table>" << endl;
	ofs << "<tr>" << endl;
	ofs << "<th>Num</th>" << endl;
	ofs << "<th>BO</th>" << endl;
	ofs << "<th>BoOrd</th>" << endl;
	ofs << "<th>Flags</th>" << endl;
	ofs << "<th>BegAt</th>" << endl;
	ofs << "<th>EndAt</th>" << endl;
	ofs << "<th>EqLen</th>" << endl;
	ofs << "<th>Len</th>" << endl;
	ofs << "<th>Arom</th>" << endl;
	ofs << "<th>InRg</th>" << endl;
	ofs << "<th>Rot</th>" << endl;
	ofs << "<th>Ami</th>" << endl;
	ofs << "<th>PriAm</th>" << endl;
	ofs << "<th>SecAm</th>" << endl;
	ofs << "<th>Est</th>" << endl;
	ofs << "<th>Carb</th>" << endl;
	ofs << "<th>Sin</th>" << endl;
	ofs << "<th>Dou</th>" << endl;
	ofs << "<th>KSing</th>" << endl;
	ofs << "<th>KDoub</th>" << endl;
	ofs << "<th>KTrip</th>" << endl;
	ofs << "<th>Clos</th>" << endl;
	ofs << "<th>Up</th>" << endl;
	ofs << "<th>Down</th>" << endl;
	ofs << "<th>Wedg</th>" << endl;
	ofs << "<th>Hash</th>" << endl;
	ofs << "</tr>" << endl;

	for (unsigned int i = 0; i < molPtr->NumBonds(); ++i) {
		OBBond* boPtr = molPtr->GetBond(i);
		ofs << "<tr>" << endl;
		ofs << "<td>" << i+1 << "</td>" << endl;
		ofs << "<td>" << boPtr->GetBO() << "</td>" << endl;
		ofs << "<td>" << boPtr->GetBondOrder() << "</td>" << endl;
		ofs << "<td>" << boPtr->GetFlags() << "</td>" << endl;
		ofs << "<td>" << boPtr->GetBeginAtomIdx() << "</td>" << endl;
		ofs << "<td>" << boPtr->GetEndAtomIdx() << "</td>" << endl;
		ofs << "<td>" << boPtr->GetEquibLength() << "</td>" << endl;
		ofs << "<td>" << boPtr->GetLength() << "</td>" << endl;
		ofs << "<td>" << boPtr->IsAromatic() << "</td>" << endl;
		ofs << "<td>" << boPtr->IsInRing() << "</td>" << endl;
		ofs << "<td>" << boPtr->IsRotor() << "</td>" << endl;
		ofs << "<td>" << boPtr->IsAmide() << "</td>" << endl;
		ofs << "<td>" << boPtr->IsPrimaryAmide() << "</td>" << endl;
		ofs << "<td>" << boPtr->IsSecondaryAmide() << "</td>" << endl;
		ofs << "<td>" << boPtr->IsEster() << "</td>" << endl;
		ofs << "<td>" << boPtr->IsCarbonyl() << "</td>" << endl;
		ofs << "<td>" << boPtr->IsSingle() << "</td>" << endl;
		ofs << "<td>" << boPtr->IsDouble() << "</td>" << endl;
		ofs << "<td>" << boPtr->IsKSingle() << "</td>" << endl;
		ofs << "<td>" << boPtr->IsKDouble() << "</td>" << endl;
		ofs << "<td>" << boPtr->IsKTriple() << "</td>" << endl;
		ofs << "<td>" << boPtr->IsClosure() << "</td>" << endl;
		ofs << "<td>" << boPtr->IsUp() << "</td>" << endl;
		ofs << "<td>" << boPtr->IsDown() << "</td>" << endl;
		ofs << "<td>" << boPtr->IsWedge() << "</td>" << endl;
		ofs << "<td>" << boPtr->IsHash() << "</td>" << endl;
		ofs << "</tr>" << endl;
	}
	ofs << "</table>" << endl;
}

}
