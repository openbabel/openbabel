/**********************************************************************
Copyright (C) 2002- Peter Murray-Rust.

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

The following conditions apply to anyone or organization making changes to this code
who has not contacted me.

The code itself is OpenSource. You may make any changes under
the terms of the GPL above.

HOWEVER: ANY alteration to the code means that your system cannot be called
a CML system. CML has been trademarked precisely to protect its implementation.
(The International Union of Crystallography has trademarked their CIF
specification, so this policy is in keeping with accepted practice in
Open Specifications).

YOU MUST therefore specifically announce that your modified code does not process CML,
does not obey CML semantics and does not emit CML. You may not generate XML
which uses the CML namespaces posted on http://www.xml-cml.org and you may not use
the <cml> tag, with or without a namespace declaration.

This notice must be included with any modified software.

Peter Murray-Rust, 2002
*/

#include "mol.h"

using namespace std;
namespace OpenBabel {

// default molecule size
int ATOM_SIZE = 100;
string _EMPTY = "";

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
void processFloatTokens(vector <double> &v, int n, string att);
// split a list of ints
void processIntegerTokens(vector <int> &v, int n, string att);
// split a list of ints
void processStringTokens(vector <int> &v, int n, string att);

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
vector<pair<string,string> > namespaceVector;
// has rootElement been read?
bool readRoot;

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

// ----------------CML elements----------------
const string CML1_ELEMENT_NAMES = "angle atom atomArray bond bondArray coordinate2 coordinate3 crystal electron feature float floatArray floatMatrix integer integerArray link list molecule reaction sequence string stringArray torsion";

const string CML1_NAMESPACE = "http://www.xml-cml.org/dtd/cml_1_0_1.dtd";
const string CML2_NAMESPACE = "http://www.xml-cml.org/schema/cml2/core";
const string STMML_NAMESPACE = "http://www.xml-cml.org/schema/stmml";
string cmlType = "";

vector <string> CML_ELEMENT_VECTOR;
string CML_ELEMENT_NAMES = "amount angle atom atomArray atomParity bond bondArray cml crystal electron formula length molecule stereo substance substanceList symmetry torsion ";

vector <string> STMML_ELEMENT_VECTOR;
string STMML_ELEMENT_NAMES = "action actionList alternative annotation array appinfo definition description dictionary dimension documentation entry enumeration link list matrix metadataList metadata object observation relatedEntry scalar stmml table unit unitList unitType ";


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
vector <string> ATOM_ATTRIBUTE_VECTOR;
string ATOM_ATTRIBUTES = "id title convention dictRef";
string ATOM_BUILTINS =
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
string SCALAR_ATTRIBUTES =  "id title convention dataType size units";
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
vector <string> STRING_ATTRIBUTE_VECTOR;
string STRING_ATTRIBUTES =  "id title count convention builtin";

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
void getAtomRefs(int size, vector <OBAtom*> &atomRef, string atomRefs);

// ----------------babel stuff---------------

// numeric value (babel) from CML order; unknown returns -1
int getBabelBondOrder(string o);
// numeric value (babel) from CML stereo; unknown returns -1
int getBabelBondFlag(string s);
// get serial number of atom (starting at 1); 0 for not found
OBAtom *getAtomPtr(string s);

// debug
void debug(ostream &ofs);

// ----------------variables-----------------

// all atoms need to be unique IDs in CML; this links IDs to AtomPtrs
vector <pair <string, OBAtom*> > atomIdVector;
// current Molecule
OBMol *molPtr;
// current Atom
OBAtom *atomPtr;
// current Bond
OBBond *bondPtr;
// dimensionality of coordinates
char *dimension;
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
int nbonds;
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

bool inputCML1;
bool inputCML2;
bool inputNamespace;
char *inputPrefix;
bool inputArray;

bool outputCML1;
bool outputCML2;
bool outputDoctype;
bool outputPretty;
bool outputNamespace;
char *outputPrefix;
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
vector <double[12]> rotTransVector;
vector <double[9]> rotVector;

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
vector <OBInternalCoord*> internalVector;

/** --------------------- initialization-----------------------*/
void makeAllowedElementLists() {
	tokenize(CML_ELEMENT_VECTOR, CML_ELEMENT_NAMES, " \n");
	tokenize(STMML_ELEMENT_VECTOR, STMML_ELEMENT_NAMES, " \n");
}

void makeAllowedAttributeLists() {
	tokenize(ANGLE_ATTRIBUTE_VECTOR, ANGLE_ATTRIBUTES, " \n");
	string s = ATOM_ATTRIBUTES;
	s.append(" ");
	s.append(ATOM_BUILTINS);
	tokenize(ATOM_ATTRIBUTE_VECTOR, s, " \n");
	tokenize(ATOMARRAY_ATTRIBUTE_VECTOR, ATOMARRAY_ATTRIBUTES, " \n");
	s = BOND_ATTRIBUTES;
	s.append(" ");
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
bool ReadXML(istream &ifs) {
	char buffer[BUFF_SIZE];
	size_t lt;
	size_t rt;

	currentElem = "";
	string token = "";
	bool lookForOpenTag = true;

	makeAllowedElementLists();
	makeAllowedAttributeLists();
	startDocument();
	while (ifs.getline(buffer,BUFF_SIZE)) {
		string buff(buffer);
// omit whitespace lines
		if (trim(buff) == "") continue;
		if (readRoot) {
			cerr << "no nonWhitespace allowed after root element: " << buff << endl;
			break;
		}
// normalize Newlines to " "
		if (token != "") token += " ";
		for (;;) {
			if (lookForOpenTag) {
				lt = buff.find("<");
// not found, more input...
				if (lt > buff.size()) {
					token += buff;
					buff = "";
					break;
				} else {
// found start of tag
					token += buff.substr(0,lt);
// process inter-tag characters
					characters(token);
					buff = buff.substr(lt);
					token = "";
					lookForOpenTag = false;
				}
			} else {
				rt = buff.find(">");
// not found, more input...
				if (rt > buff.size()) {
					token += buff;
					buff = "";
					break;
				} else {
// found end
					string ss = buff.substr(0,rt+1);
					token += ss;
					tag(token);
					buff = buff.substr(rt+1);
					token = "";
					lookForOpenTag = true;
				}
			}
			if (buff == "") break;
		}
	}
	endDocument();
}

// process anything in balanced <...>
void tag(string s) {
	// attributes
	vector <pair <string, string> > atts;

	string name;
	int l = s.length();
	string sl = toLowerCase(s);
// XML declaration
	if (sl.substr(0, 5) == "<?xml") {
		if (s.substr(l-2, 2) == "?>") {
//			cout << "xml declaration" << s << endl;
			string ss = s.substr(5, l-7);
			splitAttributes(ss, atts);
			string standalone = getAttribute(atts, "standalone");
			if (standalone == "no") {
				cerr << "cannot process standalone='no' yet" << endl;
			}
			string version = getAttribute(atts, "version");
			if (version != "1.0") {
				cerr << "XML version must be 1.0" << endl;
			}
			string encoding = toLowerCase(getAttribute(atts, "encoding"));
			if (encoding != "utf-8" && encoding != "") {
				cerr << "Cannot support encoding: " << encoding << endl;
			}
		} else {
			cerr << "Bad XML declaration: " << s << endl;
		}
// DOCTYPE is not processed
	} else if (s.substr(0,9) == "<!DOCTYPE") {
		if (s.find("[") <= s.size()) {
			cout << "cannot process internal subset of DOCTYPE " << s << endl;
		} else {
			cout << "DOCTYPE info ignored" << endl;
		}
// comments are ignored
	} else if (s.substr(0,4) == "<!--") {
		if (s.substr(l-3, 3) == "-->") {
//			cout << "Comment ignored: " << s << endl;
		} else {
			cerr << "Bad comment: " << s << endl;
		}
// Processing instructions
	} else if (s.substr(0,2) == "<?") {
		if (s.substr(l-2, 2) == "?>") {
			s = s.substr(2, l-4);
			int idx = s.find(" ");
			string target = (idx < s.size()) ? s.substr(0, idx) : s;
			string data = (idx < s.size()) ? trim(s.substr(idx)) : "";
			processingInstruction(target, data);
		} else {
			cerr << "Bad PI: " << s << endl;
		}
// CDATA sections
	} else if (s.substr(0,9) == "<![CDATA[") {
		if (s.substr(l-3, 3) == "]]>") {
			pcdata += s.substr(9, l-12);
		} else {
			cerr << "Bad CDATA: " << s << endl;
		}
// end tag
	} else if (s.substr(1,1) == "/") {
		endElement(s.substr(2, l-3));
// empty tag
	} else if (s.substr(l-2, 1) == "/") {
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
	for (int i = 0; i < s.length(); i++) {
		ii = (int) cc[i];
		if (cc[i] == '&') {
			ss.append("&amp;");
		} else if (cc[i] == '"') {
			ss.append("&quot;");
		} else if (cc[i] == '\'') {
			ss.append("&apos;");
		} else if (cc[i] == '<') {
			ss.append("&lt;");
		} else if (cc[i] == '>') {
			ss.append("&gt;");
// characters above 255
		} else if (ii > 255) {
			cerr << "characters above 255 not supported in CML" << ii <<  endl;
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
			cerr << "non-printing characters not suported: " << (int)cc[i] << endl;
		}
	}
	return ss;
}

string processXMLEntities(string s) {
	string s0(s);
	string ss;
	for (;;) {
		int idx = s.find("&");
		if (idx >= s.length()) {
			ss.append(s);
			break;
		}
		cerr << idx << endl;
		ss.append(s.substr(0, idx));
		s = s.substr(idx+1);
		idx = s.find(";");
		if (idx >= s.length()) {
			cerr << "entity without closing ; in :" << s0 << ":" << endl;
		}
		string e = s.substr(0, idx);
		if (e == "quot") {
			ss.append("\"");
		} else if (e == "apos") {
			ss.append("'");
		} else if (e == "lt") {
			ss.append("<");
		} else if (e == "gt") {
			ss.append(">");
		} else if (e == "amp") {
			ss.append("&");
		} else if (e.substr(0, 1) == "#") {
			int i = atoi((char*)e.substr(1).c_str());
			if (i >= 32 && i < 256 || i == 9 || i==10 || i==13) {
				ss.append(1, (char)i);
			} else {
				cerr << "unsupported character: #" << i << endl;
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
		cerr << "cannot process entity references..." << s << endl;
	}
	string ss = s;
	string name;
	int idx = s.find(" ");
	if (idx > s.size()) {
		name = s;
		s = "";
	} else {
		name = s.substr(0, idx);
		s = trim(s.substr(idx+1));
	}
	splitAttributes(s, atts);
	if (!isXMLName(name)) {
		cerr << "invalid XML name: " << name << endl;
	}
	startElement(name, atts);
	return name;
}

void splitAttributes(string s, vector <pair <string, string> > &atts) {
	pair<string, string> att;

	while (true) {
		int idx = s.find("=");
		if (idx > s.size()) {
			if (trim(s) != "") {
				cerr << "Bad attribute at " << s << endl;
			}
			break;
		}
		att.first = trim(s.substr(0, idx));
		s = trim(s.substr(idx+1));
		if (s.length() < 2) {
			cerr << "Bad attribute value: " << s << endl;
			break;
		}
// quote or apos
		string quoter = s.substr(0, 1);
		if (quoter != "\"" && quoter != "\'") {
			cerr << "Unquoted attribute value: " << s << endl;
			break;
		}
		s = s.substr(1);
		idx = s.find(quoter);
		if (idx > s.size()) {
			cerr << "Unbalanced quotes in attribute value: " << s << endl;
			break;
		}
		att.second = processXMLEntities(s.substr(0, idx));
		atts.push_back(att);
		s = trim(s.substr(idx+1));
		if (trim(s) == "") break;
	}
}

// check attributes against allowed list; result is unknown attributes
vector <string> getUnknownAttributes(vector <string> &allowed, vector <pair <string, string> > &atts) {
	vector <string> badAtts;
	for (int i = 0; i < atts.size(); i++) {
		string attName = atts[i].first;
		if (attName.substr(0, 5) == "xmlns") continue;
		bool ok = false;
		for (int j = 0; j < allowed.size(); j++) {
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
	for (int i = 0; i < v.size(); i++) {
		ofs << v[i] << " ";
	}
}

void noteUnusedElementName(string name, string msg) {
	if (!isInStringVector(unusedElementNameVector, name)) {
		cout << msg << name << endl;
		unusedElementNameVector.push_back(name);
	}
}


bool writeAttribute(ostream&ofs, string name, int value) {
	ofs << " " << name << "=\"" << value << "\"";
}

bool writeAttribute(ostream&ofs, string name, double value) {
	ofs << " " << name << "=\"" << value << "\"";
}

bool writeAttribute(ostream&ofs, string name, string value) {
	value = trim(value);
	if (value != "") {
		string value1 = escapeXMLEntities(value);
		ofs << " " << name << "=\"" << value1 << "\"";
	}
}

bool writeBuiltin(ostream&ofs, string name, int value) {
	ofs << "<integer builtin=\"" <<  name << "\">" << value << "</integer>" << endl;
}

bool writeBuiltin(ostream&ofs, string name, double value) {
	ofs << "<float builtin=\"" <<  name << "\">" << value << "</float>" << endl;
}

bool writeBuiltin(ostream&ofs, string name, string value) {
	value = trim(value);
	if (value != "") {
		value = escapeXMLEntities(value);
		ofs << "<string builtin=\"" <<  name << "\">" << value << "</string>" << endl;
	}
}

bool appendToArray(string &array, int value) {
	if (array != "") array.append(" ");
	char ss[20];
	sprintf(ss, "%i", value);
	string s(ss);
	array.append(trim(ss));
}

bool appendToArray(string &array, double value) {
	if (array != "") array.append(" ");
	char ss[20];
	sprintf(ss, "%f", value);
	string s(ss);
	array.append(trim(ss));
}

bool appendToArray(string &array, string value) {
	value = escapeXMLEntities(value);
	if (array != "") array.append(" ");
	array.append(trim(value));
}

bool writePCDATA(ostream&ofs, string value) {
	ofs << escapeXMLEntities(value);
}

// ------------------------ SAX events -------------------
// SAX-like call back
void startDocument() {
//  cout << "starting CML document; crude XML parser. Assumes well-formed; ignores DTDs and entities" << endl;
  readRoot = false;
}

// SAX-like call back
void endDocument() {
	cout << "read CML document" << endl;
	for (int i = 0; i < namespaceVector.size(); i++) {
		cout << "namespace :" << namespaceVector[i].first << ":" << namespaceVector[i].second << endl;
	}
}

void startElement(string name, vector<pair<string,string> > &atts) {
	processAttributes(atts);
	pair <string, string> nsPair = getNamespacePair(name);
	name = (nsPair.first == "") ? name : name.substr(nsPair.first.length()+1);
	startElement(nsPair.second, name, nsPair.first, atts);
}

pair <string, string> getNamespacePair(string name) {
	pair <string, string> nsPair;
	nsPair.first = "";
	nsPair.second = "";
	int idx = name.find(":");
	if (idx < name.length()) {
		nsPair.first = name.substr(0, idx);
		name = name.substr(idx+1);
	}
	for (int i = 0; i < namespaceVector.size(); i++) {
		if (namespaceVector[i].first == nsPair.first) {
			nsPair.second = namespaceVector[i].second;
			break;
		}
	}
	return nsPair;
}

// records whether this is CMLV1.0 or CML2 (Schema)
void setCMLType(string ct) {
	if (cmlType == "") {
		 cmlType = ct;
	} else if (cmlType != ct) {
		cerr << "Cannot mix CML namespaces" << ct << "/" << cmlType << endl;
	}
}

void startElement(string namespaceURI, string localName, string prefix, vector<pair<string,string> > &atts) {
//	if (prefix != "") cout << "PREF : " << prefix << endl;
//	if (namespaceURI != "") cout << "NSURI : " << namespaceURI << endl;
	if (currentElem != "") elementStack.push_back(currentElem);
	parent = currentElem;
	currentElem = trim(localName);
	currentAtts = atts;
	pcdata = "";
	useBuiltin = false;
    if (false) {
		startMolecule(atts);
	} else if (localName == "atom") {
        startAtom(atts);
	} else if (localName == "atomArray") {
        startAtomArray(atts);
	} else if (localName == "atomParity") {
		setCMLType("CML2");
        startAtomParity(atts);
	} else if (localName == "bond") {
        startBond(atts);
	} else if (localName == "bondArray") {
        startBondArray(atts);
    } else if (localName == "cml") {
		startCML(atts);
    } else if (localName == "crystal") {
		startCrystal(atts);
    } else if (localName == "electron") {
		startElectron(atts);
    } else if (localName == "feature") {
		startFeature(atts);
		setCMLType("CML1");
    } else if (localName == "formula") {
		startFormula(atts);
    } else if (localName == "molecule") {
		startMolecule(atts);
	} else if (
		localName == "coordinate2" ||
		localName == "coordinate3" ||
		localName == "float" ||
		localName == "float" ||
		localName == "integer" ||
	    localName == "string") {
			setCMLType("CML1");
	// delay processing till endElement as we need the pcdata
	} else if (
		localName == "floatMatrix" ||
		localName == "floatArray" ||
		localName == "integerArray" ||
	    localName == "stringArray") {
			setCMLType("CML1");
	// delay processing till endElement as we need the pcdata
	} else if (localName == "length") {
		setCMLType("CML2");
		startLength(atts);
	} else if (localName == "angle") {
		startAngle(atts);
	} else if (localName == "torsion") {
		startTorsion(atts);
	} else if (localName == "scalar") {
		startScalar(atts);
		setCMLType("CML2");
	} else if (localName == "stereo") {
		startStereo(atts);
		setCMLType("CML2");
	} else if (localName == "array") {
		setCMLType("CML2");
	} else if (localName == "matrix") {
		setCMLType("CML2");
// other CML2 elements neglected in babel
	} else if (
		localName == "substance" ||
		localName == "substanceList" ||
		localName == "amount") {
		setCMLType("CML2");
		noteUnusedElementName(localName, "CML2 element not relevant to babel: ");
	} else if (isInStringVector(CML_ELEMENT_VECTOR, localName)) {
		cout << "[debug] CML element not relevant to babel: " << localName << endl;
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
		setCMLType("CML2");
		noteUnusedElementName(localName, "STMML element not relevant to babel: ");
	} else if (isInStringVector(STMML_ELEMENT_VECTOR, localName)) {
		setCMLType("CML2");
		cout << "[debug] STMML element not relevant to babel: " << localName << endl;
	} else {
		cout << "start element ignored: " << localName << endl;
	}
}

void processAttributes(vector<pair<string,string> > &atts) {
	for (int i = 0; i < atts.size(); i++) {
		string name = atts[i].first;
		if (!isXMLName(name)) {
			cerr << "invalid XML name: " << name << endl;
		} else if (name.substr(0,5) == "xmlns") {
			processNamespace(name.substr(5), atts[i].second);
		}
	}
}

void processNamespace(string name, string value) {
	pair <string, string> ns;

	int idx = name.find(":");
	ns.first = (idx < name.size()) ? name.substr(idx) : "";
	ns.second = value;
	bool nsExists = false;
	for (int i = 0; i < namespaceVector.size(); i++) {
		if (ns.first == namespaceVector[i].first) {
			nsExists = true;
			if (namespaceVector[i].second != value) {
				cerr << "redefinition of namespace: " <<
					namespaceVector[i].second << " => " << value << endl;
			}
			break;
		}
	}
	if (!nsExists) {
		namespaceVector.push_back(ns);
		if (ns.second == STMML_NAMESPACE) {
			setCMLType("CML2");
		} else if (ns.second == CML2_NAMESPACE) {
			setCMLType("CML2");
		} else if (ns.second == CML1_NAMESPACE) {
		}
	}
}

void endElement(string name) {
	pair <string, string> nsPair = getNamespacePair(name);
	name = (nsPair.first == "") ? name : name.substr(nsPair.first.length()+1);
	endElement(nsPair.second, name, nsPair.first);
}

void endElement(string namespaceURI, string localName, string prefix) {
	vector <string> strings;

	string name = trim(localName);
	if (name != currentElem) {
		cerr << "unbalanced tags at: " << name << endl;
	}
    if (name == "molecule") {
		endMolecule();
	} else if (name == "atom") {
		endAtom();
	} else if (name == "atomArray") {
		endAtomArray();
	} else if (name == "atomParity") {
		endAtomArray();
	} else if (name == "bond") {
		endBond();
	} else if (name == "bondArray") {
		endBondArray();
	} else if (name == "crystal") {
		endCrystal();
	} else if (name == "electron") {
		endElectron();
	} else if (name == "formula") {
		endFormula();
	} else if (name == "feature") {
//		endFormula();
	} else if (
		name == "coordinate2" ||
		name == "coordinate3") {
		if ("atom" == parent) {
			processAtomBuiltin();
		} else {
			cout << "IGNORED <" << name << "> as not child of <atom>" << endl;
		}
	} else if (
		name == "float" ||
		name == "integer" ||
	    name == "string") {
		if ("atom" == parent) {
			processAtomBuiltin();
		} else if ("bond" == parent) {
			processBondBuiltin();
		} else if ("molecule" == parent) {
			addString();
		} else {
			cout << "IGNORED <" << name << "> as not child of <molecule>, <atom> or <bond>" << endl;
		}
	} else if (localName == "array") {
		setCMLType("CML2");
	} else if (localName == "matrix") {
		setCMLType("CML2");
	} else if (
		localName == "floatMatrix") {
		setCMLType("CML1");
		// no action
	} else if (
		localName == "floatArray" ||
		localName == "integerArray" ||
	    localName == "stringArray") {
		if (parent == "atomArray") {
			processAtomArrayChild();
		} else if (parent == "bondArray") {
			processBondArrayChild();
		}
	} else if (name == "length") {
		endLength();
	} else if (name == "angle") {
		endAngle();
	} else if (name == "torsion") {
		endTorsion();
	} else if (name == "reaction") {
		endReaction();
	} else if (name == "sequence") {
		endSequence();
	} else {
//		cout << "end element ignored: " << name << endl;
	}
// I assume there is a vector<> function which is neater...
    int ns = elementStack.size();
    if (ns > 0) {
	    currentElem = elementStack[ns-1];
	    parent = (ns <= 1) ? "" : elementStack[ns-2];
	    elementStack.pop_back();
	} else {
	}
	if (ns == 0) {
		readRoot = true;
	}
	pcdata = "";
}

void characters(string s) {
	pcdata = processXMLEntities(s);
}

void processingInstruction(string target, string data) {
	cout << "PI: " << target << " " << data << endl;
}

void skippedEntity(string name) {
	cout << "skipped entity: " << name << endl;
}

// gets attribute of given name; "" if not found
string getAttribute(vector <pair<string, string> > &atts, string name) {
	string s;
    for (int i = 0; i < atts.size(); i++) {
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
	while (c = *str++) {
		if (c >= '0' && c <= '9') {
		} else if (c >= 'a' && c <= 'z') {
		} else if (c >= 'A' && c <= 'Z') {
		} else if (c == '_' || c == ':' || c == '-' || c == '.') {
		} else {
			ok = false;
		}
	}
	if (!ok) {
		cerr << "invalid XML name: " << n << endl;
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
	for (i = 0;i < s.size();i++) {
	    ss[i] = tolower(s[i]);
	}
	return ss;
}

string toUpperCase(string s) {
	string ss(s);
	unsigned int i;
	for (i = 0;i < s.size();i++) {
	    ss[i] = toupper(s[i]);
	}
	return ss;
}

bool isInStringVector(vector <string> v, string s) {
	for (int i = 0; i < v.size(); i++) {
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
		cerr << "must give value for builtin" << endl;
		pcdata = "";
	} else {
		pcdata = strings[0];
	}
}

// get pointer to atom
OBAtom *getAtomPtr(string s) {
// this is crude... I expect vector<> has a map utility
	for (int j = 0; j < atomIdVector.size(); j++) {
		if (s == atomIdVector[j].first) {
			return atomIdVector[j].second;
		}
	}
	return 0;
}

// process atomRefs2="a1 a2", etc.
// add results to vector
void getAtomRefs(int size, vector <OBAtom*> &v, string atomRefString) {
	vector <string> sv;
	atomRefString += " ";
	tokenize(sv, atomRefString, " \n");
	if (sv.size() != size) {
		cerr << "unexpected size for atomRefs attribute: " << sv.size() << "/" << size << endl;
		return;
	}
	for (int i = 0; i < size; i++) {
		OBAtom* atPtr = getAtomPtr(sv[i]);
		if (atPtr == 0) {
			cerr << "cannot find atom: " << sv[i] << endl;
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
		cerr << "unknown attributes on <cml>: ";
		printVector(badAtts, cerr);
		cerr << endl;
	}

}

bool endCML() {
}
// --------------------<angle>/<angleArray>----------------

bool startAngle(vector <pair<string,string> > &atts) {
	vector <string> badAtts = getUnknownAttributes(ANGLE_ATTRIBUTE_VECTOR, atts);
	if (badAtts.size() > 0) {
		cerr << "unknown attributes on <angle>: ";
		printVector(badAtts, cerr);
		cerr << endl;
	}
	angleUnits = "degrees";
	atomRefs3Vector.clear();
// check other attributes
	for (int i = 0; i < atts.size(); i++) {
		if (false) {
	    } else if (atts[i].first == "id") {
	    } else if (atts[i].first == "title") {
	    } else if (atts[i].first == "convention") {
	    } else if (atts[i].first == "atomRefs") {
			setCMLType("CML1");
			getAtomRefs(3, atomRefs3Vector, atts[i].second);
	    } else if (atts[i].first == "atomRefs3") {
			setCMLType("CML2");
			getAtomRefs(3, atomRefs3Vector, atts[i].second);
	    } else if (atts[i].first == "units") {
			angleUnits = atts[i].second;
	    } else {
			cout << "IGNORED angle attribute: " << atts[i].first << endl;
		}
	}
}

bool endAngle() {
	pair <vector<OBAtom *>, double> angle;
	if (atomRefs3Vector.size() != 3) {
		cerr << "must have defined 3 atoms for angle" << endl;
	}
	for (int i = 0; i < 3; i++) {
		angle.first.push_back(atomRefs3Vector[i]);
	}
	angle.second = atof((char*)pcdata.c_str());
	angleVector.push_back(angle);
}

// not yet finished
bool WriteAngle(ostream &ofs, pair <vector<OBAtom*>, double> angle) {
	ofs << "<angle";
	ofs << " atomRefs3=\"a" << angle.first[0]->GetIdx() << " a" << angle.first[1]->GetIdx() << " a" << angle.first[2]->GetIdx() << "\">";
	ofs << angle.second;
	ofs << "</angle>" << endl;
}

// --------------------<atom>/<atomArray>----------------

bool startAtom(vector <pair<string,string> > &atts) {
	vector <string> badAtts = getUnknownAttributes(ATOM_ATTRIBUTE_VECTOR, atts);
	if (badAtts.size() > 0) {
		cerr << "unknown attributes on <atom>: ";
		printVector(badAtts, cerr);
		cerr << endl;
	}
	currentX = currentY = currentZ = 0;
	formalCharge = 0;

	atomicNum = etab.GetAtomicNum((char*)getAttribute(atts, "elementType").c_str());
	atomId = getAttribute(atts, "id");
	formalCharge = atoi(getAttribute(atts, "formalCharge").c_str());
	string x2 = getAttribute(atts, "x2");
	string y2 = getAttribute(atts, "y2");
	string x3 = getAttribute(atts, "x3");
	string y3 = getAttribute(atts, "y3");
	string z3 = getAttribute(atts, "z3");
	if (x3 != "") {
		currentX = atof(x3.c_str());
		setCMLType("CML2");
	} else if (x2 != "") {
		currentX = atof(x2.c_str());
		setCMLType("CML2");
	}
	if (y3 != "") {
		currentY = atof(y3.c_str());
		setCMLType("CML2");
	} else if (y2 != "") {
		currentY = atof(y2.c_str());
		setCMLType("CML2");
	}
	if (z3 != "") {
		currentZ = atof(z3.c_str());
		setCMLType("CML2");
	}
// check other attributes
	for (int i = 0; i < atts.size(); i++) {
		if (atts[i].first == "elementType") {
	    } else if (atts[i].first == "id") {
	    } else if (atts[i].first == "formalCharge") {
	    } else if (atts[i].first == "x2") {
	    } else if (atts[i].first == "y2") {
	    } else if (atts[i].first == "x3") {
	    } else if (atts[i].first == "y3") {
	    } else if (atts[i].first == "z3") {
	    } else {
			cout << "IGNORED atom attribute: " << atts[i].first << endl;
		}
	}
}

bool processAtomArrayChild() {
	vector <string> strings;

	string builtin = getAttribute(currentAtts, "builtin");
	if (builtin == "") {
		cerr << "must have builtin attribute on: " << currentElem << endl;
	}
	pcdata += "\n";
	tokenize(strings, pcdata, " \n\t");
	if (natoms == 0) {
		natoms = strings.size();
		if (natoms == 0) {
			cerr << "no atoms in array: " << pcdata << endl;
		}
	}
	if (natoms != strings.size()) {
		cerr << "inconsistent atoms in arrays: " << pcdata << endl;
	}
	for (int i = 0; i < natoms; i++) {
		if (builtin == "elementType") {
			atomicNumVector.push_back(etab.GetAtomicNum((char*)strings[i].c_str()));
		} else if (builtin == "atomId") {
			idVector.push_back(strings[i]);
		} else if (builtin == "formalCharge") {
			formalChargeVector.push_back(atoi((char*)strings[i].c_str()));
		} else if (builtin == "x2") {
			x2Vector.push_back(atof((char*)strings[i].c_str()));
		} else if (builtin == "y2") {
			y2Vector.push_back(atof((char*)strings[i].c_str()));
		} else if (builtin == "x3") {
			x3Vector.push_back(atof((char*)strings[i].c_str()));
		} else if (builtin == "y3") {
			y3Vector.push_back(atof((char*)strings[i].c_str()));
		} else if (builtin == "z3") {
			z3Vector.push_back(atof((char*)strings[i].c_str()));
		}
	}
}

// adds builtin attributes for atom
bool processAtomBuiltin() {
	Vector v;


	string builtin = getAttribute(currentAtts, "builtin");
	if (builtin == "") {
		cerr << "No builtin attribute for <atom><" << currentElem << ">" << endl;
		return false;
	}
	setCMLType("CML1");
	processBuiltinPCDATA();
	if (false) {
	} else if (currentElem == "coordinate2") {
		vector <double> fv;
		processFloatTokens(fv, 2, pcdata);
// 3d takes precedence over 2D
	    if (builtin == "xy2") {
			if (cmlDimension != "3") {
				currentX = fv[0];
				currentY = fv[1];
			}
	    } else {
			cerr << "IGNORED coordinate2 builtin: " << builtin << endl;
			return false;
	    }
	} else if (currentElem == "coordinate3") {
		vector <double> fv;
		processFloatTokens(fv, 3, pcdata);
// 3d takes precedence over 2D
	    if (builtin == "xyz3") {
			currentX = fv[0];
			currentY = fv[1];
			currentZ = fv[2];
	    } else if (builtin == "xyzFract") {
			currentX = fv[0];
			currentY = fv[1];
			currentZ = fv[2];
			fractional = true;
	    } else {
			cerr << "IGNORED coordinate2 builtin: " << builtin << endl;
			return false;
	    }
	} else if (currentElem == "float") {
		double value = atof(pcdata.c_str());
// 3d takes precedence over 2D
	    if (builtin == "x2") {
			if (cmlDimension != "3") currentX = value;
	    } else if (builtin == "y2") {
			if (cmlDimension != "3") currentY = value;
	    } else if (builtin == "x3") {
			cmlDimension = "3";
			currentX = value;
	    } else if (builtin == "y3") {
			cmlDimension = "3";
			currentY = value;
	    } else if (builtin == "z3") {
			cmlDimension = "3";
			currentZ = value;
	    } else {
			cerr << "IGNORED float builtin: " << builtin << endl;
			return false;
	    }
	} else if (currentElem == "integer") {
		int ival = atoi(pcdata.c_str());
	    if (builtin == "formalCharge") {
			formalCharge = ival;
	    } else {
			cerr << "IGNORED integer builtin: " << builtin << endl;
			return false;
		}
	} else if (currentElem == "string") {
	    if (builtin == "elementType") {
  	        atomicNum = etab.GetAtomicNum((char*)pcdata.c_str());
	    } else if (builtin == "atomId") {
			atomId = pcdata;
	    } else {
			cerr << "IGNORED string builtin: " << builtin << endl;
			return false;
		}
	}
	return true;
}

bool endAtom() {
	OBAtom atom;
	pair<string, OBAtom*> at;

	atom.SetAtomicNum(atomicNum);
	atom.SetFormalCharge(formalCharge);
	atom.SetVector(currentX, currentY, currentZ);

    molPtr->AddAtom(atom);
    int nat = molPtr->NumAtoms();
    OBAtom* atPtr = molPtr->GetAtom(nat);	// counts from 1?
// store atoms
	at.first = atomId;
	at.second = atPtr;
	atomIdVector.push_back(at);
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
		ofs << "    <atom";
		writeAttribute(ofs, "id", id);
		if (outputCML2) {
			writeAttribute(ofs, "elementType", elementType);
			if (charge != 0) writeAttribute(ofs, "formalCharge", charge);
			if (molPtr->HasNonZeroCoords()) {
				if (strcmp(dimension, "2D")) {
					writeAttribute(ofs, "x2", x);
					writeAttribute(ofs, "y2", y);
				} else if (strcmp(dimension, "3D")) {
					writeAttribute(ofs, "x3", x);
					writeAttribute(ofs, "y3", y);
					writeAttribute(ofs, "z3", z);
				}
			}
			ofs << ">" << endl;
		} else {
			ofs << ">" << endl;
			writeBuiltin(ofs, "elementType", elementType);
			if (charge != 0) writeBuiltin(ofs, "formalCharge", charge);
			if (molPtr->HasNonZeroCoords()) {
				if (strcmp(dimension, "2D")) {
					writeBuiltin(ofs, "x2", x);
					writeBuiltin(ofs, "y2", y);
				} else if (strcmp(dimension, "3D")) {
					writeBuiltin(ofs, "x3", x);
					writeBuiltin(ofs, "y3", y);
					writeBuiltin(ofs, "z3", z);
				}
			}
		}
		ofs << "    </atom>" << endl;
	} else {
		appendToArray(idArray, id);
		appendToArray(elementArray, elementType);
		appendToArray(chargeArray, charge);
		if (molPtr->HasNonZeroCoords()) {
			if (strcmp(dimension, "2D")) {
				appendToArray(x2Array, x);
				appendToArray(y2Array, y);
			} else if (strcmp(dimension, "3D")) {
				appendToArray(x3Array, x);
				appendToArray(y3Array, y);
				appendToArray(z3Array, z);
			}
		}
	}
}

void processStringTokens(vector <string> &v, int n, string att) {
	if (att == "") return;
	vector <string> sv;
	att += " ";
	tokenize(sv, att, " \n");
	if (sv.size() != n) {
		cerr << "inconsistent array attribute sizes: " << sv.size() << "/" << n << endl;
		return;
	}
	for (int i = 0; i < n; i++) v[i] = sv[i];
}

void processIntTokens(vector <int> &v, int n, string att) {
	if (att == "") return;
	vector <string> sv;
	att += " ";
	tokenize(sv, att, " \n");
	if (sv.size() != n) {
		cerr << "inconsistent array attribute sizes: " << sv.size() << "/" << n << endl;
		return;
	}
	for (int i = 0; i < n; i++) v[i] = atoi((char*)sv[i].c_str());
}

void processFloatTokens(vector <double> &v, int n, string att) {
	if (att == "") return;
	vector <string> sv;
	att += " ";
	tokenize(sv, att, " \n");
	if (sv.size() != n) {
		cerr << "inconsistent array attribute sizes: " << sv.size() << "/" << n << endl;
		return;
	}
	for (int i = 0; i < n; i++) v[i] = atof((char*)sv[i].c_str());
}

bool startAtomArray(vector <pair<string,string> > &atts) {
	vector <string> sv;
	string atomID = getAttribute(atts, "atomID");
	if (atomID == "") return false;
	setCMLType("CML2");
	atomId += " ";
	tokenize(sv, atomID, " \n");
	int natoms = sv.size();
	processStringTokens(idVector, natoms, atomID);
	processStringTokens(elementTypeVector, natoms, getAttribute(atts, "elementType"));
	processIntTokens(formalChargeVector, natoms, getAttribute(atts, "formalCharge"));
	processIntTokens(hydrogenCountVector, natoms, getAttribute(atts, "hydrogenCount"));
	processFloatTokens(x2Vector, natoms, getAttribute(atts, "x2"));
	processFloatTokens(y2Vector, natoms, getAttribute(atts, "y2"));
	processFloatTokens(x3Vector, natoms, getAttribute(atts, "x3"));
	processFloatTokens(y3Vector, natoms, getAttribute(atts, "y3"));
	processFloatTokens(z3Vector, natoms, getAttribute(atts, "z3"));
}

bool endAtomArray() {
	pair<string, OBAtom*> at;
	if (inputCML2) {
		for (int i = 0; i < natoms; i++) {
			OBAtom atom;
			atom.SetAtomicNum(atomicNumVector[i]);
			if (elementTypeVector.size() > 0) {
				atom.SetAtomicNum(etab.GetAtomicNum((char*)elementTypeVector[i].c_str()));
			}
			if (formalChargeVector.size() > 0) atom.SetFormalCharge(formalChargeVector[i]);
			Vector v;
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
	for (atom = molPtr->BeginAtom(i);atom;atom = molPtr->NextAtom(i)) {
		WriteAtom(ofs, atom, ++count);
	}
	if (outputArray) {
		ofs << "<stringArray builtin=\"atomId\">" << idArray << "</stringArray>" << endl;
		ofs << "<stringArray builtin=\"elementType\">" << elementArray << "</stringArray>" << endl;
		ofs << "<integerArray builtin=\"formalCharge\">" << chargeArray << "</integerArray>" << endl;
		if (molPtr->HasNonZeroCoords()) {
			if (strcmp(dimension, "2D")) {
				ofs << "<floatArray builtin=\"x2\">" << x2Array << "</floatArray>" << endl;
				ofs << "<floatArray builtin=\"y2\">" << y2Array << "</floatArray>" << endl;
			} else if (strcmp(dimension, "3D")) {
				ofs << "<floatArray builtin=\"x3\">" << x3Array << "</floatArray>" << endl;
				ofs << "<floatArray builtin=\"y3\">" << y3Array << "</floatArray>" << endl;
				ofs << "<floatArray builtin=\"z3\">" << z3Array << "</floatArray>" << endl;
			}
		}
	}
}

bool startAtomParity(vector <pair<string,string> > &atts) {
	atomRefs4 = getAttribute(atts, "atomRefs4");
}

bool endAtomParity(vector <pair<string,string> > &atts) {
	pair <vector<OBAtom*>, double> ap;
	vector <OBAtom*> atomRef;
	getAtomRefs(4, atomRef, atomRefs4);
	if (atomRef.size() != 4) {
		cerr << "atomRefs4 must referemce 4 atoms" << endl;
		return false;
	}
	for (int i = 0; i < 4; i++) ap.first.push_back(atomRef[i]);
	setCMLType("CML2");
	ap.second = atof((char*)pcdata.c_str());
	atomParityVector.push_back(ap);
}

bool WriteAtomParity(ostream &ofs) {
	cout << "WriteAtomParity NYI" << endl;
}

// --------------------<bond><bondArray>----------------
bool startBond(vector <pair<string,string> > &atts) {
	vector <string> badAtts = getUnknownAttributes(BOND_ATTRIBUTE_VECTOR, atts);
	if (badAtts.size() > 0) {
		cerr << "unknown attributes on <bond>: ";
		printVector(badAtts, cerr);
		cerr << endl;
	}

	vector <string> atomRefs;

	bondBeginAtom = "";
	bondEndAtom = "";
    orderString = getAttribute(currentAtts, "order");
    stereoString = getAttribute(currentAtts, "stereo");
	tokenize(atomRefs, (char*)getAttribute(currentAtts, "atomRefs2").c_str(), " \n\t,");
	if (atomRefs.size() == 0) {
		return false;
	} else if (atomRefs.size() != 2) {
		cerr << "must have 2 atom Refs per bond" << endl;
		return false;
	} else {
		setCMLType("CML2");
	}
	bondBeginAtom = atomRefs[0];
	bondEndAtom = atomRefs[1];
}

// adds builtin attributes for bond
bool processBondBuiltin() {

	string builtin = getAttribute(currentAtts, "builtin");
	if (builtin == "") {
		cerr << "No builtin attribute for <bond><" << currentElem << ">" << endl;
		return false;
	}
	setCMLType("CML1");
	if (currentElem == "float") {
		double value = atof(pcdata.c_str());
	    if (false) {
		} else if (builtin == "length") {
			length = value;
	    } else {
			cerr << "IGNORED float builtin for bond: " << builtin << endl;
			return false;
		}
	} else if (currentElem == "integer") {
		int ival = atoi(pcdata.c_str());
	    if (false) {
	    } else {
			cerr << "IGNORED integer builtin: " << builtin << endl;
			return false;
		}
	} else if (currentElem == "string") {
	    if (builtin == "atomRef") {
			if (bondBeginAtom == "") {
				bondBeginAtom = pcdata;
			} else {
				if (bondEndAtom != "") {
					cerr << "too many atomRef builtins" << endl;
					return false;
				}
				bondEndAtom = pcdata;
			}
		} else if (builtin == "order") {
			orderString = pcdata;
	    } else if (builtin == "stereo") {
			stereoString = pcdata;
	    } else {
			cerr << "IGNORED integer builtin: " << builtin << endl;
			return false;
		}
	}
	return true;
}


bool processBondArrayChild() {
	vector <string> strings;

	string builtin = getAttribute(currentAtts, "builtin");
	if (builtin == "") {
		cerr << "must have builtin attribute on: " << currentElem << endl;
	}
	pcdata += "\n";
	tokenize(strings, pcdata, " \n\t");
	if (nbonds == 0) {
		nbonds = strings.size();
		if (nbonds == 0) {
			cerr << "no bonds in array: " << pcdata << endl;
		}
	}
	if (nbonds != strings.size()) {
		cerr << "inconsistent bonds in arrays: " << pcdata << endl;
	}
	bool atomRef1 = (atomRef1Vector.size() == 0);
	for (int i = 0; i < nbonds; i++) {
		if (builtin == "atomRef") {
			if (atomRef1) {
				atomRef1Vector.push_back(strings[i]);
			} else {
				atomRef2Vector.push_back(strings[i]);
			}
		} else if (builtin == "order") {
			orderVector.push_back(strings[i]);
		} else if (builtin == "stereo") {
			stereoVector.push_back(strings[i]);
		}
	}
}

bool endBond() {
	pair <vector<OBAtom*>, double> len;
	OBBond bond;

    bondPtr = &bond;
	OBAtom* beginAtomPtr = getAtomPtr(bondBeginAtom);
	OBAtom* endAtomPtr = getAtomPtr(bondEndAtom);
	if (beginAtomPtr == 0 || endAtomPtr == 0) {
		cerr << "could not find atom refs in bond" << endl;
		return false;
	}
    bondPtr->SetBegin(beginAtomPtr);
    bondPtr->SetEnd(endAtomPtr);
    if (orderString != "") bondPtr->SetBO(getBabelBondOrder(orderString));
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
		ofs << "    <bond";
//		writeAttribute(ofs, "id", id);
		if (outputCML2) {
			string aa(atomRef1+" "+atomRef2);
			writeAttribute(ofs, "atomRefs2", aa);
			if (strcmp(boChar, "") != 0) writeAttribute(ofs, "order", boChar);
			ofs << "/>" << endl;
		} else {
			ofs << ">" << endl;
			writeBuiltin(ofs, "atomRef", atomRef1);
			writeBuiltin(ofs, "atomRef", atomRef2);
			if (strcmp(boChar, "") != 0) writeBuiltin(ofs, "order", boChar);
			ofs << "    </bond>" << endl;
		}
	} else {
		appendToArray(atomRef1Array, atomRef1);
		appendToArray(atomRef2Array, atomRef2);
		appendToArray(orderArray, boChar);
	}
}

bool startBondArray(vector <pair<string,string> > &atts) {
	vector <string> sv;
	string atomRef1 = getAttribute(atts, "atomRef1");
	if (atomRef1 == "") return false;
	setCMLType("CML2");
	atomRef1 += " ";
	tokenize(sv, atomRef1, " \n");
	int nbonds = sv.size();
	processStringTokens(atomRef1Vector, nbonds, atomRef1);
	processStringTokens(atomRef2Vector, nbonds, getAttribute(atts, "atomRef2"));
	processStringTokens(orderVector, nbonds, getAttribute(atts, "order"));
	processStringTokens(stereoVector, nbonds, getAttribute(atts, "stereo"));
}

bool endBondArray() {
	if (inputCML2) {
		if (atomRef1Vector.size() == 0 ||
			atomRef2Vector.size() == 0) {
			cerr << "atomRef arrays must be given for bonds" << endl;
		}
		for (int i = 0; i < nbonds; i++) {
			OBBond bond;
			bondPtr = &bond;
			OBAtom* beginAtomPtr = getAtomPtr(atomRef1Vector[i]);
			OBAtom* endAtomPtr = getAtomPtr(atomRef2Vector[i]);
			if (beginAtomPtr == 0 || endAtomPtr == 0) {
				cerr << "could not find atom refs in bond" << endl;
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

	ofs << "  <bondArray";
	ofs << ">" << endl;
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
		ofs << "<stringArray builtin=\"atomRef\">" << atomRef1Array << "</stringArray>" << endl;
		ofs << "<stringArray builtin=\"atomRef\">" << atomRef2Array << "</stringArray>" << endl;
		ofs << "<stringArray builtin=\"order\">" << orderArray << "</stringArray>" << endl;
	}
	ofs << "  </bondArray>" << endl;
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
		cerr << "unknown attributes on <crystal>: ";
		printVector(badAtts, cerr);
		cerr << endl;
	}
// check other attributes
	for (int i = 0; i < atts.size(); i++) {
		if (false) {
	    } else if (atts[i].first == "id") {
	    } else if (atts[i].first == "title") {
	    } else if (atts[i].first == "convention") {
	    } else if (atts[i].first == "spaceGroup") {
	    } else if (atts[i].first == "pointGroup") {
	    } else {
			cout << "IGNORED crystal attribute: " << atts[i].first << endl;
		}
	}
}

// adds builtin attributes for cryst
bool processCrystalBuiltin() {
	Vector v;

	string builtin = getAttribute(currentAtts, "builtin");
	if (builtin == "") {
		cerr << "No builtin attribute for <cryst><" << currentElem << ">" << endl;
		return false;
	}
	setCMLType("CML1");
	processBuiltinPCDATA();
	if (currentElem == "float") {
		double f = atof((char*)pcdata.c_str());
		if (false) {
		} else if (currentElem == "acell") {
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
			cerr << "IGNORED float builtin: " << builtin << endl;
			return false;
		}
	} else {
		cerr << "IGNORED builtin for " << currentElem << " in crystal; " << builtin << endl;
	}
	return true;
}

bool endCrystal() {
}

// not yet finished
bool WriteCrystal(ostream &ofs) {
	ofs << "<crystal";
	ofs << ">";
	ofs << "</crystal>" << endl;
}

// --------------------<electron>----------------

bool startElectron(vector <pair<string,string> > &atts) {
	vector <string> badAtts = getUnknownAttributes(ELECTRON_ATTRIBUTE_VECTOR, atts);
	if (badAtts.size() > 0) {
		cerr << "unknown attributes on <electron>: ";
		printVector(badAtts, cerr);
		cerr << endl;
	}
// check other attributes
	for (int i = 0; i < atts.size(); i++) {
		if (false) {
	    } else if (atts[i].first == "id") {
	    } else if (atts[i].first == "title") {
	    } else if (atts[i].first == "convention") {
	    } else {
			cout << "IGNORED electron attribute: " << atts[i].first << endl;
		}
	}
}

bool endElectron() {
//	pair <vector<OBAtom*>, double> electron;
//	electronVector.push_back(electron);
}

// not yet finished
bool WriteElectron(ostream &ofs) {
	ofs << "<electron";
	ofs << ">";
	ofs << "</electron>" << endl;
}

// --------------------<feature>----------------

bool startFeature(vector <pair<string,string> > &atts) {
	vector <string> badAtts = getUnknownAttributes(FEATURE_ATTRIBUTE_VECTOR, atts);
	if (badAtts.size() > 0) {
		cerr << "unknown attributes on <feature>: ";
		printVector(badAtts, cerr);
		cerr << endl;
	}
// check other attributes
	for (int i = 0; i < atts.size(); i++) {
		if (false) {
	    } else if (atts[i].first == "id") {
	    } else if (atts[i].first == "title") {
	    } else if (atts[i].first == "convention") {
	    } else {
			cout << "IGNORED feature attribute: " << atts[i].first << endl;
		}
	}
}

bool endFeature() {
}

bool WriteFeature(ostream &ofs) {
	ofs << "<feature";
	ofs << ">";
	ofs << "</feature>" << endl;
}

// --------------------<formula>----------------

bool startFormula(vector <pair<string,string> > &atts) {
	vector <string> badAtts = getUnknownAttributes(FORMULA_ATTRIBUTE_VECTOR, atts);
	if (badAtts.size() > 0) {
		cerr << "unknown attributes on <formula>: ";
		printVector(badAtts, cerr);
		cerr << endl;
	}
// check other attributes
	for (int i = 0; i < atts.size(); i++) {
		if (false) {
	    } else if (atts[i].first == "id") {
	    } else if (atts[i].first == "title") {
	    } else if (atts[i].first == "convention") {
	    } else {
			cout << "IGNORED formula attribute: " << atts[i].first << endl;
		}
	}
}

bool endFormula() {
}

bool WriteFormula(ostream &ofs) {
	ofs << "<formula";
	ofs << ">";
	ofs << "</formula>" << endl;
}

// --------------------<length>----------------

bool startLength(vector <pair<string,string> > &atts) {
	vector <string> badAtts = getUnknownAttributes(LENGTH_ATTRIBUTE_VECTOR, atts);
	if (badAtts.size() > 0) {
		cerr << "unknown attributes on <length>: ";
		printVector(badAtts, cerr);
		cerr << endl;
	}
	lengthUnits = "angstrom";
	atomRefs2Vector.clear();
// check other attributes
	for (int i = 0; i < atts.size(); i++) {
		if (false) {
	    } else if (atts[i].first == "id") {
	    } else if (atts[i].first == "title") {
	    } else if (atts[i].first == "convention") {
	    } else if (atts[i].first == "atomRefs2") {
			getAtomRefs(2, atomRefs2Vector, atts[i].second);
	    } else if (atts[i].first == "units") {
			lengthUnits = atts[i].second;
	    } else {
			cout << "IGNORED length attribute: " << atts[i].first << endl;
		}
	}
}

bool endLength() {
	pair <vector<OBAtom*>, double> length;
	if (atomRefs2Vector.size() != 2) {
		cerr << "must have defined 2 atoms for length" << endl;
	}
	for (int i = 0; i < 2; i++) {
		length.first.push_back(atomRefs2Vector[i]);
	}
	length.second = atof((char*)pcdata.c_str());
	lengthVector.push_back(length);
}

// not yet finished
bool WriteLength(ostream &ofs, pair <vector<OBAtom*>, double> length) {
	ofs << "<length";
	ofs << " atomRefs2=\"a" << length.first[0]->GetIdx() << " a" << length.first[1]->GetIdx() << "\">";
	ofs << length.second;
	ofs << "</length>" << endl;
}

// ------------------ <molecule> ----------------

bool startMolecule(vector <pair<string,string> > &atts) {
	vector <string> badAtts = getUnknownAttributes(MOLECULE_ATTRIBUTE_VECTOR, atts);
	if (badAtts.size() > 0) {
		cerr << "unknown attributes on <molecule>: ";
		printVector(badAtts, cerr);
		cerr << endl;
	}

	molPtr->BeginModify();
	molPtr->ReserveAtoms(ATOM_SIZE);
	molPtr->SetTitle((char*)getAttribute(atts, "title").c_str());
}

bool debugMolecule(ostream &ofs) {
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
		for (int i = 0; i < lengthVector.size(); i++) {
			pair <vector<OBAtom*>, double> length = lengthVector[i];
			WriteLength(ofs, length);
		}
	}
	if (angleVector.size() > 0) {
		ofs << "Angles: " << endl;
		for (int i = 0; i < angleVector.size(); i++) {
			pair <vector<OBAtom*>, double> angle = angleVector[i];
			WriteAngle(ofs, angle);
		}
	}
	if (torsionVector.size() > 0) {
		ofs << "Torsions: " << endl;
		for (int i = 0; i < torsionVector.size(); i++) {
			pair <vector<OBAtom*>, double> torsion = torsionVector[i];
			WriteTorsion(ofs, torsion);
		}
	}
}

// returns index of length between atoms ( starts at 0; -1 = not found)
int getLengthIndex(OBAtom* a0, OBAtom* a1) {
	for (int i = 0; i < lengthVector.size(); i++) {
		if (a0 == lengthVector[i].first[0] &&
			a1 == lengthVector[i].first[1]) return i;
		if (a0 == lengthVector[i].first[1] &&
			a1 == lengthVector[i].first[0]) return i;
	}
	return -1;
}

// returns index of angle between atoms ( starts at 0; -1 = not found)
int getAngleIndex(OBAtom* a0, OBAtom* a1, OBAtom* a2) {
	for (int i = 0; i < angleVector.size(); i++) {
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
	for (int i = 0; i < torsionVector.size(); i++) {
		if (a0 == torsionVector[i].first[0] &&
			a1 == torsionVector[i].first[1] &&
			a2 == torsionVector[i].first[2] &&
			a3 == torsionVector[i].first[3]) return i+1;
		if (a0 == torsionVector[i].first[3] &&
			a1 == torsionVector[i].first[2] &&
			a2 == torsionVector[i].first[1] &&
			a3 == torsionVector[i].first[0]) return -(i+1);
	}
	return 0;
}

//attempts to get torsions in right order for internals
// returns serial of torsion starting at 1. If not found
// returns zero. If torsion is "wrong way round", returns
// negative serial
int getFirstTorsionIndexForAtom(OBAtom* a0) {
	int k = a0->GetIdx();
	for (int i = 0; i < torsionVector.size(); i++) {
		if (a0 == torsionVector[i].first[0]) {
			bool ok = true;
			for (int j = 1; j <= 3; j++) {
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
			for (int j = 0; j < 3; j++) {
				OBAtom* atPtr = torsionVector[i].first[j];
				if (atPtr->GetIdx() > k) {
					ok = true;
					break;
				}
			}
			if (ok) return -(i+1);
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
	int nTors = torsionVector.size();
	if (nTors == 0) return;
	int nAng = angleVector.size();
	if (nAng == 0) return;
	int nLen = lengthVector.size();
	if (nLen == 0) return;
	if (nLen + 1 < molPtr->NumAtoms()) {
		if (nLen > 0) cout << "Not enough lengths to generate all internals" << endl;
	}
	if (nAng + 2 < molPtr->NumAtoms()) {
		if (nAng > 0) cout << "Not enough angles to generate all internals" << endl;
	}
	if (nTors + 3 < molPtr->NumAtoms()) {
		if (nTors > 0) cout << "Not enough torsions to generate all internals" << endl;
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
	int idx = getLengthIndex(at0, at1);
	if (idx == -1) {
		cerr << "cannot find length: " << at0->GetIdx() << "/" << at1->GetIdx() << endl;
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
		cerr << "cannot find length: " << at1 << "/" << at2 << endl;
		return;
	}
	coord->_dst = lengthVector[idx].second;
	idx = getAngleIndex(at0, at1, at2);
	if (idx == -1) {
		cerr << "cannot find angle: " << at0 << "/" << at1 << "/" << at2 << endl;
		return;
	}
	coord->_ang = angleVector[idx].second;
	internalVector.push_back(coord);

	for (int i = 3; i < molPtr->NumAtoms(); i++) {
		OBAtom* at0 = molPtr->GetAtom(i+1);
		idx = getFirstTorsionIndexForAtom(at0);
		if (idx == 0) {
			cerr << "cannot find torsion... " << endl;
			return;
		}
		int iTor = (idx > 0) ? idx-1 : -idx - 1;
		at0 = torsionVector[iTor].first[0];
		OBAtom* at1 = torsionVector[iTor].first[1];
		OBAtom* at2 = torsionVector[iTor].first[2];
		OBAtom* at3 = torsionVector[iTor].first[3];
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
			cerr << "cannot find length: " << at2 << "/" << at3 << endl;
			return;
		}
		coord->_dst = lengthVector[idx].second;
		int idx = getAngleIndex(at1, at2, at3);
		if (idx == -1) {
			cerr << "cannot find angle: " << at1 << "/" << at2 << "/" << at3 << endl;
			return;
		}
		coord->_ang = angleVector[idx].second;
		coord->_tor = torsionVector[iTor].second;
		internalVector.push_back(coord);
	}
	for (int i = 0; i < internalVector.size(); i++) {
		OBInternalCoord* coord = internalVector[i];
		int aa = (coord->_a != 0) ? coord->_a->GetIdx() : 0;
		int bb = (coord->_b != 0) ? coord->_b->GetIdx() : 0;
		int cc = (coord->_c != 0) ? coord->_c->GetIdx() : 0;
		cout << "a" << cc << " ";
		cout << "a" << bb << ":";
		cout << "a" << aa << ":";
		cout << "a" << (i+1) << ":";
		cout << coord->_dst << " ";
		cout << coord->_ang << " ";
		cout << coord->_tor << endl;
	}
}

bool endMolecule() {
	debugMolecule(cout);
	generateInternals();
    InternalToCartesian(internalVector, *molPtr);

	molPtr->EndModify();

	molPtr->ConnectTheDots();

	if (outputDebug) {
		debug(cout);
	}
}

bool WriteMolecule(ostream &ofs) {

	ofs << "<molecule";
	writeAttribute(ofs, "title", molPtr->GetTitle());
// a mechanism is needed for IDs on elements
	writeAttribute(ofs, "id", "m1");
	ofs << ">" << endl;

	if (molPtr->HasData(obCommentData)) {
		OBCommentData *cd = (OBCommentData*)molPtr->GetData(obCommentData);
		ofs << "<string title=\"comment\">" << cd->GetData() << "</comment>" << endl;
	}

	if (outputDebug) debug(ofs);

	WriteAtomArray(ofs);
	WriteBondArray(ofs);

	vector<OBGenericData*>::iterator k;
	vector<OBGenericData*> vdata = molPtr->GetData();
	for (k = vdata.begin();k != vdata.end();k++) {
		if ((*k)->GetDataType() == obPairData) {
			ofs << "<string title=\"" << (*k)->GetAttribute() << "\">"
				<< ((OBPairData*)(*k))->GetValue() << "</string>" << endl;
		}
	}

	ofs << "</molecule>" << endl;

	return(true);
}

// --------------------<reaction>----------------

bool startReaction(vector <pair<string,string> > &atts) {
	vector <string> badAtts = getUnknownAttributes(REACTION_ATTRIBUTE_VECTOR, atts);
	if (badAtts.size() > 0) {
		cerr << "unknown attributes on <reaction>: ";
		printVector(badAtts, cerr);
		cerr << endl;
	}
// check other attributes
	for (int i = 0; i < atts.size(); i++) {
		if (false) {
	    } else if (atts[i].first == "id") {
	    } else if (atts[i].first == "title") {
	    } else if (atts[i].first == "convention") {
	    } else {
			cout << "IGNORED reaction attribute: " << atts[i].first << endl;
		}
	}
}

bool endReaction() {
}

bool WriteReaction(ostream &ofs) {
	ofs << "<reaction";
	ofs << ">";
	ofs << "</reaction>" << endl;
}

// --------------------<scalar>----------------

bool startScalar(vector <pair<string,string> > &atts) {
	vector <string> badAtts = getUnknownAttributes(SCALAR_ATTRIBUTE_VECTOR, atts);
	if (badAtts.size() > 0) {
		cerr << "unknown attributes on <scalar>: ";
		printVector(badAtts, cerr);
		cerr << endl;
	}
// check other attributes
	for (int i = 0; i < atts.size(); i++) {
		if (false) {
	    } else if (atts[i].first == "id") {
	    } else if (atts[i].first == "title") {
	    } else if (atts[i].first == "convention") {
	    } else if (atts[i].first == "dataType") {
			scalarDataType = atts[i].second;
	    } else if (atts[i].first == "units") {
			scalarUnits = atts[i].second;
	    } else {
			cout << "IGNORED scalar attribute: " << atts[i].first << endl;
		}
	}
}

bool endScalar() {
	string title = getAttribute(currentAtts, "title");
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
}

// not yet finished
bool WriteScalar(ostream &ofs) {
	ofs << "<scalar/>" << endl;
}

// --------------------<sequence>----------------

bool startSequence(vector <pair<string,string> > &atts) {
	vector <string> badAtts = getUnknownAttributes(SEQUENCE_ATTRIBUTE_VECTOR, atts);
	if (badAtts.size() > 0) {
		cerr << "unknown attributes on <sequence>: ";
		printVector(badAtts, cerr);
		cerr << endl;
	}
// check other attributes
	for (int i = 0; i < atts.size(); i++) {
		if (false) {
	    } else if (atts[i].first == "id") {
	    } else if (atts[i].first == "title") {
	    } else if (atts[i].first == "convention") {
	    } else {
			cout << "IGNORED sequence attribute: " << atts[i].first << endl;
		}
	}
}

bool endSequence() {
}

// not yet finished
bool WriteSequence(ostream &ofs) {
	ofs << "<sequence";
	ofs << ">";
	ofs << "</sequence>" << endl;
}

// -------------------------<stereo>-------------------

bool startStereo(vector <pair<string,string> > &atts) {
	atomRefs4 = getAttribute(atts, "atomRefs4");
}

bool endStereo(vector <pair<string,string> > &atts) {
	pair <vector<OBAtom*>, string> st;
	vector <OBAtom*> atomRef;
	getAtomRefs(4, atomRef, atomRefs4);
	if (atomRef.size() != 4) {
		cerr << "atomRefs4 must referemce 4 atoms" << endl;
		return false;
	}
	for (int i = 0; i < 4; i++) st.first.push_back(atomRef[i]);
	setCMLType("CML2");
	st.second = pcdata;
	stereoSVector.push_back(st);
}

bool WriteStereo(ostream &ofs) {
	cout << "WriteStereo NYI" << endl;
}


// -------------------------<symmetry>-------------------

bool startSymmetry(vector <pair<string,string> > &atts) {
	vector <string> badAtts = getUnknownAttributes(SYMMETRY_ATTRIBUTE_VECTOR, atts);
	if (badAtts.size() > 0) {
		cerr << "unknown attributes on <symmetry>: ";
		printVector(badAtts, cerr);
		cerr << endl;
	}
	spacegroup = getAttribute(atts, "spacegroup");
	pointgroup = getAttribute(atts, "pointgroup");
// check other attributes
	for (int i = 0; i < atts.size(); i++) {
		if (false) {
	    } else if (atts[i].first == "id") {
	    } else if (atts[i].first == "title") {
	    } else if (atts[i].first == "convention") {
	    } else {
			cout << "IGNORED symmetry attribute: " << atts[i].first << endl;
		}
	}
}

bool endSymmetry(vector <pair<string,string> > &atts) {
}

bool WriteSymmetry(ostream &ofs) {
	ofs << "<symmetry";
	ofs << " spacegroup=\"" << spacegroup << "\"";
	ofs << " pointgroup=\"" << pointgroup << "\"";
	ofs << ">";
	ofs << "</symmetry>" << endl;
}


// -------------------------<string>-------------------
bool addString() {
	string title = getAttribute(currentAtts, "title");
	if (title != "") {
		OBPairData *dp = new OBPairData;
		dp->SetAttribute(title);
		dp->SetValue(pcdata);
		molPtr->SetData(dp);
	}
}

// --------------------<torsion>/<torsionArray>----------------

bool startTorsion(vector <pair<string,string> > &atts) {
	vector <string> badAtts = getUnknownAttributes(TORSION_ATTRIBUTE_VECTOR, atts);
	if (badAtts.size() > 0) {
		cerr << "unknown attributes on <torsion>: ";
		printVector(badAtts, cerr);
		cerr << endl;
	}
	torsionUnits = "degrees";
	atomRefs4Vector.clear();
// check other attributes
	for (int i = 0; i < atts.size(); i++) {
		if (false) {
	    } else if (atts[i].first == "id") {
	    } else if (atts[i].first == "title") {
	    } else if (atts[i].first == "convention") {
	    } else if (atts[i].first == "atomRefs") {
			setCMLType("CML1");
			getAtomRefs(4, atomRefs4Vector, atts[i].second);
	    } else if (atts[i].first == "atomRefs4") {
			setCMLType("CML2");
			getAtomRefs(4, atomRefs4Vector, atts[i].second);
	    } else if (atts[i].first == "units") {
			torsionUnits = atts[i].second;
	    } else {
			cout << "IGNORED torsion attribute: " << atts[i].first << endl;
		}
	}
}

bool endTorsion() {
	pair <vector<OBAtom*>, double> torsion;
	if (atomRefs4Vector.size() != 4) {
		cerr << "must have defined 4 atoms for torsion" << endl;
	}
	for (int i = 0; i < 4; i++) {
		torsion.first.push_back(atomRefs4Vector[i]);
	}
	torsion.second = atof((char*)pcdata.c_str());
	torsionVector.push_back(torsion);
}

// not yet finished
bool WriteTorsion(ostream &ofs, pair <vector<OBAtom*>, double> torsion) {
	ofs << "<torsion";
	ofs << " atomRefs4=\"a" << torsion.first[0]->GetIdx() << " a" << torsion.first[1]->GetIdx() << " a" << torsion.first[2]->GetIdx() << " a" << torsion.first[3]->GetIdx() << "\">";
	ofs << torsion.second;
	ofs << "</torsion>" << endl;
}


/**------------------Babel INPUT OUTPUT --------------------------*/

bool ReadCML(istream &ifs,OBMol &mol, char *title) {
	molPtr = &mol;
	ReadXML(ifs);
	return true;
}

// output routines

bool WriteCML(ostream &ofs,OBMol &mol,char *dim, char* xmlOptions) {

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
		cout << "output options: ";
		if (outputCML2) cout << "CML2 ";
		if (outputDoctype) cout << "doctype ";
		if (outputPretty) cout << "pretty ";
		if (outputNamespace) cout <<  "namespace ";
		cout << outputPrefix << " ";
		if (outputArray) cout << "arrays ";
		if (outputDebug) cout << "debug ";
		cout << endl;
	}

	molPtr = &mol;
	dimension = dim;
	WriteMolecule(ofs);
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

	for (int i = 0; i < molPtr->NumAtoms(); i++) {
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

	for (int i = 0; i < molPtr->NumBonds(); i++) {
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
