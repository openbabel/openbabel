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

// ----------------CML events----------------
// <atom>
bool startAtom(vector <pair<string,string> > &atts);
// process builtin children of <atom>
bool processAtomBuiltin();
// </atom>
bool endAtom();
bool WriteAtom(ostream &ofs, OBAtom* atom, int count);
vector <string> ATOM_ATTRIBUTE_VECTOR;
string ATOM_ATTRIBUTES = "id title convention dictRef";
string ATOM_BUILTINS =
	"x2 y2 x3 y3 z3 xy2 xyz3 xFract yFract zFract elementType formalCharge hydrogenCount";

// <atomArray>
bool startAtomArray(vector <pair<string,string> > &atts);
// </atomArray>
bool endAtomArray();
bool WriteAtomArray(ostream &ofs);
vector <string> ATOMARRAY_ATTRIBUTE_VECTOR;
string ATOMARRAY_ATTRIBUTES = "id title convention dictRef";

// <bond>
bool startBond(vector <pair<string,string> > &atts);
// process builtin children of <bond>
bool processBondBuiltin();
// </bond>
bool endBond();
bool WriteBond(ostream &ofs, OBBond* bond);
vector <string> BOND_ATTRIBUTE_VECTOR;
string BOND_ATTRIBUTES = "id title convention dictRef";
string BOND_BUILTINS =  "atomRef atomRefs2 order stereo";

// <bondArray>
bool startBondArray(vector <pair<string,string> > &atts);
// </bondArray>
bool endBondArray();
bool WriteBondArray(ostream &ofs);
vector <string> BONDARRAY_ATTRIBUTE_VECTOR;
string BONDARRAY_ATTRIBUTES = "id title convention dictRef";

// <cml>
bool startCML(vector <pair<string,string> > &atts);
// </cml>
bool endCML();
string CML_ATTRIBUTES =  "id title convention";
vector <string> CML_ATTRIBUTE_VECTOR;

// <molecule>
bool startMolecule(vector <pair<string,string> > &atts);
// </molecule>
bool endMolecule();
bool WriteMolecule(ostream &ofs);
vector <string> MOLECULE_ATTRIBUTE_VECTOR;
string MOLECULE_ATTRIBUTES =  "id title count convention";

// <string>
bool startString(vector <pair<string,string> > &atts);
// </string>
bool endString();
vector <string> STRING_ATTRIBUTE_VECTOR;
string STRING_ATTRIBUTES =  "id title count convention builtin";

// adds  <string>
bool addString();

// normalize PCDATA for builtins
void processBuiltinPCDATA();
// <fooArray> children of atomArray
bool processAtomArrayChild();
// <fooArray> children of bondArray
bool processBondArrayChild();

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

vector <string> idVector;
vector <int> atomicNumVector;
vector <int> formalChargeVector;
vector <float> x2Vector;
vector <float> y2Vector;
vector <float> x3Vector;
vector <float> y3Vector;
vector <float> z3Vector;

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

bool ReadXML(istream &ifs) {
	char buffer[BUFF_SIZE];
	size_t lt;
	size_t rt;

	currentElem = "";
	string token = "";
	bool lookForOpenTag = true;

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

void makeAllowedAttributeLists() {
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
	tokenize(MOLECULE_ATTRIBUTE_VECTOR, MOLECULE_ATTRIBUTES, " \n");
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

// ------------------------ SAX events -------------------
// SAX-like call back
void startDocument() {
  cout << "starting CML document; crude XML parser. Assumes well-formed; ignores DTDs and entities" << endl;
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
	startElement("", name, "", atts);
}

void startElement(string namespaceURI, string localName, string qName, vector<pair<string,string> > &atts) {
	if (currentElem != "") elementStack.push_back(currentElem);
	parent = currentElem;
	currentElem = trim(localName);
	currentAtts = atts;
	pcdata = "";
	useBuiltin = false;
	processAttributes(atts);
    if (localName == "molecule") {
		startMolecule(atts);
	} else if (localName == "atom") {
        startAtom(atts);
	} else if (localName == "atomArray") {
	} else if (localName == "bond") {
        startBond(atts);
	} else if (localName == "bondArray") {
	} else if (
		localName == "float" ||
		localName == "integer" ||
	    localName == "string") {
		if (parent == "atom" || parent == "bond") {
			if (inputCML2) {
				cerr << "conflict between CML1 and CML2" << endl;
			} else {
				inputCML1 = true;
			}
		}
	// delay processing till endElement as we need the pcdata
	} else if (
		localName == "floatArray" ||
		localName == "integerArray" ||
	    localName == "stringArray") {
		if (parent == "atomArray" || parent == "bondArray") {
			if (inputCML1) {
				cerr << "conflict between CML1 and CML2" << endl;
			} else {
				if (!inputCML2) cout << "reading CML as arrays" << endl;
				inputCML2 = true;
			}
		}
	// delay processing till endElement as we need the pcdata
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
			break;
		}
	}
	if (!nsExists) {
		namespaceVector.push_back(ns);
	}
}

void endElement(string s) {
	endElement("", s, "");
}

void endElement(string namespaceURI, string localName, string qName) {
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
		// test
	} else if (name == "bond") {
		endBond();
	} else if (name == "bondArray") {
		endBondArray();
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
	} else if (
		localName == "floatArray" ||
		localName == "integerArray" ||
	    localName == "stringArray") {
		if (parent == "atomArray") {
			processAtomArrayChild();
		} else if (parent == "bondArray") {
			processBondArrayChild();
		}
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

// ---------------------- CML events -----------------
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

bool endMolecule() {
	molPtr->EndModify();
	if (outputDebug) {
		debug(cout);
	}
}

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
	} else if (x2 != "") {
		currentX = atof(x2.c_str());
	}
	if (y3 != "") {
		currentY = atof(y3.c_str());
	} else if (y2 != "") {
		currentY = atof(y2.c_str());
	}
	if (z3 != "") {
		currentZ = atof(z3.c_str());
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

// adds builtin attributes for atom
bool processAtomBuiltin() {
	Vector v;


	string builtin = getAttribute(currentAtts, "builtin");
	if (builtin == "") {
		cerr << "No builtin attribute for <atom><" << currentElem << ">" << endl;
		return false;
	}
	processBuiltinPCDATA();
	if (currentElem == "float") {
		double value = atof(pcdata.c_str());
	    if (builtin == "x2") {
			currentX = value;
	    } else if (builtin == "y2") {
			currentY = value;
	    } else if (builtin == "x3") {
			currentX = value;
	    } else if (builtin == "y3") {
			currentY = value;
	    } else if (builtin == "z3") {
			currentZ = value;
	    } else {
			cerr << "IGNORED float builtin: " << builtin << endl;
			return false;
	    }
	} else if (currentElem == "integer") {
		int ival = atoi(pcdata.c_str());
	    if (builtin == "formalCharge") {
	        atomPtr->SetFormalCharge(ival);
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

bool endAtomArray() {
	pair<string, OBAtom*> at;
	if (inputCML2) {
		for (int i = 0; i < natoms; i++) {
			OBAtom atom;
			atom.SetAtomicNum(atomicNumVector[i]);
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
//		cerr << "using builtin protocol for bonds..." << endl;
		return false;
	} else if (atomRefs.size() != 2) {
		cerr << "must have 2 atom Refs per bond" << endl;
		return false;
	}
	bondBeginAtom = atomRefs[0];
	bondEndAtom = atomRefs[1];
}

// adds builtin attributes for atom
bool processBondBuiltin() {

	string builtin = getAttribute(currentAtts, "builtin");
	if (builtin == "") {
		cerr << "No builtin attribute for <bond><" << currentElem << ">" << endl;
		return false;
	}
	if (currentElem == "float") {
	    if (false) {
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

bool endBond() {
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
    molPtr->AddBond(*bondPtr);
	return true;
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

bool addString() {
	string title = getAttribute(currentAtts, "title");
	if (title != "") {
		OBPairData *dp = new OBPairData;
		dp->SetAttribute(title);
		dp->SetValue(pcdata);
		molPtr->SetData(dp);
	}
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

	molPtr = &mol;
	dimension = dim;
	WriteMolecule(ofs);
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

bool WriteAtomArray(ostream &ofs) {
	OBAtom *atom;
	vector<OBNodeBase*>::iterator i;

	ofs << "  <atomArray>" << endl;
	int count = 0;
	for (atom = molPtr->BeginAtom(i);atom;atom = molPtr->NextAtom(i)) {
		WriteAtom(ofs, atom, ++count);
	}
	if (outputArray) {
		ofs << "<stringArray builtin=\"atomId\">" << idArray << "</stringArray>" << endl;
		ofs << "<stringArray builtin=\"elementType\">" << elementArray << "</stringArray>" << endl;
		ofs << "<integerArray builtin=\"formalCharge\">" << chargeArray << "</integerArray>" << endl;
		if (strcmp(dimension, "2D")) {
			ofs << "<floatArray builtin=\"x2\">" << x2Array << "</floatArray>" << endl;
			ofs << "<floatArray builtin=\"y2\">" << y2Array << "</floatArray>" << endl;
		} else if (strcmp(dimension, "3D")) {
			ofs << "<floatArray builtin=\"x3\">" << x3Array << "</floatArray>" << endl;
			ofs << "<floatArray builtin=\"y3\">" << y3Array << "</floatArray>" << endl;
			ofs << "<floatArray builtin=\"z3\">" << z3Array << "</floatArray>" << endl;
		}
	}
	ofs << "  </atomArray>" << endl;
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
			if (strcmp(dimension, "2D")) {
				writeAttribute(ofs, "x2", x);
				writeAttribute(ofs, "y2", y);
			} else if (strcmp(dimension, "3D")) {
				writeAttribute(ofs, "x3", x);
				writeAttribute(ofs, "y3", y);
				writeAttribute(ofs, "z3", z);
			}
			ofs << "/>" << endl;
		} else {
			ofs << ">" << endl;
			writeBuiltin(ofs, "elementType", elementType);
			if (charge != 0) writeBuiltin(ofs, "formalCharge", charge);
			if (strcmp(dimension, "2D")) {
				writeBuiltin(ofs, "x2", x);
				writeBuiltin(ofs, "y2", y);
			} else if (strcmp(dimension, "3D")) {
				writeBuiltin(ofs, "x3", x);
				writeBuiltin(ofs, "y3", y);
				writeBuiltin(ofs, "z3", z);
			}
			ofs << "    </atom>" << endl;
		}
	} else {
		appendToArray(idArray, id);
		appendToArray(elementArray, elementType);
		appendToArray(chargeArray, charge);
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
		case 0: break;
		case 1: boChar = "1"; break;
		case 2: boChar = "2"; break;
		case 3: boChar = "3"; break;
		case 5: boChar = "A"; break;
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
			writeAttribute(ofs, "order", boChar);
			ofs << "/>" << endl;
		} else {
			ofs << ">" << endl;
			writeBuiltin(ofs, "atomRef", atomRef1);
			writeBuiltin(ofs, "atomRef", atomRef2);
			writeBuiltin(ofs, "order", boChar);
			ofs << "    </bond>" << endl;
		}
	} else {
		appendToArray(atomRef1Array, atomRef1);
		appendToArray(atomRef2Array, atomRef2);
		appendToArray(orderArray, boChar);
	}

}

}
