/*
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.07
 * April 30, 2024
 *
 * MIT License
 *
 * Copyright (c) 2024 IUPAC and InChI Trust
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * The InChI library and programs are free software developed under the
 * auspices of the International Union of Pure and Applied Chemistry (IUPAC).
 * Originally developed at NIST.
 * Modifications and additions by IUPAC and the InChI Trust.
 * Some portions of code were developed/changed by external contributors
 * (either contractor or volunteer) which are listed in the file
 * 'External-contributors' included in this distribution.
 *
 * info@inchi-trust.org
 *
 */

#ifndef _UTIL_H_
#define _UTIL_H_

#include "inpdef.h"

#define EL_NUMBER_H  ((U_CHAR)1)
#define EL_NUMBER_B  ((U_CHAR)5)
#define EL_NUMBER_C  ((U_CHAR)6)
#define EL_NUMBER_N  ((U_CHAR)7)
#define EL_NUMBER_O  ((U_CHAR)8)
#define EL_NUMBER_F  ((U_CHAR)9)
#define EL_NUMBER_SI ((U_CHAR)14)
#define EL_NUMBER_P  ((U_CHAR)15)
#define EL_NUMBER_S  ((U_CHAR)16)
#define EL_NUMBER_CL ((U_CHAR)17)
#define EL_NUMBER_GE ((U_CHAR)32)
#define EL_NUMBER_AS ((U_CHAR)33)
#define EL_NUMBER_SE ((U_CHAR)34)
#define EL_NUMBER_BR ((U_CHAR)35)
#define EL_NUMBER_SB ((U_CHAR)51)
#define EL_NUMBER_TE ((U_CHAR)52)
#define EL_NUMBER_I  ((U_CHAR)53)
#define EL_NUMBER_PO ((U_CHAR)84)
#define EL_NUMBER_AT ((U_CHAR)85)

#define EL_NUMBER_ZY ((U_CHAR)119)
#define EL_NUMBER_ZZ ((U_CHAR)120)

#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

/**
 * @brief Get the atomic mass object
 *
 * @param elname Element name
 * @return Returns the atomic mass as integer
 */
int get_atomic_mass(const char *elname);

/**
 * @brief Get the atomic mass from elnum object
 *
 * @param nAtNum Element number
 * @return Returns the atomic mass as integer
 */
int get_atomic_mass_from_elnum(int nAtNum);

/**
 * @brief Get the number of attached hydrogens
 *
 * @param elname Element name
 * @param inp_num_H Input number of H
 * @param num_iso_H Number of isotopic H
 * @param charge Charge
 * @param radical Radical
 * @param chem_bonds_valence Number of chemical valence bonds
 * @param atom_input_valence Number of atom valence
 * @param bAliased Flag indicating whether the atom is aliased??
 * @param bDoNotAddH Flag indicating whether to add hydrogens
 * @param bHasMetalNeighbor Flag indicating whether there is a metal neighbor
 * @return int Number of attached hydrogens
 */
int get_num_H(const char *elname, int inp_num_H, S_CHAR num_iso_H[], int charge, int radical, int chem_bonds_valence, int atom_input_valence, int bAliased, int bDoNotAddH, int bHasMetalNeighbor);

/**
 * @brief Extract charges and radicals from element name
 *
 * @param elname Element name
 * @param pnRadical Number of radicals
 * @param pnCharge Number of charges
 * @return int Number of charges or radicals extracted
 */
int extract_charges_and_radicals(char *elname, int *pnRadical, int *pnCharge);

/**
 * @brief Extract H atoms from element name
 *
 * @param elname Element name
 * @param num_iso_H Number of isotopic H
 * @return int Number of H atoms extracted
 */
int extract_H_atoms(char *elname, S_CHAR num_iso_H[]);

/**
 * @brief Normalize string (remove leading & trailing spaces, replace consecutive spaces with a single space, remove tabs.)
 *
 * @param name Input string to normalize
 * @return int Number of characters in normalized string
 */
int normalize_string(char *name);

/**
 * @brief Read up to any delimiter from the string
 *
 * @param pstring Pointer to input string
 * @param field Pointer to output field
 * @param maxlen Maximum length of the field
 * @param delims Delimiters
 * @return int
 */
int read_upto_delim(char **pstring, char *field, int maxlen, char *delims);

/**
 * @brief Check if a character is in the list of possible delimiters
 * @note  same as isspace if delims is " \t\n\v\f\r" (0x20 and 0x09-0x0D)
 * @param c Character to check
 * @param delims String of delimiter characters
 * @return int 1 if c is a delimiter, 0 otherwise
 */
int is_matching_any_delim(char c, char *delims);

/**
 * @brief Replace non-ASCII characters with '.',
 *
 * @param line Pointer to the input line
 * @return int Number of replacements
 */
int dotify_non_printable_chars(char *line);

/**
 * @brief Trim leading and trailing spaces from a string
 *
 * @param p Pointer to the input string
 * @param nLen Length of the trimmed string
 * @return char* Pointer to the trimmed string
 */
char *lrtrim(char *p, int *nLen);

/**
 * @brief Remove trailing spaces from a string
 *
 * @param p Input string
 */
void remove_trailing_spaces(char *p);

/**
 * @brief Remove one line feed character from the end of a string
 *
 * @param p Input string
 */
void remove_one_lf(char *p);

/**
 * @brief Copies up to maxlen characters INCLUDING end null from source to target.
 *        Fills out the rest of the target with null bytes.
 *        Protected from non-zero-terminated source and overlapped target/source.
 * @note  If source is NULL or maxlen is 0, nothing is copied and 0 is returned.
 * @param target Target string
 * @param source Source string
 * @param maxlen Maximum length to copy
 * @return 1 on success, 0 if target is NULL or maxlen is 0 or source is NULL
 */
int mystrncpy(char *target, const char *source, unsigned maxlen);

/**
 * @brief Reverse a string in place
 *
 * @param p Pointer to the string to reverse
 */
void mystrrev(char *p);

#define ALPHA_BASE 27

/**
 * @brief Compare two memory blocks in a case-insensitive manner.
 *
 * @param p1 Pointer to the first memory block.
 * @param p2 Pointer to the second memory block.
 * @param length Number of bytes to compare.
 * @return int < 0 if p1 < p2, 0 if p1 == p2, > 0 if p1 > p2
 */
int inchi_memicmp(const void *p1, const void *p2, size_t length);

/**
 * @brief Case-insensitive string comparison
 *
 * @param s1 Pointer to first string
 * @param s2 Pointer to second string
 * @return int < 0 if p1 < p2, 0 if p1 == p2, > 0 if p1 > p2
 */
int inchi_stricmp(const char *s1, const char *s2);

/**
 * @brief Set a string to a specified value for a given length
 *
 * @param s Pointer to the string
 * @param val Value to set
 * @param length Length of the string to set
 * @return char* Pointer to the modified string
 */
char *inchi__strnset(char *s, int val, size_t length);

/**
 * @brief Duplicate a string
 *
 * @param string Pointer to the input string
 * @return char* Pointer to the duplicated string
 */
char *inchi__strdup(const char *string);

/**
 * @brief Convert string to long integer
 *
 * @param str Pointer to the input string
 * @param p Pointer to the position where conversion stopped
 * @param base Base for conversion
 * @return long Converted value
 */
long inchi_strtol(const char *str, const char **p, int base);

/**
 * @brief Convert string to double
 *
 * @param str Pointer to the input string
 * @param p Pointer to the position where conversion stopped
 * @return double Converted value
 */
double inchi_strtod(const char *str, const char **p);

/**
 * @brief Checks if an atom is in the list/path
 *
 * @param pathAtom Pointer to the list/path of atoms
 * @param nNextAtom Atom to find
 * @param nPathLen List length
 * @return AT_NUMB* Pointer to the found atom in the list/path, NULL if not found
 */
AT_NUMB *is_in_the_list(AT_NUMB *pathAtom, AT_NUMB nNextAtom, int nPathLen);

/**
 * @brief Checks if an integer is in the list/path
 *
 * @param pathAtom Pointer to the list/path of integers
 * @param nNextAtom Integer to find
 * @param nPathLen Length of the list/path
 * @return int* Pointer to the found integer in the list/path, NULL if not found
 */
int *is_in_the_ilist(int *pathAtom, int nNextAtom, int nPathLen);

/**
 * @brief Checks if one list of integers is inside another list (ilist in ilist2)
 *
 * @param ilist Pointer to the first list of integers
 * @param nlist Length of ilist
 * @param ilist2 Pointer to the second list of integers
 * @param nlist2 Length of ilist2
 * @return int Returns 1 if ilist is inside ilist2, 0 otherwise
 */
int is_ilist_inside(int *ilist, int nlist, int *ilist2, int nlist2);

/**
 * @brief Extract InChI substring embedded into a longer string
 *
 * @param buf Pointer to output buffer
 * @param str Pointer to input string
 * @param slen Length of the input string
 */
void extract_inchi_substring(char **buf, const char *str, size_t slen);

/**
 * @brief Extract AuxInfo substring embedded into a longer string
 *
 * @param buf Pointer to output buffer
 * @param str Pointer to input string
 * @param slen Length of the input string
 */
void extract_auxinfo_substring(char **buf, const char *str, size_t slen);

/**
 * @brief Parse AuxInfostring and get a list of original atom numbers orig[cano_num]
 *
 * @param saux Pointer to AuxInfo string
 * @param orig Pointer to output array of original atom numbers
 * @return int Status code
 */
int extract_orig_nums_from_auxinfo_string(char *saux, int *orig);

/**
 * @brief Parse AuxInfostring and get non-stereo equivalence classes
 *
 * @param saux Pointer to AuxInfo string
 * @param nat Number of atoms
 * @param orig Pointer to array of original atom numbers
 * @param have_eclass_info Pointer to output flag indicating if equivalence class info is present
 * @param eclass Pointer to output array of equivalence classes by canonical atom number
 * @param eclass_by_origs Pointer to output array of equivalence classes by original atom number
 * @return int Status code
 */
int extract_nonstereo_eq_classes_from_auxinfo_string(char *saux, int nat, int *orig, int *have_eclass_info, int *eclass, int *eclass_by_origs);

/**
 * @brief Extract stereo information from InChI string
 *
 * @param sinchi Pointer to InChI string
 * @param nat Number of atoms
 * @param orig Pointer to array of original atom numbers
 * @param at_stereo_mark Pointer to output array of stereo marks for atoms
 * @return int Status code
 */
int extract_stereo_info_from_inchi_string(char *sinchi, int nat, int *orig, int *at_stereo_mark);

/**
 * @brief Extract all backbone bonds from InChI string
 *
 * @param sinchi Pointer to InChI string
 * @param n_all_bkb_orig Pointer to number of all backbone bonds
 * @param orig Pointer to array of original atom numbers
 * @param all_bkb_orig Pointer to output array of all backbone bonds by original atom numbers
 * @return int Status code
 */
int extract_all_backbone_bonds_from_inchi_string(char *sinchi, int *n_all_bkb_orig, int *orig, int *all_bkb_orig);

/**
 * @brief Get the periodic table number object
 *
 * @param elname Pointer to element name
 * @return int Periodic table number
 */
int get_periodic_table_number(const char *elname);

/**
 * @brief Check if an element is a metal
 *
 * @param nPeriodicNum Periodic table number
 * @return int 1 if the element is a metal, 0 otherwise
 */
int is_el_a_metal(int nPeriodicNum);

/**
 * @brief Get reference value of atom valence at given charge
 *
 * @param nPeriodicNum Periodic table number
 * @param charge Charge
 * @param val_num Valence number
 * @return int Reference valence value
 */
int get_el_valence(int nPeriodicNum, int charge, int val_num);

/**
 * @brief Output valence needed to unambiguosly reconstruct bonds
 *
 * @param nPeriodicNum Periodic table number
 * @param charge Charge
 * @param radical Radical
 * @param bonds_valence Bonds valence
 * @param num_H Number of hydrogens
 * @param num_bonds Number of bonds
 * @return int Valence needed to unambiguously reconstruct bonds
 */
int get_unusual_el_valence(int nPeriodicNum, int charge, int radical, int bonds_valence, int num_H, int num_bonds);

/*  Output valence that does not fit any known valences */

/**
 * @brief Output valence that does not fit any known valences
 *
 * @param nPeriodicNum Periodic table number
 * @param charge Charge
 * @param radical Radical
 * @param bonds_valence Bonds valence
 * @param num_H Number of hydrogens
 * @param num_bonds Number of bonds
 * @return int Chemical valence if unusual, 0 otherwise
 */
int detect_unusual_el_valence(int nPeriodicNum, int charge, int radical, int bonds_valence, int num_H, int num_bonds);

/**
 * @brief Output valence needed to unambiguosly reconstruct number of H
 *
 * @param nPeriodicNum Periodic table number
 * @param charge Charge
 * @param radical Radical
 * @param bonds_valence Bonds valence
 * @param actual_bonds_val Actual bonds valence
 * @param num_H Number of hydrogens
 * @param num_bonds Number of bonds
 * @return int Chemical valence needed to unambiguously reconstruct number of H
 */
int needed_unusual_el_valence(int nPeriodicNum, int charge, int radical, int bonds_valence, int actual_bonds_val, int num_H, int num_bonds);

/**
 * @brief Get the el type object
 *
 * @param nPeriodicNum Periodic table number
 * @return int Element type
 */
int get_el_type(int nPeriodicNum);

/**
 * @brief Check if no H addition allowed
 *
 * @param nPeriodicNum Periodic table number
 * @return int Non-zero if no H addition is allowed, zero otherwise
 */
int if_skip_add_H(int nPeriodicNum);

/**
 * @brief Finds chemical symbol for element of given number.
 *
 * @param nAtNum Atom number
 * @param szElement Pointer to chemical symbol
 * @return int Returns 0 if OK and -1 if element was not found.
 */
int get_element_chemical_symbol(int nAtNum, char *szElement);

/**
 * @brief Finds symbol for element of given number. Accounts for (translates)pseudoelements.
 *
 * @param nAtNum Atom number
 * @param szElement Pointer to chemical symbol
 * @return int Returns 0 if OK and -1 if element was not found.
 */
int get_element_or_pseudoelement_symbol(int nAtNum, char *szElement);

/**
 * @brief Creates a string from the removed protons (?)
 *
 * @param nNumRemovedProtons Number of removed protons
 * @param nNumExchgIsotopicH Pointer to number of exchangeable isotopic hydrogens
 * @param nNumRemovedProtonsIsotopic Pointer to array of removed protons
 * @param bIsotopic Flag to handle isotopes
 * @param szRemovedProtons Pointer to the removed protons
 * @param num_removed_iso_H Pointer to the number of removed isotopic hydrogens
 * @return int
 */
int MakeRemovedProtonsString(int nNumRemovedProtons, NUM_H *nNumExchgIsotopicH, NUM_H *nNumRemovedProtonsIsotopic, int bIsotopic, char *szRemovedProtons, int *num_removed_iso_H);

/*
    Ion pairs and fixing bonds
*/

/**
 * @brief Get the number of hydrogens
 *
 * @param at Pointer to atom list
 * @param iat Input atom id
 * @return int Number of total hydrogens atoms
 */
int num_of_H(inp_ATOM *at, int iat);

/**
 * @brief Get the element group of an element. The base element rather than the periodic group is used to aid readability.
 *
 * @param el Element number
 * @return U_CHAR Returns element group  (NitrogenGroup = 7 (EL_NUMBER_N), OxygenGroup = 8 (EL_NUMBER_O), Carbon = 6 (EL_NUMBER_C))
 */
U_CHAR ion_el_group(int el);

/**
 * @brief Check whether an atom has ion neighbors
 *
 * @param at Pointer to atom array
 * @param iat Atom id
 * @param iat_ion_neigh Atom id of ion neighbor
 * @return int Returns 1 for having an ion neighbor, 0 for none
 */
int has_other_ion_neigh(inp_ATOM *at, int iat, int iat_ion_neigh);

/**
 * @brief Check whether an atom has ion neighbors within sphere (breath first search (BFS) up to  r=2)
 *
 * @param at Pointer to atom array
 * @param iat Atom id
 * @param iat_ion_neigh Atom id of ion neighbor
 * @return int Number of ions in sphere
 */
int has_other_ion_in_sphere_2(inp_ATOM *at, int iat, int iat_ion_neigh);

/**
 * @brief Returns the number of non-metal bonds
 *
 * @param at Pointer to atom list
 * @param at_no Atom number
 * @return int Number of number of no metal bonds
 */
int nNoMetalNumBonds(inp_ATOM *at, int at_no);

/**
 * @brief Returns the number of non-metal bond valences
 *
 * @param at Pointer to atom list
 * @param at_no Atom number
 * @return int Number of non-metal bond valences
 */
int nNoMetalBondsValence(inp_ATOM *at, int at_no);

/**
 * @brief Get the index of the first element that is not a metal
 *
 * @param at Pointer to atom list
 * @param at_no Atom number
 * @return int Atom index, -1 if nothing found
 */
int nNoMetalNeighIndex(inp_ATOM *at, int at_no);

/**
 * @brief Get the index of an element that is not a metal excluding a given index
 *
 * @param at Pointer to atom list
 * @param at_no Atom number
 * @param cur_neigh Excluding atom number
 * @return int Atom index, -1 if nothing found
 */
int nNoMetalOtherNeighIndex(inp_ATOM *at, int at_no, int cur_neigh);

/**
 * @brief Get the index of an element that is not a metal excluding a 2 given indexes
 *
 * @param at Pointer to atom list
 * @param at_no Atom number
 * @param cur_neigh Excluding atom number 1
 * @param cur_neigh2 Excluding atom number 2
 * @return int Atom index, -1 if nothing found
 */
int nNoMetalOtherNeighIndex2(inp_ATOM *at, int at_no, int cur_neigh, int cur_neigh2);

/**
 * @brief Gets the number of bond valences to a metal atom
 *
 * @param at Pointer to atom list
 * @param iat Atom number
 * @return int Number of bond valences, -1 if bond to metal order is not well defined
 */
int nBondsValToMetal(inp_ATOM *at, int iat);

/**
 * @brief Gets the number of bond valences
 *
 * @param at Pointer to atom list
 * @param nNumAltBonds Pointer to the number of alternative bonds (4)
 * @param nNumWrongBonds Pointer to the number of wrong bonds
 * @return int The number of bond valences
 */
int nBondsValenceInpAt(const inp_ATOM *at, int *nNumAltBonds, int *nNumWrongBonds);

/**
 * @brief Checks whether a hetero atom may have exchangeable isotopic hydrogens
 *
 * @param atom Pointer to atom list
 * @param iat Atom number
 * @return int 2 if atom is hydrogen (?) ,1 if yes, 0 if no
 */
int bHeteroAtomMayHaveXchgIsoH(inp_ATOM *atom, int iat);

/**
 * @brief Get the endpoint valence object
 *
 * @param el_number Element number
 * @return int 3 if No, 2 if O, S, SE, TE or 0 otherwise
 */
int get_endpoint_valence(U_CHAR el_number);
#if (KETO_ENOL_TAUT == 1)

/**
 * @brief Get the endpoint valence KET object
 *
 * @param el_number Element number
 * @return int 4 if C, 2 if O or 0 otherwise
 */
int get_endpoint_valence_KET(U_CHAR el_number);
#endif

/* Forward declaration */
struct tagCANON_GLOBALS;

/**
 * @brief Frees bit string in canonicalisation data structure
 *
 * @param pCG Canonicalisation data structure
 * @return int 1 if success, 0 if failed
 */
int SetBitFree(struct tagCANON_GLOBALS *pCG);

/**
 * @brief Write coordinate (double) to string
 *
 * @param str Pointer to output string
 * @param x Input double
 */
void WriteCoord(char *str, double x);
extern const int ERR_ELEM;
extern const int nElDataLen;

#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif

#ifndef INCHI_BUILD_PLATFORM

#if defined(_WIN32)

#if defined(_WIN64)
#define INCHI_BUILD_PLATFORM "Windows 64-bit"
#else
#define INCHI_BUILD_PLATFORM "Windows 32-bit"
#endif

#elif defined(__linux__)

#if defined(__x86_64__) || defined(__ppc64__) || defined(__aarch64__) /* djb-rwth: macro added for 64-bit ARM CPUs -- GH issue #10, thanks to Vincent F. Scalfani */
#define INCHI_BUILD_PLATFORM "Linux 64-bit"
#else
#define INCHI_BUILD_PLATFORM "Linux 32-bit"
#endif

#elif defined(__APPLE__)
#define INCHI_BUILD_PLATFORM "OSX"

#else
#define INCHI_BUILD_PLATFORM ""
#endif
#endif

#ifndef INCHI_BUILD_DEBUG
#if defined(_DEBUG)
#define INCHI_BUILD_DEBUG " Debug"
#else
#define INCHI_BUILD_DEBUG ""
#endif
#endif

#ifndef INCHI_SRC_REV
#if defined(_DEBUG)
#define INCHI_SRC_REV "rev. 9b6f1414ebf3+"
#else
#define INCHI_SRC_REV ""
#endif
#endif

#ifndef INCHI_BUILD_COMPILER

#if defined(_MSC_VER)

#if _MSC_VER > 1900
#define INCHI_BUILD_COMPILER "MS VS 2017 or later"
#elif _MSC_VER == 1900
#define INCHI_BUILD_COMPILER "MS VS 2015"
#elif _MSC_VER == 1800
#define INCHI_BUILD_COMPILER "MS VS 2013"
#elif _MSC_VER == 1700
#define INCHI_BUILD_COMPILER "MS VS 2012"
#elif _MSC_VER == 1600
#define INCHI_BUILD_COMPILER "MS VS 2010"
#elif _MSC_VER == 1500
#define INCHI_BUILD_COMPILER "MS VS 2008"
#elif _MSC_VER == 1400
#define INCHI_BUILD_COMPILER "MS VS 2005"
#elif _MSC_VER == 1310
#define INCHI_BUILD_COMPILER "MS VS 2003"
#elif _MSC_VER == 1300
#define INCHI_BUILD_COMPILER "MS VS 2002"
#elif _MSC_VER == 1200
#define INCHI_BUILD_COMPILER "MS VS 6.0"
#else
#define INCHI_BUILD_COMPILER "MS VC++ 5.0 or earlier"
#endif

#else

#if defined(__GNUC__)
#define INCHI_BUILD_COMPILER "gcc " __VERSION__ ""
#else
#define INCHI_BUILD_COMPILER ""
#endif
#endif

#endif

#endif /* _UTIL_H_ */
