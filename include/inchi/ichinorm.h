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


 /**
  * Underivatization, ring-chain tautomerism, OriGAtData edits, etc.
  */

#ifndef _ICHINORM_H_
#define _ICHINORM_H_

#include "mode.h"
#include "ichitaut.h"

#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

/* main normalization procedure */
int mark_alt_bonds_and_taut_groups( struct tagINCHI_CLOCK *ic,
                                    struct tagCANON_GLOBALS *pCG,
                                    inp_ATOM *at,
                                    inp_ATOM *at_fixed_bonds_out,
                                    int num_atoms,
                                    struct tagInchiTime *ulTimeOutTime,
                                    T_GROUP_INFO *t_group_info,
                                    INCHI_MODE *inpbTautFlags,
                                    INCHI_MODE *inpbTautFlagsDone,
                                    int nebend, int *ebend );

int MarkTautomerGroups( struct tagCANON_GLOBALS *pCG, inp_ATOM *at,
                        int num_atoms,
                        T_GROUP_INFO *t_group_info,
                        C_GROUP_INFO *c_group_info,
                        struct BalancedNetworkStructure *pBNS,
                        struct BalancedNetworkData *pBD );

int MarkChargeGroups( struct tagCANON_GLOBALS *pCG,
                      inp_ATOM *at,
                      int num_atoms,
                      C_GROUP_INFO *c_group_info,
                      T_GROUP_INFO *t_group_info,
                      struct BalancedNetworkStructure *pBNS,
                      struct BalancedNetworkData *pBD );
int MarkSaltChargeGroups( struct tagCANON_GLOBALS *pCG,
                          inp_ATOM *at,
                          int num_atoms,
                          S_GROUP_INFO *s_group_info,
                          T_GROUP_INFO *t_group_info,
                          C_GROUP_INFO *c_group_info,
                          struct BalancedNetworkStructure *pBNS,
                          struct BalancedNetworkData *pBD );
int MarkSaltChargeGroups2( struct tagCANON_GLOBALS *pCG,
                           inp_ATOM *at,
                           int num_atoms,
                           S_GROUP_INFO *s_group_info,
                           T_GROUP_INFO *t_group_info,
                           C_GROUP_INFO *c_group_info,
                           struct BalancedNetworkStructure *pBNS,
                           struct BalancedNetworkData *pBD );
int MergeSaltTautGroups( struct tagCANON_GLOBALS *pCG,
                         inp_ATOM *at, int num_atoms,
                         S_GROUP_INFO *s_group_info,
                         T_GROUP_INFO *t_group_info,
                         C_GROUP_INFO *c_group_info,
                         struct BalancedNetworkStructure *pBNS );
int MakeIsotopicHGroup( inp_ATOM *at,
                        int num_atoms,
                        S_GROUP_INFO *s_group_info,
                        T_GROUP_INFO *t_group_info );

int remove_terminal_HDT( int num_atoms,
                         inp_ATOM *at,
                         int bFixTermHChrg );
int RemoveExcessiveImplicitH( int num_atoms,
                              int num_removed_H,
                              inp_ATOM *at );
int add_DT_to_num_H( int num_atoms,
                     inp_ATOM *at );
int MarkRingSystemsInp( inp_ATOM *at,
                        int num_atoms,
                        int start );
int free_t_group_info( T_GROUP_INFO *t_group_info );
int make_a_copy_of_t_group_info( T_GROUP_INFO *t_group_info,
                                 T_GROUP_INFO *t_group_info_orig );
int set_tautomer_iso_sort_keys( T_GROUP_INFO *t_group_info );
int CountTautomerGroups( sp_ATOM *at,
                         int num_atoms,
                         T_GROUP_INFO *t_group_info );
int CountTautomerGroupsInpAt( inp_ATOM *at,
                              int num_atoms,
                              T_GROUP_INFO *t_group_info );
int SortTautomerGroupsAndEndpoints( struct tagCANON_GLOBALS *pCG,
                                    T_GROUP_INFO *t_group_info,
                                    int num_atoms,
                                    int num_at_tg,
                                    AT_RANK *nRank );
int FillIsotopicAtLinearCT( int num_atoms,
                            sp_ATOM* at,
                            const AT_RANK *nAtomNumber,
                            AT_ISOTOPIC *LinearCTIsotopic,
                            int nMaxLenLinearCTIsotopic,
                            int *pnLenLinearCTIsotopic );
int FillTautLinearCT2( struct tagCANON_GLOBALS *pCG,
                       int num_atoms,
                       int num_at_tg,
                       int bIsoTaut,
                       const AT_RANK *nRank,
                       const AT_RANK *nAtomNumber,
                       const AT_RANK *nSymmRank,
                       const AT_RANK *nRankIso,
                       const AT_RANK *nAtomNumberIso,
                       const AT_RANK *nSymmRankIso,
                       AT_TAUTOMER   *LinearCTTautomer,
                       int nMaxLenLinearCTTautomer,
                       int *pnLenLinearCTTautomer,
                       AT_ISO_TGROUP *LinearCTIsotopicTautomer,
                       int nMaxLenLinearCTIsotopicTautomer,
                       int *pnLenLinearCTIsotopicTautomer,
                       T_GROUP_INFO *t_group_info );


#ifndef COMPILE_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif



#endif    /* _ICHINORM_H_ */
