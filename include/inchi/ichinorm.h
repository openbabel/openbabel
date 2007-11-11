/*
 * International Chemical Identifier (InChI)
 * Version 1
 * Software version 1.02-beta
 * August 23, 2007
 * Developed at NIST
 *
 * The InChI library and programs are free software developed under the
 * auspices of the International Union of Pure and Applied Chemistry (IUPAC);
 * you can redistribute this software and/or modify it under the terms of 
 * the GNU Lesser General Public License as published by the Free Software 
 * Foundation:
 * http://www.opensource.org/licenses/lgpl-license.php
 */


#ifndef __INCHINORM_H__
#define __INCHINORM_H__


#include "mode.h"
#include "ichi_bns.h"


#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

/* main normalization procedure */
int mark_alt_bonds_and_taut_groups ( inp_ATOM *at, inp_ATOM *at_fixed_bonds_out, int num_atoms,
                                     T_GROUP_INFO *t_group_info, INCHI_MODE *inpbTautFlags, INCHI_MODE *inpbTautFlagsDone );

int MarkTautomerGroups( inp_ATOM *at, int num_atoms, T_GROUP_INFO *t_group_info, C_GROUP_INFO *c_group_info,
                        struct BalancedNetworkStructure *pBNS, struct BalancedNetworkData *pBD);
int MarkChargeGroups ( inp_ATOM *at, int num_atoms, C_GROUP_INFO *c_group_info, T_GROUP_INFO *t_group_info,
                      struct BalancedNetworkStructure *pBNS, struct BalancedNetworkData *pBD );
int MarkSaltChargeGroups ( inp_ATOM *at, int num_atoms, S_GROUP_INFO *s_group_info,
                          T_GROUP_INFO *t_group_info, C_GROUP_INFO *c_group_info,
                          struct BalancedNetworkStructure *pBNS, struct BalancedNetworkData *pBD );
int MarkSaltChargeGroups2 ( inp_ATOM *at, int num_atoms, S_GROUP_INFO *s_group_info,
                          T_GROUP_INFO *t_group_info, C_GROUP_INFO *c_group_info,
                          struct BalancedNetworkStructure *pBNS, struct BalancedNetworkData *pBD );
int MergeSaltTautGroups( inp_ATOM *at, int num_atoms, S_GROUP_INFO *s_group_info,
                          T_GROUP_INFO *t_group_info, C_GROUP_INFO *c_group_info,
                          struct BalancedNetworkStructure *pBNS );
int MakeIsotopicHGroup( inp_ATOM *at, int num_atoms, S_GROUP_INFO *s_group_info,
                          T_GROUP_INFO *t_group_info );

int remove_terminal_HDT( int num_atoms, inp_ATOM *at );
int RemoveExcessiveImplicitH( int num_atoms, int num_removed_H, inp_ATOM *at );
int add_DT_to_num_H( int num_atoms, inp_ATOM *at );
int MarkRingSystemsInp( inp_ATOM *at, int num_atoms, int start );
int free_t_group_info( T_GROUP_INFO *t_group_info );
int make_a_copy_of_t_group_info( T_GROUP_INFO *t_group_info, T_GROUP_INFO *t_group_info_orig );
int set_tautomer_iso_sort_keys( T_GROUP_INFO *t_group_info );
int CountTautomerGroups( sp_ATOM *at, int num_atoms, T_GROUP_INFO *t_group_info );
int CountTautomerGroupsInpAt( inp_ATOM *at, int num_atoms, T_GROUP_INFO *t_group_info );
int SortTautomerGroupsAndEndpoints( T_GROUP_INFO *t_group_info, int num_atoms, int num_at_tg, AT_RANK *nRank );
int FillIsotopicAtLinearCT( int num_atoms, sp_ATOM* at, const AT_RANK *nAtomNumber,
                            AT_ISOTOPIC *LinearCTIsotopic, int nMaxLenLinearCTIsotopic, int *pnLenLinearCTIsotopic );

int FillTautLinearCT2( int num_atoms, int num_at_tg, int bIsoTaut,
        const AT_RANK *nRank, const AT_RANK *nAtomNumber, const AT_RANK *nSymmRank,
        const AT_RANK *nRankIso, const AT_RANK *nAtomNumberIso, const AT_RANK *nSymmRankIso,
        AT_TAUTOMER   *LinearCTTautomer, int nMaxLenLinearCTTautomer, int *pnLenLinearCTTautomer,
        AT_ISO_TGROUP *LinearCTIsotopicTautomer, int nMaxLenLinearCTIsotopicTautomer, int *pnLenLinearCTIsotopicTautomer,
        T_GROUP_INFO *t_group_info );


#ifndef INCHI_ALL_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif /* __INCHINORM_H__ */
