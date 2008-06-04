 /***************************************************************************
 KFc.h - header file for KFReader: library for handling ADF binary file
         format, the so-called KF-files. The current version supports only 
         reading KF files. 
         Even though the binary files are system-dependent, they are 
         converted on the fly and the caller routine gets the data in the
         native format.

 Copyright (C) 2006-2008 by Scientific Computing and Modelling NV.
 For support, contact SCM support (support@scm.com)

 This file is part of the ADF software
 For more information, see <http://www.scm.com>

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation version 2 or 3 of the License.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 SCM owns the intellectual property right for this file and reserves the
 right to distribute it under a license other than LGPL
 ****************************************************************************/
#ifndef _KF_H_
#define _KF_H_


#include "ArrayList.h"



/* variable types */

#define KF_T_INTEGER       1
#define KF_T_DOUBLE        2
#define KF_T_STRING        3
#define KF_T_LOGICAL       4

/* block types */
#define KF_BT_EMPTY        0
#define KF_BT_FREE         1
#define KF_BT_SUPERINDEX   2
#define KF_BT_INDEX        3
#define KF_BT_DATA         4

/* byte order */
#define KF_BIG_ENDIAN      1
#define KF_LITTLE_ENDIAN   2

/* KF file signature */
#define KF_SIG             "SUPERINDEX                      "
#define KF_EMPTY_SIG       "EMPTY                           "
#define KF_SIG_LENGTH      32

/* KF block length */
#define KF_BLOCKLENGTH     4096
#define KF_N_DATATYPES     4


#define KF_SECTION_NAME_LENGTH 32

#define KF_NUM_SUPERINDEX_BLOCK_ENTRIES32 ((KF_BLOCKLENGTH-sizeof(KFSuperIndexBlockHeader32))/sizeof(KFSuperIndexBlockEntry32))
#define KF_NUM_INDEX_BLOCK_ENTRIES32 ((KF_BLOCKLENGTH-sizeof(KFIndexBlockHeader32))/sizeof(KFIndexBlockEntry32))
#define KF_NUM_SUPERINDEX_BLOCK_ENTRIES64 ((KF_BLOCKLENGTH-sizeof(KFSuperIndexBlockHeader64))/sizeof(KFSuperIndexBlockEntry64))
#define KF_NUM_INDEX_BLOCK_ENTRIES64 ((KF_BLOCKLENGTH-sizeof(KFIndexBlockHeader64))/sizeof(KFIndexBlockEntry64))

#ifdef _MSC_VER
   typedef int          INT32;
   typedef _int64       INT64;
#else
   typedef int          INT32;
   typedef long long    INT64;
#endif

/** --------------------------
 *  SuperIndex block header 
 *  -------------------------- */ 
typedef struct _KFSuperIndexBlockHeader32 {
   char name[KF_SIG_LENGTH];
   INT32 mxbl;     /* The total number of blocks in the file */
   INT32 mxsbl;    /* The total number of superblocks in the file */
   INT32 mxsec;    /* The number of sections in the file */
   INT32 nextPhys; /* Physical block number of the next superindex block in the chain */
} KFSuperIndexBlockHeader32;
typedef struct _KFSuperIndexBlockHeader64 {
   char name[KF_SIG_LENGTH];
   INT64 mxbl;     /* The total number of blocks in the file */
   INT64 mxsbl;    /* The total number of superblocks in the file */
   INT64 mxsec;    /* The number of sections in the file */
   INT64 nextPhys; /* Physical block number of the next superindex block in the chain */
} KFSuperIndexBlockHeader64;

/** --------------------------
 *  SuperIndex block entry
 *  -------------------------- */ 
typedef struct _KFSuperIndexBlockEntry32 {
   char name[KF_SECTION_NAME_LENGTH];
   INT32 physBlk;  /* Physical block number */
   INT32 logicBlk; /* Logical block number within the section, starts from one. */
   INT32 numBlks;  /* Number of blocks in this block run */
   INT32 type;     /* Type of this entry, one of the BT_* macros */
} KFSuperIndexBlockEntry32;
typedef struct _KFSuperIndexBlockEntry64 {
   char name[KF_SECTION_NAME_LENGTH];
   INT64 physBlk;  /* Physical block number */
   INT64 logicBlk; /* Logical block number within the section, starts from one. */
   INT64 numBlks;  /* Number of blocks in this block run */
   INT64 type;     /* Type of this entry, one of the BT_* macros */
} KFSuperIndexBlockEntry64;

/** --------------------------
 *  Index block header 
 *  -------------------------- */ 
typedef struct _KFIndexBlockHeader32 {
   char name[KF_SECTION_NAME_LENGTH];
   INT32 mxibl;    /* The number of index blocks in this section */
   INT32 mxdbl;    /* The number of data  blocks in this section */
   INT32 mxbyte;   /* The number of used bytes in the last data block of this section */
   INT32 mxtype[KF_N_DATATYPES];/* The number of used enttries of each type in the last data block of this section */
                   /* One per data type in the order: int, double, char, logical */
} KFIndexBlockHeader32;
typedef struct _KFIndexBlockHeader64 {
   char name[KF_SECTION_NAME_LENGTH];
   INT64 mxibl;    /* The number of index blocks in this section */
   INT64 mxdbl;    /* The number of data  blocks in this section */
   INT64 mxbyte;   /* The number of used bytes in the last data block of this section */
   INT64 mxtype[KF_N_DATATYPES];/* The number of used enttries of each type in the last data block of this section */
                   /* One per data type in the order: int, double, char, logical */
} KFIndexBlockHeader64;

/** --------------------------
 *  Index block entry
 *  -------------------------- */ 
typedef struct _KFIndexBlockEntry32 {
   char name[KF_SECTION_NAME_LENGTH];
   INT32 firstLBlk;     /* The logical block number of the data block where the varialbe starts */
   INT32 firstBlkIndx;  /* The index into the data block where the variable starts. 
                            Each index refers to the part of the data block containing data of its own type. */
   INT32 length;        /* Length of the variable in elements (1 for scalars) */
   INT32 firstBlkLen;   /* The number of elements of the variable contained in its 1st data block */
   INT32 usedLen;       /* The number of elements of the variable that has actually been written in the file */
   INT32 type;          /* The type of the variable */
} KFIndexBlockEntry32;
typedef struct _KFIndexBlockEntry64 {
   char name[KF_SECTION_NAME_LENGTH];
   INT64 firstLBlk;     /* The logical block number of the data block where the varialbe starts */
   INT64 firstBlkIndx;  /* The index into the data block where the variable starts. 
                            Each index refers to the part of the data block containing data of its own type. */
   INT64 length;        /* Length of the variable in elements (1 for scalars) */
   INT64 firstBlkLen;   /* The number of elements of the variable contained in its 1st data block */
   INT64 usedLen;       /* The number of elements of the variable that has actually been written in the file */
   INT64 type;          /* The type of the variable */
} KFIndexBlockEntry64;

/** --------------------------
 *  SuperIndex block 
 *  -------------------------- */ 
typedef struct _KFSuperIndexBlock32 {
   KFSuperIndexBlockHeader32 header;
   KFSuperIndexBlockEntry32 entries[KF_NUM_SUPERINDEX_BLOCK_ENTRIES32];
} KFSuperIndexBlock32;
typedef struct _KFSuperIndexBlock64 {
   KFSuperIndexBlockHeader64 header;
   KFSuperIndexBlockEntry64 entries[KF_NUM_SUPERINDEX_BLOCK_ENTRIES64];
} KFSuperIndexBlock64;

/** --------------------------
 *  Data block header 
 *  -------------------------- */ 
typedef struct _KFDataBlockHeader32 {
   INT32 index[KF_N_DATATYPES]; /* Number of elements in this block, one per data type in the order: int, double, char, logical */
} KFDataBlockHeader32;
typedef struct _KFDataBlockHeader64 {
   INT64 index[KF_N_DATATYPES]; /* Number of elements in this block, one per data type in the order: int, double, char, logical */
} KFDataBlockHeader64;

typedef struct _KFBlockRun {
   int physBlock;
   int logBlock;
   int count;
} KFBlockRun;

typedef struct _KFFile {
   char *name;
   ArrayList sections; /* ArrayList containing KFSection elements */
   int fd;            /* File descriptor if the file is open */
   int byteOrder;     /* BIG_ENDIAN or LITTLE_ENDIAN */
   int integerSize;   /* 4 or 8 */
   int open;          /* 1 if open, otherwise closed */
   int indexHeaderLength;
   int indexEntryLength;
   int superIndexHeaderLength;
   int superIndexEntryLength;
   int typeSize[5];
} KFFile;

typedef struct _KFSection {
   KFFile *file;
   ArrayList indexBlockRuns;  /* Array of KFBlockRun objects containing index */
   ArrayList dataBlockRuns; /* Array of KFBlockRun objects containing data */
   ArrayList variables; /* Array of KFVariable objects for this section */
   int totalDataBlocks;
   int totalIndexBlocks;
   char name[KF_SECTION_NAME_LENGTH+1];
} KFSection;

typedef struct _KFVariable {
   KFSection *section;
   int firstLogBlk;     /* The logical block number of the data block where the varialbe starts */
   int firstBlkIndex;  /* The index into the data block where the variable starts. 
                            Each index refers to the part of the data block containing data of its own type. */
   int length;        /* Length of the variable in elements (1 for scalars) */
   int firstBlkLen;   /* The number of elements of the variable contained in its 1st data block */
   int usedLen;       /* The number of elements of the variable that has actually been written in the file */
   int type;          /* The type of the variable */
   int multiBlock;    /* 0 if this variable occupies more than one block */
   char name[KF_SECTION_NAME_LENGTH+1];
} KFVariable;

#ifdef __cplusplus
extern "C" {
#endif

extern char *KFTypeNames[];

/* Support functions */
int getHostByteOrder();

/* Primary methods for accessing KF file data */

int   openKFFile(KFFile *kf, char *name); /* Return file descriptor or -1 if there was an error opening file */
void  closeKFFile (KFFile *kf);           /* Free memory taken by the file structures and close the file */

int   getKFVariableLength(KFFile *kf, const char *name); /* returns variable length in units 
                                                            of the corresponding type */
int   getKFVariableType(KFFile *kf, const char *name); /* returns variable type (one of the 
                                                          T_* macros) or 0 (zero) */

/* Fills in the memory pointed to by buf with data from the KF file.
   The memory space pointed to by buf must be large enough to hold all data.
   IMPORTANT: the data is converted to a native for this platform type

   kf - pointer to an open KF file
   name - variable name in the form section%variable
   buf - pointer to the receiving buffer that must be large enough to hold all data

   returns the number of elements read or -1 if an error has occured
*/
int  getKFData(KFFile *kf, const char *name, void *buf);


#ifdef __cplusplus
}
#endif



#endif
