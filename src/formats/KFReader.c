 /***************************************************************************
 KFReader.h - Implementation of KFReader library for handling ADF binary file
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
#include "KFc.h"

#include <stdlib.h>
#ifdef _MSC_VER
#include <io.h>
#else
#include <unistd.h>
#endif
#include <stdio.h>
#include <fcntl.h>
#include <string.h>

/* (de-)initializers */

static KFVariable *createKFVariable (KFSection *sec, char* name, int lbl, int dex, int len, int fln, int usd, int typ);
static void deleteKFVariable(KFVariable *var);
static KFSection *createKFSection (KFFile *kf, char *name, int phBlk, int logBlk, int numBlks);
static void deleteKFSection(KFSection *sec);
static KFBlockRun *createKFBlockRun(int phBlk, int logBlk, int numBlks);
static void deleteKFBlockRun(KFBlockRun *br);

/* Functions for KFVariable */
static int getPhysicalBlockNumber(KFVariable *var);
static int getPhysicalBlockNumberForLogical(KFVariable *var, int logBlock);
static int calculateDataOffset(KFVariable *var, void *data, int firstBlock);
static int calculateDataSize(KFVariable *var, void *data, int firstBlock);
static int variablesNamesComparator(const void *varName, const void *var);
static KFVariable *findVariableInSection(const KFSection *sec, const char *name);
int getKFVariableData(KFVariable *var, void *buf);

/* Functions for KFSection */
static void addIndexBlockRun (KFSection *sec, int phBlk, int logBlk, int numBlks);
static void addDataBlockRun (KFSection *sec, int phBlk, int logBlk, int numBlks);
static void parseIndexBlock(KFSection *sec, KFFile *kf, void *buf);
static void parseIndexEntry(KFSection *sec, KFFile *kf, void *buf);
static int sectionsComparator(const void *secName, const void *sec);
static KFSection *findSection(KFFile *kf, const char *name);

/* Functions for KFfile */
static void initialize(KFFile *kf, void *buf);
static void getSectionFromSuperEntry(KFFile *kf, void *buf);
static void getDataBlockFromSuperEntry(KFFile *kf, void *buf);
static KFVariable *findVariable(KFFile *kf, const char *name);

/* Utility functions */
static int readBlock(int fd, void *buf, int block);
static void getSectionFromName(const char *complexName, char *secName);
static void getVariableFromName(const char *complexName, char *varName);
static int guessIntegerSize(void *buf);
static int guessByteOrder(void *buf, int size);
static int verifySuperIndex(void *buf);
static void swapBytesSuperIndexBlock(KFFile *kf, char *buf);
static void swapBytesIndexBlock(KFFile *kf, char *buf);
static void swapNBytes(char *p, int size, int count);

char *KFTypeNames[5] = {"Unknown","Integer","Real","String","Logical"};

/* --------------------------- *
 * Function definitions follow *
 * --------------------------- */

int openKFFile(KFFile *kf, char *name){
   char   buf[KF_BLOCKLENGTH];
   KFSuperIndexBlock32 *buf32 = (KFSuperIndexBlock32 *)buf;
   KFSuperIndexBlock64 *buf64 = (KFSuperIndexBlock64 *)buf;
   int off;

   int lastSBlock = 0;
   int firstSBlock = 1;
   int recordNum = 1;

   memset(kf, 0, sizeof(KFFile));
   kf->name = name;
#ifdef _MSC_VER
   kf->fd = open(kf->name, O_BINARY|O_RDONLY);
#else
   kf->fd = open(kf->name, O_RDONLY);
#endif
   if (kf->fd == -1){
      perror(kf->name);
      return -1;
   }
   kf->open = 1;

   /* First pass to collect and parse only section header blocks */
   do {

      if (readBlock(kf->fd,(void *)buf,recordNum) != KF_BLOCKLENGTH) return -1; 
   
      if (firstSBlock) {
         if (verifySuperIndex(buf)) {
            /* Figure out basic file characteristics - byte order and int size. */
            initialize(kf, buf);
            firstSBlock = 0;
         }
         else {
            close (kf->fd);
            return -1;
         }
      }
      if (kf->byteOrder != getHostByteOrder())
         swapBytesSuperIndexBlock(kf, buf);
      if (kf->integerSize == 4)
         recordNum = (int)buf32->header.nextPhys;
      else if (kf->integerSize == 8)
         recordNum = (int)buf64->header.nextPhys;
      if (recordNum == 1){
         lastSBlock = 1;
      }
      for (off = kf->superIndexHeaderLength; off <= KF_BLOCKLENGTH - kf->superIndexEntryLength; \
                     off += kf->superIndexEntryLength) {
         getSectionFromSuperEntry(kf, (char*)buf+off);
      }
      
   } while (! lastSBlock);

   /* Second pass to collect data blocks for all sections */ 
   lastSBlock = 0;
   firstSBlock = 1;
   recordNum = 1;

   do {

      if (readBlock(kf->fd,(void *)buf,recordNum) != KF_BLOCKLENGTH) return -1; 
   
      if (kf->byteOrder != getHostByteOrder())
         swapBytesSuperIndexBlock(kf, buf);
      if (kf->integerSize == 4)
         recordNum = (int)buf32->header.nextPhys;
      else if (kf->integerSize == 8)
         recordNum = (int)buf64->header.nextPhys;
      if (recordNum == 1){
         lastSBlock = 1;
      }
      for (off = kf->superIndexHeaderLength; off <= KF_BLOCKLENGTH - kf->superIndexEntryLength; \
                     off += kf->superIndexEntryLength) {
         getDataBlockFromSuperEntry(kf, (char*)buf+off);
      }
      
   } while (! lastSBlock);

   return kf->fd;
}

void closeKFFile (KFFile *kf){
   int i;
   for( i=0; i < kf->sections.length; i++){
      deleteKFSection(getArrayListElement(&(kf->sections),i));
   }
   clearArrayList(&(kf->sections));
   if (kf->fd >= 0)
      close(kf->fd);
}

/*

   returns variable length in units of the corresponding type 

*/

int getKFVariableLength(KFFile *kf, const char *name){
   KFVariable *var = findVariable(kf, name);
   if (var != NULL)
      return var->length;
   else
      return -1;
}

/* The following function is a wrapper around getKFVariableData()
   kf - pointer to open KF file
   name - variable name in the form section%variable
   buf - pointer to the receiving buffer that must be large enough to hold all data

   returns the number of bytes read or -1 if an error has occured
*/

int getKFData(KFFile *kf, const char *name, void *buf){
   KFVariable *var = findVariable(kf, name);
   if (var != NULL)
      return getKFVariableData(var, buf);
   else
      return -1;
}

/* 

   returns variable type (one of the T_* macros) or 0 (zero) 

*/

int getKFVariableType(KFFile *kf, const char *name){
   KFVariable *var = findVariable(kf, name);
   if (var != NULL)
      return var->type;
   else
      return 0;
}

/* 
   the memory space pointed to by buf must be
   large enough to hold all data.
   IMPORTANT: the data is converted to native for this platform types

   return the number of bytes read

*/

int getKFVariableData(KFVariable *var, void *buf){
   char blk[KF_BLOCKLENGTH];
   int retValue = 0;
   KFFile *kf = var->section->file;
   int needInversion = (getHostByteOrder() != kf->byteOrder);
   int off, size, i, 
      logicalBlk = var->firstLogBlk,
      firstBlk = 1,
      done = 0,
      targetIndex = 0;

   while (!done) {
      int physBlock = getPhysicalBlockNumberForLogical(var, logicalBlk);
      if (readBlock(kf->fd, blk, physBlock) == KF_BLOCKLENGTH) {
         if (needInversion)
            swapNBytes(blk, kf->integerSize, KF_N_DATATYPES);
           /* where in the data block the data for this variable starts*/
           off = calculateDataOffset(var, blk, firstBlk); 
           size = calculateDataSize(var, blk, firstBlk);  
           if (size > var->usedLen - targetIndex)
               size = var->usedLen - targetIndex;
           /* size = number of elements of this variable in this block*/
         if (needInversion && var->type != KF_T_STRING) 
            swapNBytes(blk+off, kf->typeSize[var->type], size);
         switch (var->type){
            case KF_T_INTEGER:
            case KF_T_LOGICAL:
               if (kf->integerSize == 4){
                  INT32 *base = (INT32 *)(blk+off);
                  for (i = 0; i < size; i++) {
                     ((int*)buf)[targetIndex] = (int)base[i];
                     targetIndex++;
                  }
               }
               else if (kf->integerSize == 8){
                  INT64 *base = (INT64 *)(blk+off);
                  for (i = 0; i < size; i++) {
                     ((int*)buf)[targetIndex] = (int)base[i];
                     targetIndex++;
                  }
               }
               break;
            case KF_T_DOUBLE:
               memcpy((double*)buf+targetIndex, blk+off, size*sizeof(double));
               targetIndex += size;
               break;
            case KF_T_STRING:
               memcpy((char*)buf+targetIndex, blk+off, size);
               targetIndex += size;
               break;
         }
         if (firstBlk) {
            firstBlk = 0;
         }
         if (targetIndex > var->usedLen - 1) {
            done = 1;
            retValue = targetIndex;
         }
         logicalBlk++;
      }
      else {
         fprintf (stderr, "Error reading %s: incomplete block read. Corrupted file?\n",kf->name);
         retValue = -1;
         done = 1;
      }
   }
   return retValue;
}

int getHostByteOrder(){
   int i = 1;
   char *p = (char *)&i;
   if (*p == 0x01)
      return KF_LITTLE_ENDIAN;
   else
      return KF_BIG_ENDIAN;
}

/* Return 1 if the file has a valid superblock. */
int isValidKFFile(char* fname) {

   int n;
   char blk[KF_BLOCKLENGTH];

#ifdef _MSC_VER
   int fd = open(fname, O_BINARY|O_RDONLY);
#else
   int fd = open(fname, O_RDONLY);
#endif

    
   if (fd < 0) return 0;

   n = readBlock(fd, blk, 1);

   if (n != KF_BLOCKLENGTH) return 0;
   
   n = verifySuperIndex(blk);
   close(fd);
   
   return n;
}


/***********************
 *  Static functions:  *
 ***********************/

static int calculateDataOffset(KFVariable *var, void* data, int firstBlock) {
    /*  First read data block header ( 4 integers )
      NOTE: blkInds[0] is left empty because T_* constants have
      values from 1 to 4 */
   int blkInds[5] = {0,0,0,0,0}, i;
   KFFile *kf = var->section->file;
   int off = 0;

   if (kf->integerSize == 4) {
      KFDataBlockHeader32 *header = (KFDataBlockHeader32 *)data;
      off = sizeof(KFDataBlockHeader32);
      blkInds[1] = (int)header->index[0];
      blkInds[2] = (int)header->index[1];
      blkInds[3] = (int)header->index[2];
      blkInds[4] = (int)header->index[3];
   }
   else if (var->section->file->integerSize == 8) {
      KFDataBlockHeader64 *header = (KFDataBlockHeader64 *)data;
      off = sizeof(KFDataBlockHeader64);
      blkInds[1] = (int)header->index[0];
      blkInds[2] = (int)header->index[1];
      blkInds[3] = (int)header->index[2];
      blkInds[4] = (int)header->index[3];
   }
    for (i = 1; i < var->type; i++) {
      off += blkInds[i] * kf->typeSize[i];
    }

    if (firstBlock) {
      off += kf->typeSize[var->type] * (var->firstBlkIndex - 1);
    }
    return off;
}
static int calculateDataSize(KFVariable *var, void *data, int firstBlock) {
   if (firstBlock) {
      return var->firstBlkLen;
   }
   else {
      KFFile *kf = var->section->file;
      if (kf->integerSize == 4) {
         KFDataBlockHeader32 *header = (KFDataBlockHeader32 *)data;
         return (int)header->index[var->type-1];
      }
      else if (var->section->file->integerSize == 8) {
         KFDataBlockHeader64 *header = (KFDataBlockHeader64 *)data;
         return (int)header->index[var->type-1];
      }
   }
   return 0;
}

/* return the number of read bytes */
static int readBlock(int fd, void *buf, int block){
   int lread;
#ifdef _MSC_VER
   _int64 off = (_int64) (block - 1) * KF_BLOCKLENGTH;
   _int64 off2 = _lseeki64(fd, off, SEEK_SET);
   if (off != off2) {
      fprintf (stderr, "Lseek %d, %I64d", fd, off);
      perror("lseek failed");
   }
#else
   off_t off = (off_t) (block - 1) * KF_BLOCKLENGTH;
   off_t off2 = lseek(fd, off, SEEK_SET);
   if (off != off2) {
      fprintf (stderr, "lseek (%d, %ld) returned %ld", fd, (long)off, (long)off2);
      perror("");
   }
#endif
   lread = read(fd, buf, KF_BLOCKLENGTH);
   if (lread != KF_BLOCKLENGTH){
      perror("read failed");
   }
   return lread;
}

static int guessIntegerSize(void *buf){
   KFSuperIndexBlock32 *bl32;
   KFSuperIndexBlock64 *bl64;

   bl32 = (KFSuperIndexBlock32 *)buf;
   bl64 = (KFSuperIndexBlock64 *)buf;
   if (strncmp(bl32->entries[0].name,KF_SIG,KF_SIG_LENGTH) == 0) {
      return 4;
   }
   else if (strncmp(bl64->entries[0].name,KF_SIG,KF_SIG_LENGTH) == 0) {
      return 8;
   }
   else {
      return 0;
   }
}
static int guessByteOrder(void *buf, int size){
   KFSuperIndexBlock32 *bl32 = (KFSuperIndexBlock32 *)buf;
   KFSuperIndexBlock64 *bl64 = (KFSuperIndexBlock64 *)buf;
   char *c;
   if (size == 4) 
      c = (char *)&(bl32->entries[0].physBlk);
   else if (size == 8)
      c = (char *)&(bl64->entries[0].physBlk);

   if (*c == '\0') 
      return KF_BIG_ENDIAN;
   else if (*c == '\001') 
      return KF_LITTLE_ENDIAN;
   else {
      fprintf (stderr, "Corrupted file header\n");
      return 0;
   }
}


#ifdef __SUNPRO_C
static void swapBytes(char *d, int size){
#else
static void __inline swapBytes(char *d, int size){
#endif

   int i,j, halfSize = size>>1;
   for (i = 0, j = size - 1; i < halfSize ; i++, j--){
      char b = *(d+i);
      *(d+i) = *(d+j);
      *(d+j) = b;
   }
}
static void swapNBytes(char *p, int size, int count){
   int cnt;
   for (cnt = 0; cnt < count; cnt++){
      swapBytes(p+cnt*size, size);
   }
}
static void swapBytesSuperIndexBlock(KFFile *kf, char *buf){
   int size = kf->integerSize;
   int headerSize = kf->superIndexHeaderLength;
   int entrySize = kf->superIndexEntryLength;
   char *p;

   /* swap bytes in the header */
   swapNBytes (buf+KF_SIG_LENGTH, size, 4);
   /* swap bytes in the superindex entries */
   for (p = buf+headerSize; p <= buf+(KF_BLOCKLENGTH-entrySize); p+=entrySize)
      swapNBytes (p+KF_SECTION_NAME_LENGTH, size, 4);
}
static void swapBytesIndexBlock(KFFile *kf, char *buf){
   int size = kf->integerSize;
   int headerSize = kf->indexHeaderLength;
   int entrySize = kf->indexEntryLength;
   char *p;

   /* swap bytes in the header */
   swapNBytes (buf+KF_SECTION_NAME_LENGTH, size, 3+KF_N_DATATYPES);
   /* swap bytes in the superindex entries */
   for (p = buf+headerSize; p <= buf+(KF_BLOCKLENGTH-entrySize); p+=entrySize)
      swapNBytes (p+KF_SECTION_NAME_LENGTH, size, 6);
}




/*************** Functions for KFVariable **************/
static KFVariable *createKFVariable (KFSection *sec, char *name, int lbl, int dex, int len, int fln, int usd, int typ){
   KFVariable *var = malloc (sizeof(KFVariable));
   if (var != NULL) {
      var->section = sec;
      strncpy(var->name, name, KF_SECTION_NAME_LENGTH);
      var->name[KF_SECTION_NAME_LENGTH] = '\0';
      var->firstLogBlk = lbl;
      var->firstBlkIndex = dex;
      var->length = len;
      var->firstBlkLen = fln;
      var->usedLen = usd;
      var->type = typ;
      var->multiBlock = ((usd > fln)? 1 : 0);
   }
   return var;
}
static void deleteKFVariable (KFVariable *var){
   free (var);
}
/** Calculates physical block number (record) where the variable starts
 * First block number has number 1.
 *
 * @return block (record) number, one-based.
 */

static int getPhysicalBlockNumber(KFVariable *var){
   int record = 0, i;
   for (i = var->section->dataBlockRuns.length - 1; i >= 0; i--) {
      KFBlockRun *br = getArrayListElement(&(var->section->dataBlockRuns),i);
      if (var->firstLogBlk >= br->logBlock) {
         record = br->physBlock + (var->firstLogBlk - br->logBlock);
         break;
      }
   }
   return record;
}
static int getPhysicalBlockNumberForLogical(KFVariable *var, int logBlock){
   int i;
   for (i = 0; i < var->section->dataBlockRuns.length; i++) {
      KFBlockRun *br = getArrayListElement(&(var->section->dataBlockRuns),i);
      int delta = logBlock - br->logBlock;
      if ( delta >=0 && delta < br->count ) {
         return br->physBlock + delta;
      }
   }
   return 0;
}

/* Functions for KFSection */

static KFSection *createKFSection (KFFile *kf, char *name, int phBlk, int logBlk, int numBlks){
   KFSection *sec = malloc (sizeof(KFSection));
   memset(sec, 0, sizeof(KFSection));
   if (sec != NULL) {
      strncpy(sec->name, name, KF_SECTION_NAME_LENGTH);
      sec->name[KF_SECTION_NAME_LENGTH] = '\0';
      addIndexBlockRun (sec, phBlk, logBlk, numBlks);
   }
   sec->file = kf;
   return sec;
}

static void deleteKFSection(KFSection *sec){
   int i;
   for (i=0; i < sec->dataBlockRuns.length; i++){
      deleteKFBlockRun(getArrayListElement(&(sec->dataBlockRuns),i));
   }
   clearArrayList(&(sec->dataBlockRuns));
   for (i=0; i < sec->indexBlockRuns.length; i++){
      deleteKFBlockRun(getArrayListElement(&(sec->indexBlockRuns),i));
   }
   clearArrayList(&(sec->indexBlockRuns));
   for (i=0; i < sec->variables.length; i++){
      deleteKFVariable(getArrayListElement(&(sec->variables),i));
   }
   clearArrayList(&(sec->variables));
   free(sec);
}

static void addIndexBlockRun (KFSection *sec, int phBlk, int logBlk, int numBlks){
   int insIndex = sec->indexBlockRuns.length;
   int i;
   KFBlockRun *br;
   for (i = 0; i < sec->indexBlockRuns.length; i++) {
      int iLogBlk = ((KFBlockRun *)getArrayListElement(&(sec->indexBlockRuns),i))->logBlock;
      if (logBlk < iLogBlk) {
         insIndex = i;
      }
   }
   br = createKFBlockRun(phBlk, logBlk, numBlks);
   insertArrayListElement(&(sec->indexBlockRuns), br, insIndex);
   sec->totalIndexBlocks += numBlks;
}

static void addDataBlockRun (KFSection *sec, int phBlk, int logBlk, int numBlks){
   int insIndex = sec->dataBlockRuns.length;
   int i, iLogBlk;
   KFBlockRun *br;
   for (i = 0; i < sec->dataBlockRuns.length; i++) {
      br = (KFBlockRun *)getArrayListElement(&(sec->dataBlockRuns),i);
      iLogBlk = br->logBlock;
      if (logBlk < iLogBlk) {
         insIndex = i;
      }
   }
   br = createKFBlockRun(phBlk, logBlk, numBlks);
   insertArrayListElement(&(sec->dataBlockRuns), br, insIndex);
   sec->totalDataBlocks += numBlks;
}

static KFBlockRun *createKFBlockRun(int phBlk, int logBlk, int numBlks){
   KFBlockRun *br = malloc (sizeof(KFBlockRun));
   br->count = numBlks;
   br->logBlock = logBlk;
   br->physBlock = phBlk;
   return br;
}
static void deleteKFBlockRun(KFBlockRun *br){
   free(br);
}

static void parseIndexBlock(KFSection *sec, KFFile *kf, void *buf){
   if (strncmp(sec->name, (char *)buf, KF_SECTION_NAME_LENGTH) == 0){
      int off;
      for (off = kf->indexHeaderLength; off <= KF_BLOCKLENGTH - kf->indexEntryLength; off += kf->indexEntryLength) {
         if (strncmp((char*)buf+off, KF_EMPTY_SIG, KF_SECTION_NAME_LENGTH) != 0)
            parseIndexEntry(sec, kf, (char *)buf+off);
      }
   }
}

static void parseIndexEntry(KFSection *sec, KFFile *kf, void *buf){
   /* Data in the index block must already have correct byte order */
   /* Here buf is a pointer to an entry in the index block */
   KFIndexBlockEntry32 *ent32 = buf;
   KFIndexBlockEntry64 *ent64 = buf;
   int firstLogBlk, firstBlkIndex, length, firstBlkLen, usedLen, type;
   /* Check that the variable does not already exist */
   KFVariable *var = findVariableInSection(sec, buf); 

   if (kf->integerSize == 4) {
      firstLogBlk = (int)ent32->firstLBlk;
      firstBlkIndex = (int)ent32->firstBlkIndx;
      length = (int)ent32->length;
      firstBlkLen = (int)ent32->firstBlkLen;
      usedLen = (int)ent32->usedLen;
      type = (int)ent32->type;
   }
   else if (kf->integerSize == 8) {
      firstLogBlk = (int)ent64->firstLBlk;
      firstBlkIndex = (int)ent64->firstBlkIndx;
      length = (int)ent64->length;
      firstBlkLen = (int)ent64->firstBlkLen;
      usedLen = (int)ent64->usedLen;
      type = (int)ent64->type;
   }
   if (var == NULL) {
      var = createKFVariable(sec, buf, firstLogBlk, firstBlkIndex, length, firstBlkLen, usedLen, type);
      addArrayListElement(&(sec->variables), var);
   }
}

/* Functions for KFFile */
static int verifySuperIndex(void *buf){
   /* return 1 if superindex signature is valid */
   int ret = strncmp((char *)buf, KF_SIG, KF_SIG_LENGTH);

   if ( ret == 0){
      return 1;
   }
   else {
      return 0;
   }
}

static void initialize(KFFile *kf, void *buf){
   kf->integerSize = guessIntegerSize(buf);
   kf->typeSize[KF_T_INTEGER] = kf->typeSize[KF_T_LOGICAL] = kf->integerSize;
   kf->typeSize[KF_T_DOUBLE] = 8;
   kf->typeSize[KF_T_STRING] = 1;
   if (kf->integerSize == 4){
      kf->indexHeaderLength = sizeof(KFIndexBlockHeader32);
      kf->indexEntryLength = sizeof(KFIndexBlockEntry32);
      kf->superIndexHeaderLength = sizeof(KFSuperIndexBlockHeader32);
      kf->superIndexEntryLength = sizeof(KFSuperIndexBlockEntry32);
   }
   else if (kf->integerSize == 8){
      kf->indexHeaderLength = sizeof(KFIndexBlockHeader64);
      kf->indexEntryLength = sizeof(KFIndexBlockEntry64);
      kf->superIndexHeaderLength = sizeof(KFSuperIndexBlockHeader64);
      kf->superIndexEntryLength = sizeof(KFSuperIndexBlockEntry64);
   }
   kf->byteOrder = guessByteOrder(buf, kf->integerSize);

}

static void getSectionFromSuperEntry(KFFile *kf, void *buf){
   char name[KF_SECTION_NAME_LENGTH+1];
   int phBlk, logBlk, numBlks, type;
   int iblk;
   KFSuperIndexBlockEntry32 *ent32 = buf;
   KFSuperIndexBlockEntry64 *ent64 = buf;
   char buf2[KF_BLOCKLENGTH];
   KFSection *section;

   strncpy (name, buf, KF_SECTION_NAME_LENGTH);
   name[KF_SECTION_NAME_LENGTH] = '\0';
   if (kf->integerSize == 4){
      phBlk = (int)ent32->physBlk;
      logBlk = (int)ent32->logicBlk;
      numBlks = (int)ent32->numBlks;
      type = (int)ent32->type;
   }
   else if (kf->integerSize == 8){
      phBlk = (int)ent64->physBlk;
      logBlk = (int)ent64->logicBlk;
      numBlks = (int)ent64->numBlks;
      type = (int)ent64->type;
   }
   switch (type){
      case KF_BT_EMPTY:
      case KF_BT_FREE:
      case KF_BT_SUPERINDEX:
      case KF_BT_DATA:
         break;
      case KF_BT_INDEX:
         /* The entry points to an section index block run
            Add the block run to the section. */
         if (logBlk == 1) {
            /* Create a new section since it's a first index block run */
            section = createKFSection(kf, name, phBlk, logBlk, numBlks);
            addArrayListElement(&(kf->sections), section);
         }
         else {
            /* Add a new index block run to an existing section */
            section = findSection(kf, name);
            addIndexBlockRun(section, phBlk, logBlk, numBlks);
         }

         /* Now we want to read and parse the index block run
             so it will be cached in memory */
         for (iblk = 0; iblk < numBlks; iblk++) {
            int curPhysBlk = phBlk + iblk;
            readBlock(kf->fd, buf2, curPhysBlk);
            if (kf->byteOrder != getHostByteOrder())
               swapBytesIndexBlock(kf, buf2);
            parseIndexBlock(section, kf, buf2);
         }
         break;
      default:
         fprintf (stderr, "Unknown block type flag: %d\n", type);
         break;
   }
}

static void getDataBlockFromSuperEntry(KFFile *kf, void *buf){
   char name[KF_SECTION_NAME_LENGTH+1];
   int phBlk, logBlk, numBlks, type;
   int iblk;
   KFSuperIndexBlockEntry32 *ent32 = buf;
   KFSuperIndexBlockEntry64 *ent64 = buf;
   KFSection *section;

   strncpy (name, buf, KF_SECTION_NAME_LENGTH);
   name[KF_SECTION_NAME_LENGTH] = '\0';
   if (kf->integerSize == 4){
      phBlk = (int)ent32->physBlk;
      logBlk = (int)ent32->logicBlk;
      numBlks = (int)ent32->numBlks;
      type = (int)ent32->type;
   }
   else if (kf->integerSize == 8){
      phBlk = (int)ent64->physBlk;
      logBlk = (int)ent64->logicBlk;
      numBlks = (int)ent64->numBlks;
      type = (int)ent64->type;
   }
   switch (type){
      case KF_BT_EMPTY:
      case KF_BT_FREE:
      case KF_BT_SUPERINDEX:
      case KF_BT_INDEX:
         break;
      case KF_BT_DATA:
         section = findSection(kf, name);
         if (section != NULL) {
            addDataBlockRun(section, phBlk, logBlk, numBlks);
         }
         else {
            fprintf (stderr, "Could not find section %s: data block found before index block\n", name);
         }
         break;
      default:
         fprintf (stderr, "Unknown block type flag: %d\n", type);
         break;
   }
}

static void getSectionFromName(const char *complexName, char *secName){

   char *percent = strchr(complexName, '%');
   if (percent != NULL){
      int len = percent-complexName;
      strncpy(secName, complexName, len);
      secName[len] = '\0';
   }
   else {
      *secName = '\0';
   }
}

static void getVariableFromName(const char *complexName, char *varName){
   char *percent = strchr(complexName, '%');
   if (percent != NULL)
      strncpy(varName, percent+1, KF_SECTION_NAME_LENGTH);
   else
      strncpy(varName, complexName, KF_SECTION_NAME_LENGTH);
   
}

static KFVariable *findVariable(KFFile *kf, const char *name){
   char secName[KF_SECTION_NAME_LENGTH+1];
   char varName[KF_SECTION_NAME_LENGTH+1];
   KFSection *sec;
   KFVariable *var;

   if (strchr(name, '%') == NULL){
      fprintf (stderr, "Incomplete variable specification \"%s\" in findVariable()\n", name);
      return NULL;
   }
   getSectionFromName(name, secName);
   getVariableFromName(name, varName);
   sec = findSection (kf, secName);
   if (sec == NULL)
      return NULL;
   var = findVariableInSection(sec, varName);
   return var;
}

static int sectionsComparator(const void *secName, const void *sec){
   int i;
   char str[KF_SECTION_NAME_LENGTH+1];
   strncpy (str, secName, KF_SECTION_NAME_LENGTH);
   for (i = strlen(str); i < KF_SECTION_NAME_LENGTH; i++)
      str[i] = ' ';
   str[KF_SECTION_NAME_LENGTH] = '\0';
   return strncmp(str, ((KFSection *)sec)->name, KF_SECTION_NAME_LENGTH);
}
static int variablesNamesComparator(const void *varName, const void *var){
   int i;
   char str[KF_SECTION_NAME_LENGTH+1];
   strncpy (str, varName, KF_SECTION_NAME_LENGTH);
   for (i = strlen(str); i < KF_SECTION_NAME_LENGTH; i++)
      str[i] = ' ';
   str[KF_SECTION_NAME_LENGTH] = '\0';
   return strncmp(str, ((KFVariable *)var)->name, KF_SECTION_NAME_LENGTH);
}

static KFVariable *findVariableInSection(const KFSection *sec, const char *name){
   return (KFVariable *)findArrayListElement (&(sec->variables), (void *)name, variablesNamesComparator);
}

static KFSection *findSection(KFFile *kf, const char *name){
   return (KFSection *) findArrayListElement(&(kf->sections),name, sectionsComparator);
}
