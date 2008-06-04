 /***************************************************************************
 ArrayList.c - Implementation of a simple searchable and extendable list type
 
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
#include "ArrayList.h"

#include <stdlib.h>
#include <string.h>

static void extend (ArrayList *array);

void addArrayListElement (ArrayList *array, void *elem){
   if (array->length >= array->allocatedSize)
      extend (array);
   array->data[array->length] = elem;
   array->length++;
}
void insertArrayListElement (ArrayList *array, void *elem, int position){
   int i;
   while (array->length >= array->allocatedSize || position >= array->allocatedSize)
      extend (array);
   for (i = array->length; i > position; i--){
      array->data[i] = array->data[i-1];
   }
   array->data[position] = elem;
   array->length = (array->length <= position)? position+1 : array->length+1 ;
}
void *getArrayListElement (const ArrayList *array, int index){
   if (index < array->length && index >= 0)
      return array->data[index];
   else
      return NULL;
}
void *removeArrayListElement (ArrayList *array, int index){
   if (index < array->length && index >= 0){
        void *elem = array->data[index];
      int i;
      for (i = index; i < array->length; i++){
         array->data[i] = array->data[i+1];
      }
      array->data[array->length] = NULL;
      array->length--;
      return elem;
   }
   return NULL;
}

void *findArrayListElement(const ArrayList *array, const void *item, int (*comparator)(const void *item1, const void *item2)){
   int i;
   for (i = 0; i < array->length; i++){
      int diff = (*comparator)((char*)item, (char*)array->data[i]);
      if (diff == 0){
         return array->data[i];
      }
   }
   return NULL;
}
int findIndexOfArrayListElement(const ArrayList *array, const void *item, int (*comparator)(const void *item1, const void *item2)){
   int i;
   for (i = 0; i < array->length; i++){
      int diff = (*comparator)((char*)item, (char*)array->data[i]);
      if (diff == 0){
         return i;
      }
   }
   return -1;
}
void clearArrayList(ArrayList *array){
   if (array->data != NULL)
      free(array->data);
   array->length = 0;
   array->allocatedSize = 0;
   array->data = NULL;
}

static void extend (ArrayList *array){
   int oldN = array->allocatedSize;
   void **old = array->data;
   if (array->allocatedSize == 0){
      array->allocatedSize = 16;
   }
   else {
      array->allocatedSize *= 2;
   }
   array->data = malloc ((size_t)(array->allocatedSize*sizeof(void *)));
   if (oldN > 0){
      memcpy(array->data, old, (size_t)(oldN*sizeof(void *)));
      free(old);
   }
   memset(array->data + oldN, 0, (size_t)((array->allocatedSize-oldN)*sizeof(void *)));
}


