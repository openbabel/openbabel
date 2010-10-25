 /***************************************************************************
 ArrayList.h - header file for ArrayList: a simple searchable and extendable
               list type

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

#ifndef _ARRAY_LIST_H_
#define _ARRAY_LIST_H_

typedef struct _ArrayList {
   void **data;
   int allocatedSize;
   int length;
} ArrayList;


void addArrayListElement (ArrayList *array, void *elem);
void insertArrayListElement (ArrayList *array, void *elem, int position);
void *getArrayListElement (const ArrayList *array, int index);
void *removeArrayListElement (ArrayList *array, int index);
void clearArrayList(ArrayList *array);
void *findArrayListElement(const ArrayList *array, const void *searchItem, \
                           int (*comparator)(const void *searchItem, const void *arrayItem));
int findIndexOfArrayListElement(const ArrayList *array, const void *searchItem, \
                           int (*comparator)(const void *searchItem, const void *arrayItem));

#endif

