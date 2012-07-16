/*
 * cuckoo.h
 *
 * Header file for cuckoo hashing.
 *
 * Permission to use, copy, modify and distribute this software and
 * its documentation for any purpose is hereby granted without fee,
 * provided that due acknoweledgement to the authors are provided and
 * this permission notice appears in all copies of the software.
 * The software is provided "as is". There is no warranty of any kind.
 *
 *
 * Authors:
 *     Rasmus Pagh and Flemming Friche Rodler
 *     BRICS (Basic Research In Computer Science}
 *     Department of Computer Science
 *     University of Aarhus, Denmark
 *     {pagh,ffr}@brics.dk
 * 
 * Date: June 27, 2001.  
*/

#ifndef CUCKOO_H
#define CUCKOO_H 

typedef struct cell {       /* hash table cell type */ 
  int key; 
} celltype;

typedef struct {            /* dictionary type */ 
  int size;                 /* current size */
  int shift;                /* value used for hash function */
  int tablesize;            /* size of hash tables */
  int minsize,meansize;     /* rehash trigger sizes */
  int maxchain;             /* max. iterations in insert */
  struct cell *T1;          /* point to hash table 1*/
  struct cell *T2;          /* point to hash table 2*/
  int a1[3];                /* hash function 1 */
  int a2[3];                /* hash function 2 */
} dict;

typedef dict *dict_ptr;
//#include "ddriver.h"

typedef int boolean;
dict_ptr    construct_dict(int min_size); 
boolean     insertD       (dict_ptr D, int key);
boolean     lookup       (dict_ptr D, int key); 
boolean     delete_key   (dict_ptr D, int key); 
int         keyval       (dict_ptr D, int key);
int         size        (dict_ptr D); 
void        clear        (dict_ptr D, int min_size); 
dict_ptr    destruct_dict(dict_ptr D); 

dict_ptr alloc_dict(int tablesize);
boolean rehash_insert (dict_ptr D, int key);
void rehash(dict_ptr D, int new_size); 

/* The below hash function was found to work well in practice */
/* There is no proof that this is always the case, and there  */
/* may be a better choice of function.                        */

#define inithashcuckoo(a) \
{\
  a[0] = ((int)rand() << 1) + 1;\
  a[1] = ((int)rand() << 1) + 1;\
  a[2] = ((int)rand() << 1) + 1;\
}

#define hashcuckoo(h,a,shift,key)\
{\
  h = (a[0]*key) ^ (a[1]*key) ^ (a[2]*key);\
  h = (unsigned int)h >> shift;\
}

#define DETAIL 5 /* 0 = no output trace, 1 = output trace */ 

#define TRUE 1
#define FALSE 0



#endif 
