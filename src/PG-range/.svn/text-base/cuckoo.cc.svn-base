#ifndef CUCKOO
#define CUCKOO


/*
 * cuckoo.c
 *
 * Dictionary code for cuckoo hashing.
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


#include<stdio.h>
#include<stdlib.h>
#include<strings.h>
#include<time.h>
#include<math.h>
#include "cuckoo.h" 
 

//struct CuckooHash{


dict_ptr alloc_dict(int tablesize) {

  dict_ptr D;
  
  D = (dict_ptr) calloc(1,sizeof(dict));
  D->size = 0;
  D->tablesize = tablesize;
  D->meansize = 5*(2*tablesize)/12;
  D->minsize =  (2*tablesize)/5;  
  D->shift = 32 - (int)(log(tablesize)/log(2)+0.5);
  D->maxchain = 4+(int)(4*log(tablesize)/log(2)+0.5);
  if((D->T1 = (celltype *)calloc(tablesize,sizeof(celltype))) == NULL) {
    fprintf(stderr,"Error while allocating mem for T1\n");
    exit(0);
  }
  if((D->T2 = (celltype *)calloc(tablesize,sizeof(celltype))) == NULL) {
    fprintf(stderr,"Error while allocating mem for T2\n");
    exit(0);
  }
  inithashcuckoo(D->a1);
  inithashcuckoo(D->a2);

  return(D);
}

 
/*------insert taylored to rehash-------------------------------*/
boolean rehash_insert (dict_ptr D, int key) 
{ 
  unsigned long hkey;
  int j;
  celltype x,temp;

  x.key = key; 
  for(j = 0; j < D->maxchain; j++) {
    hashcuckoo(hkey,D->a1,D->shift,x.key);
    temp = D->T1[hkey];
    D->T1[hkey] = x;
    if (!temp.key) return TRUE;
    x = temp;

    hashcuckoo(hkey,D->a2,D->shift,x.key);
    temp = D->T2[hkey];
    D->T2[hkey] = x;
    if (!temp.key) return TRUE;
    x = temp;
  }

  bzero(D->T1,D->tablesize * sizeof(celltype));
  bzero(D->T2,D->tablesize * sizeof(celltype));

  inithashcuckoo(D->a1);
  inithashcuckoo(D->a2);

  return FALSE;
} /*rehash_insert*/ 


/*------rehash--------------------------------------------*/
void rehash(dict_ptr D, int new_size) 
{
  dict_ptr D_new;
  int k;

  D_new = alloc_dict(new_size);

  for(k = 0; k < D->tablesize; k++) {
    if ((D->T1[k].key) && (!rehash_insert(D_new,D->T1[k].key))) 
      { k=-1; continue; }
    if ((D->T2[k].key) && (!rehash_insert(D_new,D->T2[k].key)))
      k=-1;
  }
  free(D->T1);
  free(D->T2);

  D_new->size = D->size;
  *D = *D_new;
  free(D_new);
}/*rehash*/


/*------construct_dict---------------------------------*/
dict_ptr construct_dict(int min_size) 
{
  srand(time(NULL));
  return alloc_dict(min_size);
}/*construct_dict*/ 


/*------insert-----------------------------------------*/
boolean insertD (dict_ptr D, int key) 
{ 

  unsigned long h1,h2;
  int j;
  celltype x,temp;
  
  /*If element already in D then replace and return*/
  hashcuckoo(h1,D->a1,D->shift,key);
  if(D->T1[h1].key==key) {
    return FALSE;
  }
  hashcuckoo(h2,D->a2,D->shift,key);
  if(D->T2[h2].key==key) {
    return FALSE;
  }

  /*else insert new element in D*/
  x.key = key; 
  for(j = 0; j < D->maxchain; j++) {
    temp = D->T1[h1];
    D->T1[h1] = x;
    if (!temp.key) {
      D->size++;
      if(D->tablesize < D->size) rehash(D, 2*D->tablesize);
      return TRUE;
    }
    x = temp;
    hashcuckoo(h2,D->a2,D->shift,x.key);
    
    temp = D->T2[h2];
    D->T2[h2] = x;
    if (!temp.key) {
      D->size++;
      if(D->tablesize < D->size) rehash(D, 2*D->tablesize);
      return TRUE;
    }
    x = temp;
    hashcuckoo(h1,D->a1,D->shift,x.key);
  }
 
  /* Forced rehash */
  if(D->size < D->meansize) 
    rehash(D, D->tablesize);
  else {
    rehash(D, 2*D->tablesize);
  }
  insertD(D,x.key);
  return TRUE;
} /*insert*/ 

/*-------lookup--------------------------------------*/
boolean lookup (dict_ptr D, int key) 
{
  unsigned long hkey;

  hashcuckoo(hkey,D->a1,D->shift,key);
  if(D->T1[hkey].key==key) 
    return TRUE;

  hashcuckoo(hkey,D->a2,D->shift,key);
  return (D->T2[hkey].key==key);

}/*lookup*/ 

/*-------delete---------------------------------------*/ 
boolean delete_key (dict_ptr D, int key) 
{ 
  unsigned long hkey;

  hashcuckoo(hkey,D->a1,D->shift,key);
  if(D->T1[hkey].key==key) {
    D->T1[hkey].key = 0;
    D->size--;
    if(D->size < D->minsize) rehash(D,D->tablesize/2);
    return TRUE;
  }
  else {
    hashcuckoo(hkey,D->a2,D->shift,key);
    if(D->T2[hkey].key==key) {
      D->T2[hkey].key = 0;
      D->size--;
      if(D->size < D->minsize) rehash(D,D->tablesize/2);
      return TRUE;
    }
  }
  return FALSE;
}/*delete*/ 


/*-------size-------------------------------------------*/
int size (dict_ptr D) 
{
  return (D->size);
} /*size*/ 


/*-------clear------------------------------------------*/
void clear (dict_ptr D, int min_size)
{
  dict_ptr D_new;

  D_new = construct_dict(min_size);
  free(D->T1);
  free(D->T2);
  *D = *D_new;
}/*clear*/

/*--------destruct_dict-----------------------------------*/
dict_ptr destruct_dict (dict_ptr D)
{
  free(D->T1);
  free(D->T2);
  free(D);
  return(NULL);
} /*destruct_dict*/
//};

#endif
