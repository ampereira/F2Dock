/**-----------------------------------------------------------------------------                                                                                                  
**                                                                                                                                                                               
**  Copyright (C) : Structural Bioinformatics Laboratory, Boston University.                                                                                                                        
**                                                                                                                                                                               
**  This software was developed at the Boston University 2006-2011, by                                                                                                      
**  Structural Bioinformatics Laboratory, as part of NIH funded research.                                                                                                                      
**                                                                                                                                                                               
**  Explicit permission is hereby granted to US Universities and US                                                                                                     
**  Government supported Research Institutions to copy and modify this                                                                                                           
**  software for educational and research purposes, provided copies include                                                                                                      
**  this notice. This software (or modified copies thereof) may not be                                                                                                           
**  distributed to any other institution without express permission from the                                                                                                     
**  Structural Bioinformatics Laboratory and  Boston University. Requests to use this software (or                                                                                 **  modified copies therof) in any other way should be sent to Dima Kozakov,                                                                                                     
**  Department of Biomedical Engineering: "midas@bu.edu".                                                                                                                  
**                                                                                                                                                                               
**---------------------------------------------------------------------------*/
#ifndef _MOL_MYHELPERS_H_
#define _MOL_MYHELPERS_H_

namespace LibMol{
#ifdef freeMem
   #undef freeMem
#endif
#define freeMem( ptr ) { if ( ptr != NULL ) free( ptr ); }

#ifdef _DEBUG_
#define mymalloc(size) malloc(size)
#define myrealloc(ptr, size) realloc(ptr, size)
#define mycalloc(nmemb, size) calloc(nmemb, size)
#endif

/** \file myhelpers.h
	This file contains wrapper functions
	and supporting string functions.
*/

/**
	This is a wrapper to malloc that will print
	an error message if malloc returns NULL.
	\param size allocate this many bytes of memory
	\return a pointer to the allocated memory
*/
#ifndef _DEBUG_
void* mymalloc(size_t size);
#endif

/**
	This is a wrapper to calloc that will print
	an error message if malloc returns NULL.
	\param nmemb size of array
	\param size size of each element in the array
	\return a pointer to the allocated memory
*/
#ifndef _DEBUG_
void* mycalloc(size_t nmemb, size_t size);
#endif

/**
	This is a wrapper to realloc that will print
	an error message accordingly.
	\param ptr pointer to the original memory
	\param size allocate this many bytes of memory
	\return a pointer to the allocated memory
*/
#ifndef _DEBUG_
void* myrealloc(void* ptr, size_t size);
#endif

#ifdef _DARWIN_

int  getline (char **lineptr, size_t *n, FILE *stream);

#endif /* _DARWIN_ */


void myexit (int status);

/**
	This is a wrapper to fopen that will print
	an error message and exit if fopen returns NULL.
*/
FILE* myfopen (const char* path, const char* mode);

/**
	This is a wrapper to fclose that will print
	an error message and exit if an error occurs.
*/
void myfclose (FILE* fp);

/**
	Returns 1 if line consists of whitespace
	characters, else 0;
*/
int iswhiteline (char* line);

/**
	Returns 1 if line consists of whitespace
	characters, else 0;
*/
char* rstrip (char* string);

/**
	Prints a formatted error message.
	[ Rezaul Chowdhury, Nov. 17, 2010 ]	
*/
void print_error( char *format, ... );


struct list
{
int n;
int *list;
};
struct pointlist
{
int n;
int** list;
};

};
#endif
