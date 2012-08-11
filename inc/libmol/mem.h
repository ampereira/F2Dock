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
#ifndef _MOL_MEM_H_
#define _MOL_MEM_H_

namespace LibMol{

#ifdef _DEBUG_
#define _mol_malloc(size) malloc(size)
#define _mol_realloc(ptr, size) realloc(ptr, size)
#define _mol_calloc(nmemb, size) calloc(nmemb, size)
#endif

/**
	This is a wrapper to malloc that will print
	an error message if malloc returns NULL.
	\param size allocate this many bytes of memory
	\return a pointer to the allocated memory
*/
#ifndef _DEBUG_
void*
_mol_malloc (size_t size);
#endif

/**
	This is a wrapper to calloc that will print
	an error message if malloc returns NULL.
	\param nmemb size of array
	\param size size of each element in the array
	\return a pointer to the allocated memory
*/
#ifndef _DEBUG_
void*
_mol_calloc (size_t nmemb, size_t size);
#endif

/**
	This is a wrapper to realloc that will print
	an error message accordingly.
	\param ptr pointer to the original memory
	\param size allocate this many bytes of memory
	\return a pointer to the allocated memory
*/
#ifndef _DEBUG_
void*
_mol_realloc (void* ptr, size_t size);
#endif

};
#endif
