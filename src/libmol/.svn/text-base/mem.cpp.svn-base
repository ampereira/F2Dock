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
#include <stdlib.h>
#include <stdio.h>

#include <libmol/libmol.h>

namespace LibMol{

#ifndef _DEBUG_
void*
_mol_malloc (size_t size)
{
    if (size == 0) {
        return NULL;
    }
	void* v = (void*) malloc (size);
	if (v == NULL)
	{
		fprintf (stderr, "insufficient memory request for %zd\n", size);
		exit (EXIT_FAILURE);
	}
	return v;
}

void*
_mol_calloc (size_t nmemb, size_t size)
{
    if (nmemb == 0 || size == 0) {
        return NULL;
    }
	void* v = (void*) calloc (nmemb, size);
	if (v == NULL)
	{
		perror ("calloc"), exit (EXIT_FAILURE);
	}
	return v;
}

void*
_mol_realloc (void* ptr, size_t size)
{
	if (size < 1)
	{
		fprintf (stderr, "warning: _mol_realloc called with size 0, no realloc will occur\n");
		return ptr;
		
	}
	void* v = (void*) realloc (ptr, size);
	if (v == NULL && ptr != NULL)
	{
		perror ("realloc"), exit (EXIT_FAILURE);
	}
	return v;
}
#endif /* _DEBUG_ */

};
