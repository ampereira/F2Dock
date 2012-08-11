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
#include <string.h>

#include <libmol/libmol.h>

namespace LibMol{

void
mol_version (char** lineptr, size_t *n)
{
    size_t length = strlen(_MOL_VERSION_);

    if (*lineptr == NULL) {
        *lineptr = (char*)_mol_malloc( sizeof(char) * length );
        *n = length;
    } else if( *n < length) {
        *lineptr = (char*)_mol_realloc(*lineptr, sizeof(char) * length );
        *n = length;
    }

	sprintf(*lineptr, "%s", _MOL_VERSION_);
}

void
mol_svn_version (char** lineptr, size_t *n)
{
    size_t length = strlen(_SVN_VERSION_);

    if (*lineptr == NULL) {
        *lineptr = (char*)_mol_malloc( sizeof(char) * length );
        *n = length;
    } else if( *n < length) {
        *lineptr = (char*)_mol_realloc(*lineptr, sizeof(char) * length );
        *n = length;
    }

	sprintf(*lineptr, "%s", _SVN_VERSION_);
}

};
