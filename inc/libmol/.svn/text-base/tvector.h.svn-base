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
#ifndef _MOL_TVECTOR_H_
#define _MOL_TVECTOR_H_

namespace LibMol{
/**
	Holds a translation of
	each dimension.
*/

typedef struct tvector mol_tvector;

struct tvector
{
	float X; /**< translation in the X dimension */
	float Y; /**< translation in the Y dimension */
	float Z; /**< translation in the Z dimension */
};

mol_tvector *
mol_tvector_create ();

void
mol_tvector_destroy (mol_tvector *v);

/**
  initializes elements of v to x, y, z
*/
void
mol_tvector_init (mol_tvector *v, float x, float y, float z);

/**
  copies from src to previously created dst
*/
void
mol_tvector_copy (mol_tvector* src, mol_tvector* dst);

/**
  returns a copy of a
*/
mol_tvector *
mol_tvector_clone (mol_tvector *a);

/**
  a = a+b
*/
void
mol_tvector_add_inplace (mol_tvector *a, mol_tvector *b);

/**
  returns a+b
*/
mol_tvector *
mol_tvector_add (mol_tvector *a, mol_tvector *b);

/**
  normalize in-place (returns 0 if normalization fails, 1 otherwise)
*/
int
mol_tvector_normalize (mol_tvector *a);

/**
	Print the translation vector.
*/
void
mol_tvector_print (mol_tvector* t);

};
#endif // _MOL_TVECTOR_H_
