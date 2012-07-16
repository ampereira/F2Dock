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
#ifndef _MOL_BOND_H_
#define _MOL_BOND_H_

namespace LibMol{

typedef struct atombond mol_bond;

/**
  	bond struct
	  a0---a1

	e = k * (l - l0)^2
*/

// deprecate: rename atombond to mol_bond
struct atombond
{
	// deprecated
	struct atom *a0, *a1; /**< atoms involved in the bond */
	// end deprecated

	// indices of atoms that form this bond
	int ai, aj;

	float l; /**< bond length */

	float l0; /**< equilibrium length */
	float k; /**< spring constant */
};


void
mol_bond_create (mol_bond* b);

void
mol_bond_destroy (mol_bond* b);

void
mol_bond_copy (mol_bond* bsrc, mol_bond* bdst);

};
#endif
