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
#ifndef _MOL_INTCRDS_H_
#define _MOL_INTCRDS_H_

namespace LibMol{

/**
	\file intcrds.h
	This file contains structures and functions
	that are used to represent a molecule as
	a list of bond lengths, bond angles, and
	torsional angles (an internal coordinates
	representation).
*/

/**
	struct zmat represents the
	molecule as a Z-matrix
*/
struct zmat
{
	int natoms; /**< number of atomic units in the molecule */
	struct zdat* zdats; /**< pointer to list of the data in the zmatrix */
}

/**
  data of the zmat
 */
struct zdat
{
	int atom_typen; /**< atom type number */
	int r0; /**< first atom reference */
	int r0val; /**< r0 value (bond length) */
	int r1; /**< second atom reference */
	int r1val; /**< r1 value (bond angle) */
	int r2; /**< third atom reference */
	int r2val; /**< r2 value (torsional angle) */
}

};

#endif
