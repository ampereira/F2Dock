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
#ifndef _MOL_POTENTIAL_H_
#define _MOL_POTENTIAL_H_

namespace LibMol{
/** \file potential.h
	This file contains structures and functions
	that are used for generating interaction
	potentials
*/

/**
  finds all atom pairs in agA and agB that are within the distance r1 and r2
  and increments the corresponding atom_typen index in matrix2df
  \param only_sab only consider sa atoms
*/
struct matrix2df* potential_matrix2df_ncontacts_bin (struct atomgrp* agA, struct atomgrp* agB, struct prm* prm, float r1, float r2, int only_sab);

/**
  creates and returns a matrix corresponding to the subatom
  mapping of input matrix A
*/
struct matrix2df* potential_matrix2df_rowcol_subatom_join (struct matrix2df* A, struct prm* prm);
};
#endif
