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
#ifndef _MOL_PROTEIN_H_
#define _MOL_PROTEIN_H_


namespace LibMol{
/**
  extracts the residues from ag (as separated by nitrogens)
  and stores the result in an array of atomgrps
*/
struct atomgrp** extract_nitrogen_residues (struct atomgrp* ag, struct prm* prm);

/**
  printf the contents of ress to stdout
*/
void print_residues (struct atomgrp** ress, struct prm* prm);

/**
  free the memory allocated to the array of atomgrps
*/
void free_atomgrps (struct atomgrp** ress);

/**
  create a new atomgroup around the nlist atoms of ag listed in list (excluded) 
  up to a distance distup
*/
struct atomgrp *around(int nlist, int *list, struct atomgrp* ag, double distup);

};
#endif
