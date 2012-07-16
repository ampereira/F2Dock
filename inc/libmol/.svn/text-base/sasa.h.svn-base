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
#ifndef _MOL_SASA_H_
#define _MOL_SASA_H_

namespace LibMol{
/*
   wrapper for baccs: fill sa term in ag
   using a common set of parameters
*/
//void msur (struct atomgrp* ag, struct prm* prm);
// msur_k == -1 => don't msur
void msur (struct atomgrp* ag, struct prm* prm, float msur_k);

/* convert float atomic surface arrays into 0/1 representation (sa in the structure atom) */
/* atoms with zero radii in prm structure are excluded from calculation */
/* r_solv -solvent radius */
/* cont_acc =1 - contact surface
                        =0 - accessible surface */
/* rpr      =1 - use radii from the parameter structure prm
                        =0 - use preset radii */
/* sthresh   if as > sthresh sa=1
             if as<= sthresh sa=0   */
void baccs(struct atomgrp* ag, struct prm* prm,
                float r_solv, short cont_acc, short rpr, float sthresh);
void baccs1(struct atomgrp* ag, int n_at, int* restat,
                double r_solv, short cont_acc, float sthresh);

/* atoms with zero radii in prm structure are excluded from calculation */
/* r_solv -solvent radius */
/* cont_acc =1 - contact surface
			=0 - accessible surface */
/* rpr      =1 - use radii from the parameter structure prm
			=0 - use preset radii */				
/* as - atomic surface (output) */
void accs (struct atomgrp* ag, struct prm* prm, float r_solv, short cont_acc, short rpr, float* as);
void accs1 (struct atomgrp* ag, int n_at, int* restat, double r_solv, short cont_acc, double* as);
void mark_sasa (struct atomgrp* ag, int* sasas);

int* read_sasa (const char* path);

int numbersasas (int* sasas);
};
#endif
