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
#ifndef _MOL_ENERGY_H_
#define _MOL_ENERGY_H_

namespace LibMol{

/** \file energy.h
	This file contains structures and functions
	that are used for calculating offgrid energies
*/

/**
  find the pairwise potential on the two structures and return it
*/
float pairwise_potential_energy (struct atomgrp* agA, struct atomgrp* agB, struct prm* prm, int only_sab);

/**
  find the coulombic energy the two structures and return it
*/
float coulombic_elec_energy (struct atomgrp* agA, struct atomgrp* agB, struct prm* prm);

/**
  numerical module energy
*/
float nummod_energy (struct atomgrp* agA, struct atomgrp* agB, struct prm* prm);
};
#endif
