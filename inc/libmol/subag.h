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
#ifndef _MOL_SUBAG_H_
#define _MOL_SUBAG_H_

namespace LibMol{
/**
 * Replaces atoms in a larger atomgrp with the atoms
 * from another atomgrp.  The first item in list should
 * correspond to the first atom in new_actives, and so on.
 * @param[out] ag          The atomgrp that has its atoms replaced
 * @param[in]  new_actives The source of atoms to be placed
 * @param[in]  nlist       The number of atoms to be placed
 * @param[in]  list        The indices of the replaced atoms in ag
 */
void setup_subag(struct atomgrp *ag, struct atomgrp *new_actives, int nlist, int *list);

/**
 * The inverse of seutp_subag.  Takes atoms from a larger atomgrp
 * and places them in another atomgrp.  The first item in list should
 * correspond to the first atom in new_actives, and so on.
 * @param[in]  ag          The source of atoms to be placed
 * @param[out] new_actives The atomgrp that hs its atoms replaced
 * @param[in]  nlist       The number of atoms to be placed
 * @param[in]  list        The indices of the atoms in ag to be copied to new_actives
 */
void unsetup_subag(struct atomgrp *ag, struct atomgrp *new_actives, int nlist, int *list);
};
#endif
