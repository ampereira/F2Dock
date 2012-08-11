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
#ifndef _MOL_ROTAMER_H_
#define _MOL_ROTAMER_H_

namespace LibMol{
/**
 * Gets the rotation matrix and translation vector for moving
 * a rotamer from the canonical position into the
 * backbone of residue in atomgrp
 * \param[in]  ag            atomgrp containing residue
 * \param[in]  residue_index index inside iares of residue
 * \param[out] rot_to_bb     rmatrix to rotate residue to
 *                           backbone
 * \param[out] trans_to_bb   tvector to translate residue to
 *                           backbone
 */
void get_detransform(const struct atomgrp *ag, const int residue_index,
		     struct rmatrix *rot_to_bb, struct tvector *trans_to_bb);
};
#endif
