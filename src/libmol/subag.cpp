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
#include <libmol/libmol.h>

namespace LibMol{

void setup_subag(struct atomgrp *ag, struct atomgrp *new_actives, int nlist,
		 int *list)
{
	for (int i = 0; i < nlist; i++) {
		int j = list[i];
		ag->atoms[j].X = new_actives->atoms[i].X;
		ag->atoms[j].Y = new_actives->atoms[i].Y;
		ag->atoms[j].Z = new_actives->atoms[i].Z;
	}
}

void unsetup_subag(struct atomgrp *ag, struct atomgrp *new_actives, int nlist,
		   int *list)
{
	for (int i = 0; i < nlist; i++) {
		int j = list[i];
		new_actives->atoms[i].X = ag->atoms[j].X;
		new_actives->atoms[i].Y = ag->atoms[j].Y;
		new_actives->atoms[i].Z = ag->atoms[j].Z;
	}
}

};
