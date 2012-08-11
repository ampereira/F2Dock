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
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <libmol/libmol.h>

namespace LibMol{

void mask_atomgrp (struct atomgrp* target, struct atomgrp* mask, float maskr)
{
	int Aatomi, Batomi; // loop iters

	// squared vals for euclidean dist
	float maskrsq = powf (maskr, 2.0);

	// loop through every atom in target
	for (Aatomi = 0; Aatomi < target->natoms; Aatomi++)
	{
		float AX = target->atoms[Aatomi].X;
		float AY = target->atoms[Aatomi].Y;
		float AZ = target->atoms[Aatomi].Z;

		// loop through every atom in mask
		for (Batomi = 0; Batomi < mask->natoms; Batomi++)
		{
				float BX = mask->atoms[Batomi].X;
				float BY = mask->atoms[Batomi].Y;
				float BZ = mask->atoms[Batomi].Z;

			// calculate euclidean distance
			float rsq = (powf ((AX - BX), 2.0) +
					powf ((AY - BY), 2.0) +
					powf ((AZ - BZ), 2.0));

			if (rsq < maskrsq) // atom distance is within the mask radius
			{
				//target->atoms[Aatomi].sa = 0; // mask it with not sa
				target->atoms[Aatomi].mask = 1; // mask it with not sa
			}
		}
	}
}

};
