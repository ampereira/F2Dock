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
#include <stdio.h>
#include <math.h>
#include <libmol/libmol.h>

namespace LibMol{

// approximately equal (for good epsilon try 5e-7)
int appxeq (double d1, double d2, double epsilon)
{
	if (fabs(d1 - d2) <= (fabs(d1 + d2) * epsilon))
		return 1; // d1 is approximately equal to d2

	return 0;
}
// definitely greater than (for good epsilon try 1e-6)
int defgt (double d1, double d2, double epsilon)
{
	if ((d1 - d2) > (fabs(d1) * epsilon))
		return 1; // d1 is definitely greater than d2

	return 0;
}
// definitely less than (for good epsilon try 1e-6)
int deflt (double d1, double d2, double epsilon)
{
	if ((d2 - d1) > (fabs(d1) * epsilon))
		return 1; // d1 is definitely less than d2

	return 0;
}
};
