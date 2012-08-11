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
#ifndef _QUARTERNIONS_H_
#define _QUARTERNIONS_H_

namespace LibMol{

struct quaternion {
	double q0;
	double q1;
	double q2;
	double q3;
};

void quaternion_to_rmatrix(struct quaternion *q, struct rmatrix *U);
void rmatrix_to_quaternion(struct rmatrix *U, struct quaternion *q);
void quaternion_sum(struct quaternion *qa, struct quaternion *qb,
		     struct quaternion *sum);
void quaternion_scaled(struct quaternion *q, double a,
			struct quaternion *q_scaled);
void quaternion_product(struct quaternion *qa, struct quaternion *qb,
			 struct quaternion *product);
void quaternion_conjugate(struct quaternion *q, struct quaternion *qc);
void quaternion_inverse(struct quaternion *q, struct quaternion *q_inv);

};
#endif
