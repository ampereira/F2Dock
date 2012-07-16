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
#include <math.h>
#include <libmol/libmol.h>

namespace LibMol{

void quaternion_to_rmatrix(struct quaternion *q, struct rmatrix *U)
{				//This function builds the rotation matrix from the given quartenion.
	double Norm =
	    pow(q->q0, 2) + pow(q->q1, 2) + pow(q->q2, 2) + pow(q->q3, 2);

	U->a11 =
	    (pow(q->q0, 2) + pow(q->q1, 2) - pow(q->q2, 2) -
	     pow(q->q3, 2)) / Norm;
	U->a12 = 2 * ((q->q1 * q->q2 - q->q3 * q->q0)) / Norm;
	U->a13 = 2 * ((q->q1 * q->q3 + q->q2 * q->q0)) / Norm;

	U->a21 = 2 * ((q->q1 * q->q2 + q->q3 * q->q0)) / Norm;
	U->a22 =
	    (pow(q->q0, 2) - pow(q->q1, 2) + pow(q->q2, 2) -
	     pow(q->q3, 2)) / Norm;
	U->a23 = 2 * ((q->q2 * q->q3 - q->q1 * q->q0)) / Norm;

	U->a31 = 2 * ((q->q1 * q->q3 - q->q2 * q->q0)) / Norm;
	U->a32 = 2 * ((q->q2 * q->q3 + q->q1 * q->q0)) / Norm;
	U->a33 =
	    (pow(q->q0, 2) - pow(q->q1, 2) - pow(q->q2, 2) +
	     pow(q->q3, 2)) / Norm;
}

void rmatrix_to_quaternion(struct rmatrix *U, struct quaternion *q)
{				//This function builds the quaternion from a given rotation matrix.
	//Determine q0:
	double S12 = ((U->a12) + (U->a21)) / 2;
	double S13 = ((U->a13) + (U->a31)) / 2;
	double S23 = ((U->a23) + (U->a32)) / 2;
	double A12 = ((U->a12) - (U->a21)) / 2;
	double A13 = ((U->a13) - (U->a31)) / 2;
	double A23 = ((U->a23) - (U->a32)) / 2;
	q->q0 =
	    cbrt(sqrt
		 (pow(A12, 2) * pow(A13, 2) * pow(A23, 2) / S12 / S13 / S23 /
		  8));

	//Find the remaining qi:
	if (q->q0 != 0) {
		q->q1 = -A23 / (2 * (q->q0));
		q->q2 = A13 / (2 * (q->q0));
		q->q3 = -A12 / (2 * (q->q0));
	} else {
		q->q1 = sqrt(S12 * S13 / S23 / 2);
		q->q2 = sqrt(S12 * S23 / S13 / 2);
		q->q3 = sqrt(S13 * S23 / S12 / 2);
	}
}

void quaternion_sum(struct quaternion *qa, struct quaternion *qb,
		     struct quaternion *sum)
{
	sum->q0 = qa->q0 + qb->q0;
	sum->q1 = qa->q1 + qb->q1;
	sum->q2 = qa->q2 + qb->q2;
	sum->q3 = qa->q3 + qb->q3;
}

void quaternion_scaled(struct quaternion *q, double a,
			struct quaternion *q_scaled)
{
	q_scaled->q0 = a * q->q0;
	q_scaled->q1 = a * q->q1;
	q_scaled->q2 = a * q->q2;
	q_scaled->q3 = a * q->q3;
}

void quaternion_product(struct quaternion *qa, struct quaternion *qb,
			 struct quaternion *product)
{
	product->q0 =
	    qa->q0 * qb->q0 - qa->q1 * qb->q1 - qa->q2 * qb->q2 -
	    qa->q3 * qb->q3;
	product->q1 =
	    qa->q0 * qb->q1 + qa->q1 * qb->q0 + qa->q2 * qb->q3 -
	    qa->q3 * qb->q2;
	product->q2 =
	    qa->q0 * qb->q2 + qa->q2 * qb->q0 + qa->q3 * qb->q1 -
	    qa->q1 * qb->q3;
	product->q3 =
	    qa->q0 * qb->q3 + qa->q3 * qb->q0 + qa->q1 * qb->q2 -
	    qa->q2 * qb->q1;
}

void quaternion_conjugate(struct quaternion *q, struct quaternion *qc)
{
	qc->q0 = q->q0;
	qc->q1 = -q->q1;
	qc->q2 = -q->q2;
	qc->q3 = -q->q3;
}

void quaternion_inverse(struct quaternion *q, struct quaternion *q_inv)
{
	struct quaternion qc;
	struct quaternion Norm;
	quaternion_conjugate(q, &qc);
	quaternion_product(q, &qc, &Norm);
	quaternion_scaled(&qc, 1 / (Norm.q0), q_inv);
}

};
