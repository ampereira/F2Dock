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

void get_detransform(const struct atomgrp *ag, const int residue_index,
		     struct rmatrix *rot_to_bb, struct tvector *trans_to_bb)
{
	struct tvector e1, e2, e3;
	float d;

	struct atom *nl = &(ag->atoms[ag->iares[residue_index]]);
	struct atom *cal = &(ag->atoms[ag->iares[residue_index] + 2]);
	struct atom *cl = &(ag->atoms[ag->iares[residue_index + 1] - 2]);	//second to last atom in res 

	e1.X = nl->X - cal->X;
	e1.Y = nl->Y - cal->Y;
	e1.Z = nl->Z - cal->Z;

	//normalize e1
	d = sqrtf((e1.X * e1.X) + (e1.Y * e1.Y) + (e1.Z * e1.Z));
	e1.X /= d;
	e1.Y /= d;
	e1.Z /= d;

	e2.X = cl->X - cal->X;
	e2.Y = cl->Y - cal->Y;
	e2.Z = cl->Z - cal->Z;

	//Orthogonalize:
	//dot product to get component of e1 along e2
	d = (e1.X * e2.X) + (e1.Y * e2.Y) + (e1.Z * e2.Z);
	printf("d = %f\n", d);
	//subtract component
	e2.X -= e1.X * d;
	e2.Y -= e1.Y * d;
	e2.Z -= e1.Z * d;

	//normalize e2
	d = sqrtf((e2.X * e2.X) + (e2.Y * e2.Y) + (e2.Z * e2.Z));
	e2.X /= d;
	e2.Y /= d;
	e2.Z /= d;

	//generate third orthonormal component as cross product
	e3.X = (e1.Y * e2.Z) - (e1.Z * e2.Y);
	e3.Y = (e1.Z * e2.X) - (e1.X * e2.Z);
	e3.Z = (e1.X * e2.Y) - (e1.Y * e2.X);

	//assign rotation matrix and translation vector
	rot_to_bb->a11 = e1.X;
	rot_to_bb->a12 = e2.X;
	rot_to_bb->a13 = e3.X;
	rot_to_bb->a21 = e1.Y;
	rot_to_bb->a22 = e2.Y;
	rot_to_bb->a23 = e3.Y;
	rot_to_bb->a31 = e1.Z;
	rot_to_bb->a32 = e2.Z;
	rot_to_bb->a33 = e3.Z;

	trans_to_bb->X = cal->X;
	trans_to_bb->Y = cal->Y;
	trans_to_bb->Z = cal->Z;
}

};
