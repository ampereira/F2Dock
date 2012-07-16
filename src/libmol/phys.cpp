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

float* moment_of_inertia (struct atomgrp* ag)
{
	struct tvector* com = center_of_mass (ag);

	float sq_x_com = powf (com->X, 2.0);
	float sq_y_com = powf (com->Y, 2.0);
	float sq_z_com = powf (com->Z, 2.0);

	float sum_sq_x = 0;
	float sum_sq_y = 0;
	float sum_sq_z = 0;

	float sum_mult_xy = 0;
	float sum_mult_xz = 0;
	float sum_mult_yz = 0;

	int i;
	for (i = 0; i < ag->natoms; i++)
	{
		sum_sq_x += powf (ag->atoms[i].X, 2.0);
		sum_sq_y += powf (ag->atoms[i].Y, 2.0);
		sum_sq_z += powf (ag->atoms[i].Z, 2.0);

		sum_mult_xy += ag->atoms[i].X * ag->atoms[i].Y;
		sum_mult_xz += ag->atoms[i].X * ag->atoms[i].Z;
		sum_mult_yz += ag->atoms[i].Y * ag->atoms[i].Z;
	}

	float* moi_matrix = (float*) _mol_malloc (sizeof (float) * 9);

	moi_matrix[0] = sum_sq_y + sum_sq_z - (sq_y_com + sq_z_com) * ag->natoms;
	moi_matrix[1] = -sum_mult_xy + com->X * com->Y * ag->natoms;
	moi_matrix[2] = -sum_mult_xz + com->X * com->Z * ag->natoms;

	moi_matrix[3] = -sum_mult_xy + com->X * com->Y * ag->natoms;
	moi_matrix[4] = sum_sq_x + sum_sq_z - (sq_x_com + sq_z_com) * ag->natoms;
	moi_matrix[5] = -sum_mult_yz + com->Y * com->Z * ag->natoms;

	moi_matrix[6] = -sum_mult_xz + com->X * com->Z * ag->natoms;
	moi_matrix[7] = -sum_mult_yz + com->Y * com->Z * ag->natoms;
	moi_matrix[8] = sum_sq_x + sum_sq_y - (sq_x_com + sq_y_com) * ag->natoms;

	free (com);

	return moi_matrix;
}

void print_moment_of_inertia_matrix (float* moimat)
{
	printf ("%.3f\t%.3f\t%.3f\n", moimat[0], moimat[1], moimat[2]);
	printf ("%.3f\t%.3f\t%.3f\n", moimat[3], moimat[4], moimat[5]);
	printf ("%.3f\t%.3f\t%.3f\n", moimat[6], moimat[7], moimat[8]);
}

};
