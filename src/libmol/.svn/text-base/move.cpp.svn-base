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
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <libmol/libmol.h>

namespace LibMol{

float calc_distance (struct tvector* va, struct tvector* vb)
{
	float d = sqrtf (
			(va->X - vb->X) * (va->X - vb->X) +
			(va->Y - vb->Y) * (va->Y - vb->Y) +
			(va->Z - vb->Z) * (va->Z - vb->Z) );

	if ( errno )
	{
		perror("sqrtf");
		exit (EXIT_FAILURE);
	}

	return d;
}

struct tvector* center_of_mass (struct atomgrp* ag)
{
	// sums of atom positions
	float sumX = 0;
	float sumY = 0;
	float sumZ = 0;

	int i;
	for (i = 0; i < ag->natoms; i++) // sum the atoms
	{
		sumX += ag->atoms[i].X;
		sumY += ag->atoms[i].Y;
		sumZ += ag->atoms[i].Z;
	}

	struct tvector* tvec = (struct tvector*) _mol_malloc (sizeof (struct tvector));

	tvec->X = (sumX / ag->natoms);
	tvec->Y = (sumY / ag->natoms);
	tvec->Z = (sumZ / ag->natoms);

	return tvec;
}

struct tvector* negative_center_of_mass (struct atomgrp* ag)
{
	struct tvector* tvec = center_of_mass (ag);

	// negate the tvec
	tvec->X = -tvec->X;
	tvec->Y = -tvec->Y;
	tvec->Z = -tvec->Z;

	return tvec;
}

struct tvector* center_of_extrema (struct atomgrp* ag)
{
	if (ag->natoms == 0)
		return NULL;

	// extrema of atom positions
	float minX = ag->atoms[0].X, maxX = ag->atoms[0].X;
	float minY = ag->atoms[0].Y, maxY = ag->atoms[0].Y;
	float minZ = ag->atoms[0].Z, maxZ = ag->atoms[0].Z;

	int i;
	for (i = 1; i < ag->natoms; i++)
	{
		if (ag->atoms[i].X < minX)
			minX = ag->atoms[i].X;
		if (ag->atoms[i].X > maxX)
			maxX = ag->atoms[i].X;

		if (ag->atoms[i].Y < minY)
			minY = ag->atoms[i].Y;
		if (ag->atoms[i].Y > maxY)
			maxY = ag->atoms[i].Y;

		if (ag->atoms[i].Z < minZ)
			minZ = ag->atoms[i].Z;
		if (ag->atoms[i].Z > maxZ)
			maxZ = ag->atoms[i].Z;
	}

	struct tvector* tvec = (struct tvector*) _mol_malloc (sizeof (struct tvector));

	tvec->X = ((maxX - minX) / 2) + minX;
	tvec->Y = ((maxY - minY) / 2) + minY;
	tvec->Z = ((maxZ - minZ) / 2) + minZ;

	return tvec;
}

struct tvector* negative_center_of_extrema (struct atomgrp* ag)
{
	struct tvector* tvec = center_of_extrema (ag);

	// negate the tvec
	tvec->X = -tvec->X;
	tvec->Y = -tvec->Y;
	tvec->Z = -tvec->Z;

	return tvec;
}

struct tvector* negative_tvector (struct tvector* tv)
{
	struct tvector* neg_tv = (struct tvector*) _mol_malloc (sizeof (struct tvector));

	neg_tv->X = -tv->X;
	neg_tv->Y = -tv->Y;
	neg_tv->Z = -tv->Z;

	return neg_tv;
}

void center_atomgrp_on_center_of_extrema (struct atomgrp* ag, struct atomgrp* centered_ag)
{
	struct tvector* neg_com = negative_center_of_extrema (ag);
	translate_atomgrp (ag, centered_ag, neg_com);
	free (neg_com);
}

void center_atomgrp_on_center_of_mass (struct atomgrp* ag, struct atomgrp* centered_ag)
{
	struct tvector* neg_com = negative_center_of_mass (ag);
	translate_atomgrp (ag, centered_ag, neg_com);
	free (neg_com);
}

void translate_atomgrp (struct atomgrp* src, struct atomgrp* dest, struct tvector* tvec)
{
	if (src == NULL || dest == NULL)
	{
		fprintf (stderr, "src and dest pointers must be malloced\n");
		exit (EXIT_FAILURE);
	}

	int i;
	for (i = 0; i < src->natoms; i++)
	{
		dest->atoms[i].X = src->atoms[i].X + tvec->X;
		dest->atoms[i].Y = src->atoms[i].Y + tvec->Y;
		dest->atoms[i].Z = src->atoms[i].Z + tvec->Z;
	}
}

struct rotset* read_rotset (const char* path)
{
	FILE* fp = myfopen (path, "r"); // open parms file

	char* line = NULL;
	size_t len = 0;

	int nrots = 0; // number of sets of rots
	while (getline (&line, &len, fp) != -1)
	{
		if (! iswhiteline (line))
		{
			nrots++;
		}
	}
	rewind (fp); // reset the file pointer

	struct rotset* rotset = (struct rotset*) _mol_malloc (sizeof (struct rotset));
	rotset->nrots = nrots; // save the value
	rotset->rmats = (struct rmatrix*) _mol_malloc (rotset->nrots * sizeof (struct rmatrix));

	int i = 0;
	while (getline (&line, &len, fp) != -1)
	{
		int tmpi; // tmp roti
		if (iswhiteline (line))
		{
			continue;
		}
		if (sscanf (line,	"%d %f %f %f %f %f %f %f %f %f",
							&tmpi,
							&rotset->rmats[i].a11, &rotset->rmats[i].a12, &rotset->rmats[i].a13,
							&rotset->rmats[i].a21, &rotset->rmats[i].a22, &rotset->rmats[i].a23,
							&rotset->rmats[i].a31, &rotset->rmats[i].a32, &rotset->rmats[i].a33
							) < 10)
		{
			fprintf (stderr, "read error on file %s\n", path);
			fprintf (stderr, "expecting one index and 9 floats per line %d\n", i);
			exit (EXIT_FAILURE);
		}
		if (tmpi != i)
		{
			fprintf (stderr, "read error on file %s\n", path);
			fprintf (stderr, "expecting index %d, found index %d\n", i, tmpi);
			exit (EXIT_FAILURE);
		}
		i++;
	}
	if (line)
		free (line);
	myfclose (fp);

	return rotset;
}

struct axisset* read_axisset (const char* path)
{
	FILE* fp = myfopen (path, "r"); // open parms file

	char* line = NULL;
	size_t len = 0;

	int naxes = 0; // number of sets of rots
	while (getline (&line, &len, fp) != -1)
	{
		if (! iswhiteline (line))
		{
			naxes++;
		}
	}
	rewind (fp); // reset the file pointer

	struct axisset* axisset = (struct axisset*) _mol_malloc (sizeof (struct axisset));
	axisset->naxes = naxes; // save the value
	axisset->axes = (struct axis*) _mol_malloc (axisset->naxes * sizeof (struct axis));

	int i = 0;
	while (getline (&line, &len, fp) != -1)
	{
		int tmpi; // tmp roti
		if (iswhiteline (line))
		{
			continue;
		}
		if (sscanf (line,	"%d %f %f %f",
							&tmpi,
							&axisset->axes[i].a,
							&axisset->axes[i].b,
							&axisset->axes[i].c
							) < 4)
		{
			fprintf (stderr, "read error on file %s\n", path);
			fprintf (stderr, "expecting one index and 3 floats per line %d\n", i);
			exit (EXIT_FAILURE);
		}
		if (tmpi != i)
		{
			fprintf (stderr, "read error on file %s\n", path);
			fprintf (stderr, "expecting index %d, found index %d\n", i, tmpi);
			exit (EXIT_FAILURE);
		}
		i++;
	}
	if (line)
		free (line);
	myfclose (fp);

	return axisset;
}

void print_axisset (struct axisset* axisset)
{
	int i;
	printf ("number of axes: %d\n", axisset->naxes);
	for (i = 0; i < axisset->naxes; i++)
	{
		printf ("%d:\t", i);
		printf ("%8.4f\t%8.4f\t%8.4f\n",
					axisset->axes[i].a,
					axisset->axes[i].b,
					axisset->axes[i].c);
	}
}

struct rotset* create_sym_rotset (int nmonomers)
{
	struct rotset* rotset = (struct rotset*) _mol_malloc (sizeof (struct rotset));
	rotset->nrots = 2; // 2 rotations: 1 for cyclical symmetry and 1 for dihedral symmetry
	rotset->rmats = (struct rmatrix*) _mol_malloc (rotset->nrots * sizeof (struct rmatrix));

	float angle = 2.0 * 3.14159 / (float) nmonomers;

	// for cyclic symmetry:
	rotset->rmats[0].a11 = cos (angle);
	rotset->rmats[0].a12 = sin (angle);
	rotset->rmats[0].a13 = 0.0;
	rotset->rmats[0].a21 = -sin (angle);
	rotset->rmats[0].a22 = cos (angle);
	rotset->rmats[0].a23 = 0.0;
	rotset->rmats[0].a31 = 0.0;
	rotset->rmats[0].a32 = 0.0;
	rotset->rmats[0].a33 = 1.0;

	// for dihedral symmetry:

	return rotset;
}

void free_rotset (struct rotset* rotset)
{
	free (rotset->rmats);
	free (rotset);
}

void rotate_atomgrp (struct atomgrp* src, struct atomgrp* dest, struct rmatrix* rmatrix)
{
	if (src == NULL || dest == NULL)
	{
		fprintf (stderr, "src and dest pointers must be malloced\n");
		exit (EXIT_FAILURE);
	}

	int atomi;
	for (atomi = 0; atomi < src->natoms; atomi++)
	{
		double X = src->atoms[atomi].X;
		double Y = src->atoms[atomi].Y;
		double Z = src->atoms[atomi].Z;
		dest->atoms[atomi].X = rmatrix->a11*X + rmatrix->a12*Y + rmatrix->a13*Z;
		dest->atoms[atomi].Y = rmatrix->a21*X + rmatrix->a22*Y + rmatrix->a23*Z;
		dest->atoms[atomi].Z = rmatrix->a31*X + rmatrix->a32*Y + rmatrix->a33*Z;
	}
}

void move_atomgrp (struct atomgrp* ag, struct atomgrp* moved_ag, struct rmatrix* rmatrix, struct tvector* tv)
{
	rotate_atomgrp (ag, moved_ag, rmatrix);
	translate_atomgrp (moved_ag, moved_ag, tv);
}

void print_rotset (struct rotset* rotset)
{
	int i;
	printf ("number of rotations: %d\n", rotset->nrots);
	for (i = 0; i < rotset->nrots; i++)
	{
		printf ("%d:\t", i);
		printf ("%8.4f\t%8.4f\t%8.4f\t%8.4f\t%8.4f\t%8.4f\t%8.4f\t%8.4f\t%8.4f\n",
					rotset->rmats[i].a11,
					rotset->rmats[i].a12,
					rotset->rmats[i].a13,
					rotset->rmats[i].a21,
					rotset->rmats[i].a22,
					rotset->rmats[i].a23,
					rotset->rmats[i].a31,
					rotset->rmats[i].a32,
					rotset->rmats[i].a33);
	}
}

void print_rmatrix (struct rmatrix* rmatrix)
{
	printf ("rmatrix:");
	printf ("%8.4f\t%8.4f\t%8.4f\t%8.4f\t%8.4f\t%8.4f\t%8.4f\t%8.4f\t%8.4f\n",
				rmatrix->a11,
				rmatrix->a12,
				rmatrix->a13,
				rmatrix->a21,
				rmatrix->a22,
				rmatrix->a23,
				rmatrix->a31,
				rmatrix->a32,
				rmatrix->a33);
}

void print_tvector (struct tvector* t)
{
	printf ("tvector: %.4f, %.4f, %.4f\n", t->X, t->Y, t->Z);
}

};
