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
#include <errno.h>

#include <libmol/libmol.h>

namespace LibMol{

struct matrix2df* potential_matrix2df_ncontacts_bin (struct atomgrp* agA, struct atomgrp* agB, struct prm* prm, float r1, float r2, int only_sab)
{
	if (r1 < 0.0 || r2 < 0.0)
	{
		fprintf (stderr, "begin error\n");
		fprintf (stderr, "in function potential_matrix2df_ncontacts_bin\n");
		fprintf (stderr, "at least one of the bin limits is less than 0\n");
		fprintf (stderr, "end error\n");
		exit (EXIT_FAILURE);
	}

	int Aatomi, Batomi; // loop iters

	// squared vals for euclidean dist
	float r1sq = powf (r1, 2.0);
	float r2sq = powf (r2, 2.0);

	struct matrix2df* M = matrix2df_create (prm->natoms, prm->natoms);
	matrix2df_init (M, 0); // init all matrix vals to 0

	// loop through every atom in agA
	for (Aatomi = 0; Aatomi < agA->natoms; Aatomi++)
	{
		if (only_sab && ! agA->atoms[Aatomi].sa)
			continue;

		int Atypen = agA->atoms[Aatomi].atom_typen;

		if (Atypen < 0 || Atypen > prm->natoms-1)
		{
			fprintf (stderr, "begin error\n");
			fprintf (stderr, "in function potential_matrix2df_ncontacts_bin\n");
			fprintf (stderr, "in the first atom group\n");
			fprintf (stderr, "atom type number of atom index %d is not defined in the argument atom prm\n", Aatomi);
			fprintf (stderr, "end error\n");
			exit (EXIT_FAILURE);
		}

		float AX = agA->atoms[Aatomi].X;
		float AY = agA->atoms[Aatomi].Y;
		float AZ = agA->atoms[Aatomi].Z;

		// loop through every atom in agB
		for (Batomi = 0; Batomi < agB->natoms; Batomi++)
		{
			if (only_sab && ! agB->atoms[Batomi].sa)
				continue;

			int Btypen = agB->atoms[Batomi].atom_typen;

			if (Atypen < 0 || Atypen > prm->natoms-1)
			{
				fprintf (stderr, "begin error\n");
				fprintf (stderr, "in function potential_matrix2df_ncontacts_bin\n");
				fprintf (stderr, "in the second atom group\n");
				fprintf (stderr, "atom type number of atom index %d is not defined in the argument atom prm\n", Batomi);
				fprintf (stderr, "end error\n");
				exit (EXIT_FAILURE);
			}

			float BX = agB->atoms[Batomi].X;
			float BY = agB->atoms[Batomi].Y;
			float BZ = agB->atoms[Batomi].Z;

			// calculate euclidean distance
			float rsq = (powf ((AX - BX), 2.0) +
					powf ((AY - BY), 2.0) +
					powf ((AZ - BZ), 2.0));

			if (rsq >= r1sq && rsq < r2sq) // atom distance is within the bin
			{
				// (symmetric matrix)
				M->vals[Atypen][Btypen]++;
				M->vals[Btypen][Atypen]++;
			}
		}
	}

	return M;
}

struct matrix2df* potential_matrix2df_rowcol_subatom_join (struct matrix2df* A, struct prm* prm)
{
	// create subatom matrix
	struct matrix2df* B = matrix2df_create (prm->nsubatoms, prm->nsubatoms);
	matrix2df_init (B, 0); // init all matrix vals to 0

	int i, j;
	for (i = 0; i < A->ni; i++)
	{
		if (i < 0 || i > prm->natoms-1)
		{
			fprintf (stderr, "begin error\n");
			fprintf (stderr, "in function potential_matrix2df_rowcol_subatom_join\n");
			fprintf (stderr, "the row %d does not have a corresponding atom type", i);
			fprintf (stderr, "in the argument atom prm\n");
			fprintf (stderr, "end error\n");
			exit (EXIT_FAILURE);
		}

		int subi = prm->atoms[i].subid; // get subatom mapping

		if (subi < 0)
			continue; // ignore this subatom type

		if (subi > prm->nsubatoms-1)
		{
			fprintf (stderr, "begin error\n");
			fprintf (stderr, "in function potential_matrix2df_rowcol_subatom_join\n");
			fprintf (stderr, "in the first atom group\n");
			fprintf (stderr, "the argument atom prm subatom mapping of atom (row) %d", i);
			fprintf (stderr, "is greater than the maximum subatom type index\n");
			fprintf (stderr, "end error\n");
			exit (EXIT_FAILURE);
		}

		for (j = 0; j < A->nj; j++)
		{
			if (j < 0 || j > prm->natoms-1)
			{
				fprintf (stderr, "begin error\n");
				fprintf (stderr, "in function potential_matrix2df_rowcol_subatom_join\n");
				fprintf (stderr, "the column %d does not have a corresponding atom type", j);
				fprintf (stderr, "in the argument atom prm\n");
				fprintf (stderr, "end error\n");
				exit (EXIT_FAILURE);
			}

			int subj = prm->atoms[j].subid; // get subatom mapping

			if (subj < 0)
				continue; // ignore this subatom type

			if (subj > prm->nsubatoms-1)
			{
				fprintf (stderr, "begin error\n");
				fprintf (stderr, "in function potential_matrix2df_rowcol_subatom_join\n");
				fprintf (stderr, "in the first atom group\n");
				fprintf (stderr, "the argument atom prm subatom mapping of atom (column) %d", j);
				fprintf (stderr, "is greater than the maximum subatom type index\n");
				fprintf (stderr, "end error\n");
				exit (EXIT_FAILURE);
			}

			B->vals[subi][subj] += A->vals[i][j];
		}
	}

	return B;
}

};
