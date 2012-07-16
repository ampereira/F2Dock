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

float pairwise_potential_energy (struct atomgrp* agA, struct atomgrp* agB, struct prm* prm, int only_sab)
{
	if (prm->pwpot->r1 < 0.0 || prm->pwpot->r2 < 0.0)
	{
		fprintf (stderr, "begin error\n");
		fprintf (stderr, "at least one of the potential's bin limits is less than 0\n");
		fprintf (stderr, "end error\n");
		exit (EXIT_FAILURE);
	}

	int Aatomi, Batomi; // loop iters

	// squared vals for euclidean dist
	float r1sq = powf (prm->pwpot->r1, 2.0);
	float r2sq = powf (prm->pwpot->r2, 2.0);

	float E = 0;

	// loop through every atom in agA
	for (Aatomi = 0; Aatomi < agA->natoms; Aatomi++)
	{
		if (only_sab && ! agA->atoms[Aatomi].sa)
			continue;

		int Atypen = agA->atoms[Aatomi].atom_typen;
		if (Atypen < 0 || Atypen > prm->natoms-1)
		{
			fprintf (stderr, "begin error\n");
			fprintf (stderr, "atom type number of atom index %d is not defined in the argument atom prm\n", Aatomi);
			fprintf (stderr, "end error\n");
			exit (EXIT_FAILURE);
		}

		int subA = prm->atoms[Atypen].subid; // get subatom mapping
		if (subA < 0)
			continue; // ignore this subatom type
		if (subA > prm->nsubatoms-1)
		{
			fprintf (stderr, "begin error\n");
			fprintf (stderr, "the argument atom prm subatom mapping of atom %d", Atypen);
			fprintf (stderr, "is greater than the maximum subatom type index\n");
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
				fprintf (stderr, "atom type number of atom index %d is not defined in the argument atom prm\n", Batomi);
				fprintf (stderr, "end error\n");
				exit (EXIT_FAILURE);
			}

			int subB = prm->atoms[Btypen].subid; // get subatom mapping
			if (subB < 0)
				continue; // ignore this subatom type
			if (subB > prm->nsubatoms-1)
			{
				fprintf (stderr, "begin error\n");
				fprintf (stderr, "the argument atom prm subatom mapping of atom %d", Btypen);
				fprintf (stderr, "is greater than the maximum subatom type index\n");
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
				for (int ek = 0; ek < prm->pwpot->k; ek++)
				{
					E += prm->pwpot->lambdas[ek] * prm->pwpot->Xs[(ek*prm->nsubatoms)+subA] * prm->pwpot->Xs[(ek*prm->nsubatoms)+subB];
				}
			}
		}
	}

	return E;
}

float coulombic_elec_energy (struct atomgrp* agA, struct atomgrp* agB, struct prm* prm)
{
	int Aatomi, Batomi; // loop iters

	float E = 0;

	float maxr = 20.0;
	float maxrsq = maxr * maxr;
	float minr = 2.0;
	float minrsq = minr * minr;

	// loop through every atom in agA
	for (Aatomi = 0; Aatomi < agA->natoms; Aatomi++)
	{
		int Atypen = agA->atoms[Aatomi].atom_typen;
		if (Atypen < 0 || Atypen > prm->natoms-1)
		{
			fprintf (stderr, "atom type number of atom index %d is not defined in atom prm\n", Aatomi);
			exit (EXIT_FAILURE);
		}

		float q1 = prm->atoms[Atypen].q;

		float AX = agA->atoms[Aatomi].X;
		float AY = agA->atoms[Aatomi].Y;
		float AZ = agA->atoms[Aatomi].Z;

		// loop through every atom in agB
		for (Batomi = 0; Batomi < agB->natoms; Batomi++)
		{
			int Btypen = agB->atoms[Batomi].atom_typen;
			if (Btypen < 0 || Btypen > prm->natoms-1)
			{
				fprintf (stderr, "atom type number of atom index %d is not defined in the atom prm\n", Batomi);
				exit (EXIT_FAILURE);
			}

			float q2 = prm->atoms[Btypen].q;

			float BX = agB->atoms[Batomi].X;
			float BY = agB->atoms[Batomi].Y;
			float BZ = agB->atoms[Batomi].Z;

			// calculate euclidean distance
			float rsq = (powf ((AX - BX), 2.0) +
					powf ((AY - BY), 2.0) +
					powf ((AZ - BZ), 2.0));

			//E += (q1 * q2) / rsq;
			if (rsq <= maxrsq)
			{
				if (rsq < minrsq) // set a limit
					rsq = minrsq;

				E += (q1 * q2) * ( (1/rsq) - ( (1/maxrsq) * (2-(rsq/maxrsq)) ) );
			}
		}
	}

	return E;
}

float nummod_energy (struct atomgrp* agA, struct atomgrp* agB, struct prm* prm)
{
	int Aatomi, Batomi; // loop iters

	float vrE = 0;
	float vaE = 0;
	float eE = 0;

	float r2sq = powf (6.5, 2.0);

	// loop through every atom in agA
	for (Aatomi = 0; Aatomi < agA->natoms; Aatomi++)
	{
		int Atypen = agA->atoms[Aatomi].atom_typen;
		if (Atypen < 0 || Atypen > prm->natoms-1)
		{
			fprintf (stderr, "atom type number of atom index %d is not defined in atom prm\n", Aatomi);
			exit (EXIT_FAILURE);
		}

		float r1sq = powf (prm->atoms[Atypen].r, 2.0);
		float q1 = prm->atoms[Atypen].q;

		float AX = agA->atoms[Aatomi].X;
		float AY = agA->atoms[Aatomi].Y;
		float AZ = agA->atoms[Aatomi].Z;

		// loop through every atom in agB
		for (Batomi = 0; Batomi < agB->natoms; Batomi++)
		{
			int Btypen = agB->atoms[Batomi].atom_typen;
			if (Btypen < 0 || Btypen > prm->natoms-1)
			{
				fprintf (stderr, "atom type number of atom index %d is not defined in the atom prm\n", Batomi);
				exit (EXIT_FAILURE);
			}

			//float r2 = prm->atoms[Btypen].r;
			float q2 = prm->atoms[Btypen].q;

			float BX = agB->atoms[Batomi].X;
			float BY = agB->atoms[Batomi].Y;
			float BZ = agB->atoms[Batomi].Z;

			// calculate euclidean distance
			float rsq = (powf ((AX - BX), 2.0) +
					powf ((AY - BY), 2.0) +
					powf ((AZ - BZ), 2.0));

			if (rsq < r1sq)
			{
				//vrE += vw * 1.0;
				vrE += 1.0;
			}
			else if (rsq < r2sq)
			{
				//vaE -= vw * 1.0;
				vaE -= 1.0;
			}

			//E += ew * ((q1 * q2) / rsq);
			eE += ((q1 * q2) / rsq);
		}
	}

	printf ("vrE: %.3f\n", vrE);
	printf ("vaE: %.3f\n", vaE);
	printf ("eE: %.3f\n", eE);

	float E = 0.0;

	return E;
}

};
