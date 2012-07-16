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
#include <string.h>
#include <libmol/libmol.h>

namespace LibMol{

struct atomgrp** extract_nitrogen_residues (struct atomgrp* ag, struct prm* prm)
{
	char nitrogen[2] = {'N','\0'};
	struct atomgrp** ress; // array of residues result

	// if no atoms or first atom is not N, return ag as index 0 and NULL as index 1 of ress
	if (ag->natoms < 1 || strcmp (nitrogen, prm->atoms[ag->atoms[0].atom_typen].typemin) != 0)
	{
		ress = (struct atomgrp**) _mol_malloc (2 * sizeof (struct atomgrp*));
		ress[0] = copy_atomgrp (ag);
		ress[1] = NULL; // flag the end with NULL
		return ress;
	}

	int nresidues = 0;
	int atomi; // atom index
	int natoms = 0; // natoms in the residue
	int maxnatoms = 0; // maximum number of atoms in a residue
	for (atomi = 0; atomi < ag->natoms; atomi++) // first count the number of residues
	{
		natoms++;
		if (strcmp (nitrogen, prm->atoms[ag->atoms[atomi].atom_typen].typemin) == 0)
		{
			nresidues++;
			if (natoms > maxnatoms)
			{
				maxnatoms = natoms-1;
			}
			natoms = 1;
		}
	}
	if (natoms > maxnatoms) // check natoms in last residue
		maxnatoms = natoms;

	if (nresidues == 0) // no res, return ag as index 0 and NULL as index 1 of ress
	{
		ress = (struct atomgrp**) _mol_malloc (2 * sizeof (struct atomgrp*));
		ress[0] = copy_atomgrp (ag);
		ress[1] = NULL; // flag the end with NULL
	}
	else
	{
		ress = (struct atomgrp**) _mol_malloc ((nresidues+1) * sizeof (struct atomgrp*));
		ress[nresidues] = NULL; // flag the end with NULL

		int resn = -1; // residue number
		int resatomn = 0; // residue atom number
		for (atomi = 0; atomi < ag->natoms; atomi++) // separate and store the ress
		{
			if (strcmp (nitrogen, prm->atoms[ag->atoms[atomi].atom_typen].typemin) == 0)
			{
				resn++;
				if (resn != 0)
					ress[resn-1]->natoms = resatomn; // attach the natoms for the prev res
				resatomn = 0; // reset

				ress[resn] = (struct atomgrp*) _mol_calloc (1, sizeof (struct atomgrp));
				ress[resn]->atoms = (struct atom*) _mol_malloc (maxnatoms * sizeof (struct atom));
			}

			copy_atom (&ag->atoms[atomi], &ress[resn]->atoms[resatomn]);

			resatomn++;
		}
		ress[resn]->natoms = resatomn; // attach the natoms for the last res
	}

	return ress;
}

void free_atomgrps (struct atomgrp** ress)
{
	int resn = 0;
	while (ress[resn] != NULL)
	{
		free_atomgrp (ress[resn]);
		resn++;
	}
	free (ress);
}

void print_residues (struct atomgrp** ress, struct prm* prm)
{
	int resi = 0;
	while (ress[resi] != NULL)
	{
		printf ("residue index: %d\n", resi);
		print_atomgrp (ress[resi], prm);
		resi++;
	}
}


struct atomgrp *around(int nlist, int *list, struct atomgrp* ag, double distup)
{
        int i, j;
        int natoms=ag->natoms;
        double x1,y1,z1,x2,y2,z2,d2;
        double du2=distup*distup;
        if(natoms==0)
        {
           fprintf (stderr,
                "around> ERROR: no atoms initially");
                exit (EXIT_FAILURE); 
        }
        int *tmplist=(int*)_mol_malloc(natoms*sizeof(int));

        for(i=0; i<natoms; i++)tmplist[i]=1;
        for(i=0; i<nlist; i++)tmplist[list[i]]=0;

        int nout=0;
        for(i=0; i<nlist; i++)
        {
           x1=ag->atoms[list[i]].X;
           y1=ag->atoms[list[i]].Y;
           z1=ag->atoms[list[i]].Z;
           for(j=0; j<natoms; j++)
           {
              if(tmplist[j]!=1)continue;
              x2=ag->atoms[j].X-x1;
              y2=ag->atoms[j].Y-y1;
              z2=ag->atoms[j].Z-z1;
              d2=x2*x2+y2*y2+z2*z2;
              if(d2<=du2)
              {
                   tmplist[j]=2;
                   nout++;
              }
           }
        }
        struct atomgrp* destag = (struct atomgrp*) _mol_calloc (1, sizeof (struct atomgrp));
        destag->natoms = nout;
        destag->atoms = (struct atom*) _mol_malloc (sizeof (struct atom) * destag->natoms);
        j=0;
        for (i = 0; i < natoms; i++)
        {
                if(tmplist[i]==2)copy_atom (&ag->atoms[i], &destag->atoms[j++]);
        }
        free(tmplist);
        return destag;
}

};
