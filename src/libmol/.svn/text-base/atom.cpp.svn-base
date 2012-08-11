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
**  Structural Bioinformatics Laboratory and  Boston University. Requests to use this software (or                                                                                 
**  modified copies therof) in any other way should be sent to Dima Kozakov,                                                                                                     
**  Department of Biomedical Engineering: "midas@bu.edu".                                                                                                                  
**                                                                                                                                                                               
**---------------------------------------------------------------------------*/

/**-----------------------------------------------------------------------------
/** Modified at Computational Visualization Center (CVC) at University of Texas
**  at Austin to be used as part of F2Dock.
	Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>
**---------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <libmol/libmol.h>

namespace LibMol{

// warning: this function does not handle all of mol_atom
void
mol_atom_create (mol_atom* a, int nbondis)
{
	assert (a != NULL);

	a->nbondis = nbondis;
	if (a->nbondis > 0)
		a->bondis = (int*) _mol_malloc (a->nbondis * sizeof (int));
	else
		a->bondis = NULL;

    a->nbonds = 0;
    a->nangs = 0;
    a->ntors = 0;
    a->nimps = 0;

	return;
}

// warning: this function does not handle all of mol_atom
void
mol_atom_destroy (mol_atom* a)
{
	assert (a != NULL);

	if (a->nbondis > 0)
		free (a->bondis);
	if (a->nbonds > 0)
		free (a->bonds);
	if (a->nangs > 0)
		free (a->angs);
	if (a->ntors > 0)
		free (a->tors);
	if (a->nimps > 0)
		free (a->imps);

	return;
}

// warning: this function does not handle all of mol_atom
void
mol_atom_copy (mol_atom* asrc, mol_atom* adst)
{
	int i;

	assert (asrc != NULL);
	assert (adst != NULL);
	assert (asrc->nbondis == adst->nbondis);
	if (adst->nbondis > 0)
		assert (adst->bondis != NULL);

	adst->ingrp = asrc->ingrp;
	adst->atom_typen = asrc->atom_typen;
	adst->atom_ftypen = asrc->atom_ftypen;

	adst->sa = asrc->sa;
	adst->fixed = asrc->fixed;
	adst->mask = asrc->mask;
	adst->attl = asrc->attl;

	adst->X = asrc->X;
	adst->Y = asrc->Y;
	adst->Z = asrc->Z;

	for (i = 0; i < adst->nbondis; i++)
		adst->bondis[i] = asrc->bondis[i];

	return;
}

void free_atom(struct atom* atm)
{
   free(atm->bonds);
   free(atm->angs);
   free(atm->tors);
   free(atm->imps);
}

void fullcopy_atom (struct atom* src, struct atom* dest)
{
	if (src == NULL || dest == NULL)
	{
		fprintf (stderr, "err in function copy_atom: src or dest arg is NULL\n");
		exit (EXIT_FAILURE);
	}

        dest->ingrp=src->ingrp;
	dest->atom_typen = src->atom_typen;
        dest->atom_ftypen = src->atom_ftypen;
	dest->sa = src->sa;
	dest->fixed = src->fixed;
	dest->attl = src->attl;
	dest->mask = src->mask;

	dest->X = src->X;
	dest->Y = src->Y;
	dest->Z = src->Z;

        dest->GX = src->GX;
        dest->GY = src->GY;
        dest->GZ = src->GZ;

        dest->eps = src->eps;
        dest->rminh = src->rminh;
        dest->eps03 = src->eps03;
        dest->rminh03 = src->rminh03;
        dest->acevolume = src->acevolume;
        dest->chrg = src->chrg;
        dest->B = src->B;
        int i;
        dest->nbonds = src->nbonds;
        dest->bonds=(struct atombond** )_mol_malloc(dest->nbonds*sizeof(struct atombond* ));
        for(i=0; i<dest->nbonds; i++)
                 dest->bonds[i] = src->bonds[i];
        dest->nangs = src->nangs;
        dest->angs=(struct atomangle** )_mol_malloc(dest->nangs*sizeof(struct atomangle*));
        for(i=0; i<dest->nangs; i++)
                 dest->angs[i] = src->angs[i];
        dest->ntors = src->ntors;
        dest->tors=(struct atomtorsion** )_mol_malloc(dest->ntors*sizeof(struct atomtorsion*));
        for(i=0; i<dest->ntors; i++)
                 dest->tors[i] = src->tors[i];
        dest->nimps = src->nimps;
        dest->imps=(struct atomimproper** )_mol_malloc(dest->nimps*sizeof(struct atomimproper*));
        for(i=0; i<dest->nimps; i++)
                 dest->imps[i] = src->imps[i]; 
}

void copy_atom (struct atom* src, struct atom* dest)
{
        if (src == NULL || dest == NULL)
        {
                fprintf (stderr, "err in function copy_atom: src or dest arg is NULL\n");
                exit (EXIT_FAILURE);
        }

        dest->atom_typen = src->atom_typen;
        dest->sa = src->sa;
		dest->attl = src->attl;
		dest->mask = src->mask;
		dest->fixed = src->fixed;

        dest->X = src->X;
        dest->Y = src->Y;
        dest->Z = src->Z;

        /*
        dest->bonds[0] = src->bonds[0];
        dest->bonds[1] = src->bonds[1];
        dest->bonds[2] = src->bonds[2];
        dest->bonds[3] = src->bonds[3];
        */
}

void fprintf_stderr_atom (struct atom* atom, struct prm* prm)
{
	fprintf (stderr, "\tatom type number: %d\n", atom->atom_typen);
	fprintf (stderr, "\tatom type name prefix: %s\n", prm->atoms[atom->atom_typen].typemaj);
	fprintf (stderr, "\tatom type name suffix: %s\n", prm->atoms[atom->atom_typen].typemin);
	fprintf (stderr, "\tcharge: %.2f\n", prm->atoms[atom->atom_typen].q);
	fprintf (stderr, "\tsa: %d\n", atom->sa);
	fprintf (stderr, "\tX: %8.3f\n", atom->X);
	fprintf (stderr, "\tY: %8.3f\n", atom->Y);
	fprintf (stderr, "\tZ: %8.3f\n", atom->Z);
	fprintf (stderr, "\n");
}

void free_springset(struct springset *sprst)
{
        int i;
        for(i=0; i<sprst->nsprings; i++)
              free(sprst->springs[i].laspr);
        if(sprst->springs != NULL)
               free(sprst->springs);
        if(sprst!= NULL)
               free(sprst);
}

};
