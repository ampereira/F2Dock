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
#include <limits.h>

#include <libmol/libmol.h>

namespace LibMol{

// identifiers in the parameter file:
#define VERSION_KEY   "version"
#define ATOM_KEY      "atom"
#define VDW_KEY       "vdw"
#define PWPOT_KEY "pwpot"
#define PWPOT_R_KEY "pwpot.r"
#define BOND_KEY "bond"

#define MAXLEN 100 // max string length

// atype comparison function
static int comp_prmatom (const void* a1, const void *a2)
{
	struct prmatom* atom1 = (struct prmatom*) a1; 
	struct prmatom* atom2 = (struct prmatom*) a2;

	if (strcmp (atom1->typemaj, atom2->typemaj) < 0)
		return -1;
	if (strcmp (atom1->typemaj, atom2->typemaj) > 0)
		return 1;
	return (strcmp (atom1->typemin, atom2->typemin));
}

static int comp_prmbond (const void* b1, const void *b2)
{
	struct prmbond* bond1 = (struct prmbond*) b1; 
	struct prmbond* bond2 = (struct prmbond*) b2;

	if (bond1->i < bond2->i)
		return -1;
	if (bond1->i > bond2->i)
		return 1;
	// (i's are equal)
	if (bond1->j < bond2->j)
		return -1;
	if (bond1->j > bond2->j)
		return 1;

	fprintf (stderr, "error: identical %s subatom mappings\n", BOND_KEY);
	exit (EXIT_FAILURE);

	return 0;
}

struct prm* read_prm (const char* path, const char* bin_version)
{
	struct prm* prm = (struct prm*) _mol_malloc (sizeof (struct prm));

	read_prmversion (path, bin_version);
	read_prmatom (prm, path);
	read_prmpwpot (prm, path);
	read_prmbond (prm, path);

	return prm;
}

void read_prmversion (const char* path, const char* bin_version)
{
	char idstr[MAXLEN];
	char prm_version[MAXLEN];
	int versionf = 0;
	char* line = NULL;
	size_t len = 0;

	FILE* fp = myfopen (path, "r"); // open prm file

	while (getline (&line, &len, fp) != -1)
	{
		if (sscanf (line, "%80s", idstr) < 1) // get id
			continue; // whitespace line

		if (strncmp (idstr, VERSION_KEY, MAXLEN) == 0) // cmp id
		{
			if (sscanf (line, "%*s %s", prm_version) < 1)
			{
				fprintf (stderr, "in file %s:\n", path);
				fprintf (stderr, "wrong format for %s line\n", VERSION_KEY);
				exit (EXIT_FAILURE);
			}

			if (strncmp (prm_version, bin_version, MAXLEN) != 0)
			{
				fprintf (stderr, "in file %s:\n", path);
				fprintf (stderr, "prm version (%s) inconsistent with binary version (%s)\n", prm_version, bin_version);
				exit (EXIT_FAILURE);
			}

			versionf = 1;
		}
	}
	if (line)
		free (line);
	myfclose (fp);

	if (! versionf)
	{
		fprintf (stderr, "in file %s:\n", path);
		fprintf (stderr, "version specifier not found\n");
		exit (EXIT_FAILURE);
	}
}

void read_prmatom (struct prm* prm, const char* path)
{
	prm->atoms = NULL;
	prm->natoms = 0;
	prm->nsubatoms = -1;
	int atomsi = 0;
	char idstr[MAXLEN];

	char* line = NULL;
	size_t len = 0;

	FILE* fp = myfopen (path, "r"); // open prm file
	
	// preparse
	while (getline (&line, &len, fp) != -1)
	{
		if (sscanf (line, "%80s", idstr) < 1) // get id
			continue; // whitespace line

		if (strncmp (idstr, ATOM_KEY, MAXLEN) == 0) // cmp id
		{
			prm->natoms++;
		}
	}

	// reset the fp
	rewind (fp);

	// malloc
	prm->atoms = (struct prmatom*) _mol_malloc (sizeof (struct prmatom) * prm->natoms);

	// fill vals
	while (getline (&line, &len, fp) != -1)
	{
		if (sscanf (line, "%80s", idstr) < 1) // get id
			continue; // whitespace line

		if (strncmp (idstr, ATOM_KEY, MAXLEN) == 0) // cmp id
		{
			prm->atoms[atomsi].typemaj = (char*) _mol_malloc (20 * sizeof (char));
			prm->atoms[atomsi].typemin = (char*) _mol_malloc (20 * sizeof (char));
			if (sscanf (line, "%*s %d %19s %19s %d %f %f",
							&prm->atoms[atomsi].id, // this is reset below after qsort
							prm->atoms[atomsi].typemaj,
							prm->atoms[atomsi].typemin,
							&prm->atoms[atomsi].subid,
							&prm->atoms[atomsi].r,
							&prm->atoms[atomsi].q
							) < 6)
			{
				fprintf (stderr, "error in file %s:\n", path);
				fprintf (stderr, "wrong format for line:\n%s", line);
				exit (EXIT_FAILURE);
			}
			// find nsubatoms
			if (prm->atoms[atomsi].subid+1 > prm->nsubatoms)
				prm->nsubatoms = prm->atoms[atomsi].subid+1;

			atomsi++;
		}
	}
	if (line)
		free (line);

	myfclose (fp);

	if (prm->nsubatoms == 0)
	{
		fprintf (stderr, "error in file %s:\n", path);
		fprintf (stderr, "there are no subatom types >= 0");
		exit (EXIT_FAILURE);
	}

	// qsort atoms for bsearch
	qsort (prm->atoms, prm->natoms, sizeof (struct prmatom), comp_prmatom);
	for (int i=0; i<prm->natoms; i++) // assign the id after sort
	{
		prm->atoms[i].id=i;
	}
}


/**
	Read the PDB file, and populate only the typemin and typemaj fields 
	of prm->atoms with only a single entry for each unique <typemaj, typemin> 
	tuple in the PDB. This is mainly for compatibility with other parts
	of the code already written. However, for large PDB's this will also 
	save some space.
	[ Rezaul Chowdhury, Nov. 17, 2010 ]
*/
void read_typeinfo_from_pdb (struct prm* prm, const char* path)
{
	prm->atoms = NULL;
	prm->natoms = 0;
	prm->nsubatoms = -1;

	char* line = NULL;
	size_t len = 0;

	FILE* fp = myfopen (path, "r"); // open pdb file
	
	printf( "\n\n##%s##\n\n", path );
	fflush( stdout );

	// preparse
	while (getline (&line, &len, fp) != -1)
	{
		if (strncmp (line, "ATOM  ", 6) == 0 ) continue; // check for ATOM line

		prm->natoms++;
        } 

	// reset the fp
	rewind (fp);

	// malloc
	prm->atoms = (struct prmatom*) _mol_malloc (sizeof (struct prmatom) * prm->natoms);

	// read every line of the pdb file
	int atomsi = 0;
	while (getline (&line, &len, fp) != -1)
	{
		if (strncmp (line, "ATOM  ", 6) != 0 ) // check for ATOM line
			continue;

      		prm->atoms[atomsi].typemaj = (char*) _mol_malloc (20 * sizeof (char));
      		prm->atoms[atomsi].typemin = (char*) _mol_malloc (20 * sizeof (char));

		if (sscanf (line, "%*s %*d %4s %4s", prm->atoms[atomsi].typemin, prm->atoms[atomsi].typemaj) < 2)
		{
			fprintf (stderr, "error: in file %s line %s: incorrect atom line\n", path, line);
			exit (EXIT_FAILURE);
		}

      		atomsi++;
      	}
      	
	if (line)
		free (line);

	myfclose (fp);

	// qsort atoms for bsearch
	qsort (prm->atoms, prm->natoms, sizeof (struct prmatom), comp_prmatom);

        atomsi = 0;
	prm->atoms[atomsi].id=atomsi;
	
        int atomi = 1;
        while (atomi < prm->natoms)
        {
        	if ( ( strcmp( prm->atoms[atomsi].typemaj, prm->atoms[atomi].typemaj ) == 0 ) 
        		&& ( strcmp( prm->atoms[atomsi].typemin, prm->atoms[atomi].typemin ) == 0 ) )
        	{
        		free( prm->atoms[atomi].typemaj );
        		free( prm->atoms[atomi].typemin );        		
        	}
        	else
        	{
        		atomsi++;
        		prm->atoms[atomsi].typemaj = prm->atoms[atomi].typemaj;
        		prm->atoms[atomsi].typemin = prm->atoms[atomi].typemin;        		
			prm->atoms[atomsi].id=atomsi;        		
        		atomi++;
        	}
        }      		
                
        prm->natoms = atomsi + 1;

	prm->atoms = (struct prmatom*) _mol_realloc ((void *) prm->atoms, sizeof (struct prmatom) * prm->natoms);
}



void read_prmpwpot (struct prm* prm, const char* path)
{
	int k = 0; // number of eigenvalues
	int rflag = 0; // radius flag
	int lambdasi=0;
	int Xsi=0,j;
	char idstr[MAXLEN];
	prm->pwpot = NULL;
	float tmpscanfloat = 0.0;

	char* line = NULL;
	size_t len = 0;
	char* tmpline; // for incrementing the scanned line
	int n; // n chars read in sscanf


	FILE* fp = myfopen (path, "r"); // open prm file

	// preparse
	while (getline (&line, &len, fp) != -1)
	{
		if (sscanf (line, "%80s", idstr) < 1) // get id
			continue; // whitespace line

		if (strncmp (idstr, PWPOT_KEY, MAXLEN) == 0) // cmp id
		{
			k++;
		}
	}

	if (k == 0)
		return; // no pwpot

	// malloc
	prm->pwpot = (struct prmpwpot*) _mol_malloc (sizeof (struct prmpwpot));
	prm->pwpot->k = k;
	prm->pwpot->lambdas = (float*) _mol_malloc (sizeof (float) * prm->pwpot->k);
	prm->pwpot->Xs = (float*) _mol_malloc (sizeof (float) * prm->pwpot->k * prm->nsubatoms);

	// reset the fp
	rewind (fp);

	// fill vals
	while (getline (&line, &len, fp) != -1)
	{
		if (sscanf (line, "%80s", idstr) < 1) // get id
			continue; // whitespace line

		if (strncmp (idstr, PWPOT_KEY, MAXLEN) == 0) // cmp id
		{
			tmpline = line; // for incrementing the scan

			if (sscanf (tmpline, "%*s %f %n", &prm->pwpot->lambdas[lambdasi], &n) < 1)
			{
				fprintf (stderr, "error in file %s:\n", path);
				fprintf (stderr, "wrong format for %s line\n", PWPOT_KEY);
				exit (EXIT_FAILURE);
			}

			lambdasi++;
			tmpline += n;

			for (j = 0; j < prm->nsubatoms; j++)
			{
				if (sscanf (tmpline, "%f %n", &prm->pwpot->Xs[(Xsi*prm->nsubatoms)+j], &n) < 1)
				{
					fprintf (stderr, "error in file %s:\n", path);
					fprintf (stderr, "wrong format for potential line,\n");
					fprintf (stderr, "expecting %d components of eigenvector %d\n", prm->nsubatoms, Xsi);
					exit (EXIT_FAILURE);
				}
				tmpline += n;
			}
			// check to see if there are any more floats
			if (sscanf (tmpline, "%f", &tmpscanfloat) == 1)
			{
				fprintf (stderr, "error in file %s:\n", path);
				fprintf (stderr, "wrong format for potential line,\n");
				fprintf (stderr, "length of eigenvector not equal to number of subatoms\n");
				exit (EXIT_FAILURE);
			}
			Xsi++;
		}
	}

	// reset the fp
	rewind (fp);

	// get r1,r2
	while (getline (&line, &len, fp) != -1)
	{
		if (sscanf (line, "%80s", idstr) < 1) // get id
			continue; // whitespace line

		if (strncmp (idstr, PWPOT_R_KEY, MAXLEN) == 0) // cmp id
		{
			if (sscanf (line, "%*s %f %f", &prm->pwpot->r1, &prm->pwpot->r2) < 2)
				print_readerr (ERR_RADIUS, path, line);
			rflag = 1;
		}
	}
	if (line)
		free (line);

	myfclose (fp);

	if (rflag == 0)
	{
		fprintf (stderr, "error in file %s:\n", path);
		fprintf (stderr, "there are no pwpot.r defintions");
		exit (EXIT_FAILURE);
	}
}

void read_prmbond (struct prm* prm, const char* path)
{
	prm->bonds = NULL;
	prm->nbonds = 0;
	int bondsi = 0;
	char idstr[MAXLEN];

	char* line = NULL;
	size_t len = 0;

	FILE* fp = myfopen (path, "r"); // open prm file

	// preparse
	while (getline (&line, &len, fp) != -1)
	{
		if (sscanf (line, "%80s", idstr) < 1) // get id
			continue; // whitespace line

		if (strncmp (idstr, BOND_KEY, MAXLEN) == 0) // cmp id
		{
			prm->nbonds++;
		}
	}

	// malloc
	prm->bonds = (struct prmbond*) _mol_malloc (sizeof (struct prmbond) * prm->nbonds);

	// reset the fp
	rewind (fp);

	// fill vals
	while (getline (&line, &len, fp) != -1)
	{
		if (sscanf (line, "%80s", idstr) < 1) // get id
			continue; // whitespace line

		if (strncmp (idstr, BOND_KEY, MAXLEN) == 0) // cmp id
		{
			if (sscanf (line, "%*s %d %d %f %f",
							&prm->bonds[bondsi].i,
							&prm->bonds[bondsi].j,
							&prm->bonds[bondsi].k,
							&prm->bonds[bondsi].l0
							) < 4)
			{
				fprintf (stderr, "error in file %s:\n", path);
				fprintf (stderr, "wrong format for line:\n%s", line);
				exit (EXIT_FAILURE);
			}
			// check subatom types for consistency
			if (prm->bonds[bondsi].i+1 > prm->nsubatoms ||
				prm->bonds[bondsi].j+1 > prm->nsubatoms)
			{
				fprintf (stderr, "error in file %s:\n", path);
				fprintf (stderr, "at line:\n%s", line);
				fprintf (stderr, "%s subatom mapping is > max atom subatom mapping\n", BOND_KEY);
				exit (EXIT_FAILURE);
			}

			bondsi++;
		}
	}
	if (line)
		free (line);

	myfclose (fp);

	// qsort bonds
	qsort (prm->bonds, prm->nbonds, sizeof (struct prmbond), comp_prmbond);
}

int atomid (struct prm* prm, const char* typemaj, const char* typemin)
{
	struct prmatom atomkey;
	atomkey.typemaj = (char*) _mol_malloc (20 * sizeof (char));
	atomkey.typemin = (char*) _mol_malloc (20 * sizeof (char));

	atomkey.typemaj = strncpy (atomkey.typemaj, typemaj, 20);
	atomkey.typemin = strncpy (atomkey.typemin, typemin, 20);

	struct prmatom* atomres = (prmatom*) bsearch (&atomkey, prm->atoms, prm->natoms, sizeof (struct prmatom), comp_prmatom);

	free (atomkey.typemaj);
	free (atomkey.typemin);

	if (atomres == NULL)
		return -1;
	else
		return atomres->id;
}

//currently only copies prmatom
struct prm* copy_prm (struct prm* iprm)
{
	int i;
	int natoms = iprm->natoms;

	struct prm* oprm = (struct prm*) _mol_malloc (sizeof (struct prm));

	oprm->natoms = iprm->natoms;
	oprm->atoms = (struct prmatom*) _mol_malloc (natoms * sizeof (struct prmatom));

	oprm->nsubatoms = iprm->nsubatoms;

	for (i = 0; i < natoms; i++)
	{
		copy_prmatom (&iprm->atoms[i], &oprm->atoms[i]);
	}

	/*
	oprm->hydrogen_typen = iprm->hydrogen_typen;

	oprm->pwpot->r1 = iprm->pwpot->r1;
	oprm->pwpot->r2 = iprm->pwpot->r2;

	oprm->pwpot->k = iprm->pwpot->k;
	oprm->lambdas = (float*) _mol_malloc (sizeof (float) * oprm->pwpot->k);
	oprm->Xs = (float*) _mol_malloc (sizeof (float) * oprm->pwpot->k * oprm->pwpot->k);

	for (i = 0; i < oprm->pwpot->k; i++)
	{
		oprm->pwpot->lambdas[i] = iprm->pwpot->lambdas[i];
		int j;
		for (j = 0; j < oprm->pwpot->k; j++)
		{
			int ai = (i*oprm->pwpot->k)+j;
			oprm->pwpot->Xs[ai] = iprm->Xs[ai];
		}
	}
	*/

	return oprm;
}

void copy_prmatom (struct prmatom* iatom, struct prmatom* oatom)
{
	size_t satypemaj = 20 * sizeof (char); //size of atypemaj for malloc and strncpy
	size_t satypemin = 20 * sizeof (char); //size of atypemin for malloc and strncpy

	oatom->id = iatom->id;

	oatom->typemaj = (char*) _mol_malloc (satypemaj);
	strncpy (oatom->typemaj, iatom->typemaj, satypemaj);

	oatom->typemin = (char*) _mol_malloc (satypemin);
	strncpy (oatom->typemin, iatom->typemin, satypemin);

	oatom->subid = iatom->subid;
	oatom->r = iatom->r;
	oatom->q = iatom->q;
}

void modify_prms_radii (struct prm* prm, float k)
{
	int i;
	for (i = 0; i < prm->natoms; i++)
	{
		prm->atoms[i].r *= k;
	}
}

/*
void free_prms (struct prm* prm)
{
	// free all components of prm
	int i;
	for (i = 0; i < prm->natoms; i++)
	{
		//free (prm->atom_type_names[i]);
		free (prm->atoms[i].typemaj);
		free (prm->atoms[i].typemin);
	}
	//free (prm->atom_type_names);

	free (prm->subatoms);
	free (prm->rs);
	free (prm->chrgs);
	free (prm->Xs);
	free (prm->lambdas);

	free (prm); // free prm itself
}
*/

void print_readerr (enum ereaderr readerr, const char* path, char* line)
{
	fprintf (stderr, "error in file %s:\n", path);
	switch (readerr)
	{
		case ERR_VERSION:
			fprintf (stderr, "wrong format for VERSION line\n");
			break;

		case ERR_ATOM:
			fprintf (stderr, "wrong format for atom line\n");
			break;

		case ERR_NAMELEN:
			fprintf (stderr, "length of atom type name is too long (it should be less than 40 chars)\n");
			break;

		case ERR_HYDROGEN:
			fprintf (stderr, "wrong format for hydrogen line\n");
			break;

		case ERR_RADIUS:
			fprintf (stderr, "wrong format for radius line\n");
			break;

		case ERR_POTENTIAL:
			fprintf (stderr, "wrong format for potential line\n");
			break;

			/*
			   case ERR_EIGENVALUE:
			   printf ("not enough eigenvalues\n");
			   break;

			   case ERR_EIGENVECTOR:
			   printf ("not enough eigenvector values\n");
			   break;
			 */
	}
	fprintf (stderr, "line: %s\n", line);

	exit (EXIT_FAILURE);
}

void print_prm (struct prm* prm)
{
	print_prmatom (prm);
	print_prmpwpot (prm);
	print_prmbond (prm);

	/*
	printf ("atypeh:\n");
	for (i = 0; i < prm->natoms; i++)
	{
		printf ("atypeh.key: %s  atypeh.val: %d\n", prm->atypeh[i].key, prm->atypeh[i].val);
	}

				printf ("hydrogen typen:\n");
				printf ("%d\n", prm->hydrogen_typen);
				*/
}

void print_prmatom (struct prm* prm)
{
	int i;

	if (prm->atoms == NULL)
		return;

	printf ("number of atoms: %d\n", prm->natoms);
	printf ("number of subatoms: %d\n", prm->nsubatoms);
	for (i = 0; i < prm->natoms; i++)
	{
		printf ("%d\t", prm->atoms[i].id);
		printf ("%s\t", prm->atoms[i].typemaj);
		printf ("%s\t", prm->atoms[i].typemin);
		printf ("%d\t", prm->atoms[i].subid);
		printf ("%8.3f\t", prm->atoms[i].r);
		printf ("%8.3f\t", prm->atoms[i].q);
		printf ("\n");
	}
}

void print_prmpwpot (struct prm* prm)
{
	int i,j;

	if (prm->pwpot == NULL)
		return;

	printf ("eigenvalues:\n");
	for (i = 0; i < prm->pwpot->k; i++)
	{
		printf ("%2d%10.4f\n", i, prm->pwpot->lambdas[i]);
	}

	printf ("eigenvectors:\n");
	for (i = 0; i < prm->pwpot->k; i++)
	{
		for (j = 0; j < prm->nsubatoms; j++)
		{
			//printf ("%02d,%02d: %08.4f  ", i, j, prm->pwpot->Xs[(i*prm->pwpot->k)+j]);
			printf ("%8.4f", prm->pwpot->Xs[(i*prm->nsubatoms)+j]);
		}
		printf ("\n");
	}
}

void print_prmbond (struct prm* prm)
{
	int i;

	if (prm->bonds == NULL)
		return;

	printf ("number of bonds: %d\n", prm->nbonds);
	for (i = 0; i < prm->nbonds; i++)
	{
		printf ("%d\t", prm->bonds[i].i);
		printf ("%d\t", prm->bonds[i].j);
		printf ("%.3f\t", prm->bonds[i].k);
		printf ("%.3f\t", prm->bonds[i].l0);
		printf ("\n");
	}
}

void
destroy_prmatom (struct prmatom* atom)
{
    free (atom->typemaj);
    free (atom->typemin);
}

void
destroy_prmpwpot (struct prmpwpot* pwpot)
{
    free (pwpot->lambdas);
    free (pwpot->Xs);
}

void
free_prmpwpot (struct prmpwpot* pwpot)
{
    destroy_prmpwpot (pwpot);
    free(pwpot);
}

void
destroy_prm (struct prm* prm)
{
    for (int i=0; i<prm->natoms; i++)
    {
        destroy_prmatom (&(prm->atoms[i]));
    }
    free (prm->atoms);

    free_prmpwpot ( prm->pwpot );

    free (prm->bonds);
}

void
free_prm (struct prm* prm)
{
    destroy_prm (prm);
    free(prm);
}

};
