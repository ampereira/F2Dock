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
#include <stdbool.h>
#include <string.h>
#include <libmol/libmol.h>

namespace LibMol{


void assign_combined_residue_sequence( struct atomgrp *ag )
{
   int k = 1; 
   for ( int i = 0; i < ag->natoms; i++ )
     {
       if ( i && ( ag->atoms[ i ].res_seq != ag->atoms[ i - 1 ].res_seq ) )
         {
           k++;
           if ( ag->atoms[ i ].res_seq < ag->atoms[ i - 1 ].res_seq ) k += 99;
         }
         
       ag->atoms[ i ].comb_res_seq = k;  
     }
}


struct atomgrp* read_pdb (const char* path, struct prm* prm)
{
	FILE* fp = myfopen (path, "r");

	struct atomgrp* ag = (struct atomgrp*) _mol_calloc (1, sizeof (struct atomgrp));
	ag->natoms = 100; // just a guess, realloc as necessary
	ag->atoms = (struct atom*) _mol_malloc (sizeof (struct atom) * ag->natoms);

	char* line = NULL;
	size_t len = 0;

	char atypemaj[5];
	char atypemin[5];

	// read every line of the pdb file
	int atomi = 0;
	while (getline (&line, &len, fp) != -1)
	{
		if (strncmp (line, "ATOM  ", 6) != 0 ) // check for ATOM line
			continue;

		if (atomi+1 > ag->natoms)
		{
			ag->natoms *= 2;
			ag->atoms = (struct atom*) _mol_realloc (ag->atoms, sizeof (struct atom) * ag->natoms);
		}

		/*
		// init bonds
		ag->atoms[atomi].bonds[0] = -1;
		ag->atoms[atomi].bonds[1] = -1;
		ag->atoms[atomi].bonds[2] = -1;
		ag->atoms[atomi].bonds[3] = -1;
		*/

		// init sa
		ag->atoms[atomi].sa = -1;
		// init mask
		ag->atoms[atomi].mask = 0;
		// init attl
		ag->atoms[atomi].attl = 0.0;

		if (sscanf (line, "%*s %*d %4s %4s", atypemin, atypemaj) < 2)
		{
			fprintf (stderr, "error: in file %s line %s: incorrect atom line\n", path, line);
			exit (EXIT_FAILURE);
		}

		if (strncmp (atypemin, "HT", 2) == 0)
			continue; // ignore terminal hydrogens
		if (strncmp (atypemin, "OCT1", 4) == 0)
		{
			atypemin[1] = '\0';
		}
		if (strncmp (atypemin, "OCT2", 4) == 0)
			continue;


		ag->atoms[atomi].atom_typen = atomid (prm, atypemaj, atypemin);
		if (ag->atoms[atomi].atom_typen == -1) // val not found
		{
			if (atypemin[0] == 'H') // try one more time for hydrogen
			{
				atypemin[1] = '\0';
				ag->atoms[atomi].atom_typen = atomid (prm, atypemaj, atypemin);
			}
			if (ag->atoms[atomi].atom_typen == -1) // val still not found
			{
				fprintf (stderr, "error: in file %s line %s: atom type number of %s %s not defined in prm\n", path, line, atypemaj, atypemin);
				exit (EXIT_FAILURE);
			}
		}				
		
		if ( !strcmp( atypemin, "C" ) || !strcmp( atypemin, "CA" ) || !strcmp( atypemin, "N" ) || !strcmp( atypemin, "O" ) || !strcmp( atypemin, "H" ) )
			ag->atoms[atomi].backbone = 1;
		else	
			ag->atoms[atomi].backbone = 0;
		
		ag->atoms[atomi].X = atof (&line[30]);
		ag->atoms[atomi].Y = atof (&line[38]);
		ag->atoms[atomi].Z = atof (&line[46]);
                ag->atoms[atomi].B = atof (&line[60]);
                
                ag->atoms[atomi].icode = line[26];                
                line[26] = 0;
                ag->atoms[atomi].res_seq = atoi(&line[21]);                

                ag->atoms[atomi].hprop = UNKNOWN_HPROP;
                ag->atoms[atomi].base = ag->atoms[atomi].base2 = -1;                

		atomi++;
	}
	if (line)
		free (line);
	myfclose (fp);

	// final realloc of the arrays to make them tight
	ag->natoms = atomi;
	ag->atoms = (struct atom*) _mol_realloc (ag->atoms, sizeof (struct atom) * ag->natoms);
	
        ag->is_psf_read = false;
        ag->prm = prm;

//	check_atomgrp (ag, prm);

	return ag;
}

/**
	Read the PDB file, but assign an integer id to each atom
	based on <residue-name, atom-name> instead of explicitly 
	storing the name of the atom and the residue it is part of.
	This is mainly for compatibility with other parts of the 
	code already written. However, for large PDB's this will also 
	save some space.
	[ Rezaul Chowdhury, Nov. 17, 2010 ]
*/
struct atomgrp* read_pdb_with_compressed_typeinfo (const char* path, struct prm* prm)
{
        read_typeinfo_from_pdb (prm, path);

	FILE* fp = myfopen (path, "r");

	struct atomgrp* ag = (struct atomgrp*) _mol_calloc (1, sizeof (struct atomgrp));
	ag->natoms = 100; // just a guess, realloc as necessary
	ag->atoms = (struct atom*) _mol_malloc (sizeof (struct atom) * ag->natoms);

	char* line = NULL;
	size_t len = 0;

	char atypemaj[5];
	char atypemin[5];

	// read every line of the pdb file
	int atomi = 0;
	while (getline (&line, &len, fp) != -1)
	{
		if (strncmp (line, "ATOM  ", 6) != 0 ) // check for ATOM line
			continue;

		if (atomi+1 > ag->natoms)
		{
			ag->natoms *= 2;
			ag->atoms = (struct atom*) _mol_realloc (ag->atoms, sizeof (struct atom) * ag->natoms);
		}

		/*
		// init bonds
		ag->atoms[atomi].bonds[0] = -1;
		ag->atoms[atomi].bonds[1] = -1;
		ag->atoms[atomi].bonds[2] = -1;
		ag->atoms[atomi].bonds[3] = -1;
		*/

		// init sa
		ag->atoms[atomi].sa = -1;

		if (sscanf (line, "%*s %*d %4s %4s", atypemin, atypemaj) < 2)
		{
			fprintf (stderr, "error: in file %s line %s: incorrect atom line\n", path, line);
			exit (EXIT_FAILURE);
		}

		ag->atoms[atomi].atom_typen = atomid (prm, atypemaj, atypemin);

		if ( !strcmp( atypemin, "C" ) || !strcmp( atypemin, "CA" ) || !strcmp( atypemin, "N" ) || !strcmp( atypemin, "O" ) || !strcmp( atypemin, "H" ) )
			ag->atoms[atomi].backbone = 1;
		else	
			ag->atoms[atomi].backbone = 0;

		ag->atoms[atomi].X = atof (&line[30]);
		ag->atoms[atomi].Y = atof (&line[38]);
		ag->atoms[atomi].Z = atof (&line[46]);
                ag->atoms[atomi].B = atof (&line[60]);
                
                ag->atoms[atomi].icode = line[26];                
                line[26] = 0;
                ag->atoms[atomi].res_seq = atoi(&line[21]);                

                ag->atoms[atomi].hprop = UNKNOWN_HPROP;
                ag->atoms[atomi].base = ag->atoms[atomi].base2 = -1;                
                                
		atomi++;
	}
	if (line)
		free (line);
	myfclose (fp);

	// final realloc of the arrays to make them tight
	ag->natoms = atomi;
	ag->atoms = (struct atom*) _mol_realloc (ag->atoms, sizeof (struct atom) * ag->natoms);

	ag->is_psf_read = false;
        ag->prm = prm;

	return ag;
}


struct atomgrp* read_pdb_nopar (const char* path)
{
        FILE* fp = myfopen (path, "r");

        struct atomgrp* ag = (struct atomgrp*) _mol_calloc (1, sizeof (struct atomgrp));
        ag->natoms = 100; // just a guess, realloc as necessary
        ag->atoms = (struct atom*) _mol_malloc (sizeof (struct atom) * ag->natoms);

        char* line = NULL;
        size_t len = 0;


        // read every line of the pdb file
        int atomi = 0;
        while (getline (&line, &len, fp) != -1)
        {
                if (strncmp (line, "ATOM  ", 6) != 0 && strncmp (line, "HETATM", 6) != 0)
                        continue;
                if (atomi+1 > ag->natoms)
                {
                        ag->natoms *= 2;
                        ag->atoms = (struct atom*) _mol_realloc (ag->atoms, sizeof (struct atom) * ag->natoms);
                }

                // init bonds
/*
                ag->atoms[atomi].bonds[0] = -1;
                ag->atoms[atomi].bonds[1] = -1;
                ag->atoms[atomi].bonds[2] = -1;
                ag->atoms[atomi].bonds[3] = -1;
*/

                // init sa
                ag->atoms[atomi].sa = -1;
                // init mask
//                ag->atoms[atomi].mask = 0;
                // init attl
//                ag->atoms[atomi].attl = 0.0;


                ag->atoms[atomi].atom_typen = 1;

		if ( !strncmp( line + 13, "C ", 2 ) || !strncmp( line + 13, "CA ", 3 ) || !strncmp( line + 13, "N ", 2 ) || !strncmp( line + 13, "O ", 2 ) || !strncmp( line + 13, "H ", 2 ) )
			ag->atoms[atomi].backbone = 1;
		else	
			ag->atoms[atomi].backbone = 0;
                
                ag->atoms[atomi].X = atof (&line[30]);
                ag->atoms[atomi].Y = atof (&line[38]);
                ag->atoms[atomi].Z = atof (&line[46]);

                ag->atoms[atomi].icode = line[26];                
                line[26] = 0;
                ag->atoms[atomi].res_seq = atoi(&line[21]);

                ag->atoms[atomi].hprop = UNKNOWN_HPROP;
                ag->atoms[atomi].base = ag->atoms[atomi].base2 = -1;                
                                
                atomi++;
        }
        if (line)
                free (line);
        myfclose (fp);

        // final realloc of the arrays to make them tight
        ag->natoms = atomi;
        ag->atoms = (struct atom*) _mol_realloc (ag->atoms, sizeof (struct atom) * ag->natoms);

        ag->is_psf_read = false;
        ag->prm = NULL;

        return ag;
}

struct atomgrp** read_pdb_models (const char* path, struct prm* prm, int* rmodels)
{
	FILE* fp = myfopen (path, "r");

    int nmodels=100;//Initial guess
	struct atomgrp** ag_models = (struct atomgrp**) _mol_malloc (nmodels*sizeof (struct atomgrp*));
	char* line = NULL;
	size_t len = 0;
	char atypemaj[5];
	char atypemin[5];
        int  modeli=-1;
        while (getline (&line, &len, fp) != -1)
	{   
            if (strncmp (line, "MODEL", 5) != 0) // check for MODEL line
			continue;
            //We are now in model loop
            if ((modeli+1 > (*rmodels-1)) && ((*rmodels)!=-1))//If were over requested models and we don't read everyhting
		break;
            if (modeli+1 > nmodels)
		{
		        printf("realloc\n");
		        nmodels *= 2;
			ag_models = (struct atomgrp**) _mol_realloc (ag_models, nmodels*sizeof (struct atomgrp*));
		}
            modeli++;
	    ag_models[modeli] = (struct atomgrp*) _mol_calloc (1, sizeof (struct atomgrp));
	    ag_models[modeli]->natoms = 100; // just a guess, realloc as necessary
            ag_models[modeli]->atoms = (struct atom*) _mol_malloc (sizeof (struct atom) * ag_models[modeli]->natoms);
	    int atomi = 0;
	    while (getline (&line, &len, fp) != -1)
	    {       
	    	if (strncmp (line, "ENDMDL", 6) == 0) //if ENDMDL we finished group reading
		{    
		     ag_models[modeli]->natoms = atomi;
		     ag_models[modeli]->atoms = (struct atom*) _mol_realloc (ag_models[modeli]->atoms, sizeof (struct atom) * ag_models[modeli]->natoms);
		     		     
		     check_atomgrp (ag_models[modeli], prm);
		    break;
		}
	        if (strncmp (line, "ATOM  ", 6) != 0) // check for ATOM line
			continue;
                if (atomi+1 > ag_models[modeli]->natoms)
		{
			ag_models[modeli]->natoms *= 2;
			ag_models[modeli]->atoms = (struct atom*) _mol_realloc (ag_models[modeli]->atoms, sizeof (struct atom) * ag_models[modeli]->natoms);
		}

		// init sa
		ag_models[modeli]->atoms[atomi].sa = -1;
		// init mask
		ag_models[modeli]->atoms[atomi].mask = 0;

		if (sscanf (line, "%*s %*d %4s %4s", atypemin, atypemaj) < 2)
		{
			fprintf (stderr, "error: in file %s line %s: incorrect atom line\n", path, line);
			exit (EXIT_FAILURE);
		}

		if (strncmp (atypemin, "HT", 2) == 0)
			continue; // ignore terminal hydrogens
		if (strncmp (atypemin, "OCT1", 4) == 0)
		{
			atypemin[1] = '\0';
		}
		if (strncmp (atypemin, "OCT2", 4) == 0)
			continue;


		ag_models[modeli]->atoms[atomi].atom_typen = atomid (prm, atypemaj, atypemin);
		if (ag_models[modeli]->atoms[atomi].atom_typen == -1) // val not found
		{
			if (atypemin[0] == 'H') // try one more time for hydrogen
			{
				atypemin[1] = '\0';
				ag_models[modeli]->atoms[atomi].atom_typen = atomid (prm, atypemaj, atypemin);
			}
			if (ag_models[modeli]->atoms[atomi].atom_typen == -1) // val still not found
			{
				fprintf (stderr, "error: in file %s line %s: atom type number of %s %s not defined in prm\n", path, line, atypemaj, atypemin);
				exit (EXIT_FAILURE);
			}
		}


		if ( !strcmp( atypemin, "C" ) || !strcmp( atypemin, "CA" ) || !strcmp( atypemin, "N" ) || !strcmp( atypemin, "O" ) || !strcmp( atypemin, "H" ) )
			ag_models[modeli]->atoms[atomi].backbone = 1;
		else	
			ag_models[modeli]->atoms[atomi].backbone = 0;

		ag_models[modeli]->atoms[atomi].X = atof (&line[30]);
		ag_models[modeli]->atoms[atomi].Y = atof (&line[38]);
		ag_models[modeli]->atoms[atomi].Z = atof (&line[46]);
		
                ag_models[modeli]->atoms[atomi].icode = line[26];                
                line[26] = 0;
                ag_models[modeli]->atoms[atomi].res_seq = atoi(&line[21]);		
		
                ag_models[modeli]->atoms[atomi].hprop = UNKNOWN_HPROP;		
                ag_models[modeli]->atoms[atomi].base = ag_models[modeli]->atoms[atomi].base2 = -1;                
		
                ag_models[modeli]->is_psf_read = false;
                ag_models[modeli]->prm = prm;

		atomi++;
	}
	}
	if (line)
		free (line);
	myfclose (fp);
	nmodels = modeli+1;
	ag_models = (struct atomgrp**) _mol_realloc (ag_models, nmodels*sizeof (struct atomgrp));
	// final realloc of the arrays to make them tight

        *rmodels=modeli+1;
	return ag_models;
}

struct atomgrp **read_pdb_modelsnopar(const char *path, int *rmodels)
{
	FILE *fp = myfopen(path, "r");
	int nmodels = 100;	//Initial guess
	struct atomgrp **ag_models = (struct atomgrp **)_mol_malloc(nmodels * sizeof(struct atomgrp *));
	char *line = NULL;
	size_t len = 0;
	int modeli = -1;
	while (getline(&line, &len, fp) != -1) {
		if (strncmp(line, "MODEL", 5) != 0)	// check for MODEL line
			continue;
		//We are now in model loop
		if ((modeli + 1 > (*rmodels - 1)) && ((*rmodels) != -1))	//If were over requested models and we don't read everyhting
			break;
		if (modeli + 1 > nmodels) {
			printf("realloc\n");
			nmodels *= 2;
			ag_models = (struct atomgrp **)
			    _mol_realloc(ag_models,
					 nmodels * sizeof(struct atomgrp *));
		}
		modeli++;
		ag_models[modeli] = (struct atomgrp *)_mol_calloc(1, sizeof(struct atomgrp));
		ag_models[modeli]->natoms = 16;	// just a guess, realloc as necessary
		ag_models[modeli]->atoms = (struct atom*)
		    _mol_malloc(sizeof(struct atom) *
				ag_models[modeli]->natoms);
                ag_models[modeli]->prm = NULL;				
		int atomi = 0;
		while (getline(&line, &len, fp) != -1) {
			if (strncmp(line, "ENDMDL", 6) == 0)	//if ENDMDL we finished group reading
			{
				ag_models[modeli]->natoms = atomi;
				ag_models[modeli]->atoms = (struct atom*)
				    _mol_realloc(ag_models[modeli]->atoms,
						 sizeof(struct atom) *
						 ag_models[modeli]->natoms);
						 
				break;
			}
			if (strncmp(line, "ATOM  ", 6) != 0)	// check for ATOM line
				continue;
			if (atomi + 1 > ag_models[modeli]->natoms) {
				ag_models[modeli]->natoms *= 2;
				ag_models[modeli]->atoms = (struct atom*)
				    _mol_realloc(ag_models[modeli]->atoms,
						 sizeof(struct atom) *
						 ag_models[modeli]->natoms);
			}

                	if ( !strncmp( line + 13, "C ", 2 ) || !strncmp( line + 13, "CA ", 3 ) || !strncmp( line + 13, "N ", 2 ) || !strncmp( line + 13, "O ", 2 ) || !strncmp( line + 13, "H ", 2 ) )
                		ag_models[modeli]->atoms[atomi].backbone = 1;
                	else	
                		ag_models[modeli]->atoms[atomi].backbone = 0;

			ag_models[modeli]->atoms[atomi].X = atof(&line[30]);
			ag_models[modeli]->atoms[atomi].Y = atof(&line[38]);
			ag_models[modeli]->atoms[atomi].Z = atof(&line[46]);

                        ag_models[modeli]->atoms[atomi].icode = line[26];                
                        line[26] = 0;
                        ag_models[modeli]->atoms[atomi].res_seq = atoi(&line[21]);		

                        ag_models[modeli]->atoms[atomi].hprop = UNKNOWN_HPROP;		                        
                        ag_models[modeli]->atoms[atomi].base = ag_models[modeli]->atoms[atomi].base2 = -1;                
                	
			atomi++;
		}
	}
	if (line)
		free(line);
	myfclose(fp);
	nmodels = modeli + 1;
	ag_models = (struct atomgrp**)_mol_realloc(ag_models, nmodels * sizeof(struct atomgrp));
	// final realloc of the arrays to make them tight

	*rmodels = modeli + 1;
	return ag_models;
}

void fprint_pdb (struct atomgrp* ag, struct prm* prm, const char* path)
{
	FILE* fp = myfopen (path, "w");

	//struct atomgrp** ags = (struct atomgrp**) _mol_malloc (2 * sizeof (struct atomgrp*));
	struct atomgrp** ags = extract_nitrogen_residues (ag, prm);
	//ags[0] = copy_atomgrp (ag);
	//ags[1] = NULL; // flag the end with NULL

	int sum_atomi = 0;
	int agi = 0;
	while (ags[agi] != NULL)
	{
		int atomi;
		for (atomi = 0; atomi < ags[agi]->natoms; atomi++, sum_atomi++)
		{
			fprintf (fp, "%-6s", "ATOM"); // atom number
			fprintf (fp, "%5d", sum_atomi+1); // atom number
			//fprintf (fp, "  "); // 2 spaces
			fprintf (fp, " "); // 1 space

            char atomname[5];
            if (strlen(prm->atoms[ags[agi]->atoms[atomi].atom_typen].typemin) == 4) {
                strcpy(atomname, prm->atoms[ags[agi]->atoms[atomi].atom_typen].typemin);
            } else {
                atomname[0] = ' ';
                strcpy(atomname+1, prm->atoms[ags[agi]->atoms[atomi].atom_typen].typemin);
            }
			fprintf (fp, "%-4.4s", atomname); // atom typemin
			fprintf (fp, " "); // Alternate location indicator
			fprintf (fp, "%-3s", prm->atoms[ags[agi]->atoms[atomi].atom_typen].typemaj); // atom typemaj

			fprintf (fp, " "); // 1 space
			fprintf (fp, "%1s", "A"); // chain id
			fprintf (fp, "%4d", agi+1); // residue number

			fprintf (fp, " "); // code for insertion of residues
			fprintf (fp, "   "); // 3 spaces

			fprintf (fp, "%8.3f", ags[agi]->atoms[atomi].X); // X coordinate
			fprintf (fp, "%8.3f", ags[agi]->atoms[atomi].Y); // Y coordinate
			fprintf (fp, "%8.3f", ags[agi]->atoms[atomi].Z); // Z coordinate
                        fprintf (fp, "  1.00  1.00      AAAA"); // segid
			fprintf (fp, "\n"); // new line
		}
		agi++;
	}
        fprintf (fp,"END\n");//END

	//free (ags);
	free_atomgrps (ags);

	myfclose (fp);
}

void write_pdb_nopar (struct atomgrp* ag, const char* inf, const char* ouf)
{
        FILE* fp = myfopen (inf, "r");
        FILE* fop = myfopen (ouf, "w");

        char* line = NULL;
        size_t len = 0;
        char c;


        // read every line of the pdb file
        int atomi = 0;
        while (getline (&line, &len, fp) != -1)
        {
                if (strncmp (line, "ATOM  ", 6) != 0 && strncmp (line, "HETATM", 6) != 0)
                {
                         fprintf (fop,"%s",line);
                        continue;
                }
                c=line[54];
                sprintf(line+30,"%8.3f",ag->atoms[atomi].X );
                sprintf(line+38,"%8.3f",ag->atoms[atomi].Y );
                sprintf(line+46,"%8.3f",ag->atoms[atomi].Z );
                atomi++;
                line[54]=c;
                fprintf (fop,"%s",line);
        }
        if (line)
                free (line);
        myfclose (fp);
        myfclose (fop);
}

};
