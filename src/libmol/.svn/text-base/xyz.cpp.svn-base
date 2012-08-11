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

struct atomgrp* read_xyz (const char* path, struct prm* prm)
{
	FILE* fp = myfopen (path, "r");

	struct atomgrp* xyz = (struct atomgrp*) _mol_calloc (1, sizeof (struct atomgrp));

	char* line = NULL;
	size_t len = 0;
	ssize_t read;
	read = getline (&line, &len, fp); // get the header (the number of atoms in xyz file)
	if (read == -1)
	{
		fprintf (stderr, "read error on %s:\n", path);
		exit (EXIT_FAILURE);
	}
	if (sscanf (line, "%d", &xyz->natoms) < 1)
	{
		fprintf (stderr, "format read error on %s:\n", path);
		exit (EXIT_FAILURE);
	}

	xyz->atoms = (struct atom*) _mol_malloc (sizeof (struct atom) * xyz->natoms);

	// read every remaining line of the xyz file
	int atomn = 0;
	while (getline (&line, &len, fp) != -1)
	{
		if (line[0] == '\n' || line[0] == '\r') // this is a blank line
			continue;

		// init bonds
		/*
		xyz->atoms[atomn].bonds[0] = -1;
		xyz->atoms[atomn].bonds[1] = -1;
		xyz->atoms[atomn].bonds[2] = -1;
		xyz->atoms[atomn].bonds[3] = -1;
		*/

		// init sa
		xyz->atoms[atomn].sa = -1;

		//if (sscanf (line, "%*d %*s %f %f %f %d %d %d %d %d",
		if (sscanf (line, "%*d %*s %lf %lf %lf %d",
					&xyz->atoms[atomn].X,
					&xyz->atoms[atomn].Y,
					&xyz->atoms[atomn].Z,
					&xyz->atoms[atomn].atom_typen
					/*
					&xyz->atoms[atomn].bonds[0],
					&xyz->atoms[atomn].bonds[1],
					&xyz->atoms[atomn].bonds[2],
					&xyz->atoms[atomn].bonds[3]
					*/
		   			) < 4)
		{
			fprintf (stderr, "format read error on %s:\n", path);
			exit (EXIT_FAILURE);
		}

		atomn++;
	}
	if (line)
		free (line);
	myfclose (fp);

	check_atomgrp (xyz, prm);

	return xyz;
}

struct atomgrp* read_oldxyz (const char* path, struct prm* prm)
{
	FILE* fp = myfopen (path, "r");

	struct atomgrp* xyz = (struct atomgrp*) _mol_calloc (1, sizeof (struct atomgrp));

	char* line = NULL;
	size_t len = 0;
	ssize_t read;
	read = getline (&line, &len, fp); // get the header (the number of atoms in xyz file)
	if (read == -1)
	{
		fprintf (stderr, "read error on %s:\n", path);
		exit (EXIT_FAILURE);
	}
	if (sscanf (line, "%d", &xyz->natoms) < 1)
	{
		fprintf (stderr, "format read error on %s:\n", path);
		exit (EXIT_FAILURE);
	}

	xyz->atoms = (struct atom*) _mol_malloc (sizeof (struct atom) * xyz->natoms);

	// read every remaining line of the xyz file
	int atomn = 0;
	while (getline (&line, &len, fp) != -1)
	{
		if (line[0] == '\n' || line[0] == '\r') // this is a blank line
			continue;

		// init bonds
		/*
		xyz->atoms[atomn].bonds[0] = -1;
		xyz->atoms[atomn].bonds[1] = -1;
		xyz->atoms[atomn].bonds[2] = -1;
		xyz->atoms[atomn].bonds[3] = -1;
		*/

		// init sa
		xyz->atoms[atomn].sa = -1;
		float sa; // old xyz file has floating point sa

		//if (sscanf (line, "%*d %*s %f %f %f %f %d %d %d %d %d",
		if (sscanf (line, "%*d %*s %lf %lf %lf %f %d",
					&xyz->atoms[atomn].X,
					&xyz->atoms[atomn].Y,
					&xyz->atoms[atomn].Z,
					&sa,
					&xyz->atoms[atomn].atom_typen
					/*
					&xyz->atoms[atomn].bonds[0],
					&xyz->atoms[atomn].bonds[1],
					&xyz->atoms[atomn].bonds[2],
					&xyz->atoms[atomn].bonds[3]
					*/
		   			) < 5)
		{
			fprintf (stderr, "format read error on %s:\n", path);
			exit (EXIT_FAILURE);
		}
		xyz->atoms[atomn].sa = (int) sa;

		atomn++;
	}
	if (line)
		free (line);
	myfclose (fp);

	check_atomgrp (xyz, prm);

	return xyz;
}

void print_xyz (struct atomgrp* xyz, struct prm* prm)
{
	printf ("%6d\n", xyz->natoms); // number of atoms
	int atomn;
	for (atomn = 0; atomn < xyz->natoms; atomn++)
	{
		printf ("%6d", atomn); // atom number
		printf ("  "); // 2 spaces

		printf ("%-3.3s", prm->atoms[xyz->atoms[atomn].atom_typen].typemin); // atom type (left aligned), only print first three chars

		printf ("%12.6f", xyz->atoms[atomn].X); // X coordinate
		printf ("%12.6f", xyz->atoms[atomn].Y); // Y coordinate
		printf ("%12.6f", xyz->atoms[atomn].Z); // Z coordinate
		printf ("   "); // 3 spaces
		printf ("%3d", xyz->atoms[atomn].atom_typen); // atom type number

		/*
		int bi = 0;
		while (bi < 4 && xyz->atoms[atomn].bonds[bi] != -1)
		{
			printf ("%6d", xyz->atoms[atomn].bonds[bi]+1); // atom type number
			bi++;
		}
		*/

		printf ("\n"); // new line
	}
}

void fprint_xyz (struct atomgrp* xyz, struct prm* prm, const char* path)
{
	// open xyz ofile
	FILE* fp = myfopen (path, "w");

	fprintf (fp, "%6d\n", xyz->natoms); // number of atoms
	int atomn;
	for (atomn = 0; atomn < xyz->natoms; atomn++)
	{
		fprintf (fp, "%6d", atomn); // atom number
		fprintf (fp, "  "); // 2 spaces

		fprintf (fp, "%-3.3s", prm->atoms[xyz->atoms[atomn].atom_typen].typemin); // atom type (left aligned), only print first three chars

		fprintf (fp, "%12.6f", xyz->atoms[atomn].X); // X coordinate
		fprintf (fp, "%12.6f", xyz->atoms[atomn].Y); // Y coordinate
		fprintf (fp, "%12.6f", xyz->atoms[atomn].Z); // Z coordinate
		fprintf (fp, "   "); // 3 spaces
		fprintf (fp, "%3d", xyz->atoms[atomn].atom_typen); // atom type number

		/*
		int bi = 0;
		while (bi < 4 && xyz->atoms[atomn].bonds[bi] != -1)
		{
			fprintf (fp, "%6d", xyz->atoms[atomn].bonds[bi]+1); // bonding
			bi++;
		}
		*/

		fprintf (fp, "\n"); // new line
	}
	myfclose (fp);
}

void fprint_oldxyz (struct atomgrp* xyz, struct prm* prm, const char* path)
{
	// open xyz ofile
	FILE* fp = myfopen (path, "w");

	fprintf (fp, "%6d\n", xyz->natoms); // number of atoms
	int atomn;
	for (atomn = 0; atomn < xyz->natoms; atomn++)
	{
		fprintf (fp, "%6d", atomn); // atom number
		fprintf (fp, "  "); // 2 spaces

		fprintf (fp, "%-3.3s", prm->atoms[xyz->atoms[atomn].atom_typen].typemin); // atom type (left aligned), only print first three chars

		fprintf (fp, "%12.6f", xyz->atoms[atomn].X); // X coordinate
		fprintf (fp, "%12.6f", xyz->atoms[atomn].Y); // Y coordinate
		fprintf (fp, "%12.6f", xyz->atoms[atomn].Z); // Z coordinate
		fprintf (fp, "   "); // 3 spaces
		fprintf (fp, "%4.2f", (float) xyz->atoms[atomn].sa);
		fprintf (fp, "%6d", xyz->atoms[atomn].atom_typen); // atom type number

		/*
		int bi = 0;
		while (bi < 4 && xyz->atoms[atomn].bonds[bi] != -1)
		{
			fprintf (fp, "%6d", xyz->atoms[atomn].bonds[bi]+1); // atom type number
			bi++;
		}
		*/

		fprintf (fp, "\n"); // new line
	}
	myfclose (fp);
}

};
