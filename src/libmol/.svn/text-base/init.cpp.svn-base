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

#define MAXSLEN 200

char* main_dir ()
{
	fprintf (stderr, "%s: Deprecated function\n", __func__);

	size_t slen; // string length

	char here[] = ".";

	char* mdir = (char*) _mol_malloc (MAXSLEN * sizeof (char));

	slen = strlen (here);
	if (slen > MAXSLEN-20)
	{
		fprintf (stderr, "string %s is too long\n", here);
		exit (EXIT_FAILURE);
	}

	mdir = strcpy (mdir, here);
	//mdir = strcat (mdir, "/.mol");
	mdir = strcat (mdir, "");

	return mdir;
}

char* prms_dir ()
{
	fprintf (stderr, "%s: Deprecated function\n", __func__);

	size_t slen; // string length

	char* mdir = main_dir ();
	char* pdir = (char*) _mol_malloc (MAXSLEN * sizeof (char));

	slen = strlen (mdir);
	if (slen > MAXSLEN-20)
	{
		fprintf (stderr, "string %s is too long\n", mdir);
		exit (EXIT_FAILURE);
	}

	pdir = strcpy (pdir, mdir);
	free (mdir);
	pdir = strcat (pdir, "/prm");

	return pdir;
}

char* atom_prm_file (const char* atom_prm)
{
	fprintf (stderr, "%s: Deprecated function\n", __func__);

	size_t slen; // string length

	char* pdir = prms_dir ();
	char* atom_pfile = (char*) _mol_malloc (MAXSLEN * sizeof (char));

	slen = strlen (pdir);
	if (slen > MAXSLEN-20)
	{
		fprintf (stderr, "string %s is too long\n", pdir);
		exit (EXIT_FAILURE);
	}
	slen = strlen (pdir) + strlen (atom_prm);
	if (slen > MAXSLEN-20)
	{
		fprintf (stderr, "string %s is too long\n", atom_prm);
		exit (EXIT_FAILURE);
	}

	atom_pfile = strcpy (atom_pfile, pdir);
	free (pdir);
	atom_pfile = strcat (atom_pfile, "/");
	atom_pfile = strcat (atom_pfile, atom_prm);

	return atom_pfile;
}

char* coeffs_prm_file (const char* coeffs_prm)
{
	fprintf (stderr, "%s: Deprecated function\n", __func__);

	size_t slen; // string length

	char* pdir = prms_dir ();
	char* coeffs_pfile = (char*) _mol_malloc (MAXSLEN * sizeof (char));

	slen = strlen (pdir);
	if (slen > MAXSLEN-20)
	{
		fprintf (stderr, "string %s is too long\n", pdir);
		exit (EXIT_FAILURE);
	}
	slen = strlen (pdir) + strlen (coeffs_prm);
	if (slen > MAXSLEN-20)
	{
		fprintf (stderr, "string %s is too long\n", coeffs_prm);
		exit (EXIT_FAILURE);
	}

	coeffs_pfile = strcpy (coeffs_pfile, pdir);
	free (pdir);
	coeffs_pfile = strcat (coeffs_pfile, "/");
	coeffs_pfile = strcat (coeffs_pfile, coeffs_prm);

	return coeffs_pfile;
}

char* rots_prm_file (const char* rots_prm)
{
	fprintf (stderr, "%s: Deprecated function\n", __func__);

	size_t slen; // string length

	char* pdir = prms_dir ();
	char* rots_pfile = (char*) _mol_malloc (MAXSLEN * sizeof (char));

	slen = strlen (pdir);
	if (slen > MAXSLEN-20)
	{
		fprintf (stderr, "string %s is too long\n", pdir);
		exit (EXIT_FAILURE);
	}
	slen = strlen (pdir) + strlen (rots_prm);
	if (slen > MAXSLEN-20)
	{
		fprintf (stderr, "string %s is too long\n", rots_prm);
		exit (EXIT_FAILURE);
	}

	rots_pfile = strcpy (rots_pfile, pdir);
	free (pdir);
	rots_pfile = strcat (rots_pfile, "/");
	rots_pfile = strcat (rots_pfile, rots_prm);

	return rots_pfile;
}

};
