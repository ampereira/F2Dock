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

#include <libmol/libmol.h>

namespace LibMol{


struct matrix2di* matrix2di_create (int ni, int nj)
{
	// error checking
	if (ni < 1 || nj < 1)
	{
		fprintf (stderr, "error: in function matrix2di_create, argument matrix size(s) is less than 1\n");
		exit (EXIT_FAILURE);
	}

	// allocate memory for the matrix
	struct matrix2di* A = (struct matrix2di*)_mol_malloc (sizeof (struct matrix2di));

	// connect the matrix sizes
	A->ni = ni;
	A->nj = nj;

	// allocate memory for the pointers to the rows
	A->vals = (int**)_mol_malloc (ni * sizeof (int*));

	// allocate contiguous block of memory for the matrix vals
	A->vals[0] = (int*)_mol_malloc (ni * nj * sizeof (int));

	// point to the remaining rows
	int i;
	for (i = 0; i < A->ni; i++)
	{
		A->vals[i] = A->vals[0] + (i * A->nj);
	}

	// return the newly created array
	return A;
}
struct matrix2df* matrix2df_create (int ni, int nj)
{
	// error checking
	if (ni < 1 || nj < 1)
	{
		fprintf (stderr, "error: in function matrix2df_create, argument matrix size(s) is less than 1\n");
		exit (EXIT_FAILURE);
	}

	// allocate memory for the matrix
	struct matrix2df* A = (struct matrix2df*)_mol_malloc (sizeof (struct matrix2df));

	// connect the matrix sizes
	A->ni = ni;
	A->nj = nj;

	// allocate memory for the pointers to the rows
	A->vals = (float**) _mol_malloc (ni * sizeof (float*));

	// allocate contiguous block of memory for the matrix vals
	A->vals[0] = (float*)_mol_malloc (ni * nj * sizeof (float));

	// point to the remaining rows
	int i;
	for (i = 0; i < A->ni; i++)
	{
		A->vals[i] = A->vals[0] + (i * A->nj);
	}

	// return the newly created array
	return A;
}


void matrix2di_destroy (struct matrix2di* A)
{
	free (A->vals[0]); // free the contiguous block of vals
	free (A->vals); // free the pointers to the rows
	free (A); // free the matrix itself
}
void matrix2df_destroy (struct matrix2df* A)
{
	free (A->vals[0]); // free the contiguous block of vals
	free (A->vals); // free the pointers to the rows
	free (A); // free the matrix itself
}


void matrix2di_init (struct matrix2di* A, int initval)
{
	int i, j;
	for (i = 0; i < A->ni; i++)
	{
		for (j = 0; j < A->nj; j++)
		{
			A->vals[i][j] = initval;
		}
	}
}
void matrix2df_init (struct matrix2df* A, float initval)
{
	int i, j;
	for (i = 0; i < A->ni; i++)
	{
		for (j = 0; j < A->nj; j++)
		{
			A->vals[i][j] = initval;
		}
	}
}


void matrix2di_print (struct matrix2di* A)
{
	int i, j;
	for (i = 0; i < A->ni; i++)
	{
		for (j = 0; j < A->nj-1; j++)
		{
			printf ("%d\t", A->vals[i][j]);
		}

		// print the last value of the row followed by a newline
		printf ("%d\n", A->vals[i][j]);
	}
}
void matrix2df_print (struct matrix2df* A)
{
	int i, j;
	for (i = 0; i < A->ni; i++)
	{
		for (j = 0; j < A->nj-1; j++)
		{
			printf ("%f\t", A->vals[i][j]);
		}

		// print the last value of the row followed by a newline
		printf ("%f\n", A->vals[i][j]);
	}
}


void matrix2di_fprint (struct matrix2di* A, const char* path)
{
	FILE* fp = myfopen (path, "w");

	int i, j;
	for (i = 0; i < A->ni; i++)
	{
		for (j = 0; j < A->nj-1; j++)
		{
			fprintf (fp, "%d\t", A->vals[i][j]);
		}

		// print the last value of the row followed by a newline
		fprintf (fp, "%d\n", A->vals[i][j]);
	}

	myfclose (fp);
}
void matrix2df_fprint (struct matrix2df* A, const char* path)
{
	FILE* fp = myfopen (path, "w");

	int i, j;
	for (i = 0; i < A->ni; i++)
	{
		for (j = 0; j < A->nj-1; j++)
		{
			fprintf (fp, "%f\t", A->vals[i][j]);
		}

		// print the last value of the row followed by a newline
		fprintf (fp, "%f\n", A->vals[i][j]);
	}

	myfclose (fp);
}


struct matrix2di* matrix2di_read (const char* path)
{
	char* line = NULL;
	size_t len = 0;
	int val; // tmp var for scanf
	char* tmpline; // for incrementing the scanned line

	int ni=0, nj=0; // matrix sizes
	int i, j; // loop index
	int n; // n chars read in sscanf

	FILE* fp = myfopen (path, "r");

	// count the number of newlines (number of rows)
	while (getline (&line, &len, fp) != -1)
	{
		ni++;
	}

	// count the number of words on the first line (number of cols)
	tmpline = line; // for incrementing the scan
	while (sscanf (tmpline, "%d %n", &val, &n) != EOF)
	{
		nj++;
		tmpline += n;
	}

	struct matrix2di* A = matrix2di_create (ni, nj); // create the new matrix of size ni x nj

	rewind (fp); // reset the fp
	// now read the matrix
	i = 0;
	while (getline (&line, &len, fp) != -1)
	{
		tmpline = line; // for incrementing the scan

		for (j = 0; j < nj; j++)
		{
			if (sscanf (tmpline, "%d %n", &val, &n) < 1)
			{
				fprintf (stderr, "begin error\n");
				fprintf (stderr, "in function matrix2di_read\n");
				fprintf (stderr, "inconsistent number of columns in matrix\n");
				fprintf (stderr, "at file:\n");
				fprintf (stderr, "%s\n", path);
				fprintf (stderr, "at line:\n");
				fprintf (stderr, "%s\n", line);
				fprintf (stderr, "end error\n");
				exit (EXIT_FAILURE);
			}
			A->vals[i][j] = val;
			tmpline += n;
		}
		i++;
	}
	if (line)
		free (line);
	myfclose (fp);

	return A;
}
struct matrix2df* matrix2df_read (const char* path)
{
	char* line = NULL;
	size_t len = 0;
	float val; // tmp var for scanf
	char* tmpline; // for incrementing the scanned line

	int ni=0, nj=0; // matrix sizes
	int i, j; // loop index
	int n; // n chars read in sscanf

	FILE* fp = myfopen (path, "r");

	// count the number of newlines (number of rows)
	while (getline (&line, &len, fp) != -1)
	{
		ni++;
	}

	// count the number of words on the first line (number of cols)
	tmpline = line; // for incrementing the scan
	while (sscanf (tmpline, "%f %n", &val, &n) != EOF)
	{
		nj++;
		tmpline += n;
	}

	struct matrix2df* A = matrix2df_create (ni, nj); // create the new matrix of size ni x nj

	rewind (fp); // reset the fp
	// now read the matrix
	i = 0;
	while (getline (&line, &len, fp) != -1)
	{
		tmpline = line; // for incrementing the scan

		for (j = 0; j < nj; j++)
		{
			if (sscanf (tmpline, "%f %n", &val, &n) < 1)
			{
				fprintf (stderr, "begin error\n");
				fprintf (stderr, "in function matrix2df_read\n");
				fprintf (stderr, "inconsistent number of columns in matrix\n");
				fprintf (stderr, "at file:\n");
				fprintf (stderr, "%s\n", path);
				fprintf (stderr, "at line:\n");
				fprintf (stderr, "%s\n", line);
				fprintf (stderr, "end error\n");
				exit (EXIT_FAILURE);
			}
			A->vals[i][j] = val;
			tmpline += n;
		}
		i++;
	}
	if (line)
		free (line);
	myfclose (fp);

	return A;
}


int matrix2di_size_diff (struct matrix2di* A, struct matrix2di* B)
{
	if (A->ni == B->ni && A->nj == B->nj)
		return 0; // matrix sizes are equal
	return 1; // matrix sizes differ
}
int matrix2df_size_diff (struct matrix2df* A, struct matrix2df* B)
{
	if (A->ni == B->ni && A->nj == B->nj)
		return 0; // matrix sizes are equal
	return 1; // matrix sizes differ
}


void matrix2di_add (struct matrix2di* A, struct matrix2di* B, struct matrix2di* C)
{
	matrix2di_pairwise_arith (A, B, C, 0);
}
void matrix2di_sub (struct matrix2di* A, struct matrix2di* B, struct matrix2di* C)
{
	matrix2di_pairwise_arith (A, B, C, 1);
}
void matrix2di_pairwise_mult (struct matrix2di* A, struct matrix2di* B, struct matrix2di* C)
{
	matrix2di_pairwise_arith (A, B, C, 2);
}
void matrix2di_pairwise_div (struct matrix2di* A, struct matrix2di* B, struct matrix2di* C)
{
	matrix2di_pairwise_arith (A, B, C, 3);
}
void matrix2di_pairwise_arith (struct matrix2di* A, struct matrix2di* B, struct matrix2di* C, int op)
{
	if (matrix2di_size_diff (A, B) ||
		matrix2di_size_diff (A, C))
	{
		fprintf (stderr, "begin error\n");
		fprintf (stderr, "at least one pair of argument matrices differ in size\n");
		fprintf (stderr, "end error\n");
		exit (EXIT_FAILURE);
	}

	int i, j;
	for (i = 0; i < A->ni; i++)
	{
		for (j = 0; j < A->nj; j++)
		{
			if (op == 0) // add
			{
				C->vals[i][j] = A->vals[i][j] + B->vals[i][j];
			}
			else if (op == 1) // subtract
			{
				C->vals[i][j] = A->vals[i][j] - B->vals[i][j];
			}
			else if (op == 2) // multiply
			{
				C->vals[i][j] = A->vals[i][j] * B->vals[i][j];
			}
			else if (op == 3) // divide
			{
				C->vals[i][j] = A->vals[i][j] / B->vals[i][j];
			}
			else
			{
				fprintf (stderr, "begin error\n");
				fprintf (stderr, "requested matrix operation unknown\n");
				fprintf (stderr, "end error\n");
				exit (EXIT_FAILURE);
			}
		}
	}
}


void matrix2df_add (struct matrix2df* A, struct matrix2df* B, struct matrix2df* C)
{
	matrix2df_pairwise_arith (A, B, C, 0);
}
void matrix2df_sub (struct matrix2df* A, struct matrix2df* B, struct matrix2df* C)
{
	matrix2df_pairwise_arith (A, B, C, 1);
}
void matrix2df_pairwise_mult (struct matrix2df* A, struct matrix2df* B, struct matrix2df* C)
{
	matrix2df_pairwise_arith (A, B, C, 2);
}
void matrix2df_pairwise_div (struct matrix2df* A, struct matrix2df* B, struct matrix2df* C)
{
	matrix2df_pairwise_arith (A, B, C, 3);
}
void matrix2df_pairwise_arith (struct matrix2df* A, struct matrix2df* B, struct matrix2df* C, int op)
{
	if (matrix2df_size_diff (A, B) ||
		matrix2df_size_diff (A, C))
	{
		fprintf (stderr, "begin error\n");
		fprintf (stderr, "at least one pair of argument matrices differ in size\n");
		fprintf (stderr, "end error\n");
		exit (EXIT_FAILURE);
	}

	int i, j;
	for (i = 0; i < A->ni; i++)
	{
		for (j = 0; j < A->nj; j++)
		{
			if (op == 0) // add
			{
				C->vals[i][j] = A->vals[i][j] + B->vals[i][j];
			}
			else if (op == 1) // subtract
			{
				C->vals[i][j] = A->vals[i][j] - B->vals[i][j];
			}
			else if (op == 2) // multiply
			{
				C->vals[i][j] = A->vals[i][j] * B->vals[i][j];
			}
			else if (op == 3) // divide
			{
				C->vals[i][j] = A->vals[i][j] / B->vals[i][j];
			}
			else
			{
				fprintf (stderr, "begin error\n");
				fprintf (stderr, "requested matrix operation unknown\n");
				fprintf (stderr, "end error\n");
				exit (EXIT_FAILURE);
			}
		}
	}
}


void matrix2df_log (struct matrix2df* A, struct matrix2df* B)
{
	if (matrix2df_size_diff (A, B))
	{
		fprintf (stderr, "begin error\n");
		fprintf (stderr, "in function matrix log\n");
		fprintf (stderr, "at least one pair of argument matrices differ in size\n");
		fprintf (stderr, "end error\n");
		exit (EXIT_FAILURE);
	}

	int i, j;
	for (i = 0; i < A->ni; i++)
	{
		for (j = 0; j < A->nj; j++)
		{
			B->vals[i][j] = logf (A->vals[i][j]);
		}
	}
}

float matrix2df_element_sum (struct matrix2df* A)
{ float sum = 0.0;
	int i, j;
	for (i = 0; i < A->ni; i++)
	{
		for (j = 0; j < A->nj; j++)
		{
			sum += A->vals[i][j];
		}
	}
	return sum;
}

void matrix2df_scalar_mult (struct matrix2df* A, float val, struct matrix2df* B)
{
	if (matrix2df_size_diff (A, B))
	{
		fprintf (stderr, "begin error\n");
		fprintf (stderr, "in function matrix scalar multiplication\n");
		fprintf (stderr, "at least one pair of argument matrices differ in size\n");
		fprintf (stderr, "end error\n");
		exit (EXIT_FAILURE);
	}

	int i, j;
	for (i = 0; i < A->ni; i++)
	{
		for (j = 0; j < A->nj; j++)
		{
			B->vals[i][j] = A->vals[i][j] * val;
		}
	}
}

};
