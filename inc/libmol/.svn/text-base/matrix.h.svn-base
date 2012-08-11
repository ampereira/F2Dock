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
#ifndef _MOL_MATRIX_H_
#define _MOL_MATRIX_H_

namespace LibMol{
/** \file matrix.h
	structures and functions for matrix operations
*/

struct matrix2di // 2d matrix of integers
{
	int** vals; // pointer to array of rows
	int ni, nj; // number of rows, number of cols
};

struct matrix2df // 2d matrix of floats
{
	float** vals; // pointer to array of rows
	int ni, nj; // number of rows, number of cols
};

/**
  creates a matrix2di of size ni x nj,
  with vals (vals[0]) allocated as a contiguous block of memory
*/
struct matrix2di* matrix2di_create (int ni, int nj);
struct matrix2df* matrix2df_create (int ni, int nj);

/**
  frees the matrix A and all of its allocated memory
*/
void matrix2di_destroy (struct matrix2di* A);
void matrix2df_destroy (struct matrix2df* A);

/**
  initializes the vals of A to initval
*/
void matrix2di_init (struct matrix2di* A, int initval);
void matrix2df_init (struct matrix2df* A, float initval);

/**
  print the matrix A to stdout
*/
void matrix2di_print (struct matrix2di* A);
void matrix2df_print (struct matrix2df* A);

/**
  fprint the matrix A to the file whose name is the string pointed to by path
*/
void matrix2di_fprint (struct matrix2di* A, const char* path);
void matrix2df_fprint (struct matrix2df* A, const char* path);

/**
  reads the matrix in the file at path
*/
struct matrix2di* matrix2di_read (const char* path);
struct matrix2df* matrix2df_read (const char* path);

/**
  compares ni of A and B, and nj of A and B
  \return 0 if the sizes are equal, 1 if the sizes differ
*/
int matrix2di_size_diff (struct matrix2di* A, struct matrix2di* B);
int matrix2df_size_diff (struct matrix2df* A, struct matrix2df* B);

/**
  perform operation on matrices A and B and store the result in C
*/
void matrix2di_add (struct matrix2di* A, struct matrix2di* B, struct matrix2di* C);
void matrix2di_sub (struct matrix2di* A, struct matrix2di* B, struct matrix2di* C);
void matrix2di_pairwise_mult (struct matrix2di* A, struct matrix2di* B, struct matrix2di* C);
void matrix2di_pairwise_div (struct matrix2di* A, struct matrix2di* B, struct matrix2di* C);
void matrix2di_pairwise_arith (struct matrix2di* A, struct matrix2di* B, struct matrix2di* C, int op);

void matrix2df_add (struct matrix2df* A, struct matrix2df* B, struct matrix2df* C);
void matrix2df_sub (struct matrix2df* A, struct matrix2df* B, struct matrix2df* C);
void matrix2df_pairwise_mult (struct matrix2df* A, struct matrix2df* B, struct matrix2df* C);
void matrix2df_pairwise_div (struct matrix2df* A, struct matrix2df* B, struct matrix2df* C);
void matrix2df_pairwise_arith (struct matrix2df* A, struct matrix2df* B, struct matrix2df* C, int op);

void matrix2df_log (struct matrix2df* A, struct matrix2df* B);
float matrix2df_element_sum (struct matrix2df* A);
void matrix2df_scalar_mult (struct matrix2df* A, float val, struct matrix2df* B);
};
#endif
