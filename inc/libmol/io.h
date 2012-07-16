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
#ifndef _MOL_IO_H_
#define _MOL_IO_H_


namespace LibMol{

/** \file io.h
	This file contains structures and functions
	that are used for reading from and
	writing to files
*/

typedef enum
{
	FILE_PDB,
	FILE_XYZ,
	FILE_MS, // marksur pdb
	FILE_DX, // dx grid
	FILE_POT, // potential grid
	FILE_UNKNOWN
} File_Type;

/**
  determines the file ext of path
*/
File_Type file_ext (const char* path);

/**
  prints to stderr a list of known file types
*/
void fprintf_stderr_atomgrp_file_exts ();

/**
  reads an atomgrp file based on its ext
*/
struct atomgrp* read_file_atomgrp (const char* path, struct prm* prm, float msur_k);

//struct grid* read_file_grid (const char* path);

/**
  prints an atomgrp file based on its ext
*/
void fprint_file_atomgrp (const char* path, struct atomgrp* ag, struct prm* prm);

/**
  reads modified vdw parameters from pdb file 
*/
void read_mod_vdw(char *mfile, int *nmod, int **mod, double **modeps, double **modrminh);
};
#endif
