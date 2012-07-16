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
#ifndef _MOL_PDB_H_
#define _MOL_PDB_H_

namespace LibMol{
/**
  generate a single residue sequence combining all chains
*/
void assign_combined_residue_sequence( struct atomgrp *ag );

/**
  read a pdb file
*/
struct atomgrp* read_pdb (const char* path, struct prm* prm);

/**
  read a pdb file
*/
struct atomgrp* read_pdb_nopar (const char* path);

/**
	Read the PDB file, but assign an integer id to each atom
	based on <residue-name, atom-name> instead of explicitly 
	storing the name of the atom and the residue it is part of.
	This is mainly for compatibility with other parts of the 
	code already written. However, for large PDB's this will also 
	save some space.
	[ Rezaul Chowdhury, Nov. 17, 2010 ]
*/
struct atomgrp* read_pdb_with_compressed_typeinfo (const char* path, struct prm* prm);

/**
  read pdb model files
  rmodels specifies how much models is required, if set to -1    
*/
struct atomgrp** read_pdb_models (const char* path, struct prm* prm, int* rmodels);

/**
 * read pdb model files without requiring the atoms prms file
 * \return an array of pointers to atomgrp with models
 * @param[in] path location of the models file
 * @param[in,out] rmodels number of models required from file, if -1, read all, 
 * updated to number of models read in
 */
struct atomgrp **read_pdb_modelsnopar(const char *path, int *rmodels);

/**
  print a pdb file
*/
void fprint_pdb (struct atomgrp* ag, struct prm* prm, const char* path);

/**
  print a pdb file
*/
void write_pdb_nopar (struct atomgrp* ag, const char* inf, const char* ouf);
};
#endif
