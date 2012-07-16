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
#ifndef _MOL_XYZ_H_
#define _MOL_XYZ_H_

namespace LibMol{
/**
	\file xyz.h
	This file contains structures and functions
	that are used to initialize the atomgrp
	structure from the xyz file.
*/

/**
	Returns a pointer to an xyz struct,
	read from xyz file at path.
	Exits if file does not exist.
*/
struct atomgrp* read_xyz (const char* path, struct prm* prm);
/**
  reads xyz file with extra sa column.
*/
struct atomgrp* read_oldxyz (const char* path, struct prm* prm);

/**
	Prints an xyz file from the xyz coordinates.
*/
void print_xyz (struct atomgrp* xyz, struct prm* prm);
void fprint_xyz (struct atomgrp* xyz, struct prm* prm, const char* path);
void fprint_oldxyz (struct atomgrp* xyz, struct prm* prm, const char* path);

};
#endif
