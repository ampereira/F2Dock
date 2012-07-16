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
#ifndef _MOL_MOVE_H_
#define _MOL_MOVE_H_

namespace LibMol{
/** \file move.h
	This file contains structures and functions
	that are used to translate the coordinates
	of the xyz to its center of mass.
	This file also contains structures and functions
	that are used to do an euler rotation of
	the xyz.
*/

struct rmatrix
{
	float a11; float a12; float a13;
	float a21; float a22; float a23;
	float a31; float a32; float a33;
};

/**
	This holds the desired rotations that
	are defined in the rotation parameter
	file.
*/
struct rotset
{
	int nrots;
	struct rmatrix* rmats;
};


struct axis
{
	float a,b,c;
};

struct axisset
{
	int naxes;
	struct axis* axes;
};

/**
	Calculates the distance between endpoints of 2 vectors,
	otherwise known as distance between 2 points
	\param va the first vector
	\param vb the second vector
	\return float of the distance
*/
float calc_distance (struct tvector* va, struct tvector* vb);

/**
	Calculates the center of mass of each
	dimension and returns the result in the
	translate_vector struct.
	\param xyz the xyz whose center of mass is to be calculated
	\return the tvector center of mass result
*/
struct tvector* center_of_mass (struct atomgrp* xyz);

/**
	Calculates the center of mass of each
	dimension and returns the negative of the
	result in the translate_vector struct.
	This tvector can then be sent to the
	translate_coords function to center the
	protein.
	\param xyz the xyz whose center of mass is to be calculated
	\return the tvector center of mass result
*/
struct tvector* negative_center_of_mass (struct atomgrp* xyz);


struct tvector* center_of_extrema (struct atomgrp* ag);
struct tvector* negative_center_of_extrema (struct atomgrp* ag);

/**
  Returns a malloced negative of tv
*/
struct tvector* negative_tvector (struct tvector* tv);

/**
	Translates by the tvec, all coords in xyz.
	\param xyz xyz to be translated
	\param tvec the tvector defining the translation
*/
void translate_atomgrp (struct atomgrp* src, struct atomgrp* dest, struct tvector* tvec);

/**
	Reads the contents of the rotation
	matrix parameter file and fills the
	struct rmatrix accordingly.
	\return a pointer to the filled rmatrix
*/
struct rotset* read_rotset (const char* path);

/**
  axis functions for symmetry constraints
*/
struct axisset* read_axisset (const char* path);
void print_axisset (struct axisset* axisset);

/**
  create a rotset for ligand monomer
*/
struct rotset* create_sym_rotset (int nmonomers);

/**
	Frees the rmatrix and its contents.
*/
void free_rotset (struct rotset* rotset);

/**
	Prints the contents of struct rmatrix.
	This is a debugging function.
	\param rmatrix the rmatrix to be printed.
*/
void print_rotset (struct rotset* rotset);
void print_rmatrix (struct rmatrix* rmatrix);

/**
	Rotate xyz prot and store the result
	in xyz rotprot.
	\param prot	the xyz to be rotated
	\param rotprot the results of the rotation are stored here
	\param rmatrix the matrix defining the rotations
	\param rotn index into rmatrix, this is the specific rotation to be made
*/
void rotate_atomgrp (struct atomgrp* src, struct atomgrp* dest, struct rmatrix* rmatrix);

/**
  Rotate and translate atomgrp.
*/
void move_atomgrp (struct atomgrp* ag, struct atomgrp* moved_ag, struct rmatrix* rmatrix, struct tvector* tv);

/**
  Center atomgrp on center of extrema.
*/
void center_atomgrp_on_center_of_extrema (struct atomgrp* ag, struct atomgrp* centered_ag);
/**
  Center atomgrp on center of mass.
*/
void center_atomgrp_on_center_of_mass (struct atomgrp* ag, struct atomgrp* centered_ag);

/**
	Print the translation vector.
*/
void print_tvector (struct tvector* t);
};
#endif
