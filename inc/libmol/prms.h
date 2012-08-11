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
#ifndef _MOL_PRMS_H_
#define _MOL_PRMS_H_

namespace LibMol{
/**
	\file prm.h
	This file contains structures and functions
	that are used to initialize prm from
	the atom parameter file.
*/

struct prm
{
	// atom
	int natoms; /**< number of atoms */
	struct prmatom* atoms;
	int nsubatoms; /**< number of subatom types */

	//float* rs;
	//float* chrgs;

	// pairwise potential
	struct prmpwpot* pwpot;

	// bond
	int nbonds;
	struct prmbond* bonds;

};

struct prmatom
{
	int id; // atom id
	char* typemaj; // e.g. GLY
	char* typemin; // e.g. CA
	int subid; // subatom id

	float r; // vdw radius
	float q; // partial charge
};

/**
  pairwise potential
*/
struct prmpwpot
{
	float r1,r2; /**< r1,r2 for potential. */
	int k; /**< number of eigenvalues */
	float* lambdas; /**< pointer to array of eigenvalues */
	float* Xs; /**< pointer to matrix of eigenvectors */
};

struct prmbond
{
	int i,j; // subatom i,j
	float k; // spring constant
	float l0; // equilibrium length
};

/**
	Current line type being read.
*/
enum ereadstate
{
	VERSION,
	ATOM,
	HYDROGEN,
	RADIUS,
	POTENTIAL,
	BOND,
};

enum ereaderr
{
	ERR_VERSION,
	ERR_ATOM,
	ERR_NAMELEN,
	ERR_HYDROGEN,
	ERR_RADIUS,
	ERR_POTENTIAL,
};

/**
	Defines atoms types, pairwise energies,
	and charges.
*/
struct prm* read_prm (const char* path, const char* bin_version);
void read_prmversion (const char* path, const char* bin_version);
void read_prmatom (struct prm* prm, const char* path);
void read_prmpwpot (struct prm* prm, const char* path);
void read_prmbond (struct prm* prm, const char* path);

/**
	Read the PDB file, and populate only the typemin and typemaj fields 
	of prm->atoms with only a single entry for each unique <typemaj, typemin> 
	tuple in the PDB. This is mainly for compatibility with other parts
	of the code already written. However, for large PDB's this will also 
	save some space.
	[ Rezaul Chowdhury, Nov. 17, 2010 ]
*/
void read_typeinfo_from_pdb (struct prm* prm, const char* path);

/**
  Return the atom id of the atypemaj, atypemin key
*/
int atomid (struct prm* prm, const char* atypemaj, const char* atypemin);

/**
	Returns a malloced copy of iprm,
	currently only copies prmatom

	\param iprms parameter struct to be copied
	\return pointer to copy of iprms
*/
struct prm* copy_prm (struct prm* iprm);
void copy_prmatom (struct prmatom* iatom, struct prmatom* oatom);

/**
	Multiply radii in prm by k.
	\param prm parameter struct to be modified
	\param k value with which to multiply radii
	\return void
*/
void modify_prms_radii (struct prm* prm, float k);

/**
	Frees the contents of prm and prm itself.
*/
/*
void free_prms (struct prm* prm);
*/

/**
	Prints an error message based on readerr.
*/
void print_readerr (enum ereaderr readerr, const char* path, char* line);

/**
	Prints the prm struct to stdout
*/
void print_prm (struct prm* prm);
void print_prmatom (struct prm* prm);
void print_prmpwpot (struct prm* prm);
void print_prmbond (struct prm* prm);

void free_prm (struct prm* prm);

};
#endif
