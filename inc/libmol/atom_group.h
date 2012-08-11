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
#ifndef _MOL_ATOM_GROUP_H_
#define _MOL_ATOM_GROUP_H_
#include <stdbool.h>

namespace LibMol{

typedef struct atomgrp mol_atom_group;

enum mol_res_type {
	UNK,
	ALA,
	ARG,
	ARGN,
	ASN,
	ASP,
	ASPH,
	CYS,
	GLN,
	GLU,
	GLUH,
	GLY,
	HIS,
	ILE,
	LEU,
	LYS,
	LYSN,
	MET,
	PHE,
	PRO,
	SER,
	SEP,
	THR,
	TRP,
	TYR,
	VAL,
	HSC,
	HSD,
	ACE,
	PCA,
	HYL,
	HYP,
	HSE,
	ORN,
	PEN,
	ALB,
	ABU,
	ST2,
	TIP3,
	OH2,
	HOH,
	WAT,
};

/**
	Holds list of atoms and
	the number of atoms in
	the list.
*/
struct atomgrp
{
	int natoms; /**< number of atoms in the group */
	mol_atom* atoms; /**< pointer to the array of atoms in the group */

        int nactives;/**< number of active atoms in the group */
        int *activelist;/**< list of active atoms in the group */

	int nbonds;
	mol_bond* bonds; /**< all first level bonds of this atomgrp */
        
        int nbact;
        struct atombond** bact; /**< array of pointers to the bonds of active atom/s */

	int nangs;
	struct atomangle* angs; /**< all angs of this atomgrp */

        int nangact;
        struct atomangle** angact; /**< all active angs of this atomgrp */

	int ntors;
	struct atomtorsion* tors; /**< all tors of this atomgrp */

        int ntoract;
        struct atomtorsion** toract; /**< all active tors of this atomgrp */

	int nimps;
	struct atomimproper* imps; /**< all imps of this atomgrp */
        
        int nimpact;
        struct atomimproper** impact; /**< all active imps of this atomgrp */

        int nres;
        int *iares;
        char **idres;
        enum mol_res_type *res_type;

        bool is_psf_read; //psf has been read in
        struct prm *prm;        
        
        void *flow_struct; // for netfork-flow based hydrogen bonding
};


// warning: this function does not handle all of mol_atom_group
// - this function mallocs
void
mol_atom_group_create_from_template (
		mol_atom_group* ag, mol_atom_group* agtmplt);


// warning: this function does not handle all of mol_atom_group
// - this function frees
void
mol_atom_group_destroy (mol_atom_group* ag);


// warning: this function does not handle all of mol_atom_group
// - this function neither mallocs nor frees
void
mol_atom_group_copy (mol_atom_group* agsrc, mol_atom_group* agdst);


// input req: minmaxs should be sizeof 6 floats
// output: minmaxs filled with min and max coords
//         of ag for every dimension
void
mol_atom_group_minmaxs (mol_atom_group* ag, float *minmaxs);


/**
	Frees all atoms in the atomgrp and the atomgrp itself.
*/
void full_free_atomgrp (struct atomgrp* ag);

/**
        Frees atomgroup light way 
*/
void free_atomgrp (struct atomgrp* ag);

/**
	Creates a copy of srcag and returns it.
*/
struct atomgrp* copy_atomgrp (struct atomgrp* srcag);

/**
        Creates a copy of srcag and returns it.
*/
struct atomgrp* fullcopy_atomgrp (struct atomgrp* srcag);

/**
  extract all atoms of type and return them in atomgrp
*/
struct atomgrp* extract_type (struct atomgrp* ag, const char* type, struct prm* prm);

/**
  remove all atoms of type and return the remaining atoms in atomgrp
*/
struct atomgrp* rm_type (struct atomgrp* ag, const char* type, struct prm* prm);

struct atomgrp* exrm_type (struct atomgrp* ag, const char* type, struct prm* prm, int direction);

/**
  typemaj extraction and removal
*/
struct atomgrp* extract_typemaj (struct atomgrp* ag, const char* typemaj, struct prm* prm);
struct atomgrp* rm_typemaj (struct atomgrp* ag, const char* typemaj, struct prm* prm);
struct atomgrp* exrm_typemaj (struct atomgrp* ag, const char* typemaj, struct prm* prm, int direction);

/**
  make all ag atoms sa
*/
void full_sa (struct atomgrp* ag);

/**
  join ags into a single atom group ag and return ag
*/
struct atomgrp* join_atomgrps (struct atomgrp** ags);
/**
  join ag1 & ag2 into a single atom group ag and return ag
*/
struct atomgrp* join_2atomgrps (struct atomgrp* ag1, struct atomgrp* ag2);

/**
	Prints the contents of the atomgrp struct to stdout
*/
void print_atomgrp (struct atomgrp* ag, struct prm* prm);

/**
  check the properties of atomgrp for consistency with prm
*/
void check_atomgrp (struct atomgrp* ag, struct prm* prm);

/**
  fill the atomgroup index in atoms
*/
void fill_ingrp (struct atomgrp* ag);

/**
  Copy coordinates from atom group to an array
*/

void ag2array( double* array, struct atomgrp* ag);

/**
  Copy coordinates from an array to the atom group
*/

void array2ag ( double* array, struct atomgrp* ag);

};

#endif
