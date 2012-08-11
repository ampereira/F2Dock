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
#ifndef _MOL_ICHARMM_H_
#define _MOL_ICHARMM_H_

namespace LibMol{

/**
	routines for reading charmm forcefield parameters for the
	particular atom structure 
*/

/**
	return an array of charmm atom type names for the
	particular atom structure
*/
char* get_atnam_rtf(char* rtffile, int* atoms, int nat);
/**
        read charmm psf file: 
	o_nat-number of atoms, 		o_atoms-atom type numbers,
	o_nb-number of bonds, 		o_bonds-2atom numbers for each bond,
	o_nang-number of angles, 	o_ang-3atom numbers for each angle,
	o_ndih-number of dihedrals, 	o_dihs-4atom numbers for each dihedral,
	o_nimp-number of impropers, 	o_imps-4atom numbers for each improper,
	o_cha-charges for atoms
*/
void read_psf(const char* psffile, int* o_nat,  int** o_atoms, char*** o_atom_names, int* o_nb,   int** o_bonds,
                                   int* o_nang, int** o_angs,  int* o_ndih, int** o_dihs,
                                   int* o_nimp, int** o_imps, float** o_cha, 
                                   int* o_nres, int** o_ires, char** o_dres,
                                   int* o_ndon, int** o_dons, int* o_nacc, int** o_accs);
/**
	read forcefield constants from charmm papameter file
	**kbond, **lbond:-bonds
	**o_kangle, **o_fangle:-angles
	**o_kdih, **o_fdih, **o_pdih:-dihedrals (pdih-periodicity)
	**o_kimp,**o_fimp:-impropers
	**o_epa, **o_siga:-vdw parameters
*/
void read_par(char* prmfile, int nat,  char* atnam,
                int nb,   int *bonds, float **kbond, float **lbond,
                int nang, int *angles, float **o_kangle, float **o_fangle,
                int *ndih, int *dihs, float **o_kdih, float **o_fdih, int **o_pdih, int **o_tdih,
                int nimp, int *imps, float **o_kimp, float **o_fimp,
                float **o_epa, float **o_siga, float **o_acevolumes);
/**
        return *pptext, a pointer to the n-th word (starting from 0) in a string ptext
	with a delimiter ch
*/
int nthv(char **pptext, char *ptext, char ch, int n);
/**
	get bond parameters from a string buffer
*/
void chkbond(char* buffer, int nb, int* bonds,
        char* atnam, float* kbond, float* lbond);
/**
	based on existing bond i generate a string tag
	to be matched by a parameter record
*/
void bondtag(int i, int* bonds, char* atnam, char* tag);
/**
        get angle parameters from a string buffer
*/
void chkangle(char* buffer, int nang, int* angles,
        char* atnam, float* kangle, float* fangle);
/**
        based on existing angle i generate a string tag
        to be matched by a parameter record
*/
void angletag(int i, int* angles, char* atnam, char* tag);
/**
        get torsion parameters from a string buffer
*/
void chkdih(char* buffer, int *ndih, int ndih0, int* dihs,
                char* atnam, float* kdih, float* fdih, int* pdih, int *tdih, int *wdih);
/**
        based on existing torsion or improper i generate a string tag
        to be matched by a parameter record (mode is a way of dealing
	with wild cards)
*/
void dihtag(int i, int* dihs, char* atnam, int mode, char* tag);
/**
        get improper parameters from a string buffer
*/
void chkimp(char* buffer, int nimp, int* imps,
                char* atnam, float* kimp, float* fimp);
/** 
	put "X   " replacement into string p starting with character n
*/
void putXtoP(char* p, int n);
/**
        get vdw parameters from a string buffer
*/
void chkvdw(char* buffer, int nat, char* atnam, float* epa, float* siga);
/**
        get ACE volumes parmaeter from string
*/
void chkacevolumes(char* buffer, int nat, char* atnam, float* acevolumes);
/**
        main routine to read all charmm parameters for an atom group ag
*/

void read_ff_charmm(const char* psffile, char* prmfile, char* rtffile, struct atomgrp* ag);
/**
       update fixed atom index and bonded active structures. list contains indices of fixed atoms
*/
void fixed_update(struct atomgrp* ag, int nlist, int* list);
/**
	initialize active structures
*/
void fixed_init(struct atomgrp* ag);
/**
         routine to read bond parameters for an atom group ag from psf
*/
void read_ff_charmm_light(const char* psffile, struct atomgrp* ag);

};
#endif
