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
#ifndef _MOL_ATOM_H_
#define _MOL_ATOM_H_


namespace LibMol{

typedef struct atom mol_atom;

//enum HBondProp;  /* properties related to hydrogen bonding -> defined in hbond.c */

// deprecate: rename atom to mol_atom
struct atom
{
	int ingrp ; /**< atom index in the atomgroup */
	int atom_typen; /**< atom type number */
	int atom_ftypen;/**< atom type number in the forcefield */
	char *name;

	int sa; /**< solvent accessible: 1 => solvent accessible, 0 => !1, -1 => undefined */
	int fixed; /**< =1 if atom is immovable , 0 otherwise*/
	int mask; /**< mask this atom: 1=>mask */
	float attl; /**< attraction level: 0 => repel,
				  1 => standard sa atom,
				  >1 => attract at this level (up to 9)
				  -1 => undefined
				 */

	double X,Y,Z; /**< X,Y,Z coordinates of atom */

	double GX,GY,GZ; /**< forces acting on an atom */

	double eps, rminh, eps03, rminh03; /**< vdw parameters of atom */
    
        double acevolume;

	double chrg; /**< electrostatic charge of atom */

        double B; /**< b-factor of an atom */
        
        int res_seq; /**< residue sequence number */
        int comb_res_seq; /**< single sequence of residue numbers combining all chains */
        char icode;  /**< insertion code */
        
        int backbone; /**< 1 if this atom is part of the backbone, 0 otherwise */

	// deprecated
	int nbonds;
	struct atombond** bonds; /**< first level bonds of this atom */
	// end deprecated

	// bond indices of this atom's bonds
	int nbondis;
	int* bondis;

	int nangs;
	struct atomangle** angs; /**< angles this atom is involved in */

	int ntors;
	struct atomtorsion** tors; /**< torsions this atom is involved in */

	int nimps;
	struct atomimproper** imps; /**< impropers this atom is involved in */
	
	int hprop; /**< properties related to hydrogen bonding */
	           /**< this should have been "enum HBondProp hprop" (HbondProp is defined in hbond.h) */
                   /**< but C/C++ standard does not seem to allow forward decelarations */
	int hybridization; /**< hybridization state of the hydrogen bond acceptor atom */
		           /**< this should have been "enum HybridizationState hybridization" (HybridizationState is defined in hbond.h) */
                           /**< but C/C++ standard does not seem to allow forward decelarations */
        int base, base2; /**< indices to the base atoms of the hbond acceptor and the donor atom (in base) for hydrogens */                   
        int octree_ptr; /**< index (ptr) to octree leaf node to which this atoms belongs */           
};

// warning: this function does not handle all of mol_atom
// - this function mallocs
void
mol_atom_create (mol_atom* a, int nbondis);

// warning: this function does not handle all of mol_atom
// - this function frees
void
mol_atom_destroy (mol_atom* a);

// warning: this function does not handle all of mol_atom
// - this function neither mallocs nor frees
void
mol_atom_copy (mol_atom* asrc, mol_atom* adst);


/**********************************************************
the functions and types below should be deprecated or
moved to separate files
**********************************************************/



/**
  	angle struct
	a0
	 \
	  \
	   a1---a2

	e = k * (th - th0)^2
*/
struct atomangle
{
	struct atom *a0, *a1, *a2; /**< atoms involved in the angle */
	float th; /**< angle theta */

	float th0; /**< equilibrium angle theta */
	float k; /**< spring constant */
};

/**
  	torsion struct
	a0         a3
	 \         /
	  \       /
	   a1---a2

	e = k * (1 + cos (n*chi - d))
*/
struct atomtorsion
{
	struct atom *a0, *a1, *a2, *a3; /**< atoms involved in the torsion */
	float chi; /**< angle chi: angle between the planes a0,a1,a2 and a1,a2,a3 */

	float k; /**< k constant */
	float d; /**< delta constant */
	int n; /**< number of bonds made by atoms a1 and a2 */
};

/**
  	improper struct
	 a1---a0---a3
	      |
		  |
		  a2

	e = k * (psi - psi0)^2
*/
struct atomimproper
{
	struct atom *a0, *a1, *a2, *a3; /**< atoms involved in the improper */
	float psi; /**< angle psi: angle between the planes a0,a1,a2 and a1,a2,a3 */

	float k; /**< spring constant */
	float psi0; /**< equilibrium angle psi */
};

/**
        spring
*/
struct spring
{
        struct atomgrp *agp;  /**< affected atomgroup */
        int naspr ; /**< number of affected atoms */
        int *laspr; /**< list of atoms */
        double fkspr; /**< force constant */
        double X0,Y0,Z0; /**< anchor point */
};

struct springset
{
        int nsprings; /**< number of springs */
        struct spring *springs;  /**< array of springs */
};

/**
        frees a spring set
*/
void free_springset(struct springset *sprst);

/**
  Copies contents of atom struct src to dest.
  The pointers must be to previously allocated memory.
*/
void copy_atom (struct atom* src, struct atom* dest);

/**
  Copies contents of atom struct src to dest.
  The pointers must be to previously allocated memory.
*/
void fullcopy_atom (struct atom* src, struct atom* dest);

/**
	Prints to stderr the contents of the atom struct
*/
void fprintf_stderr_atom (struct atom* atom, struct prm* prm);

/**
 Free bonds, angles, tors, impropers pointers
*/
void free_atom(struct atom* atm);

};

#endif
