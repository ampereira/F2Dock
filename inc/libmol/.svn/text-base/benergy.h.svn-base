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
#ifndef _MOL_BENERGY_H_
#define _MOL_BENERGY_H_

namespace LibMol{

/** \file benergy.h
        This file contains  functions
        for calculating bonded energies and forces
	(bonds, angles, dihedrals, impropers)
*/

/**
  find the bond energy and forces
*/
void beng(struct atomgrp *ag, double* en);

/**
  find the angle energy and forces
*/
void aeng(struct atomgrp *ag, double* en);

/**
  find the improper energy and forces
*/
void ieng(struct atomgrp *ag, double* en);

/**
  find the dihedral energy and forces
*/
void teng(struct atomgrp *ag, double* en);

void zero_grads(struct atomgrp *ag);

void check_b_grads(struct atomgrp *ag, double d,
                void (*efun)(struct atomgrp *, double*));
void check_epeng_grads(struct atomgrp *ag, struct grid* pot, double d,
                void (*efun)(double, struct atomgrp *, double*, struct grid*));
void check_speng_grads(int nstart, int nend,
                struct atomgrp *ag, double d, double stens,
                double* hx0, double* hy0, double* hz0, 
                int nx, int ny, int nz,
                double dcx, double dcy, double dcz,
                double cx, double cy, double cz, double w,
                void (*efun)(double, double, struct atomgrp *, double*,
                double*, double*, double*, 
                int,int,int,
                double,double,double,
                double,double,double,double));

};
#endif
