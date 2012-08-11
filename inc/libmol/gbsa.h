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
#ifndef _MOL_GBSA_H_
#define _MOL_GBSA_H_

namespace LibMol{

#define L_PI           3.141592653589793  /* CHARMM pi */
struct acesetup {
    int ntypes;
    double efac;
    int nbsize;//Size of nblist
    int* list0123;
    int n0123;
    double* eself;
    double* rborn;
    double* swarr;
    double* dswarr;
    double* darr;
    //d(rb)/d(eself)
    double* dbrdes;
    double* xsf;
    double* ysf;
    double* zsf;
    double* xf;
    double* yf;
    double* zf;
    double* diarr;
    double  *lwace,*rsolv,*vsolv,*s2ace,*uace,*wace,*hydr;
   
};
//Initialize ace data types
void ace_ini(struct atomgrp* ag,struct acesetup* ac_s);
//Update ace lists once fixedlist was updated
void ace_fixedupdate(struct atomgrp* ag,struct agsetup* ags ,struct acesetup* ac_s);
//Update nblst once nblist is updated
void ace_updatenblst(struct agsetup* ags, struct acesetup* ac_s);
//Calculate ace energy and gradients
void aceeng(struct atomgrp* ag,double *en,struct acesetup* ac_s,struct agsetup* ags);
//Free ace data
void destroy_acesetup(struct acesetup* ac_s);
void free_acesetup(struct acesetup* ac_s);
//Test ace gradients
void test_acegrads(struct atomgrp *ag,struct agsetup* ags, struct acesetup* acs,double d);
};
#endif
