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
#ifndef _MOL_MINIMIZE_H_
#define _MOL_MINIMIZE_H_

/** \file minimize.h
        This file contains  functions
        for local structure minimization
*/

namespace LibMol{

typedef struct {
    int     n;
    int     m;
    int     niter;       /* number of iterations so far                        */
    int     nfuns;       /* number of function evaluations so far              */
    int     iflag;
    int     diagco;
    int     iprint[2];   /* see the comment in lbfgs.f for usage of this field */
    double  eps;
    double  xtol;
    double *diag;
    double *w;
} lbfgs_t;

typedef enum
{
    MOL_LBFGS,
    MOL_CONJUGATE_GRADIENTS,
    MOL_POWELL,
} mol_min_method;


double max (double x, double y);

lbfgs_t* lbfgs_create(int n, int m, double eps);

int lbfgs_run(lbfgs_t* obj, double* x, double f, double* g);

void lbfgs_destroy(lbfgs_t* obj);

void linestep ( double* stx , double* fx , double* dx , double* sty , double* fy , double* dy , double* stp , double fp , double dp , short* brackt , double stpmin , double stpmax , int* info );

double ddot ( int n, double* dx, int ix0, int incx, double* dy, int iy0, int incy );

void lb1 ( int* iprint , int* iter , int* nfun , double gnorm , int n , int m , double* x , double f , double* g , double* stp , short finish );

void daxpy ( int n , double da , double* dx , int ix0, int incx , double* dy , int iy0, int incy);

//LBFGS minimizer for arbitrary function 
void lbfgsMin(double* orig, unsigned int maxIt, double tol, int ndim,  void* prms,void (*egfun)(int , double* , void* , double* , double*),double* min, double* fmim );
//Wrapper for atom group minimizers;
void minimize_ag(mol_min_method min_type, unsigned int maxIt, double tol, struct atomgrp* ag,void* minprms, void (*egfun)(int , double* , void* , double* , double*));

// Powell/dirPowell below

void bracket(double* orig, double* dir, double step,
             int ndim, void* prms,
             void (*egfun)(int , double* , void* , double* , double*),
             double *fb,  double *la, double *lb, double *lc);

void brent(double* orig, double* dir,
           double fb, double la, double lb, double lc,
           int ndim, void* prms,
           void (*egfun)(int , double* , void* , double* , double*),
           double tol, unsigned int maxtimes,
           double*  min, double* fmim);

void dirbrent(double* orig, double* dir,
           double fb, double la, double lb, double lc,
           int ndim, void* prms,
           void (*egfun)(int , double* , void* , double* , double*),
           double tol, int maxtimes,
           double*  min, double* fmim, double*  grad);

void powell(double* orig, double* directions, unsigned int maxIt, double tol,
            int ndim, void* prms,
            void (*egfun)(int , double* , void* , double*, double* ),
            double* min, double* fmim);

void dirMin(double* orig, unsigned int maxIt, double tol,
           int ndim, void* prms,
           void (*egfun)(int , double* , void* , double* , double*),
           double* min, double* fmim);

void limin(double* orig, double* dir, unsigned int maxIt, double tol,
           int ndim, void* prms,
           void (*egfun)(int , double* , void* , double* , double*),
           double* min, double* fmim);

void old_brent(double* orig, double* dir,
           double fb, double la, double lb, double lc,
           int ndim, void* prms,
           void (*egfun)(int , double* , void* , double*, double* ),
           double tol, int maxtimes,
           double*  min, double* fmim);
// Powell/dirPowell above
};
#endif
