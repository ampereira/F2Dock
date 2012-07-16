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
#ifndef _HBOND_H_
#define _HBOND_H_

#include <cmath>
#include <math.h>
//#include <gmp.h>
#include <limits.h>

namespace LibMol{

//#define USE_LONG_DOUBLE
#define NON_POSITIVE_POLY

#ifdef USE_LONG_DOUBLE
   #define FLOAT long double
#else
   #define FLOAT double
#endif   

#define LINEAR_FADE  1
#define LOG_FADE     2
#define BSPLINE_FADE 3

#ifndef FADING_FUNCTION 
  #define FADING_FUNCTION BSPLINE_FADE
#endif  

// FROM ROSETTA -- FROM ROSETTA -- FROM ROSETTA -- FROM ROSETTA -- FROM ROSETTA
//
#define SHORT_RANGE_CUTOFF            4
#define GENERIC_SHORTRANGE_SEQ_CUTOFF 1

#define MAX_R             3.0
#define MIN_R             1.4 /* AH distance */
#define MIN_xH            0.0 /* cos( radians( 180.0 - 90.0 ) ): psi cutoff */
#define MIN_xD            0.0 /* cos( radians( 180.0 - 90.0 ) ): theta cutoff */
#define MAX_xH            1.0 /* cos( radians( 180.0 - 180.0 ) ): psi cutoff */
#define MAX_xD            1.0 /* cos( radians( 180.0 - 180.0 ) ): theta cutoff */
#define R_INTERP_EDGE     2.1
#define ANGLE_INTERP_EDGE 0.05

#define MIN_R2 ( MIN_R * MIN_R )
#define MAX_R2 ( MAX_R * MAX_R )

#define MAX_HB_ENERGY     0.0 /* at and above this cutoff, not considered a hbond */

#define SWITCH_DIS        2.1 /* distance cutoff for sidechain hbonds tables */
#define INTERP_RANGE      0.2
 // interpolate if AH dist falls between switch +/- range
#define INTERP_MIN ( SWITCH_DIS - INTERP_RANGE )
#define INTERP_MAX ( SWITCH_DIS + INTERP_RANGE )
//
// FROM ROSETTA -- FROM ROSETTA -- FROM ROSETTA -- FROM ROSETTA -- FROM ROSETTA


// FROM ROSETTA -- FROM ROSETTA -- FROM ROSETTA -- FROM ROSETTA -- FROM ROSETTA
//
// JSS Useful #defines to set up hbonds_geom scoring polynomials of the form
// name(/*in*/ double x, /*out*/ double &value, double &deriv);
// You shouldn't have to look at these; The idea is to auto-generate
// the calls to create_poly#() using matlab scripts to fit polynomials.
// POLY_CLIPVAL is the polynomial value outside of clip interval [xmin xmax]
//
// FROM ROSETTA -- FROM ROSETTA -- FROM ROSETTA -- FROM ROSETTA -- FROM ROSETTA
//
#define POLY_CLIPVAL 0.0
#define create_poly8(name, xmin, xmax, c7, c6, c5, c4, c3, c2, c1, c0) inline void name(FLOAT x, FLOAT *value, FLOAT *deriv) { if ( (x<=xmin) || (x>=xmax) ) {(*value) = POLY_CLIPVAL; (*deriv) = 0.0; return; } (*value) = (*deriv) = c7; (*value) = (*value)*x+c6; (*deriv) = (*deriv)*x+(*value); (*value) = (*value)*x+c5; (*deriv) = (*deriv)*x+(*value); (*value) = (*value)*x+c4; (*deriv) = (*deriv)*x+(*value); (*value) = (*value)*x+c3; (*deriv) = (*deriv)*x+(*value); (*value) = (*value)*x+c2; (*deriv) = (*deriv)*x+(*value); (*value) = (*value)*x+c1; (*deriv) = (*deriv)*x+(*value); (*value) = (*value)*x+c0; }
#define create_poly7(name, xmin, xmax, c6, c5, c4, c3, c2, c1, c0) inline void name(FLOAT x, FLOAT *value, FLOAT *deriv) { if ( (x<=xmin) || (x>=xmax) ) {(*value) = POLY_CLIPVAL; (*deriv) = 0.0; return; } (*value) = (*deriv) = c6; (*value) = (*value)*x+c5; (*deriv) = (*deriv)*x+(*value); (*value) = (*value)*x+c4; (*deriv) = (*deriv)*x+(*value); (*value) = (*value)*x+c3; (*deriv) = (*deriv)*x+(*value); (*value) = (*value)*x+c2; (*deriv) = (*deriv)*x+(*value); (*value) = (*value)*x+c1; (*deriv) = (*deriv)*x+(*value); (*value) = (*value)*x+c0; }
#define create_poly6(name, xmin, xmax, c5, c4, c3, c2, c1, c0) inline void name(FLOAT x, FLOAT *value, FLOAT *deriv) { if ( (x<=xmin) || (x>=xmax) ) {(*value) = POLY_CLIPVAL; (*deriv) = 0.0; return; } (*value) = (*deriv) = c5; (*value) = (*value)*x+c4; (*deriv) = (*deriv)*x+(*value); (*value) = (*value)*x+c3; (*deriv) = (*deriv)*x+(*value); (*value) = (*value)*x+c2; (*deriv) = (*deriv)*x+(*value); (*value) = (*value)*x+c1; (*deriv) = (*deriv)*x+(*value); (*value) = (*value)*x+c0; }
#define create_poly5(name, xmin, xmax, c4, c3, c2, c1, c0) inline void name(FLOAT x, FLOAT *value, FLOAT *deriv) { if ( (x<=xmin) || (x>=xmax) ) {(*value) = POLY_CLIPVAL; (*deriv) = 0.0; return; } (*value) = (*deriv) = c4; (*value) = (*value)*x+c3; (*deriv) = (*deriv)*x+(*value); (*value) = (*value)*x+c2; (*deriv) = (*deriv)*x+(*value); (*value) = (*value)*x+c1; (*deriv) = (*deriv)*x+(*value); (*value) = (*value)*x+c0; }
#define create_poly4(name, xmin, xmax, c3, c2, c1, c0) inline void name(FLOAT x, FLOAT *value, FLOAT *deriv) { if ( (x<=xmin) || (x>=xmax) ) {(*value) = POLY_CLIPVAL; (*deriv) = 0.0; return; } (*value) = (*deriv) = c3; (*value) = (*value)*x+c2; (*deriv) = (*deriv)*x+(*value); (*value) = (*value)*x+c1; (*deriv) = (*deriv)*x+(*value); (*value) = (*value)*x+c0; }
#define create_poly3(name, xmin, xmax, c2, c1, c0) inline void name(FLOAT x, FLOAT *value, FLOAT *deriv) { if ( (x<=xmin) || (x>=xmax) ) {(*value) = POLY_CLIPVAL; (*deriv) = 0.0; return; } (*value) = (*deriv) = c2; (*value) = (*value)*x+c1; (*deriv) = (*deriv)*x+(*value); (*value) = (*value)*x+c0; }
#define create_poly2(name, xmin, xmax, c1, c0) inline void name(FLOAT x, FLOAT *value, FLOAT *deriv) { if ( (x<=xmin) || (x>=xmax) ) {(*value) = POLY_CLIPVAL; (*deriv) = 0.0; return; } (*value) = (*deriv) = c1; (*value) = (*value)*x+c0; }
//
// FROM ROSETTA -- FROM ROSETTA -- FROM ROSETTA -- FROM ROSETTA -- FROM ROSETTA


// FROM ROSETTA -- FROM ROSETTA -- FROM ROSETTA -- FROM ROSETTA -- FROM ROSETTA
//
// define polynomials for use in hbond score calculation
//
#ifdef NON_POSITIVE_POLY

create_poly8(POLY_AHdisBBHelix, 1.78218633, 2.6757, // these are the roots to eliminate positive values.
             12.93768086,-221.0155722,1604.391304,-6409.335773,15200.86425,-21375.00216,16475.98811,-5361.55644)
create_poly8(POLY_AHdisBBOther, 1.6971523, 2.679339, // non-positive interval
             13.58980244,-224.0452428,1568.933094,-6044.257847,13820.1498,-18730.96076,13912.92238,-4361.995425)
create_poly5(POLY_AHdisSP2, 1.6941, 2.5, // non-positive interval
             10.98727738,-100.2401419,340.9733405,-511.6111233,285.0061262)
create_poly5(POLY_AHdisSP3, 1.755, 2.521385, // non-positive interval
             7.011735538,-68.99968829,251.820931,-403.3593133,238.7378958)
create_poly8(POLY_xDBBHelix, 0.3746, 1.04,
             223.5268153,-757.7254095,1019.593508,-689.2232431,240.1436064,-37.84119583,0.85868904,0.278181985)
create_poly8(POLY_xDBBOther, 0.76, 1.09,
             111.9877946,-380.3066184,514.7650204,-352.4092342,124.6219703,-19.94401946,0.149314979,0.635771774)
create_poly3(POLY_xDSP2short, 0.7071, 1.01,
             -0.562582503,-0.746682668,0.809265171)
create_poly3(POLY_xDSP2long, 0.0, 1.01,
             0.094962885,-0.254313172,0.0)
create_poly3(POLY_xDSP3short, 0.61566, 1.01,
             -0.100140144,-1.139139041,0.739279186)
create_poly3(POLY_xDSP3long, 0.0, 1.01,
             0.089380221,-0.207503776,0.0)
create_poly8(POLY_xHBBHelix, 0.156, 1.03,
             54.80664331,-196.8196655,295.9418886,-232.105602,96.99124565,-20.60918361,1.573169816,0.000745458)
create_poly8(POLY_xHBBOther, 0.61566, 1.07,
             43.94483847,-144.3836033,193.5865176,-132.4469355,47.28137288,-8.945888012,-0.227035135,0.791902995)
create_poly3(POLY_xHSP2short, 0.0, 1.08,
             1.720984644,-1.855254573,0.0)
create_poly3(POLY_xHSP2long, 0.0, 1.01,
             0.439598249,-0.444673076,0.0)
create_poly3(POLY_xHSP3, 0.0, 1.06,
             1.761487842,-1.876959406,0.0)
create_poly7(POLY_xHRing, 0.7608, 1.089,
             37.744316,-117.731674,143.0759275,-86.2258835,26.7448175,-4.4699705,0.6458455)

#else

create_poly8(POLY_AHdisBBHelix, MIN_R, 2.8, // 1.78218633, 2.6757, // these are the roots to eliminate positive values.
             12.93768086,-221.0155722,1604.391304,-6409.335773,15200.86425,-21375.00216,16475.98811,-5361.55644)
create_poly8(POLY_AHdisBBOther, MIN_R, 2.745, // 1.6971523, 2.679339, // non-positive interval
             13.58980244,-224.0452428,1568.933094,-6044.257847,13820.1498,-18730.96076,13912.92238,-4361.995425)
create_poly5(POLY_AHdisSP2, MIN_R, 2.5, // 1.6941, 2.5, // non-positive interval
             10.98727738,-100.2401419,340.9733405,-511.6111233,285.0061262)
create_poly5(POLY_AHdisSP3, MIN_R, 2.5, // 1.755, 2.521385, // non-positive interval
             7.011735538,-68.99968829,251.820931,-403.3593133,238.7378958)
create_poly8(POLY_xDBBHelix, MIN_xD, MAX_xD, // 0.3746, 1.04,
             223.5268153,-757.7254095,1019.593508,-689.2232431,240.1436064,-37.84119583,0.85868904,0.278181985)
create_poly8(POLY_xDBBOther, MIN_xD, MAX_xD, // 0.76, 1.09,
             111.9877946,-380.3066184,514.7650204,-352.4092342,124.6219703,-19.94401946,0.149314979,0.635771774)
create_poly3(POLY_xDSP2short, MIN_xD, MAX_xD, // 0.7071, 1.01,
             -0.562582503,-0.746682668,0.809265171)
create_poly3(POLY_xDSP2long, MIN_xD, MAX_xD, // 0.0, 1.01,
             0.094962885,-0.254313172,0.0)
create_poly3(POLY_xDSP3short, MIN_xD, MAX_xD, // 0.61566, 1.01,
             -0.100140144,-1.139139041,0.739279186)
create_poly3(POLY_xDSP3long, MIN_xD, MAX_xD, // 0.0, 1.01,
             0.089380221,-0.207503776,0.0)
create_poly8(POLY_xHBBHelix, MIN_xH, MAX_xH, // 0.156, 1.03,
             54.80664331,-196.8196655,295.9418886,-232.105602,96.99124565,-20.60918361,1.573169816,0.000745458)
create_poly8(POLY_xHBBOther, MIN_xH, MAX_xH, // 0.61566, 1.07,
             43.94483847,-144.3836033,193.5865176,-132.4469355,47.28137288,-8.945888012,-0.227035135,0.791902995)
create_poly3(POLY_xHSP2short, MIN_xH, MAX_xH, // 0.0, 1.08,
             1.720984644,-1.855254573,0.0)
create_poly3(POLY_xHSP2long, MIN_xH, MAX_xH, // 0.0, 1.01,
             0.439598249,-0.444673076,0.0)
create_poly3(POLY_xHSP3, MIN_xH, MAX_xH, // 0.0, 1.06,
             1.761487842,-1.876959406,0.0)
create_poly7(POLY_xHRing, MIN_xH, MAX_xH, // 0.7608, 1.089,
             37.744316,-117.731674,143.0759275,-86.2258835,26.7448175,-4.4699705,0.6458455)

#endif             
//
// FROM ROSETTA -- FROM ROSETTA -- FROM ROSETTA -- FROM ROSETTA -- FROM ROSETTA


// FOLLOWING ROSETTA -- FOLLOWING ROSETTA -- FOLLOWING ROSETTA -- FOLLOWING ROSETTA
//
#define create_fade_interval( name, min0, fmin, fmax, max0 ) \
inline void name##_value_deriv(FLOAT x, FLOAT *val, FLOAT *deriv) { \
     const FLOAT name##_dfade_min = ( fmin == min0 ) ? 0.0 : ( 1.0 / ( fmin - min0 ) ), name##_dfade_max = ( max0 == fmax ) ? 0.0 : ( 1.0 / ( max0 - fmax ) ); \
     if ( x < fmax ) { if ( x >= fmin ) { *val = 1.0; *deriv = 0.0; return; } \
                       if ( x <= min0 ) { *val = *deriv = 0.0; return; } \
                       *deriv = name##_dfade_min; *val = ( x - min0 ) * name##_dfade_min; \
                     } \
     else { if ( x >= max0 ) { *val = *deriv = 0.0; return; } \
 	    *deriv = -name##_dfade_max; *val = ( max0 - x ) * name##_dfade_max; } \
}
//
// FOLLOWING ROSETTA -- FOLLOWING ROSETTA -- FOLLOWING ROSETTA -- FOLLOWING ROSETTA


// FOLLOWING ROSETTA -- FOLLOWING ROSETTA -- FOLLOWING ROSETTA -- FOLLOWING ROSETTA
//

// used to adjust xD,xH
create_fade_interval( fade_rBB, MIN_R, MIN_R /*+ 0.1*/, R_INTERP_EDGE, MAX_R )

create_fade_interval( fade_rshort, MIN_R, MIN_R /*+ 0.1*/, INTERP_MIN, INTERP_MAX )
create_fade_interval( fade_rlong, INTERP_MIN, INTERP_MAX, INTERP_MAX, MAX_R )

// xD=theta should fade r,xH sooner!
create_fade_interval( fade_xD, MIN_xD, ANGLE_INTERP_EDGE, 1, 1 )

// fades r,xD
create_fade_interval( fade_xH, MIN_xH, ANGLE_INTERP_EDGE, 1, 1 )
//
// FOLLOWING ROSETTA -- FOLLOWING ROSETTA -- FOLLOWING ROSETTA -- FOLLOWING ROSETTA


#ifdef USE_LONG_DOUBLE
  #define create_log_fade_interval( name, a, b, c, d, e, f ) \
  inline void name##_log_value_deriv( FLOAT x, FLOAT *val, FLOAT *deriv ) { \
       FLOAT lnx = logl( x ); \
       *val = ( ( ( ( f * lnx + e ) * lnx + d ) * lnx + c ) * lnx + b ) * lnx + a; \
       *deriv = ( ( ( ( 5 * f * lnx + 4 * e ) * lnx + 3 * d ) * lnx + 2 * c ) * lnx + b ) / x; \
  }
#else
  #define create_log_fade_interval( name, a, b, c, d, e, f ) \
  inline void name##_log_value_deriv( FLOAT x, FLOAT *val, FLOAT *deriv ) { \
       FLOAT lnx = log( x ); \
       *val = ( ( ( ( f * lnx + e ) * lnx + d ) * lnx + c ) * lnx + b ) * lnx + a; \
       *deriv = ( ( ( ( 5 * f * lnx + 4 * e ) * lnx + 3 * d ) * lnx + 2 * c ) * lnx + b ) / x; \
  }     
#endif

// used to adjust xD,xH
create_log_fade_interval( fade_rBB, -37.3180547817287, 285.741141861363, -834.307606917721, 1189.6798830102, -826.409671561707, 222.919001825369 )

create_log_fade_interval( fade_rshort, 5.17967058862565, -73.0669137757406, 315.57016260663, -559.689100793323, 438.586611628675, -126.548509998553 )
create_log_fade_interval( fade_rlong, -14.0813542571684, 130.05181593123, -453.153559518528, 740.025843886255, -562.47374254511, 160.168566150279 )

// xD=theta should fade r,xH sooner!
create_log_fade_interval( fade_xD, 0.133459218358311, 13.0816902916604, -57.2491258523768, 86.8476548565786, -3.78050184877689, -58.7848249871094 )

// fades r,xD
create_log_fade_interval( fade_xH, 0.133459218358311, 13.0816902916604, -57.2491258523768, 86.8476548565786, -3.78050184877689, -58.7848249871094 )


enum HBondProp { UNKNOWN_HPROP = 0, HBOND_DONOR = 1, HBOND_ACCEPTOR = 2, DONATABLE_HYDROGEN = 4, ROTATABLE_HYDROGEN = 8, FIXED_HYDROGEN = 16 };

// FOLLOWING ROSETTA -- FOLLOWING ROSETTA -- FOLLOWING ROSETTA -- FOLLOWING ROSETTA
//
enum HybridizationState { UNKNOWN_HYBRID = 0, SP1_HYBRID, SP2_HYBRID, SP3_HYBRID, RING_HYBRID };

// donor classes
enum HBDonChemType
{
	hbdon_NO = 0,
	hbdon_BB,
	hbdon_SC,
	hbdon_MAX
};


// acceptor classes
enum HBAccChemType
{
	hbacc_NO = 0,
	hbacc_BB,
	hbacc_SP2,
	hbacc_SP3,
	hbacc_RING,
	hbacc_MAX
};

enum HBEvalType
{
	hbe_NONE = 0,  /* no hbond */
	/* backbone/backbone (BB/BB) -- order must be preserved */
	hbe_BB,        /* generic BB/BB (only refined are used) */
	hbe_BBTURN,    /* BB/BB short range */
	hbe_BBHELIX,   /* BB/BB short range */
	hbe_BBOTHER,   /* keep BB/BB consecutive for hbond_compute_energy */
	/* backbone (BB) donor, sidechain (SC) acceptor */
	hbe_SP2B,
	hbe_SP3B,
	hbe_RINGB,
	/* BB acceptor, SC donor */
	hbe_BSC,
	/* sidechain/sidechain (SC/SC) */
	hbe_SP2SC,
	hbe_SP3SC,
	hbe_RINGSC,
	/* hbe_MAX should always be # of last */
	hbe_MAX
};


enum HBondWeightType 
{
	hbw_NONE = 0,
	hbw_SR_BB,
	hbw_LR_BB,
	hbw_BB_SC,
	hbw_SC,
	hbw_H2O
};
//
// FOLLOWING ROSETTA -- FOLLOWING ROSETTA -- FOLLOWING ROSETTA -- FOLLOWING ROSETTA


struct dvector
{
   FLOAT X; /**< translation in the X dimension */
   FLOAT Y; /**< translation in the Y dimension */
   FLOAT Z; /**< translation in the Z dimension */
};


#define ARY( a, n, i, j ) a[ ( i ) * ( n ) + ( j ) ]

typedef struct
{
    int hydro_id, acc_id;
    double en;
} FLOW_EDGE;


typedef struct
{
    int id, cap;
} FLOW_MAP;


typedef struct
{
    int max_n_fedge, max_n;
    FLOW_EDGE *flow_edge;
    FLOW_MAP *flow_map;
    int *cap;
    FLOAT *cost;
    int *deg;
    int *par;
    int *q;
    int *inq;
    FLOAT *pi;
    FLOAT *d;
    int *adj;
    int *fnet;
} FLOW_STRUCT;


int init_flow_struct( FLOW_STRUCT *flow_struct );
void free_flow_struct( FLOW_STRUCT *flow_struct );
int reinit_flow_struct( FLOW_STRUCT *flow_struct, int max_n_fedge, int max_n );

void mark_hbond_donors (
		 struct atomgrp* ag,
		 struct prm* prm
		);

void mark_hbond_acceptors (
		 struct atomgrp* ag,
		 struct prm* prm
		);

void fix_acceptor_bases( struct atomgrp *ag, struct prm *prm );

void hbondeng( struct atomgrp *ag, double *energy, struct nblist *nblst );

void hbondengcat( struct atomgrp *ag, double *energy, struct nblist *nblst );

void hbondeng_all( struct atomgrp *ag, double *energy, struct nblist *nblst );

void flow_hbondeng( struct atomgrp *ag, double *energy, struct nblist *nblst );

void hbondeng_octree_single_mol( OCTREE_PARAMS *octpar, double *energy );

double *alloc_categorized_hbondeng( void );

void init_categorized_hbondeng( double *engcat );

void print_categorized_hbondeng( double *engcat );

void get_categorized_hbondeng( double *engcat, double *bb_bb_sr, double *bb_bb_lr, double *bb_sc, double *sc_sc );

void set_categorized_hbondeng( double *engcat, double bb_bb_sr, double bb_bb_lr, double bb_sc, double sc_sc );

void print_acceptor_hybridization_states( struct atomgrp *ag, struct prm *prm );

void check_fading_funcs( void );

};

#endif
