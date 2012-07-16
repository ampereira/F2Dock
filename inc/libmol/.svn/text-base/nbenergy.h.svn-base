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
#ifndef _MOL_NBENERGY_H_
#define _MOL_NBENERGY_H_

namespace LibMol{
/** \file nbenergy.h
        This file contains  functions
        for calculating nonbonded energies and forces
	and mantaining nonbonded lists
*/
#define small 0.0000001

struct cluster
{
	int natoms; /**<number of atoms in the cluster */
        double gcent[3]; /**<coordinates of the geom center */
        double mdc;      /**<maximun distance to the center */
        int nextincube; /**<index of the next cluster in cube */
        int* iatom;/**<atom index of the cluster */
};

struct clusterset
{
        int nclusters; /**< number of clusters in the coor set */
        struct cluster* clusters; /**< pointer to the array of clusters */
        double marg; /**< sum of two max distances from the geom center */
};

struct cube
{
        int hstincube; /**< index of the highest numbered cluster in cube */
	int nncubes; /**< number of neighboring cubes */
        int *icubes; /**< neighboring cubes indices */
        int ix, iy, iz; /**< indices of the cube */
};

struct cubeset
{
        int ncubes; /**< number of cubes for the cluster set */
        int nfcubes; /**<number of filled cubes for the cluster set */
        int* ifcubes; /**<index of filled cubes */
        struct cube* cubes; /**< pointer to the array of cubes */
        double cubel; /**< cube length */
};

struct nblist
{
        int npairs;   /**< number of pairs in nblist */
        int nfat;     /**< number of first atoms in pairs */
        int *ifat;    /**< index of first atom */
        int *nsat;    /**< number of second atoms for each first */
        int **isat;   /**< pointer to the array of second atoms */
	double nbcut; /**< nonbond list cutoff length */
        double nbcof; /**< forcefield cutoff length */
        float *crds;  /**< coordinates at nblist generation */
};
//Wrapper for nonbonded list
struct agsetup{
    int n02, n03,ndm,nf03;
    int *list02;
    int *list03;
    int *listf03;//Fixed 03 list
    int *excl_list;
    int** pd1;
    int** pd2;
    struct nblist *nblst;
    struct clusterset *clst;
};

void destroy_agsetup(struct agsetup* ags);
void free_agsetup(struct agsetup* ags);

void test_nbgrads(struct atomgrp *ag, double d, struct nblist *nblst,
                  int n03, int* list03);

void vdweng03(double f, struct atomgrp *ag, double* ven, int n03, int* list03);

void vdwengs03(double f, double rc, struct atomgrp *ag, double* ven,
               int n03, int* list03);

void vdweng(struct atomgrp *ag, double* ven, struct nblist *nblst);

void eleng03(double f, struct atomgrp *ag, double eps, double* een,
             int n03, int* list03);

void elengs03(double f, double rc, struct atomgrp *ag, double eps, double* een,
              int n03, int* list03);

void eleng(struct atomgrp *ag, double eps, double* een, struct nblist *nblst);

void destroy_nblist(struct nblist *nblst);
void free_nblist(struct nblist *nblst);

void gen_nblist(struct atomgrp *ag, struct cubeset *cust, struct clusterset *clst, 
                int *excl_list, int **pd1, int **pd2, int ndm, struct nblist *nblst);

void free_cubeset(struct cubeset *cust);

void gen_cubeset(double nbcut, struct clusterset *clst, struct cubeset *cust);

void findmarg(struct atomgrp *ag, struct clusterset *clst);

void free_clset(struct clusterset *clst);

void gen_clset(int natoms, struct clusterset *clst, 
               int nclust, int *clust);

void clust14(struct atomgrp* ag, int *nclust, int *clust);

void addpair(int i, int k, int *npair, int* pairs);

int natpair(struct atomgrp* ag, int *pair, int *clust);

void addclust(struct atomgrp* ag, int *pair, int *clust, int nclust);

void comp_list01(struct atomgrp* ag, int *list01, int *na01, int **pna01);

void comp_n23(int natoms, int *na01, int **pna01, int *n02, int *n03);

void comp_list02(int natoms, int *na01, int **pna01, 
                 int *na02, int **pna02, int *n02, int *list02); 

void comp_list03(int natoms, int *na01, int **pna01, int *list03);

void trim_list03(int natoms, 
                 int *na01, int **pna01, 
                 int *na02, int **pna02, 
                 int *n03, int *list03);

int trim_comp(void *s1, void *s2);

void excl_dims(int natoms, int *na01, int **pna01,
               int n02, int *list02, int n03, int *list03,
               int *nd1, int *nd2, int *ndm, int *atmind);

int exta(int a1, int a2, int *excl_list, int **pd1, int **pd2, int ndm);

// assuming ndm <= 20
inline int exta2( int a1, int a2, int *excl_list, int **pd1, int ndm ) 
{
    int i = a2 - a1 - 1;
        
    if ( i >= ndm ) return 0;
    
    if ( i < 10 ) return excl_list[ ( a1 << 3 ) + ( a1 << 1 ) + i ];
    else return pd1[ a1 ][ i - 10 ];
}


void excl_tab(int natoms, int *na01, int **pna01,
               int n02, int *list02, int n03, int *list03,
               int nd1, int nd2, int ndm, int *atmind,
               int *excl_list, int **pd1, int **pd2);
/* extract non fixed 03 list */
void fix_list03(struct atomgrp *ag, int n03, int* list03, 
                                    int *nf03, int *listf03);
//Wrapper for nblist initialisation
void init_nblst(struct atomgrp* ag, struct agsetup* ags);
//Wrapper for nblist update
void update_nblst(struct atomgrp* ag, struct agsetup* ags);
//Updating nblist if atoms moved enough
int check_clusterupdate(struct atomgrp* ag,struct agsetup* ags);
/*create 012 fixed/active lists for ace*/
void give_012(struct atomgrp* ag, struct agsetup* ags,
              int* na012, int** la012, int* nf012, int** lf012);
/*springs energy*/
void springeng(struct springset *sprst, double* spren);

/* modify vdw parameters and save old values */
void modify_vdw_save(struct atomgrp* ag, int nmod, int* mod, double* enew, double* rnew,
		                                             double* eold, double* rold);
/* modify vdw parameters back to saved old values */
void modify_vdw_back(struct atomgrp* ag, int nmod, int* mod, double* eold, double* rold);

/* get 1D vdw epsilon and/or sigma arrays scaled for further vdw modification */
void get_mod_vdw_all(double lambe, double lambr, struct atomgrp* ag,
                               int *nmod, int **mod, double **modeps, double **modrminh);

/* ####################################################################### */
/* #                       START: ADDED BY REZAUL                        # */
/* ####################################################################### */


//void vdweng_octree_single_mol( OCTREE *octree_static, int node_static, OCTREE *octree_moving, int node_moving,
//                               void *params, double dist_cutoff, double *trans_mat, double *energy );

void vdweng_octree_single_mol( OCTREE_PARAMS *octpar, double *energy );
                                                                
//void eleng_octree_single_mol( OCTREE *octree_static, int node_static, OCTREE *octree_moving, int node_moving,
//                              void *params, double dist_cutoff, double *trans_mat, double *energy );
                              
void eleng_octree_single_mol( OCTREE_PARAMS *octpar, double *energy );                              
                              
/* ####################################################################### */
/* #                        END: ADDED BY REZAUL                         # */
/* ####################################################################### */

};
#endif
