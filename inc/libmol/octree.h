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
**  Structural Bioinformatics Laboratory and  Boston University. Requests to use this software (or                                                                                 
**  modified copies therof) in any other way should be sent to Dima Kozakov,                                                                                                     
**  Department of Biomedical Engineering: "midas@bu.edu".                                                                                                                  
**                                                                                                                                                                               
**---------------------------------------------------------------------------*/
#ifndef _MOL_DYN_OCTREE_H_
#define _MOL_DYN_OCTREE_H_

#define HALF_SQRT_THREE 0.86602540378443864676372317075294

#define INIT_NUM_OCTREE_NODES 100
#define LOW_BITS  14
#define LOW_MASK  0x3FFF

#define INIT_MIGRATION_ARRAY_SIZE 8

#define create_octree_ptr( a, b ) ( ( ( a ) << LOW_BITS ) + b )
#define get_node_id( c ) ( ( c ) >> LOW_BITS )
#define get_index_in_node( c ) ( ( c ) & LOW_MASK )

#ifdef zeroIfLess
   #undef zeroIfLess
#endif
#define zeroIfLess( a, b ) ( ( ( a ) < ( b ) ) ? 0 : 1 )

#ifndef RECURSION_DEPTH
   #define RECURSION_DEPTH 50
#endif


//#define ADD_ATTR

namespace LibMol{
typedef struct
{
        double lx, ly, lz;     /* the corner of the cell cube with the minimum x, y, z co-ordinates */
        double dim;            /* dimension of the cube */

#ifdef ADD_ATTR
	double sx, sy, sz;     /* sums of the x, y and z co-ordinates of the atoms under this node */
	double sq;             /* sum of charges of atoms under this node */
#endif

	int n;                 /* number of items/atoms under this node */
	int nfixed;            /* number of fixed items/atomd under this node */

	int p_ptr;             /* index of (ptr to) the parent node in the octree (-1 for the parent of the root) */
	int c_ptr[ 8 ];        /* indices of (ptrs to) child nodes */

	int leaf;              /* 1 if this node is a leaf (i.e., no further subdivision), 0 otherwise */
	int id_cap;            /* capacity of the indices array */
	int id_num;            /* number of items in the indices array (used for temporary storage in non-leafs) */
	int *indices;          /* indices to the items/atoms in this node (fixed items are listed first) */

} OCTREE_NODE;


typedef struct
{
	OCTREE_NODE *nodes;    /* nodes of the octree with nodes[ 0 ] being the root */

	int max_leaf_size   ;  /* maximum number of nodes in a leaf */
 	double max_leaf_dim;   /* maximum dimension of a leaf node */

 	mol_atom *atoms;       /* pointer to the list of atoms */
 	int natoms;            /* number of atoms in the list */

	int num_nodes;         /* number of nodes in the nodes array */
	int free_node_ptr;     /* index of the first free node in the nodes array, -1 if array if full */
} OCTREE;


typedef struct OPAR
{
   struct agsetup *ags;
   double eps;

   OCTREE *octree_static;
   int node_static;
   OCTREE *octree_moving;
   int node_moving;
   double dist_cutoff, approx_cutoff;
   double hdist_cutoff;
   int fixed_cull;
   double *trans;
   double *engcat;
   void *proc_func_params;
   void ( * processing_function )( struct OPAR *, double * );
} OCTREE_PARAMS;


typedef struct
{
   int node_static, static_cid;
   int node_moving, moving_cid;
} REC_PARAMS;


inline void transform_point( double x, double y, double z, double *trans_mat, double *nx, double *ny, double *nz )
{
   *nx = trans_mat[  0 ] * x + trans_mat[  1 ] * y + trans_mat[  2 ] * z + trans_mat[  3 ];
   *ny = trans_mat[  4 ] * x + trans_mat[  5 ] * y + trans_mat[  6 ] * z + trans_mat[  7 ];
   *nz = trans_mat[  8 ] * x + trans_mat[  9 ] * y + trans_mat[ 10 ] * z + trans_mat[ 11 ];
}



inline double min_pt2bx_dist2( double bx, double by, double bz, double dim, double x, double y, double z )
{
   double dx, dy, dz;
   double dl, dr;
   
   dl = bx - x; 
   dr = - dim - dl;
   dx = ( dl > 0 ) ? dl : ( ( dr > 0 ) ? dr : 0 );   

   dl = by - y; 
   dr = - dim - dl;
   dy = ( dl > 0 ) ? dl : ( ( dr > 0 ) ? dr : 0 );      
   
   dl = bz - z; 
   dr = - dim - dl;
   dz = ( dl > 0 ) ? dl : ( ( dr > 0 ) ? dr : 0 );   
      
   return ( dx * dx + dy * dy + dz * dz );
}


inline int within_distance_cutoff( double x1, double y1, double z1, double dim1, double x2, double y2, double z2, double dim2, double ext )
{
   double d1 = dim1 + ext, d2 = - ( dim2 + ext );

   if ( ( x2 - x1 >= d1 ) || ( x2 - x1 <= d2 )
     || ( y2 - y1 >= d1 ) || ( y2 - y1 <= d2 )
     || ( z2 - z1 >= d1 ) || ( z2 - z1 <= d2 ) ) return 0;   
   
   d1 = 0.5 * ( dim2 - dim1 );
  
   d2 = ( x2 - x1 + d1 ) * ( x2 - x1 + d1 ) 
      + ( y2 - y1 + d1 ) * ( y2 - y1 + d1 )
      + ( z2 - z1 + d1 ) * ( z2 - z1 + d1 );            

   d1 = ( HALF_SQRT_THREE * ( dim1 + dim2 ) + ext );
   d1 *= d1;
   
   return ( d2 < d1 );   
}


inline int inside_node( OCTREE_NODE *node, mol_atom *atom )
{
   return ( ( atom->X - node->lx >= 0 ) && ( atom->X - node->lx < node->dim )
         && ( atom->Y - node->ly >= 0 ) && ( atom->Y - node->ly < node->dim ) 
         && ( atom->Z - node->lz >= 0 ) && ( atom->Z - node->lz < node->dim ) );
}




int build_octree( OCTREE *octree, int max_leaf_size, double max_leaf_dim, double slack_factor, mol_atom_group *ag );

int build_octree_excluding_fixed_atoms( OCTREE *octree, int max_leaf_size, double max_leaf_dim, double slack_factor, mol_atom_group *ag );

void print_octree( OCTREE *octree );

int get_octree_size( OCTREE *octree );

void update_octree( OCTREE *octree, mol_atom* atom );

int reorganize_octree( OCTREE *octree, int batch_update );

double octree_accumulation_excluding_far( OCTREE *octree_static, OCTREE *octree_moving,
                                          double dist_cutoff, double approx_cutoff, int fixed_cull, double *trans,
                                          void *proc_func_params,
                                          void ( * processing_function )( OCTREE_PARAMS *, double * ) );

void destroy_octree( OCTREE *octree );

};
#endif
