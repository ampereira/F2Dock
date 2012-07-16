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
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <libmol/libmol.h>

namespace LibMol{

//#define ALIGNMENT_BLOCK_SIZE 8

int init_free_node_server( OCTREE *octree )
{
   octree->num_nodes = 0;

   octree->free_node_ptr = -1;

   octree->nodes = ( OCTREE_NODE * ) _mol_malloc( INIT_NUM_OCTREE_NODES * sizeof( OCTREE_NODE ) );
//   octree->nodes = ( OCTREE_NODE * ) memalign( ALIGNMENT_BLOCK_SIZE, INIT_NUM_OCTREE_NODES * sizeof( OCTREE_NODE ) );   

   if ( octree->nodes == NULL ) return 0;

   octree->num_nodes = INIT_NUM_OCTREE_NODES;

   for ( int i = octree->num_nodes - 1; i >= 0; i-- )
     {
       octree->nodes[ i ].p_ptr = octree->free_node_ptr;
       octree->free_node_ptr = i;
       
       octree->nodes[ i ].indices = NULL;
      
       octree->nodes[ i ].n = 0;
       octree->nodes[ i ].nfixed = 0;       
       octree->nodes[ i ].id_num = 0;
       octree->nodes[ i ].id_cap = 0;       
     }

   return 1;
}


int next_free_node( OCTREE *octree )
{
   if ( octree->free_node_ptr == -1 )
     {
       int new_num_nodes = 2 * octree->num_nodes;

       if ( new_num_nodes <= 0 ) new_num_nodes = INIT_NUM_OCTREE_NODES;

       octree->nodes = ( OCTREE_NODE * ) _mol_realloc( octree->nodes, new_num_nodes * sizeof( OCTREE_NODE ) );

       if ( octree->nodes == NULL ) return -1;              

       for ( int i = new_num_nodes - 1; i >= octree->num_nodes; i-- )
         {
           octree->nodes[ i ].p_ptr = octree->free_node_ptr;
           octree->free_node_ptr = i;

           octree->nodes[ i ].indices = NULL;
          
           octree->nodes[ i ].n = 0;
           octree->nodes[ i ].nfixed = 0;           
           octree->nodes[ i ].id_num = 0;
           octree->nodes[ i ].id_cap = 0;       
         }
         
        octree->num_nodes = new_num_nodes; 
     }

   int next_node = octree->free_node_ptr;

   octree->free_node_ptr = octree->nodes[ next_node ].p_ptr;

   return next_node;
}


void free_node( OCTREE *octree, int node_id )
{
   freeMem( octree->nodes[ node_id ].indices );
   octree->nodes[ node_id ].indices = NULL;

   octree->nodes[ node_id ].n = 0;
   octree->nodes[ node_id ].nfixed = 0;   
   octree->nodes[ node_id ].id_num = 0;
   octree->nodes[ node_id ].id_cap = 0;

   octree->nodes[ node_id ].p_ptr = octree->free_node_ptr;
   octree->free_node_ptr = node_id;
}



inline void compute_root_bounding_box( int node_id, OCTREE *octree, double slack_factor,
                                       int *indices, int start_id, int end_id )
{
   OCTREE_NODE *node = &( octree->nodes[ node_id ] );
   mol_atom *atoms = octree->atoms;
         
   int s = indices[ start_id ];

   double minX = atoms[ s ].X, minY = atoms[ s ].Y, minZ = atoms[ s ].Z;
   double maxX = atoms[ s ].X, maxY = atoms[ s ].Y, maxZ = atoms[ s ].Z;

   for ( int i = start_id + 1; i <= end_id; i++ )
     {
      int j = indices[ i ];

      if ( atoms[ j ].X < minX ) minX = atoms[ j ].X;
      if ( atoms[ j ].X > maxX ) maxX = atoms[ j ].X;

      if ( atoms[ j ].Y < minY ) minY = atoms[ j ].Y;
      if ( atoms[ j ].Y > maxY ) maxY = atoms[ j ].Y;

      if ( atoms[ j ].Z < minZ ) minZ = atoms[ j ].Z;
      if ( atoms[ j ].Z > maxZ ) maxZ = atoms[ j ].Z;
     }

   double cx = ( minX + maxX ) / 2,
   	 cy = ( minY + maxY ) / 2,
   	 cz = ( minZ + maxZ ) / 2;

   double dim = max( maxX - minX, maxY - minY );

   dim = max( dim, maxZ - minZ );
   dim *= slack_factor;
   
   node->lx = cx - dim * 0.5;   
   node->ly = cy - dim * 0.5;
   node->lz = cz - dim * 0.5;      
   
   node->dim = dim;
}



inline void compute_non_root_bounding_box( int node_id, OCTREE *octree, int child_id )
{
   OCTREE_NODE *node = &( octree->nodes[ node_id ] );
   
   if ( node->p_ptr < 0 ) return;
   
   OCTREE_NODE *pnode = &( octree->nodes[ node->p_ptr ] );
   
   if ( child_id < 0 )
     {
       for ( int i = 0; i < 8; i++ )
         if ( pnode->c_ptr[ i ] == node_id )
           {
             child_id = i;
             break;
           }
           
       if ( child_id < 0 ) return;   
     }
     
   double lx = pnode->lx,
   	 ly = pnode->ly,
   	 lz = pnode->lz;
   double dim = pnode->dim;
   
   dim *= 0.5;
   
   if ( child_id & 1 ) lx += dim;	       
   if ( child_id & 2 ) ly += dim;
   if ( child_id & 4 ) lz += dim;      
         
   node->lx = lx;
   node->ly = ly;
   node->lz = lz;
   
   node->dim = dim;      
}



inline void compute_non_leaf_attributes( int node_id, OCTREE *octree )
{
#ifdef ADD_ATTR         
   OCTREE_NODE *node = &( octree->nodes[ node_id ] );

   double sumX = 0, sumY = 0, sumZ = 0;	
   double sumQ = 0;

   for ( int i = 0; i < 8; i++ )
     if ( node->c_ptr[ i ] >= 0 )
       {
         int j = node->c_ptr[ i ];
         OCTREE_NODE *cnode = &( octree->nodes[ j ] );

         sumX += cnode->sx;
         sumY += cnode->sy;
         sumZ += cnode->sz;
         
         sumQ += cnode->sq;
       }        

   node->sx = sumX;
   node->sy = sumY;
   node->sz = sumZ;      
   
   node->sq = sumQ;
#endif   
}


inline void compute_leaf_attributes( int node_id, OCTREE *octree,
                                     int *indices, int start_id, int end_id )
{
#ifdef ADD_ATTR
   OCTREE_NODE *node = &( octree->nodes[ node_id ] );        
   mol_atom *atoms = octree->atoms;

   double sumX = 0, sumY = 0, sumZ = 0, sumQ = 0;
   
   for ( int i = start_id; i <= end_id; i++ )
     {
      int j = indices[ i ];

      sumX += atoms[ j ].X;
      sumY += atoms[ j ].Y;
      sumZ += atoms[ j ].Z;
      
      sumQ += atoms[ j ].chrg;                  
     }

   node->sx = sumX;   
   node->sy = sumY;
   node->sz = sumZ;      
   
   node->sq = sumQ;   
#endif   
}


inline int get_child_id( OCTREE_NODE *node, mol_atom *atom )
{
   double dim = 0.5 * node->dim;
   double cx = node->lx + dim, cy = node->ly + dim, cz = node->lz + dim;
 
   int k = ( zeroIfLess( atom->Z, cz ) << 2 )
         + ( zeroIfLess( atom->Y, cy ) << 1 )
         + ( zeroIfLess( atom->X, cx ) );

   return k;
}



int expand_octree_node( int node_id, OCTREE *octree,
                        int *indices, int *indices_temp, int start_id, int end_id )
{
   OCTREE_NODE *node = &( octree->nodes[ node_id ] );
   mol_atom *atoms = octree->atoms;

   int n = node->n = end_id - start_id + 1;
   double dim = node->dim;

   if ( ( n <= octree->max_leaf_size ) || ( dim <= octree->max_leaf_dim ) )
     {     
      node->leaf = 1;

      node->indices = ( int * ) _mol_malloc( 2 * n * sizeof ( int ) );
//      node->indices = ( int * ) memalign( ALIGNMENT_BLOCK_SIZE, 2 * n * sizeof ( int ) );

      if ( node->indices == NULL )
        {
          print_error( "Failed to allocate leaf node memory for octree!" );
          return 0;
        }

      node->id_cap = 2 * n;
      
      int nfixed = 0;
      
      for ( int i = start_id; i <= end_id; i++ )
        {
          int j = indices[ i ];
          
          if ( atoms[ j ].fixed ) 
            {           
              node->indices[ nfixed ] = j;
              atoms[ j ].octree_ptr = create_octree_ptr( node_id, nfixed );
              nfixed++;
            }  
        }
                
      node->nfixed = nfixed;        

      if ( nfixed < n ) 
         {
           for ( int i = start_id, k = nfixed; i <= end_id; i++ )
             {
               int j = indices[ i ];
               
               if ( atoms[ j ].fixed ) continue;

               node->indices[ k ] = j;
               atoms[ j ].octree_ptr = create_octree_ptr( node_id, k );
               k++;
             }
         }                     
         
#ifdef ADD_ATTR        
       compute_leaf_attributes( node_id, octree, indices, start_id, end_id );        
#endif       
     }
   else
     {
      node->leaf = 0;

      int count[ 8 ] = { 0, 0, 0, 0, 0, 0, 0, 0 };

      for ( int i = start_id; i <= end_id; i++ )
        {
         int j = indices[ i ];
         int k = get_child_id( node, &( atoms[ j ] ) );
         count[ k ]++;
        }

      int start_index[ 8 ];
      int cur_index[ 8 ];

      cur_index[ 0 ] = start_index[ 0 ] = start_id;
      for ( int i = 1; i < 8; i++ )
        cur_index[ i ] = start_index[ i ] = start_index[ i - 1 ] + count[ i - 1 ];

      for ( int i = start_id; i <= end_id; i++ )
        {
         int j = indices[ i ];
         int k = get_child_id( node, &( atoms[ j ] ) );

         indices_temp[ cur_index[ k ] ] = j;
         cur_index[ k ]++;
        }

      node->nfixed = 0;
            
      for ( int i = 0; i < 8; i++ )
       { 
        if ( count[ i ] > 0 )
          {
           int j = next_free_node( octree );
           
           node = &( octree->nodes[ node_id ] );
           
           node->c_ptr[ i ] = j;
           
           octree->nodes[ j ].p_ptr = node_id;                      

           compute_non_root_bounding_box( j, octree, i );           
           
           if ( !expand_octree_node( j, octree, indices_temp, indices, start_index[ i ], start_index[ i ] + count[ i ] - 1 ) ) return 0;

           node = &( octree->nodes[ node_id ] );           
           
           node->nfixed += octree->nodes[ j ].nfixed;
          }
        else node->c_ptr[ i ] = -1;        
       } 
      
#ifdef ADD_ATTR        
       compute_non_leaf_attributes( node_id, octree );        
#endif       
     }
     
   return 1;
}


void collect_atoms_from_leaves( int node_id, OCTREE *octree, int *indices, int start_id )
{
   OCTREE_NODE *node = &( octree->nodes[ node_id ] );
   
   if ( node->leaf )
     {
       for ( int i = 0; i < node->n; i++ )
          indices[ start_id + i ] = node->indices[ i ];
     }   
   else
     {
       for ( int i = 0; i < 8; i++ )
         if ( node->c_ptr[ i ] >= 0 )
           {
             int j = node->c_ptr[ i ];
             
             collect_atoms_from_leaves( j, octree, indices, start_id );
             
             start_id += octree->nodes[ j ].n;
             
             free_node( octree, j );
           }        
     }  
}


int contract_octree_node( int node_id, OCTREE *octree )
{
   OCTREE_NODE *node = &( octree->nodes[ node_id ] );
   mol_atom *atoms = octree->atoms;

   int n = node->n;
   
   if ( ( node->leaf ) || ( n > octree->max_leaf_size ) ) return 1;

   node->indices = ( int * ) _mol_malloc( sizeof ( int ) * 2 * n );
//   node->indices = ( int * ) memalign( ALIGNMENT_BLOCK_SIZE, sizeof ( int ) * 2 * n );

   if ( node->indices == NULL )
     {
       print_error( "Failed to allocate leaf node memory for octree!" );
       return 0;
     }

   collect_atoms_from_leaves( node_id, octree, node->indices, 0 );

   node->leaf = 1;
   node->id_cap = 2 * n;

   for ( int i = 0, k = 0; i < n; i++ )
     {
       int j = node->indices[ i ];
       
       if ( atoms[ j ].fixed )
         {
           int l = node->indices[ k ];
           node->indices[ k++ ] = j;
           node->indices[ i ] = l;
         }
     }
   
   for ( int i = 0; i < n; i++ )
     {
       int j = node->indices[ i ];
       atoms[ j ].octree_ptr = create_octree_ptr( node_id, i );       
     }

   return 1;
}




int build_octree( OCTREE *octree, int max_leaf_size, double max_leaf_dim, double slack_factor, mol_atom_group *ag )
{
   int *indices, *indices_temp;

   indices = ( int * ) _mol_malloc( sizeof ( int ) * ag->natoms );
//   indices = ( int * ) memalign( ALIGNMENT_BLOCK_SIZE, sizeof ( int ) * ag->natoms );
   indices_temp = ( int * ) _mol_malloc( sizeof ( int ) * ag->natoms );
//   indices_temp = ( int * ) memalign( ALIGNMENT_BLOCK_SIZE, sizeof ( int ) * ag->natoms );

   if ( ( indices == NULL ) || ( indices_temp == NULL ) )
     {
      print_error( "Failed to allocate temporary storage for octree!" );
      freeMem( indices );
      freeMem( indices_temp );
      return 0;
     }

   for ( int i = 0; i < ag->natoms; i++ )
     indices[ i ] = i;

   octree->max_leaf_size = max_leaf_size;
   octree->max_leaf_dim = max_leaf_dim;
   octree->atoms = ag->atoms;
   octree->natoms = ag->natoms;   

   init_free_node_server( octree );

   int octree_root = next_free_node( octree );

   compute_root_bounding_box( octree_root, octree, slack_factor, indices, 0, ag->natoms - 1 );
   
   octree->nodes[ octree_root ].p_ptr = -1;

   int built = expand_octree_node( octree_root, octree, indices, indices_temp, 0, ag->natoms - 1 );

   freeMem( indices );
   freeMem( indices_temp );

   return built;
}



int build_octree_excluding_fixed_atoms( OCTREE *octree, int max_leaf_size, double max_leaf_dim, double slack_factor, mol_atom_group *ag )
{
   int *indices, *indices_temp;
   int nflex = 0;
   
   for ( int i = 0; i < ag->natoms; i++ )
     if ( !ag->atoms[ i ].fixed ) nflex++;
     
   printf( "ag->natoms = %d, nflex = %d\n", ag->natoms, nflex );  
     
   indices = ( int * ) _mol_malloc( sizeof ( int ) * nflex );
//   indices = ( int * ) memalign( ALIGNMENT_BLOCK_SIZE, sizeof ( int ) * nflex );
   indices_temp = ( int * ) _mol_malloc( sizeof ( int ) * nflex );
//   indices_temp = ( int * ) memalign( ALIGNMENT_BLOCK_SIZE, sizeof ( int ) * nflex );

   if ( ( indices == NULL ) || ( indices_temp == NULL ) )
     {
      print_error( "Failed to allocate temporary storage for octree!" );
      freeMem( indices );
      freeMem( indices_temp );
      return 0;
     }

   for ( int i = 0, k = 0; i < ag->natoms; i++ )
      if ( !ag->atoms[ i ].fixed ) indices[ k++ ] = i;

   octree->max_leaf_size = max_leaf_size;
   octree->max_leaf_dim = max_leaf_dim;
   octree->atoms = ag->atoms;

   init_free_node_server( octree );

   int octree_root = next_free_node( octree );

   compute_root_bounding_box( octree_root, octree, slack_factor, indices, 0, nflex - 1 );
   
   octree->nodes[ octree_root ].p_ptr = -1;

   int built = expand_octree_node( octree_root, octree, indices, indices_temp, 0, nflex - 1 );

   freeMem( indices );
   freeMem( indices_temp );

   return built;
}



void traverse_octree( int node_id, OCTREE *octree )
{   
   OCTREE_NODE *node = &( octree->nodes[ node_id ] );
 
   printf( "%d ( %d, %d, %lf ): ", node_id, node->n, node->nfixed, node->dim );
   
   if ( !node->leaf )
      for ( int i = 0; i < 8; i++ )
        if ( node->c_ptr[ i ] >= 0 ) printf( "%d ", node->c_ptr[ i ] );
     
   printf( "\n" );  

   if ( !node->leaf )
      for ( int i = 0; i < 8; i++ )
        if ( node->c_ptr[ i ] >= 0 ) 
           traverse_octree( node->c_ptr[ i ], octree );
}


void print_octree( OCTREE *octree )
{
   traverse_octree( 0, octree );
}


int get_subtree_size( int node_id, OCTREE *octree )
{   
   OCTREE_NODE *node = &( octree->nodes[ node_id ] );
   
   int s = sizeof( OCTREE_NODE );
   
   s += node->id_cap * sizeof( int );
 
   if ( !node->leaf )
      for ( int i = 0; i < 8; i++ )
        if ( node->c_ptr[ i ] >= 0 ) 
           s += get_subtree_size( node->c_ptr[ i ], octree );
           
   return s;        
}


int get_octree_size( OCTREE *octree )
{
   return get_subtree_size( 0, octree );
}



inline int remove_atom_from_leaf( int node_id, OCTREE *octree, int atom_id )
{
   mol_atom *atom = &( octree->atoms[ atom_id ] );

   int j = get_index_in_node( atom->octree_ptr );
   
   OCTREE_NODE *node = &( octree->nodes[ node_id ] );
   
   int n = node->n;
   
   if ( atom->fixed )
     {
       int nf = node->nfixed - 1;
       
       node->indices[ j ] = node->indices[ nf ];   
       node->indices[ nf ] = node->indices[ n - 1 ];
       
       octree->atoms[ node->indices[ nf ] ].octree_ptr = octree->atoms[ node->indices[ j ] ].octree_ptr;                    
       octree->atoms[ node->indices[ j ] ].octree_ptr = atom->octree_ptr;   

       node->nfixed = nf;
     }
   else
     {     
       node->indices[ j ] = node->indices[ n - 1 ];   
       octree->atoms[ node->indices[ j ] ].octree_ptr = atom->octree_ptr;   
     }
       
   node->n = n - 1;
   
   if ( n <= ( node->id_cap >> 2 ) )
     {
       node->indices = ( int * ) _mol_realloc( node->indices, ( node->id_cap >> 1 ) * sizeof( int ) );
       
       if ( node->indices == NULL ) 
         {
           print_error( "Failed to contract leaf storage for octree!" );         
           return 0;
         }   
       
       node->id_cap >>= 1;
     }

#ifdef ADD_ATTR        
   node->sx -= atom->X;
   node->sy -= atom->Y;
   node->sz -= atom->Z;
    
   node->sq -= atom->chrg;            
#endif       

   atom->octree_ptr = -1;

   return 1;
}


inline int remove_atom_from_non_leaf( int node_id, OCTREE *octree, int atom_id )
{
   OCTREE_NODE *node = &( octree->nodes[ node_id ] );
   mol_atom *atom = &( octree->atoms[ atom_id ] );

   int n = node->n;   
   node->n = n - 1;
   if ( atom->fixed ) node->nfixed--;

#ifdef ADD_ATTR        
   node->sx -= atom->X;
   node->sy -= atom->Y;
   node->sz -= atom->Z;
    
   node->sq -= atom->chrg;            
#endif       

   if ( node->n < ( octree->max_leaf_size >> 1 ) ) return contract_octree_node( node_id, octree );
   else return 1;
}



inline void add_atom_to_non_leaf( int node_id, OCTREE *octree, int atom_id )
{
   OCTREE_NODE *node = &( octree->nodes[ node_id ] );
   mol_atom *atom = &( octree->atoms[ atom_id ] );

   int n = node->n;   
   node->n = n + 1;
   if ( atom->fixed ) node->nfixed++;

#ifdef ADD_ATTR        
   node->sx += atom->X;
   node->sy += atom->Y;
   node->sz += atom->Z;
    
   node->sq += atom->chrg;
#endif       
}


inline int add_upward_migrating_atom_to_non_leaf( int node_id, OCTREE *octree, int atom_id )
{
   OCTREE_NODE *node = &( octree->nodes[ node_id ] );

   int m = node->id_num;   

   if ( !m )
     {
       node->indices = ( int * ) _mol_malloc( INIT_MIGRATION_ARRAY_SIZE * sizeof( int ) );
//       node->indices = ( int * ) memalign( ALIGNMENT_BLOCK_SIZE, INIT_MIGRATION_ARRAY_SIZE * sizeof( int ) );
       
       if ( node->indices == NULL )
         {
           print_error( "Failed to allocate temporary migration array for octree!" );
	   return 0;
         }                  
         
       node->id_cap = INIT_MIGRATION_ARRAY_SIZE;  
     }
     
   if ( m == node->id_cap )
     {
       node->indices = ( int * ) _mol_realloc( node->indices, ( node->id_cap << 1 ) * sizeof( int ) );
       
       if ( node->indices == NULL )
         {
           print_error( "Failed to reallocate temporary migration array for octree!" );
	   return 0;
         }                  
         
       node->id_cap <<= 1;  
     }
     
   node->indices[ m ] = atom_id;  
   node->id_num = m + 1;
   
   return 1;  
}


inline int remove_upward_migrating_atom_from_non_leaf( int node_id, OCTREE *octree, int index )
{
   OCTREE_NODE *node = &( octree->nodes[ node_id ] );

   int m = node->id_num - 1;   

   int atom_id = node->indices[ index ];
   mol_atom *atom = &( octree->atoms[ atom_id ] );

   node->indices[ index ] = node->indices[ m ]; 
   node->id_num = m;
   ( node->n )--;
   if ( atom->fixed ) ( node->nfixed )--;     

#ifdef ADD_ATTR        
   node->sx -= atom->X;
   node->sy -= atom->Y;
   node->sz -= atom->Z;
    
   node->sq -= atom->chrg;            
#endif       

   if ( ( m <= ( node->id_cap >> 2 ) ) && ( m >= ( INIT_MIGRATION_ARRAY_SIZE >> 1 ) ) )
     {
       node->indices = ( int * ) _mol_realloc( node->indices, ( m << 1 ) * sizeof( int ) );
       
       if ( node->indices == NULL )
         {
           print_error( "Failed to allocate temporary migration array for octree!" );
	   return 0;
         }                  
         
       node->id_cap = ( m << 1 );  
     }
        
   return 1;  
}


inline int add_downward_migrating_atom_to_non_leaf( int node_id, OCTREE *octree, int atom_id )
{
   OCTREE_NODE *node = &( octree->nodes[ node_id ] );
   mol_atom *atom = &( octree->atoms[ atom_id ] );

   int m = node->id_num;   

   if ( !m )
     {
       node->indices = ( int * ) _mol_malloc( INIT_MIGRATION_ARRAY_SIZE * sizeof( int ) );
//       node->indices = ( int * ) memalign( ALIGNMENT_BLOCK_SIZE, INIT_MIGRATION_ARRAY_SIZE * sizeof( int ) );
       
       if ( node->indices == NULL )
         {
           print_error( "Failed to allocate temporary migration array for octree!" );
	   return 0;
         }                  
         
       node->id_cap = INIT_MIGRATION_ARRAY_SIZE;  
     }
     
   if ( m == node->id_cap )
     {
       node->indices = ( int * ) _mol_realloc( node->indices, ( node->id_cap << 1 ) * sizeof( int ) );
       
       if ( node->indices == NULL )
         {
           print_error( "Failed to reallocate temporary migration array for octree!" );
	   return 0;
         }                  
         
       node->id_cap <<= 1;  
     }
     
   node->indices[ m ] = atom_id;  
   node->id_num = m + 1;
   ( node->n )++;
   if ( atom->fixed ) ( node->nfixed )++;

#ifdef ADD_ATTR        
   node->sx += atom->X;
   node->sy += atom->Y;
   node->sz += atom->Z;
    
   node->sq += atom->chrg;            
#endif       
        
   return 1;  
}


inline int remove_downward_migrating_atom_from_non_leaf( int node_id, OCTREE *octree, int index )
{
   OCTREE_NODE *node = &( octree->nodes[ node_id ] );

   int m = node->id_num - 1;   

   if ( index != m ) node->indices[ index ] = node->indices[ m ]; 

   node->id_num = m;

   if ( ( m <= ( node->id_cap >> 2 ) ) && ( m >= ( INIT_MIGRATION_ARRAY_SIZE >> 1 ) ) )
     {
       node->indices = ( int * ) _mol_realloc( node->indices, ( m << 1 ) * sizeof( int ) );
       
       if ( node->indices == NULL )
         {
           print_error( "Failed to allocate temporary migration array for octree!" );
	   return 0;
         }                  
         
       node->id_cap = ( m << 1 );  
     }
        
   return 1;  
}



inline int add_atom_to_leaf( int node_id, OCTREE *octree, int atom_id )
{
   OCTREE_NODE *node = &( octree->nodes[ node_id ] );
   mol_atom *atom = &( octree->atoms[ atom_id ] );
   
   int n = node->n;

   if ( n == node->id_cap )
     {
       if ( node->id_cap == 0 )
         {
           node->id_cap = 1;
           node->indices = ( int * ) _mol_malloc( ( node->id_cap << 1 ) * sizeof( int ) ); 
//           node->indices = ( int * ) memalign( ALIGNMENT_BLOCK_SIZE, ( node->id_cap << 1 ) * sizeof( int ) ); 
         }  
       else node->indices = ( int * ) _mol_realloc( node->indices, ( node->id_cap << 1 ) * sizeof( int ) );
       
       if ( node->indices == NULL ) 
         {
           print_error( "Failed to expand leaf storage for octree!" );         
           return 0;
         }   
       
       node->id_cap <<= 1;
     }

   if ( atom->fixed )
     {
       int nf = node->nfixed;
       
       if ( n > 0 )
         {
           node->indices[ n ] = node->indices[ nf ];   
           octree->atoms[ node->indices[ n ] ].octree_ptr = create_octree_ptr( node_id, n );
         }

       node->indices[ nf ] = atom_id;   
       atom->octree_ptr = create_octree_ptr( node_id, nf );
       
       node->nfixed = nf + 1;
     } 
   else
     {    
       node->indices[ n ] = atom_id;   
       atom->octree_ptr = create_octree_ptr( node_id, n );
     }
       
   node->n = n + 1;

#ifdef ADD_ATTR        
   node->sx += atom->X;
   node->sy += atom->Y;
   node->sz += atom->Z;
    
   node->sq += atom->chrg;    
#endif       

   if ( node->n > ( octree->max_leaf_size << 1 ) )
     {
       int *temp = ( int * ) _mol_malloc( node->n * sizeof( int ) );
//       int *temp = ( int * ) memalign( ALIGNMENT_BLOCK_SIZE, node->n * sizeof( int ) );
       
       if ( temp == NULL ) 
         {
           print_error( "Failed to allocate temporary storage for octree!" );         
           return 0;
         }   
     
       int done = expand_octree_node( node_id, octree, node->indices, temp, 0, node->n - 1 );
       
       freeMem( temp );
       
       return done;
     }
   else return 1; 
}



int push_down( int node_id, OCTREE *octree, int atom_id )
{
   OCTREE_NODE *node = &( octree->nodes[ node_id ] );
   mol_atom *atom = &( octree->atoms[ atom_id ] );
   
   if ( node->leaf ) return add_atom_to_leaf( node_id, octree, atom_id ); 
   else 
     {
       add_atom_to_non_leaf( node_id, octree, atom_id );             

       for ( int i = 0; i < 8; i++ )
         { 
          if ( node->c_ptr[ i ] >= 0 )
            {
              OCTREE_NODE *cnode = &( octree->nodes[ node->c_ptr[ i ] ] );           
                            
              if ( inside_node( cnode, atom ) ) return push_down( node->c_ptr[ i ], octree, atom_id );
            } 
          else
            {
              double lx = node->lx,
   	            ly = node->ly,
   	            lz = node->lz;
   	             
              double hdim = 0.5 * node->dim;
      
              if ( i & 1 ) lx += hdim;	       
              if ( i & 2 ) ly += hdim;
              if ( i & 4 ) lz += hdim;
              
              if ( !( ( atom->X >= lx ) && ( atom->X < lx + hdim )
                   && ( atom->Y >= ly ) && ( atom->Y < ly + hdim ) 
                   && ( atom->Z >= lz ) && ( atom->Z < lz + hdim ) ) ) continue;
                                
              int j = next_free_node( octree );
               
              node = &( octree->nodes[ node_id ] );
               
              node->c_ptr[ i ] = j;
              
              OCTREE_NODE *cnode = &( octree->nodes[ j ] );
                             
              cnode->p_ptr = node_id;                     
              
              cnode->lx = lx; 
              cnode->ly = ly;
              cnode->lz = lz;                            
              
              cnode->dim = hdim;              
      
              cnode->leaf = 1;
              
              return push_down( j, octree, atom_id );
            }  
         }
         
       return 0;
     }  
}



int pull_up( int node_id, OCTREE *octree, int atom_id )
{
   OCTREE_NODE *node = &( octree->nodes[ node_id ] );
   mol_atom *atom = &( octree->atoms[ atom_id ] );
   
   if ( !inside_node( node, atom ) )
     {
       if ( node->p_ptr < 0 )
         {
           print_error( "Atom has moved outside the root bounding box!" );         
           return 0;
         } 
       
       if ( node->leaf ) remove_atom_from_leaf( node_id, octree, atom_id ); 
       else remove_atom_from_non_leaf( node_id, octree, atom_id );
       
       return pull_up( node->p_ptr, octree, atom_id );
     }
   else return push_down( node_id, octree, atom_id );
}



int batch_pull_up( int node_id, OCTREE *octree, int *empty )
{
   OCTREE_NODE *node = &( octree->nodes[ node_id ] );
   
   if ( node->n == node->nfixed )
     {
       *empty = ( node->n == 0 );
       return 1;
     }

   if ( node->leaf )
     {
       if ( node->p_ptr >= 0 )
         {
           for ( int i = node->nfixed; i < node->n; i++ )
             {
               int atom_id = node->indices[ i ];
               mol_atom *atom = &( octree->atoms[ atom_id ] );
               
               if ( !inside_node( node, atom ) )
                 {
                   if ( !remove_atom_from_leaf( node_id, octree, atom_id ) ) return 0;
                   if ( !add_upward_migrating_atom_to_non_leaf( node->p_ptr, octree, atom_id ) ) return 0;
                   i--;
                 }
             }
         }   
     }
   else
     {
       int m = node->id_num;
       int nc = 0;
       
       for ( int i = 0; i < 8; i++ )
         if ( node->c_ptr[ i ] >= 0 ) 
           {
             int emp;
           
             if ( !batch_pull_up( node->c_ptr[ i ], octree, &emp ) ) return 0;
             
             if ( emp ) node->c_ptr[ i ] = -1;
             else nc++;  
           }  

       if ( node->p_ptr >= 0 )
         {           
           for ( int i = m; i < node->id_num; i++ )
             {
               int atom_id = node->indices[ i ];
               mol_atom *atom = &( octree->atoms[ atom_id ] );
               
               if ( !inside_node( node, atom ) )
                 {
                   if ( !remove_upward_migrating_atom_from_non_leaf( node_id, octree, i ) ) return 0;
                   if ( !add_upward_migrating_atom_to_non_leaf( node->p_ptr, octree, atom_id ) ) return 0;
                   i--;
                 }           
             }
         }
         
       if ( nc == 0 )
         {
           node->leaf = 1;
           node->n = node->id_num;
           node->id_num = 0;
           
           int nf = 0;
           
           for ( int i = 0; i < node->n; i++ )
             {
               int j = node->indices[ i ];
               
               if ( octree->atoms[ j ].fixed )
                 {
                   int l = node->indices[ nf ];
                   node->indices[ nf++ ] = j;
                   node->indices[ i ] = l;
                 }
             }
             
           node->nfixed = nf;             
           
           for ( int i = 0; i < node->n; i++ )
             {
               int j = node->indices[ i ];
               octree->atoms[ j ].octree_ptr = create_octree_ptr( node_id, i );       
             }           
         }               
     }
     
   if ( node->n + node->id_num == 0 )
     {
       free_node( octree, node_id );
       *empty = 1;
     }            
   else *empty = 0;     
       
   return 1;   
}


int batch_push_down( int node_id, OCTREE *octree )
{
   OCTREE_NODE *node = &( octree->nodes[ node_id ] );

   if ( node->leaf ) 
     {
       if ( node->n > ( octree->max_leaf_size << 1 ) )
         {
           int *temp = ( int * ) _mol_malloc( node->n * sizeof( int ) );
//           int *temp = ( int * ) memalign( ALIGNMENT_BLOCK_SIZE, node->n * sizeof( int ) );
           
           if ( temp == NULL ) 
             {
               print_error( "Failed to allocate temporary storage for octree!" );         
               return 0;
             }   
         
           int done = expand_octree_node( node_id, octree, node->indices, temp, 0, node->n - 1 );
           
           freeMem( temp );
           
           return done;
         }
       else return 1;     
     } 

   for ( int i = node->id_num - 1; i >= 0; i-- )
     {
       int atom_id = node->indices[ i ];
       mol_atom *atom = &( octree->atoms[ atom_id ] );
       
       int k = get_child_id( node, atom );
       
       if ( node->c_ptr[ k ] < 0 )
         {
           int j = next_free_node( octree );

           node = &( octree->nodes[ node_id ] );

           octree->nodes[ j ].p_ptr = node_id;           

           compute_non_root_bounding_box( j, octree, k );           

           node->c_ptr[ k ] = j;
           
           octree->nodes[ j ].n = 0;
           octree->nodes[ j ].nfixed = 0;           
           octree->nodes[ j ].leaf = 1;
         }
       
       if ( !remove_downward_migrating_atom_from_non_leaf( node_id, octree, i ) ) return 0;
       
       OCTREE_NODE *cnode = &( octree->nodes[ node->c_ptr[ k ] ] );

       if ( cnode->leaf )
         {
           if ( !add_atom_to_leaf( node->c_ptr[ k ], octree, atom_id ) ) return 0;
         }
       else 
         {
           if ( !add_downward_migrating_atom_to_non_leaf( node->c_ptr[ k ], octree, atom_id ) ) return 0;
         }  
     }
     
   node->id_num = 0;      
     
   for ( int i = 0; i < 8; i++ )
     if ( node->c_ptr[ i ] >= 0 )
       { 
         if ( !batch_push_down( node->c_ptr[ i ], octree ) ) return 0;
         node = &( octree->nodes[ node_id ] );
       }  

   if ( node->n < ( octree->max_leaf_size >> 1 ) ) return contract_octree_node( node_id, octree );

   return 1;   
}



void update_octree( OCTREE *octree, mol_atom* atom )
{
   int node_id = get_node_id( atom->octree_ptr );
   OCTREE_NODE *node = &( octree->nodes[ node_id ] );
   int j = get_index_in_node( atom->octree_ptr );
   int atom_id = node->indices[ j ];

   if ( !inside_node( node, atom ) ) pull_up( node_id, octree, atom_id );
}



int reorganize_octree( OCTREE *octree, int batch_update )
{
   if ( batch_update )
     {
       int emp;
       
       if ( !batch_pull_up( 0, octree, &emp ) ) return 0;
       if ( !batch_push_down( 0, octree ) ) return 0;   
     }
   else
     {
       for ( int i = 0; i < octree->natoms; i++ )
          {
            mol_atom *atom = &( octree->atoms[ i ] );
            if ( !atom->fixed ) update_octree( octree, atom );
          }       
     }  
     
   return 1;
}



void free_subtree_nodes( OCTREE *octree, int node_id )
{
   OCTREE_NODE *node = &( octree->nodes[ node_id ] );
   
   if ( !node->leaf )
     {
       for ( int i = 0; i < 8; i++ )
         if ( node->c_ptr[ i ] >= 0 )
            free_subtree_nodes( octree, node->c_ptr[ i ] );
     }
     
   free_node( octree, node_id );    
}


void destroy_octree( OCTREE *octree )
{
   free_subtree_nodes( octree, 0 );
   freeMem( octree->nodes );
}




double accumulate_excluding_far( OCTREE_PARAMS *octpar )                                 
{
   OCTREE_NODE *static_node = &( octpar->octree_static->nodes[ octpar->node_static ] );
   OCTREE_NODE *moving_node = &( octpar->octree_moving->nodes[ octpar->node_moving ] );

   if ( moving_node->nfixed == moving_node->n ) return 0;     

   if ( octpar->trans == NULL )
     {
       if ( !within_distance_cutoff( static_node->lx, static_node->ly, static_node->lz, static_node->dim,
                                     moving_node->lx, moving_node->ly, moving_node->lz, moving_node->dim,
                                     octpar->approx_cutoff ) ) return 0;     
     }
   else
     {
       double sumRad = HALF_SQRT_THREE * ( static_node->dim + moving_node->dim );
       double maxD2 = ( sumRad + octpar->dist_cutoff );
       maxD2 *= maxD2;
      
       double half_dim = 0.5 * moving_node->dim;
       double mx = moving_node->lx + half_dim,
              my = moving_node->ly + half_dim,
              mz = moving_node->lz + half_dim;
      
       transform_point( mx, my, mz, octpar->trans, &mx, &my, &mz );
      
       half_dim = 0.5 * static_node->dim;
      
       double dx = static_node->lx + half_dim - mx,
              dy = static_node->ly + half_dim - my,
              dz = static_node->lz + half_dim - mz;      
       
       double d2 = dx * dx + dy * dy + dz * dz;            
       
       if ( d2 >= maxD2 ) return 0;                                                     
     }
            
   double energy = 0;

   if ( static_node->leaf )  
     {
       if ( moving_node->leaf ) octpar->processing_function( octpar, &energy );
       else
         {
           for ( int j = 0; j < 8; j++ )
             if ( moving_node->c_ptr[ j ] >= 0 )
               {
                 octpar->node_moving = moving_node->c_ptr[ j ];                           
                 energy += accumulate_excluding_far( octpar );
               }                                                   
         }
     }
   else
     {
       if ( moving_node->leaf )
         {
           for ( int i = 0; i < 8; i++ )
             if ( static_node->c_ptr[ i ] >= 0 )
               {
                 octpar->node_static = static_node->c_ptr[ i ];
                 energy += accumulate_excluding_far( octpar );
               }                                                   
         }
       else
         {
           for ( int i = 0; i < 8; i++ )
             if ( static_node->c_ptr[ i ] >= 0 )
                for ( int j = 0; j < 8; j++ )
                  if ( moving_node->c_ptr[ j ] >= 0 )
                    {
                      octpar->node_static = static_node->c_ptr[ i ];                                            
                      octpar->node_moving = moving_node->c_ptr[ j ];                           
                      energy += accumulate_excluding_far( octpar );
                    }
         }  
     }  
      
   return energy;
}



double octree_accumulation_excluding_far( OCTREE *octree_static, OCTREE *octree_moving,
                                          double dist_cutoff, double approx_cutoff, int fixed_cull, double *trans,
                                          void *proc_func_params,
                                          void ( * processing_function )( OCTREE_PARAMS *, double * ) )
{
   OCTREE_PARAMS octpar;
   
   octpar.octree_static = octree_static;
   octpar.node_static = 0;

   octpar.octree_moving = octree_moving;
   octpar.node_moving = 0;
      
   octpar.dist_cutoff = dist_cutoff;
   octpar.approx_cutoff = approx_cutoff;   
   octpar.fixed_cull = fixed_cull;   
   octpar.trans = trans;
   octpar.proc_func_params = proc_func_params;
   octpar.processing_function = processing_function;
                
   return accumulate_excluding_far( &octpar );
}

};
