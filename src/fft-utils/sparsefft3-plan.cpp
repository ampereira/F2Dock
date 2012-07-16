/* Copyright (C) 2000 Massachusetts Institute of Technology.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
   Modified at CVC Lab, UT Austin for use in F2Dock with FFTW3 

	Authors: Rezaul Alam Chowdhury <shaikat@cs.utexas.edu>
        Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of F2Dock.
*/

#include <stdlib.h>
#include <stdio.h>

#include "sparsefft3-timers.h"
#include "sparsefft3.h"

#define MIN2( a, b ) ( ( a ) < ( b ) ? ( a ) : ( b ) )
#define MAX2( a, b ) ( ( a ) > ( b ) ? ( a ) : ( b ) )
#define MAX3( a, b, c ) MAX2( a, MAX2( b, c ) )

/***************************************************************************/

/* Update r, which points to a linked list of [min,max]
   ranges/intervals, to include the new integers value x (but not any
   other new integers).  Creates a new interval [x,x] if necessary, or
   extends an old one that was adjacent to x.  Upon entry, r may not
   contain any overlapping or adjacent intervals (e.g. [1,2] [3,4]),
   and this property is preserved upon exit. */
static void update_range( sparse3DFFT_range *r, int x )
{
  sparse3DFFT_range last_r = NULL; /* the last element of r's linked list */
  sparse3DFFT_range cur_r = *r;

  while ( cur_r ) 
    {
      if ( ( x >= cur_r->min ) && ( x <= cur_r->max ) ) return;  /* x is already in an element of r */
      else if ( x == cur_r->min - 1 ) 
             {
	       cur_r->min--;
	       break;
	     }
	   else if ( x == cur_r->max + 1 ) 
	          {
	            cur_r->max++;
	            break;
	          }  
	          
      last_r = cur_r;
      cur_r = cur_r->next;
    }
     
  if ( cur_r ) 
    { 
      /* we extended an existing interval */
      sparse3DFFT_range cur2_r = cur_r->next, prev_r = cur_r;

      /* Need to check whether, by extending the existing interval,
         we have caused this interval to be adjacent with some other
	 interval.  In that case, we can merge the two intervals.
	 Note that we don't have to search intervals prior to
	 cur_r, since we already examined whether x was in them
	 or adjacent to them. */

      while ( cur2_r ) 
        {
	  if ( ( cur2_r->min == cur_r->max + 1 ) || ( cur2_r->max == cur_r->min - 1 ) ) 
	    {
	      /* merge cur2_r into cur_r: */
	      cur_r->max = MAX2( cur_r->max, cur2_r->max );
	      cur_r->min = MIN2( cur_r->min, cur2_r->min );

	      /* Delete cur2_r: */
	      prev_r->next = cur2_r->next;
	      FFTW_free( cur2_r );
	      cur2_r = prev_r->next;
	    }
	  else 
	    {
	      prev_r = cur2_r;
	      cur2_r = cur2_r->next;
	    }
	}
    }
  else 
    {  
      /* need a new interval */
      sparse3DFFT_range new_r = ( sparse3DFFT_range) FFTW_malloc( sizeof( sparse3DFFT_range_struct ) );
      
      new_r->min = new_r->max = x;
      new_r->next = NULL;
	
      if ( last_r ) last_r->next = new_r;
      else *r = new_r;
     }     
}

/***************************************************************************/

#ifdef USING_GETTIMEOFDAY

sparse3DFFT_time sparse3DFFT_gettimeofday_get_time( void )
{
  struct timeval tv;

  gettimeofday( &tv, 0 );

  return tv;
}

sparse3DFFT_time sparse3DFFT_gettimeofday_time_diff( sparse3DFFT_time t1, sparse3DFFT_time t2 )
{
  sparse3DFFT_time diff;

  diff.tv_sec = t1.tv_sec - t2.tv_sec;
  diff.tv_usec = t1.tv_usec - t2.tv_usec;
  
  /* normalize */
  while ( diff.tv_usec < 0 ) 
    {
      diff.tv_usec += 1000000L;
      diff.tv_sec -= 1;
    }

  return diff;
}

#endif

double time_sparse3DFFT( sparse3DFFT_plan p, FFTW_complex *data_in, FFTW_complex *data_out )
{
  double t;
  
  TIME_FFT( sparse3DFFT( p, data_in, data_out ), t, SPARSE3DFFT_TIME_MIN );
  
  return t;
}

/***************************************************************************/

sparse3DFFT_plan sparse3DFFT_create_plan( int nx, int ny, int nz,
				          int dir, int flags,
				          sparse3DFFT_sparsedir sparsedir,
				          sparse3DFFT_nonzero_func nonzero,
				          void *nonzero_data,
				          FFTW_complex *data_in,
				          FFTW_complex *data_out )
{
  if ( ( nx < 1 ) || ( ny < 1 ) || ( nz < 1 ) || ( data_in == NULL ) || ( data_out == NULL ) ) return NULL;

  sparse3DFFT_plan p = ( sparse3DFFT_plan ) FFTW_malloc( sizeof( sparse3DFFT_plan_struct ) );

  p->sparsedir = sparsedir;
  p->n[ 0 ] = nx;
  p->n[ 1 ] = ny;
  p->n[ 2 ] = nz;

  for ( int i = 0; i < 3; ++i ) 
    {
      p->nafter[i] = 1;

      for ( int j = 0; j < 3; ++j ) 
        {
	  if ( j > i ) p->nafter[ i ] *= p->n[ j ];
	  
	  if ( j != i ) 
	    {
	      p->range1[ i ][ j ] = ( sparse3DFFT_range * ) FFTW_malloc( sizeof( sparse3DFFT_range ) * p->n[ i ] );
	      
	      for ( int x = 0; x < p->n[ i ]; ++x )
   	         p->range1[ i ][ j ][ x ] = NULL;
	    }
	  else p->range1[ i ][ j ] = NULL;
	}
	
      p->range2[i] = NULL;
    }

  for ( int i = 0; i < 3; i++ ) 
    for ( int j = 0; j < 3; j++ )
      {           
  	for ( int k = 0, d = 1; k <= MAX_LOG_HOWMANY; k++, d *= 2 )
  	  {  	       
  	    if ( ( d - 1 ) * p->nafter[ j ] + ( p->n[ i ] - 1 ) * p->nafter[ i ] >= p->n[ 0 ] * p->n[ 1 ] * p->n[ 2 ] ) 
  	       p->p[ i ][ j ][ k ] = NULL;
  	    else
               p->p[ i ][ j ][ k ] = FFTW_plan_many_dft( 1, &( p->n[ i ] ), d,
                                                         data_out, NULL, p->nafter[ i ], p->nafter[ j ],
                                                         data_out, NULL, p->nafter[ i ], p->nafter[ j ], 
                                                         dir, flags );
  	  }                                            
      }                                                 

  for ( int ix = 0; ix < nx; ++ix )
    for ( int iy = 0; iy < ny; ++iy )
      for ( int iz = 0; iz < nz; ++iz ) 
        {
	  int x[ 3 ];
	  
	  x[ 0 ] = ix; x[ 1 ] = iy; x[ 2 ] = iz;
	  
	  if ( nonzero( x, nonzero_data ) ) 
	    {
	      for ( int i = 0; i < 3; ++i ) 
	        {
		  update_range( &p->range2[ i ], x[ i ] );
		  
		  for ( int j = 0; j < 3; ++j )
  		    if ( j != i ) update_range( &p->range1[ i ][ j ][ x[ i ] ], x[ j ] );
		}
	    }
        }

  /* Determine the fastest way in which to do the FFT,
     preferably by trying and timing all of the possibilities.
     Or, just guess. */
  if ( ( flags & OPT_FFTW_SEARCH_TYPE ) != FFTW_ESTIMATE ) 
    {
      int best_fft2_first = 0, best_first_dim = 0, best_second_dim = 1, best_third_dim = 2;
      double best_time = 1e20;

      for ( int dim2 = 0; dim2 < nx * ny * nz; ++dim2 )
	c_re( data_in[ dim2 ] ) = c_im( data_in[ dim2 ] ) = c_re( data_out[ dim2 ] ) = c_im( data_out[ dim2 ] ) = 0.0;
	
      for ( p->fft2_first = 0; p->fft2_first <= 1; ++p->fft2_first )
        for ( p->first_dim = 0; p->first_dim < 3; ++p->first_dim )
          for ( int dim2 = 1; dim2 <= 2; ++dim2 ) 
            {
	      p->second_dim = ( p->first_dim + dim2 ) % 3;
	      p->third_dim = ( p->first_dim + 3 - dim2 ) % 3;
	     
       	      double t = time_sparse3DFFT( p, data_in, data_out );
       	     
	      if ( t < best_time ) 
	        {
		  best_time = t;
		  best_fft2_first = p->fft2_first;
		  best_first_dim = p->first_dim;
		  best_second_dim = p->second_dim;
		  best_third_dim = p->third_dim;
	        }
	    }
	    
      p->fft2_first = best_fft2_first;
      p->first_dim = best_first_dim;
      p->second_dim = best_second_dim;
      p->third_dim = best_third_dim;
    }
  else 
    {
      p->fft2_first = 0;
      p->first_dim = 0;
      p->second_dim = 1;
      p->third_dim = 2;
    }

  return p;
}

/***************************************************************************/

static void destroy_range( sparse3DFFT_range r )
{
  while ( r ) 
    {
      sparse3DFFT_range cur_r = r;
      r = r->next;
      FFTW_free(cur_r);
    }
}

void sparse3DFFT_destroy_plan( sparse3DFFT_plan p )
{
  if ( p ) 
    {
      for ( int i = 0; i < 3; ++i ) 
        {
	  for ( int j = 0; j < 3; ++j )
            for ( int k = 0; k <= MAX_LOG_HOWMANY; k++ )
              if ( p->p[ i ][ j ][ k ] != NULL ) FFTW_destroy_plan( p->p[ i ][ j ][ k ] );  
	    
	  destroy_range( p->range2[ i ] );
	  
	  for ( int j = 0; j < 3; ++j ) 
	    if ( p->range1[ i ][ j ] ) 
	      {
		for ( int x = 0; x < p->n[ i ]; ++x )
		   destroy_range( p->range1[ i ][ j ][ x ] );
		     
		FFTW_free( p->range1[ i ][ j ] );
	      }
	}
	
      FFTW_free( p );
    }
}

/***************************************************************************/
