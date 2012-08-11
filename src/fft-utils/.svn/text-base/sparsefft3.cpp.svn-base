
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

#include <string.h>
#include "sparsefft3.h"

/***************************************************************************/

/* generic FFT: */

static void performFFT( sparse3DFFT_plan p, int dim1, int dim2, int howmany, FFTW_complex *data )
{
   int dist = 0;
   
   while ( howmany > 0 )
     {
      int k = 0, d = 1;
      
      while ( ( k < MAX_LOG_HOWMANY ) && ( ( d << 1 ) <= howmany ) )
        {
         k++; d <<= 1;
        }
        
      FFTW_execute_dft( p->p[ dim1 ][ dim2 ][ k ], data + dist, data + dist );
        
      howmany -= d;  
      dist += d * ( p->nafter[ dim2 ] );
     }
}


/***************************************************************************/
/* The following are various routines to transform one or two dimensions
   of a 3d array at a time, taking advantage of the sparsity pattern that
   was passed to sparsefft3_create_plan.   For some of them (the ones
   that do two FFTs at a time) we also supply an inverse routine. */

/* Given the 3d array data conforming to the sparsity pattern,
   transform first_dim and then second_dim, arranged in planes
   along plane_dim. */
static void fft2_planes( sparse3DFFT_plan p, int plane_dim, int first_dim, int second_dim, FFTW_complex *data )
{
  for ( int iplane = 0; iplane < p->n[ plane_dim ]; ++iplane ) 
    {
      sparse3DFFT_range r = p->range1[ plane_dim ][ second_dim ][ iplane ];
      
      while ( r ) 
        {
	  performFFT( p, first_dim, second_dim, r->max - r->min + 1, 
	               data + iplane * p->nafter[ plane_dim ] + r->min * p->nafter[ second_dim ] );	  	  	  	    	  
	  r = r->next;
	}
	  
      if ( p->range1[ plane_dim ][ second_dim ][ iplane ] )
	  performFFT( p, second_dim, first_dim, p->n[ first_dim ], data + iplane * p->nafter[ plane_dim ] );	  	  	  	    
     }
}

/* inverse of fft2_planes. */
static void ifft2_planes( sparse3DFFT_plan p, int plane_dim, int first_dim, int second_dim, FFTW_complex *data )
{
  for ( int iplane = 0; iplane < p->n[ plane_dim ]; ++iplane ) 
    {
      sparse3DFFT_range r = p->range1[ plane_dim ][ second_dim ][ iplane ];

      if ( r ) performFFT( p, second_dim, first_dim, p->n[ first_dim ], data + iplane * p->nafter[ plane_dim ] );

      while ( r ) 
        {
	  performFFT( p, first_dim, second_dim, r->max - r->min + 1, 
	               data + iplane * p->nafter[ plane_dim ] + r->min * p->nafter[ second_dim ] );	  	  	  
          r = r->next;
        }
    }
}

/* As fft2_planes, except that we only transform first_dim and leave
   second_dim untransformed. */
static void fft1_planes( sparse3DFFT_plan p, int plane_dim, int first_dim, int second_dim, FFTW_complex *data )
{
  for ( int iplane = 0; iplane < p->n[ plane_dim ]; ++iplane ) 
    {
      sparse3DFFT_range r = p->range1[ plane_dim ][ second_dim ][ iplane ];
      
      while ( r ) 
        {
	  performFFT( p, first_dim, second_dim, r->max - r->min + 1, 
	               data + iplane * p->nafter[ plane_dim ] + r->min * p->nafter[ second_dim ] );	  	    
	  r = r->next;
	}
    }
}

/* Given the 3d array data conforming to the sparsity pattern, but
   after the first dim (plane_dim) has been transformed, transform
   second_dim then third_dim, arranged in planes along plane_dim. */
static void fft2_planes2( sparse3DFFT_plan p, int plane_dim, int second_dim, int third_dim, FFTW_complex *data )
{
  for ( int iplane = 0; iplane < p->n[ plane_dim ]; ++iplane ) 
    {
      sparse3DFFT_range r = p->range2[ third_dim ];
      
      while ( r ) 
        {
	  performFFT( p, second_dim, third_dim, r->max - r->min + 1, 
	              data + iplane * p->nafter[ plane_dim ] + r->min * p->nafter[ third_dim ] );	  	  
          r = r->next;
        }

      performFFT( p, third_dim, second_dim, p->n[ second_dim ], data + iplane * p->nafter[ plane_dim ] );	  	  
    }
}

/* inverse of fft2_planes2 */
static void ifft2_planes2( sparse3DFFT_plan p, int plane_dim, int second_dim, int third_dim, FFTW_complex *data )
{
  for ( int iplane = 0; iplane < p->n[ plane_dim ]; ++iplane ) 
    {
      sparse3DFFT_range r = p->range2[ third_dim ];

      performFFT( p, third_dim, second_dim, p->n[ second_dim ], data + iplane * p->nafter[ plane_dim ] );	  	  

      while ( r ) 
        {
	  performFFT( p, second_dim, third_dim, r->max - r->min + 1, 
	               data + iplane * p->nafter[ plane_dim ] + r->min * p->nafter[ third_dim ] );	  
          r = r->next;
        }
    }
}

/* Given the 3d array data, completely transform the dimension which_dim,
   assuming no sparsity (i.e. all sparsity has been destroyed already). */
static void fft1_all( sparse3DFFT_plan p, int which_dim, FFTW_complex *data )
{
  /* If we are doing the first or the last dimension, we can
	process the whole thing in one do_fft call. */
  if ( which_dim == 0 ) performFFT( p, 0, 2, p->nafter[ 0 ], data );     
  else if ( which_dim == 2 ) performFFT( p, 2, 1, p->n[ 0 ] * p->n[ 1 ], data );         
       else
          {
	    int first_dim, second_dim;
	    
	    if ( ( which_dim + 1 ) % 3 < ( which_dim + 2 ) % 3 ) 
	      {
	        first_dim = ( which_dim + 1 ) % 3;
	        second_dim = ( which_dim + 2 ) % 3;
	      }
	    else 
	      {
	        first_dim = ( which_dim + 2 ) % 3;
	        second_dim = ( which_dim + 1 ) % 3;
	      }
	      
	    for ( int ifirst = 0; ifirst < p->n[ first_dim ]; ++ifirst )
	       performFFT( p, which_dim, second_dim, p->n[ second_dim ], data + ifirst * p->nafter[ first_dim ] );
          }
}

/***************************************************************************/

void sparse3DFFT( sparse3DFFT_plan p, FFTW_complex *data_in, FFTW_complex *data_out )
{
  if ( ( data_in == NULL ) || ( data_out == NULL ) ) return;
  
  if ( data_out != data_in ) memcpy( data_out, data_in, ( p->n[ 0 ] ) * ( p->n[ 1 ] ) * ( p->n[ 2 ] ) * sizeof( FFTW_complex ) );
  
  if ( p->sparsedir == SPARSE3DFFT_SPARSEINPUT ) 
    {    
      if ( p->fft2_first ) 
        {       
	  fft2_planes( p, p->third_dim, p->first_dim, p->second_dim, data_out );
	  fft1_all( p, p->third_dim, data_out );
	}
      else 
        {
	  fft1_planes( p, p->third_dim, p->first_dim, p->second_dim, data_out );
	  fft2_planes2( p, p->first_dim, p->second_dim, p->third_dim, data_out );
	}
    }
  else 
    {
      if ( p->fft2_first ) 
        {
	  fft1_all( p, p->third_dim, data_out );
	  ifft2_planes( p, p->third_dim, p->first_dim, p->second_dim, data_out );
	}
      else 
        {
	  ifft2_planes2( p, p->first_dim, p->second_dim, p->third_dim, data_out );
	  fft1_planes( p, p->third_dim, p->first_dim, p->second_dim, data_out );
	}
    }
}
