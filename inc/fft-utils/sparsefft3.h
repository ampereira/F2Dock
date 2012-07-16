/*
   Modified at CVC Lab, UT Austin for use in F2Dock with FFTW3 
*/

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

#ifndef SPARSE3DFFT_H
#define SPARSE3DFFT_H

#include "fftw3.h"
#include "fftwPrecision.h"

#define c_re( c ) ( ( c )[ 0 ] )
#define c_im( c ) ( ( c )[ 1 ] )

#define MAX_LOG_HOWMANY 7

struct sparse3DFFT_range_struct_ 
  {
    int min, max;
    struct sparse3DFFT_range_struct_ *next;
  };
  
typedef struct sparse3DFFT_range_struct_ sparse3DFFT_range_struct;
typedef sparse3DFFT_range_struct *sparse3DFFT_range;

typedef enum 
  {
    SPARSE3DFFT_SPARSEINPUT,
    SPARSE3DFFT_SPARSEOUTPUT
  } sparse3DFFT_sparsedir;

typedef struct 
  {
    sparse3DFFT_sparsedir sparsedir;
    int n[ 3 ], nafter[ 3 ];
    FFTW_plan p[ 3 ][ 3 ][ MAX_LOG_HOWMANY + 1 ];
    sparse3DFFT_range *range1[ 3 ][ 3 ], range2[ 3 ];
    FFTW_complex *scratch;
    int fft2_first, first_dim, second_dim, third_dim;
  } sparse3DFFT_plan_struct;
  
typedef sparse3DFFT_plan_struct *sparse3DFFT_plan;

typedef int ( *sparse3DFFT_nonzero_func ) ( int x[ 3 ], void *data );

extern sparse3DFFT_plan sparse3DFFT_create_plan( int nx, int ny, int nz,
					         int dir, int flags,
					         sparse3DFFT_sparsedir sparsedir,
					         sparse3DFFT_nonzero_func nonzero,
					         void *nonzero_data,
					         FFTW_complex *data_in,
					         FFTW_complex *data_out );
					      
extern void sparse3DFFT_destroy_plan( sparse3DFFT_plan p );

extern void sparse3DFFT( sparse3DFFT_plan p, FFTW_complex *data_in, FFTW_complex *data_out );

#endif /* SPARSE3DFFT_H */
