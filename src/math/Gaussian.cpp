/*
  Copyright 2011 The University of Texas at Austin

        Authors: Rezaul Alam Chowdhury <shaikat@cs.utexas.edu>
        Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of F2Dock.

  F2Dock is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  F2Dock is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
*/


/*****************************************************************************/
// Gaussian.cpp: implementation of the Gaussian class.
/*****************************************************************************/

/*

We precompute the dilated, periodized Gaussian bell function 

*/

#include <math.h>
#include <stdio.h>
#include "Gaussian.h"

#ifndef M_PI
#define 	M_PI   3.14159265358979323846
#endif

Gaussian::Gaussian(double alpha, int m, int n, int N) : SmoothingFunction(alpha, m, n, N)
{
	precompute();
}

Gaussian::~Gaussian()
{

}

/*

Precomputation of the dilated Gaussian mentioned on page 21 (eqn 5.9) of the
following paper:

Daniel Potts and Gabriele Steidl
Fast Summation at Nonequispaced Knots by NFFTs
SIAM Journal on Scientific Computing, vol. 24, pp. 2013 - 2037, 2003.

The m_Phi values are equivalent to the Psi values in the
following paper:

Potts, D., Steidl G., and Tasche M.
Fast Fourier transforms for nonequispaced data: A tutorial.
in: Modern Sampling Theory: Mathematics and Applications, J.J. Benedetto and P. Ferreira (Eds.), 
Chapter 12, pages 249-274, 1998.

The Psi values in the paper above are truncated to [-m, m], and so
only r = 0 has any nonzero contribution to the sum.

*/

bool Gaussian::precompute( )
{
    if ( ( m_PhiLength < 1 ) || ( m_PhiHatLength < 1 ) ) return false;

    m_Phi = new double[ m_PhiLength ];        // an array of size N * ( m + 1 )
    m_PhiHat = new double[ m_PhiHatLength ];  // an array of size n = alpha * N

    double b = 2 * alpha * m / ( M_PI * ( 2.0 * alpha - 1 ) );  // b value suggested in the paper

    for ( int i = 0; i < m_PhiLength; i++ ) 
       {
        double j = i / ( ( double ) N );	
	
        m_Phi[ i ] = pow( M_PI * b, -0.5 ) * exp( - j * j / b );
		     // m_Phi[ i ] = ( pi * b )^(- 1/2) * e^( - ( (i / N)^2 ) / b )
		     // getPhi( j ) returns m_Phi[ | j * N | ] = ( pi * b )^(- 1/2) * e^( - ( j^2 ) / b ),
		     // where -m <= j <= m
       }

    for ( int i = 0; i < m_PhiHatLength; i++ )
       {
        // this j value is different from that in the paper ( ( M_PI * i ) / n );
        // we get better docking results with this value
//        double j = ( M_PI * i ) / ( 1.25 * n );
        double j = ( M_PI * i ) / n;
	
        m_PhiHat[ i ] = ( 1.0 / n ) * exp( - j * j * b );
      	               // m_PhiHat[ i ] = ( 1 / n ) * e^( - ( ( pi * i / n )^2 ) * b )
      	               // getPhiHat( i ) returns m_PhiHat[ i ]      	    
       }


    m_Initialized = true;
    
    return true;
}
