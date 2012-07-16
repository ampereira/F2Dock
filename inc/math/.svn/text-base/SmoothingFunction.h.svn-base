/*
  Copyright 2011 The University of Texas at Austin

        Authors: Rezaul Alam Chowdhury <shaikat@cs.utexas.edu>, Vinay Siddavanahalli <skvinay@cs.utexas.edu>
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
// SmoothingFunction.h: interface for the SmoothingFunction class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_SMOOTHINGFUNCTION_H__29D14D8B_2356_4D85_86BF_AAD558BB98E5__INCLUDED_)
#define AFX_SMOOTHINGFUNCTION_H__29D14D8B_2356_4D85_86BF_AAD558BB98E5__INCLUDED_

#include <math.h>
#include <stdio.h>

//using namespace std;

class SmoothingFunction  
{
public:
	SmoothingFunction(double alpha, int m, int n, int N);
	virtual ~SmoothingFunction();

	inline double getPhi( double index );
	inline double getPhiHat( int index );
	virtual inline double getPhiAtRealPos( double pos ) = 0;
	virtual inline double getPhiHatAtRealPos( double pos ) = 0;

protected:

	virtual bool precompute() = 0;

	double* m_Phi;
	double* m_PhiHat;
	int m_PhiLength;
	int m_PhiHatLength;

	double alpha;
	int m;
	int n;
	int N; // used to compare with large FFT.
	bool m_Initialized;
};

// this index is a real number between 0 and m+1.
// When multiplied by N, this should be an integer. 
// This is useful to compare with a FFT method.
double SmoothingFunction::getPhi( double index )
{
	if( !m_Initialized ) precompute();

	if( !m_Phi ) return 0;
	if( (index<-(m+1)) || (index >= m+1) ) 
	  {
//	   fprintf( stderr, "( index = %lf, m = %d )\n", index, m );
	   return 0;
	  } 

	//printf(" AAAA %d\n", (int) fabs((double)(index*N)));
	return m_Phi[(int) fabs((double)(index*N))];
}

double SmoothingFunction::getPhiHat( int index )
{
	if( !m_Initialized ) precompute();

	if( !m_PhiHat ) return 0;
	
	if( (index<0) || (index >= m_PhiHatLength) ) return 0;

	return m_PhiHat[index];
}

#endif // !defined(AFX_SMOOTHINGFUNCTION_H__29D14D8B_2356_4D85_86BF_AAD558BB98E5__INCLUDED_)
