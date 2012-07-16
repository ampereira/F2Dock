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
// Gaussian.h: interface for the Gaussian class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_GAUSSIAN_H__C3E0B6A2_14F1_4872_9F56_D70142C89D07__INCLUDED_)
#define AFX_GAUSSIAN_H__C3E0B6A2_14F1_4872_9F56_D70142C89D07__INCLUDED_

#include "SmoothingFunction.h"

class Gaussian : public SmoothingFunction  
{
public:
	Gaussian(double alpha, int m, int n, int N);
	virtual ~Gaussian();
	
	inline double getPhiAtRealPos( double pos );
	inline double getPhiHatAtRealPos( double pos );
        
protected:
	bool precompute();
};

double Gaussian::getPhiAtRealPos( double pos )
{
	return 0;
}

double Gaussian::getPhiHatAtRealPos( double pos )
{
	return 0;
}

#endif // !defined(AFX_GAUSSIAN_H__C3E0B6A2_14F1_4872_9F56_D70142C89D07__INCLUDED_)
