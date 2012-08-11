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

// SmoothingFunction.cpp: implementation of the SmoothingFunction class.
//
//////////////////////////////////////////////////////////////////////

#include <math.h>
#include "SmoothingFunction.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

SmoothingFunction::SmoothingFunction(double _alpha, int _m, int _n, int _N)
{
	alpha = _alpha;
	m = _m;
	n = _n ;
	N = _N;

	m_Phi = 0;
	m_PhiHat = 0;
	m_PhiLength = _N * ( _m + 1 );
	m_PhiHatLength = _n;

	m_Initialized = false;
}

SmoothingFunction::~SmoothingFunction()
{
	if( m_Phi )
	{
		delete [] m_Phi; m_Phi = 0;
	}
	if( m_PhiHat )
	{
		delete [] m_PhiHat; m_PhiHat = 0;
	}
}
