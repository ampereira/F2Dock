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


#include <stdio.h>

#include "ValuePosition3D.h"

/********************************************************************/
/*                                                                  */
/*  Constructor                                                     */
/*                                                                  */
/********************************************************************/
ValuePosition3D::ValuePosition3D()
{
  reset( );
}

/********************************************************************/
/*                                                                  */
/*  Destructor                                                      */
/*                                                                  */
/********************************************************************/
ValuePosition3D::~ValuePosition3D()
{

}


ValuePosition3D& ValuePosition3D::operator=(const ValuePosition3D& copy)
{
	return set(copy);
}


ValuePosition3D& ValuePosition3D::set(const ValuePosition3D& copy)
{
	if (this!=&copy) {
	         m_Value = copy.m_Value;
                 m_RealValue = copy.m_RealValue;
                 m_SkinSkinRealValue = copy.m_SkinSkinRealValue;
                 m_CoreCoreRealValue = copy.m_CoreCoreRealValue;
                 m_SkinCoreRealValue = copy.m_SkinCoreRealValue;
                 m_ImaginaryValue = copy.m_ImaginaryValue;
                 m_SkinSkinImaginaryValue = copy.m_SkinSkinImaginaryValue;
                 m_CoreCoreImaginaryValue = copy.m_CoreCoreImaginaryValue;
                 m_SkinCoreImaginaryValue = copy.m_SkinCoreImaginaryValue;
                 m_elecValue = copy.m_elecValue;
                 m_hbondValue = copy.m_hbondValue;
                 m_hydrophobicityValue = copy.m_hydrophobicityValue;
                 m_vdWPotential = copy.m_vdWPotential;
                 m_nClashes = copy.m_nClashes;
                 m_simpComp = copy.m_simpComp;
                 m_pGsol = copy.m_pGsol;
                 m_pGsolH = copy.m_pGsolH;
                 m_delDispE = copy.m_delDispE;
                 m_Translation[ 0 ] = copy.m_Translation[ 0 ];
                 m_Translation[ 1 ] = copy.m_Translation[ 1 ];
                 m_Translation[ 2 ] = copy.m_Translation[ 2 ];
                 m_origScore = copy.m_origScore; 
                 m_rerankerScore = copy.m_rerankerScore;
                 m_clusterPenalty = copy.m_clusterPenalty;
                 m_origRank = copy.m_origRank;
                 m_rerankerRank = copy.m_rerankerRank;
                 m_RotationIndex = copy.m_RotationIndex;
                 m_FineRotationIndex = copy.m_FineRotationIndex;
                 m_ConformationIndex = copy.m_ConformationIndex;
	}
	
	return *this;
}


void ValuePosition3D::reset(void)
{
  m_Value = -100000000.0;
  m_SkinSkinRealValue = -100000000.0, m_CoreCoreRealValue = -100000000.0, m_SkinCoreRealValue = -100000000.0;
  m_SkinSkinImaginaryValue = -100000000.0, m_CoreCoreImaginaryValue = -100000000.0, m_SkinCoreImaginaryValue = -100000000.0;
  m_RealValue = -100000000.0;  
  m_ImaginaryValue = -100000000.0;
  m_elecValue = 0;
  m_hbondValue = 0;  
  m_hydrophobicityValue = 0;
  m_vdWPotential = 0.0;  
  m_simpComp = 0.0;
  m_pGsol = 0.0;
  m_pGsolH = 0.0;  
  m_delDispE = 0.0;
  m_nClashes = 0;
  m_clusterPenalty = 0.0;
  m_Translation[0] = m_Translation[1] = m_Translation[2] = 0.0;
  m_origScore = m_rerankerScore = m_Value;
  m_origRank = m_rerankerRank = -1;  
  m_RotationIndex = m_FineRotationIndex = 0;
  m_ConformationIndex = -1;
}
