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
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/


#ifndef CCV_VALUE_POSITION_3D_H
#define CCV_VALUE_POSITION_3D_H

class ValuePosition3D
{
 public:
  
  ValuePosition3D();
  virtual ~ValuePosition3D();

  ValuePosition3D& operator=(const ValuePosition3D& copy);
  ValuePosition3D& set(const ValuePosition3D& copy);

  void reset(void);

  double m_Value;
  double m_SkinSkinRealValue, m_CoreCoreRealValue, m_SkinCoreRealValue;
  double m_SkinSkinImaginaryValue, m_CoreCoreImaginaryValue, m_SkinCoreImaginaryValue;
  double m_RealValue;  
  double m_ImaginaryValue;
  double m_elecValue;
  double m_hbondValue;  
  double m_hydrophobicityValue;  
  double m_vdWPotential;
  double m_simpComp;
  double m_pGsol;
  double m_pGsolH;
  double m_delDispE;
  double m_origScore, m_rerankerScore;
  double m_clusterPenalty;
  int m_origRank, m_rerankerRank;
  int m_nClashes;
  int m_RotationIndex;
  int m_FineRotationIndex;  
  int m_ConformationIndex;
  double m_Translation[3];
};

#endif
