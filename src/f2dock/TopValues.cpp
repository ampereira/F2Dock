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
#include <math.h>

#include "TopValues.h"

using CCVOpenGLMath::Matrix;
using CCVOpenGLMath::Vector;

#define gridId( v, nB, nF ) ( ( ( ( v ) >> ( ( nB ) + ( nB ) ) ) * ( nF ) + ( ( ( v ) >> ( nB ) ) & ( ( 1 << ( nB ) ) - 1 ) ) ) * ( nF ) + ( ( v ) & ( ( 1 << ( nB ) ) - 1 ) ) )

/********************************************************************/
/*                                                                  */
/*  Constructor                                                     */
/*  The constructor tells us how many of the top positions needs to */
/*  be searched for.                                                */
/*                                                                  */
/********************************************************************/
TopValues::TopValues( int numberOfPositions, int gridSize )
{
   m_NumberOfPositions = numberOfPositions;
   m_CurrentNumberOfPositions = 0;
   m_GridSize = gridSize;
   m_ValuePosition3D = new ValuePosition3D[ numberOfPositions ];
   m_Heap = new PairingHeap( true, true );

   int i;
   m_Heap->Find_Min( i, m_CurMin );
}

/********************************************************************/
/*                                                                  */
/*  Destructor                                                      */
/*                                                                  */
/********************************************************************/
TopValues::~TopValues( )
{
   delete [ ] m_ValuePosition3D; m_ValuePosition3D = 0;
   delete m_Heap; m_Heap = 0;
}


int TopValues::getCurrentNumberOfPositions( )
{
   return m_CurrentNumberOfPositions;
}


int TopValues::getGridSize( )
{
   return m_GridSize;
}


int TopValues::getCurMin( )
{
   return m_CurMin;
}


inline void TopValues::deleteMin( )
{
   if ( m_CurrentNumberOfPositions > 0 )
     {
      int i;
      VAL_TYPE v;

      m_Heap->Delete_Min( i, v );
      m_Heap->Find_Min( i, m_CurMin );
//      m_CurMin -= ( INF >> 1 );

      m_CurrentNumberOfPositions--;
     }
}


bool TopValues::extractMin( ValuePosition3D &sol )
{
   if ( m_CurrentNumberOfPositions > 0 )
     {
      int i;
      VAL_TYPE v;

      m_Heap->Delete_Min( i, v );

      sol = m_ValuePosition3D[ i ];

      m_CurrentNumberOfPositions--;
      
      return true;
     }
   else
     {
      sol.reset( );
      
      return false;      
     }
}



void TopValues::extractMin( double* value,
                            double* realValue, double* skinSkinRealValue, double* coreCoreRealValue, double* skinCoreRealValue,
                            double* unrealValue, double* skinSkinUnrealValue, double* coreCoreUnrealValue, double* skinCoreUnrealValue,
                            double* elecValue,
                            double* hbondValue,
                            double* hydrophobicityValue,
                            double* vdWPotential,
                            int* nClashes,
                            double* simpComp,
                            double* pGsol,
                            double* pGsolH,
                            double* delDispE,
                            double* x, double* y, double* z, int* rotationIndex, int* fineRotationIndex, int* conformationIndex )
{
   if ( m_CurrentNumberOfPositions > 0 )
     {
      int i;
      VAL_TYPE v;

      m_Heap->Delete_Min( i, v );

      *value = m_ValuePosition3D[ i ].m_Value;
      *realValue = m_ValuePosition3D[ i ].m_RealValue;
      *skinSkinRealValue = m_ValuePosition3D[ i ].m_SkinSkinRealValue;
      *coreCoreRealValue = m_ValuePosition3D[ i ].m_CoreCoreRealValue;
      *skinCoreRealValue = m_ValuePosition3D[ i ].m_SkinCoreRealValue;
      *unrealValue = m_ValuePosition3D[ i ].m_ImaginaryValue;
      *skinSkinUnrealValue = m_ValuePosition3D[ i ].m_SkinSkinImaginaryValue;
      *coreCoreUnrealValue = m_ValuePosition3D[ i ].m_CoreCoreImaginaryValue;
      *skinCoreUnrealValue = m_ValuePosition3D[ i ].m_SkinCoreImaginaryValue;
      *elecValue = m_ValuePosition3D[ i ].m_elecValue;
      *hbondValue = m_ValuePosition3D[ i ].m_hbondValue;
      *hydrophobicityValue = m_ValuePosition3D[ i ].m_hydrophobicityValue;
      *vdWPotential = m_ValuePosition3D[ i ].m_vdWPotential;
      *nClashes = m_ValuePosition3D[ i ].m_nClashes;
      *simpComp = m_ValuePosition3D[ i ].m_simpComp;
      *pGsol = m_ValuePosition3D[ i ].m_pGsol;
      *pGsolH = m_ValuePosition3D[ i ].m_pGsolH;
      *delDispE = m_ValuePosition3D[ i ].m_delDispE;
      *x = m_ValuePosition3D[ i ].m_Translation[ 0 ];
      *y = m_ValuePosition3D[ i ].m_Translation[ 1 ];
      *z = m_ValuePosition3D[ i ].m_Translation[ 2 ];
      *rotationIndex = m_ValuePosition3D[ i ].m_RotationIndex;
      *fineRotationIndex = m_ValuePosition3D[ i ].m_FineRotationIndex;
      *conformationIndex = m_ValuePosition3D[ i ].m_ConformationIndex;

      m_CurrentNumberOfPositions--;
     }
   else
     {
      *value = -1000000000;
     }
}


void TopValues::extractMin( double* value,
                            double *realValue, double* unrealValue,
                            double* elecValue,
                            double* hbondValue,
                            double *hydrophobicityValue,
                            double* vdWPotential,
                            int *nClashes,
                            double* simpComp,
                            double* pGsol,
                            double* pGsolH,
                            double* delDispE,
                            double* x, double* y, double* z, int* rotationIndex, int* fineRotationIndex, int* conformationIndex )
{
   double skinSkinRealValue, coreCoreRealValue, skinCoreRealValue;
   double skinSkinUnrealValue, coreCoreUnrealValue, skinCoreUnrealValue;

   extractMin( value, realValue, &skinSkinRealValue, &coreCoreRealValue, &skinCoreRealValue,
                      unrealValue, &skinSkinUnrealValue, &coreCoreUnrealValue, &skinCoreUnrealValue,
                      elecValue,
                      hbondValue,
                      hydrophobicityValue,
                      vdWPotential,
                      nClashes,
                      simpComp,
                      pGsol,
                      pGsolH,
                      delDispE,
                      x, y, z, rotationIndex, fineRotationIndex, conformationIndex );
}


inline bool TopValues::replaceMin( ValuePosition3D &sol )
{
   if ( m_CurrentNumberOfPositions > 0 )
     {
      int i;
      VAL_TYPE v;

      m_Heap->Delete_Min( i, v );

      m_ValuePosition3D[ i ] = sol;

      m_Heap->Insert( i, ( VAL_TYPE ) sol.m_Value );

      m_Heap->Find_Min( i, m_CurMin );
      
      return true;
     }
   else return false;  
}



inline void TopValues::replaceMin( double value,
                                   double realValue, double skinSkinRealValue, double coreCoreRealValue, double skinCoreRealValue,
                                   double unrealValue, double skinSkinUnrealValue, double coreCoreUnrealValue, double skinCoreUnrealValue,
                                   double elecValue,
                                   double hbondValue,
                                   double hydrophobicityValue,
                                   double vdWPotential,
                                   int nClashes,
                                   double simpComp,
                                   double pGsol,
                                   double pGsolH,
                                   double delDispE,
                                   double x, double y, double z, int rotationIndex, int fineRotationIndex, int conformationIndex )
{
   if ( m_CurrentNumberOfPositions > 0 )
     {
      int i;
      VAL_TYPE v;

      m_Heap->Delete_Min( i, v );

      m_ValuePosition3D[ i ].m_Value = value;
      m_ValuePosition3D[ i ].m_RealValue = realValue;
      m_ValuePosition3D[ i ].m_SkinSkinRealValue = skinSkinRealValue;
      m_ValuePosition3D[ i ].m_CoreCoreRealValue = coreCoreRealValue;
      m_ValuePosition3D[ i ].m_SkinCoreRealValue = skinCoreRealValue;
      m_ValuePosition3D[ i ].m_ImaginaryValue = unrealValue;
      m_ValuePosition3D[ i ].m_SkinSkinImaginaryValue = skinSkinUnrealValue;
      m_ValuePosition3D[ i ].m_CoreCoreImaginaryValue = coreCoreUnrealValue;
      m_ValuePosition3D[ i ].m_SkinCoreImaginaryValue = skinCoreUnrealValue;
      m_ValuePosition3D[ i ].m_elecValue = elecValue;
      m_ValuePosition3D[ i ].m_hbondValue = hbondValue;
      m_ValuePosition3D[ i ].m_hydrophobicityValue = hydrophobicityValue;
      m_ValuePosition3D[ i ].m_vdWPotential = vdWPotential;
      m_ValuePosition3D[ i ].m_nClashes = nClashes;
      m_ValuePosition3D[ i ].m_simpComp = simpComp;
      m_ValuePosition3D[ i ].m_pGsol = pGsol;
      m_ValuePosition3D[ i ].m_pGsolH = pGsolH;
      m_ValuePosition3D[ i ].m_delDispE = delDispE;
      m_ValuePosition3D[ i ].m_Translation[ 0 ] = x;
      m_ValuePosition3D[ i ].m_Translation[ 1 ] = y;
      m_ValuePosition3D[ i ].m_Translation[ 2 ] = z;
      m_ValuePosition3D[ i ].m_RotationIndex = rotationIndex;
      m_ValuePosition3D[ i ].m_FineRotationIndex = fineRotationIndex;
      m_ValuePosition3D[ i ].m_ConformationIndex = conformationIndex;

      m_Heap->Insert( i, ( VAL_TYPE ) value );

      m_Heap->Find_Min( i, m_CurMin );
     }
}


inline void TopValues::replaceMin( double value,
                                   double realValue, double skinSkinRealValue, double coreCoreRealValue, double skinCoreRealValue,
                                   double unrealValue, double skinSkinUnrealValue, double coreCoreUnrealValue, double skinCoreUnrealValue,
                                   double elecValue,
                                   double hbondValue,
                                   double hydrophobicityValue,
                                   double vdWPotential,
                                   int nClashes,
                                   double simpComp,
                                   double pGsol,
                                   double pGsolH,
                                   double delDispE,
                                   double x, double y, double z, int rotationIndex, int fineRotationIndex )
{
   replaceMin( value, realValue, skinSkinRealValue, coreCoreRealValue, skinCoreRealValue,
               unrealValue, skinSkinUnrealValue, coreCoreUnrealValue, skinCoreUnrealValue,
               elecValue,
               hbondValue,
               hydrophobicityValue,
               vdWPotential,
               nClashes,
               simpComp,
               pGsol,
               pGsolH,
               delDispE,
               x, y, z, rotationIndex, fineRotationIndex, 0 );
}


inline void TopValues::replaceMin( double value,
                                   double realValue, double unrealValue,
                                   double elecValue,
                                   double hbondValue,
                                   double hydrophobicityValue,
                                   double vdWPotential,
                                   int nClashes,
                                   double simpComp,
                                   double pGsol,
                                   double pGsolH,
                                   double delDispE,
                                   double x, double y, double z, int rotationIndex, int fineRotationIndex, int conformationIndex )
{
   replaceMin( value, realValue, 0, 0, 0, unrealValue, 0, 0, 0, elecValue, hbondValue, hydrophobicityValue, vdWPotential, nClashes, simpComp,
               pGsol, pGsolH, delDispE, x, y, z, rotationIndex, fineRotationIndex, conformationIndex );
}


inline void TopValues::replaceMin( double value,
                                   double realValue, double unrealValue,
                                   double elecValue,
                                   double hbondValue,
                                   double hydrophobicityValue,
                                   double vdWPotential,
                                   int nClashes,
                                   double simpComp,
                                   double pGsol,
                                   double pGsolH,
                                   double delDispE,
                                   double x, double y, double z, int rotationIndex, int fineRotationIndex )
{
   replaceMin( value, realValue, 0, 0, 0, unrealValue, 0, 0, 0, elecValue, hbondValue, hydrophobicityValue, vdWPotential, nClashes, simpComp,
               pGsol, pGsolH, delDispE, x, y, z, rotationIndex, fineRotationIndex, 0 );
}


inline bool TopValues::insert( ValuePosition3D &sol )
{
   if ( m_CurrentNumberOfPositions < m_NumberOfPositions )
     {
      int i = m_CurrentNumberOfPositions++;

      m_ValuePosition3D[ i ] = sol;

      m_Heap->Insert( i, ( VAL_TYPE ) sol.m_Value );
      
      return true;
     }
   else return false;  
}


inline void TopValues::insert( double value,
                               double realValue, double skinSkinRealValue, double coreCoreRealValue, double skinCoreRealValue,
                               double unrealValue, double skinSkinUnrealValue, double coreCoreUnrealValue, double skinCoreUnrealValue,
                               double elecValue,
                               double hbondValue,
                               double hydrophobicityValue,
                               double vdWPotential,
                               int nClashes,
                               double simpComp,
                               double pGsol,
                               double pGsolH,
                               double delDispE,
                               double x, double y, double z, int rotationIndex, int fineRotationIndex, int conformationIndex )
{
   if ( m_CurrentNumberOfPositions < m_NumberOfPositions )
     {
      int i = m_CurrentNumberOfPositions++;

      m_ValuePosition3D[ i ].m_Value = value;
      m_ValuePosition3D[ i ].m_RealValue = realValue;
      m_ValuePosition3D[ i ].m_SkinSkinRealValue = skinSkinRealValue;
      m_ValuePosition3D[ i ].m_CoreCoreRealValue = coreCoreRealValue;
      m_ValuePosition3D[ i ].m_SkinCoreRealValue = skinCoreRealValue;
      m_ValuePosition3D[ i ].m_ImaginaryValue = unrealValue;
      m_ValuePosition3D[ i ].m_SkinSkinImaginaryValue = skinSkinUnrealValue;
      m_ValuePosition3D[ i ].m_CoreCoreImaginaryValue = coreCoreUnrealValue;
      m_ValuePosition3D[ i ].m_SkinCoreImaginaryValue = skinCoreUnrealValue;
      m_ValuePosition3D[ i ].m_elecValue = elecValue;
      m_ValuePosition3D[ i ].m_hbondValue = hbondValue;
      m_ValuePosition3D[ i ].m_hydrophobicityValue = hydrophobicityValue;
      m_ValuePosition3D[ i ].m_vdWPotential = vdWPotential;
      m_ValuePosition3D[ i ].m_nClashes = nClashes;
      m_ValuePosition3D[ i ].m_simpComp = simpComp;
      m_ValuePosition3D[ i ].m_pGsol = pGsol;
      m_ValuePosition3D[ i ].m_pGsolH = pGsolH;
      m_ValuePosition3D[ i ].m_delDispE = delDispE;
      m_ValuePosition3D[ i ].m_Translation[ 0 ] = x;
      m_ValuePosition3D[ i ].m_Translation[ 1 ] = y;
      m_ValuePosition3D[ i ].m_Translation[ 2 ] = z;
      m_ValuePosition3D[ i ].m_RotationIndex = rotationIndex;
      m_ValuePosition3D[ i ].m_FineRotationIndex = fineRotationIndex;
      m_ValuePosition3D[ i ].m_ConformationIndex = conformationIndex;

      m_Heap->Insert( i, ( VAL_TYPE ) value );
     }
}


inline void TopValues::insert( double value,
                               double realValue, double skinSkinRealValue, double coreCoreRealValue, double skinCoreRealValue,
                               double unrealValue, double skinSkinUnrealValue, double coreCoreUnrealValue, double skinCoreUnrealValue,
                               double elecValue,
                               double hbondValue,
                               double hydrophobicityValue,
                               double vdWPotential,
                               int nClashes,
                               double simpComp,
                               double pGsol,
                               double pGsolH,
                               double delDispE,
                               double x, double y, double z, int rotationIndex, int fineRotationIndex )
{
   insert( value, realValue, skinSkinRealValue, coreCoreRealValue, skinCoreRealValue,
                  unrealValue, skinSkinUnrealValue, coreCoreUnrealValue, skinCoreUnrealValue,
                  elecValue,
                  hbondValue,
                  hydrophobicityValue,
                  vdWPotential,
                  nClashes,
                  simpComp,
                  pGsol,
                  pGsolH,
                  delDispE,
                  x, y, z, rotationIndex, fineRotationIndex, 0 );
}


inline void TopValues::insert( double value, double realValue, double unrealValue,
                               double elecValue,
                               double hbondValue,
                               double hydrophobicityValue,
                               double vdWPotential,
                               int nClashes,
                               double simpComp,
                               double pGsol,
                               double pGsolH,
                               double delDispE,
                               double x, double y, double z, int rotationIndex, int fineRotationIndex, int conformationIndex )
{
   insert( value, realValue, 0, 0, 0, unrealValue, 0, 0, 0, elecValue, hbondValue, hydrophobicityValue, vdWPotential, nClashes, simpComp, 
           pGsol, pGsolH, delDispE, x, y, z, rotationIndex, fineRotationIndex, conformationIndex );
}


inline void TopValues::insert( double value, double realValue, double unrealValue,
                               double elecValue,
                               double hbondValue,
                               double hydrophobicityValue,
                               double vdWPotential,
                               int nClashes,
                               double simpComp,
                               double pGsol,
                               double pGsolH,
                               double delDispE,
                               double x, double y, double z, int rotationIndex, int fineRotationIndex )
{
   insert( value, realValue, 0, 0, 0, unrealValue, elecValue, hbondValue, hydrophobicityValue, vdWPotential, nClashes, simpComp, 
           pGsol, pGsolH, delDispE, 0, 0, 0, x, y, z, rotationIndex, fineRotationIndex, 0 );
}


inline void TopValues::insert( int valuePosIndex )
{
   if ( m_CurrentNumberOfPositions < m_NumberOfPositions )
     {
      VAL_TYPE v = ( VAL_TYPE ) m_ValuePosition3D[ valuePosIndex ].m_Value;

      m_Heap->Insert( valuePosIndex, v );

      m_CurrentNumberOfPositions++;
     }
}


inline bool TopValues::updateList( ValuePosition3D &sol )
{
   if ( m_CurrentNumberOfPositions < m_NumberOfPositions )
     {
      insert( sol );

      if ( sol.m_Value < m_CurMin ) m_CurMin = ( VAL_TYPE ) sol.m_Value;

      return true;
     }
   else
     {
      if ( sol.m_Value > m_CurMin )
        {
         replaceMin( sol );

         return true;
        }
      else return false;
     }
}


inline bool TopValues::updateList( double value,
                                   double realValue, double skinSkinRealValue, double coreCoreRealValue, double skinCoreRealValue,
                                   double unrealValue, double skinSkinUnrealValue, double coreCoreUnrealValue, double skinCoreUnrealValue,
                                   double elecValue,
                                   double hbondValue,
                                   double hydrophobicityValue,
                                   double vdWPotential,
                                   int nClashes,
                                   double simpComp,
                                   double pGsol,
                                   double pGsolH,
                                   double delDispE,
                                   double x, double y, double z, int rotationIndex, int fineRotationIndex, int conformationIndex )
{
   if ( m_CurrentNumberOfPositions < m_NumberOfPositions )
     {
      insert( value, realValue, skinSkinRealValue, coreCoreRealValue, skinCoreRealValue,
              unrealValue, skinSkinUnrealValue, coreCoreUnrealValue, skinCoreUnrealValue,
              elecValue,
              hbondValue,
              hydrophobicityValue,
              vdWPotential,
              nClashes,
              simpComp,
              pGsol,
              pGsolH,
              delDispE,
              x, y, z, rotationIndex, fineRotationIndex, conformationIndex );

      if ( value < m_CurMin ) m_CurMin = ( VAL_TYPE ) value;

      return true;
     }
   else
     {
      if ( value > m_CurMin )
        {
         replaceMin( value, realValue, skinSkinRealValue, coreCoreRealValue, skinCoreRealValue,
                     unrealValue, skinSkinUnrealValue, coreCoreUnrealValue, skinCoreUnrealValue,
                     elecValue,
                     hbondValue,
                     hydrophobicityValue,
                     vdWPotential,
                     nClashes,
                     simpComp,
                     pGsol,
                     pGsolH,
                     delDispE,
                     x, y, z, rotationIndex, fineRotationIndex, conformationIndex );

         return true;
        }
      else return false;
     }
}



inline bool TopValues::updateList( double value,
                                   double realValue, double skinSkinRealValue, double coreCoreRealValue, double skinCoreRealValue,
                                   double unrealValue, double skinSkinUnrealValue, double coreCoreUnrealValue, double skinCoreUnrealValue,
                                   double elecValue,
                                   double hbondValue,
                                   double hydrophobicityValue,
                                   double vdWPotential,
                                   int nClashes,
                                   double simpComp,
                                   double pGsol,
                                   double pGsolH,
                                   double delDispE,
                                   double x, double y, double z, int rotationIndex, int fineRotationIndex )
{
   return updateList( value, realValue, skinSkinRealValue, coreCoreRealValue, skinCoreRealValue,
                      unrealValue, skinSkinUnrealValue, coreCoreUnrealValue, skinCoreUnrealValue,
                      elecValue,
                      hbondValue,
                      hydrophobicityValue,
                      vdWPotential,
                      nClashes,
                      simpComp,
                      pGsol,
                      pGsolH,
                      delDispE,
                      x, y, z, rotationIndex, fineRotationIndex, 0 );
}


inline bool TopValues::updateList( double value, double realValue, double unrealValue,
                                   double elecValue,
                                   double hbondValue,
                                   double hydrophobicityValue,
                                   double vdWPotential,
                                   int nClashes,
                                   double simpComp,
                                   double pGsol,
                                   double pGsolH,
                                   double delDispE,
                                   double x, double y, double z, int rotationIndex, int fineRotationIndex, int conformationIndex )
{
   return updateList( value, realValue, 0, 0, 0, unrealValue, 0, 0, 0, elecValue, hbondValue, hydrophobicityValue, vdWPotential, nClashes, simpComp, 
                      pGsol, pGsolH, delDispE, x, y, z, rotationIndex, fineRotationIndex, conformationIndex );
}


inline bool TopValues::updateList( double value, double realValue, double unrealValue,
                                   double elecValue,
                                   double hbondValue,
                                   double hydrophobicityValue,
                                   double vdWPotential,
                                   int nClashes,
                                   double simpComp,
                                   double pGsol,
                                   double pGsolH,
                                   double delDispE,
                                   double x, double y, double z, int rotationIndex, int fineRotationIndex )
{
   return updateList( value, realValue, 0, 0, 0, unrealValue, 0, 0, 0, elecValue, hbondValue, hydrophobicityValue, vdWPotential, nClashes, simpComp, 
                      pGsol, pGsolH, delDispE, x, y, z, rotationIndex, fineRotationIndex, 0 );
}


bool TopValues::updateTopValues( ValuePosition3D &sol )
{
   return updateList( sol );
}


bool TopValues::updateTopValues( double value,
                                 double realValue, double skinSkinRealValue, double coreCoreRealValue, double skinCoreRealValue,
                                 double unrealValue, double skinSkinUnrealValue, double coreCoreUnrealValue, double skinCoreUnrealValue,
                                 double elecValue,
                                 double hbondValue,
                                 double hydrophobicityValue,
                                 double vdWPotential,
                                 int nClashes,
                                 double simpComp,
                                 double pGsol,
                                 double pGsolH,
                                 double delDispE,
                                 double x, double y, double z, int rotationIndex, int fineRotationIndex, int conformationIndex )
{
   return updateList( value, realValue, skinSkinRealValue, coreCoreRealValue, skinCoreRealValue,
                      unrealValue, skinSkinUnrealValue, coreCoreUnrealValue, skinCoreUnrealValue,
                      elecValue,
                      hbondValue,
                      hydrophobicityValue,
                      vdWPotential,
                      nClashes,
                      simpComp,
                      pGsol,
                      pGsolH,
                      delDispE,
                      x, y, z, rotationIndex, fineRotationIndex, conformationIndex );
}

bool TopValues::updateTopValues( double value,
                                 double realValue, double skinSkinRealValue, double coreCoreRealValue, double skinCoreRealValue,
                                 double unrealValue, double skinSkinUnrealValue, double coreCoreUnrealValue, double skinCoreUnrealValue,
                                 double elecValue,
                                 double hbondValue,
                                 double hydrophobicityValue,
                                 double vdWPotential,
                                 int nClashes,
                                 double simpComp,
                                 double pGsol,
                                 double pGsolH,
                                 double delDispE,
                                 double x, double y, double z, int rotationIndex, int fineRotationIndex )
{
   return updateList( value, realValue, skinSkinRealValue, coreCoreRealValue, skinCoreRealValue,
                      unrealValue, skinSkinUnrealValue, coreCoreUnrealValue, skinCoreUnrealValue,
                      elecValue,
                      hbondValue,
                      hydrophobicityValue,
                      vdWPotential,
                      nClashes,
                      simpComp,
                      pGsol,
                      pGsolH,
                      delDispE,
                      x, y, z, rotationIndex, fineRotationIndex, 0 );
}


bool TopValues::updateTopValues( double value, double realValue, double unrealValue,
                                 double elecValue,
                                 double hbondValue,
                                 double hydrophobicityValue,
                                 double vdWPotential,
                                 int nClashes,
                                 double simpComp,
                                 double pGsol,
                                 double pGsolH,
                                 double delDispE,
                                 double x, double y, double z, int rotationIndex, int fineRotationIndex, int conformationIndex )
{
   return updateList( value, realValue, 0, 0, 0, unrealValue, 0, 0, 0, elecValue, hbondValue, hydrophobicityValue, vdWPotential, nClashes, simpComp, 
                      pGsol, pGsolH, delDispE, x, y, z, rotationIndex, fineRotationIndex, conformationIndex );
}


bool TopValues::updateTopValues( double value, double realValue, double unrealValue,
                                 double elecValue,
                                 double hbondValue,
                                 double hydrophobicityValue,
                                 double vdWPotential,
                                 int nClashes,
                                 double simpComp,
                                 double pGsol,
                                 double pGsolH,
                                 double delDispE,
                                 double x, double y, double z, int rotationIndex, int fineRotationIndex )
{
   return updateList( value, realValue, 0, 0, 0, unrealValue, 0, 0, 0, elecValue, hbondValue, hydrophobicityValue, vdWPotential, nClashes, simpComp, 
                      pGsol, pGsolH, delDispE, x, y, z, rotationIndex, fineRotationIndex, 0 );
}



inline void TopValues::computeAndPrintRMSD( int i, int j, int k, int c,
                                            FFTW_complex* grid, double realSCWeight, double imaginarySCWeight,
                                            int rotationIndex, int fineRotationIndex, int conformationIndex,
                                            PARAMS *pr )
{
   Matrix transformation;
   double rx, ry, rz, rmsd;

   retrieveTransformation( k, j, i, conformationIndex, rotationIndex, fineRotationIndex,
                           pr,
                           rx, ry, rz,
                           transformation );

   // compute square of translation magnitude
   double trans2 = rx * rx + ry * ry + rz * rz;

   if ( trans2 <= 75.0 )  // translation is less than 5. ==> compute RMSD and print results if RMSD < 5.
     {
       rmsd = getRMSD( pr->pri, transformation );

       if ( rmsd < 5.0 )
	 {
	   double score = realSCWeight * grid[ c ][ 0 ] + imaginarySCWeight * grid[ c ][ 1 ];
	   double rss = 0.0, rcc = 0.0, rsc = 0.0, iss = 0.0, icc = 0.0, isc = 0.0;

           printf( "%16.5f %16.5f %16.5f %16.5f %16.5f %16.5f %16.5f ",
           	   score / pr->functionScaleFactor,
        	   rss / pr->functionScaleFactor,
        	   rcc / pr->functionScaleFactor,
        	   rsc / pr->functionScaleFactor,
        	   iss / pr->functionScaleFactor,
        	   icc / pr->functionScaleFactor,
        	   isc / pr->functionScaleFactor );

	   for ( int s = 0; s < 4; s++ )
              for ( int t = 0; t < 4; t++ )
                 printf( "%9.3f ", transformation.get( s, t ) );

           printf( "%2d %9.3f\n", conformationIndex, rmsd );
         }
     }
}


/********************************************************************/
/*                                                                  */
/*  Given a grid of complex values , we                             */
/*  update the current top positions with the top values from the   */
/*  grid if any are larger. We store both value and position        */
/*                                                                  */
/*  We need to make sure that the same grid size is given as input  */
/*  each time.                                                      */
/*                                                                  */
/********************************************************************/
bool TopValues::updateTopPositions( int *validOutputMap,
   FFTW_complex* grid, FFTW_complex* skinSkinGrid,
   FFTW_complex* coreCoreGrid, FFTW_complex* skinCoreGrid,
   FFTW_complex* elecGrid, double elecWeight,
   FFTW_complex* hbondGrid, double hbondWeight,
   FFTW_complex* hydrophobicityGrid, double realHydrophobicityWeight, double imaginaryHydrophobicityWeight,
   FFTW_complex* simpleComplementarityGrid,
   int rotationIndex, int fineRotationIndex, int conformationIndex, PARAMS *pr )
//
// updateTopPositions used when breakDownScores is True
//
{
//   if ( !grid ) return false;
   bool goodNode = false;

   double eps = 0.00001;
   double functionScaleFactor = pr->functionScaleFactor;
   PairingHeap* pairH = NULL;

   // same conformation of moving molecule so conf always zero
   int i, j, k, c, numGoodPeaks, conf = 0;
   double score, realScore, unrealScore, rmsd;
   bool sortPeaks = false;
   FFTW_DATA_TYPE *elecGridReal = ( FFTW_DATA_TYPE * ) elecGrid;

   if ( ( pr->pri->clusterTransRad > 0 ) || ( pr->pri->peaksPerRotation < pr->numFreq * pr->numFreq * pr->numFreq ) )
     {
      sortPeaks = true;
      pairH = new PairingHeap( true, false );
     }
     
   bool filterOn = false;
   
   if ( pr->pri->applyVdWFilter || pr->pri->applyClashFilter || pr->pri->applyPseudoGsolFilter || pr->pri->applyDispersionFilter || pr->pri->applyMiscFilter )  
      filterOn = true;

   Matrix transformation;
   double trans[ 16 ];

   for ( i = 0, c = 0, numGoodPeaks = 0; i < (pr->numFreq); i++ )
      for ( j = 0; j < (pr->numFreq); j++ )
         for ( k = 0; k < (pr->numFreq); k++ )
         {
#ifdef WITH_ALL_RMSD
// 	      // if rotation index is < 1000 compute RMSD
// 	      if (rotationIndex<1000 && pr->reDocking) {

 		double rx, ry, rz, rmsd;

                retrieveTransformation( k, j, i, conf, rotationIndex, fineRotationIndex,
                                        pr,
                                        rx, ry, rz,
                                        transformation );

		// compute square of translation magnitude
		double trans2 = rx*rx + ry*ry + rz*rz;

		if (trans2 <= 75.){  // translation is less than 5. ==> compute RMSD and print results if RMSD < 5.

		    // find final transformation by adding the translation
//		    transformation = transformation.preMultiplication( Matrix::translation( rx, ry, rz ) );

		    rmsd = getRMSD( pr->pri, transformation);
		    if (rmsd < 5.0)
		      {
			double score = grid[ c ][ 0 ];
			double rss = skinSkinGrid[ c ][ 0 ];
			double rcc = coreCoreGrid[ c ][ 0 ];
			double rsc = skinCoreGrid[ c ][ 0 ];
			double iss = skinSkinGrid[ c ][ 1 ];
			double icc = coreCoreGrid[ c ][ 1 ];
			double isc = skinCoreGrid[ c ][ 1 ];

			printf( "%16.5f %16.5f %16.5f %16.5f %16.5f %16.5f %16.5f ",
				score / functionScaleFactor,
				rss / functionScaleFactor,
				rcc / functionScaleFactor,
				rsc / functionScaleFactor,
				iss / functionScaleFactor,
				icc / functionScaleFactor,
				isc / functionScaleFactor );

			for ( int i = 0; i < 4; i++ )
			  for ( int j = 0; j < 4; j++ )
			    printf( "%9.3f ", transformation.get( i, j ) );
		        printf( "%2d %9.3f\n", conf, rmsd );

		      }
		}
#endif
              if ( validOutputMap[ c ] && ( ( grid[ c ][ 0 ] > eps ) || ( - grid[ c ][ 0 ] > eps ) ) )
                {
                  score = grid[ c ][ 0 ];
                  if ( elecWeight != 0 ) score -= elecWeight * elecGridReal[ c ];
                  if ( hbondWeight != 0 ) score += hbondWeight * hbondGrid[ c ][ 0 ];
                  double hydrophobicityScore = 0, simpleComplementarityScore = 0;

                  if ( pr->pri->hydrophobicityWeight != 0 )
                    {
                      double hval = 1.0;

                      if ( hydrophobicityGrid[ c ][ 0 ] > 0.0 ) hval = ( hydrophobicityGrid[ c ][ 1 ] / hydrophobicityGrid[ c ][ 0 ] );

                      if ( hval > 3.0 ) hval = -1.0;

                      hydrophobicityScore = pr->pri->hydrophobicityWeight * hval * pr->functionScaleFactor;
                    }
                  else if ( ( realHydrophobicityWeight != 0 ) || ( imaginaryHydrophobicityWeight != 0 ) )
                    {
                      hydrophobicityScore = ( realHydrophobicityWeight * hydrophobicityGrid[ c ][ 0 ]
                               + imaginaryHydrophobicityWeight * hydrophobicityGrid[ c ][ 1 ] );
                    }

                  score += hydrophobicityScore;

                  if ( ( pr->pri->simpleShapeWeight > 0 ) || ( pr->pri->simpleShapeWeight > 0 ) )
                    {
                      simpleComplementarityScore = simpleComplementarityGrid[ c ][ 0 ];
                      score += simpleComplementarityScore;
                    }

                  if ( sortPeaks )
                    {
                     if ( ( m_CurrentNumberOfPositions < m_NumberOfPositions ) || ( score > m_CurMin ) )
                       {
                        pairH->Insert( c, -score );
                        numGoodPeaks++;
                       }
                    }
                  else
                    {
                      realScore = skinSkinGrid[ c ][ 0 ] + coreCoreGrid[ c ][ 0 ] + skinCoreGrid[ c ][ 0 ];
                      unrealScore = skinSkinGrid[ c ][ 1 ] + coreCoreGrid[ c ][ 1 ] + skinCoreGrid[ c ][ 1 ];

                      if ( ( m_CurrentNumberOfPositions < m_NumberOfPositions ) || ( score > m_CurMin ) )
                        {
//                         if ( pr->pri->applyVdWFilter || pr->pri->applyClashFilter )
//                           {
                            double vdwp = 0, vdwpQ;
                            int nclashes = 0, nclose;

                            retrieveTransformation( k, j, i, conf, rotationIndex, fineRotationIndex,
                                                    pr,
                                                    transformation, trans );

                            bool filtered = false;
                            
//                            if ( pr->pri->applyVdWFilter || pr->pri->applyClashFilter )
//                               {
//                                 computeSingleVdWPotential( pr->pri->vdWParams, transformation, &vdwp, &vdwpQ, &nclashes, &nclose );
//
//                                 if ( !pr->pri->applyVdWFilter ) vdwp = vdwpQ = 0;
//                                 if ( !pr->pri->applyClashFilter ) nclashes = 0;
//
//                                 if ( pr->pri->compQuadVdW ) vdwp = vdwpQ;
//
//
//                                 if ( ( pr->pri->applyVdWFilter && ( vdwp > pr->pri->vdWCutoffLow ) )
//                                   || ( pr->pri->applyClashFilter && ( nclashes > pr->pri->clashTolerance ) ) )
//                                     filtered = true;
//                               }      

                            double pseudoGsol = 0, pseudoGsolH = 0;

                            if ( !filtered && pr->pri->applyPseudoGsolFilter )
                               {
                                 pr->pri->pGsol->getPseudoGsol( pr->threadID - 1, trans, &pseudoGsol, &pseudoGsolH );
                                 filtered = updatePseudoGsolCutoff( pr->pri, pseudoGsol );
                                 //if ( pseudoGsolH < 0 ) filtered = true;
                               }

                            double delDispE = 0;

                            if ( !filtered && pr->pri->applyDispersionFilter )
                               {
                                 delDispE = pr->pri->dispF->getDelDispE( pr->threadID - 1, trans );
                                 if ( delDispE > pr->pri->dispersionCutoff ) filtered = true;
                               }

                            if ( !filtered )
                               {
                                 updateList( score, realScore, skinSkinGrid[ c ][ 0 ], coreCoreGrid[ c ][ 0 ], skinCoreGrid[ c ][ 0 ],
                	                     unrealScore, skinSkinGrid[ c ][ 1 ], coreCoreGrid[ c ][ 1 ], skinCoreGrid[ c ][ 1 ],
                	                     ( elecWeight != 0 ) ? elecGridReal[ c ] : 0,
                	                     ( hbondWeight != 0 ) ? hbondGrid[ c ][ 0 ] : 0,
                	                     hydrophobicityScore,
                		             vdwp,
                		             nclashes,
                		             simpleComplementarityScore,
                		             pseudoGsol,
                		             pseudoGsolH,
                		             delDispE,
                		             k, j, i, rotationIndex, fineRotationIndex, conformationIndex );
                		             
                		 goodNode = true;            
                	       }             
//                           }
//                         else
//                           {
//                              updateList( score, realScore, skinSkinGrid[ c ][ 0 ], coreCoreGrid[ c ][ 0 ], skinCoreGrid[ c ][ 0 ],
//                	                  unrealScore, skinSkinGrid[ c ][ 1 ], coreCoreGrid[ c ][ 1 ], skinCoreGrid[ c ][ 1 ],
//                	                  ( elecWeight != 0 ) ? elecGridReal[ c ] : 0,
//                	                  ( hbondWeight != 0 ) ? hbondGrid[ c ][ 0 ] : 0,
//                	                  hydrophobicityScore,
//                		          0, 0,
//                		          simpleComplementarityScore,
//                		          0, 0, 0, k, j, i, rotationIndex, fineRotationIndex, conformationIndex );
//
//                              goodNode = true;
//                	   }
                        }
    	            }
	        }
	      c++;
	 }

   if ( sortPeaks )
      {
//       quickSort( grid, 1.0, 0.0, elecGridReal, elecWeight, pr->sortedPeaks, 0, numGoodPeaks - 1 );

       for ( c = 0/*numGoodPeaks*/; c < pr->numFreq * pr->numFreq * pr->numFreq; c++ )
          pr->sortedPeaks[ c ] = 0;//1;

       double gridSpacing = 1.0 / ( ( pr->numFreq - 1 ) * pr->scaleB );

       int checkedForFiltering = 0;

       for ( int p = 0, q = 0; p < numGoodPeaks; p++ )
         {
          pairH->Delete_Min( c, score );
          score = -score;

          if ( pr->sortedPeaks[ c ] ) continue;

          if ( ( m_CurrentNumberOfPositions < m_NumberOfPositions ) || ( score > m_CurMin ) )
             {
              realScore = skinSkinGrid[ c ][ 0 ] + coreCoreGrid[ c ][ 0 ] + skinCoreGrid[ c ][ 0 ];
              unrealScore = skinSkinGrid[ c ][ 1 ] + coreCoreGrid[ c ][ 1 ] + skinCoreGrid[ c ][ 1 ];

              double hydrophobicityScore = 0, simpleComplementarityScore = 0;

              if ( pr->pri->hydrophobicityWeight != 0 )
                 {
                   double hval = 1.0;

                   if ( hydrophobicityGrid[ c ][ 0 ] > 0.0 ) hval = ( hydrophobicityGrid[ c ][ 1 ] / hydrophobicityGrid[ c ][ 0 ] );

                   if ( hval > 3.0 ) hval = -1.0;

                   hydrophobicityScore = pr->pri->hydrophobicityWeight * hval * pr->functionScaleFactor;
                 }
              else if ( ( realHydrophobicityWeight != 0 ) || ( imaginaryHydrophobicityWeight != 0 ) )
                 {
                   hydrophobicityScore = realHydrophobicityWeight * hydrophobicityGrid[ c ][ 0 ]
                                       + imaginaryHydrophobicityWeight * hydrophobicityGrid[ c ][ 1 ];
                 }

              if ( ( pr->pri->simpleShapeWeight > 0 ) || ( pr->pri->simpleShapeWeight > 0 ) )
                 {
                   simpleComplementarityScore = simpleComplementarityGrid[ c ][ 0 ];
                   score += simpleComplementarityScore;
                 }

              int cc = c;

              k = cc % pr->numFreq;
              cc /= pr->numFreq;
              j = cc % pr->numFreq;
              cc /= pr->numFreq;
              i = cc;
              
              if ( filterOn && ( checkedForFiltering++ > pr->pri->filterDepth ) ) break;

//              if ( pr->pri->applyVdWFilter || pr->pri->applyClashFilter )
//                {
                 if ( checkedForFiltering++ > pr->pri->filterDepth ) break;

                 double vdwp = 0, vdwpQ;
                 int nclashes = 0, nclose;

                 retrieveTransformation( k, j, i, conf, rotationIndex, fineRotationIndex,
                                         pr,
                                         transformation, trans );

                 bool filtered = false;

//                 if ( pr->pri->applyVdWFilter || pr->pri->applyClashFilter )
//                    {
//                       computeSingleVdWPotential( pr->pri->vdWParams, transformation, &vdwp, &vdwpQ, &nclashes, &nclose );
//      
//                       if ( !pr->pri->applyVdWFilter ) vdwp = vdwpQ = 0;
//                       if ( !pr->pri->applyClashFilter ) nclashes = 0;
//      
//                       if ( pr->pri->compQuadVdW ) vdwp = vdwpQ;
//      
//                                   
//                       if ( ( pr->pri->applyVdWFilter && ( vdwp > pr->pri->vdWCutoffLow ) )
//                         || ( pr->pri->applyClashFilter && ( nclashes > pr->pri->clashTolerance ) ) )
//                              filtered = true;
//                    }    
//
                 double pseudoGsol = 0, pseudoGsolH = 0;

                 if ( !filtered && pr->pri->applyPseudoGsolFilter )
                    {
                      pr->pri->pGsol->getPseudoGsol( pr->threadID - 1, trans, &pseudoGsol, &pseudoGsolH );
                      filtered = updatePseudoGsolCutoff( pr->pri, pseudoGsol );
                      //if ( pseudoGsolH < 0 ) filtered = true;
                    }
                    
                 double delDispE = 0;

                 if ( !filtered && pr->pri->applyDispersionFilter )
                    {
                      delDispE = pr->pri->dispF->getDelDispE( pr->threadID - 1, trans );
                      if ( delDispE > pr->pri->dispersionCutoff ) filtered = true;
                    }
                    
                 if ( !filtered )
                   {
                     updateList( score, realScore, skinSkinGrid[ c ][ 0 ], coreCoreGrid[ c ][ 0 ], skinCoreGrid[ c ][ 0 ],
    	                         unrealScore, skinSkinGrid[ c ][ 1 ], coreCoreGrid[ c ][ 1 ], skinCoreGrid[ c ][ 1 ],
    	                         ( elecWeight != 0 ) ? elecGridReal[ c ] : 0,
       	                         ( hbondWeight != 0 ) ? hbondGrid[ c ][ 0 ] : 0,
    	                         hydrophobicityScore,
    		                 vdwp,
    		                 nclashes,
       		                 simpleComplementarityScore,
       		                 pseudoGsol,
       		                 pseudoGsolH,
       		                 delDispE,
     		                 k, j, i, rotationIndex, fineRotationIndex, conformationIndex );
     		                 
     		     goodNode = true;            
     		   }              
    		 else continue;
//                }
//              else
//                {
//                  updateList( score, realScore, skinSkinGrid[ c ][ 0 ], coreCoreGrid[ c ][ 0 ], skinCoreGrid[ c ][ 0 ],
//    	                      unrealScore, skinSkinGrid[ c ][ 1 ], coreCoreGrid[ c ][ 1 ], skinCoreGrid[ c ][ 1 ],
//    	                      ( elecWeight != 0 ) ? elecGridReal[ c ] : 0,
//    	                      ( hbondWeight != 0 ) ? hbondGrid[ c ][ 0 ] : 0,
//    	                      hydrophobicityScore,
//    		              0, 0,
//                	      simpleComplementarityScore,
//                	      0, 0, k, j, i, rotationIndex, fineRotationIndex, conformationIndex );
//
//                  goodNode = true;
//    		}

              if ( ++q == pr->pri->peaksPerRotation ) break;

              markCluster( i, j, k, pr->sortedPeaks, pr->numFreq, gridSpacing, pr->pri->clusterTransRad, 1 );
             }
         }
      }

   return goodNode;
}


bool TopValues::updateTopPositions( int *validOutputMap,
   FFTW_complex* grid, FFTW_complex* skinSkinGrid,
   FFTW_complex* coreCoreGrid, FFTW_complex* skinCoreGrid,
   FFTW_complex* elecGrid, double elecWeight,
   FFTW_complex* hbondGrid, double hbondWeight,
   FFTW_complex* hydrophobicityGrid, double realHydrophobicityWeight, double imaginaryHydrophobicityWeight,
   FFTW_complex* simpleComplementarityGrid,
   int rotationIndex, int fineRotationIndex, PARAMS *pr )
{
  return updateTopPositions( validOutputMap, grid, skinSkinGrid, coreCoreGrid, skinCoreGrid, elecGrid, elecWeight, hbondGrid, hbondWeight,
                             hydrophobicityGrid, realHydrophobicityWeight, imaginaryHydrophobicityWeight,
                             simpleComplementarityGrid, rotationIndex, fineRotationIndex, 0, pr );
}



bool TopValues::updateTopPositions( int *validOutputMap, FFTW_complex* grid, double realSCWeight, double imaginarySCWeight,
                                    FFTW_complex* elecGrid, double elecWeight, FFTW_complex* hbondGrid, double hbondWeight,
                                    FFTW_complex* hydrophobicityGrid, double realHydrophobicityWeight, double imaginaryHydrophobicityWeight,
                                    FFTW_complex* simpleComplementarityGrid,
                                    int rotationIndex, int fineRotationIndex, int conformationIndex, PARAMS *pr )
//
// updateTopPositions used when breakDownScores is False
//
{
   bool goodNode = false;

   PairingHeap* pairH = new PairingHeap( true, false );
   PG *clustPG = pr->clustPG; //( pr->pri->clusterTransRad > 0 ) ? ( new PG( 10.0, - pr->numFreq * pr->pri->gridSpacing, 5.0 ) ) : NULL;

   // same conformation of moving molecule so conf always zero
   int conf = 0;
   FFTW_DATA_TYPE *elecGridReal = ( FFTW_DATA_TYPE * ) elecGrid;

   Matrix transformation;
   double trans[ 16 ];

   double nf3 = pr->numFreq * pr->numFreq * pr->numFreq;

   double hvalPenalty = -10;   

   bool filterOn = false;
   
   if ( pr->pri->applyVdWFilter || pr->pri->applyClashFilter || pr->pri->applyPseudoGsolFilter || pr->pri->applyDispersionFilter || pr->pri->applyMiscFilter )  
      filterOn = true;
      
   int numGoodPeaks = 0;

   for ( int i = 0, c = 0; i < pr->numFreq; i++ )
      for ( int j = 0; j < pr->numFreq; j++ )
   	 for ( int k = 0; k < pr->numFreq; k++, c++ )
            {

#ifdef WITH_ALL_RMSD
              computeAndPrintRMSD( i, j, k, c,
                                   grid, realSCWeight, imaginarySCWeight,
                                   rotationIndex, fineRotationIndex, conf,
                                   pr );
#endif

              if ( validOutputMap[ c ] )
                {
                  double score = realSCWeight * grid[ c ][ 0 ] + imaginarySCWeight * grid[ c ][ 1 ];

                  if ( elecWeight != 0 ) score -= elecWeight * elecGridReal[ c ];
                  if ( hbondWeight != 0 ) score += hbondWeight * hbondGrid[ c ][ 0 ];

                  if ( pr->pri->hydrophobicityWeight != 0 )
                    {
                      double hval = hvalPenalty;

                      if ( pr->pri->useInterfacePropensity )
                        {
                          if ( ( hydrophobicityGrid[ c ][ 1 ] > pr->pri->hydroRatioNumeratorLow * pr->functionScaleFactor ) 
                            && ( hydrophobicityGrid[ c ][ 1 ] < pr->pri->hydroRatioNumeratorHigh * pr->functionScaleFactor ) 
                            && ( hydrophobicityGrid[ c ][ 0 ] > pr->pri->hydroRatioDenominatorLow * pr->functionScaleFactor )
                            && ( hydrophobicityGrid[ c ][ 0 ] < pr->pri->hydroRatioDenominatorHigh * pr->functionScaleFactor ) ) 
                               hval = ( hydrophobicityGrid[ c ][ 1 ] / hydrophobicityGrid[ c ][ 0 ] );

//                          if ( hydrophobicityGrid[ c ][ 0 ] > 0.1 * pr->functionScaleFactor ) hval = ( hydrophobicityGrid[ c ][ 1 ] / hydrophobicityGrid[ c ][ 0 ] );
                                                     
                          if ( ( hval < pr->pri->hydroMinRatio ) || ( hval > pr->pri->hydroRatioTolerance ) ) hval = hvalPenalty;
                          score += ( pr->pri->hydrophobicityWeight * hval * pr->functionScaleFactor
                                   + pr->pri->hydrophobicityProductWeight * hval * hydrophobicityGrid[ c ][ 1 ] );
                        }
                      else 
                        { 
                          if ( hydrophobicityGrid[ c ][ 0 ] > 0 ) hval = ( hydrophobicityGrid[ c ][ 1 ] / hydrophobicityGrid[ c ][ 0 ] );
                        
                          if ( hval > 3.0 ) hval = hvalPenalty;    
                          score += pr->pri->hydrophobicityWeight * hval * pr->functionScaleFactor;
                        }  
                    }
                  else if ( ( realHydrophobicityWeight != 0 ) || ( imaginaryHydrophobicityWeight != 0 ) )
                    {
                      score += ( realHydrophobicityWeight * hydrophobicityGrid[ c ][ 0 ]
                               + imaginaryHydrophobicityWeight * hydrophobicityGrid[ c ][ 1 ] );
                    }

                  if ( ( pr->pri->simpleShapeWeight > 0 ) || ( pr->pri->simpleChargeWeight > 0 ) )
                      score += simpleComplementarityGrid[ c ][ 0 ];

                  if ( ( m_CurrentNumberOfPositions < m_NumberOfPositions ) || ( score > m_CurMin ) )
                    {
                     pairH->Insert( c, -score );
                     numGoodPeaks++;
                    }
                }
	    }


//   int sI[ pr->pri->peaksPerRotation ], sJ[ pr->pri->peaksPerRotation ], sK[ pr->pri->peaksPerRotation ];

   int numPassed = 0, inserted = 0;
   for ( int p = 0, checkedForFiltering = 0; p < numGoodPeaks; p++ )
    {
      int c;
      double score;
      
      double hvalT = 0;

      pairH->Delete_Min( c, score );
      score = -score;

      int i = c / pr->numFreq, j, k = c % pr->numFreq;

      j = i % pr->numFreq;
      i /= pr->numFreq;

      int clusterPenalty = 0;

      if ( pr->pri->clusterTransRad > 0 )
        {
          int cVal1 = getClusterValue( i, j, k, pr->sortedPeaks, pr->numFreq, pr->pri->gridSpacing, pr->pri->clusterTransRad * 1, score );

          if ( cVal1 >= pr->pri->clusterTransSize * 1 + 0 ) clusterPenalty = 3;
          else
             {
               int cVal2 = getClusterValue( i, j, k, pr->sortedPeaks, pr->numFreq, pr->pri->gridSpacing, pr->pri->clusterTransRad * 2, score );

               if ( cVal2 >= pr->pri->clusterTransSize * 2 + 1 ) clusterPenalty = 2;
               else
                  {
                    int cVal3 = getClusterValue( i, j, k, pr->sortedPeaks, pr->numFreq, pr->pri->gridSpacing, pr->pri->clusterTransRad * 3, score );

                    if ( cVal3 >= pr->pri->clusterTransSize * 3 + 3 ) clusterPenalty = 1;
                  }
             }
        }
        
//      if ( clusterPenalty > 0 ) continue;  

      if ( ( m_CurrentNumberOfPositions < m_NumberOfPositions ) || ( score > m_CurMin ) )
        {
         double hydrophobicityScore = 0, simpleComplementarityScore = 0;

         if ( pr->pri->hydrophobicityWeight != 0 )
            {
              double hval = hvalPenalty;

              if ( pr->pri->useInterfacePropensity )
                {
                  if ( ( hydrophobicityGrid[ c ][ 1 ] > pr->pri->hydroRatioNumeratorLow * pr->functionScaleFactor )
                    && ( hydrophobicityGrid[ c ][ 1 ] < pr->pri->hydroRatioNumeratorHigh * pr->functionScaleFactor ) 
                    && ( hydrophobicityGrid[ c ][ 0 ] > pr->pri->hydroRatioDenominatorLow * pr->functionScaleFactor )
                    && ( hydrophobicityGrid[ c ][ 0 ] < pr->pri->hydroRatioDenominatorHigh * pr->functionScaleFactor ) ) 
                       hval = ( hydrophobicityGrid[ c ][ 1 ] / hydrophobicityGrid[ c ][ 0 ] );

                  if ( ( hval < pr->pri->hydroMinRatio ) || ( hval > pr->pri->hydroRatioTolerance ) ) hval = hvalPenalty;
                  
                  hydrophobicityScore = pr->pri->hydrophobicityWeight * hval * pr->functionScaleFactor
                                      + pr->pri->hydrophobicityProductWeight * hval * hydrophobicityGrid[ c ][ 1 ];
                }
              else 
                { 
                  if ( hydrophobicityGrid[ c ][ 0 ] > 0 ) hval = ( hydrophobicityGrid[ c ][ 1 ] / hydrophobicityGrid[ c ][ 0 ] );                                  
                  if ( hval > 3.0 ) hval = hvalPenalty;    
                  hydrophobicityScore = pr->pri->hydrophobicityWeight * hval * pr->functionScaleFactor;
                } 
                
              hvalT = hval;   
            }
         else if ( ( realHydrophobicityWeight != 0 ) || ( imaginaryHydrophobicityWeight != 0 ) )
            {
              hydrophobicityScore = realHydrophobicityWeight * hydrophobicityGrid[ c ][ 0 ]
                                  + imaginaryHydrophobicityWeight * hydrophobicityGrid[ c ][ 1 ];
            }


         if ( ( pr->pri->simpleShapeWeight > 0 ) || ( pr->pri->simpleChargeWeight > 0 ) )
            {
              simpleComplementarityScore = simpleComplementarityGrid[ c ][ 0 ];
            }
        
#ifdef RERANK_DEBUG        
         updateList( score, grid[ c ][ 0 ], 0, 0, 0,
	                    grid[ c ][ 1 ], 0, 0, 0,
	                    ( elecWeight != 0 ) ? elecGridReal[ c ] : 0,
	                    ( hbondWeight != 0 ) ? hbondGrid[ c ][ 0 ] : ( numPassed * pr->functionScaleFactor ), //( pseudoGsolS * pr->functionScaleFactor ),
	                    hydrophobicityScore,
		            0,
		            0,
		            simpleComplementarityScore,
		            0,
		            0,
		            0,
		            k, j, i, rotationIndex, fineRotationIndex, 0 );
#else
         updateList( score, grid[ c ][ 0 ], 0, 0, 0,
	                    grid[ c ][ 1 ], 0, 0, 0,
	                    ( elecWeight != 0 ) ? elecGridReal[ c ] : 0,
	                    ( hbondWeight != 0 ) ? hbondGrid[ c ][ 0 ] : ( hvalT * pr->functionScaleFactor ), //( pseudoGsolS * pr->functionScaleFactor ),
	                    hydrophobicityScore,
		            0,
		            0,
		            simpleComplementarityScore,
		            hydrophobicityGrid[ c ][ 1 ] / pr->functionScaleFactor,
		            hydrophobicityGrid[ c ][ 0 ] / pr->functionScaleFactor,
		            0,
		            k, j, i, rotationIndex, fineRotationIndex, numPassed );
#endif		            

         goodNode = true;

         if ( clusterPenalty == 0 ) markCluster( i, j, k, pr->sortedPeaks, pr->numFreq, pr->pri->gridSpacing, 0, score );

         if ( ++numPassed == pr->pri->peaksPerRotation ) break;
        }
    }

   delete pairH;

   return goodNode;
}



bool TopValues::updateTopPositions( int *validOutputMap, FFTW_complex* grid, double realSCWeight, double imaginarySCWeight,
                                    FFTW_complex* elecGrid, double elecWeight, FFTW_complex* hbondGrid, double hbondWeight,
                                    FFTW_complex* hydrophobicityGrid, double realHydrophobicityWeight, double imaginaryHydrophobicityWeight,
                                    FFTW_complex* simpleComplementarityGrid,
                                    int rotationIndex, int fineRotationIndex, PARAMS *pr )
{
   return updateTopPositions( validOutputMap, grid, realSCWeight, imaginarySCWeight, elecGrid, elecWeight,
                              hbondGrid, hbondWeight, hydrophobicityGrid, realHydrophobicityWeight, imaginaryHydrophobicityWeight,
                              simpleComplementarityGrid,
                              rotationIndex, fineRotationIndex, 0, pr );
}


