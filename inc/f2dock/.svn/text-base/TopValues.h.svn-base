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


#ifndef CCV_TOP_VALUES_H
#define CCV_TOP_VALUES_H

#include <stdio.h>

#include "fftw3.h"

#include "fast-PQ/PairingHeap.h"
#include "Docking.h"
#include "ValuePosition3D.h"
#include "math/Matrix.h"
#include "math/Vector.h"
#include "PG-range/PG.h"

using CCVOpenGLMath::Matrix;
using CCVOpenGLMath::Vector;

#ifndef VAL_TYPE
   #define VAL_TYPE double
#endif

class ValuePosition3D;

class TopValues
{
public:

	TopValues( int numberOfPositions, int gridSize );
	virtual ~TopValues( );

        bool updateTopValues( ValuePosition3D &sol );
        bool updateTopValues( double value, double realValue, double skinSkinRealValue, double coreCoreRealValue, double skinCoreRealValue,
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
                                 double x, double y, double z, int rotationIndex, int fineRotationIndex, int conformationIndex );
        bool updateTopValues( double value, double realValue, double skinSkinRealValue, double coreCoreRealValue, double skinCoreRealValue,
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
                                 double x, double y, double z, int rotationIndex, int fineRotationIndex );
	bool updateTopValues( double value, double realValue, double unrealValue,
	                         double elecValue,
	                         double hbondValue,
	                         double hydrophobicityValue,
	                         double vdWPotential,
	                         int nClashes,
                                 double simpComp,
                                 double pGsol,
                                 double pGsolH,
                                 double delDispE,
	                         double x, double y, double z, int rotIndex, int fineRotIndex );
	bool updateTopValues( double value, double realValue, double unrealValue,
	                         double elecValue,
	                         double hbondValue,
	                         double hydrophobicityValue,
	                         double vdWPotential,
	                         int nClashes,
                                 double simpComp,
                                 double pGsol,
                                 double pGsolH,
                                 double delDispE,
	                         double x, double y, double z, int rotIndex, int fineRotIndex, int confIndex );

        inline void computeAndPrintRMSD( int i, int j, int k, int c,
                                         FFTW_complex* grid, double realSCWeight, double imaginarySCWeight,
                                         int rotationIndex, int fineRotationIndex, int conformationIndex,
                                         PARAMS *pr );

        bool updateTopPositions( int *validOutputMap,
		   FFTW_complex* grid, FFTW_complex* skinSkinGrid,
		   FFTW_complex* coreCoreGrid, FFTW_complex* skinCoreGrid,
		   FFTW_complex* elecGrid, double elecWeight,
		   FFTW_complex* hbondGrid, double hbondWeight,
		   FFTW_complex* hydrophobicityGrid, double realHydrophobicityWeight, double imaginaryHydrophobicityWeight,
                   FFTW_complex* simpleComplementarityGrid,
		   int rotationIndex, int fineRotationIndex, int conformationIndex, PARAMS *pr );

	bool updateTopPositions( int *validOutputMap,
    		   FFTW_complex* grid, FFTW_complex* skinSkinGrid,
    		   FFTW_complex* coreCoreGrid, FFTW_complex* skinCoreGrid,
    		   FFTW_complex* elecGrid, double elecWeight,
    		   FFTW_complex* hbondGrid, double hbondWeight,
    		   FFTW_complex* hydrophobicityGrid, double realHydrophobicityWeight, double imaginaryHydrophobicityWeight,
 		   FFTW_complex* simpleComplementarityGrid,
    		   int rotationIndex, int fineRotationIndex, PARAMS *pr );

	bool updateTopPositions( int *validOutputMap, FFTW_complex* grid, double realSCWeight, double imaginarySCWeight,
	                         FFTW_complex* elecGrid, double elecWeight,
	                         FFTW_complex* hbondGrid, double hbondWeight,
	                         FFTW_complex* hydrophobicityGrid, double realHydrophobicityWeight, double imaginaryHydrophobicityWeight,
	                         FFTW_complex* simpleComplementarityGrid,
	                         int rotationIndex, int fineRotationIndex, PARAMS *pr );

	bool updateTopPositions( int *validOutputMap, FFTW_complex* grid, double realSCWeight, double imaginarySCWeight,
	                         FFTW_complex* elecGrid, double elecWeight,
	                         FFTW_complex* hbondGrid, double hbondWeight,
	                         FFTW_complex* hydrophobicityGrid, double realHydrophobicityWeight, double imaginaryHydrophobicityWeight,
	                         FFTW_complex* simpleComplementarityGrid,
	                         int rotIndex, int fineRotIndex, int confIndex, PARAMS *pr );

	int getCurrentNumberOfPositions( );
	int getGridSize( );
	int getCurMin( );
	
        bool extractMin( ValuePosition3D &sol );	
        void extractMin( double* value, double* realValue, double* skinSkinRealValue, double* coreCoreRealValue, double* skinCoreRealValue,
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
                            double* x, double* y, double* z, int* rotationIndex, int* fineRotationIndex, int* conformationIndex );
        void extractMin( double* value, double* realValue, double* unrealValue,
                            double* elecValue,
                            double* hbondValue,
                            double* hydrophobicityValue,
                            double *vdWPotential,
                            int *nClashes,
                            double* simpComp,
                            double* pGsol,
                            double* pGsolH,
                            double* delDispE,
                            double* x, double* y, double* z, int* rotationIndex, int* fineRotationIndex, int* conformationIndex );


protected:
        inline bool updateList( ValuePosition3D &sol );
        inline bool updateList( double value, double realValue, double skinSkinRealValue, double coreCoreRealValue, double skinCoreRealValue,
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
                                   double x, double y, double z, int rotationIndex, int fineRotationIndex, int conformationIndex );
        inline bool updateList( double value, double realValue, double skinSkinRealValue, double coreCoreRealValue, double skinCoreRealValue,
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
                                   double x, double y, double z, int rotationIndex, int fineRotationIndex );
	inline bool updateList( double value, double realValue, double unrealValue,
	                           double elecValue,
	                           double hbondValue,
	                           double hydrophobicityValue,
	                           double vdWPotential,
	                           int nClashes,
                                   double simpComp,
                                   double pGsol,
                                   double pGsolH,
                                   double delDispE,
	                           double x, double y, double z, int rotIndex, int fineRotIndex );
	inline bool updateList( double value, double realValue, double unrealValue,
	                           double elecValue,
	                           double hbondValue,
	                           double hydrophobicityValue,
	                           double vdWPotential,
	                           int nClashes,
                                   double simpComp,
                                   double pGsol,
                                   double pGsolH,
                                   double delDispE,
	                           double x, double y, double z, int rotIndex, int fineRotIndex, int confIndex );

private:
	int m_NumberOfPositions; // number of top positions to capture
	int m_CurrentNumberOfPositions; // number of top positions captured
	int m_GridSize; // The size of the input grid.
	VAL_TYPE m_CurMin;
	ValuePosition3D* m_ValuePosition3D; // top positions in unsorted order
	PairingHeap* m_Heap; // an auxiliary buffer heap on valuePositions

        inline void deleteMin( );
        inline bool replaceMin( ValuePosition3D &sol );        
        inline void replaceMin( double value, double realValue, double skinSkinRealValue, double coreCoreRealValue, double skinCoreRealValue,
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
                                   double x, double y, double z, int rotationIndex, int fineRotationIndex, int conformationIndex );
        inline void replaceMin( double value, double realValue, double skinSkinRealValue, double coreCoreRealValue, double skinCoreRealValue,
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
                                   double x, double y, double z, int rotationIndex, int fineRotationIndex );
        inline void replaceMin( double value, double realValue, double unrealValue,
                                   double elecValue,
                                   double hbondValue,
                                   double hydrophobicityValue,
                                   double vdWPotential,
                                   int nClashes,
                                   double simpComp,
                                   double pGsol,
                                   double pGsolH,
                                   double delDispE,
                                   double x, double y, double z, int rotationIndex, int fineRotationIndex, int conformationIndex );
        inline void replaceMin( double value, double realValue, double unrealValue,
                                   double elecValue,
                                   double hbondValue,
                                   double hydrophobicityValue,
                                   double vdWPotential,
                                   int nClashes,
                                   double simpComp,
                                   double pGsol,
                                   double pGsolH,
                                   double delDispE,
                                   double x, double y, double z, int rotationIndex, int fineRotationIndex );
        inline bool insert( ValuePosition3D &sol );                                   
        inline void insert( double value, double realValue, double skinSkinRealValue, double coreCoreRealValue, double skinCoreRealValue,
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
                               double x, double y, double z, int rotationIndex, int fineRotationIndex, int conformationIndex );
        inline void insert( double value, double realValue, double skinSkinRealValue, double coreCoreRealValue, double skinCoreRealValue,
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
                               double x, double y, double z, int rotationIndex, int fineRotationIndex );
        inline void insert( double value, double realValue, double unrealValue,
                               double elecValue,
                               double hbondValue,
                               double hydrophobicityValue,
                               double vdWPotential,
                               int nClashes,
                               double simpComp,
                               double pGsol,
                               double pGsolH,
                               double delDispE,
                               double x, double y, double z, int rotationIndex, int fineRotationIndex, int conformationIndex );
        inline void insert( double value, double realValue, double unrealValue,
                               double elecValue,
                               double hbondValue,
                               double hydrophobicityValue,
                               double vdWPotential,
                               int nClashes,
                               double simpComp,
                               double pGsol,
                               double pGsolH,
                               double delDispE,
                               double x, double y, double z, int rotationIndex, int fineRotationIndex );
        inline void insert( int valuePosIndex );

//        inline void markCluster( int x, int y, int z, int *sortedPeaks, int numFreq, double gridSpacing, double clusterRadius );
};
extern "C"
{
  float getRMSD(PARAMS_IN *pr, Matrix transformation );
}
#endif
