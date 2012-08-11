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


#ifndef __DOCKINGMAIN_H__
#define __DOCKINGMAIN_H__

#define RERANK_DEBUG

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#if ! defined(__APPLE__)
#include <malloc.h>
#endif

#include <vector>

#include "math/Matrix.h"
#include "math/Gaussian.h"
#include "math/Vector.h"
#include "math/SmoothingFunction.h"
#include "fft-utils/fftw3.h"
#include "fft-utils/fftwPrecision.h"
#include "fft-utils/rank-fftw.h"
#include "fft-utils/sparsefft3.h"
#include "fft-utils/sparsefft3-timers.h"
#include "fft-utils/fastfft.h"
#include "fft-utils/rank-fftw.h"
#include "fast-clash/clashFilter.h"
#include "PG-range/PG.h"
#include "fast-hydro/pseudoGsol.h"
#include "fast-GB/fastDispE.h"
#include "fast-LJ/fastLJ.h"
#include "fast-resCont/resContFilter.h"

#ifdef LIBMOL_FOUND
#include <hbondFilter/hbondFilter.h>
#endif

#include "misc-ident/miscIdent.h"
#include "utils/utils.h"
#include "vol/RAWIV.h"
#include "ValuePosition3D.h"

#include <pthread.h>
#include <stdarg.h>


using CCVOpenGLMath::Matrix;
//using CCVOpenGLMath::Vector;

//int computeScoreForBoundComplex( int argc, char* argv[] );

struct PARAMS_VDW;

typedef struct
{
  int x, y, z;
  FFTW_complex v;
} NONZERO_GRIDCELLS;

typedef struct 
{
  double cx, cy, cz, r2;
} INSPHERE_DATA;

typedef struct 
{
  int n;
  int *validOutputMap;
} VALID_OUTPUT_DATA;


typedef struct 
{
  bool badNode;
  int startNeighbors, endNeighbors;
} ROTATION_GRAPH_NODE;


typedef struct
  {
    int *grid;
    int dimX, dimY, dimZ;
    double minX, minY, minZ;
    double gridSpacing;
    double cRad;
    int pw;
    double maxRatio, upScale;
  } CURVATURE_GRID;


// this structure is used in main to store arguments
typedef struct
  {
// computevdw and vdwSmoothWidth were commented out out since the vdw scores computed were not believeable and did not compare to a python implementation that was using the same libraries.
    char *id;
    
    char *paramFile;
    
    int performDocking;
    int numThreads;

    // the following members are the SAME for each thread
    int breakDownScores;  // when True skin-skin core-core and skin core
                           // components get computed

    //int computevdw;  // when True vdw interation are evaluated with AutoDockScorer

    int numberOfPositions; // number of topranking docking to save

    int gridSize;  // size of the grid (has to be a power of 2???)
    bool gridSizeSpecified; // true if gridSize is specified in the parameter file    
    double gridSpacing;  // spacing (in Angstroms) of the spatial grid
    bool enforceExactGridSpacing; // use the exact grid spacing specified by the user
    bool gridSpacingSpecified; // true if gridSpacing is specified in the parameter file
    double interpFuncExtentInAngstroms;    
    int numFreq;   // number of frequencies used in FFT
    bool numFreqSpecified; // true if numFreq is specified in the parameter file
    bool smoothSkin; // whether to use a smoothing function or not
    bool singleLayerLigandSkin;  // if set to 'true', only the surface atoms of the ligand are considered as skin atoms
    double pseudoAtomRadius; // if set to a non-negative value, this value overrides the receptor pseudo atoms radii read from the F2D file
    double pseudoAtomDistance; // distnace of the pseudo atom centers from the vdW surface of the static molecule    
    bool rotateVolume; // rotate the volume instead of the atoms of molecule B
    bool dockVolume; // dock volumes instead of molecules with explicit atoms
    
    bool useSparseFFT; // speed up FFT using the sparsity of input and/or output matrices
    bool narrowBand;   // if 'true', consider only the solutions within a narrow band in the output grid
    
//    char *gridFile; // name of the file containing efficient grid sizes for FFTW computation
    int numEfficientGridSizes; // number of efficient grid sizes
    int *efficientGridSizes;   // array of efficient grid sizes
    int minEffGridSize;  // minimum grid size cubed for efficient grid size computation
    int maxEffGridSize;  // maximum grid size cubed for efficient grid size computation
    
//    int interpFuncExtent; // ??
    int numCentersA;  // number of centers (core + skin)
    int numCentersB;  // number of centers (core + skin)
    int numberOfRotations; // set automatically when rotations is set
    
    Matrix initRot;  // initial rotation matrix ( default is < 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 > unless either explicitly set or randomRotate is 'true' )

    double distanceCutoff;
    double alpha;   // ??
    double blobbiness;
    bool blobbinessSpecified; // true if blobbiness is specified in the parameter file    
    double skinSkinWeight; // ??
    double coreCoreWeight; // ??
    double skinCoreWeight; // ??
    double realSCWeight; // ??
    double imaginarySCWeight; // ??
    
    double elecKernelVoidRad;
    double elecKernelDistLow, elecKernelDistHigh;
    double elecKernelValLow, elecKernelValHigh;  
    
    double elecScale; // ??
    double elecRadiusInGrids;  // the radius of the sphere inside which a charge is diffused using a Gaussian
    double hbondWeight;    
    double hbondDistanceCutoff;  // the distance cutoff ( in angstroms ) for hydrogen bonds
    double hydrophobicityWeight;    
    double hydrophobicityProductWeight;        
    double hydroRatioTolerance;
    double hydroMinRatio;
    double hydroRatioNumeratorLow, hydroRatioNumeratorHigh;
    double hydroRatioDenominatorLow, hydroRatioDenominatorHigh;
    bool twoWayHydrophobicity;
    double hydroPhobicPhobicWeight;    
    double hydroPhilicPhilicWeight;    
    double hydroPhobicPhilicWeight;    
    double hydroRadExt; // in Angstroms
    bool useInterfacePropensity;    // 'true' by default (instead of just hydrophobic residues)   
    bool perResidueHydrophobicity;  // 'false' by default
    int numTopHydrophobicResidues;  // only the residues with the topmost `numTopHydrophobicResidues' hydrophobicity
                                    // values will be considered 
    double staticMolHydroDistCutoff;  // distance cutoff from the skin atom centers for the atoms of the static 
                                      // molecule which contribute to hydrophobicity computation
                                      // (this is a center-to-surface distance cutoff)

    double simpleShapeWeight;      // simplified shape complementarity
    double simpleChargeWeight;     // simplified charge complementarity
    double simpleRadExt; // in Angstroms

    double clashWeight;          // weight given to each clash when added to total score
    double scoreScaleUpFactor;

    double bandwidth;     // suggested width of each band of core atoms having the same weight value 
    double gradFactor;    // the multiplicative factor for creating the gradiant of weights across bands
    
    bool curvatureWeightedStaticMol;   // construct curvature-weighted receptor skin      
    bool curvatureWeightedMovingMol;   // construct curvature-weighted ligand skin          
    double curvatureWeightingRadius;   // radius (in angstroms) of the influence zone sphere for curvature weighting
    CURVATURE_GRID cGrid;              // curvature grid for the static molecule
    bool spreadReceptorSkin;   // spread the receptor skin using a VDW-type weight function
    bool randomRotate;         // apply an initial random rotation on the ligand molecule

    char *outputFilename;
    char *staticMoleculePdb;
    char *movingMoleculePdb;
    char *staticMoleculeF2d;
    char *movingMoleculeF2d;
    
    char *staticMoleculePQR;
    char *movingMoleculePQR;    

    char *staticMoleculeSCReRaw;
    char *staticMoleculeSCImRaw;
    char *staticMoleculeElecReRaw;        
    
    char *movingMoleculeSCReRaw;
    char *movingMoleculeSCImRaw;
    char *movingMoleculeElecReRaw;        

    // needed for hbond computation. must be specified in the .inp file
    char *staticMoleculePSF;	// staticMoleculePSF <filename>
    char *movingMoleculePSF;	// movingMoleculePSF <filename>
    char *staticMoleculeMol2;	// staticMoleculeMol2 <filename>
    char *movingMoleculeMol2;	// movingMoleculeMol2 <filename>

    // files needed for hbond computation using libmol. same files can be used for all complexes. but these should still be specified in the .inp file.
    char *prmFile;	// prmFile parm.prm
    char *aprmFile;	// aprmFile atoms.0.0.6.prm.ms.3cap+0.5ace.Hr0rec
    char *rtfFile;	// rtfFile pdbamino.rtf

    char *typeA;      // string of "IIIIIIIIIIIIIEEEEE"
    char *typeB;
    char **atNamesA;
    char **atNamesB;
    char **resTypesA;
    char **resTypesB;

    int *resNumsA;
    int *resNumsB;

    float *chargesA;
    float *hydrophobicityA;
    float *hydrophobicity2A;
    float *radiiA;
    float *chargesB;
    float *hydrophobicityB;    
    float *hydrophobicity2B;    
    float *radiiB;
    
    char *hbondTypeA;      // string of "AAAAANNNNBBBBBDDDDD", A = acceptor, D = donor, B = both, N = none
    char *hbondTypeB;
    
    float *rotations;  // array of (numberOfRotations, 9) float
    ROTATION_GRAPH_NODE *rotGraph;
    std::vector< int > rotNeighbors;
    double pruneAngle;    // pruning angle in degrees for rotational search
   
    // static molecule
    double *xkAOrig;
    double *ykAOrig;
    double *zkAOrig;

    // moving molecule
    double *xkBOrig;
    double *ykBOrig;
    double *zkBOrig;

    // RMSD calculations
    int nbRMSDAtoms;  // # of atoms used to compute RMSD
    int *atNums;      // indices of atoms in moving molecules
    float *xRef;  // X-coord of atoms in reference position
    float *yRef;  // Y-coord of atoms in reference position
    float *zRef;  // Z-coord of atoms in reference position

    //double vdwSmoothWidth;  // variable used to adjust the well widthof the vdw function

    char *transformationFilename; // contains 3 x 3 transformation matrices for the moving molecule for postprocessing (e.g., vdW calculation)
    int vdWGridSize;              // grid size for vdW potential computation
    bool compQuadVdW;             // if set to true, the vdW potential is also computed using the quadratic time algorithm
    bool surfaceBasedVdW;         // if set to true, only the surface atoms of the moving molecule are used for vdW computation

    bool applyVdWFilter;          // if set to true, on-the-fly filtering based on vdW potential is performed
    double vdWCutoffLow;          // when applyVdWFilter is true and #clashes < clashTolerance / 2, poses with vdW potential > vdWCutoffLow are penalized
    double vdWCutoffHigh;         // when applyVdWFilter is true and #clashes >= clashTolerance / 2, poses with vdW potential > vdWCutoffHigh are penalized    
    double vdWEqmRadScale;        // all r_eqm values are multiplied by this factor

    double vdWWellWidth;          // smooth vdW energy well width

    double vdWTolerance;          // when pruneAngle is positive and applyVdWFilter is set, 
                                  // a node is considered good it has a neighbor with
                                  // vdWScore <= vdWCutoff + vdWTolerance
                                  
    fastLJ *ljFilter;              // octree-based Lennard-Jones potential estimator                              
                                  
    bool applyClashFilter;        // if set to true, on-the-fly filtering based on number of atomic clashes is performed
    double eqmDistFrac;           // two atoms clash if distance between atom centers < eqmDistFrac * r_eqm_XY
    int clashTolerance;           // maximum number of clashes tolerated
    clashFilter *cFilter;

    bool applyMiscFilter;         // when set to 'true' various minor filters are applied
    
    int filterDepth;              // for a given rotation, all configurations that are ranked beyond filterDepth
                                  // when sorted by decreasing order by score are automatically filtered

    struct PARAMS_VDW *vdWParams;

    bool applyPseudoGsolFilter;   // if set to true, on-the-fly filtering based on pseudo solvation energy is performed
    double pseudoGsolCutoff;      // when applyPseudoGsolFilter is set to true, all solutions with pseudo solvation energy above pseudoGsolCutoff are discarded
    double pseudoGsolWeight;      // after cutoff surviving docking poses have their scores increased (additive) by weighted pseudoGsol value
    int pseudoGsolFilterLowestRank;

    pseudoGsol *pGsol;    
    
    pthread_mutex_t pseudoGsolLock;

#ifdef LIBMOL_FOUND
    bool applyHbondFilter;        // if set to true, on-the-fly filtering based on hbond energy is performed
    double hBondFilterWeight;	  // score is updated by hBondFilterWeight*deltaHBondEnergy
    hbondFilter *hbFilter;
#endif
    
    bool applyDispersionFilter;     // if set to true, on-the-fly filtering based on solute-solvent dispersion energy is performed
    double dispersionCutoff;        // when applyDispersionFilter is set to true, all solutions with dispersion energy change below dispersionCutoff are discarded
    double dispersionWeight;        // after cutoff surviving docking poses have their scores increased (additive) by weighted disperions energy change
    double dispersionEnergyLimit;   // dispersion energy will be trancated to remain within [ dispersionEnergyLimit, dispersionEnergyLimit ]
    double dispersionMinAtomRadius; // smallest atom radius for dispersion energy calculation 

    fastGB::fastDispE *dispF;    

    double filterScaleDownFactor; // if a solution is filtered scale down its score by this factor instead of discarding it
        
    double clusterTransRad;       // radius (in angstroms) of a cluster in translational space;
                                  // if a solution is kept, no solution within clusterTransRad of it can be retained
    int clusterTransSize;         // maximum number of solutions retained within a cluster in translational space                               
    double clusterRotRad;         // radius (in degrees) of a cluster in rotational space;
                                  // if a solution is kept, no solution within clusterTransRad and clusterRotRad of it can be retained                                  
    int peaksPerRotation;         // at most how many solutions per rotation should be retained

    char complexType;             // complex type: antibody-antigen (A), enzyme-inhibitor (E), generic (G), unknown (U)

    bool rerank;

    int antibody;                 // 1 if static molecule is an antibody, -1 if moving molecule is an antibody, and 0 otherwise
    bool applyAntibodyFilter;     

    int enzyme;                   // 1 if static molecule is an enzyme, -1 if moving molecule is an enzyme, and 0 otherwise    
    bool applyEnzymeFilter;     

    bool applyResidueContactFilter;     
    
    char *spectrum;
    int numBands;
    int *bands;

    int numRerank;                // number of top positions to rerank
    double rerankerPseudoGsolWeight;
    double rerankerDispersionWeightHigh;
    double rerankerDispersionWeightLow;    
    double rerankerMinF2DockRank;
    double rerankerF2DockScoreWeight;    
   
    int control;

  } PARAMS_IN;
  
  class TopValues;


// this structure is used to pass arguments to threads
typedef struct
  {
    int threadID; // unique integer ID of the thread
    
    int confID; // deprecated and set to 0

    // the following members are the SAME for each thread
    int breakDownScores;  // when True skin-skin core-core and skin core
                           // components get computed

    int numberOfPositions; // number of topranking docking to save

    int gridSize;  // size of the grid (has to be a power of 2???)
    int numFreq;   // number of frequencies used in FFT
    int interpFuncExtent; // ??

//    double largestEdge;  // length of largest edge across the boxes containing
                         // the moving molecules and the one containg the 
                         // static molecule
    double alpha;   // ??
    double blobbiness;  // ??
    double skinSkinWeight; // ??
    double coreCoreWeight; // ??
    double skinCoreWeight; // ??
    double realSCWeight; // ?? 
    double imaginarySCWeight; // ??
    double elecScale; // ??
    double elecRadiusInGrids;  // the radius of the sphere inside which a charge is diffused using a Gaussian    
    double hbondWeight;
    double hbondDistanceCutoff;  // the distance cutoff ( in angstroms ) for hydrogen bonds
    double hydrophobicityWeight;
    double hydroPhobicPhobicWeight;    
    double hydroPhilicPhilicWeight;    
    double hydroPhobicPhilicWeight;
    double simpleShapeWeight;    
    double simpleChargeWeight;    
    double scaleA;  // ??
    double functionScaleFactor;
    double gridFactor;

    int numberOfRotations; // set automatically when rotations is set
    float *rotations;  // array of (numberOfRotations, 9) float
    
    Matrix randRot; // random rotation matrix for the initial rotation of the moving molecule
   
    // static molecule
    double *xkAOrig, *ykAOrig, *zkAOrig; // coordinates from file
    double *xkA, *ykA, *zkA; // place holder for coordinates scaled into
                             // unit cube
    int numCentersA;  // number of centers (core + skin)
    char *typeA;      // string of "IIIIIIIIIIIIIEEEEE"
    float *chargesA;
    float *radiiA;
    float *rkA;
    float *translate_A;

    // moving molecule
    double *xkBOrig, *ykBOrig, *zkBOrig; // coordinates from file
    double *xkB, *ykB, *zkB; // place holder for transformed coordinates 
                             // scaled into unit cube
    int numCentersB;  // number of centers (core + skin)
    char *typeB;
    float *chargesB;
    float *radiiB;
    float *rkB;    
    float *translate_B;
    float scaleB;

    char *hbondTypeA;      // string of "AAAAANNNNBBBBBDDDDD", A = acceptor, D = donor, B = both, N = none
    char *hbondTypeB;
    
    bool rotateVolume;
    NONZERO_GRIDCELLS *gridBCells;
    int numNonzeroGridBCells;

    NONZERO_GRIDCELLS *gridBCells_01;
    int numNonzeroGridBCells_01;
    NONZERO_GRIDCELLS *gridBCells_10;
    int numNonzeroGridBCells_10;
    NONZERO_GRIDCELLS *gridBCells_11;
    int numNonzeroGridBCells_11;

    FFTW_complex *elecGridB;   
    NONZERO_GRIDCELLS *elecGridBCells;
    int numNonzeroElecGridBCells;

    FFTW_complex *hbondGridB;   
    NONZERO_GRIDCELLS *hbondGridBCells;
    int numNonzeroHbondGridBCells;

    FFTW_complex *hydrophobicityGridB;   
    NONZERO_GRIDCELLS *hydrophobicityGridBCells;
    int numNonzeroHydrophobicityGridBCells;

    FFTW_complex *hydrophobicityTwoGridB;   
    NONZERO_GRIDCELLS *hydrophobicityTwoGridBCells;
    int numNonzeroHydrophobicityTwoGridBCells;

    FFTW_complex *simpleComplementarityGridB;   
    NONZERO_GRIDCELLS *simpleComplementarityGridBCells;
    int numNonzeroSimpleComplementarityGridBCells;
 
    int *validOutputMap;
 
    FFTW_complex *gridA, *gridB;
 
    FFTW_complex *centerFrequenciesA, *centerElecFrequenciesA, *centerHbondFrequenciesA, 
                 *centerHydrophobicityFrequenciesA, *centerHydrophobicityTwoFrequenciesA, *centerSimpleComplementarityFrequenciesA;
    FFTW_complex *centerFrequenciesB, *centerElecFrequenciesB;   
    FFTW_complex *centerFrequenciesProduct, *centerFrequenciesElecProduct;   
    FFTW_complex *sparseProfile, *sparseShapeProfile, *sparseElecProfile, *sparseHbondProfile, 
                 *sparseHydrophobicityProfile, *sparseHydrophobicityTwoProfile, *sparseSimpleComplementarityProfile;
    FFTW_plan freqPlan, elecFreqPlan;   
    sparse3DFFT_plan sparseFreqPlanBackward;//, sparseHbondFreqPlanBackward;          
    FFTW_complex *fkB, *fkBElec, *fkBHbond, *fkBHydrophobicity, *fkBHydrophobicityTwo, *fkBSimpleComplementarity;
    FFTW_complex* freqHat;
    FFTW_plan freqHatPlan;   
    FFTW_plan moreFreqPlan, moreElecFreqPlan;
    sparse3DFFT_plan sparseFreqPlanForward;//, sparseHbondFreqPlanForward;      
    FFTW_complex *smallElectrostaticsKernel;
    bool smoothSkin;    
    SmoothingFunction *smoothingFunction;
    TopValues *localTopValues;
    FFTW_complex *centerFrequenciesA_01, *centerFrequenciesA_10, *centerFrequenciesA_11;
    FFTW_complex *fkB_01, *fkB_10, *fkB_11;   
    FFTW_complex *sparseProfile_01, *sparseProfile_10, *sparseProfile_11;
    FFTW_complex *ourMoreFrequencies;   
    FFTW_complex *ourMoreFrequenciesOut;   
    
    double *sortedPeaks;    // an array used to sort peaks during translational search (used when clusterRadius > 0 and/or peaksPerRotation < numFreq^3)
    PG *clustPG;         // Packing Grid used for clustering
       
    pseudoGsol *pGsol;    
    double pGsolAvg;
    
    // RMSD calculations
    int nbRMSDAtoms;  // # of atoms used to compute RMSD
    int *atNums;      // indices of atoms in moving molecules
    float *xRef;  // X-coord of atoms in reference position
    float *yRef;  // Y-coord of atoms in reference position
    float *zRef;  // Z-coord of atoms in reference position
    //RMSDcalculator *RMSDcalc;
    //int startFineRotation, endFineRotation, numberOfFineRotations;
    //float *fineRotations;   
    PARAMS_IN *pri;

  } PARAMS;


// this structure is used to pass arguments to filtering threads
typedef struct
  {
    int threadID; // unique integer ID of the thread
    
    TopValues *TopValuesIn, *TopValuesOut;

    int numFreq;   // number of frequencies used in FFT

    double functionScaleFactor;

    float *rotations;  // array of (numberOfRotations, 9) float
    
    Matrix randRot; // random rotation matrix for the initial rotation of the moving molecule

    float *translate_A;
    float *translate_B;

    float scaleB;
    
    double pgsolSum;
    
    PARAMS_IN *pri;

  } FILTER_PARAMS;


extern "C"
{
  int dock( PARAMS_IN *pr );
  int scoreUntransformed( PARAMS_IN *pr);
  int saveGrid( PARAMS_IN *pr);
  float getRMSD(PARAMS_IN *pr, Matrix transformation );
//  int computeEffGrid( PARAMS_IN *pr );  
}

void shapeScoreSingleConfiguration( FFTW_complex *gridA, FFTW_complex *gridB, int n, int xt, int yt, int zt, 
                               double *scoreSS, double *scoreCC, double *scoreSC );
void elecScoreSingleConfiguration( int nA, double *xA, double *yA, double *zA, float *rA, float *qA, char *typeA,
                                   int nB, double *xB, double *yB, double *zB, float *rB, float *qB, char *typeB, Matrix transMat,
                                   double *elecPos, double *elecNeg );                               
void retrieveTransformation( int ix, int iy, int iz, int c, int r, int f,
                             PARAMS *pr,
                             double &rx, double &ry, double &rz,
                             Matrix &trans );
void retrieveTransformation( int ix, int iy, int iz, int c, int r, int f,
                             PARAMS *pr,
                             double &rx, double &ry, double &rz,
                             Matrix &trans, double *dtrans );
void retrieveTransformation( int ix, int iy, int iz, int c, int r, int f,
                             PARAMS *pr,
                             Matrix &trans );
void retrieveTransformation( int ix, int iy, int iz, int c, int r, int f,
                             PARAMS *pr,
                             Matrix &trans, double *dtrans );
void markCluster( int x, int y, int z, double *markedPeaks, int numFreq, double gridSpacing, double clusterTransRad, double val );                                                          
int getClusterValue( int x, int y, int z, double *markedPeaks, int numFreq, double gridSpacing, double clusterTransRad, double val );
void getMinMax( double* x, double* minVal, double* maxVal, int n );
bool getLargestSize( double* x, double* y, double* z, int n, double* maxLength );
bool getLargestPairwiseDistance( double* x, double* y, double* z, int n, double* maxLength );
bool getMaxRadius( float* r1, int n1, float* r2, int n2, float *maxRadius );
// double getTime( );
bool updatePseudoGsolCutoff( PARAMS_IN *pr, double pseudoGsol );


#endif


