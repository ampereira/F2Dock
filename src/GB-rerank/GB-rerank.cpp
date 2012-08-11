
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

#ifdef _MSC_VER
#pragma warning(disable:4786)
#endif

#include "GB-rerank/GB-rerank.h"


void getAtomsQuadsAndGpol( PARAMS_IN *params, char *pqrFile, char *quadFile, int *numAtoms, double **atomsPQR, int *numQPoints, double **qPoints, double *Gpol )
{
  fastGB::fastBornRadius fastBR( quadFile, pqrFile, false );
  
  fastBR.getAtomsPQR( numAtoms, atomsPQR );
  fastBR.getQPoints( numQPoints, qPoints );  
  
  fastBR.setNumThreads( params->numThreadsBR );
  fastBR.setEpsilon( params->epsilonBR );

  fastBR.computeBornRadii( );
  
  int nAtoms;
  double *atomsPQRR = NULL;
  
  fastBR.getAtomsPQRR( &nAtoms, &atomsPQRR );
  
  fastGB::fastGpol fastGp( nAtoms, atomsPQRR, false );

  fastGp.setNumThreads( params->numThreadsGpol );
  fastGp.setEpsilon( params->epsilonGpol );
  fastGp.useApproxMathFunctions( params->useApproxMath );

  fastGp.computeFastGpol( Gpol );
  
  freeMem( atomsPQRR );  
}



void getGpol( PARAMS_IN *params, int numAtoms, double *atomsPQR, int numQPoints, double *qPoints, double *Gpol )
{
  fastGB::fastBornRadius fastBR( numQPoints, qPoints, numAtoms, atomsPQR, false );

  fastBR.setNumThreads( params->numThreadsBR );
  fastBR.setEpsilon( params->epsilonBR );
  
  fastBR.computeBornRadii( );
  
  int nAtoms;
  double *atomsPQRR = NULL;
  
  fastBR.getAtomsPQRR( &nAtoms, &atomsPQRR );
  
  fastGB::fastGpol fastGp( nAtoms, atomsPQRR, false );

  fastGp.setNumThreads( params->numThreadsGpol );
  fastGp.setEpsilon( params->epsilonGpol );
  fastGp.useApproxMathFunctions( params->useApproxMath );
  
  fastGp.computeFastGpol( Gpol );
}


void transformAndCopyMovingMolecule( double *atomsPQR, double *qPoints, double *transMat, MOLECULE_INFO *staticMol, MOLECULE_INFO *movingMol )
{
   for ( int i = 0; i < movingMol->numAtoms; i++ )
     {
       int j = staticMol->numAtoms + i;
       
       double x = movingMol->atoms[ 5 * i + 0 ],
              y = movingMol->atoms[ 5 * i + 1 ],
              z = movingMol->atoms[ 5 * i + 2 ],
              q = movingMol->atoms[ 5 * i + 3 ],
              r = movingMol->atoms[ 5 * i + 4 ];
              
       double xx = transMat[  0 ] * x + transMat[  1 ] * y + transMat[  2 ] * z + transMat[  3 ],
              yy = transMat[  4 ] * x + transMat[  5 ] * y + transMat[  6 ] * z + transMat[  7 ],
              zz = transMat[  8 ] * x + transMat[  9 ] * y + transMat[ 10 ] * z + transMat[ 11 ];       
              
       atomsPQR[ 5 * j + 0 ] = xx;       
       atomsPQR[ 5 * j + 1 ] = yy;
       atomsPQR[ 5 * j + 2 ] = zz;
       atomsPQR[ 5 * j + 3 ] = q;
       atomsPQR[ 5 * j + 4 ] = r;                            
     }

   for ( int i = 0; i < movingMol->numQPoints; i++ )
     {
       double x  = movingMol->qPoints[ 7 * i + 0 ],
              y  = movingMol->qPoints[ 7 * i + 1 ],
              z  = movingMol->qPoints[ 7 * i + 2 ],
              nx = movingMol->qPoints[ 7 * i + 3 ],
              ny = movingMol->qPoints[ 7 * i + 4 ],
              nz = movingMol->qPoints[ 7 * i + 5 ],
              w  = movingMol->qPoints[ 7 * i + 6 ];
              
       double xx = transMat[  0 ] * x + transMat[  1 ] * y + transMat[  2 ] * z + transMat[  3 ],
              yy = transMat[  4 ] * x + transMat[  5 ] * y + transMat[  6 ] * z + transMat[  7 ],
              zz = transMat[  8 ] * x + transMat[  9 ] * y + transMat[ 10 ] * z + transMat[ 11 ];       
              
       double nxx = ( transMat[  0 ] * ( x + nx ) + transMat[  1 ] * ( y + ny ) + transMat[  2 ] * ( z + nz ) + transMat[  3 ] ) - xx,
              nyy = ( transMat[  4 ] * ( x + nx ) + transMat[  5 ] * ( y + ny ) + transMat[  6 ] * ( z + nz ) + transMat[  7 ] ) - yy,
              nzz = ( transMat[  8 ] * ( x + nx ) + transMat[  9 ] * ( y + ny ) + transMat[ 10 ] * ( z + nz ) + transMat[ 11 ] ) - zz;              
              
       qPoints[ 7 * i + 0 ] = xx;       
       qPoints[ 7 * i + 1 ] = yy;
       qPoints[ 7 * i + 2 ] = zz;
       qPoints[ 7 * i + 3 ] = nxx;
       qPoints[ 7 * i + 4 ] = nyy;
       qPoints[ 7 * i + 5 ] = nzz;
       qPoints[ 7 * i + 6 ] = w;                                   
     }
}


double computeXlateForPG( MOLECULE_INFO *staticMol, MOLECULE_INFO *movingMol )
{
   double minXYZ;
   
   for ( int i = 0; i < staticMol->numQPoints; i++ )
     {
       if ( i == 0 ) minXYZ = staticMol->qPoints[ 5 * i + 0 ];
       
       for ( int j = 0; j < 3; j++ )
         if ( staticMol->qPoints[ 5 * i + j ] < minXYZ )
            minXYZ = staticMol->qPoints[ 5 * i + j ];
     }

   double min_X_Y_Z[ 3 ];
   double max_X_Y_Z[ 3 ];

   for ( int i = 0; i < movingMol->numQPoints; i++ )
     {
       if ( i == 0 ) 
         {
           for ( int j = 0; j < 3; j++ )
             min_X_Y_Z[ j ] = max_X_Y_Z[ j ] = movingMol->qPoints[ 5 * i + j ];
         }  
       else  
         {
           for ( int j = 0; j < 3; j++ )
             {
               if ( movingMol->qPoints[ 5 * i + j ] < min_X_Y_Z[ j ] )
                  min_X_Y_Z[ j ] = movingMol->qPoints[ 5 * i + j ];
                  
               if ( movingMol->qPoints[ 5 * i + j ] > max_X_Y_Z[ j ] )
                  max_X_Y_Z[ j ] = movingMol->qPoints[ 5 * i + j ];
             }
         }         
     }

   double maxD = 0;
   
   for ( int i = 0; i < 3; i++ )
     {   
       double d = fabs( max_X_Y_Z[ i ] - min_X_Y_Z[ i ] );
      
       if ( d > maxD ) maxD = d; 
     }  
     
   minXYZ -= ( 2 * maxD );
   
   if ( minXYZ < 0 ) return -minXYZ;
   else return 0;  
}



void printRerankingParamters( PARAMS_IN *p, FILE* fp )
{
  fprintf( fp, "# RERANKING PARAMETERS:\n#\n" );	

  fprintf( fp, "# \t staticMoleculePQR = %s\n", p->staticMoleculePQR );
  fprintf( fp, "# \t movingMoleculePQR = %s\n", p->movingMoleculePQR );  

  fprintf( fp, "# \t staticMoleculeQUAD = %s\n", p->staticMoleculeQUAD );
  fprintf( fp, "# \t movingMoleculeQUAD = %s\n", p->movingMoleculeQUAD );  
  
  fprintf( fp, "# \t F2DockOutputFile = %s\n", p->F2DockOutputFile );
  fprintf( fp, "# \t rerankedOutputFile = %s\n", p->rerankedOutputFile );
  
  fprintf( fp, "# \t F2DockScoreWeight = %lf\n", p->F2DockScoreWeight );
  fprintf( fp, "# \t GpolWeight = %lf\n", p->GpolWeight );
  fprintf( fp, "# \t GnonpolWeight = %lf\n", p->GnonpolWeight );    
  
  fprintf( fp, "# \t distanceCutoff = %lf\n", p->distanceCutoff );    
  
  fprintf( fp, "# \t epsilonBR = %lf\n", p->epsilonBR );    
  fprintf( fp, "# \t epsilonGpol = %lf\n", p->epsilonGpol );    
  
  fprintf( fp, "# \t useApproxMath = %d\n", p->useApproxMath );    
  
  fprintf( fp, "# \t numThreadsBR = %d\n", p->numThreadsBR );    
  fprintf( fp, "# \t numThreadsGpol = %d\n", p->numThreadsGpol );    
  
  fprintf( fp, "# \t numSol = %d\n", p->numSol );    
  
  fprintf( fp, "# \t spectrum = " );      

  for ( int i = 0; i < p->numBands; i++ )
      fprintf( fp, "%d:%d%s", ( i > 0 ) ? ( p->bands[ 2 * i ] - p->bands[ 2 * i - 2 ] ) : p->bands[ 2 * i ], p->bands[ 2 * i + 1 ], ( i < p->numBands - 1 ) ? "-" : "" );

  fprintf( fp, "\n#\n" );	
  fprintf( fp, "#\n" );
}


#define freeClose1( ) { fclose( ifp );  }
#define freeClose2( ) { freeClose1( ); fclose( ofp ); }
#define freeClose3( ) { freeClose2( ); freeMem( atomsPQR ); }
#define freeClose4( ) { freeClose3( ); freeMem( qPoints ); }
#define freeClose5( ) { freeClose4( ); freeMem( staticQPoints ); freeMem( movingQPoints ); }
#define freeClose6( ) { freeClose5( ); delete staticPG; delete pairH; }


bool rerankF2DockOutput( PARAMS_IN *params, MOLECULE_INFO *staticMol, MOLECULE_INFO *movingMol )
{
   int range[ ] = { 1, 10, 100, 1000, 10000, 100000 };
   int nRange = sizeof( range ) / sizeof( range[ 0 ] );
   int hitsInRange[ nRange ];
   double rmsdGood = 5.0;
   int numGoodPeaks = 0;
   SOLUTION_INFO *solInfo = NULL;
   int rankMin, indexMinRMSD, rankMinRMSD;
   double minRMSD;
   double timeGpol;
   int numSol = 0;
   
   FILE *ifp = fopen( params->F2DockOutputFile, "rt" );

   if ( ifp == NULL ) 
     {
       printError( (char *)"Failed to open F2Dock output file ( %s )!", params->F2DockOutputFile );
       return false;
     }  
     
   FILE *ofp = fopen( params->rerankedOutputFile, "wt" );
     
   if ( ofp == NULL ) 
     {
       printError( (char *)"Failed to create reranked output file ( %s )!", params->rerankedOutputFile );
       freeClose1( );
       return false;
     }  

   char s[ 2000 ], t[ 2000 ];
   int matLoc = 0, rmsdLoc = 0, scoreLoc = 0;
                  
   while ( fgets( s, 1999, ifp ) != NULL )
     {
       if ( skipInitial( s, (char *)"# OUTPUT FORMAT:", t ) ) 
         {
           printRerankingParamters( params, ofp );
                    
           int m;
           
           if ( sscanf( t, "%d", &m ) != 1 )
             {
               printError( (char *)"Error reading F2Dock output file ( %s )!", params->F2DockOutputFile );  
               freeClose2( );
               return false;
             }
             
           fprintf( ofp, "# OUTPUT FORMAT: %d\n", m + 4 );  
           fprintf( ofp, "#	 COLNAME new_rank int rank of this docking result after GB based reranking\n" );  
           fprintf( ofp, "#	 COLNAME new_score float new score computed from delgpol and original F2Dock score\n" );             
           fprintf( ofp, "#	 COLNAME delgpol float reduction in polarization energy in kcal/mol\n" );  
           fprintf( ofp, "#	 COLNAME areaprop float approximation to constant * surfaceArea \n" );  
           
           m = 1;
           
           while ( fgets( s, 1999, ifp ) != NULL )
             {
               if ( ( m == 1 ) && !skipInitial( s, (char *)"#	 COLNAME rank int", t ) )
                 {
                   printError( (char *)"First column of F2Dock output must be the rank ( int )!" );  
                   freeClose2( );
                   return false;
                 }
                 
               if ( skipInitial( s, (char *)"#	 COLNAME mat1 ", t ) ) matLoc = m;
               else if ( skipInitial( s, (char *)"#	 COLNAME rmsd ", t ) ) rmsdLoc = m;
	       else if ( skipInitial( s, (char *)"#	 COLNAME score ", t ) ) scoreLoc = m;
               
               fprintf( ofp, "%s", s );  
               
               m++;
               if ( ( matLoc > 0 ) && ( rmsdLoc > 0 ) && ( scoreLoc > 0 ) ) break;
             }             
         }
       else if ( skipInitial( s, (char *)"# START PEAKS", t ) ) 
             {
               if ( ( matLoc <= 0 ) || ( rmsdLoc <= 0 ) )
                 {
                   printError( (char *)"Missing format data in F2Dock output file ( %s )!", params->F2DockOutputFile );  
                   freeClose2( );
                   return false;
                 }
                 
               if ( ( matLoc <= rmsdLoc ) && ( matLoc + 12 > rmsdLoc ) )
                 {
                   printError( (char *)"Invalid format data in F2Dock output file ( %s )!", params->F2DockOutputFile );  
                   freeClose2( );
                   return false;
                 }  
             
               fprintf( ofp, "%s", s );
               
               timeGpol = getTime( );
                              
               double *atomsPQR = ( double * ) malloc( 5 * ( staticMol->numAtoms + movingMol->numAtoms ) * sizeof( double ) );
               double *qPoints = ( double * ) malloc( 7 * ( staticMol->numQPoints + movingMol->numQPoints ) * sizeof( double ) );
               
               if ( ( atomsPQR == NULL ) || ( qPoints == NULL ) )
                 {
                   printError( (char *)"Failed to allocate memory!" );  
                   freeClose4( );
                   return false;                 
                 }

               Point *staticQPoints = ( Point * ) malloc( staticMol->numQPoints * sizeof( Point ) );
               Point *movingQPoints = ( Point * ) malloc( movingMol->numQPoints * sizeof( Point ) );
               
               if ( ( staticQPoints == NULL ) || ( movingQPoints == NULL ) )
                 {
                   printError( (char *)"Failed to allocate memory!" );  
                   freeClose5( );
                   return false;                 
                 }
                 
               for ( int i = 0; i < 5 * staticMol->numAtoms; i++ )
                 atomsPQR[ i ] = staticMol->atoms[ i ];  
                
               double xlatePG = computeXlateForPG( staticMol, movingMol );  

               PG *staticPG = new PG( 10.0, xlatePG, 5.0 );  

               for ( int i = 0; i < staticMol->numQPoints; i++ )
                 {
                   staticQPoints[ i ].x = staticMol->qPoints[ 7 * i + 0 ];  
                   staticQPoints[ i ].y = staticMol->qPoints[ 7 * i + 1 ];  
                   staticQPoints[ i ].z = staticMol->qPoints[ 7 * i + 2 ];  
                   
                   staticPG->addPoint( &staticQPoints[ i ] );
                 }  

               PairingHeap *pairH = new PairingHeap( true, true );
               int i = 0, solIndex;
               bool gotNumSol = false;
             
               while ( fgets( s, 1999, ifp ) != NULL )
                 {
                   if ( skipInitial( s, (char *)"# END PEAKS", t ) ) break;
                   
                   if ( !gotNumSol )
                     {
                       double d;

                       if ( !getDoublesInRange( s, 1, 1, &d ) )
                         {
                           printError( (char *)"Error reading F2Dock output file ( %s )!", params->F2DockOutputFile );  
                           freeClose6( );
                           return false;
                         }
                                                 
                       numSol = ( int ) floor( d + 0.5 );  
                       solIndex = numSol;
                       
                       if ( ( params->numSol > 0 ) && ( numSol > params->numSol ) ) numSol = params->numSol;
                       
                       solInfo = ( SOLUTION_INFO * ) malloc( ( numSol + 1 ) * sizeof( SOLUTION_INFO ) );
                       
                       if ( solInfo == NULL )
                         {
                           printError( (char *)"Failed to allocate memory!" );  
                           freeClose6( );
                           return false;
                         }    
                         
                       for ( int j = 0; j < nRange; j++ )
                           hitsInRange[ j ] = 0;
      
                       rankMin = solIndex + 1;
                       indexMinRMSD = -1;
                       numGoodPeaks = 0;                               
                     }
                     
                   if ( solIndex-- > numSol ) 
                     {
                       fprintf( ofp, "%6d %16.5lf %16.5lf %16.5lf %s", solIndex + 1, ( double ) 0.0, ( double ) 0.0, ( double ) 0.0, s );                     
                       
                       double rmsd;
                       
                       if ( !getDoublesInRange( s, rmsdLoc, rmsdLoc, &rmsd ) )
                         {
                           printError( (char *)"Error reading F2Dock output file ( %s )!", params->F2DockOutputFile );  
                           freeClose6( );
                           return false;
                         }
                       
                       indexMinRMSD = -2;
                         
                       if ( !gotNumSol || ( rmsd < minRMSD ) )
                         {
                           minRMSD = rmsd;
                           rankMinRMSD = solIndex + 1;
                         }
                       
                       if ( rmsd < rmsdGood )
                         {
                           numGoodPeaks++;
                           rankMin = solIndex;
                           for ( int j = 0; j < nRange; j++ )
                              if ( solIndex < range[ j ] ) hitsInRange[ j ]++;
                         }   
                         
                       gotNumSol = true;                 
  
                       continue;  
                     }  

                   gotNumSol = true;                 

                   printf( "i = %d\n", i );
                   fflush( stdout );
                                  
                   double trans[ 12 ], rmsd;
                   
                   if ( !getDoublesInRange( s, matLoc, matLoc + 11, trans ) || !getDoublesInRange( s, rmsdLoc, rmsdLoc, &rmsd ) )
                     {
                       printError( (char *)"Error reading F2Dock output file ( %s )!", params->F2DockOutputFile );  
                       freeClose6( );
                       return false;
                     }  

                   transformAndCopyMovingMolecule( atomsPQR, qPoints, trans, staticMol, movingMol );

                   PG *movingPG = new PG( 10.0, xlatePG, 5.0 );  
    
                   for ( int j = 0; j < movingMol->numQPoints; j++ )
                     {
                       movingQPoints[ j ].x = qPoints[ 7 * j + 0 ];  
                       movingQPoints[ j ].y = qPoints[ 7 * j + 1 ];  
                       movingQPoints[ j ].z = qPoints[ 7 * j + 2 ];  
                       
                       movingPG->addPoint( &movingQPoints[ j ] );
                     }  
                     
                   double totalQPointsWeight = 0;
                   int totalQPoints = 0;

                   for ( int i = 0; i < movingMol->numQPoints; i++ )
                     {
                       if ( !staticPG->pointsWithinRange( &movingQPoints[ i ], params->distanceCutoff ) )
                         {
                           //totalQPointsWeight += movingMol->qPoints[ 7 * i + 6 ];                             
                           
                           for ( int j = 0; j < 7; j++ )
                             qPoints[ 7 * totalQPoints + j ] = qPoints[ 7 * i + j ];  
                             
                           totalQPoints++;  
                         }  
                       else totalQPointsWeight += qPoints[ 7 * i + 6 ];                               
                     }  

                   for ( int i = 0; i < staticMol->numQPoints; i++ )
                     {
                       if ( !movingPG->pointsWithinRange( &staticQPoints[ i ], params->distanceCutoff ) )
                         {
                           //totalQPointsWeight += staticMol->qPoints[ 7 * i + 6 ];                             
                           
                           for ( int j = 0; j < 7; j++ )
                             qPoints[ 7 * totalQPoints + j ] = staticMol->qPoints[ 7 * i + j ];  
                             
                           totalQPoints++;  
                         }  
                       else totalQPointsWeight += staticMol->qPoints[ 7 * i + 6 ];                               
                     }  

//                   for ( int i = 0; i < movingMol->numQPoints; i++ )
//                     {
//                       if ( !staticPG->pointsWithinRange( &movingQPoints[ i ], params->distanceCutoff ) )
//                         {
//   //                           totalQPointsWeight += movingMol->qPoints[ 7 * i + 6 ];                             
//                           
//                           for ( int j = 0; j < 7; j++ )
//                             qPoints[ 7 * totalQPoints + j ] = movingMol->qPoints[ 7 * i + j ];  
//                             
//                           totalQPoints++;  
//                         }  
//                       else totalQPointsWeight += movingMol->qPoints[ 7 * i + 6 ];                               
//                     }  
                   
                   delete movingPG;
                                          
                   double Gpol = 0;
                   
                   if ( params->GpolWeight != 0 ) 
                      getGpol( params, staticMol->numAtoms + movingMol->numAtoms, atomsPQR, totalQPoints, qPoints, &Gpol );

                   double origScore;
                       
                   if ( !getDoublesInRange( s, scoreLoc, scoreLoc, &origScore ) )
                     {
                       printError( (char *)"Error reading F2Dock output file ( %s )!", params->F2DockOutputFile );  
                       freeClose6( );
                       return false;
                     }
                     
                   solInfo[ i ].sol = strdup( s );
                   solInfo[ i ].F2DockRank = solIndex + 1;                   
                   solInfo[ i ].rmsd = rmsd;
                   solInfo[ i ].DelGpol = Gpol - ( staticMol->Gpol + movingMol->Gpol );
                   solInfo[ i ].areaProp = totalQPointsWeight; 
                   solInfo[ i ].newScore = params->GnonpolWeight * solInfo[ i ].areaProp + params->GpolWeight * solInfo[ i ].DelGpol + params->F2DockScoreWeight * origScore;                   
                   
                   pairH->Insert( i, solInfo[ i ].newScore );
                   
                   i++;
                 }

               delete staticPG;
               
               freeMem( atomsPQR );
               freeMem( qPoints );                 
               
               freeMem( staticQPoints );
               freeMem( movingQPoints );               
               
               int curIndex = -1, prevIndex, nextIndex;

               while ( i-- )
                 {
                   int k;
                   double d;
                   
                   pairH->Delete_Min( k, d );                   
                   
                   solInfo[ k ].nextIndex = curIndex;
                   curIndex = k;                   
                 }
                 
               solInfo[ numSol ].nextIndex = curIndex;  
                                               
               delete pairH;                 
               
               i = numSol;               
               curIndex = numSol;
               for ( int j = 0, k = 0; j < params->numBands; j++ )
                 {
                   prevIndex = curIndex;
                   nextIndex = solInfo[ curIndex ].nextIndex;
                   int maxK = ( ( params->bands[ 2 * j ] > numSol ) ? numSol : params->bands[ 2 * j ] );
                   
//                   printf( "\n j = %d, params->bands[ 2 * j ] = %d, params->bands[ 2 * j + 1 ] = %d, maxK = %d\n", 
//                               j, params->bands[ 2 * j ], params->bands[ 2 * j + 1 ], maxK );
//                   fflush( stdout );
//                    
//                   printf( "curIndex = %d, prevIndex = %d, nextIndex = %d\n", curIndex, prevIndex, nextIndex ); 
//                   fflush( stdout );
                    
                   while ( k < maxK )
                     {
                       while ( solInfo[ nextIndex ].F2DockRank > params->bands[ 2 * j + 1 ] )
                         {
                           prevIndex = nextIndex;
                           nextIndex = solInfo[ nextIndex ].nextIndex;
//                           printf( "prevIndex = %d, nextIndex = %d\n", prevIndex, nextIndex );                           
//                           fflush( stdout );
                         }  
                       
                       int tempIndex = solInfo[ curIndex ].nextIndex;
                       
                       solInfo[ prevIndex ].nextIndex = solInfo[ nextIndex ].nextIndex;
                       solInfo[ nextIndex ].nextIndex = solInfo[ curIndex ].nextIndex;
                       solInfo[ curIndex ].nextIndex = nextIndex;
                       curIndex = nextIndex;
                       prevIndex = solInfo[ prevIndex ].nextIndex;
                       nextIndex = solInfo[ prevIndex ].nextIndex;

//                       printf( "k = %d, curIndex = %d, prevIndex = %d, nextIndex = %d\n", k, curIndex, prevIndex, nextIndex ); 
//                       fflush( stdout );
                                                   
                       k++;
                     }                   
                     
//                   printf( "(outside) curIndex = %d, prevIndex = %d, nextIndex = %d\n", curIndex, prevIndex, nextIndex ); 
//                   fflush( stdout );                     
                 }
                                  
               curIndex = solInfo[ numSol ].nextIndex;           
               prevIndex = numSol;
                     
               while ( curIndex != -1 )
                 {
                   nextIndex = solInfo[ curIndex ].nextIndex;
                   solInfo[ curIndex ].nextIndex = prevIndex;
                   prevIndex = curIndex;
                   curIndex = nextIndex;
                 }               
                 
               i = numSol;  
               curIndex = prevIndex;
               while ( curIndex != numSol )
                 {
                   int k = curIndex;
                   
                   fprintf( ofp, "%6d %16.5lf %16.5lf %16.5lf %s", i, solInfo[ k ].newScore, solInfo[ k ].DelGpol, solInfo[ k ].areaProp, solInfo[ k ].sol );
                  
                   if ( ( indexMinRMSD == -1 ) || ( solInfo[ k ].rmsd < minRMSD ) )
                     {
                       indexMinRMSD = k;
                       minRMSD = solInfo[ k ].rmsd;
                       rankMinRMSD = i;
                     }
                   
                   if ( solInfo[ k ].rmsd < rmsdGood )
                     {
                       numGoodPeaks++;
                       rankMin = i - 1;
                       for ( int j = 0; j < nRange; j++ )
                          if ( i <= range[ j ] ) hitsInRange[ j ]++;
                     }                 
                 
                   curIndex = solInfo[ curIndex ].nextIndex;
                   i--;
                 }
                 
                                             
/*               while ( i-- )
                 {
                   int k;
                   double d;
                   
                   pairH->Delete_Min( k, d );                   
                   fprintf( ofp, "%6d %16.5lf %16.5lf %16.5lf %s", i + 1, solInfo[ k ].newScore, solInfo[ k ].DelGpol, solInfo[ k ].areaProp, solInfo[ k ].sol );
                  
                   if ( ( indexMinRMSD == -1 ) || ( solInfo[ k ].rmsd < minRMSD ) )
                     {
                       indexMinRMSD = k;
                       minRMSD = solInfo[ k ].rmsd;
                       rankMinRMSD = i + 1;
                     }
                   
                   if ( solInfo[ k ].rmsd < rmsdGood )
                     {
                       numGoodPeaks++;
                       rankMin = i;
                       for ( int j = 0; j < nRange; j++ )
                          if ( i < range[ j ] ) hitsInRange[ j ]++;
                     }
                 }
*/
                 
               fprintf( ofp, "# END PEAKS\n" );  
                                             
               timeGpol = getTime( ) - timeGpol;
             }
             
      else if ( skipInitial( s, (char *)"#Hits in Range:", t ) ) 
             {
               f_printf( ofp, (char *)"%s", s );
               for ( int i = 0; i < nRange; i++ )
		 f_printf( ofp, (char *)"# [%d, %d] --> %d\n", 1, range[ i ], hitsInRange[ i ] );

               while ( fgets( s, 1999, ifp ) != NULL )
                 {
                   if ( !skipInitial( s, (char *)"# [", t ) ) break;
                 }
                 
               f_printf( ofp, (char *)"%s", s );                 
             }
      else if ( skipInitial( s, (char *)"# good peaks under ", t ) ) 
             {           
               f_printf( ofp, (char *)"# good peaks under %.0lf A: count = %d highest rank = %d min RMSD = %lf\n", rmsdGood, numGoodPeaks, rankMin, minRMSD );  
             }  
      else if ( skipInitial( s, (char *)"# best peak: rmsd = ", t ) ) 
             {
                double d;
                int i = 0;
                
                i = getDouble( t, i, &d );  // rmsd value
                i = getString( t, i, s );  // "rank"
                i = getString( t, i, s );  // "="
                i = getDouble( t, i, &d );  // rank value
                
                f_printf( ofp, (char *)"# best peak: rmsd = %lf rank = %d delgpol = %lf %s", minRMSD, rankMinRMSD, 
                         ( indexMinRMSD < 0 ) ? ( ( double ) 0.0 ) : solInfo[ indexMinRMSD ].DelGpol, t + i );                  
             }             
      else if ( skipInitial( s, (char *)"# total time = ", t ) )       
             {
                double d;
                int i = 0;
                
                i = getDouble( t, i, &d );  // time in seconds

                fprintf( ofp, "# reranking time = %lf sec\n#\n", timeGpol );

                fprintf( ofp, "# total time = %lf sec\n", d + timeGpol );
             }
      else fprintf( ofp, "%s", s );       
     }
          
   for ( int i = 0; i < numSol; i++ ) 
      freeMem( solInfo[ i ].sol ); 
      
   freeMem( solInfo ); 
          
   fclose( ifp );
   fclose( ofp );  
     
   return true;  
}


#define returnSpectrumError1( ) { printError( (char *)"Invalid spectrum ( %s )!", spectrum ); return false; }
#define returnSpectrumError2( ) { freeMem( p->bands ); returnSpectrumError1( ); }

bool decodeSpectrum( PARAMS_IN *p, char *spectrum )
{
   if ( spectrum == NULL ) return true;
   
   int numBands = 0;
   int i = 0;
  
   while ( spectrum[ i ] )
     {
       if ( spectrum[ i ] == ':' ) numBands++;
       i++;
     }
     
   p->bands = ( int * ) malloc( 2 * ( numBands + 1 ) * sizeof( int ) ); 
   
   if ( p->bands == NULL )  returnSpectrumError1( );
     
   int k = 0;
   i = 0;
   numBands = 0;
   
   while ( spectrum[ i ] )  
     {
       i = skipWhiteSpaces( spectrum, i );
       
       if ( !spectrum[ i ] ) break;
       
       int v1, v2;
       
       int j = getInt( spectrum, i, &v1 );
       
       if ( ( i == j ) || ( v1 < 0 ) ) returnSpectrumError2( );
         
       i = skipWhiteSpaces( spectrum, j );  
       
       if ( spectrum[ i ] != ':' ) returnSpectrumError2( );
       
       i = skipWhiteSpaces( spectrum, i + 1 );
       
       if ( !spectrum[ i ] ) returnSpectrumError2( );
       
       j = getInt( spectrum, i, &v2 );
       
       if ( ( i == j ) || ( v2 < 0 ) ) returnSpectrumError2( );       
       
       if ( !k )
         {
           p->bands[ k++ ] = v1;
           p->bands[ k++ ] = v2;
           numBands = 1;
         }
       else
         {  
           if ( v2 < p->bands[ k - 1 ] ) returnSpectrumError2( );
           
           if ( v2 > p->bands[ k - 1 ] )
             {
               p->bands[ k ] = p->bands[ k - 2 ] + v1;
               p->bands[ k + 1 ] = v2;
               if ( p->bands[ k ] > p->bands[ k + 1 ] ) returnSpectrumError2( );               
               k += 2;
               numBands++;
             }
           else 
             { 
               p->bands[ k - 2 ] += v1;
               if ( p->bands[ k - 2 ] > p->bands[ k - 1 ] ) returnSpectrumError2( );
             }  
         }
                     
       i = skipWhiteSpaces( spectrum, j );         
       
       if ( spectrum[ i ] == '-' ) i++;
     }
     
   p->numBands = numBands;  
     
   return true;  
}


bool getParamsFromFile( PARAMS_IN *p, char *paramFile )
{
  char s[ 2000 ];
  char key[ 500 ], val[ 500 ];
  char *spectrum;
  FILE *fp;

  fp = fopen( paramFile, "r" );

  if ( fp == NULL )
    {
      printError( (char *)"Failed to open parameter file %s!", paramFile );
      return false;
    }

  p->staticMoleculePQR = p->movingMoleculePQR = NULL;
  p->staticMoleculeQUAD = p->staticMoleculeQUAD = NULL;
  p->F2DockOutputFile = p->rerankedOutputFile = NULL;
  p->F2DockScoreWeight = 0.0;
  p->GpolWeight = -1.0;
  p->GnonpolWeight = 0.005;  
  p->distanceCutoff = 1.5;
  p->numThreadsBR = p->numThreadsGpol = 4;
  p->epsilonBR = 0.5;
  p->epsilonGpol = 0.7;
  p->useApproxMath = false;
  p->numSol = 100;
  p->numBands = 0;
  p->bands = NULL;
  
  spectrum = NULL;

  while ( fgets( s, 1999, fp ) != NULL )
    {
      if ( sscanf( s, "%s %s", key, val ) != 2 ) continue;
    
      if ( !strcasecmp( key, "staticMoleculePQR" ) ) p->staticMoleculePQR = strdup( val );
      else if ( !strcasecmp( key, "movingMoleculePQR" ) ) p->movingMoleculePQR = strdup( val );
      else if ( !strcasecmp( key, "staticMoleculeQUAD" ) ) p->staticMoleculeQUAD = strdup( val );
      else if ( !strcasecmp( key, "movingMoleculeQUAD" ) ) p->movingMoleculeQUAD = strdup( val );
      else if ( !strcasecmp( key, "F2DockOutputFile" ) ) p->F2DockOutputFile = strdup( val );
      else if ( !strcasecmp( key, "rerankedOutputFile" ) ) p->rerankedOutputFile = strdup( val );
      else if ( !strcasecmp( key, "F2DockScoreWeight" ) ) p->F2DockScoreWeight = atof( val );
      else if ( !strcasecmp( key, "GpolWeight" ) ) p->GpolWeight = atof( val );      
      else if ( !strcasecmp( key, "GnonpolWeight" ) ) p->GnonpolWeight = atof( val );      
      else if ( !strcasecmp( key, "distanceCutoff" ) )
             {
               double v = atof( val );
               
               if ( v < 0 )
                 {
                   printError( (char *)"%s must be a non-negative float!", key );
                   fclose( fp );
                   return false;
                 }
                 
               p->distanceCutoff = v;                 
             }                  
      else if ( !strcasecmp( key, "numThreadsBR" ) )
             {
               int v = atoi( val );
               
               if ( v < 1 )
                 {
                   printError( (char *)"%s must be a positive integer!", key );
                   fclose( fp );
                   return false;
                 }
                 
               p->numThreadsBR = v;  
             }
      else if ( !strcasecmp( key, "numThreadsGpol" ) )
             {
               int v = atoi( val );
               
               if ( v < 1 )
                 {
                   printError( (char *)"%s must be a positive integer!", key );
                   fclose( fp );
                   return false;
                 }
                 
               p->numThreadsGpol = v;                 
             }
      else if ( !strcasecmp( key, "epsilonBR" ) )
             {
               double v = atof( val );
               
               if ( v < 0 )
                 {
                   printError( (char *)"%s must be a non-negative float!", key );
                   fclose( fp );
                   return false;
                 }
                 
               p->epsilonBR = v;                 
             }
      else if ( !strcasecmp( key, "epsilonGpol" ) )
             {
               double v = atof( val );
               
               if ( v < 0 )
                 {
                   printError( (char *)"%s must be a non-negative float!", key );
                   fclose( fp );
                   return false;
                 }
                 
               p->epsilonGpol = v;                 
             }                    
      else if ( !strcasecmp( key, "useApproxMath" ) ) 
             {
	       if ( !strcasecmp( val, "true" ) ) p->useApproxMath = true;
	       else if ( !strcasecmp( val, "false" ) ) p->useApproxMath = false;
	            else  
	              {
		        printError( (char *)"%s must be a Boolean value!", key );
		        fclose( fp );
		        return false;
	              }	    
	     }             
      else if ( !strcasecmp( key, "numSol" ) )
             {
               int v = atoi( val );
               
               if ( v < 1 )
                 {
                   printError( (char *)"%s must be a positive integer!", key );
                   fclose( fp );
                   return false;
                 }
                 
               p->numSol = v;                 
             }            
      else if ( !strcasecmp( key, "spectrum" ) ) spectrum = strdup( val );
    }

  fclose( fp );
    
  if ( p->staticMoleculePQR == NULL )  
    { 
      printError( (char *)"Missing PQR file name for the static molecule!" );
      return false;
    }

  if ( p->movingMoleculePQR == NULL )  
    { 
      printError( (char *)"Missing PQR file name for the moving molecule!" );
      return false;
    }

  if ( p->staticMoleculeQUAD == NULL )  
    { 
      printError( (char *)"Missing QUAD file name for the static molecule!" );
      return false;
    }

  if ( p->movingMoleculeQUAD == NULL )  
    { 
      printError( (char *)"Missing QUAD file name for the moving molecule!" );
      return false;
    }
    
  if ( p->F2DockOutputFile == NULL )  
    { 
      printError( (char *)"Missing F2Dock output file name!" );
      return false;
    }

  if ( p->rerankedOutputFile == NULL )  
    { 
      printError( (char *)"Missing file name for reranked output!" );
      return false;
    }
        
  if ( !decodeSpectrum( p, spectrum ) ) return false;    

  freeMem( spectrum );
        
  return true;
}



void getMoleculeInformation( PARAMS_IN *params, MOLECULE_INFO *staticMol, MOLECULE_INFO *movingMol )
{
  getAtomsQuadsAndGpol( params, params->staticMoleculePQR, params->staticMoleculeQUAD, 
                        &( staticMol->numAtoms ), &( staticMol->atoms ),
                        &( staticMol->numQPoints ), &( staticMol->qPoints ), 
                        &( staticMol->Gpol ) );  
                        
  getAtomsQuadsAndGpol( params, params->movingMoleculePQR, params->movingMoleculeQUAD, 
                        &( movingMol->numAtoms ), &( movingMol->atoms ),
                        &( movingMol->numQPoints ), &( movingMol->qPoints ), 
                        &( movingMol->Gpol ) );                                                                        
}


int main( int argc, char *argv[ ] )
{
  if ( argc < 2 )
    {
      printError( (char *)"Input text file not specified!" );
      return 1;
    }
     
  PARAMS_IN params;
   
  if ( !getParamsFromFile( &params, argv[ 1 ] ) ) return 1; 
  
  MOLECULE_INFO staticMol, movingMol;
   
  getMoleculeInformation( &params, &staticMol, &movingMol );   

  if ( !rerankF2DockOutput( &params, &staticMol, &movingMol ) ) return 1;
    
  return 0;
}
