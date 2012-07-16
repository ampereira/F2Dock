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

#include "fastGpol.h"

using fastGB::fastGpol;


bool fastGpol::setMinRadius( double minRadius )
{
   if ( minRadius < 0 )
     {
      printError( (char *)"minRadius must be a non-negative real number!" );
      return false;     
     }
     
   minRad = minRadius;

   if ( printStatus ) printf( (char *)"\nminRadius is set to %lf\n", minRadius );
   
   return true;
}


bool fastGpol::setMaxLeafSize( int maxLfSize )
{
   if ( maxLfSize <= 0 )
     {
      printError( (char *)"maxLeafSize must be a positive integer!" );
      return false;     
     }
     
   maxLeafAtoms = maxLfSize;

   if ( printStatus ) printf( (char *)"\nmaxLeafSize is set to %d\n", maxLfSize );

   return true;
}


bool fastGpol::setEpsilon( double eps )
{
   if ( eps < 0 )
     {
      printError( (char *)"epsilon must be a non-negative real number!" );
      return false;     
     }
     
   epsilon = eps;

   if ( printStatus ) printf( (char *)"\nepsilon is set to %lf\n", eps );

   return true;
}


bool fastGpol::setNumThreads( int nThreads )
{
   if ( nThreads < 1 )
     {
      printError( (char *)"numThreads must be a positive integer!" );
      return false;     
     }

   numThreads = nThreads;

   if ( printStatus ) printf( (char *)"\nnumThreads is set to %d\n", nThreads );
   
   return true;
}


bool fastGpol::setIonDielectric( double e_ion )
{
   if ( e_ion < 0 )
     {
       printError( (char *)"Ion dielectric must be a positive float!" );
       return false;     
     }

   e_in = e_ion;

   if ( printStatus ) printf( (char *)"\nIon dielectric is set to %lf\n", e_ion );
   
   return true;
}


bool fastGpol::setSolventDielectric( double e_sol )
{
   if ( e_sol < 0 )
     {
       printError( (char *)"Solvent dielectric must be a positive float!" );
       return false;     
     }

   e_out = e_sol;

   if ( printStatus ) printf( (char *)"\nSolvent dielectric is set to %lf\n", e_sol );
   
   return true;
}


void fastGpol::useApproxMathFunctions( bool approxMath )
{
   useApproxMath = approxMath;
   
   if ( printStatus ) 
     {
       if ( useApproxMath ) printf( (char *)"\nApproximate math functions will be used\n" );
       else printf( (char *)"\nExact math functions will be used\n" );
     }  
}


void fastGpol::setPrintStatus( bool printStat )
{
   printStatus = printStat;
   
   if ( printStatus ) printf( (char *)"\nprintStatus is set to true\n" );
}


void fastGpol::printCurrentSettings( void )
{
   printf( (char *)"\nCurrent Parameter Settings:\n" );
   printf( (char *)"\n\tminRadius = %lf, maxLeafSize = %d, epsilon = %lf, numThreads = %d, numAtoms = %d\n", 
              minRad, maxLeafAtoms, epsilon, numThreads, numAtoms );
   printf( (char *)"\tDielectric constants: ion = %lf, solvent = %lf\n", e_in, e_out );              

   if ( useApproxMath ) printf( (char *)"\tUse approximate math functions: YES\n" );
   else printf( (char *)"\tUse approximate math functions: NO\n" );

#ifdef USE_SSE
   printf( (char *)"\tUse intrinsic SSE functions: YES\n\n" );
#else   
   printf( (char *)"\tUse intrinsic SSE functions: NO\n\n" );
#endif
}


void fastGpol::setDefaults( void )
{
   numThreads = 8;
   epsilon = 0.3;
   minRad = 2.0;
   maxLeafAtoms = 10;
   e_in = 1;
   e_out = 80.0;
   useApproxMath = false;
   printStatus = true;
}


void fastGpol::freeMemory( void )
{
   freeMem( atoms );   
   freeMem( octree );   
   freeMem( approxR );
   freeMem( qSum );
}


fastGpol::fastGpol( char *pqrRFile, bool printStat )
{
   setDefaults( );
   
   printStatus = printStat;
   
   if ( printStatus ) printCurrentSettings( );              
      
   if ( !readFromPQRRFile( pqrRFile ) )
      {
       freeMemory( );
       exit( 1 );
      }
      
   if ( !initDataStructures( ) )
      {
       freeMemory( );
       exit( 1 );
      }              
}


fastGpol::fastGpol( char *pqrRFile )
{
   fastGpol( pqrRFile, true );
}


fastGpol::fastGpol( int nAtoms, double *pqrR, bool printStat )
{
   setDefaults( );   
   
   printStatus = printStat;
   
   if ( printStatus ) printCurrentSettings( );              
   
   if ( !readFromPQRRArray( nAtoms, pqrR ) )
      {
       freeMemory( );
       exit( 1 );
      }

   if ( !initDataStructures( ) )
      {
       freeMemory( );
       exit( 1 );
      }                            
}


fastGpol::fastGpol( int nAtoms, double *pqrR )
{
   fastGpol( nAtoms, pqrR, true );
}


fastGpol::~fastGpol( )
{
   freeMemory( );
   if ( printStatus ) printf( (char *)"\n" );
}


bool fastGpol::readFromPQRRFile( char *fname )
{
   FILE *fp;

   if ( printStatus ) printf( (char *)"\nreading from PQRR file %s... ", fname );   

   double startT = getTime( );
   
   fp = fopen( fname, "rt" );
   
   if ( fp == NULL )
     {
      printError( (char *)"Failed to open input file %s!", fname );
      return false;
     }
     
   if ( fscanf( fp, (char *)"%d", &numAtoms ) != 1 )
     {
      printError( (char *)"Failed to read input file %s!", fname );
      fclose( fp );
      return false;
     }      

   if ( numAtoms <= 0 )
     {
      printError( (char *)"Invalid number of atoms (%d)!", numAtoms );
      fclose( fp );
      return false;
     }      
          
   atoms = ( ATOM * ) malloc( ( numAtoms ) * sizeof( ATOM ) );
   
   if ( atoms == NULL )
     {
      printError( (char *)"Memory allocation failed!" );
      fclose( fp );
      return false;
     }

   int i;
     
   for ( i = 0; i < numAtoms; i++ )
     {
      double x, y, z, r, q, R;
      
      if ( fscanf( fp, (char *)"%lf %lf %lf %lf %lf %lf", &x, &y, &z, &q, &r, &R ) != 6 ) break;
      
      atoms[ i ].x = x;
      atoms[ i ].y = y;
      atoms[ i ].z = z;
      atoms[ i ].r = r;
      atoms[ i ].q = q;
      atoms[ i ].R = R;      
     }  

   if ( i < numAtoms )
     {
      printError( (char *)"Not enough data in input file %s (expected %d, found %d)!", fname, numAtoms, i );
      numAtoms = 0;
      fclose( fp );
      return false;      
     }

   fclose( fp );

   double endT = getTime( );

   if ( printStatus ) printf( (char *)"done ( %lf sec, read %d atoms )\n", endT - startT, numAtoms );
     
   return true;  
}


bool fastGpol::readFromPQRRArray( int nAtoms, double *pqrR )
{   
   if ( printStatus ) printf( (char *)"\ncopying data from input array... " );   

   double startT = getTime( );

   if ( nAtoms <= 0 )
     {
      printError( (char *)"Invalid number of atoms (%d)!", nAtoms );
      return false;
     }      
     
   numAtoms = nAtoms;     
          
   atoms = ( ATOM * ) malloc( numAtoms * sizeof( ATOM ) );
   
   if ( atoms == NULL )
     {
      printError( (char *)"Memory allocation failed!" );
      return false;
     }

   for ( int i = 0; i < numAtoms; i++ )
     {
       atoms[ i ].x = pqrR[ 6 * i + 0 ];
       atoms[ i ].y = pqrR[ 6 * i + 1 ];
       atoms[ i ].z = pqrR[ 6 * i + 2 ];
       atoms[ i ].q = pqrR[ 6 * i + 3 ];
       atoms[ i ].r = pqrR[ 6 * i + 4 ];
       atoms[ i ].R = pqrR[ 6 * i + 5 ];      
     }  

   double endT = getTime( );
     
   if ( printStatus ) printf( (char *)"done ( %lf sec, copied %d atoms )\n", endT - startT, numAtoms );
     
   return true;  
}


bool fastGpol::initDataStructures( void )
{
   if ( printStatus ) printf( (char *)"\ninitializing data structures... " );   

   double startT = getTime( );
   
   double minR = atoms[ 0 ].R;
   for ( int i = 1; i < numAtoms; i++ )
      if ( atoms[ i ].R < minR ) minR = atoms[ i ].R;
           
   maxGroupR = 0;
   for ( int i = 0; i < numAtoms; i++ )
     {
      atoms[ i ].groupR = ( int ) floor( log( atoms[ i ].R / minR ) / log( 1 + epsilon ) );
      if ( atoms[ i ].groupR > maxGroupR ) maxGroupR = atoms[ i ].groupR;
     }
      
   approxR = ( double * ) malloc( ( maxGroupR + 1 ) * sizeof( double ) );
   if ( approxR == NULL )
     {
      printError( (char *)"Memory allocation failed!" );
      return false;
     }
     
   approxR[ 0 ] = minR;  
   for ( int i = 1; i <= maxGroupR; i++ )
     approxR[ i ] = ( 1 + epsilon ) * approxR[ i - 1 ];
   
   ATOM *atomsT = ( ATOM * ) malloc( numAtoms * sizeof( ATOM ) );

   if ( atomsT == NULL )
     {
      printError( (char *)"Memory allocation failed!" );
      return false;
     }
     
   int numOctreeNodes = 0;

   countOctreeNodesAndSortAtoms( atoms, 0, numAtoms - 1, minRad * minRad, maxLeafAtoms, atomsT, &numOctreeNodes );   

   free( atomsT );
      
   octree = ( OCTREE_NODE * ) malloc( numOctreeNodes * sizeof( OCTREE_NODE ) );
   qSum = ( double * ) malloc( ( maxGroupR + 1 ) * numOctreeNodes * sizeof( double ) );

   if ( ( octree == NULL ) || ( qSum == NULL ) )
     {
      printError( (char *)"Memory allocation failed!" );
      return false;
     }

   octree[ 0 ].atomStartID = 0;   
   octree[ 0 ].atomEndID = numAtoms - 1;      
   
   int freeID = 1;

   constructOctree( octree, 0, qSum, maxGroupR + 1, minRad, maxLeafAtoms, atoms, &freeID );

   double endT = getTime( );

   if ( printStatus ) printf( (char *)"done ( %lf sec )\n", endT - startT );
    
   return true;  
}


void fastGpol::computeQuadGpol( double *Gpol )
{
   if ( printStatus ) printf( (char *)"\ncomputing G_pol naively... " );

   double startT = getTime( );

   *Gpol = 0;

   for ( int i = 0; i < numAtoms; i++ )
      for ( int j = 0; j < numAtoms; j++ )
//         if ( i != j )
           {
            double rij2 = ( atoms[ i ].x - atoms[ j ].x ) * ( atoms[ i ].x - atoms[ j ].x )
                        + ( atoms[ i ].y - atoms[ j ].y ) * ( atoms[ i ].y - atoms[ j ].y )
                        + ( atoms[ i ].z - atoms[ j ].z ) * ( atoms[ i ].z - atoms[ j ].z );
            double RiRj = atoms[ i ].R * atoms[ j ].R;
            double qiqj = atoms[ i ].q * atoms[ j ].q;            
                                                  
            *Gpol += ( qiqj / sqrt( rij2 + RiRj * exp( - rij2 / ( 4 * RiRj ) ) ) );            
           }

   double tau = 1 / e_in - 1 / e_out;
   *Gpol = - 332 * 0.5 * tau * ( *Gpol );  
              
   double endT = getTime( );

   if ( printStatus ) printf( (char *)"done ( %lf sec )\n", endT - startT );           
}



int fastGpol::nextJob( int &nodeI, int &nodeJ )
{
   int jobExists = 0;
   
   pthread_mutex_lock( &jobLock );
   if ( jobs.size( ) > 0 ) 
     {
       JOB job = jobs.back( );
       
       jobs.pop_back( );
       
       nodeI = job.nodeI;
       nodeJ = job.nodeJ;       
       
       jobExists = 1;
     }  
   pthread_mutex_unlock( &jobLock );
   
   return jobExists;
}


void fastGpol::initJobServer( std::vector< JOB > &leafJobs, std::vector< JOB > &nonLeafJobs )
{
   pthread_mutex_init( &jobLock, NULL );
  
   jobs.clear( );
  
   int k = leafJobs.size( );
   while ( k-- ) 
     {
       jobs.push_back( leafJobs.back( ) );
       leafJobs.pop_back( );
     }  
   leafJobs.clear( );

   k = nonLeafJobs.size( );
   while ( k-- ) 
     {
       jobs.push_back( nonLeafJobs.back( ) );    
       nonLeafJobs.pop_back( );
     }  
   nonLeafJobs.clear( );  
}


inline float fastGpol::invSqrt( float x )
{
   if ( x < 0 ) return INFINITY;
   
   volatile UL_F_union xx;
   
   xx.ix = 0;      
   xx.fx = x;
   xx.ix = ( 0xBE6EB50CUL - xx.ix ) >> 1;   
      
   return ( 0.5f * xx.fx * ( 3.0f - x * xx.fx * xx.fx ) );
}

inline float fastGpol::fastExp( float x )
{
   if ( x > 85 ) return INFINITY;
   if ( x < -85 ) return 0.0;   
   
   volatile UL_F_union xx;
         
   // 12102203.16156148555068672305845f = ( 2^23 ) / ln(2);      
   unsigned long i = ( unsigned long ) ( x * 12102203.16156f );
   
   // 361007 = ( 0.08607133 / 2 ) * ( 2^23 )         
   xx.ix = i + ( 127L << 23 ) - 361007;
   
   return xx.fx;
}

#ifdef USE_SSE

inline float fastGpol::vectorGB( v4sf xi, v4sf yi, v4sf zi, v4sf Ri, v4sf qi,
                                 v4sf xj, v4sf yj, v4sf zj, v4sf Rj, v4sf qj )
{
   v4sf rij2, RiRj, RiRjx4, qiqj, val;
   V4SF val2;
 
   xi = _mm_sub_ps( xi, xj );                       
   rij2 = _mm_mul_ps( xi, xi );
   
   yi = _mm_sub_ps( yi, yj );
   yi = _mm_mul_ps( yi, yi );
   rij2 = _mm_add_ps( rij2, yi );
   
   zi = _mm_sub_ps( zi, zj );
   zi = _mm_mul_ps( zi, zi );                       
   rij2 = _mm_add_ps( rij2, zi );
   
   RiRj = _mm_mul_ps( Ri, Rj );
   qiqj = _mm_mul_ps( qi, qj );
   
   val = _mm_sub_ps( ZERO, rij2 );
   RiRjx4 = _mm_mul_ps( RiRj, FOUR );
   val = _mm_div_ps( val, RiRjx4 );
   
   val = exp_ps( val );
   val = _mm_mul_ps( RiRj, val );
   val = _mm_add_ps( rij2, val );
   val = _mm_sqrt_ps( val );
   val = _mm_div_ps( qiqj, val );                       
  
   val2.v = val;
   
   return ( val2.f[ 0 ] + val2.f[ 1 ] + val2.f[ 2 ] + val2.f[ 3 ] );
}


inline float fastGpol::vectorGB( float rij2f, v4sf Ri, v4sf qi, v4sf Rj, v4sf qj )
{
   v4sf rij2, RiRj, RiRjx4, qiqj, val;
   V4SF val2;
    
   rij2 = _mm_set1_ps( rij2f );  
    
   RiRj = _mm_mul_ps( Ri, Rj );
   qiqj = _mm_mul_ps( qi, qj );
   
   val = _mm_sub_ps( ZERO, rij2 );
   RiRjx4 = _mm_mul_ps( RiRj, FOUR );
   val = _mm_div_ps( val, RiRjx4 );
   
   val = exp_ps( val );
   val = _mm_mul_ps( RiRj, val );
   val = _mm_add_ps( rij2, val );
   val = _mm_sqrt_ps( val );
   val = _mm_div_ps( qiqj, val );                       
  
   val2.v = val;
   
   return ( val2.f[ 0 ] + val2.f[ 1 ] + val2.f[ 2 ] + val2.f[ 3 ] );
}

#endif



void fastGpol::countOctreeNodesAndSortAtoms( ATOM *atom, int atomStartID, int atomEndID, 
                                             double minRad2, int maxLeafAtoms,
                                             ATOM *atomT, int *numNodes )
{
   double minX = atom[ atomStartID ].x, minY = atom[ atomStartID ].y, minZ = atom[ atomStartID ].z;
   double maxX = atom[ atomStartID ].x, maxY = atom[ atomStartID ].y, maxZ = atom[ atomStartID ].z;
   
   for ( int i = atomStartID + 1; i <= atomEndID; i++ )
     {
      if ( atom[ i ].x < minX ) minX = atom[ i ].x;      
      if ( atom[ i ].x > maxX ) maxX = atom[ i ].x;      
      
      if ( atom[ i ].y < minY ) minY = atom[ i ].y;      
      if ( atom[ i ].y > maxY ) maxY = atom[ i ].y;      

      if ( atom[ i ].z < minZ ) minZ = atom[ i ].z;      
      if ( atom[ i ].z > maxZ ) maxZ = atom[ i ].z;      
     } 
   
   double centerX = ( minX + maxX ) / 2,
          centerY = ( minY + maxY ) / 2,
          centerZ = ( minZ + maxZ ) / 2;

   double r2 = ( atom[ atomStartID ].x - centerX ) * ( atom[ atomStartID ].x - centerX )
             + ( atom[ atomStartID ].y - centerY ) * ( atom[ atomStartID ].y - centerY )
             + ( atom[ atomStartID ].z - centerZ ) * ( atom[ atomStartID ].z - centerZ );
   
   for ( int i = atomStartID + 1; i <= atomEndID; i++ )
     {
      double r2T = ( atom[ i ].x - centerX ) * ( atom[ i ].x - centerX )
                 + ( atom[ i ].y - centerY ) * ( atom[ i ].y - centerY )
                 + ( atom[ i ].z - centerZ ) * ( atom[ i ].z - centerZ );
     
      if ( r2T > r2 ) r2 = r2T;
     } 
   
   *numNodes = 1;
      
   if ( ( atomEndID - atomStartID +  1 > maxLeafAtoms ) && ( r2 > minRad2 ) )
     {
      int atomCount[ 8 ] = { 0, 0, 0, 0, 0, 0, 0, 0 };
      
      for ( int i = atomStartID; i <= atomEndID; i++ )
        {
         atomT[ i ] = atom[ i ];
         
         int j = ( zeroIfLess( atom[ i ].z, centerZ ) << 2 )
               + ( zeroIfLess( atom[ i ].y, centerY ) << 1 )
               + ( zeroIfLess( atom[ i ].x, centerX ) );
         
         atomCount[ j ]++;
        }
        
      int atomStartIndex[ 8 ];
      int atomCurIndex[ 8 ];              
      
      atomCurIndex[ 0 ] = atomStartIndex[ 0 ] = atomStartID;
      for ( int i = 1; i < 8; i++ )
        atomCurIndex[ i ] = atomStartIndex[ i ] = atomStartIndex[ i - 1 ] + atomCount[ i - 1 ];

      for ( int i = atomStartID; i <= atomEndID; i++ )
        {        
         int j = ( zeroIfLess( atomT[ i ].z, centerZ ) << 2 )
               + ( zeroIfLess( atomT[ i ].y, centerY ) << 1 )
               + ( zeroIfLess( atomT[ i ].x, centerX ) );
           
         atom[ atomCurIndex[ j ] ] = atomT[ i ];
         atomCurIndex[ j ]++;  
        }        
        
      for ( int i = 0; i < 8; i++ ) 
        if ( atomCount[ i ] > 0 )
          {
           int numNodesT = 0;
           
           countOctreeNodesAndSortAtoms( atom, atomStartIndex[ i ], atomStartIndex[ i ] + atomCount[ i ] - 1, minRad2, maxLeafAtoms, atomT, &numNodesT );
         
           *numNodes += numNodesT;
          }
     }  
}



void fastGpol::constructOctree( OCTREE_NODE *octree, int nodeID, double *qSum, 
                                int numGroupR, double minRad, int maxLeafAtoms, 
                                ATOM *atom, int *freeID )
{
   int atomStartID = octree[ nodeID ].atomStartID,
       atomEndID = octree[ nodeID ].atomEndID;   

   int qPtr = numGroupR * nodeID;
   for ( int i = qPtr; i < qPtr + numGroupR; i++ )
     qSum[ i ] = 0;

   double minX = atom[ atomStartID ].x, minY = atom[ atomStartID ].y, minZ = atom[ atomStartID ].z;
   double maxX = atom[ atomStartID ].x, maxY = atom[ atomStartID ].y, maxZ = atom[ atomStartID ].z;
   
   for ( int i = atomStartID + 1; i <= atomEndID; i++ )
     {
      if ( atom[ i ].x < minX ) minX = atom[ i ].x;      
      if ( atom[ i ].x > maxX ) maxX = atom[ i ].x;      
      
      if ( atom[ i ].y < minY ) minY = atom[ i ].y;      
      if ( atom[ i ].y > maxY ) maxY = atom[ i ].y;      

      if ( atom[ i ].z < minZ ) minZ = atom[ i ].z;      
      if ( atom[ i ].z > maxZ ) maxZ = atom[ i ].z;      
     } 
   
   octree[ nodeID ].cx = ( minX + maxX ) / 2;
   octree[ nodeID ].cy = ( minY + maxY ) / 2;
   octree[ nodeID ].cz = ( minZ + maxZ ) / 2;

   double r2 = ( atom[ atomStartID ].x - octree[ nodeID ].cx ) * ( atom[ atomStartID ].x - octree[ nodeID ].cx )
             + ( atom[ atomStartID ].y - octree[ nodeID ].cy ) * ( atom[ atomStartID ].y - octree[ nodeID ].cy )
             + ( atom[ atomStartID ].z - octree[ nodeID ].cz ) * ( atom[ atomStartID ].z - octree[ nodeID ].cz );
   
   for ( int i = atomStartID + 1; i <= atomEndID; i++ )
     {
      double r2T = ( atom[ i ].x - octree[ nodeID ].cx ) * ( atom[ i ].x - octree[ nodeID ].cx )
                 + ( atom[ i ].y - octree[ nodeID ].cy ) * ( atom[ i ].y - octree[ nodeID ].cy )
                 + ( atom[ i ].z - octree[ nodeID ].cz ) * ( atom[ i ].z - octree[ nodeID ].cz );
     
      if ( r2T > r2 ) r2 = r2T;
     } 
   
   octree[ nodeID ].cr = sqrt( r2 );
         
   if ( ( atomEndID - atomStartID +  1 <= maxLeafAtoms ) || ( octree[ nodeID ].cr <= minRad ) )
     {
      octree[ nodeID ].leaf = true;
              
      for ( int i = atomStartID; i <= atomEndID; i++ )
         qSum[ qPtr + atom[ i ].groupR ] += atom[ i ].q;
     }
   else
     {
      octree[ nodeID ].leaf = false;     
      
      int atomCount[ 8 ] = { 0, 0, 0, 0, 0, 0, 0, 0 };
      
      for ( int i = atomStartID; i <= atomEndID; i++ )
        {
         int j = ( zeroIfLess( atom[ i ].z, octree[ nodeID ].cz ) << 2 )
               + ( zeroIfLess( atom[ i ].y, octree[ nodeID ].cy ) << 1 )
               + ( zeroIfLess( atom[ i ].x, octree[ nodeID ].cx ) );
         
         atomCount[ j ]++;
        }
        
      int atomStartIndex[ 8 ];

      atomStartIndex[ 0 ] = atomStartID;
      for ( int i = 1; i < 8; i++ )
        atomStartIndex[ i ] = atomStartIndex[ i - 1 ] + atomCount[ i - 1 ];

      for ( int i = 0; i < 8; i++ ) 
        if ( atomCount[ i ] > 0 )
          {
           int j = ( *freeID )++;
           octree[ nodeID ].cPtr[ i ] = j; 
         
           octree[ j ].atomStartID = atomStartIndex[ i ];
           octree[ j ].atomEndID = atomStartIndex[ i ] + atomCount[ i ] - 1;         
         
           constructOctree( octree, j, qSum, numGroupR, minRad, maxLeafAtoms, atom, freeID );
         
           for ( int k = 0; k < numGroupR; k++ )
             qSum[ qPtr + k ] += qSum[ numGroupR * j + k ];
          }
        else octree[ nodeID ].cPtr[ i ] = -1;           
     }  
}


void fastGpol::generateJobs( OCTREE_NODE *octree, int nodeI, int nodeJ, double eps, std::vector< JOB > &leafJobs, std::vector< JOB > &nonLeafJobs )
{   
   if ( octree[ nodeI ].leaf && octree[ nodeJ ].leaf )
     {            
       JOB job;
       
       job.nodeI = nodeI;
       job.nodeJ = nodeJ;       
       
       leafJobs.push_back( job );
     }
   else
     {
      double sumRad = octree[ nodeI ].cr + octree[ nodeJ ].cr;
      double rij2 = ( octree[ nodeI ].cx - octree[ nodeJ ].cx ) * ( octree[ nodeI ].cx - octree[ nodeJ ].cx )
                  + ( octree[ nodeI ].cy - octree[ nodeJ ].cy ) * ( octree[ nodeI ].cy - octree[ nodeJ ].cy )
                  + ( octree[ nodeI ].cz - octree[ nodeJ ].cz ) * ( octree[ nodeI ].cz - octree[ nodeJ ].cz );
      double rij = sqrt( rij2 ) - sumRad;
          
      if ( rij > ( 2 / eps ) * sumRad )
        {   
          JOB job;
       
          job.nodeI = nodeI;
          job.nodeJ = nodeJ;       

          leafJobs.push_back( job );        
        }
      else
        {
         if ( !octree[ nodeI ].leaf && !octree[ nodeJ ].leaf )  
           {
            for ( int i = 0; i < 8; i++ )
              if ( octree[ nodeI ].cPtr[ i ] >= 0 )
                 for ( int j = 0; j < 8; j++ )
                   if ( octree[ nodeJ ].cPtr[ j ] >= 0 ) 
                     {
                       JOB job;
                     
                       job.nodeI = octree[ nodeI ].cPtr[ i ];
                       job.nodeJ = octree[ nodeJ ].cPtr[ j ];       
      
                       nonLeafJobs.push_back( job );        
                     }                         
           }
         else
           {
            if ( octree[ nodeJ ].leaf )
              {
               int t = nodeI;
               nodeI = nodeJ;
               nodeJ = t;
              }

            for ( int j = 0; j < 8; j++ )
              if ( octree[ nodeJ ].cPtr[ j ] >= 0 )
                {
                  JOB job;
                 
                  job.nodeI = nodeI;
                  job.nodeJ = octree[ nodeJ ].cPtr[ j ];       
                 
                  nonLeafJobs.push_back( job );        
                }                          
           }           
        }              
     }   
}



void fastGpol::recursiveFastGpol( OCTREE_NODE *octree, int nodeI, int nodeJ, 
                                double *qSum, double *approxR, int numGroupR, ATOM *atom, 
                                double eps, double *Gpol )
{
   double sumRad = octree[ nodeI ].cr + octree[ nodeJ ].cr;
   double sumRad2 = sumRad * sumRad;    
   double rij2 = ( octree[ nodeI ].cx - octree[ nodeJ ].cx ) * ( octree[ nodeI ].cx - octree[ nodeJ ].cx )
               + ( octree[ nodeI ].cy - octree[ nodeJ ].cy ) * ( octree[ nodeI ].cy - octree[ nodeJ ].cy )
               + ( octree[ nodeI ].cz - octree[ nodeJ ].cz ) * ( octree[ nodeI ].cz - octree[ nodeJ ].cz );
   double rij = sqrt( rij2 ) - sumRad;
   rij2 = rij * rij;
   
   int qPtrI = numGroupR * nodeI, qPtrJ = numGroupR * nodeJ;   
   
   if ( octree[ nodeI ].leaf && octree[ nodeJ ].leaf )
     {            
#ifdef USE_SSE     
         V4SF xi, yi, zi, Ri, qi;
         V4SF xj, yj, zj, Rj, qj;         
         int k = 0;
#endif         
     
         for ( int i = octree[ nodeI ].atomStartID; i <= octree[ nodeI ].atomEndID; i++ )
            for ( int j = octree[ nodeJ ].atomStartID; j <= octree[ nodeJ ].atomEndID; j++ )
//              if ( i != j )
#ifdef USE_SSE
                 {
                   xi.f[ k ] = atom[ i ].x;
                   yi.f[ k ] = atom[ i ].y;
                   zi.f[ k ] = atom[ i ].z;                   
                   Ri.f[ k ] = atom[ i ].R;
                   qi.f[ k ] = atom[ i ].q;                   

                   xj.f[ k ] = atom[ j ].x;
                   yj.f[ k ] = atom[ j ].y;
                   zj.f[ k ] = atom[ j ].z;                   
                   Rj.f[ k ] = atom[ j ].R;
                   qj.f[ k ] = atom[ j ].q;                   
                                      
                   if ( ++k == 4 ) 
                     {
                       *Gpol += vectorGB( xi.v, yi.v, zi.v, Ri.v, qi.v, xj.v, yj.v, zj.v, Rj.v, qj.v );
                       k = 0;
                     }  
                 }
#else
                 {
                  double rij2 = ( atom[ i ].x - atom[ j ].x ) * ( atom[ i ].x - atom[ j ].x )
                              + ( atom[ i ].y - atom[ j ].y ) * ( atom[ i ].y - atom[ j ].y )
                              + ( atom[ i ].z - atom[ j ].z ) * ( atom[ i ].z - atom[ j ].z );                              
                              
                  if ( rij2 < 0.5 * 0.5 * ( atom[ i ].r + atom[ j ].r ) * ( atom[ i ].r + atom[ j ].r ) )            
                       rij2 = 0.5 * 0.5 * ( atom[ i ].r + atom[ j ].r ) * ( atom[ i ].r + atom[ j ].r );
                  
                  double RiRj = atom[ i ].R * atom[ j ].R;
                  double qiqj = atom[ i ].q * atom[ j ].q;            
                    
                  if ( useApproxMath )                  
                      *Gpol += ( qiqj * invSqrt( ( float ) ( rij2 + RiRj * fastExp( ( float ) ( - rij2 / ( 4 * RiRj ) ) ) ) ) );                                                  
                  else                                           
                      *Gpol += ( qiqj / sqrt( rij2 + RiRj * exp( - rij2 / ( 4 * RiRj ) ) ) );            
                 }
#endif                 

#ifdef USE_SSE
         while ( k-- )
           {
              double rij2 = ( xi.f[ k ] - xj.f[ k ] ) * ( xi.f[ k ] - xj.f[ k ] )
                          + ( yi.f[ k ] - yj.f[ k ] ) * ( yi.f[ k ] - yj.f[ k ] )
                          + ( zi.f[ k ] - zj.f[ k ] ) * ( zi.f[ k ] - zj.f[ k ] );                              
              
              double RiRj = Ri.f[ k ] * Rj.f[ k ];
              double qiqj = qi.f[ k ] * qj.f[ k ];            

              if ( useApproxMath )                  
                 *Gpol += ( qiqj * invSqrt( ( float ) ( rij2 + RiRj * fastExp( ( float ) ( - rij2 / ( 4 * RiRj ) ) ) ) ) );                                                  
              else                                          
                 *Gpol += ( qiqj / sqrt( rij2 + RiRj * exp( - rij2 / ( 4 * RiRj ) ) ) );            
           }
#endif
     }
   else
     {
      if ( rij > ( 2 / eps ) * sumRad ) //( rij2 > eps2 * sumRad2 ) 
        {   
#ifdef USE_SSE   
            V4SF Ri, qi;
            V4SF Rj, qj;         
            int k = 0;
#endif         
        
            for ( int i = qPtrI; i < qPtrI + numGroupR; i++ )
               if ( qSum[ i ] != 0 )
                 {
                  for ( int j = qPtrJ; j < qPtrJ + numGroupR; j++ )
                    if ( qSum[ j ] != 0 )
#ifdef USE_SSE
                      {
                        Ri.f[ k ] = approxR[ i - qPtrI ];
                        qi.f[ k ] = qSum[ i ];                   
    
                        Rj.f[ k ] = approxR[ j - qPtrJ ];
                        qj.f[ k ] = qSum[ j ];                   
                                          
                        if ( ++k == 4 ) 
                          {
                            *Gpol += vectorGB( rij2, Ri.v, qi.v, Rj.v, qj.v );
                            k = 0;
                          }  
                      }
#else                    
                      {
                        double RiRj = approxR[ i - qPtrI ] * approxR[ j - qPtrJ ];
                        double qiqj = qSum[ i ] * qSum[ j ];            

                        if ( useApproxMath )                  
                           *Gpol += ( qiqj * invSqrt( ( float ) ( rij2 + RiRj * fastExp( ( float ) ( - rij2 / ( 4 * RiRj ) ) ) ) ) );
                        else
                           *Gpol += ( qiqj / sqrt( rij2 + RiRj * exp( - rij2 / ( 4 * RiRj ) ) ) );
                      }
#endif                                                       
                 }
                 
#ifdef USE_SSE
         while ( k-- )
           {
              double RiRj = Ri.f[ k ] * Rj.f[ k ];
              double qiqj = qi.f[ k ] * qj.f[ k ];            

              if ( useApproxMath )                  
                 *Gpol += ( qiqj * invSqrt( ( float ) ( rij2 + RiRj * fastExp( ( float ) ( - rij2 / ( 4 * RiRj ) ) ) ) ) );                                                  
              else
                 *Gpol += ( qiqj / sqrt( rij2 + RiRj * exp( - rij2 / ( 4 * RiRj ) ) ) );            
           }
#endif                                       
        }
      else
        {
         double GpolT;
         
         if ( !octree[ nodeI ].leaf && !octree[ nodeJ ].leaf )  
           {
            for ( int i = 0; i < 8; i++ )
              if ( octree[ nodeI ].cPtr[ i ] >= 0 )
                 for ( int j = 0; j < 8; j++ )
                   if ( octree[ nodeJ ].cPtr[ j ] >= 0 ) 
                     {
                      GpolT = 0;
                      recursiveFastGpol( octree, octree[ nodeI ].cPtr[ i ], octree[ nodeJ ].cPtr[ j ], qSum, approxR, numGroupR, atom, eps, &GpolT );
                      *Gpol += GpolT;            
                     }                         
           }
         else
           {
            if ( octree[ nodeJ ].leaf )
              {
               int t = nodeI;
               nodeI = nodeJ;
               nodeJ = t;
              }

            for ( int j = 0; j < 8; j++ )
              if ( octree[ nodeJ ].cPtr[ j ] >= 0 )
                {
                 GpolT = 0;                
                 recursiveFastGpol( octree, nodeI, octree[ nodeJ ].cPtr[ j ], qSum, approxR, numGroupR, atom, eps, &GpolT );
                 *Gpol += GpolT;            
                }                          
           }           
        }              
     }   
}



void fastGpol::computeThreadedFastGpol( int threadID, OCTREE_NODE *octree,
                                        double *qSum, double *approxR, int numGroupR, ATOM *atom, 
                                        double eps, double *Gpol )	     	    
{
   int nodeI, nodeJ;
   
   while ( nextJob( nodeI, nodeJ ) )
       recursiveFastGpol( octree, nodeI, nodeJ, qSum, approxR, numGroupR, atom, eps, Gpol );
}



void fastGpol::computeFastGpol( double *Gpol )
{
   if ( printStatus ) printf( (char *)"\napproximating G_pol... " );

   double startT = getTime( );

   double eps2 = ( 1 + 1 / epsilon ) * ( 1 + 1 / epsilon );
   
   if ( eps2 < 4.0 ) eps2 = 4.0;

#ifdef USE_SSE
   ZERO = _mm_setzero_ps( );
   FOUR = _mm_set1_ps( 4.0f );   
#endif   
   
   std::vector< JOB > leafJobs;
   std::vector< JOB > nonLeafJobs;   
   std::vector< JOB > sourceJobs;      
   
   JOB job;
   
   job.nodeI = job.nodeJ = 0;
   nonLeafJobs.push_back( job );
   
   while ( leafJobs.size( ) + nonLeafJobs.size( ) < numThreads * numThreads )   
     {              
       sourceJobs.clear( );

       int k = nonLeafJobs.size( );
       
       while ( k-- ) 
         {
           sourceJobs.push_back( nonLeafJobs.back( ) );
           nonLeafJobs.pop_back( );
         }  
       
       nonLeafJobs.clear( );       
       
       k = sourceJobs.size( );
       
       while ( k-- )
         {
           job = sourceJobs.back( );           
           generateJobs( octree, job.nodeI, job.nodeJ, epsilon, leafJobs, nonLeafJobs );
           sourceJobs.pop_back( );                      
         }
     }

   PARAMS prT[ numThreads ];
   pthread_t p[ numThreads ];           

   initJobServer( leafJobs, nonLeafJobs );

   for ( int i = 0; i < numThreads; i++ )
     {              
       prT[ i ].threadID = i + 1;
       prT[ i ].octree = octree;
       prT[ i ].qSum = qSum;
       prT[ i ].approxR = approxR;
       prT[ i ].numGroupR = maxGroupR + 1;
       prT[ i ].atom = atoms;
       prT[ i ].eps = epsilon;
       prT[ i ].Gpol = 0;
       prT[ i ].thisPtr = this;
               
       pthread_create( &p[ i ], NULL, startThreadedFastGpol, ( void * ) &prT[ i ] );
     }

   for ( int i = 0; i < numThreads; i++ )
     pthread_join( p[ i ], NULL );
     
   *Gpol = 0;
   for ( int i = 0; i < numThreads; i++ )  
     *Gpol += prT[ i ].Gpol;

   double tau = 1 / e_in - 1 / e_out;
   *Gpol = - 332 * 0.5 * tau * ( *Gpol );  

   double endT = getTime( );

   if ( printStatus ) printf( (char *)"done ( %lf sec )\n", endT - startT );
}

