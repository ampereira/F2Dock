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

#include "pseudoGsol.h"

void pseudoGsol::freeMemory( void )
{
   freeMem( staticAtoms );
   freeMem( movingAtoms );   

   freeMem( staticQPoints );
   freeMem( movingQPoints );  

   freeMem( staticQPointsPG );
   freeMem( movingQPointsPG );  
   
   freeMem( staticAtomsOctree );
   freeMem( movingAtomsOctree );

   freeMem( staticQPointsOctree );
   freeMem( movingQPointsOctree );
   
   freeMem( staticQPointsOctreeFlags );
   freeMem( movingQPointsOctreeFlags );   
   
   delete staticPG;
   delete movingPG;   
   
   delete hydroIndex;
}


void pseudoGsol::setDefaults( void )
{
  numThreads = 4;  
    
  minRadius = 4.0;
  maxLeafSize = 500;

  numStaticAtoms = numMovingAtoms = 0;  
  staticAtoms = movingAtoms = NULL;

  numStaticQPoints = numMovingQPoints = 0;  
  staticQPoints = movingQPoints = NULL;

  staticQPointsPG = movingQPointsPG = NULL;

  numStaticAtomsOctreeNodes = numMovingAtomsOctreeNodes = 0;
  staticAtomsOctree = movingAtomsOctree = NULL;
  staticAtomsOctreeRoot = movingAtomsOctreeRoot = -1;

  numStaticQPointsOctreeNodes = numMovingQPointsOctreeNodes = 0;
  staticQPointsOctree = movingQPointsOctree = NULL;
  staticQPointsOctreeRoot = movingQPointsOctreeRoot = -1;
  
  staticQPointsOctreeFlags = movingQPointsOctreeFlags = NULL;
  
  preprocessElementInformationTable( &hydroIndex, &nHydroIndex );  
}


void pseudoGsol::processQPoints( void )
{
  staticQPointsPG = ( Point * ) malloc( numStaticQPoints * sizeof( Point ) );
  movingQPointsPG = ( Point * ) malloc( numMovingQPoints * sizeof( Point ) );
               
  if ( ( staticQPointsPG == NULL ) || ( movingQPointsPG == NULL ) )
     {
       printError( (char *)"Failed to allocate memory!" );  
       exit( 1 );                 
     }
                 
  double xlatePG = computeXlateForPG( numStaticQPoints, staticQPoints, numMovingQPoints, movingQPoints );

  staticPG = new PG( 10.0, xlatePG, 5.0 );  

  for ( int i = 0; i < numStaticQPoints; i++ )
     {
       staticQPointsPG[ i ].x = staticQPoints[ i ].x;  
       staticQPointsPG[ i ].y = staticQPoints[ i ].y;  
       staticQPointsPG[ i ].z = staticQPoints[ i ].z;  
                   
       staticPG->addPoint( &staticQPointsPG[ i ] );
     }    

  xlatePG = computeXlateForPG( numMovingQPoints, movingQPoints, numStaticQPoints, staticQPoints );

  movingPG = new PG( 10.0, xlatePG, 5.0 );  
     
  for ( int i = 0; i < numMovingQPoints; i++ )
     {
       movingQPointsPG[ i ].x = movingQPoints[ i ].x;  
       movingQPointsPG[ i ].y = movingQPoints[ i ].y;  
       movingQPointsPG[ i ].z = movingQPoints[ i ].z;  
       
       movingPG->addPoint( &movingQPointsPG[ i ] );
     }       
     
  if ( !buildAtomsOctree( numStaticAtoms, staticAtoms, &numStaticAtomsOctreeNodes, &staticAtomsOctree, &staticAtomsOctreeRoot ) 
    || !buildAtomsOctree( numMovingAtoms, movingAtoms, &numMovingAtomsOctreeNodes, &movingAtomsOctree, &movingAtomsOctreeRoot )   
    || !buildQPointsOctree( numStaticQPoints, staticQPoints, &numStaticQPointsOctreeNodes, &staticQPointsOctree, &staticQPointsOctreeRoot ) 
    || !buildQPointsOctree( numMovingQPoints, movingQPoints, &numMovingQPointsOctreeNodes, &movingQPointsOctree, &movingQPointsOctreeRoot ) )
      {
        exit( 1 );
      }  
       
  assignHydrophobicityToQPoints( staticAtomsOctree, staticAtomsOctreeRoot, staticAtoms,  
                                 staticQPointsOctree, staticQPointsOctreeRoot, staticQPoints, 4, 0.1 );       
  assignHydrophobicityToQPoints( movingAtomsOctree, movingAtomsOctreeRoot, movingAtoms,
                                 movingQPointsOctree, movingQPointsOctreeRoot, movingQPoints, 4, 0.1 );         
                                 
  staticQPointsOctreeFlags = ( bool * ) malloc( 2 * numThreads * numStaticQPointsOctreeNodes * sizeof( bool ) );
  movingQPointsOctreeFlags = ( bool * ) malloc( 2 * numThreads * numMovingQPointsOctreeNodes * sizeof( bool ) );
               
  if ( ( staticQPointsOctreeFlags == NULL ) || ( movingQPointsOctreeFlags == NULL ) )
     {
       printError( (char *)"Failed to allocate memory!" );  
       exit( 1 );                 
     }  
               		
}



pseudoGsol::pseudoGsol( char *paramFile, int nStAtoms, double *stAtoms, int nMvAtoms, double *mvAtoms, int nThreads )
{
  setDefaults( );
  
  if ( !getParamsFromFile( &params, paramFile, false ) ) exit( 1 );
  
  if ( nThreads <= 0 )
    {
      printError( (char *)"Invalid number of threads ( %d )!", numThreads );  
      exit( 1 );                     
    }
    
  numThreads = nThreads;  
    
  if ( !copyAtomsFromArray( nStAtoms, stAtoms, &numStaticAtoms, &staticAtoms ) 
    || !copyAtomsFromArray( nMvAtoms, mvAtoms, &numMovingAtoms, &movingAtoms ) 
    || !readQPoints( params.staticMoleculeQUAD, &numStaticQPoints, &staticQPoints )
    || !readQPoints( params.movingMoleculeQUAD, &numMovingQPoints, &movingQPoints ) )
      {
        exit( 1 );
      }  

  processQPoints( );
}


pseudoGsol::pseudoGsol( char *paramFile, int nThreads )
{
  setDefaults( );
  
  if ( !getParamsFromFile( &params, paramFile, true ) ) exit( 1 );  
  
  if ( nThreads <= 0 )
    {
      printError( (char *)"Invalid number of threads ( %d )!", numThreads );  
      exit( 1 );                     
    }
    
  numThreads = nThreads;    

  if ( !readAtoms( params.staticMoleculePQR, &numStaticAtoms, &staticAtoms ) 
    || !readAtoms( params.movingMoleculePQR, &numMovingAtoms, &movingAtoms ) 
    || !readQPoints( params.staticMoleculeQUAD, &numStaticQPoints, &staticQPoints )
    || !readQPoints( params.movingMoleculeQUAD, &numMovingQPoints, &movingQPoints ) )
      {
        exit( 1 );
      }  

  processQPoints( );
}


pseudoGsol::~pseudoGsol( )
{
  freeMemory( );
}


bool pseudoGsol::readQPoints( char *qPtsFile, int *nQPoints, QPOINT **qPoints )
{
   FILE *fp;
   
   fp = fopen( qPtsFile, "rt" );
   
   if ( fp == NULL )
     {
      printError( (char *)"Failed to open quadrature points file (%s)!", qPtsFile );
      return false;
     }
   
   *nQPoints = 0;
   QPOINT qPt;  
   double nx, ny, nz;
          
   while ( fscanf( fp, (char *)"%lf %lf %lf %lf %lf %lf %lf", &qPt.x, &qPt.y, &qPt.z, &nx, &ny, &nz, &qPt.w ) == 7 ) ( *nQPoints )++; 

   fclose( fp );  
   
   ( *qPoints ) = ( QPOINT * ) malloc( ( *nQPoints ) * sizeof( QPOINT ) );
   
   if ( *qPoints == NULL )
     {
      printError( (char *)"Failed to allocate memory for quadrature points!" );
      return false;
     }
     
   fp = fopen( qPtsFile, "rt" );
   
   if ( fp == NULL )
     {
      printError( (char *)"Failed to open quadrature points file (%s)!", qPtsFile );
      return false;
     }

   for ( int i = 0; i < *nQPoints; i++ )   
     {          
      if ( fscanf( fp, (char *)"%lf %lf %lf %lf %lf %lf %lf", &qPt.x, &qPt.y, &qPt.z, &nx, &ny, &nz, &qPt.w ) != 7 )
        {
         printError( (char *)"Failed to read the quadrature points file (%s)!", qPtsFile );
         return false;
        }

      qPt.h = 0;          
      
      ( *qPoints )[ i ] = qPt;         
     }    
     
   fclose( fp );
         
   return true;
}


bool pseudoGsol::copyAtomsFromArray( int numAtoms, double *atms, int *nAtoms, ATOM **atoms )
{
   if ( numAtoms <= 0 )
     {
      printError( (char *)"No atoms to copy!" );
      return false;
     }
   
   *nAtoms = numAtoms;
   ( *atoms ) = ( ATOM * ) malloc( ( *nAtoms ) * sizeof( ATOM ) );
   
   if ( *atoms == NULL )
     {
      printError( (char *)"Failed to allocate memory for atoms!" );
      return false;
     }

   for ( int i = 0; i < *nAtoms; i++ )   
     {          
      ( *atoms )[ i ].x = atms[ 5 * i + 0 ];
      ( *atoms )[ i ].y = atms[ 5 * i + 1 ];
      ( *atoms )[ i ].z = atms[ 5 * i + 2 ];

      ( *atoms )[ i ].r = atms[ 5 * i + 3 ];      
      ( *atoms )[ i ].h = atms[ 5 * i + 4 ];      

      ( *atoms )[ i ].id = i + 1;
     }    
              
   return true;
}


void pseudoGsol::preprocessElementInformationTable( int **index, int *nIndex )
{
  int c = 1;
  int n = sizeof( elementTable ) / sizeof( elementTable[ 0 ] );

  for ( int i = 1; i < n; i++ )
    if ( strcmp( elementTable[ i ].residueName, elementTable[ i - 1 ].residueName ) ) c++;

  ( *index ) = new int [ c ];

  *nIndex = c;

  c = 0;
  ( *index )[ 0 ] = 0;
  for ( int i = 1; i < n; i++ )
    if ( strcmp( elementTable[ i ].residueName, elementTable[ i - 1 ].residueName ) ) 
        ( *index )[ ++c ] = i;  
}



int pseudoGsol::strcmp_nospace( char *s1, char *s2 )
{
   int i = 0, j = 0;
   
   while ( s1[ i ] && s2[ j ] )
     {
       if ( !isspace( s1[ i ] ) && !isspace( s2[ j ] ) )
         {
           if ( s1[ i ] == s2[ j ] )
             {
               i++; j++;
             }
           else return ( ( int ) ( s1[ i ] - s2[ j ] ) );  
         }
         
       if ( isspace( s1[ i ] ) ) i++;
       if ( isspace( s2[ j ] ) ) j++;       
     }

   while ( isspace( s1[ i ] ) ) i++;
   while ( isspace( s2[ j ] ) ) j++;       
     
   return ( ( int ) ( s1[ i ] - s2[ j ] ) );  
}


float pseudoGsol::getHydrophobicity( char *atomName, char *residueName, int *index, int nIndex, bool useInterfacePropensity, bool perResidueHydrophobicity )
{
  for ( int i = 0; i < nIndex; i++ )
    { 
      int j = index[ i ];
      
      if ( !strcmp_nospace( residueName, elementTable[ j ].residueName ) )
        {
          int k = ( i == nIndex - 1 ) ? MAX_BIOCHEM_ELEMENTS : index[ i + 1 ];
          
          while ( j < k )
            {
              if ( !strcmp_nospace( atomName, elementTable[ j ].atomName ) )
                {
                  if ( useInterfacePropensity ) return ( ( float ) elementTable[ j ].interfacePropensity );
                  
                  if ( perResidueHydrophobicity ) return ( ( float ) elementTable[ j ].perResidueHydrophobicity );
                  else return ( ( float ) elementTable[ j ].hydrophobicity );    
                }  
                
              j++;  
            }
        }
    }
    
  return 0;  
}



bool pseudoGsol::readAtoms( char *atomsFile, int *nAtoms, ATOM **atoms )
{
   FILE *fp;
   
   fp = fopen( atomsFile, "rt" );
   
   if ( fp == NULL )
     {
      printError( (char *)"Failed to open PQR file (%s)!", atomsFile );
      return false;
     }
   
   int numAtoms = 0;
   
   char line[ 101 ];

   while ( fgets( line, 100, fp ) != NULL )
     {
       int i = 0;
       while ( isspace( line[ i ] ) ) i++;
       if ( strstr( line + i, (char *)"ATOM" ) == line + i ) numAtoms++;
     }
   
   fclose( fp );  
   
   *nAtoms = numAtoms;
   ( *atoms ) = ( ATOM * ) malloc( ( *nAtoms ) * sizeof( ATOM ) );
   
   if ( *atoms == NULL )
     {
      printError( (char *)"Failed to allocate memory for atoms!" );
      return false;
     }
     
   fp = fopen( atomsFile, "rt" );
   
   if ( fp == NULL )
     {
      printError( (char *)"Failed to open PQR file (%s)!", atomsFile );
      return false;
     }

   ATOM atm;  

   int k = 0;
   while ( fgets( line, 100, fp ) != NULL )
     {
      int i = 0, j;
      double v;
      char tmp[ 100 ], anm[ 100 ];    
    
      i = getAlphaString( line, i, tmp );   // get 'ATOM'/'HETATM', and ignore
      
      if ( strcmp( tmp, "ATOM" ) ) continue;

      i = getInt( line, i, &j );            // get atom number, and ignore

      i = getString( line, i, anm );        // get atom name

      i = getAlphaString( line, i, tmp );   // get residue name
      
      ( *atoms )[ k ].h = getHydrophobicity( anm, tmp, hydroIndex, nHydroIndex, params.useInterfacePropensity, params.perResidueHydrophobicity );

      i = getInt( line, i, &j );            // get residue number, and ignore

      i = getDouble( line, i, &v );         // get X coordinate
      ( *atoms )[ k ].x = v;    
    
      i = getDouble( line, i, &v );         // get Y coordinate
      ( *atoms )[ k ].y = v;    

      i = getDouble( line, i, &v );         // get Z coordinate
      ( *atoms )[ k ].z = v;    

      i = getDouble( line, i, &v );         // get charge

      i = getDouble( line, i, &v );         // get radius
      ( *atoms )[ k ].r = v;  
      
//      switch ( anm[ 0 ] )
//        { 
//          case 'N' : ( *atoms )[ k ].r = vdwRad[ N ];
//                     break;
//
//          case 'C' : ( *atoms )[ k ].r = vdwRad[ C ];
//                     break;
//
//          case 'O' : ( *atoms )[ k ].r = vdwRad[ O ];
//                     break;
//
//          case 'H' : ( *atoms )[ k ].r = vdwRad[ H ];
//                     break;
//
//          case 'S' : ( *atoms )[ k ].r = vdwRad[ S ];
//                     break;
//
//          case 'P' : ( *atoms )[ k ].r = vdwRad[ P ];
//                     break;
//                     
//          default  : ( *atoms )[ k ].r = v;    
//                     break;           
//        }        
      
      ( *atoms )[ k ].id = k + 1;

      k++;
     }
     
   fclose( fp );
         
   return true;
}



void pseudoGsol::initFreeNodeServer( int numNodes )
{
   pthread_mutex_init( &nodesLock, NULL );
   curNode = 0;
   maxNode = numNodes;
}


int pseudoGsol::nextFreeNode( void )
{
   int nextNode = -1;
   
   pthread_mutex_lock( &nodesLock );
   if ( curNode < maxNode ) nextNode = curNode++;
   pthread_mutex_unlock( &nodesLock );
   
   return nextNode;
}


void pseudoGsol::countAtomsOctreeNodesAndSortAtoms( ATOM *atoms, int atomsStartID, int atomsEndID, 
                                                    ATOM *atomsT, int *numNodes )
{
   double minX = atoms[ atomsStartID ].x, minY = atoms[ atomsStartID ].y, minZ = atoms[ atomsStartID ].z;
   double maxX = atoms[ atomsStartID ].x, maxY = atoms[ atomsStartID ].y, maxZ = atoms[ atomsStartID ].z;
   
   for ( int i = atomsStartID + 1; i <= atomsEndID; i++ )
     {
      if ( atoms[ i ].x < minX ) minX = atoms[ i ].x;      
      if ( atoms[ i ].x > maxX ) maxX = atoms[ i ].x;      
      
      if ( atoms[ i ].y < minY ) minY = atoms[ i ].y;      
      if ( atoms[ i ].y > maxY ) maxY = atoms[ i ].y;      

      if ( atoms[ i ].z < minZ ) minZ = atoms[ i ].z;      
      if ( atoms[ i ].z > maxZ ) maxZ = atoms[ i ].z;      
     } 
   
   double cx = ( minX + maxX ) / 2,
          cy = ( minY + maxY ) / 2,
          cz = ( minZ + maxZ ) / 2;

   double r2 = ( atoms[ atomsStartID ].x - cx ) * ( atoms[ atomsStartID ].x - cx )
             + ( atoms[ atomsStartID ].y - cy ) * ( atoms[ atomsStartID ].y - cy )
             + ( atoms[ atomsStartID ].z - cz ) * ( atoms[ atomsStartID ].z - cz );
   
   for ( int i = atomsStartID + 1; i <= atomsEndID; i++ )
     {
      double r2T = ( atoms[ i ].x - cx ) * ( atoms[ i ].x - cx )
                 + ( atoms[ i ].y - cy ) * ( atoms[ i ].y - cy )
                 + ( atoms[ i ].z - cz ) * ( atoms[ i ].z - cz );
     
      if ( r2T > r2 ) r2 = r2T;
     } 
   
   *numNodes = 1;
      
   if ( ( atomsEndID - atomsStartID +  1 > maxLeafSize ) && ( r2 > minRadius * minRadius ) )
     {
      int atomsCount[ 8 ] = { 0, 0, 0, 0, 0, 0, 0, 0 };
      
      for ( int i = atomsStartID; i <= atomsEndID; i++ )
        {
         atomsT[ i ] = atoms[ i ];
         
         int j = ( zeroIfLess( atoms[ i ].z, cz ) << 2 )
               + ( zeroIfLess( atoms[ i ].y, cy ) << 1 )
               + ( zeroIfLess( atoms[ i ].x, cx ) );
         
         atomsCount[ j ]++;
        }
        
      int atomsStartIndex[ 8 ];
      int atomsCurIndex[ 8 ];              
      
      atomsCurIndex[ 0 ] = atomsStartIndex[ 0 ] = atomsStartID;
      for ( int i = 1; i < 8; i++ )
        atomsCurIndex[ i ] = atomsStartIndex[ i ] = atomsStartIndex[ i - 1 ] + atomsCount[ i - 1 ];

      for ( int i = atomsStartID; i <= atomsEndID; i++ )
        {        
         int j = ( zeroIfLess( atomsT[ i ].z, cz ) << 2 )
               + ( zeroIfLess( atomsT[ i ].y, cy ) << 1 )
               + ( zeroIfLess( atomsT[ i ].x, cx ) );
           
         atoms[ atomsCurIndex[ j ] ] = atomsT[ i ];
         atomsCurIndex[ j ]++;  
        }        
        
      for ( int i = 0; i < 8; i++ ) 
        if ( atomsCount[ i ] > 0 )
          {
           int numNodesT = 0;
           
           countAtomsOctreeNodesAndSortAtoms( atoms, atomsStartIndex[ i ], atomsStartIndex[ i ] + atomsCount[ i ] - 1, atomsT, &numNodesT );
         
           *numNodes += numNodesT;
          }
     }  
}



int pseudoGsol::constructAtomsOctree( int atomsStartID, int atomsEndID, ATOM *atoms, ATOMS_OCTREE_NODE *atomsOctree )
{
   int nodeID = nextFreeNode( );
   
   atomsOctree[ nodeID ].atomsStartID = atomsStartID;
   atomsOctree[ nodeID ].atomsEndID = atomsEndID;
   
   double minX = atoms[ atomsStartID ].x, minY = atoms[ atomsStartID ].y, minZ = atoms[ atomsStartID ].z;
   double maxX = atoms[ atomsStartID ].x, maxY = atoms[ atomsStartID ].y, maxZ = atoms[ atomsStartID ].z;
   
   for ( int i = atomsStartID + 1; i <= atomsEndID; i++ )
     {
      if ( atoms[ i ].x < minX ) minX = atoms[ i ].x;      
      if ( atoms[ i ].x > maxX ) maxX = atoms[ i ].x;      
      
      if ( atoms[ i ].y < minY ) minY = atoms[ i ].y;      
      if ( atoms[ i ].y > maxY ) maxY = atoms[ i ].y;      

      if ( atoms[ i ].z < minZ ) minZ = atoms[ i ].z;      
      if ( atoms[ i ].z > maxZ ) maxZ = atoms[ i ].z;      
     } 
   
   double cx = atomsOctree[ nodeID ].cx = ( minX + maxX ) / 2;
   double cy = atomsOctree[ nodeID ].cy = ( minY + maxY ) / 2;
   double cz = atomsOctree[ nodeID ].cz = ( minZ + maxZ ) / 2;

   double r2 = ( atoms[ atomsStartID ].x - cx ) * ( atoms[ atomsStartID ].x - cx )
             + ( atoms[ atomsStartID ].y - cy ) * ( atoms[ atomsStartID ].y - cy )
             + ( atoms[ atomsStartID ].z - cz ) * ( atoms[ atomsStartID ].z - cz );
   
   for ( int i = atomsStartID + 1; i <= atomsEndID; i++ )
     {
      double r2T = ( atoms[ i ].x - cx ) * ( atoms[ i ].x - cx )
                 + ( atoms[ i ].y - cy ) * ( atoms[ i ].y - cy )
                 + ( atoms[ i ].z - cz ) * ( atoms[ i ].z - cz );
     
      if ( r2T > r2 ) r2 = r2T;
     } 
   
   double cr = atomsOctree[ nodeID ].cr = sqrt( r2 );
         
   if ( ( atomsEndID - atomsStartID +  1 <= maxLeafSize ) || ( cr <= minRadius ) )
      atomsOctree[ nodeID ].leaf = true;
   else
     {
      atomsOctree[ nodeID ].leaf = false;     
      
      int atomsCount[ 8 ] = { 0, 0, 0, 0, 0, 0, 0, 0 };
      
      for ( int i = atomsStartID; i <= atomsEndID; i++ )
        {
         int j = ( zeroIfLess( atoms[ i ].z, cz ) << 2 )
               + ( zeroIfLess( atoms[ i ].y, cy ) << 1 )
               + ( zeroIfLess( atoms[ i ].x, cx ) );
         
         atomsCount[ j ]++;
        }
        
      int atomsStartIndex[ 8 ];

      atomsStartIndex[ 0 ] = atomsStartID;
      for ( int i = 1; i < 8; i++ )
        atomsStartIndex[ i ] = atomsStartIndex[ i - 1 ] + atomsCount[ i - 1 ];

      for ( int i = 0; i < 8; i++ ) 
        if ( atomsCount[ i ] > 0 )
          {
           int j = constructAtomsOctree( atomsStartIndex[ i ], atomsStartIndex[ i ] + atomsCount[ i ] - 1, atoms, atomsOctree );
           atomsOctree[ nodeID ].cPtr[ i ] = j; 
          }
        else atomsOctree[ nodeID ].cPtr[ i ] = -1;           
     }  
     
   return nodeID;  
}


bool pseudoGsol::buildAtomsOctree( int nAtoms, ATOM *atoms, int *numAtomsOctreeNodes, ATOMS_OCTREE_NODE **atomsOctree, int *atomsOctreeRoot )
{  
   ATOM *atomsT;         
   atomsT = ( ATOM * ) malloc( nAtoms * sizeof( ATOM ) );
   
   if ( atomsT == NULL )
     {
      printError( (char *)"Failed to allocate temporary memory for atoms!" );
      return false;
     }

   countAtomsOctreeNodesAndSortAtoms( atoms, 0, nAtoms - 1, atomsT, numAtomsOctreeNodes );

   freeMem( atomsT );
   
   ( *atomsOctree ) = ( ATOMS_OCTREE_NODE * ) malloc( ( *numAtomsOctreeNodes ) * sizeof( ATOMS_OCTREE_NODE ) );

   if ( *atomsOctree == NULL )
     {
      printError( (char *)"Unable to build atoms octree - memory allocation failed!" );
      return false;
     }
 
   initFreeNodeServer( nAtoms );
   
   *atomsOctreeRoot = constructAtomsOctree( 0, nAtoms - 1, atoms, *atomsOctree );
   
   return true;
}



void pseudoGsol::countQPointsOctreeNodesAndSortQPoints( QPOINT *qPts, int qPtsStartID, int qPtsEndID, 
                                                        QPOINT *qPtsT, int *numNodes )
{
   double minX = qPts[ qPtsStartID ].x, minY = qPts[ qPtsStartID ].y, minZ = qPts[ qPtsStartID ].z;
   double maxX = qPts[ qPtsStartID ].x, maxY = qPts[ qPtsStartID ].y, maxZ = qPts[ qPtsStartID ].z;
   
   for ( int i = qPtsStartID + 1; i <= qPtsEndID; i++ )
     {
      if ( qPts[ i ].x < minX ) minX = qPts[ i ].x;      
      if ( qPts[ i ].x > maxX ) maxX = qPts[ i ].x;      
      
      if ( qPts[ i ].y < minY ) minY = qPts[ i ].y;      
      if ( qPts[ i ].y > maxY ) maxY = qPts[ i ].y;      

      if ( qPts[ i ].z < minZ ) minZ = qPts[ i ].z;      
      if ( qPts[ i ].z > maxZ ) maxZ = qPts[ i ].z;      
     } 
   
   double cx = ( minX + maxX ) / 2,
          cy = ( minY + maxY ) / 2,
          cz = ( minZ + maxZ ) / 2;

   double r2 = ( qPts[ qPtsStartID ].x - cx ) * ( qPts[ qPtsStartID ].x - cx )
             + ( qPts[ qPtsStartID ].y - cy ) * ( qPts[ qPtsStartID ].y - cy )
             + ( qPts[ qPtsStartID ].z - cz ) * ( qPts[ qPtsStartID ].z - cz );
   
   for ( int i = qPtsStartID + 1; i <= qPtsEndID; i++ )
     {
      double r2T = ( qPts[ i ].x - cx ) * ( qPts[ i ].x - cx )
                 + ( qPts[ i ].y - cy ) * ( qPts[ i ].y - cy )
                 + ( qPts[ i ].z - cz ) * ( qPts[ i ].z - cz );
     
      if ( r2T > r2 ) r2 = r2T;
     } 
   
   *numNodes = 1;
      
   if ( ( qPtsEndID - qPtsStartID + 1 > maxLeafSize ) && ( r2 > minRadius * minRadius ) )
     {
      int qPtsCount[ 8 ] = { 0, 0, 0, 0, 0, 0, 0, 0 };
      
      for ( int i = qPtsStartID; i <= qPtsEndID; i++ )
        {
         qPtsT[ i ] = qPts[ i ];
         
         int j = ( zeroIfLess( qPts[ i ].z, cz ) << 2 )
               + ( zeroIfLess( qPts[ i ].y, cy ) << 1 )
               + ( zeroIfLess( qPts[ i ].x, cx ) );
         
         qPtsCount[ j ]++;
        }
      
      int qPtsStartIndex[ 8 ];
      int qPtsCurIndex[ 8 ];              
      
      qPtsCurIndex[ 0 ] = qPtsStartIndex[ 0 ] = qPtsStartID;
      for ( int i = 1; i < 8; i++ )
        qPtsCurIndex[ i ] = qPtsStartIndex[ i ] = qPtsStartIndex[ i - 1 ] + qPtsCount[ i - 1 ];

      for ( int i = qPtsStartID; i <= qPtsEndID; i++ )
        {        
         int j = ( zeroIfLess( qPtsT[ i ].z, cz ) << 2 )
               + ( zeroIfLess( qPtsT[ i ].y, cy ) << 1 )
               + ( zeroIfLess( qPtsT[ i ].x, cx ) );
           
         qPts[ qPtsCurIndex[ j ] ] = qPtsT[ i ];
         qPtsCurIndex[ j ]++;  
        }        
        
      for ( int i = 0; i < 8; i++ ) 
        if ( qPtsCount[ i ] > 0 )
          {
           int numNodesT = 0;
           
           countQPointsOctreeNodesAndSortQPoints( qPts, qPtsStartIndex[ i ], qPtsStartIndex[ i ] + qPtsCount[ i ] - 1, qPtsT, &numNodesT );
         
           *numNodes += numNodesT;
          }
     }  
}



int pseudoGsol::constructQPointsOctree( int qPtsStartID, int qPtsEndID, QPOINT *qPts, QPOINTS_OCTREE_NODE *qPointsOctree )
{
   int nodeID = nextFreeNode( );
   
   qPointsOctree[ nodeID ].qPtsStartID = qPtsStartID;
   qPointsOctree[ nodeID ].qPtsEndID = qPtsEndID;

   double minX = qPts[ qPtsStartID ].x, minY = qPts[ qPtsStartID ].y, minZ = qPts[ qPtsStartID ].z;
   double maxX = qPts[ qPtsStartID ].x, maxY = qPts[ qPtsStartID ].y, maxZ = qPts[ qPtsStartID ].z;
   
   for ( int i = qPtsStartID + 1; i <= qPtsEndID; i++ )
     {
      if ( qPts[ i ].x < minX ) minX = qPts[ i ].x;      
      if ( qPts[ i ].x > maxX ) maxX = qPts[ i ].x;      
      
      if ( qPts[ i ].y < minY ) minY = qPts[ i ].y;      
      if ( qPts[ i ].y > maxY ) maxY = qPts[ i ].y;      

      if ( qPts[ i ].z < minZ ) minZ = qPts[ i ].z;      
      if ( qPts[ i ].z > maxZ ) maxZ = qPts[ i ].z;      
     } 
   
   double cx = qPointsOctree[ nodeID ].cx = ( minX + maxX ) / 2;
   double cy = qPointsOctree[ nodeID ].cy = ( minY + maxY ) / 2;
   double cz = qPointsOctree[ nodeID ].cz = ( minZ + maxZ ) / 2;

   double r2 = ( qPts[ qPtsStartID ].x - cx ) * ( qPts[ qPtsStartID ].x - cx )
             + ( qPts[ qPtsStartID ].y - cy ) * ( qPts[ qPtsStartID ].y - cy )
             + ( qPts[ qPtsStartID ].z - cz ) * ( qPts[ qPtsStartID ].z - cz );
   
   for ( int i = qPtsStartID + 1; i <= qPtsEndID; i++ )
     {
      double r2T = ( qPts[ i ].x - cx ) * ( qPts[ i ].x - cx )
                 + ( qPts[ i ].y - cy ) * ( qPts[ i ].y - cy )
                 + ( qPts[ i ].z - cz ) * ( qPts[ i ].z - cz );
     
      if ( r2T > r2 ) r2 = r2T;
     } 
   
   double cr = qPointsOctree[ nodeID ].cr = sqrt( r2 );
         
   if ( ( qPtsEndID - qPtsStartID +  1 > maxLeafSize ) && ( cr > minRadius ) )   
     {
      qPointsOctree[ nodeID ].leaf = false;     
      
      int qPtsCount[ 8 ] = { 0, 0, 0, 0, 0, 0, 0, 0 };
      
      for ( int i = qPtsStartID; i <= qPtsEndID; i++ )
        {
         int j = ( zeroIfLess( qPts[ i ].z, cz ) << 2 )
               + ( zeroIfLess( qPts[ i ].y, cy ) << 1 )
               + ( zeroIfLess( qPts[ i ].x, cx ) );
         
         qPtsCount[ j ]++;
        }

      int qPtsStartIndex[ 8 ];

      qPtsStartIndex[ 0 ] = qPtsStartID;
      for ( int i = 1; i < 8; i++ )
        qPtsStartIndex[ i ] = qPtsStartIndex[ i - 1 ] + qPtsCount[ i - 1 ];

      for ( int i = 0; i < 8; i++ ) 
        if ( qPtsCount[ i ] > 0 )
          {
           int j = constructQPointsOctree( qPtsStartIndex[ i ], qPtsStartIndex[ i ] + qPtsCount[ i ] - 1, qPts, qPointsOctree );

           qPointsOctree[ nodeID ].cPtr[ i ] = j; 
          }
        else qPointsOctree[ nodeID ].cPtr[ i ] = -1;           
     }  
   else qPointsOctree[ nodeID ].leaf = true;
     
   return nodeID;  
}



bool pseudoGsol::buildQPointsOctree( int nQPoints, QPOINT *qPoints, int *numQPointsOctreeNodes, QPOINTS_OCTREE_NODE **qPointsOctree, int *qPointsOctreeRoot )
{  
   QPOINT *qPointsT;         
   qPointsT = ( QPOINT * ) malloc( nQPoints * sizeof( QPOINT ) );
   
   if ( qPointsT == NULL )
     {
      printError( (char *)"Failed to allocate temporary memory for quadrature points!" );
      return false;
     }

   countQPointsOctreeNodesAndSortQPoints( qPoints, 0, nQPoints - 1, qPointsT, numQPointsOctreeNodes );

   freeMem( qPointsT );
   
   ( *qPointsOctree ) = ( QPOINTS_OCTREE_NODE * ) malloc( ( *numQPointsOctreeNodes ) * sizeof( QPOINTS_OCTREE_NODE ) );
   
   if ( *qPointsOctree == NULL )
     {
      printError( (char *)"Unable to build quadrature points octree - memory allocation failed!" );
      return false;
     }
 
   initFreeNodeServer( nQPoints );
   
   *qPointsOctreeRoot = constructQPointsOctree( 0, nQPoints - 1, qPoints, *qPointsOctree );
   
   return true;
}



void pseudoGsol::assignHydrophobicityToQPoints( ATOMS_OCTREE_NODE *atomsOctree, int nodeA, ATOM *atoms,
                                                QPOINTS_OCTREE_NODE *qPointsOctree, int nodeQ, QPOINT *qPoints,
                                                double farDist, double rangeExt )
{
   double sumRad = atomsOctree[ nodeA ].cr + qPointsOctree[ nodeQ ].cr;
   double maxD2 = ( sumRad + farDist ) * ( sumRad + farDist );
   double dx = qPointsOctree[ nodeQ ].cx - atomsOctree[ nodeA ].cx,
          dy = qPointsOctree[ nodeQ ].cy - atomsOctree[ nodeA ].cy,
          dz = qPointsOctree[ nodeQ ].cz - atomsOctree[ nodeA ].cz;
   double d2 = dx * dx + dy * dy + dz * dz;
   
   if ( d2 <= maxD2 )
      {
       if ( atomsOctree[ nodeA ].leaf && qPointsOctree[ nodeQ ].leaf )
         {            
           for ( int i = atomsOctree[ nodeA ].atomsStartID; i <= atomsOctree[ nodeA ].atomsEndID; i++ )
            {
              double r2 = ( atoms[ i ].r + rangeExt ) * ( atoms[ i ].r + rangeExt );
              
              for ( int j = qPointsOctree[ nodeQ ].qPtsStartID; j <= qPointsOctree[ nodeQ ].qPtsEndID; j++ )
                {                            
                 dx = qPoints[ j ].x - atoms[ i ].x;
                 dy = qPoints[ j ].y - atoms[ i ].y;
                 dz = qPoints[ j ].z - atoms[ i ].z;                     
                        
                 d2 = dx * dx + dy * dy + dz * dz;       
                 
                 if ( d2 < r2 ) qPoints[ j ].h += atoms[ i ].h;
                }                 
            }   
         }
       else if ( !atomsOctree[ nodeA ].leaf && !qPointsOctree[ nodeQ ].leaf )         
              {
                for ( int i = 0; i < 8; i++ )
                  if ( atomsOctree[ nodeA ].cPtr[ i ] >= 0 )
                     for ( int j = 0; j < 8; j++ )
                       if ( qPointsOctree[ nodeQ ].cPtr[ j ] >= 0 ) 
                          assignHydrophobicityToQPoints( atomsOctree, atomsOctree[ nodeA ].cPtr[ i ], atoms,
                                                         qPointsOctree, qPointsOctree[ nodeQ ].cPtr[ j ], qPoints,
                                                         farDist, rangeExt );                         
              }
            else if ( !atomsOctree[ nodeA ].leaf )         
                   {
                     for ( int i = 0; i < 8; i++ )
                       if ( atomsOctree[ nodeA ].cPtr[ i ] >= 0 )
                          assignHydrophobicityToQPoints( atomsOctree, atomsOctree[ nodeA ].cPtr[ i ], atoms,
                                                         qPointsOctree, nodeQ, qPoints, farDist, rangeExt );                         
                   }
                 else 
                   {
                     for ( int j = 0; j < 8; j++ )
                       if ( qPointsOctree[ nodeQ ].cPtr[ j ] >= 0 )
                          assignHydrophobicityToQPoints( atomsOctree, nodeA, atoms, 
                                                         qPointsOctree, qPointsOctree[ nodeQ ].cPtr[ j ], qPoints, farDist, rangeExt );                         
                   }              
      }      
}


inline void pseudoGsol::transformPoint( double x, double y, double z, double *transMat, double *nx, double *ny, double *nz )
{
   *nx = transMat[  0 ] * x + transMat[  1 ] * y + transMat[  2 ] * z + transMat[  3 ];
   *ny = transMat[  4 ] * x + transMat[  5 ] * y + transMat[  6 ] * z + transMat[  7 ];
   *nz = transMat[  8 ] * x + transMat[  9 ] * y + transMat[ 10 ] * z + transMat[ 11 ];       
}


inline void pseudoGsol::transformPoint( Point p, double *transMat, Point *np )
{
   np->x = transMat[  0 ] * p.x + transMat[  1 ] * p.y + transMat[  2 ] * p.z + transMat[  3 ];
   np->y = transMat[  4 ] * p.x + transMat[  5 ] * p.y + transMat[  6 ] * p.z + transMat[  7 ];
   np->z = transMat[  8 ] * p.x + transMat[  9 ] * p.y + transMat[ 10 ] * p.z + transMat[ 11 ];       
}



double pseudoGsol::computeXlateForPG( int numStQPoints, QPOINT *stQPoints, int numMvQPoints, QPOINT *mvQPoints )
{
   double minXYZ;
   
   for ( int i = 0; i < numStQPoints; i++ )
     {
       if ( i == 0 ) minXYZ = stQPoints[ i ].x;
       else if ( stQPoints[ i ].x < minXYZ ) minXYZ = stQPoints[ i ].x;
       
       if ( stQPoints[ i ].y < minXYZ ) minXYZ = stQPoints[ i ].y;
       if ( stQPoints[ i ].z < minXYZ ) minXYZ = stQPoints[ i ].z;
     }

   double minX, minY, minZ;
   double maxX, maxY, maxZ;

   for ( int i = 0; i < numMvQPoints; i++ )
     {
       if ( i == 0 ) 
         {
           minX = maxX = mvQPoints[ i ].x;
           minY = maxY = mvQPoints[ i ].y;
           minZ = maxZ = mvQPoints[ i ].z;           
         }  
       else  
         {
           if ( mvQPoints[ i ].x < minX ) minX = mvQPoints[ i ].x;
           if ( mvQPoints[ i ].y < minY ) minY = mvQPoints[ i ].y;
           if ( mvQPoints[ i ].z < minZ ) minZ = mvQPoints[ i ].z;                      

           if ( mvQPoints[ i ].x > maxX ) maxX = mvQPoints[ i ].x;
           if ( mvQPoints[ i ].y > maxY ) maxY = mvQPoints[ i ].y;
           if ( mvQPoints[ i ].z > maxZ ) maxZ = mvQPoints[ i ].z;                      
         }         
     }

   double maxD = 0, d;

   d = fabs( maxX - minX );
   if ( d > maxD ) maxD = d; 

   d = fabs( maxY - minY );
   if ( d > maxD ) maxD = d; 
   
   d = fabs( maxZ - minZ );
   if ( d > maxD ) maxD = d; 

   minXYZ -= ( 2 * maxD );
   
   if ( minXYZ < 0 ) return -minXYZ;
   else return 0;  
}



void pseudoGsol::printGsolParamters( FILE* fp )
{
  fprintf( fp, (char *)"# \t staticMoleculeQUAD = %s\n", params.staticMoleculeQUAD );
  fprintf( fp, (char *)"# \t movingMoleculeQUAD = %s\n", params.movingMoleculeQUAD );  
  
  fprintf( fp, (char *)"# \t distanceCutoff = %lf\n", params.distanceCutoff );      
}


double pseudoGsol::determinant( double *trans )
{
  return ( - trans[ LID( 0, 2 ) ] * trans[ LID( 1, 1 ) ] * trans[ LID( 2, 0 ) ]
           + trans[ LID( 0, 1 ) ] * trans[ LID( 1, 2 ) ] * trans[ LID( 2, 0 ) ]
           + trans[ LID( 0, 2 ) ] * trans[ LID( 1, 0 ) ] * trans[ LID( 2, 1 ) ]
           - trans[ LID( 0, 0 ) ] * trans[ LID( 1, 2 ) ] * trans[ LID( 2, 1 ) ]
           - trans[ LID( 0, 1 ) ] * trans[ LID( 1, 0 ) ] * trans[ LID( 2, 2 ) ]
           + trans[ LID( 0, 0 ) ] * trans[ LID( 1, 1 ) ] * trans[ LID( 2, 2 ) ] );
}


void pseudoGsol::invert( double *trans, double *transI )
{
  double d = determinant( trans );
  
  if ( d != 0.0 )
     {
       transI[ LID( 0, 0 ) ] = ( trans[ LID( 1, 1 ) ] * trans[ LID( 2, 2 ) ]
                               - trans[ LID( 1, 2 ) ] * trans[ LID( 2, 1 ) ] ) / d;
                               
       transI[ LID( 0, 1 ) ] = ( trans[ LID( 0, 2 ) ] * trans[ LID( 2, 1 ) ] 
                               - trans[ LID( 0, 1 ) ] * trans[ LID( 2, 2 ) ] ) / d;
                               
       transI[ LID( 0, 2 ) ] = ( trans[ LID( 0, 1 ) ] * trans[ LID( 1, 2 ) ]
                               - trans[ LID( 0, 2 ) ] * trans[ LID( 1, 1 ) ] ) / d;
                               
       transI[ LID( 0, 3 ) ] = ( trans[ LID( 0, 3 ) ] * trans[ LID( 1, 2 ) ] * trans[ LID( 2, 1 ) ] 
                               - trans[ LID( 0, 2 ) ] * trans[ LID( 1, 3 ) ] * trans[ LID( 2, 1 ) ] 
                               - trans[ LID( 0, 3 ) ] * trans[ LID( 1, 1 ) ] * trans[ LID( 2, 2 ) ] 
                               + trans[ LID( 0, 1 ) ] * trans[ LID( 1, 3 ) ] * trans[ LID( 2, 2 ) ] 
                               + trans[ LID( 0, 2 ) ] * trans[ LID( 1, 1 ) ] * trans[ LID( 2, 3 ) ] 
                               - trans[ LID( 0, 1 ) ] * trans[ LID( 1, 2 ) ] * trans[ LID( 2, 3 ) ] ) / d;

       transI[ LID( 1, 0 ) ] = ( trans[ LID( 1, 2 ) ] * trans[ LID( 2, 0 ) ] 
                               - trans[ LID( 1, 0 ) ] * trans[ LID( 2, 2 ) ] ) / d;
                               
       transI[ LID( 1, 1 ) ] = ( trans[ LID( 0, 0 ) ] * trans[ LID( 2, 2 ) ]
                               - trans[ LID( 0, 2 ) ] * trans[ LID( 2, 0 ) ] ) / d;
                               
       transI[ LID( 1, 2 ) ] = ( trans[ LID( 0, 2 ) ] * trans[ LID( 1, 0 ) ]
                               - trans[ LID( 0, 0 ) ] * trans[ LID( 1, 2 ) ] ) / d;
                               
       transI[ LID( 1, 3 ) ] = ( trans[ LID( 0, 2 ) ] * trans[ LID( 1, 3 ) ] * trans[ LID( 2, 0 ) ] 
                               - trans[ LID( 0, 3 ) ] * trans[ LID( 1, 2 ) ] * trans[ LID( 2, 0 ) ] 
                               + trans[ LID( 0, 3 ) ] * trans[ LID( 1, 0 ) ] * trans[ LID( 2, 2 ) ] 
                               - trans[ LID( 0, 0 ) ] * trans[ LID( 1, 3 ) ] * trans[ LID( 2, 2 ) ] 
                               - trans[ LID( 0, 2 ) ] * trans[ LID( 1, 0 ) ] * trans[ LID( 2, 3 ) ] 
                               + trans[ LID( 0, 0 ) ] * trans[ LID( 1, 2 ) ] * trans[ LID( 2, 3 ) ] ) / d;

       transI[ LID( 2, 0 ) ] = ( trans[ LID( 1, 0 ) ] * trans[ LID( 2, 1 ) ]
                               - trans[ LID( 1, 1 ) ] * trans[ LID( 2, 0 ) ] ) / d;
                               
       transI[ LID( 2, 1 ) ] = ( trans[ LID( 0, 1 ) ] * trans[ LID( 2, 0 ) ]
                               - trans[ LID( 0, 0 ) ] * trans[ LID( 2, 1 ) ] ) / d;
                               
       transI[ LID( 2, 2 ) ] = ( trans[ LID( 0, 0 ) ] * trans[ LID( 1, 1 ) ]
                               - trans[ LID( 0, 1 ) ] * trans[ LID( 1, 0 ) ] ) / d;
                               
       transI[ LID( 2, 3 ) ] = ( trans[ LID( 0, 3 ) ] * trans[ LID( 1, 1 ) ] * trans[ LID( 2, 0 ) ] 
                               - trans[ LID( 0, 1 ) ] * trans[ LID( 1, 3 ) ] * trans[ LID( 2, 0 ) ] 
                               - trans[ LID( 0, 3 ) ] * trans[ LID( 1, 0 ) ] * trans[ LID( 2, 1 ) ] 
                               + trans[ LID( 0, 0 ) ] * trans[ LID( 1, 3 ) ] * trans[ LID( 2, 1 ) ] 
                               + trans[ LID( 0, 1 ) ] * trans[ LID( 1, 0 ) ] * trans[ LID( 2, 3 ) ] 
                               - trans[ LID( 0, 0 ) ] * trans[ LID( 1, 1 ) ] * trans[ LID( 2, 3 ) ] ) / d;
     }
  else 
     {
       for ( int i = 0; i < 12; i++ )
         transI[ i ] = ( ( i % 5 ) ? 0.0 : 1.0 );
     }   
}


void pseudoGsol::initOctreeFlags( int threadID )
{
   int offset = 2 * threadID * numStaticQPointsOctreeNodes;
   
   for ( int i = 0; i < numStaticQPointsOctreeNodes; i++ )
     staticQPointsOctreeFlags[ offset + i ] = false;
     
   offset = 2 * threadID * numMovingQPointsOctreeNodes;
   
   for ( int i = 0; i < numMovingQPointsOctreeNodes; i++ )
     movingQPointsOctreeFlags[ offset + i ] = false;     
}


void pseudoGsol::collectPseudoGsol( int threadID, double *trans, double *transI, double *pGsol, 
                                    double *pGsolHStaticPos, double *pGsolHStaticNeg, double *pGsolHMovingPos, double *pGsolHMovingNeg )
{
  
 *pGsol = *pGsolHStaticPos = *pGsolHStaticNeg = *pGsolHMovingPos = *pGsolHMovingNeg = 0;

	int offsetStatic = 2 * threadID * numStaticQPointsOctreeNodes;   
	int offsetMoving = 2 * threadID * numMovingQPointsOctreeNodes;
	
		call_to_kern(offsetStatic, transI, pGsol, pGsolHStaticPos, pGsolHStaticNeg, staticQPointsOctreeFlags,staticQPointsOctree, staticQPoints, numStaticQPoints, numStaticQPointsOctreeNodes, movingPG, movingPG->TRANSLATE, movingPG->getdivsize(), movingPG->getRangeCount(), params, offsetMoving, trans, pGsolHMovingPos, pGsolHMovingNeg, movingQPointsOctreeFlags, movingQPointsOctree, movingQPoints, numMovingQPoints,  numMovingQPointsOctreeNodes, staticPG, staticPG->TRANSLATE, staticPG->getdivsize(), staticPG->getRangeCount());



/*
   //int offset = 2 * threadID * numStaticQPointsOctreeNodes;   
   
   for ( int i = 0; i < numStaticQPointsOctreeNodes; i++ )
		{
     if ( staticQPointsOctreeFlags[ offsetStatic + i ] )
        {
           for ( int j = staticQPointsOctree[ i ].qPtsStartID; j <= staticQPointsOctree[ i ].qPtsEndID; j++ )
             {
               Point p, q;
               
               p.x = staticQPoints[ j ].x;
               p.y = staticQPoints[ j ].y;
               p.z = staticQPoints[ j ].z;
               
               transformPoint( p, transI, &q );
               
               if ( movingPG->pointsWithinRange( &q, params.distanceCutoff ) )
                 {
                   ( *pGsol ) += staticQPoints[ j ].w;
                   if ( staticQPoints[ j ].h > 0 ) ( *pGsolHStaticPos ) += staticQPoints[ j ].h * staticQPoints[ j ].w;
                   else ( *pGsolHStaticNeg ) += staticQPoints[ j ].h * staticQPoints[ j ].w;
                 }  
             }            
        } 
	}
     
   //offset = 2 * threadID * numMovingQPointsOctreeNodes;
   
   for ( int i = 0; i < numMovingQPointsOctreeNodes; i++ ){
     if ( movingQPointsOctreeFlags[ offsetMoving + i ] )
        {
           for ( int j = movingQPointsOctree[ i ].qPtsStartID; j <= movingQPointsOctree[ i ].qPtsEndID; j++ )
             {
               Point p, q;
               
               p.x = movingQPoints[ j ].x;
               p.y = movingQPoints[ j ].y;
               p.z = movingQPoints[ j ].z;
               
               transformPoint( p, trans, &q );
               
               if ( staticPG->pointsWithinRange( &q, params.distanceCutoff ) )
                 {
                   ( *pGsol ) += movingQPoints[ j ].w;
                   if ( movingQPoints[ j ].h > 0 ) ( *pGsolHMovingPos ) += movingQPoints[ j ].h * movingQPoints[ j ].w;
                   else ( *pGsolHMovingNeg ) += movingQPoints[ j ].h * movingQPoints[ j ].w;
                 }  
             }         
        }   
      
}
*/
	printf("CPU *pGsol, *pGsolHStaticPos, *pGsolHStaticNeg, *pGsolHMovingPos, *pGsolHMovingNeg %lf %lf %lf %lf %lf\n", *pGsol, *pGsolHStaticPos, *pGsolHStaticNeg, *pGsolHMovingPos, *pGsolHMovingNeg);

}




void pseudoGsol::markPotentialQPoints( int threadID, int nodeS, int nodeM, double *trans )
{
   double sumRad = staticQPointsOctree[ nodeS ].cr + movingQPointsOctree[ nodeM ].cr;
   double maxD2 = ( sumRad + params.distanceCutoff ) * ( sumRad + params.distanceCutoff );
   Point P, Q;
   
   P.x = movingQPointsOctree[ nodeM ].cx;
   P.y = movingQPointsOctree[ nodeM ].cy;
   P.z = movingQPointsOctree[ nodeM ].cz;
   
   transformPoint( P, trans, &Q );
   
   double dx = staticQPointsOctree[ nodeS ].cx - Q.x,
          dy = staticQPointsOctree[ nodeS ].cy - Q.y,
          dz = staticQPointsOctree[ nodeS ].cz - Q.z;
   double d2 = dx * dx + dy * dy + dz * dz;

   if ( d2 <= maxD2 )
      {
       if ( staticQPointsOctree[ nodeS ].leaf && movingQPointsOctree[ nodeM ].leaf )
         {            
           staticQPointsOctreeFlags[ 2 * threadID * numStaticQPointsOctreeNodes + nodeS ] = true;
           movingQPointsOctreeFlags[ 2 * threadID * numMovingQPointsOctreeNodes + nodeM ] = true;
         }
       else if ( !staticQPointsOctree[ nodeS ].leaf && !movingQPointsOctree[ nodeM ].leaf )         
              {
                for ( int i = 0; i < 8; i++ )
                  if ( staticQPointsOctree[ nodeS ].cPtr[ i ] >= 0 )
                     for ( int j = 0; j < 8; j++ )
                       if ( movingQPointsOctree[ nodeM ].cPtr[ j ] >= 0 ) 
                           markPotentialQPoints( threadID, staticQPointsOctree[ nodeS ].cPtr[ i ], movingQPointsOctree[ nodeM ].cPtr[ j ], trans );                         
              }
            else if ( !staticQPointsOctree[ nodeS ].leaf )         
                   {
                     for ( int i = 0; i < 8; i++ )
                       if ( staticQPointsOctree[ nodeS ].cPtr[ i ] >= 0 )
                           markPotentialQPoints( threadID, staticQPointsOctree[ nodeS ].cPtr[ i ], nodeM, trans );                         
                   }
                 else 
                   {
                     for ( int j = 0; j < 8; j++ )
                       if ( movingQPointsOctree[ nodeM ].cPtr[ j ] >= 0 )
                           markPotentialQPoints( threadID, nodeS, movingQPointsOctree[ nodeM ].cPtr[ j ], trans );                         
                   }              
      }      
}



void pseudoGsol::getPseudoGsol( int threadID, double *trans, double *pGsol,
                                double *pGsolHStaticPos, double *pGsolHStaticNeg, double *pGsolHMovingPos, double *pGsolHMovingNeg )
{
   if ( ( threadID < 0 ) || ( threadID >= numThreads ) )
     {
       printError( (char *)"Invalid thread id!" );
       exit( 1 );
     }
   
   initOctreeFlags( threadID );   
   
   markPotentialQPoints( threadID, staticQPointsOctreeRoot, movingQPointsOctreeRoot, trans );

   double transI[ 12 ];
   
   invert( trans, transI );  
   
   collectPseudoGsol( threadID, trans, transI, pGsol, pGsolHStaticPos, pGsolHStaticNeg, pGsolHMovingPos, pGsolHMovingNeg );
}


void pseudoGsol::getPseudoGsol( int threadID, double *pGsol, double *pGsolHStaticPos, double *pGsolHStaticNeg, double *pGsolHMovingPos, double *pGsolHMovingNeg )
{
  double trans[ ] = { 1, 0, 0, 0,
                      0, 1, 0, 0,
                      0, 0, 1, 0 };

  getPseudoGsol( threadID, trans, pGsol, pGsolHStaticPos, pGsolHStaticNeg, pGsolHMovingPos, pGsolHMovingNeg );   
}


void pseudoGsol::getPseudoGsol( int threadID, double *trans, double *pGsol, double *pGsolH )
{
   double pGsolHStaticPos, pGsolHStaticNeg, pGsolHMovingPos, pGsolHMovingNeg;
   
   getPseudoGsol( threadID, trans, pGsol, &pGsolHStaticPos, &pGsolHStaticNeg, &pGsolHMovingPos, &pGsolHMovingNeg );
   
   ( *pGsolH ) = pGsolHStaticPos + pGsolHStaticNeg + pGsolHMovingPos + pGsolHMovingNeg;
}


void pseudoGsol::getPseudoGsol( int threadID, double *pGsol, double *pGsolH )
{
  double trans[ ] = { 1, 0, 0, 0,
                      0, 1, 0, 0,
                      0, 0, 1, 0 };

  getPseudoGsol( threadID, trans, pGsol, pGsolH );
}



bool pseudoGsol::getParamsFromFile( PARAMS_IN *p, char *paramFile, bool atomsFromFile )
{
  char s[ 2000 ];
  char key[ 500 ], val[ 500 ];
  FILE *fp;

  fp = fopen( paramFile, "r" );

  if ( fp == NULL )
    {
      printError( (char *)"Failed to open parameter file %s!", paramFile );
      return false;
    }

  p->staticMoleculePQR = p->staticMoleculePQR = NULL;
  p->staticMoleculeQUAD = p->staticMoleculeQUAD = NULL;
  p->distanceCutoff = 1.5;
  p->useInterfacePropensity = true;  
  p->perResidueHydrophobicity = true;
  
  while ( fgets( s, 1999, fp ) != NULL )
    {
      if ( sscanf( s, (char *)"%s %s", key, val ) != 2 ) continue;

      if ( !strcasecmp( key, (char *)"staticMoleculePQR" ) ) p->staticMoleculePQR = strdup( val );
      else if ( !strcasecmp( key, (char *)"movingMoleculePQR" ) ) p->movingMoleculePQR = strdup( val );    
      else if ( !strcasecmp( key, (char *)"staticMoleculeQUAD" ) ) p->staticMoleculeQUAD = strdup( val );
      else if ( !strcasecmp( key, (char *)"movingMoleculeQUAD" ) ) p->movingMoleculeQUAD = strdup( val );
      else if ( !strcasecmp( key, (char *)"distanceCutoff" ) )
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
      else if ( !strcasecmp(key, (char *)"useInterfacePropensity" ) ) 
             {
	       if ( !strcasecmp( val, (char *)"true" ) ) p->useInterfacePropensity = true;
	       else if ( !strcasecmp( val, (char *)"false" ) ) p->useInterfacePropensity = false;
	       else {
		      printf( (char *)"Error: %s must be a Boolean value!\n", key);
		      return false;
	            }
	     }                          
      else if ( !strcasecmp(key, (char *)"perResidueHydrophobicity" ) ) 
             {
	       if ( !strcasecmp( val, (char *)"true" ) ) p->perResidueHydrophobicity = true;
	       else if ( !strcasecmp( val, (char *)"false" ) ) p->perResidueHydrophobicity = false;
	       else {
		      printf( (char *)"Error: %s must be a Boolean value!\n", key);
		      return false;
	            }
	     }             
      else if ( !strcasecmp( key, (char *)"numThreads" ) )
             {
               int v = atoi( val );
               
               if ( v < 1 )
                 {
                   printError( (char *)"%s must be a positive integer!", key );
                   fclose( fp );
                   return false;
                 }
                 
               p->numThreads = v;  
             }	     
    }

  fclose( fp );

  if ( atomsFromFile )    
    {
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
    
  return true;
}
