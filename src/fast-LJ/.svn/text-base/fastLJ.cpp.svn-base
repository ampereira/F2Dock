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

#include "fastLJ.h"

#include <iostream>

using namespace std;


void fastLJ::freeMemory( void )
{
   for ( int i = 0; i < numAtomType; i++ )
     {
       freeMem( staticAtoms[ i ] );
       freeMem( movingAtoms[ i ] );
       freeMem( movingAtomsOrig[ i ] );
     }  
     
   freeMem( staticAtomsOctree );
   freeMem( movingAtomsOctree );  
}


void fastLJ::computeConstants( void )
{
  double r_eqm_X[ ] = { 2.00, 1.00, 1.75, 1.60, 2.10, 2.00 };
  double eps_X[ ] = { 0.15, 0.02, 0.16, 0.20, 0.20, 0.20 };
  double eps_XY[ numAtomType ][ numAtomType ];

  for ( int i = 0; i < numAtomType; i++ )
    r_eqm_X[ i ] *= vdWEqmRadScale;
  
  for ( int i = 0, k = 0; i < numAtomType; i++ )
     for ( int j = 0; j < numAtomType; j++, k++ )
        {
         eps_XY[ i ][ j ] = sqrt( eps_X[ i ] * eps_X[ j ] );         
         r_eqm_XY[ i ][ j ] = ( r_eqm_X[ i ] + r_eqm_X[ j ] ) / 2.0;         
         A[ i ][ j ] = eps_XY[ i ][ j ] * pow( r_eqm_XY[ i ][ j ], 12.0 );
         B[ i ][ j ] = 2 * eps_XY[ i ][ j ] * pow( r_eqm_XY[ i ][ j ], 6.0 );
        }                               
}


void fastLJ::setDefaults( void )
{
   for ( int i = 0; i < numAtomType; i++ )
     { 
       staticAtoms[ i ] = NULL;
       movingAtoms[ i ] = NULL;
       movingAtomsOrig[ i ] = NULL;
     }  
      
   staticAtomsOctree = NULL;
   movingAtomsOctree = NULL;
   
   minRadius = 2.0;
   minRadiusUsed = -1;      
   
   maxLeafSize = 20;
   maxLeafSizeUsed = -1; 
   
   minInterAtomDist = 0.4;  
   minInterAtomDistUsed = -1;         

   epsilon = 0.5;
   
   staticAtomsOctreeBuilt = false;
   movingAtomsOctreeBuilt = false;
   
   curNode = maxNode = 0;
   
   staticAtomsOctreeRoot = -1;
   movingAtomsOctreeRoot = -1;

   vdWEqmRadScale = 0.3;
   
   numThreads = 1;

   for ( int i = 0; i < 16; i++ )
      transMatrix[ i ] = ( i % 5 ) ? 0.0 : 1.0;
     
   printStatus = true;
   
   useSSEFunctions = false;
   
   minD2 = minInterAtomDist * minInterAtomDist;

#ifdef USE_SSE
   MIND2 = _mm_set1_ps( ( float ) minD2 );   
#endif   

   computeConstants( );   
}


void fastLJ::printCurrentSettings( void )
{
   printf( "\nCurrent Parameter Settings:\n" );
   printf( "\tminRadius = %lf, maxLeafSize = %d, minInterAtomDist = %lf, epsilon = %lf, numThreads = %d\n", minRadius, maxLeafSize, minInterAtomDist, epsilon, numThreads );
#ifdef USE_SSE
   if ( useSSEFunctions )
     printf( "\tuse intrinsic SSE functions: YES\n" );
   else  
     printf( "\tuse intrinsic SSE functions: NO\n" );   
#else   
   printf( "\tuse intrinsic SSE functions: NO\n" );
#endif   
}


fastLJ::fastLJ( int numStaticAtoms, double *stAtomsXYZ, char *stAtomsType, int numMovingAtoms, double *mvAtomsXYZ, char *mvAtomsType, int nThreads, bool printStat )
{
   setDefaults( );

   if ( nThreads < 1 )
     {
       printError( (char *)"Invalid number of threads ( %d )!", nThreads );
       nThreads = 1;
     }
   
   numThreads = nThreads;
   printStatus = printStat;
         
   if ( !copyAtomsFromArray( numStaticAtoms, stAtomsXYZ, stAtomsType, nStaticAtoms, staticAtoms ) 
     || !copyAtomsFromArray( numMovingAtoms, mvAtomsXYZ, mvAtomsType, nMovingAtoms, movingAtomsOrig ) )
      {
       freeMemory( );
       exit( 1 );
      }

   for ( int i = 0; i < numAtomType; i++ )
     {
       movingAtoms[ i ] = ( ATOM * ) malloc( numThreads * nMovingAtoms[ i ] * sizeof( ATOM ) );

       if ( movingAtoms[ i ] == NULL )
         {
          printError( (char *)"Failed to allocate memory for atoms!" );
          freeMemory( );
          exit( 1 );
         }
     }    
 
   buildOctrees( );
  
   if ( printStatus ) printCurrentSettings( );
}


fastLJ::fastLJ( char *staticAtomsFile, char *movingAtomsFile, int nThreads, bool printStat )
{
   setDefaults( );

   if ( nThreads < 1 )
     {
       printError( (char *)"Invalid number of threads ( %d )!", nThreads );
       nThreads = 1;
     }
   
   numThreads = nThreads;
   printStatus = printStat;
         
   if ( !readAtomsFromFile( staticAtomsFile, nStaticAtoms, staticAtoms ) 
     || !readAtomsFromFile( movingAtomsFile, nMovingAtoms, movingAtomsOrig ) )
      {
       freeMemory( );
       exit( 1 );
      }

   for ( int i = 0; i < numAtomType; i++ )
     {
       movingAtoms[ i ] = ( ATOM * ) malloc( numThreads * nMovingAtoms[ i ] * sizeof( ATOM ) );

       if ( movingAtoms[ i ] == NULL )
         {
          printError( (char *)"Failed to allocate memory for atoms!" );
          freeMemory( );
          exit( 1 );
         }
     }    
            
   buildOctrees( );
  
   if ( printStatus ) printCurrentSettings( );
}


fastLJ::fastLJ( char *paramFile )
{
   setDefaults( );
   
   char *staticAtomsFile;
   char *movingAtomsFile;
      
   if ( !getParamsFromFile( paramFile, &staticAtomsFile, &movingAtomsFile )
     || !readAtomsFromFile( staticAtomsFile, nStaticAtoms, staticAtoms ) 
     || !readAtomsFromFile( movingAtomsFile, nMovingAtoms, movingAtomsOrig ) )
      {
       freeMemory( );
       exit( 1 );
      }
      
   freeMem( staticAtomsFile );   
   freeMem( movingAtomsFile );      

   for ( int i = 0; i < numAtomType; i++ )
     {
       movingAtoms[ i ] = ( ATOM * ) malloc( numThreads * nMovingAtoms[ i ] * sizeof( ATOM ) );

       if ( movingAtoms[ i ] == NULL )
         {
          printError( (char *)"Failed to allocate memory for atoms!" );
          freeMemory( );
          exit( 1 );
         }
     }    
 
   buildOctrees( );
  
   if ( printStatus ) printCurrentSettings( );
}


fastLJ::~fastLJ( )
{
   freeMemory( );
   if ( printStatus ) printf( "\n" );
}



bool fastLJ::setMinRadius( double minRad )
{
   if ( minRad < 0 )
     {
      printError( (char *)"minRadius must be a non-negative real number!" );
      return false;     
     }
     
   minRadius = minRad;

   bool staticAtomsLoaded = false, movingAtomsLoaded = false;
   
   for ( int k = 0; k < numAtomType; k++ )
     {
       if ( staticAtoms[ k ] != NULL ) staticAtomsLoaded = true;
       if ( movingAtomsOrig[ k ] != NULL ) movingAtomsLoaded = true;
     }

   if ( ( staticAtomsLoaded == true ) && ( movingAtomsLoaded == true ) ) 
      buildOctrees( );

   if ( printStatus ) printf( "\nminRadius is set to %lf\n", minRad );
   
   return true;
}


bool fastLJ::setMaxLeafSize( int maxLfSize )
{
   if ( maxLfSize <= 0 )
     {
      printError( (char *)"maxLeafSize must be a positive integer!" );
      return false;     
     }
     
   maxLeafSize = maxLfSize;

   bool staticAtomsLoaded = false, movingAtomsLoaded = false;
   
   for ( int k = 0; k < numAtomType; k++ )
     {
       if ( staticAtoms[ k ] != NULL ) staticAtomsLoaded = true;
       if ( movingAtomsOrig[ k ] != NULL ) movingAtomsLoaded = true;
     }

   if ( ( staticAtomsLoaded == true ) && ( movingAtomsLoaded == true ) ) 
      buildOctrees( );

   if ( printStatus ) printf( "\nmaxLeafSize is set to %d\n", maxLfSize );

   return true;
}


bool fastLJ::setMinInterAtomDist( double minAtomDist )
{
   if ( minAtomDist <= 0 )
     {
      printError( (char *)"minInterAtomDist must be a positive real number!" );
      return false;     
     }
     
   minInterAtomDist = minAtomDist;
   minD2 =  minInterAtomDist * minInterAtomDist;
   
#ifdef USE_SSE
   MIND2 = _mm_set1_ps( ( float ) minD2 );   
#endif   

   if ( ( staticAtoms != NULL ) && ( movingAtomsOrig != NULL ) ) 
      buildOctrees( );

   if ( printStatus ) printf( "\nminInterAtomDist is set to %lf\n", minAtomDist );
   
   return true;
}


bool fastLJ::setEpsilon( double eps )
{
   if ( eps < 0 )
     {
      printError( (char *)"epsilon must be a non-negative real number!" );
      return false;     
     }
     
   epsilon = eps;

   if ( printStatus ) printf( "\nepsilon is set to %lf\n", eps );

   return true;
}



void fastLJ::useSSE( bool useSSEFuncs )
{
#ifdef USE_SSE
   useSSEFunctions = useSSEFuncs;
   
   if ( printStatus ) 
     {
       if ( useSSEFuncs ) printf( "\nSSE functions will be used\n" );
       else printf( "\nSSE functions will not be used\n" );
     }  
#else
   if ( printStatus && useSSEFuncs ) printf( "\nplease recompile with sse=true\n" );
   useSSEFunctions = false;  
#endif     
}


void fastLJ::setTransformationMatrix( double *transMat )
{
   for ( int i = 0; i < 16; i++ )
      transMatrix[ i ] = transMat[ i ];
   
   if ( printStatus ) 
     {
       printf( "\nTransformation Matrix is set to:\n" );
       for ( int i = 0; i < 16; i++ )
          printf( "%s%8.3lf%s", transMatrix[ i ], ( i % 4 ) ? "" : "\t", ( ( i + 1 ) % 4 ) ? "  " : "\n" );
       printf( "\n" );   
     }  
}


void fastLJ::setPrintStatus( bool printStat )
{
   printStatus = printStat;
   
   if ( printStatus ) printf( "\nprintStatus is set to true\n" );
}


void fastLJ::transformMovingAtoms( int threadID, double *transMat )
{
   for ( int k = 0; k < numAtomType; k++ )
     for ( int i = 0, j = threadID * nMovingAtoms[ k ]; i < nMovingAtoms[ k ]; i++ )
       {
         transform( movingAtomsOrig[ k ][ i ].x, movingAtomsOrig[ k ][ i ].y, movingAtomsOrig[ k ][ i ].z, 
                    transMat, movingAtoms[ k ][ j + i ].x, movingAtoms[ k ][ j + i ].y, movingAtoms[ k ][ j + i ].z );
         movingAtoms[ k ][ j + i ].id = movingAtomsOrig[ k ][ i ].id;           
       }             
}


int fastLJ::nextFreeNode( void )
{
   int nextNode = -1;
   
   if ( curNode < maxNode ) nextNode = curNode++;
   
   return nextNode;
}


void fastLJ::initFreeNodeServer( int numNodes )
{
   curNode = 0;
   maxNode = numNodes;
}


bool fastLJ::getParamsFromFile( char *paramFile, char **staticAtomsFile, char **movingAtomsFile )
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

  ( *staticAtomsFile ) = ( *movingAtomsFile ) = NULL;

  while ( fgets( s, 1999, fp ) != NULL )
    {
      if ( sscanf( s, "%s %s", key, val ) != 2 ) continue;
    
      if ( !strcasecmp( key, "staticMoleculePQR" ) || !strcasecmp( key, "staticMoleculePDB" ) ) ( *staticAtomsFile ) = strdup( val );
      else if ( !strcasecmp( key, "movingMoleculePQR" ) || !strcasecmp( key, "movingMoleculePDB" ) ) ( *movingAtomsFile ) = strdup( val );
      else if ( !strcasecmp( key, "numThreads" ) )
             {
               int v = atoi( val );               
               if ( v < 1 ) printError( (char *)"%s must be a positive integer!", key );
               else numThreads = v;  
             }
      else if ( !strcasecmp( key, "epsilonLJ" ) )
             {
               double v = atof( val );
               
               if ( v < 0 ) printError( (char *)"%s must be a non-negative float!", key );
               else setEpsilon( v );                 
             }
      else if ( !strcasecmp( key, "minRadius" ) )
             {
               double v = atof( val );
               
               if ( v < 0 ) printError( (char *)"%s must be a non-negative float!", key );
               else setMinRadius( v );                 
             }
      else if ( !strcasecmp( key, "maxLeafSize" ) )
             {
               int v = atoi( val );               
               if ( v < 1 ) printError( (char *)"%s must be a positive integer!", key );
               else setMaxLeafSize( v );
             }
      else if ( !strcasecmp( key, "minInterAtomDist" ) )
             {
               double v = atof( val );
               
               if ( v <= 0 ) printError( (char *)"%s must be a positive float!", key );
               else setMinInterAtomDist( v );
             }
      else if ( !strcasecmp( key, "useSSE" ) ) 
             {
	       if ( !strcasecmp( val, "true" ) ) useSSE( true );
	       else if ( !strcasecmp( val, "false" ) ) useSSE( false );
	            else printError( (char *)"%s must be a Boolean value!", key );
	     }             
      else if ( !strcasecmp( key, "vdWEqmRadScale" ) ) 
             {
	       double v = atof( val );
	       if ( v <= 0 ) printError( (char *)"%s must be a positive real value!", key );
	       else vdWEqmRadScale = v;
  	     }
      else if ( !strcasecmp( key, "printStatus" ) ) 
             {
	       if ( !strcasecmp( val, "true" ) ) setPrintStatus( true );
	       else if ( !strcasecmp( val, "false" ) ) setPrintStatus( false );
	            else printError( (char *)"%s must be a Boolean value!", key );
	     }             
    }

  fclose( fp );
    
  if ( ( ( *staticAtomsFile ) == NULL ) || ( ( *movingAtomsFile ) == NULL ) )
    { 
      printError( (char *)"Missing PDB/PQR file name for the %s molecule!", ( ( *staticAtomsFile ) == NULL ) ? "static" : "moving" );
      freeMem( ( *staticAtomsFile ) );
      freeMem( ( *movingAtomsFile ) );
      return false;
    }
        
  return true;
}


bool fastLJ::readAtomsFromPQR( char *atomsFile, int *numAtoms, ATOM **atms )
{
   if ( printStatus ) printf( "\nreading atoms from %s... ", atomsFile );

   double startT = getTime( );

   FILE *fp;
   
   fp = fopen( atomsFile, "rt" );
   
   if ( fp == NULL )
     {
      printError( (char *)"Failed to open PQR file (%s)!", atomsFile );
      return false;
     }
   
   for ( int i = 0; i < numAtomType; i++ )
     numAtoms[ i ] = 0;
   
   char line[ 101 ];

   while ( fgets( line, 100, fp ) != NULL )
     {
      int i = 0, j;
      char tmp[ 100 ];    
    
      i = getAlphaString( line, i, tmp );   // get 'ATOM'/'HETATM', and ignore
      
      if ( strcmp( tmp, "ATOM" ) ) continue;

      i = getInt( line, i, &j );            // get atom number, and ignore

      i = getString( line, i, tmp );        // get atom name
      
      switch ( tmp[ 0 ] )
        {
          case 'C': numAtoms[ C ]++; break;
          case 'H': numAtoms[ H ]++; break;
          case 'N': numAtoms[ N ]++; break;
          case 'O': numAtoms[ O ]++; break;
          case 'P': numAtoms[ P ]++; break;
          case 'S': numAtoms[ S ]++; break;                                                  
        }
     }

   int nAtoms = 0;

   for ( int i = 0; i < numAtomType; i++ )
     if ( numAtoms[ i ] > 0 ) 
        {
          atms[ i ] = ( ATOM * ) malloc( numAtoms[ i ] * sizeof( ATOM ) );
          if ( atms[ i ] == NULL )
            {
              printError( (char *)"Failed to allocate memory for atoms!" );
              fclose( fp );
              return false;
            }
          nAtoms += numAtoms[ i ];
        }  
     else atms[ i ] = NULL;
     
   rewind( fp );
   
   int curIndex[ numAtomType ];

   for ( int i = 0; i < numAtomType; i++ )
     curIndex[ i ] = 0;
   
   while ( fgets( line, 100, fp ) != NULL )
     {
      int i = 0, j, k;
      double v;
      char tmp[ 100 ];    
    
      i = getAlphaString( line, i, tmp );   // get 'ATOM'/'HETATM'
      
      if ( strcmp( tmp, "ATOM" ) ) continue;

      i = getInt( line, i, &j );            // get atom number, and ignore

      i = getString( line, i, tmp );        // get atom name

      switch ( tmp[ 0 ] )
        {
          case 'C': k = C; break;
          case 'H': k = H; break;
          case 'N': k = N; break;
          case 'O': k = O; break;
          case 'P': k = P; break;
          case 'S': k = S; break;
          default : continue;
        }      

      i = getAlphaString( line, i, tmp );   // get residue name, and ignore

      i = getInt( line, i, &j );            // get residue number, and ignore

      j = curIndex[ k ]++;

      i = getDouble( line, i, &v );         // get X coordinate
      atms[ k ][ j ].x = v;    
    
      i = getDouble( line, i, &v );         // get Y coordinate
      atms[ k ][ j ].y = v;    

      i = getDouble( line, i, &v );         // get Z coordinate
      atms[ k ][ j ].z = v;    

      i = getDouble( line, i, &v );         // get charge and ignore

      i = getDouble( line, i, &v );         // get radius and ignore
     }
     
   fclose( fp );

   double endT = getTime( );
   
   if ( printStatus ) printf( "done ( %lf sec, read %d atoms )\n", endT - startT, nAtoms );
         
   return true;
}


bool fastLJ::readAtomsFromPDB( char *atomsFile, int *numAtoms, ATOM **atms )
{
   if ( printStatus ) printf( "\nreading atoms from %s... ", atomsFile );

   double startT = getTime( );

   FILE *fp;
   
   fp = fopen( atomsFile, "rt" );
   
   if ( fp == NULL )
     {
      printError( (char *)"Failed to open PDB file (%s)!", atomsFile );
      return false;
     }
   
   for ( int i = 0; i < numAtomType; i++ )
     numAtoms[ i ] = 0;
   
   char line[ 101 ];

   while ( fgets( line, 100, fp ) != NULL )
     {
      int i = 0, j;
      char tmp[ 100 ];    
    
      i = getAlphaString( line, i, tmp );   // get 'ATOM'/'HETATM'
      
      if ( strcmp( tmp, "ATOM" ) ) continue;

      switch ( line[ 12 ] )
        {
          case 'C': numAtoms[ C ]++; break;
          case 'H': numAtoms[ H ]++; break;
          case 'N': numAtoms[ N ]++; break;
          case 'O': numAtoms[ O ]++; break;
          case 'P': numAtoms[ P ]++; break;
          case 'S': numAtoms[ S ]++; break;                                                  
        }
     }

   int nAtoms = 0;

   for ( int i = 0; i < numAtomType; i++ )
     if ( numAtoms[ i ] > 0 ) 
        {
          atms[ i ] = ( ATOM * ) malloc( numAtoms[ i ] * sizeof( ATOM ) );
          if ( atms[ i ] == NULL )
            {
              printError( (char *)"Failed to allocate memory for atoms!" );
              fclose( fp );
              return false;
            }
          nAtoms += numAtoms[ i ];
        }  
     else atms[ i ] = NULL;
     
   rewind( fp );
   
   int curIndex[ numAtomType ];

   for ( int i = 0; i < numAtomType; i++ )
     curIndex[ i ] = 0;
   
   while ( fgets( line, 100, fp ) != NULL )
     {
      int i = 0, j, k;
      double v;
      char tmp[ 100 ];    
    
      i = getAlphaString( line, i, tmp );   // get 'ATOM'/'HETATM'
      
      if ( strcmp( tmp, "ATOM" ) ) continue;

      switch ( line[ 12 ] )
        {
          case 'C': k = C; break;
          case 'H': k = H; break;
          case 'N': k = N; break;
          case 'O': k = O; break;
          case 'P': k = P; break;
          case 'S': k = S; break;
          default : continue;
        }      
      
      j = curIndex[ k ]++;

      i = 30;

      i = getDouble( line, i, &v );         // get X coordinate
      atms[ k ][ j ].x = v;    
    
      i = getDouble( line, i, &v );         // get Y coordinate
      atms[ k ][ j ].y = v;    

      i = getDouble( line, i, &v );         // get Z coordinate
      atms[ k ][ j ].z = v;    
     }
     
   fclose( fp );

   double endT = getTime( );
   
   if ( printStatus ) printf( "done ( %lf sec, read %d atoms )\n", endT - startT, nAtoms );
         
   return true;
}



bool fastLJ::readAtomsFromFile( char *atomsFile, int *numAtoms, ATOM **atms )
{
   int l = strlen( atomsFile );
   
   if ( ( l > 3 ) && ( atomsFile[ l - 4 ] == '.' ) && ( toupper( atomsFile[ l - 3 ] ) == 'P' ) 
     && ( toupper( atomsFile[ l - 2 ] ) == 'D' ) && ( toupper( atomsFile[ l - 1 ] ) == 'B' ) )
     return readAtomsFromPDB( atomsFile, numAtoms, atms );
   else if ( ( l > 3 ) && ( atomsFile[ l - 4 ] == '.' ) && ( toupper( atomsFile[ l - 3 ] ) == 'P' ) 
     && ( toupper( atomsFile[ l - 2 ] ) == 'Q' ) && ( toupper( atomsFile[ l - 1 ] ) == 'R' ) )
     return readAtomsFromPQR( atomsFile, numAtoms, atms );
   else  
     {
      printError( (char *)"Unknown file type (%s)!", atomsFile );
      return false;
     }
}


bool fastLJ::copyAtomsFromArray( int numAtomsSrc, double *atmsSrcXYZ, char *atmsSrcType, int *numAtomsDest, ATOM **atmsDest )
{
   if ( printStatus ) printf( "\ncopying atoms from array... " );

   double startT = getTime( );
   
   if ( numAtomsSrc <= 0 )
     {
      printError( (char *)"No atoms to copy!" );
      return false;
     }

   for ( int i = 0; i < numAtomType; i++ )
     numAtomsDest[ i ] = 0;
   
   for ( int i = 0; i < numAtomsSrc; i++ )
     {
      switch ( atmsSrcType[ i ] )
        {
          case 'C': numAtomsDest[ C ]++; break;
          case 'H': numAtomsDest[ H ]++; break;
          case 'N': numAtomsDest[ N ]++; break;
          case 'O': numAtomsDest[ O ]++; break;
          case 'P': numAtomsDest[ P ]++; break;
          case 'S': numAtomsDest[ S ]++; break;                                                  
        }     
     }
   
   for ( int i = 0; i < numAtomType; i++ )
     if ( numAtomsDest[ i ] > 0 ) 
        {
          atmsDest[ i ] = ( ATOM * ) malloc( numAtomsDest[ i ] * sizeof( ATOM ) );
          if ( atmsDest[ i ] == NULL )
            {
              printError( (char *)"Failed to allocate memory for atoms!" );
              return false;
            }
        }  
     else atmsDest[ i ] = NULL;

   int curIndex[ numAtomType ];

   for ( int i = 0; i < numAtomType; i++ )
     curIndex[ i ] = 0;

   for ( int i = 0; i < numAtomsSrc; i++ )   
     {          
      int j, k;

      switch ( atmsSrcType[ i ] )
        {
          case 'C': k = C; break;
          case 'H': k = H; break;
          case 'N': k = N; break;
          case 'O': k = O; break;
          case 'P': k = P; break;
          case 'S': k = S; break;
          default : continue;
        }      
      
      j = curIndex[ k ]++;
          
      atmsDest[ k ][ j ].x = atmsSrcXYZ[ 3 * i + 0 ];
      atmsDest[ k ][ j ].y = atmsSrcXYZ[ 3 * i + 1 ];
      atmsDest[ k ][ j ].z = atmsSrcXYZ[ 3 * i + 2 ];
     }    
     
   double endT = getTime( );
   
   if ( printStatus ) printf( "done ( %lf sec, copied %d atoms )\n", endT - startT, numAtomsSrc );
         
   return true;
}


void fastLJ::countAtomsOctreeNodesAndSortAtoms( ATOM **atoms, int *atomsStartID, int *atomsEndID, 
                                                ATOM **atomsT, int *numNodes )
{

  //if (*atomsStartID > *atomsEndID) return;
  //for ( int k = 0; k < numAtomType; k++ )
  //  cout << "Enter " << k << " " << atomsStartID[k] << " " << atomsEndID[k] << endl;
  //cout << "NN " << *numNodes << endl;

   //double minX = atoms[ 0 ][ atomsStartID[ 0 ] ].x, minY = atoms[ 0 ][ atomsStartID[ 0 ] ].y, minZ = atoms[ 0 ][ atomsStartID[ 0 ] ].z;
   //double maxX = atoms[ 0 ][ atomsStartID[ 0 ] ].x, maxY = atoms[ 0 ][ atomsStartID[ 0 ] ].y, maxZ = atoms[ 0 ][ atomsStartID[ 0 ] ].z;

   double minX = 100000000.0;
   double minY = 100000000.0;
   double minZ = 100000000.0;

   double maxX = -100000000.0;
   double maxY = -100000000.0;
   double maxZ = -100000000.0;

   for ( int k = 0; k < numAtomType; k++ )
     for ( int i = atomsStartID[ k ]; i <= atomsEndID[ k ]; i++ )
        { 
          if ( atoms[ k ][ i ].x < minX ) minX = atoms[ k ][ i ].x;      
          if ( atoms[ k ][ i ].x > maxX ) maxX = atoms[ k ][ i ].x;      
      
          if ( atoms[ k ][ i ].y < minY ) minY = atoms[ k ][ i ].y;      
          if ( atoms[ k ][ i ].y > maxY ) maxY = atoms[ k ][ i ].y;      

          if ( atoms[ k ][ i ].z < minZ ) minZ = atoms[ k ][ i ].z;      
          if ( atoms[ k ][ i ].z > maxZ ) maxZ = atoms[ k ][ i ].z;      
        } 

   double cx = ( minX + maxX ) / 2,
          cy = ( minY + maxY ) / 2,
          cz = ( minZ + maxZ ) / 2;
   
   //
   //double r2 = ( atoms[ 0 ][ atomsStartID[ 0 ] ].x - cx ) * ( atoms[ 0 ][ atomsStartID[ 0 ] ].x - cx )
   //          + ( atoms[ 0 ][ atomsStartID[ 0 ] ].y - cy ) * ( atoms[ 0 ][ atomsStartID[ 0 ] ].y - cy )
   //          + ( atoms[ 0 ][ atomsStartID[ 0 ] ].z - cz ) * ( atoms[ 0 ][ atomsStartID[ 0 ] ].z - cz );
   double r2 = 0.0;

   for ( int k = 0; k < numAtomType; k++ )
     for ( int i = atomsStartID[ k ]; i <= atomsEndID[ k ]; i++ )
        { 
          double r2T = ( atoms[ k ][ i ].x - cx ) * ( atoms[ k ][ i ].x - cx )
                     + ( atoms[ k ][ i ].y - cy ) * ( atoms[ k ][ i ].y - cy )
                     + ( atoms[ k ][ i ].z - cz ) * ( atoms[ k ][ i ].z - cz );
         
          if ( r2T > r2 ) r2 = r2T;        
        }

   double cr = sqrt( r2 ) + ( minInterAtomDist * 2.0 );
   
   int numAtoms = 0;

   for ( int k = 0; k < numAtomType; k++ )
     if ( atomsEndID[ k ] >= atomsStartID[ k ] ) 
        numAtoms += ( atomsEndID[ k ] - atomsStartID[ k ] + 1 ); 
      
   *numNodes = 1;
      
   //cout << "arand " << numAtoms << " " << cr << " " << maxLeafSize << " " << minRadius << endl;

   if ( ( numAtoms > maxLeafSize ) && ( cr > minRadius ) )
     {
      int atomsCount[ 8 ][ numAtomType ], atomsCountAllTypes[ 8 ] = { 0, 0, 0, 0, 0, 0, 0, 0 };
      
      for ( int i = 0; i < 8; i++ )
         for ( int j = 0; j < numAtomType; j++ )
            atomsCount[ i ][ j ] = 0;

      for ( int k = 0; k < numAtomType; k++ )
        for ( int i = atomsStartID[ k ]; i <= atomsEndID[ k ]; i++ )
          {
           atomsT[ k ][ i ] = atoms[ k ][ i ];
           
           int j = ( zeroIfLess( atoms[ k ][ i ].z, cz ) << 2 )
                 + ( zeroIfLess( atoms[ k ][ i ].y, cy ) << 1 )
                 + ( zeroIfLess( atoms[ k ][ i ].x, cx ) );
           
           atomsCount[ j ][ k ]++;
           atomsCountAllTypes[ j ]++;
          }       
     // for (int j=0; j<8; j++) {
//	cout << atomsCount[j][1] << " ";
      // cout << "done" << endl;

      int atomsStartIndex[ 8 ][ numAtomType ], atomsEndIndex[ 8 ][ numAtomType ];
      int atomsCurIndex[ 8 ][ numAtomType ];              
      
      for ( int k = 0; k < numAtomType; k++ )
        {
         atomsCurIndex[ 0 ][ k ] = atomsStartIndex[ 0 ][ k ] = atomsStartID[ k ];
         for ( int i = 1; i < 8; i++ )
            atomsCurIndex[ i ][ k ] = atomsStartIndex[ i ][ k ] = atomsStartIndex[ i - 1 ][ k ] + atomsCount[ i - 1 ][ k ];
            
         for ( int i = atomsStartID[ k ]; i <= atomsEndID[ k ]; i++ )
           {        
            int j = ( zeroIfLess( atomsT[ k ][ i ].z, cz ) << 2 )
                  + ( zeroIfLess( atomsT[ k ][ i ].y, cy ) << 1 )
                  + ( zeroIfLess( atomsT[ k ][ i ].x, cx ) );
             
            atoms[ k ][ atomsCurIndex[ j ][ k ] ] = atomsT[ k ][ i ];
            atomsCurIndex[ j ][ k ]++;  
           }                    
        }    
        
      for ( int i = 0; i < 8; i++ ) 
        if ( atomsCountAllTypes[ i ] > 0 )
          {
           int numNodesT = 0;
           
           for ( int k = 0; k < numAtomType; k++ )
	     {
	       atomsEndIndex[ i ][ k ] = atomsStartIndex[ i ][ k ] + atomsCount[ i ][ k ] - 1;
	     }	   

           countAtomsOctreeNodesAndSortAtoms( atoms, atomsStartIndex[ i ], atomsEndIndex[ i ], atomsT, &numNodesT );

           *numNodes += numNodesT;
          }
     }
}



int fastLJ::constructAtomsOctree( int *atomsStartID, int *atomsEndID, ATOM **atoms, ATOMS_OCTREE_NODE *atomsOctree )
{
   int nodeID = nextFreeNode( );
   
   if ( nodeID < 0 ) return nodeID;
      
   for ( int k = 0; k < numAtomType; k++ ) 
     {  
      atomsOctree[ nodeID ].atomsStartID[ k ] = atomsStartID[ k ];
      atomsOctree[ nodeID ].atomsEndID[ k ] = atomsEndID[ k ];
     } 

   // arand, bug fix, 08-30-2011
   //double minX = atoms[ 0 ][ atomsStartID[ 0 ] ].x, minY = atoms[ 0 ][ atomsStartID[ 0 ] ].y, minZ = atoms[ 0 ][ atomsStartID[ 0 ] ].z;
   //double maxX = atoms[ 0 ][ atomsStartID[ 0 ] ].x, maxY = atoms[ 0 ][ atomsStartID[ 0 ] ].y, maxZ = atoms[ 0 ][ atomsStartID[ 0 ] ].z;

   double minX = 1000000000.0;
   double minY = 1000000000.0;
   double minZ = 1000000000.0;
   double maxX = -1000000000.0;
   double maxY = -1000000000.0;
   double maxZ = -1000000000.0;

   for ( int k = 0; k < numAtomType; k++ )
     for ( int i = atomsStartID[ k ]; i <= atomsEndID[ k ]; i++ )
        { 
          if ( atoms[ k ][ i ].x < minX ) minX = atoms[ k ][ i ].x;      
          if ( atoms[ k ][ i ].x > maxX ) maxX = atoms[ k ][ i ].x;      
      
          if ( atoms[ k ][ i ].y < minY ) minY = atoms[ k ][ i ].y;      
          if ( atoms[ k ][ i ].y > maxY ) maxY = atoms[ k ][ i ].y;      

          if ( atoms[ k ][ i ].z < minZ ) minZ = atoms[ k ][ i ].z;      
          if ( atoms[ k ][ i ].z > maxZ ) maxZ = atoms[ k ][ i ].z;      
        } 

   double cx = atomsOctree[ nodeID ].cx = ( minX + maxX ) / 2;
   double cy = atomsOctree[ nodeID ].cy = ( minY + maxY ) / 2;
   double cz = atomsOctree[ nodeID ].cz = ( minZ + maxZ ) / 2;

   //double r2 = ( atoms[ 0 ][ atomsStartID[ 0 ] ].x - cx ) * ( atoms[ 0 ][ atomsStartID[ 0 ] ].x - cx )
   //          + ( atoms[ 0 ][ atomsStartID[ 0 ] ].y - cy ) * ( atoms[ 0 ][ atomsStartID[ 0 ] ].y - cy )
   //          + ( atoms[ 0 ][ atomsStartID[ 0 ] ].z - cz ) * ( atoms[ 0 ][ atomsStartID[ 0 ] ].z - cz );

   // arand, 08-30-2011
   double r2 = 0.0;

   for ( int k = 0; k < numAtomType; k++ )
     for ( int i = atomsStartID[ k ]; i <= atomsEndID[ k ]; i++ )
        { 
          double r2T = ( atoms[ k ][ i ].x - cx ) * ( atoms[ k ][ i ].x - cx )
                     + ( atoms[ k ][ i ].y - cy ) * ( atoms[ k ][ i ].y - cy )
                     + ( atoms[ k ][ i ].z - cz ) * ( atoms[ k ][ i ].z - cz );
         
          if ( r2T > r2 ) r2 = r2T;        
        }
   
   double cr = atomsOctree[ nodeID ].cr = sqrt( r2 ) + ( minInterAtomDist * 2.0 );

   int numAtoms = 0;

   for ( int k = 0; k < numAtomType; k++ )
     {
      if ( atomsEndID[ k ] >= atomsStartID[ k ] ) 
         atomsOctree[ nodeID ].numAtoms[ k ] = atomsEndID[ k ] - atomsStartID[ k ] + 1;
      else 
         atomsOctree[ nodeID ].numAtoms[ k ] = 0;   
      numAtoms += atomsOctree[ nodeID ].numAtoms[ k ]; 
     }

   if ( ( numAtoms <= maxLeafSize ) || ( cr <= minRadius ) )
      atomsOctree[ nodeID ].leaf = true;
   else
     {
      atomsOctree[ nodeID ].leaf = false;     

      int atomsCount[ 8 ][ numAtomType ], atomsCountAllTypes[ 8 ] = { 0, 0, 0, 0, 0, 0, 0, 0 };
      
      for ( int i = 0; i < 8; i++ )
         for ( int j = 0; j < numAtomType; j++ )
            atomsCount[ i ][ j ] = 0;

      for ( int k = 0; k < numAtomType; k++ )
        for ( int i = atomsStartID[ k ]; i <= atomsEndID[ k ]; i++ )
          {
           int j = ( zeroIfLess( atoms[ k ][ i ].z, cz ) << 2 )
                 + ( zeroIfLess( atoms[ k ][ i ].y, cy ) << 1 )
                 + ( zeroIfLess( atoms[ k ][ i ].x, cx ) );
           
           atomsCount[ j ][ k ]++;
           atomsCountAllTypes[ j ]++;
          }
        
      int atomsStartIndex[ 8 ][ numAtomType ], atomsEndIndex[ 8 ][ numAtomType ];

      for ( int k = 0; k < numAtomType; k++ )
        {
         atomsStartIndex[ 0 ][ k ] = atomsStartID[ k ];
         for ( int i = 1; i < 8; i++ )
            atomsStartIndex[ i ][ k ] = atomsStartIndex[ i - 1 ][ k ] + atomsCount[ i - 1 ][ k ];
        }    
      
      for ( int i = 0; i < 8; i++ ) 
        if ( atomsCountAllTypes[ i ] > 0 )
          {
           for ( int k = 0; k < numAtomType; k++ )
              atomsEndIndex[ i ][ k ] = atomsStartIndex[ i ][ k ] + atomsCount[ i ][ k ] - 1;

           int j = constructAtomsOctree( atomsStartIndex[ i ], atomsEndIndex[ i ], atoms, atomsOctree );
           atomsOctree[ nodeID ].cPtr[ i ] = j; 
          }
        else atomsOctree[ nodeID ].cPtr[ i ] = -1;  
     }  
     
   return nodeID;  
}


bool fastLJ::buildStaticAtomsOctree( void )
{  
   if ( printStatus ) printf( "\nbuilding static atoms octree... " );
   
   double startT = getTime( );
   
   ATOM *atomsT[ numAtomType ];
   
   for ( int k = 0; k < numAtomType; k++ )
     if ( nStaticAtoms[ k ] > 0 )
        {
         atomsT[ k ] = ( ATOM * ) malloc( nStaticAtoms[ k ] * sizeof( ATOM ) );
   
         if ( atomsT[ k ] == NULL )
           {
            printError( (char *)"Failed to allocate temporary memory for static atoms!" );
            if ( !staticAtomsOctreeBuilt ) exit( 1 );
            return false;
           }
        }
     else atomsT[ k ] = NULL;

   int atomsStartIndex[ numAtomType ], atomsEndIndex[ numAtomType ];
   
   for ( int k = 0; k < numAtomType; k++ )
     {
      atomsStartIndex[ k ] = 0;
      atomsEndIndex[ k ] = nStaticAtoms[ k ] - 1;
     }   
           
   countAtomsOctreeNodesAndSortAtoms( staticAtoms, atomsStartIndex, atomsEndIndex, atomsT, &numStaticAtomsOctreeNodes );
   
   for ( int k = 0; k < numAtomType; k++ )
     freeMem( atomsT[ k ] );
   
   ATOMS_OCTREE_NODE *atomsOctreeT;
   
   atomsOctreeT = ( ATOMS_OCTREE_NODE * ) malloc( numStaticAtomsOctreeNodes * sizeof( ATOMS_OCTREE_NODE ) );

   if ( atomsOctreeT == NULL )
     {
      printError( (char *)"Unable to %s static atoms octree - memory allocation failed!", ( staticAtomsOctreeBuilt ) ? "rebuild" : "build" );
      if ( !staticAtomsOctreeBuilt ) exit( 1 );
      return false;
     }
 
   freeMem( staticAtomsOctree );
   staticAtomsOctree = atomsOctreeT; 

   initFreeNodeServer( numStaticAtomsOctreeNodes );
   
   staticAtomsOctreeRoot = constructAtomsOctree( atomsStartIndex, atomsEndIndex, staticAtoms, staticAtomsOctree );

   staticAtomsOctreeBuilt = true;

   double endT = getTime( );
   
   if ( printStatus ) printf( "done ( %lf sec )\n", endT - startT );
   
   return true;
}



bool fastLJ::buildMovingAtomsOctree( void )
{  
  //cout << "Also here" << endl;
  printStatus = true;
   if ( printStatus ) printf( "\nbuilding moving atoms octree... " );
   //cout << endl;

   double startT = getTime( );

   ATOM *atomsT[ numAtomType ];

   //cout << "A " << numAtomType << endl;
   
   for ( int k = 0; k < numAtomType; k++ ) {
     //cout << k << " " << nMovingAtoms[k] << endl;
     if ( nMovingAtoms[ k ] > 0 )
        {
         atomsT[ k ] = ( ATOM * ) malloc( nMovingAtoms[ k ] * sizeof( ATOM ) );
   
         if ( atomsT[ k ] == NULL )
           {
            printError( (char *)"Failed to allocate temporary memory for moving atoms!" );
            if ( !movingAtomsOctreeBuilt ) exit( 1 );
            return false;
           }
        }
     else atomsT[ k ] = NULL;
   }

   int atomsStartIndex[ numAtomType ], atomsEndIndex[ numAtomType ];
   
   //cout << "AA" << endl;

   for ( int k = 0; k < numAtomType; k++ )
     {
      // cout << k << " " << nMovingAtoms[k] << endl;
      atomsStartIndex[ k ] = 0;
      atomsEndIndex[ k ] = nMovingAtoms[ k ] - 1;
     }   
   //cout << "AAA" << endl;
   
   countAtomsOctreeNodesAndSortAtoms( movingAtomsOrig, atomsStartIndex, atomsEndIndex, atomsT, &numMovingAtomsOctreeNodes );

   //cout << "B" << endl;
   
   for ( int k = 0; k < numAtomType; k++ )
     freeMem( atomsT[ k ] );
   
   ATOMS_OCTREE_NODE *atomsOctreeT;
   
   atomsOctreeT = ( ATOMS_OCTREE_NODE * ) malloc( numMovingAtomsOctreeNodes * sizeof( ATOMS_OCTREE_NODE ) );

   if ( atomsOctreeT == NULL )
     {
      printError( (char *)"Unable to %s moving atoms octree - memory allocation failed!", ( movingAtomsOctreeBuilt ) ? "rebuild" : "build" );
      if ( !movingAtomsOctreeBuilt ) exit( 1 );
      return false;
     }

  // cout << "B1" << endl;

   freeMem( movingAtomsOctree );
   movingAtomsOctree = atomsOctreeT;  

   //cout << "C" << endl;

   initFreeNodeServer( numMovingAtomsOctreeNodes );
   
   movingAtomsOctreeRoot = constructAtomsOctree( atomsStartIndex, atomsEndIndex, movingAtomsOrig, movingAtomsOctree );
   
   //cout << "D" << endl;

   movingAtomsOctreeBuilt = true;
   
   double endT = getTime( );
   
   if ( printStatus ) printf( "done ( %lf sec )\n", endT - startT );
   
   return true;
}


bool fastLJ::buildOctrees( void )
{
   if ( ( minRadius != minRadiusUsed ) || ( maxLeafSize != maxLeafSizeUsed ) 
     || ( minInterAtomDist != minInterAtomDistUsed ) || ( !staticAtomsOctreeBuilt ) ) 
        buildStaticAtomsOctree( );  
        
   //cout << "Got here" << endl;

   if ( ( minRadius != minRadiusUsed ) || ( maxLeafSize != maxLeafSizeUsed ) 
     || ( minInterAtomDist != minInterAtomDistUsed ) || ( !movingAtomsOctreeBuilt ) ) 
        buildMovingAtomsOctree( );
   
   minRadiusUsed = minRadius;
   maxLeafSizeUsed = maxLeafSize;
   minInterAtomDistUsed = minInterAtomDist; 
}


#ifdef USE_SSE

inline float fastLJ::vectorLJ( v4sf xi, v4sf yi, v4sf zi,
                               v4sf xj, v4sf yj, v4sf zj, 
                               v4sf Aij, v4sf Bij )
{
   v4sf d2, d6, v;
   V4SF val;
 
   xi = _mm_sub_ps( xi, xj );                       
   d2 = _mm_mul_ps( xi, xi );
   
   yi = _mm_sub_ps( yi, yj );
   yi = _mm_mul_ps( yi, yi );
   d2 = _mm_add_ps( d2, yi );
   
   zi = _mm_sub_ps( zi, zj );
   zi = _mm_mul_ps( zi, zi );                       
   d2 = _mm_add_ps( d2, zi );

   d2 = _mm_max_ps( d2, MIND2 );
   
   d6 = _mm_mul_ps( d2, d2 );
   d6 = _mm_mul_ps( d6, d2 );
   
   v = _mm_div_ps( Bij, d6 );
                          
   d6 = _mm_mul_ps( d6, d6 );
   d6 = _mm_div_ps( Aij, d6 );
   
   val.v = _mm_sub_ps( d6, v );
  
   return ( val.f[ 0 ] + val.f[ 1 ] + val.f[ 2 ] + val.f[ 3 ] );
}

#endif


void fastLJ::approximatePotential( int threadID, double *transMat, int nodeS, int nodeM, double *LJPot )
{
   double sumRad = staticAtomsOctree[ nodeS ].cr + movingAtomsOctree[ nodeM ].cr;
   double sumRad2 = sumRad * sumRad;   
   double cxM, cyM, czM;
   
   transform( movingAtomsOctree[ nodeM ].cx, movingAtomsOctree[ nodeM ].cy, movingAtomsOctree[ nodeM ].cz, transMat, cxM, cyM, czM );
    
   double dx = staticAtomsOctree[ nodeS ].cx - cxM,
          dy = staticAtomsOctree[ nodeS ].cy - cyM,
          dz = staticAtomsOctree[ nodeS ].cz - czM;
   double d2 = dx * dx + dy * dy + dz * dz;
   
   bool farEnough = false;
   
   if ( /*( d2 > 4 * sumRad2 ) && */( d2 > ( sumRad2 / ( epsilon * epsilon ) ) ) ) farEnough = true;
   
   register double pot = 0;
   
   if ( farEnough ) 
      {
       double d6 = 1.0 / ( d2 * d2 * d2 );
       double d12 = d6 * d6;
       
       for ( int k = 0; k < numAtomType; k++ )
         for ( int l = 0; l < numAtomType; l++ )  
             pot += movingAtomsOctree[ nodeM ].numAtoms[ k ] * staticAtomsOctree[ nodeS ].numAtoms[ l ] * ( ( A[ k ][ l ] * d12 ) - ( B[ k ][ l ] * d6 ) );       
      }  
   else
      {
       if ( staticAtomsOctree[ nodeS ].leaf && movingAtomsOctree[ nodeM ].leaf )
         {
#ifdef USE_SSE     
           V4SF xi, yi, zi;
           V4SF xj, yj, zj;         
           V4SF Aij, Bij;         
           
           int t = 0;
#endif         
           for ( int k = 0; k < numAtomType; k++ )
             for ( int l = 0; l < numAtomType; l++ )
               {
                double Akl = A[ k ][ l ], Bkl = B[ k ][ l ];
                
                int s = threadID * nMovingAtoms[ k ];
                
                for ( int i = movingAtomsOctree[ nodeM ].atomsStartID[ k ]; i <= movingAtomsOctree[ nodeM ].atomsEndID[ k ]; i++ )
                  {             
                   double mx = movingAtoms[ k ][ s + i ].x, my = movingAtoms[ k ][ s + i ].y, mz = movingAtoms[ k ][ s + i ].z;
#ifdef USE_SSE
                   if ( useSSEFunctions )
                     {               
                      for ( int j = staticAtomsOctree[ nodeS ].atomsStartID[ l ]; j <= staticAtomsOctree[ nodeS ].atomsEndID[ l ]; j++ )                
                        {    
                          xi.f[ t ] = mx;
                          yi.f[ t ] = my;
                          zi.f[ t ] = mz;                   
                        
                          xj.f[ t ] = staticAtoms[ l ][ j ].x;
                          yj.f[ t ] = staticAtoms[ l ][ j ].y;
                          zj.f[ t ] = staticAtoms[ l ][ j ].z;
                          
                          Aij.f[ t ] = Akl;
                          Bij.f[ t ] = Bkl;                   
                                            
                          if ( ++t == 4 ) 
                            {
                              pot += vectorLJ( xi.v, yi.v, zi.v, xj.v, yj.v, zj.v, Aij.v, Bij.v );
                              t = 0;
                            }  
                        }                    
                     }
                   else
                     {
                      for ( int j = staticAtomsOctree[ nodeS ].atomsStartID[ l ]; j <= staticAtomsOctree[ nodeS ].atomsEndID[ l ]; j++ )                
                        {                                            
                          dx = staticAtoms[ l ][ j ].x - mx;
                          dy = staticAtoms[ l ][ j ].y - my;
                          dz = staticAtoms[ l ][ j ].z - mz;                     
                                                    
                          d2 = dx * dx + dy * dy + dz * dz;       
                     
                          if ( d2 < minD2 ) d2 = minD2;
                          
                          double d6 = 1.0 / ( d2 * d2 * d2 );
                          double d12 = d6 * d6;
                         
                          pot += ( ( Akl * d12 ) - ( Bkl * d6 ) );       
                        }                                 
                     }  
#else
                  for ( int j = staticAtomsOctree[ nodeS ].atomsStartID[ l ]; j <= staticAtomsOctree[ nodeS ].atomsEndID[ l ]; j++ )
                    {                                            
                      dx = staticAtoms[ l ][ j ].x - mx;
                      dy = staticAtoms[ l ][ j ].y - my;
                      dz = staticAtoms[ l ][ j ].z - mz;                     
                                                    
                      d2 = dx * dx + dy * dy + dz * dz;       
                     
                      if ( d2 < minD2 ) d2 = minD2;
    
                      double d6 = 1.0 / ( d2 * d2 * d2 );
                      double d12 = d6 * d6;
                         
                      pot += ( ( Akl * d12 ) - ( Bkl * d6 ) );      
                    }                               
#endif
                 }   
               }  
             
#ifdef USE_SSE
           while ( t-- )
             {
               double d2 = ( xi.f[ t ] - xj.f[ t ] ) * ( xi.f[ t ] - xj.f[ t ] )
                         + ( yi.f[ t ] - yj.f[ t ] ) * ( yi.f[ t ] - yj.f[ t ] )
                         + ( zi.f[ t ] - zj.f[ t ] ) * ( zi.f[ t ] - zj.f[ t ] );                              
                  
               if ( d2 < minD2 ) d2 = minD2;
               
               double d6 = 1.0 / ( d2 * d2 * d2 );
               double d12 = d6 * d6;
                         
               pot += ( ( Aij.f[ t ] * d12 ) - ( Bij.f[ t ] * d6 ) );               
            }                    
#endif             
         }
       else if ( !staticAtomsOctree[ nodeS ].leaf && !movingAtomsOctree[ nodeM ].leaf )         
              {
                for ( int i = 0; i < 8; i++ )
                  if ( staticAtomsOctree[ nodeS ].cPtr[ i ] >= 0 )
                     for ( int j = 0; j < 8; j++ )
                       if ( movingAtomsOctree[ nodeM ].cPtr[ j ] >= 0 ) 
                         {
                           double ljPot;
                           
                           approximatePotential( threadID, transMat, staticAtomsOctree[ nodeS ].cPtr[ i ], movingAtomsOctree[ nodeM ].cPtr[ j ], &ljPot );
                           
                           pot += ljPot;                         
                         }  
              }
            else if ( !staticAtomsOctree[ nodeS ].leaf )         
                   {
                     for ( int i = 0; i < 8; i++ )
                       if ( staticAtomsOctree[ nodeS ].cPtr[ i ] >= 0 )
                         {
                           double ljPot;
                           
                           approximatePotential( threadID, transMat, staticAtomsOctree[ nodeS ].cPtr[ i ], nodeM, &ljPot );
                           
                           pot += ljPot;
                         }  
                   }
                 else 
                   {
                     for ( int j = 0; j < 8; j++ )
                       if ( movingAtomsOctree[ nodeM ].cPtr[ j ] >= 0 )
                         {
                           double ljPot;
                           
                           approximatePotential( threadID, transMat, nodeS, movingAtomsOctree[ nodeM ].cPtr[ j ], &ljPot );
                           
                           pot += ljPot;
                         }  
                   }              
      }      
      
   ( *LJPot ) = pot;            
}



bool fastLJ::computePotential( int threadID, double *transMat, double *LJPot )
{
   if ( !staticAtomsOctreeBuilt || !movingAtomsOctreeBuilt ) return false;

   if ( printStatus ) printf( "\napproximating Lennard-Jones potential... " );

   double startT = getTime( );

   transformMovingAtoms( threadID, transMat );
      
   approximatePotential( threadID, transMat, staticAtomsOctreeRoot, movingAtomsOctreeRoot, LJPot );   
         
   double endT = getTime( );
   
   if ( printStatus ) printf( "done ( %lf sec )\n", endT - startT );
   
   return true;
}



bool fastLJ::computePotential( int threadID, double *LJPot )
{
   return computePotential( threadID, transMatrix, LJPot );
}


void fastLJ::computePotentialNaively( int threadID, double *transMat, double *LJPot )
{
   if ( printStatus ) printf( "\ncomputing Lennard-Jones potential naively... " );

   double startT = getTime( );
 
   transformMovingAtoms( threadID, transMat );
      
   double minD2 = minInterAtomDist * minInterAtomDist;
   
   register double pot = 0;   
   
   for ( int k = 0; k < numAtomType; k++ )
     for ( int l = 0; l < numAtomType; l++ )  
        {
         double Aij = A[ k ][ l ], Bij = B[ k ][ l ];
         
         int s = threadID * nMovingAtoms[ k ];
         
         for ( int i = 0; i < nMovingAtoms[ k ]; i++ )
           {
            double mx = movingAtoms[ k ][ s + i ].x, my = movingAtoms[ k ][ s + i ].y, mz = movingAtoms[ k ][ s + i ].z;
            
            for ( int j = 0; j < nStaticAtoms[ l ]; j++ )
              {                            
               double dx = staticAtoms[ l ][ j ].x - mx;
               double dy = staticAtoms[ l ][ j ].y - my;
               double dz = staticAtoms[ l ][ j ].z - mz;                     
                      
               double d2 = dx * dx + dy * dy + dz * dz;
               
               if ( d2 < minD2 ) d2 = minD2; 
               
               double d6 = 1.0 / ( d2 * d2 * d2 );
               double d12 = d6 * d6;
               
               pot += ( ( Aij * d12 ) - ( Bij * d6 ) );       
              }                   
           }    
        }      

   ( *LJPot ) = pot;
     
   double endT = getTime( );
   
   if ( printStatus ) printf( "done ( %lf sec )\n", endT - startT );
}


void fastLJ::computePotentialNaively( int threadID, double *LJPot )
{
   computePotentialNaively( threadID, transMatrix, LJPot );
}
