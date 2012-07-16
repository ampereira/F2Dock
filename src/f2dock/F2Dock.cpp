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
#include <stdlib.h>
#include <fstream>
#include "TopValues.h"
#include "Docking.h"
#include <vector>
#include <string>
#include "ElementInformation.h"

using namespace std;

/*
int skipWhiteSpaces( char *buf, int i )
{
   int j = i;

   while ( buf[ j ] && isspace( buf[ j ] ) ) j++;

   return j;
}
*/

int countLines(char *filename)
{
  FILE* fp = fopen( filename, "r" );
  int i=0;
  char s[ 2000 ];

  while ( fgets( s, 1999, fp ) != NULL )
    {
      int j = skipWhiteSpaces( s, 0 );
      if ( s[ j ] != '#' ) i++;
    }

  fclose(fp);
  return i;
}


void generateRandomRotationMatrix( Matrix &randRot )
{
//  double alpha = ( ( ( double ) rand( ) ) / ( ( double ) ( RAND_MAX ) ) - 0.5 ) * 2 * M_PI;
//  double beta = ( ( ( double ) rand( ) ) / ( ( double ) ( RAND_MAX ) ) - 0.5 ) * M_PI;
//  double gamma = ( ( ( double ) rand( ) ) / ( ( double ) ( RAND_MAX ) ) - 0.5 ) * 2 * M_PI;;
//
//  randRot.reset( );
//  randRot.preMultiplication( Matrix::rotationZ( gamma ) );
//  randRot.preMultiplication( Matrix::rotationX( beta ) );
//  randRot.preMultiplication( Matrix::rotationZ( alpha ) );


// from: Arvo, James (1992), "Fast random rotation matrices", in David Kirk, Graphics Gems III,
//       San Diego: Academic Press Professional, pp. 117–120, ISBN 978-0-12-409671-4

    double theta = ( ( ( double ) rand( ) ) / ( ( double ) ( RAND_MAX ) ) ) * 2.0 * M_PI;
    double phi   = ( ( ( double ) rand( ) ) / ( ( double ) ( RAND_MAX ) ) ) * 2.0 * M_PI;
    double z     = ( ( ( double ) rand( ) ) / ( ( double ) ( RAND_MAX ) ) ) * 2.0;

    double r  = sqrt( z );
    double Vx = sin( phi ) * r, Vy = cos( phi ) * r, Vz = sqrt( 2.0 - z );
    double st = sin( theta ), ct = cos( theta );
    double Sx = Vx * ct - Vy * st, Sy = Vx * st + Vy * ct;

    randRot.set(  Vx * Sx - ct, Vx * Sy - st, Vx * Vz, 0.0,
                  Vy * Sx + st, Vy * Sy - ct, Vy * Vz, 0.0,
                  Vz * Sx,      Vz * Sy,      1.0 - z, 0.0,
                  0.0,          0.0,          0.0,     1.0 );
}


bool readRotations(char *rotationFile, PARAMS_IN *p)
{
  int numberOfRotations = countLines(rotationFile);

  float *rotations = new float[ ( numberOfRotations + 1 ) * 9 ];

  FILE* fpRot = fopen( rotationFile, "r" );
  if (  fpRot == NULL )
    {
      printf( "Error: Failed to open parameter file %s!\n", rotationFile);
      return false;
    }

  int i=0, off;
  while ( true ) {
    off = i*9;

    if (fscanf( fpRot, "%f %f %f %f %f %f %f %f %f\n",
		&rotations[off  ], &rotations[off+1], &rotations[off+2], \
		&rotations[off+3], &rotations[off+4], &rotations[off+5], \
		&rotations[off+6], &rotations[off+7], &rotations[off+8] ) == 9 )
         i++;
    else
      break;
  }
  fclose( fpRot );
  p->rotations = rotations;
  p->numberOfRotations = i;

  return true;
}


void createRotationsNeighborhoodGraph( PARAMS_IN *p )
{
  int nr = p->numberOfRotations;

  ROTATION_GRAPH_NODE *rotGraph = new ROTATION_GRAPH_NODE[ nr ];
  double cosTheta = cos( ( p->pruneAngle / 180.0 ) * M_PI );
  float *rot = p->rotations;

  printf( "building neighborhood graph for rotations" );

  int m = 0;
  for ( int i = 0; i < nr; i++ )
    {
      if ( !( i % 100 ) )
        {
          printf( "." );
          fflush( stdout );
        }

      rotGraph[ i ].badNode = false;
      rotGraph[ i ].startNeighbors = m;

      for ( int j = i + 1; j < nr; j++ )
        {
           double didj = 0, di2 = 0, dj2 = 0;

           for ( int k = 0; k < 3; k++ )
             {
              double di = 0, dj = 0;

              for ( int l = 0; l < 3; l++ )
                {
                 di += rot[ i * 9 + k * 3 + l ];
                 dj += rot[ j * 9 + k * 3 + l ];
                }

              di2 += di * di;
              dj2 += dj * dj;
              didj += di * dj;
             }

           double dotP = didj / sqrt( di2 * dj2 );

           if ( dotP >= cosTheta )
             {
               p->rotNeighbors.push_back( j );
               m++;
             }
        }

      rotGraph[ i ].endNeighbors = m - 1;
    }

  p->rotGraph = rotGraph;

  printf( "\n" );
}

/*
int getAlphaString( char *buf, int i, char *s )
{
   int j = skipWhiteSpaces( buf, i );

   int k = 0;

   while ( buf[ j ] && isalpha( buf[ j ] ) ) s[ k++ ] = buf[ j++ ];

   s[ k ] = 0;

   return j;
}


int getString( char *buf, int i, char *s )
{
   int j = skipWhiteSpaces( buf, i );

   int k = 0;

   while ( buf[ j ] && ( isalnum( buf[ j ] ) || ispunct( buf[ j ] ) ) ) s[ k++ ] = buf[ j++ ];

   s[ k ] = 0;

   return j;
}


int getDouble( char *buf, int i, double *v )
{
   char s[ 100 ];
   int j = skipWhiteSpaces( buf, i );

   int k = 0;

   if ( buf[ j ] == '-' ) s[ k++ ] = buf[ j++ ];

   while ( buf[ j ] && ( isdigit( buf[ j ] ) || ( buf[ j ] == '.' ) ) ) s[ k++ ] = buf[ j++ ];

   s[ k ] = 0;

   *v = atof( s );

   return j;
}


int getInt( char *buf, int i, int *v )
{
   char s[ 100 ];
   int j = skipWhiteSpaces( buf, i );

   int k = 0;

   if ( buf[ j ] == '-' ) s[ k++ ] = buf[ j++ ];

   while ( buf[ j ] && ( isdigit( buf[ j ] ) ) ) s[ k++ ] = buf[ j++ ];

   s[ k ] = 0;

   *v = atoi( s );

   return j;
}
*/


bool readF2D( char *fileName, double **xkOrig, double **ykOrig,
	      double **zkOrig, float **charges, float **radii,
	      int *numCenters, char **atype, char ***atnames,
	      char ***restype, int **resnum, char **moleculePdb, char **moleculeF2d)
{
  FILE* fp;
  int i=0;
  char buf[2000];
  char *c;
  char sep[] = " ";

  if( !fileName ) return false;

  (*numCenters) = countLines(fileName);
  printf("%d\n",(*numCenters));

  (*atype) = new char [(*numCenters)];
  (*atnames) = new char *[(*numCenters)];
  (*restype) = new char *[(*numCenters)];
  (*resnum) = new int [(*numCenters)];
  (*xkOrig) = new double[(*numCenters)];
  (*ykOrig) = new double[(*numCenters)];
  (*zkOrig) = new double[(*numCenters)];
  (*charges) = new float[(*numCenters)];
  (*radii) = new float[(*numCenters)];

  string pdbExt = "pdb";
  string f2dExt = "f2d";
  string moleculePathPdb = (string(fileName).erase((strlen(fileName)-3), 3)) + pdbExt;
  string moleculePathF2d = (string(fileName).erase((strlen(fileName)-3), 3)) + f2dExt;
  char *tmpPathPdb, *tmpPathF2d;
  tmpPathPdb = new char[moleculePathPdb.size()+1];
  tmpPathF2d = new char[moleculePathF2d.size()+1];
  strcpy(tmpPathPdb, moleculePathPdb.c_str());
  strcpy(tmpPathF2d, moleculePathF2d.c_str());
  (*moleculePdb) = tmpPathPdb;
  (*moleculeF2d) = tmpPathF2d;

  if( (fp = fopen(fileName, "r")) == NULL ) return false;

  for (i=0; i<(*numCenters); i++) {
    fgets(buf, 1999, fp);

//    printf( "%s", buf );

    int j = 0, k;
    double v;
    char tmp[ 100 ];

    j = getAlphaString( buf, j, tmp );   // get 'ATOM'/'HETATM', and ignore

    j = getInt( buf, j, &k );            // get atom number, and ignore

    j = getString( buf, j, tmp );        // get atom name
    (*atnames)[i] = strdup(tmp);

//    j = getAlphaString( buf, j, tmp );   // get residue name
    j = getString( buf, j, tmp );   // get residue name
    (*restype)[i] = strdup(tmp);

    j = getInt( buf, j, &k );            // get residue number
    (*resnum)[i] = k;

    j = getDouble( buf, j, &v );         // get X coordinate
    (*xkOrig)[i] = v;

    j = getDouble( buf, j, &v );         // get Y coordinate
    (*ykOrig)[i] = v;

    j = getDouble( buf, j, &v );         // get Z coordinate
    (*zkOrig)[i] = v;

    j = getDouble( buf, j, &v );         // get charge
    (*charges)[i] = ( float ) v;

    j = getDouble( buf, j, &v );         // get radius
    (*radii)[i] = ( float ) v;

    j = getString( buf, j, tmp );        // get atom type ('I' or 'E')
    (*atype)[i] = tmp[ 0 ];
  }
  return true;

}


void identifyHbondAcceptorsAndDonors( int numCenters, double *xk, double *yk, double *zk, float *rk,
	                              char **atnames, int *resnum, char **hbondType, char *atype, bool staticMol )
{
  ( *hbondType ) = new char [ numCenters ];
  bool *CNeg = new bool [ numCenters ];

  for ( int i = 0; i < numCenters; i++ )
   {
     CNeg[ i ] = false;

     int ii = 0;
     if ( isdigit( atnames[ i ][ ii ] ) ) ii++;

     if ( atnames[ i ][ ii ] == 'C' )
       {
         if ( !atnames[ i ][ ii + 1 ] || ( atnames[ i ][ ii + 1 ] == 'A' ) )
           {
             CNeg[ i ] = true;
             continue;
           }

         for ( int j = 0; j < numCenters; j++ )
           if ( resnum[ j ] == resnum[ i ] )
             {
               int jj = 0;
               if ( isdigit( atnames[ j ][ jj ] ) ) jj++;

               if ( ( atnames[ j ][ jj ] == 'O' ) || ( atnames[ j ][ jj ] == 'N' ) )
                 {
                   double d2 = ( xk[ i ] - xk[ j ] ) * ( xk[ i ] - xk[ j ] )
                             + ( yk[ i ] - yk[ j ] ) * ( yk[ i ] - yk[ j ] )
                             + ( zk[ i ] - zk[ j ] ) * ( zk[ i ] - zk[ j ] );
                   double r2 = ( rk[ i ] + rk[ j ] ) * ( rk[ i ] + rk[ j ] );

                   if ( d2 <= r2 )
                     {
                       CNeg[ i ] = true;
                       break;
                     }
                 }
             }
       }
   }

  for ( int i = 0; i< numCenters; i++)
   {
     ( *hbondType )[ i ] = 'N';

     if ( ( staticMol && ( atype[ i ] == 'E' ) ) || ( !staticMol && ( atype[ i ] == 'I' ) ) ) continue;

     int ii = 0;
     if ( isdigit( atnames[ i ][ ii ] ) ) ii++;

     if ( ( atnames[ i ][ ii ] == 'O' ) || ( atnames[ i ][ ii ] == 'N' ) ) ( *hbondType )[ i ] = 'A';
     else if ( atnames[ i ][ ii ] == 'H' )
            {
               for ( int j = 0; j < numCenters; j++ )
                 if ( resnum[ j ] == resnum[ i ] )
                   {
                     int jj = 0;
                     if ( isdigit( atnames[ j ][ jj ] ) ) jj++;

                     if ( ( atnames[ j ][ jj ] == 'O' ) || ( atnames[ j ][ jj ] == 'N' ) || ( ( atnames[ j ][ jj ] == 'C' ) && ( CNeg[ j ] == true ) ) )
                       {
                         double d2 = ( xk[ i ] - xk[ j ] ) * ( xk[ i ] - xk[ j ] )
                                   + ( yk[ i ] - yk[ j ] ) * ( yk[ i ] - yk[ j ] )
                                   + ( zk[ i ] - zk[ j ] ) * ( zk[ i ] - zk[ j ] );
                         double r2 = ( rk[ i ] + rk[ j ] ) * ( rk[ i ] + rk[ j ] );

                         if ( d2 <= r2 )
                           {
                             ( *hbondType )[ i ] = 'D';
                             break;
                           }
                       }
                   }
            }
   }

  delete [] CNeg;
}


void preprocessElementInformationTable( int **index, int *nIndex, double *hydrophobicityCutoff, int numTopHydrophobicResidues )
{
  int c = 1;
  int n = sizeof( elementTable ) / sizeof( elementTable[ 0 ] );

  PairingHeap* pairH = new PairingHeap( true, false );

  for ( int i = 1; i < n; i++ )
    if ( strcmp( elementTable[ i ].residueName, elementTable[ i - 1 ].residueName ) )
      {
        c++;
        if ( ( elementTable[ i ].residueIndex >= 1 ) && ( elementTable[ i ].residueIndex <= 20 ) )
           pairH->Insert( elementTable[ i ].residueIndex, elementTable[ i ].perResidueHydrophobicity );
      }

  while ( ( numTopHydrophobicResidues-- > 0 ) && !pairH->isEmpty( ) )
    {
      int k;
      double v;

      pairH->Delete_Min( k, v );

      *hydrophobicityCutoff = v;
    }

  ( *index ) = new int [ c ];

  *nIndex = c;

  c = 0;
  ( *index )[ 0 ] = 0;
  for ( int i = 1; i < n; i++ )
    if ( strcmp( elementTable[ i ].residueName, elementTable[ i - 1 ].residueName ) )
        ( *index )[ ++c ] = i;

  delete pairH;
}


int strcmp_nospace( char *s1, char *s2 )
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



float getHydrophobicity( char *atomName, char *residueName, int *index, int nIndex, bool useInterfacePropensity, bool perResidueHydrophobicity, double hydrophobicityCutoff )
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
                  if ( perResidueHydrophobicity )
                    {
                      if ( useInterfacePropensity ) return ( - ( float ) elementTable[ j ].interfacePropensity );

                      if ( elementTable[ j ].perResidueHydrophobicity > hydrophobicityCutoff ) return 0;
                      else return ( ( float ) elementTable[ j ].perResidueHydrophobicity );
                    }
                  else return ( ( float ) elementTable[ j ].hydrophobicity );
                }

              j++;
            }
        }
    }

  return 0;
}


void assignHydrophobicity( int numCenters, double *x, double *y, double *z,
	                   char **atnames, char **restype, char *atype, float *charges, bool staticMol,
	                   double distCutoff, int *index, int nIndex, bool useInterfacePropensity, bool perResidueHydrophobicity,
	                   double hydrophobicityCutoff, float **hydrophobicity, float **hydrophobicity2 )
{
  ( *hydrophobicity ) = new float [ numCenters ];
  ( *hydrophobicity2 ) = new float [ numCenters ];

  if ( staticMol )
    {
       double trans = 0;

       for ( int i = 0; i < numCenters; i++ )
         {
           if ( x[ i ] < trans ) trans = x[ i ];
           if ( y[ i ] < trans ) trans = y[ i ];
           if ( z[ i ] < trans ) trans = z[ i ];
         }

       PG *skinPG = new PG( 10.0, -trans, 5.0 );
       Point pt;

       for ( int i = 0; i < numCenters; i++ )
         if ( atype[ i ] == 'E' )
           {
             pt.x = x[ i ];
             pt.y = y[ i ];
             pt.z = z[ i ];

             skinPG->addPoint( &pt );
           }

       for ( int i = 0; i < numCenters; i++ )
         if ( atype[ i ] == 'I' )
           {
             pt.x = x[ i ];
             pt.y = y[ i ];
             pt.z = z[ i ];

//             if ( skinPG->pointsWithinRange( &pt, distCutoff ) )
             ( *hydrophobicity )[ i ] = getHydrophobicity( atnames[ i ], restype[ i ], index, nIndex, useInterfacePropensity, perResidueHydrophobicity, hydrophobicityCutoff );
             ( *hydrophobicity2 )[ i ] = getHydrophobicity( atnames[ i ], restype[ i ], index, nIndex, false, perResidueHydrophobicity, hydrophobicityCutoff );
//             else
//                ( *hydrophobicity )[ i ] = 0;
           }
         else ( *hydrophobicity )[ i ] = 0;

       delete skinPG;
    }
  else
    {
       for ( int i = 0; i < numCenters; i++ )
         {
//         if ( atype[ i ] == 'E' )
              ( *hydrophobicity )[ i ] = getHydrophobicity( atnames[ i ], restype[ i ], index, nIndex, useInterfacePropensity, perResidueHydrophobicity, hydrophobicityCutoff );
              ( *hydrophobicity2 )[ i ] = getHydrophobicity( atnames[ i ], restype[ i ], index, nIndex, false, perResidueHydrophobicity, hydrophobicityCutoff );
//         else ( *hydrophobicity )[ i ] = 0;
         }
    }
}



bool readRMSDRef(char *filename, PARAMS_IN *p)
// read a file that provides 0-based indices and corrdinates of atoms in the
// moving molecule that will be used for compting RMSD
//

// the file format is: leading comments with '#'
// number of atoms in moving molecule
// number of atoms used for RMSD
// i x y z for each atom used for RMSD calculation
{
  int numberAtomsTotal, numberAtomsRMSD;
  FILE* fp = fopen( filename, "r" );
  char s[ 2000 ];  // buffer to read lines
  char *val;
  char sep[] = " ";

  if (  fp== NULL ) {
    printf( "Error: Failed to open parameter file %s!\n", filename);
    return false;
  }

  // skip leading comments
  while ( fgets( s, 1999, fp ) != NULL ) {
    if (s[0]=='#') continue;
    break;
  }

  sscanf( s, "%d", &numberAtomsTotal);
  fgets( s, 1999, fp );
  sscanf( s, "%d", &numberAtomsRMSD);

  if (numberAtomsTotal != p->numCentersB) {
    printf( "Error: bad RMSD ref file: %s expected %d got %d!\n",
	    filename, p->numCentersB, numberAtomsTotal);
    return false;
  }

  float *xRef = new float[ numberAtomsRMSD ];
  float *yRef = new float[ numberAtomsRMSD ];
  float *zRef = new float[ numberAtomsRMSD ];
  int *atNums = new int[ numberAtomsRMSD ];

cout<<"reached here "<<numberAtomsRMSD<<endl;

  // read the index X Y Z values
  int i =0;
  while ( fgets( s, 1999, fp ) != NULL ) {
cout<<i<<endl;
    atNums[i] = atoi(strtok(s, sep));  // atom index
    xRef[i] = (float)atof(strtok(NULL, sep)); // X coordinate
    yRef[i] = (float)atof(strtok(NULL, sep)); // y coordinate
    zRef[i] = (float)atof(strtok(NULL, sep)); // z coordinate
cout<<" "<<atNums[i]<<" "<<xRef[i]<<" "<<yRef[i]<<" "<<zRef[i]<<endl;
    i++;
  }

  p->nbRMSDAtoms = numberAtomsRMSD;
  p->xRef = xRef;
  p->yRef = yRef;
  p->zRef = zRef;
  p->atNums = atNums;
  return true;
}


bool readEffGridFile( char *effGridFile, PARAMS_IN *pr )
{
  int numLines = countLines( effGridFile );

  pr->efficientGridSizes = new int[ numLines ];

  FILE* fpGrid = fopen( effGridFile, "r" );

  if (  fpGrid == NULL )
    {
      printf( "Error: Failed to open parameter file %s!\n", effGridFile );
      return false;
    }

  char buf[ 2000 ];
  int i = 0;

  while ( fgets( buf, 1999, fpGrid ) != NULL )
    {
     int j = skipWhiteSpaces( buf, 0 );

     if ( isdigit( buf[ j ] ) )
       {
        int v;
        getInt( buf, j, &v );

        printf( "%d\n", v );

        if ( ( i > 0 ) && ( v <= pr->efficientGridSizes[ i - 1 ] ) )
          {
           printf( "Error: Invalid sequence of grid sizes in parameter file %s!\n", effGridFile );
           return false;
          }

        pr->efficientGridSizes[ i++ ] = v;
       }
    }

  pr->numEfficientGridSizes = i;

  fclose( fpGrid );

  return true;
}



#define returnSpectrumError1( ) { printError( "Invalid spectrum ( %s )!", p->spectrum ); return false; }
#define returnSpectrumError2( ) { freeMem( p->bands ); returnSpectrumError1( ); }

bool decodeSpectrum( PARAMS_IN *p )
{
   if ( p->spectrum == NULL )
     {
       p->numBands = 0;
       p->bands = ( int * ) malloc( 2 * ( p->numBands + 2 ) * sizeof( int ) );

       if ( p->bands == NULL )  returnSpectrumError1( );

       return true;
     }

   int numBands = 0;
   int i = 0;

   while ( p->spectrum[ i ] )
     {
       if ( p->spectrum[ i ] == ':' ) numBands++;
       i++;
     }

   p->bands = ( int * ) malloc( 2 * ( numBands + 2 ) * sizeof( int ) );

   if ( p->bands == NULL )  returnSpectrumError1( );

   int k = 0;
   i = 0;
   numBands = 0;

   while ( p->spectrum[ i ] )
     {
       i = skipWhiteSpaces( p->spectrum, i );

       if ( !p->spectrum[ i ] ) break;

       int v1, v2;

       int j = getInt( p->spectrum, i, &v1 );

       if ( ( i == j ) || ( v1 < 0 ) ) returnSpectrumError2( );

       i = skipWhiteSpaces( p->spectrum, j );

       if ( p->spectrum[ i ] != ':' ) returnSpectrumError2( );

       i = skipWhiteSpaces( p->spectrum, i + 1 );

       if ( !p->spectrum[ i ] ) returnSpectrumError2( );

       j = getInt( p->spectrum, i, &v2 );

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

       i = skipWhiteSpaces( p->spectrum, j );

       if ( p->spectrum[ i ] == '-' ) i++;
     }

   p->numBands = numBands;

   return true;
}


bool setParamFromFile(PARAMS_IN *p, char *paramFile)
{
  char s[ 2000 ], line[ 2000 ];  // buffer to read lines
  char *key, *val;
  char sep[] = " ";
  FILE *fp;
  int nbrot=-1; // variable used to overwrite total number of rotations in file
                // if the param file provides numRot =
  int ival;
  double dval;
  bool bval;

  //printf("BBBBBBBBBBBBpr = %p %d\n", p, p->numThreads);

  p->paramFile = strdup( paramFile );

  if ( paramFile[ 0 ] )
    {
      fp = fopen( paramFile, "r" );

      if (  fp == NULL )
	{
	  printf( "\n\nError: Failed to open parameter file %s!\n\n",
		  paramFile );
	  return false;
	}

      // read a line
      while ( fgets( s, 1999, fp ) != NULL )
	{
	  if (strlen(s)<3 || s[0]=='#') continue;

	  strcpy( line, s );

	  key = strtok(s, sep);
	  val = strtok(NULL, sep);


	  if (val[strlen(val)-1]=='\n')  // remove unix new line
	    val[strlen(val)-1] = '\0';

	  if (val[strlen(val)-1]=='\r') // rmove windows new line
	    val[strlen(val)-1] = '\0';
	
	  printf("key val NUMTHREAD = %s %s %d \n", key, val, p->numThreads);
      
	  if (strcasecmp(key, "id")==0) {
	     p->id = strdup( val );

	  } else if (strcasecmp(key, "staticMolecule")==0) {
	    if ( !readF2D(val, &p->xkAOrig, &p->ykAOrig, &p->zkAOrig,
			  &p->chargesA, &p->radiiA, &p->numCentersA, &p->typeA,
			  &p->atNamesA, &p->resTypesA, &p->resNumsA, &p->staticMoleculePdb, &p->staticMoleculeF2d) ){

	      printf( "Error: failed to read %s!\n", val);
	      return false;
	    }
		else continue;

	  } else if (strcasecmp(key, "movingMolecule")==0) {
	     if ( !readF2D(val, &p->xkBOrig, &p->ykBOrig, &p->zkBOrig,
		     &p->chargesB, &p->radiiB, &p->numCentersB, &p->typeB,
		     &p->atNamesB, &p->resTypesB, &p->resNumsB, &p->movingMoleculePdb, &p->movingMoleculeF2d) ) {

	      printf( "Error: failed to read %s!\n", val);
	      return false;
	    }
		else continue;
	  } else if ( strcasecmp( key, "staticMoleculePQR" ) == 0 ) {
	     p->staticMoleculePQR = strdup( val );
		 continue;
	  } else if ( strcasecmp( key, "movingMoleculePQR" ) == 0 ) {
	     p->movingMoleculePQR = strdup( val );
		 continue;
	  } else if ( strcasecmp( key, "staticMoleculePSF" ) == 0 ) {
	     p->staticMoleculePSF = strdup( val );
		 continue;
	  } else if ( strcasecmp( key, "movingMoleculePSF" ) == 0 ) {
	     p->movingMoleculePSF = strdup( val );
		 continue;
	  } else if ( strcasecmp( key, "staticMoleculeMol2" ) == 0 ) {
	     p->staticMoleculeMol2 = strdup( val );
		 continue;
	  } else if ( strcasecmp( key, "movingMoleculeMol2" ) == 0 ) {
	     p->movingMoleculeMol2 = strdup( val );
		 continue;
	  } else if (strcasecmp(key, "staticMoleculeShapeReal")==0) {
	     p->staticMoleculeSCReRaw = strdup( val );
		 continue;
	  } else if (strcasecmp(key, "staticMoleculeShapeImaginary")==0) {
	     p->staticMoleculeSCImRaw = strdup( val );
		 continue;
	  } else if (strcasecmp(key, "staticMoleculeElecReal")==0) {
	     p->staticMoleculeElecReRaw = strdup( val );
		 continue;
	  } else if (strcasecmp(key, "movingMoleculeShapeReal")==0) {
	     p->movingMoleculeSCReRaw = strdup( val );
		 continue;
	  } else if (strcasecmp(key, "movingMoleculeShapeImaginary")==0) {
	     p->movingMoleculeSCImRaw = strdup( val );
		 continue;
	  } else if (strcasecmp(key, "movingMoleculeElecReal")==0) {
	     p->movingMoleculeElecReRaw = strdup( val );
		 continue;
	  } else if (strcasecmp(key, "prmFile")==0) {
	     p->prmFile = strdup( val );
		 continue;
	  } else if (strcasecmp(key, "aprmFile")==0) {
	     p->aprmFile = strdup( val );
		 continue;
	  } else if (strcasecmp(key, "rtfFile")==0) {
	     p->rtfFile = strdup( val );
		 continue;
	  } else if (strcasecmp(key, "rmsdAtoms")==0) {
	    if ( !readRMSDRef(val, p) ) {

	      printf( "Error: failed to read %s!\n", val);
	      return false;
	    }
		else continue;
	  } else if (strcasecmp(key, "effGridFile")==0) {
	    if ( !readEffGridFile(val, p) ) {

	      printf( "Error: failed to read %s!\n", val);
	      return false;
	    }
		else continue;
	  }
	   else if ( strcasecmp( key, "forbiddenVolFileName" ) == 0 ) {
	     p->forbiddenVolFileName = strdup( val );
	     continue;
	  }

	  else if (strcasecmp(key, "minEffGridSize")==0) {
	  	    ival = atoi(val);
	  	    if ( ival <= 0 )
	  	      {
	  		printf( "Error: minEffGridSize must be a positive integer!\n");
	  		return false;
	  	      }
	  	    p->minEffGridSize = ival;
		    continue;
	  	  }

	  else if (strcasecmp(key, "maxEffGridSize")==0) {
	  	    ival = atoi(val);
	  	    if ( ival <= 0 )
	  	    {
	  		printf( "Error: maxEffGridSize must be a positive integer!\n");
	  		return false;
	  	    }
	  	    p->maxEffGridSize = ival;
 		   continue;
   	  }
	  else if (strcasecmp(key, "numThreads")==0) {
	    ival = atoi(val);
	    if ( ival <= 0 )
	      {
		printf( "Error: %s must be a positive integer!\n", key);
                return false;
	      }
	    p->numThreads = ival;
		continue;

	  } else if (strcasecmp(key, "numFreq")==0) {
	    ival = atoi(val);
	    if ( ival <= 0 )
	      {
		printf( "Error: numFreq must be a positive integer!\n");
		return false;
	      }
	    p->numFreq = ival;
        p->numFreqSpecified = true;
		continue;

	  } else if (strcasecmp(key, "numSolutions")==0) {
	    ival = atoi(val);
	    if ( ival <= 0 && ival != -1 )
	      {
		printf( "Error: numSolutions must be a positive integer!\n");
		return false;
	      }
	    p->numberOfPositions = ival;
		continue;


	  } else if (strcasecmp(key, "gridSize")==0) {
	    ival = atoi(val);
	    if ( ival <= 0 )
	      {
		printf( "Error: gridSize must be a positive integer!\n");
		return false;
	      }
	    p->gridSize = ival;
	    p->gridSizeSpecified = true;
		continue;

	  } else if (strcasecmp(key, "blobbiness")==0) {
	    dval = atof(val);
	    p->blobbiness = dval;
	    p->blobbinessSpecified = true;
	    continue;
	  } else if (strcasecmp(key, "distCutoff")==0) {
	    dval = atof(val);
	    if ( dval <= 0 )
	      {
		printf( "Error: distCutoff must be a positive real value!\n");
		return false;
	      }
	    p->distanceCutoff = dval;
		continue;
	  } else if (strcasecmp(key, "alpha")==0) {
	    dval = atof(val);
	    if ( dval < 1. )
	      {
                printf( "Error: alpha must be a real value at least as large as 1!\n");
                return false;
	      }
	    p->alpha = dval;
		continue;
	  } else if (strcasecmp(key, "skinSkinWeight")==0) {
	    dval = atof(val);
	    p->skinSkinWeight = dval;
		continue;
	  } else if (strcasecmp(key, "coreCoreWeight")==0) {
	    dval = atof(val);
	    p->coreCoreWeight = dval;
		continue;
	  } else if (strcasecmp(key, "skinCoreWeight")==0) {
	    dval = atof(val);
	    p->skinCoreWeight = dval;
		continue;
	  } else if (strcasecmp(key, "realSCWeight")==0) {
	    dval = atof(val);
	    p->realSCWeight = dval;
		continue;
	  } else if (strcasecmp(key, "imaginarySCWeight")==0) {
	    dval = atof(val);
	    p->imaginarySCWeight = dval;

	  } else if (strcasecmp(key, "elecKernelVoidRad")==0) {
	    dval = atof(val);
	    if ( dval < 0 )
	      {
		printf( "Error: %s must be a nonnegative real value!\n", key);
		return false;
	      }
	    p->elecKernelVoidRad = dval;
		continue;
	  } else if (strcasecmp(key, "elecKernelDistLow")==0) {
	    dval = atof(val);
	    if ( dval < 0 )
	      {
		printf( "Error: %s must be a nonnegative real value!\n", key);
		return false;
	      }
	    p->elecKernelDistLow = dval;
		continue;
	  } else if (strcasecmp(key, "elecKernelDistHigh")==0) {
	    dval = atof(val);
	    if ( dval < 0 )
	      {
		printf( "Error: %s must be a nonnegative real value!\n", key);
		return false;
	      }
	    p->elecKernelDistHigh = dval;
		continue;
	  } else if (strcasecmp(key, "elecKernelValLow")==0) {
	    dval = atof(val);
	    if ( dval <= 0 )
	      {
		printf( "Error: %s must be a positive real value!\n", key);
		return false;
	      }
	    p->elecKernelValLow = dval;
		continue;
	  } else if (strcasecmp(key, "elecKernelValHigh")==0) {
	    dval = atof(val);
	    if ( dval <= 0 )
	      {
		printf( "Error: %s must be a positive real value!\n", key);
		return false;
	      }
	    p->elecKernelValHigh = dval;
		continue;
	  } else if (strcasecmp(key, "elecWeight")==0) {
	    dval = atof(val);
	    /*if ( dval < 0 )
	      {
		printf( "Error: elecWeight must be a non-negative real value!\n");
		return false;
	      }*/
	    p->elecScale = dval;
		continue;
	  } else if (strcasecmp(key, "elecRadiusInGrids")==0) {
	    dval = atof(val);
	    if ( dval <= 0 )
	      {
		printf( "Error: %s must be a positive real value!\n", key);
		return false;
	      }
	    p->elecRadiusInGrids = dval;
		continue;
	  } else if (strcasecmp(key, "hbondWeight")==0) {
	    dval = atof(val);
	    p->hbondWeight = dval;
		continue;
	  } else if (strcasecmp(key, "hbondDistanceCutoff")==0) {
	    dval = atof(val);
	    if ( dval <= 0 )
	      {
		printf( "Error: %s must be a positive real value!\n", key);
		return false;
	      }
	    p->hbondDistanceCutoff = dval;
		continue;
	  } else if (strcasecmp(key, "hydrophobicityWeight")==0) {
	    dval = atof(val);
	    p->hydrophobicityWeight = dval;
		continue;
	  } else if (strcasecmp(key, "hydrophobicityProductWeight")==0) {
	    dval = atof(val);
	    p->hydrophobicityProductWeight = dval;
		continue;
	  } else if (strcasecmp(key, "hydroRatioTolerance")==0) {
	    dval = atof(val);
	    if ( dval <= 0 )
	      {
		printf( "Error: %s must be a positive real value!\n", key);
		return false;
	      }
	    p->hydroRatioTolerance = dval;
		continue;
	  } else if (strcasecmp(key, "hydroMinRatio")==0) {
	    dval = atof(val);
	    if ( dval <= 0 )
	      {
		printf( "Error: %s must be a positive real value!\n", key);
		return false;
	      }
	    p->hydroMinRatio = dval;
		continue;
	  } else if (strcasecmp(key, "hydroRatioNumeratorLow")==0) {
	    dval = atof(val);
	    if ( dval <= 0 )
	      {
		printf( "Error: %s must be a positive real value!\n", key);
		return false;
	      }
	    p->hydroRatioNumeratorLow = dval;
		continue;
	  } else if (strcasecmp(key, "hydroRatioNumeratorHigh")==0) {
	    dval = atof(val);
	    if ( dval <= 0 )
	      {
		printf( "Error: %s must be a positive real value!\n", key);
		return false;
	      }
	    p->hydroRatioNumeratorHigh = dval;
		continue;
	  } else if (strcasecmp(key, "hydroRatioDenominatorLow")==0) {
	    dval = atof(val);
	    if ( dval <= 0 )
	      {
		printf( "Error: %s must be a positive real value!\n", key);
		return false;
	      }
	    p->hydroRatioDenominatorLow = dval;
		continue;
	  } else if (strcasecmp(key, "hydroRatioDenominatorHigh")==0) {
	    dval = atof(val);
	    if ( dval <= 0 )
	      {
		printf( "Error: %s must be a positive real value!\n", key);
		return false;
	      }
	    p->hydroRatioDenominatorHigh = dval;
		continue;
	  } else if (strcasecmp(key, "twoWayHydrophobicity")==0) {
	    if (strcasecmp(val, "true")==0 ) p->twoWayHydrophobicity = true;
	    else if (strcasecmp(val, "false")==0 ) p->twoWayHydrophobicity = false;
	    else  {
		   printf( "Error: %s must be a Boolean value!\n", key);
		   return false;
	          }
	  } else if (strcasecmp(key, "hydroPhobicPhobicWeight")==0) {
	    dval = atof(val);
	    p->hydroPhobicPhobicWeight = dval;

	  } else if (strcasecmp(key, "hydroPhilicPhilicWeight")==0) {
	    dval = atof(val);
	    p->hydroPhilicPhilicWeight = dval;

	  } else if (strcasecmp(key, "hydroPhobicPhilicWeight")==0) {
	    dval = atof(val);
	    p->hydroPhobicPhilicWeight = dval;

	  } else if (strcasecmp(key, "hydroRadExt")==0) {
	    dval = atof(val);
	    if ( dval < 0 )
	      {
		printf( "Error: %s must be a nonnegative real value!\n", key);
		return false;
	      }
	    p->hydroRadExt = dval;

	  } else if (strcasecmp(key, "useInterfacePropensity")==0) {
	    if (strcasecmp(val, "true")==0 ) p->useInterfacePropensity = true;
	    else if (strcasecmp(val, "false")==0 ) p->useInterfacePropensity = false;
	    else  {
		   printf( "Error: %s must be a Boolean value!\n", key);
		   return false;
	          }
	  } else if (strcasecmp(key, "perResidueHydrophobicity")==0) {
	    if (strcasecmp(val, "true")==0 ) p->perResidueHydrophobicity = true;
	    else if (strcasecmp(val, "false")==0 ) p->perResidueHydrophobicity = false;
	    else  {
		   printf( "Error: %s must be a Boolean value!\n", key);
		   return false;
	          }
	  } else if (strcasecmp(key, "staticMolHydroDistCutoff")==0) {
	    dval = atof(val);
	    if ( dval < 0 )
	      {
		printf( "Error: %s must be a nonnegative real value!\n", key);
		return false;
	      }
	    p->staticMolHydroDistCutoff = dval;

	  } if (strcasecmp(key, "numTopHydrophobicResidues")==0) {
	    ival = atoi(val);
	    if ( ival <= 0 )
	      {
		printf( "Error: %s must be a positive integer!\n", key);
                return false;
	      }
	    p->numTopHydrophobicResidues = ival;

	  } else if (strcasecmp(key, "simpleShapeWeight")==0) {
	    dval = atof(val);
	    if ( dval < 0 )
	      {
		printf( "Error: %s must be a nonnegative real value!\n", key);
		return false;
	      }
	    p->simpleShapeWeight = dval;

	  } else if (strcasecmp(key, "simpleChargeWeight")==0) {
	    dval = atof(val);
	    if ( dval < 0 )
	      {
		printf( "Error: %s must be a nonnegative real value!\n", key);
		return false;
	      }
	    p->simpleChargeWeight = dval;

	  } else if (strcasecmp(key, "simpleRadExt")==0) {
	    dval = atof(val);
	    if ( dval < 0 )
	      {
		printf( "Error: %s must be a nonnegative real value!\n", key);
		return false;
	      }
	    p->simpleRadExt = dval;

	  } else if (strcasecmp(key, "clashWeight")==0) {
	    dval = atof(val);
	    p->clashWeight = dval;

	  }  else if (strcasecmp(key, "pseudoAtomRadius")==0) {
	    dval = atof(val);
	    if ( dval < 0 )
	      {
		printf( "Error: %s must be a nonnegative real value!\n", key);
		return false;
	      }
	    p->pseudoAtomRadius = dval;

	  } else if (strcasecmp(key, "pseudoAtomDistance")==0) {
	    dval = atof(val);
	    if ( dval < 0 )
	      {
		printf( "Error: %s must be a nonnegative real value!\n", key);
		return false;
	      }
	    p->pseudoAtomDistance = dval;

	  } else if ( strcasecmp( key, "initialRotation" ) == 0 ) {

	    float rotMat[ 9 ];

	    if ( sscanf( "%s %f %f %f %f %f %f %f %f %f", key,
	                 &( rotMat[ 0 ] ), &( rotMat[ 1 ] ), &( rotMat[ 2 ] ),
	                 &( rotMat[ 3 ] ), &( rotMat[ 4 ] ), &( rotMat[ 5 ] ),
	                 &( rotMat[ 6 ] ), &( rotMat[ 7 ] ), &( rotMat[ 8 ] ) ) != 10 )
	      {
		printf( "Error: failed to read %s!\n", key );
		return false;
	      }

	    p->initRot.set( rotMat[ 0 ], rotMat[ 1 ], rotMat[ 2 ], 0.0,
	                    rotMat[ 3 ], rotMat[ 4 ], rotMat[ 5 ], 0.0,
	                    rotMat[ 6 ], rotMat[ 7 ], rotMat[ 8 ], 0.0,
	                            0.0,         0.0,         0.0, 1.0 );

	  } else if (strcasecmp(key, "pruneAngle")==0) {
	    dval = atof(val);
	    if ( dval < 0 )
	      {
		printf( "Error: %s must be a nonnegative real value!\n", key);
		return false;
	      }
	    p->pruneAngle = dval;

	  } else if (strcasecmp(key, "smoothSkin")==0) {
	    if (strcasecmp(val, "true")==0 ) p->smoothSkin = true;
	    else if (strcasecmp(val, "false")==0 ) p->smoothSkin = false;
	    else  {
		   printf( "Error: %s must be a Boolean value!\n", key);
		   return false;
	          }
	  } else if ( strcasecmp( key, "singleLayerLigandSkin" ) == 0 )
	           {
	             if ( !strcasecmp( val, "true" ) ) p->singleLayerLigandSkin = true;
	             else if ( !strcasecmp( val, "false" ) ) p->singleLayerLigandSkin = false;
	                  else
	                      {
		                printf( "Error: %s must be a Boolean value!\n", key );
		                return false;
	                      }
	  } else if (strcasecmp(key, "rotateVolume")==0) {
	    if (strcasecmp(val, "true")==0 ) p->rotateVolume = true;
	    else if (strcasecmp(val, "false")==0 ) p->rotateVolume = false;
	    else  {
		   printf( "Error: %s must be a Boolean value!\n", key);
		   return false;
	          }
	  } else if (strcasecmp(key, "dockVolume")==0) {
	    if (strcasecmp(val, "true")==0 ) p->dockVolume = true;
	    else if (strcasecmp(val, "false")==0 ) p->dockVolume = false;
	    else  {
		   printf( "Error: %s must be a Boolean value!\n", key);
		   return false;
	          }
	  } else if (strcasecmp(key, "useSparseFFT")==0) {
	    if (strcasecmp(val, "true")==0 ) p->useSparseFFT = true;
	    else if (strcasecmp(val, "false")==0 ) p->useSparseFFT = false;
	    else  {
		   printf( "Error: %s must be a Boolean value!\n", key);
		   return false;
	          }
	  } else if (strcasecmp(key, "narrowBand")==0) {
	    if (strcasecmp(val, "true")==0 ) p->narrowBand = true;
	    else if (strcasecmp(val, "false")==0 ) p->narrowBand = false;
	    else  {
		   printf( "Error: %s must be a Boolean value!\n", key);
		   return false;
	          }
	  } else if (strcasecmp(key, "gridSpacing")==0) {
	    dval = atof(val);
	    if ( dval < 0 )
	      {
		printf( "Error: %s must be a positive real value!\n", key );
		return false;
	      }
	    p->gridSpacing = dval;
	    p->gridSpacingSpecified = true;

	  } else if (strcasecmp(key, "enforceExactGridSpacing")==0) {
	    if (strcasecmp(val, "true")==0 ) p->enforceExactGridSpacing = true;
	    else if (strcasecmp(val, "false")==0 ) p->enforceExactGridSpacing = false;
	    else  {
		   printf( "Error: %s must be a Boolean value!\n", key);
		   return false;
	          }

	  } else if (strcasecmp(key, "interpFuncExtentInAngstroms")==0) {
	    dval = atof(val);
	    if ( dval < 0 )
	      {
		printf( "Error: %s must be a positive real value!\n", key );
		return false;
	      }
	    p->interpFuncExtentInAngstroms = dval;
	  } else if (strcasecmp(key, "scoreScaleUpFactor")==0) {
	    dval = atof(val);
	    if ( dval < 0 )
	      {
		printf( "Error: scoreScaleUpFactor must be a non-negative real value!\n");
		return false;
	      }
	    p->scoreScaleUpFactor = dval;

	  } else if (strcasecmp(key, "rotFile")==0) {
	    readRotations(val, p);

	  } else if (strcasecmp(key, "numRot")==0) {
	    ival = atoi(val);
	    if ( ival < 0 )
	      {
		printf( "Error: %s must be a non-negative integer!", key );
		return false;
	      }
	    nbrot = ival;

	  } else if (strcasecmp(key, "outFile")==0) {
	    p->outputFilename = strdup(val);

	  } else if (strcasecmp(key, "transFile")==0) {
	    p->transformationFilename = strdup(val);

	  } else if (strcasecmp(key, "vdWGridSize")==0) {
	    ival = atoi(val);
	    if ( ival < 0 )
	      {
		printf( "Error: %s must be a positive integer!", key );
		return false;
	      }
	    p->vdWGridSize = ival;

	  } else if (strcasecmp(key, "compQuadVdW")==0) {
	    if (strcasecmp(val, "true")==0 ) p->compQuadVdW = true;
	    else if (strcasecmp(val, "false")==0 ) p->compQuadVdW = false;
	    else  {
		   printf( "Error: %s must be a Boolean value!\n", key);
		   return false;
	          }

	  } else if (strcasecmp(key, "surfaceBasedVdW")==0) {
	    if (strcasecmp(val, "true")==0 ) p->surfaceBasedVdW = true;
	    else if (strcasecmp(val, "false")==0 ) p->surfaceBasedVdW = false;
	    else  {
		   printf( "Error: %s must be a Boolean value!\n", key);
		   return false;
	          }

	  } else if (strcasecmp(key, "applyVdWFilter")==0) {
	    if (strcasecmp(val, "true")==0 ) p->applyVdWFilter = true;
	    else if (strcasecmp(val, "false")==0 ) p->applyVdWFilter = false;
	    else  {
		   printf( "Error: %s must be a Boolean value!\n", key);
		   return false;
	          }

	  } else if (strcasecmp(key, "vdWCutoffLow")==0) {
	    p->vdWCutoffLow = atof(val);
	  } else if (strcasecmp(key, "vdWCutoffHigh")==0) {
	    p->vdWCutoffHigh = atof(val);
	  } else if (strcasecmp(key, "vdWWellWidth")==0) {
	    dval = atof(val);
	    if ( dval < 0 )
	      {
		printf( "Error: %s must be a non-negative real value!\n", key );
		return false;
	      }
	    p->vdWWellWidth = dval;

	  } else if (strcasecmp(key, "vdWEqmRadScale")==0) {
	    dval = atof(val);
	    if ( dval <= 0 )
	      {
		printf( "Error: %s must be a positive real value!\n", key );
		return false;
	      }
	    p->vdWEqmRadScale = dval;

	  } else if (strcasecmp(key, "vdWTolerance")==0) {
	    dval = atof(val);
	    if ( dval < 0 )
	      {
		printf( "Error: %s must be a non-negative real value!\n", key );
		return false;
	      }
	    p->vdWTolerance = dval;

	  } else if (strcasecmp(key, "applyClashFilter")==0) {
	    if (strcasecmp(val, "true")==0 ) p->applyClashFilter = true;
	    else if (strcasecmp(val, "false")==0 ) p->applyClashFilter = false;
	    else  {
		   printf( "Error: %s must be a Boolean value!\n", key);
		   return false;
	          }

	  } 
	  else if (strcasecmp(key, "applyForbiddenVolumeFilter")==0) {

	    if (strcasecmp(val, "true")==0 ) p->applyForbiddenVolumeFilter = true;
	    else if (strcasecmp(val, "false")==0 ) p->applyForbiddenVolumeFilter = false;
	    else  {
		   printf( "Error: %s must be a Boolean value!\n", key);
		   return false;
	          }

	  } 

	  else if (strcasecmp(key, "forbiddenVolumeFileType")==0) {
	    int aval = atoi(val);
	    if(aval < 0 || aval > 2) 
	      std::cout<<"Unable to read file type. Forbidden volume filter will not be applied "<<std::endl;
	      	    
	    p->forbiddenVolumeFileType = aval;

	  } 

	  else if (strcasecmp(key, "applyMiscFilter")==0) {
	    if (strcasecmp(val, "true")==0 ) p->applyMiscFilter = true;
	    else if (strcasecmp(val, "false")==0 ) p->applyMiscFilter = false;
	    else  {
		   printf( "Error: %s must be a Boolean value!\n", key);
		   return false;
	          }

	  } else if (strcasecmp(key, "eqmDistFrac")==0) {
	    dval = atof(val);
	    if ( dval < 0 )
	      {
		printf( "Error: %s must be a non-negative real value!\n", key );
		return false;
	      }
	    p->eqmDistFrac = dval;

	  } else if (strcasecmp(key, "clashTolerance")==0) {
	    ival = atoi(val);
	    if ( ival < 0 )
	      {
		printf( "Error: %s must be a non-negative integer!\n", key );
		return false;
	      }
	    p->clashTolerance = ival;

	  } else if (strcasecmp(key, "applyPseudoGsolFilter")==0) {
	    if (strcasecmp(val, "true")==0 ) p->applyPseudoGsolFilter = true;
	    else if (strcasecmp(val, "false")==0 ) p->applyPseudoGsolFilter = false;
	    else  {
		   printf( "Error: %s must be a Boolean value!\n", key);
		   return false;
	          }

	  } else if (strcasecmp(key, "pseudoGsolCutoff")==0) {
	    p->pseudoGsolCutoff = atof(val);

	  }  else if (strcasecmp(key, "pseudoGsolWeight")==0) {
	    p->pseudoGsolWeight = atof(val);

	  } else if (strcasecmp(key, "pseudoGsolFilterLowestRank")==0) {
	    ival = atoi(val);
	    if ( ival < 1 )
	      {
		printf( "Error: %s must be a positive integer!\n", key );
		return false;
	      }
	    p->pseudoGsolFilterLowestRank = ival;

	  } else if (strcasecmp(key, "applyDispersionFilter")==0) {
	    if (strcasecmp(val, "true")==0 ) p->applyDispersionFilter = true;
	    else if (strcasecmp(val, "false")==0 ) p->applyDispersionFilter = false;
	    else  {
		   printf( "Error: %s must be a Boolean value!\n", key);
		   return false;
	          }

	  } else if (strcasecmp(key, "dispersionCutoff")==0) {
	    p->dispersionCutoff = atof(val);

	  }  else if (strcasecmp(key, "dispersionWeight")==0) {
	    p->dispersionWeight = atof(val);

	  } else if (strcasecmp(key, "filterDepth")==0) {
	    ival = atoi(val);
	    if ( ival <= 0 )
	      {
		printf( "Error: %s must be a positive integer!\n", key );
		return false;
	      }
	    p->filterDepth = ival;

	  } else if (strcasecmp(key, "filterScaleDownFactor")==0) {
	    p->filterScaleDownFactor = atof(val);

	  } else if (strcasecmp(key, "dispersionEnergyLimit")==0) {
	    dval = atof(val);
	    if ( dval < 0 )
	      {
		printf( "Error: %s must be a non-negative real value!\n", key );
		return false;
	      }
	    p->dispersionEnergyLimit = dval;

	  } else if (strcasecmp(key, "dispersionMinAtomRadius")==0) {
	    dval = atof(val);
	    if ( dval < 0 )
	      {
		printf( "Error: %s must be a non-negative real value!\n", key );
		return false;
	      }
	    p->dispersionMinAtomRadius = dval;

	  } else if (strcasecmp(key, "clusterTransRad")==0) {
	    dval = atof(val);
	    if ( dval < 0 )
	      {
		printf( "Error: %s must be a non-negative real value!\n", key );
		return false;
	      }
	    p->clusterTransRad = dval;

	  } else if (strcasecmp(key, "clusterTransSize")==0) {
	    ival = atoi(val);
	    if ( ival < 1 )
	      {
		printf( "Error: %s must be a positive integer!\n", key );
		return false;
	      }
	    p->clusterTransSize = ival;

	  } else if (strcasecmp(key, "clusterRotRad")==0) {
	    dval = atof(val);
	    if ( dval < 0 )
	      {
		printf( "Error: %s must be a non-negative real value!\n", key );
		return false;
	      }
	    p->clusterRotRad = dval;

	  } else if (strcasecmp(key, "peaksPerRotation")==0) {
	    ival = atoi(val);
	    if ( ival <= 0 )
	      {
		printf( "Error: %s must be a positive integer!\n", key );
		return false;
	      }
	    p->peaksPerRotation = ival;

	  } else if (strcasecmp(key, "breakDownScores")==0) {
	    if (strcasecmp(val,"true")==0) p->breakDownScores = 1;
	    else p->breakDownScores = 0;

	  } else if (strcasecmp(key, "bandwidth")==0) {
	    dval = atof(val);
	    if ( dval < 0 )
	      {
		printf( "Error: bandwidth must be a non-negative real value!\n");
		return false;
	      }
	    p->bandwidth = dval;

	  } else if (strcasecmp(key, "gradFactor")==0) {
	    dval = atof(val);
	    if ( dval < 0 )
	      {
		printf( "Error: gradFactor must be a non-negative real value!\n");
		return false;
	      }
	    p->gradFactor = dval;

	  } else if ( strcasecmp( key, "curvatureWeightedStaticMol" ) == 0 )
	    {
	      if ( strcasecmp( val, "true" ) == 0 ) p->curvatureWeightedStaticMol = true;
	      else if ( strcasecmp( val, "false" ) == 0 ) p->curvatureWeightedStaticMol = false;
	      else  {
		     printf( "Error: %s must be a Boolean value!\n", key );
		     return false;
	            }

	  } else if ( strcasecmp( key, "curvatureWeightedMovingMol" ) == 0 )
	    {
	      if ( strcasecmp( val, "true" ) == 0 ) p->curvatureWeightedMovingMol = true;
	      else if ( strcasecmp( val, "false" ) == 0 ) p->curvatureWeightedMovingMol = false;
	      else  {
		     printf( "Error: %s must be a Boolean value!\n", key );
		     return false;
	            }

	  } else if (strcasecmp(key, "curvatureWeightingRadius")==0) {
	    dval = atof(val);
	    if ( dval < 0 )
	      {
		printf( "Error: %s must be a non-negative real value!\n", key );
		return false;
	      }
	    p->curvatureWeightingRadius = dval;

	  } else if (strcasecmp(key, "spreadReceptorSkin")==0) {
	    if (strcasecmp(val, "true")==0 ) p->spreadReceptorSkin = true;
	    else if (strcasecmp(val, "false")==0 ) p->spreadReceptorSkin = false;
	    else  {
		   printf( "Error: spreadReceptorSkin must be a Boolean value!\n");
		   return false;
	          }

	  } else if (strcasecmp(key, "randomRotate")==0) {
	    if (strcasecmp(val, "true")==0 ) p->randomRotate = true;
	    else if (strcasecmp(val, "false")==0 ) p->randomRotate = false;
	    else  {
		   printf( "Error: randomRotate must be a Boolean value!\n");
		   return false;
	          }

	  } else if (strcasecmp(key, "spectrum")==0) {
	     p->spectrum = strdup( val );
	  } else if ( strcasecmp( key, "numRerank" ) == 0 ) {
	    ival = atoi(val);
	    if ( ival < 0 )
	      {
		printf( "Error: %s must be a nonegative integer!\n", key );
		return false;
	      }
	    p->numRerank = ival;

	  } else if ( strcasecmp( key, "rerank" ) == 0 ) 
	    {
	      if ( strcasecmp( val, "true" ) == 0 ) p->rerank = true;
	      else if ( strcasecmp( val, "false" ) == 0 ) p->rerank = false;
	      else {
		     printf( "Error: %s must be a Boolean value!\n", key);
		     return false;
	            }

	  } else if ( strcasecmp( key, "applyAntibodyFilter" ) == 0 ) 
	    {
	      if ( strcasecmp( val, "true" ) == 0 ) p->applyAntibodyFilter = true;
	      else if ( strcasecmp( val, "false" ) == 0 ) p->applyAntibodyFilter = false;
	      else {
		     printf( "Error: %s must be a Boolean value!\n", key);
		     return false;
	            }

	  } else if ( strcasecmp( key, "applyEnzymeFilter" ) == 0 ) 
	    {
	      if ( strcasecmp( val, "true" ) == 0 ) p->applyEnzymeFilter = true;
	      else if ( strcasecmp( val, "false" ) == 0 ) p->applyEnzymeFilter = false;
	      else {
		     printf( "Error: %s must be a Boolean value!\n", key);
		     return false;
	            }

	  } else if ( strcasecmp( key, "applyResidueContactFilter" ) == 0 ) 
	    {
	      if ( strcasecmp( val, "true" ) == 0 ) p->applyResidueContactFilter = true;
	      else if ( strcasecmp( val, "false" ) == 0 ) p->applyResidueContactFilter = false;
	      else {
		     printf( "Error: %s must be a Boolean value!\n", key);
		     return false;
	            }

	  } 

#ifdef LIBMOL_FOUND
	    else if ( strcasecmp( key, "applyHbondFilter" ) == 0 ) 
	    {
	      if ( strcasecmp( val, "true" ) == 0 ) p->applyHbondFilter = true;
	      else if ( strcasecmp( val, "false" ) == 0 ) p->applyHbondFilter = false;
	      else {
		     printf( "Error: %s must be a Boolean value!\n", key);
		     return false;
	            }

	  } else if ( strcasecmp( key, "hBondFilterWeight" ) == 0 ) {
	    dval = atof(val);
	    p->hBondFilterWeight = dval;

	  }
#endif
	    else if ( strcasecmp( key, "rerankerPseudoGsolWeight" ) == 0 ) {
	    dval = atof(val);
	    p->rerankerPseudoGsolWeight = dval;

	  } else if ( strcasecmp( key, "rerankerDispersionWeightHigh" ) == 0 ) {
	    dval = atof(val);
	    p->rerankerDispersionWeightHigh = dval;

	  } else if ( strcasecmp( key, "rerankerDispersionWeightLow" ) == 0 ) {
	    dval = atof(val);
	    p->rerankerDispersionWeightLow = dval;

	  } else if ( strcasecmp( key, "rerankerF2DockScoreWeight" ) == 0 ) {
	    dval = atof(val);
	    p->rerankerF2DockScoreWeight = dval;

	  } else if ( strcasecmp( key, "rerankerMinF2DockRank" ) == 0 ) {
	    ival = atoi(val);
	    if ( ival < 0 )
	      {
		printf( "Error: %s must be a nonegative integer!\n", key );
		return false;
	      }
	    p->rerankerMinF2DockRank = ival;

	  } /*else if (strcasecmp(key, "vdwSmoothWidth")==0) {
	    dval = atof(val);
	    if ( dval < 0 )
	      {
		printf( "Error: vdwSmoothWidth must be a non-negative real value!\n");
		return false;
	      }
	    p->vdwSmoothWidth = dval;

	  }*/
	  else {
	    printf("WARNING: the following line from the parameter has been ignored\n  %s\n", s);
	  }

	}
      fclose( fp );
    }

  if (p->rotations == NULL) { // no user specified rotation file
    p->rotations = new float[ 9 ];
    p->rotations[ 0 ] = 1;
    p->rotations[ 1 ] = 0;
    p->rotations[ 2 ] = 0;
    p->rotations[ 3 ] = 0;
    p->rotations[ 4 ] = 1;
    p->rotations[ 5 ] = 0;
    p->rotations[ 6 ] = 0;
    p->rotations[ 7 ] = 0;
    p->rotations[ 8 ] = 1;
    p->numberOfRotations = 1;
  }

  if ( p->gridSize < p->numFreq )
    {
      printf( "Error: gridSize cannot be smaller than numFreq!\n");
      return false;
    }

  if ( p->elecKernelVoidRad > p->elecKernelDistLow )
    {
      printf( "Error: elecKernelVoidRad cannot be larger than elecKernelDistLow!\n");
      return false;
    }

  if ( p->elecKernelDistLow >= p->elecKernelDistHigh )
    {
      printf( "Error: elecKernelDistLow must be smaller than elecKernelDistHigh!\n");
      return false;
    }

//  if ( p->breakDownScores &&  p->dockVolume )
//    {
//     printf( "Error: Both breakDownScores and dockVolume cannot be set to \'true\' simultaneously!\n" );
//     return false;
//    }
//
//  if ( p->breakDownScores &&  p->rotateVolume )
//    {
//     printf( "Error: Both breakDownScores and rotateVolume cannot be set to \'true\' simultaneously!\n" );
//     return false;
//    }

//  if ( ( p->elecScale != 0 ) &&  p->rotateVolume )
//    {
//     printf( "Error: elecWeight must be zero when rotateVolume is set to \'true\'!\n" );
//     return false;
//    }

  // if user specified #rotations and
  if ( nbrot != -1 && nbrot < p->numberOfRotations )
    {
      p->numberOfRotations = nbrot;
    }

  if ( p->numThreads > p->numberOfRotations )
    {
      p->numThreads = p->numberOfRotations;
    }

//  if ( p->peaksPerRotation <= 0 ) p->peaksPerRotation = p->numFreq * p->numFreq * p->numFreq;

  if ( !p->blobbinessSpecified )
    {
     if ( p->smoothSkin ) p->blobbiness = -2.3;
     else p->blobbiness = -1.9;
    }

  if ( p->dockVolume &&  !p->rotateVolume )
    {
     printf( "Warning: Since dockVolume is set to \'true\', rotateVolume has also been set to \'true\'!\n" );
     p->rotateVolume = true;
    }

  if ( p->pseudoAtomRadius >= 0 )
    {
      for ( int j = 0; j < p->numCentersA; j++ )
        if ( p->typeA[ j ] == 'E' ) p->radiiA[ j ] = p->pseudoAtomRadius;
    }

  if ( p->pruneAngle > 0 ) createRotationsNeighborhoodGraph( p );

  if ( p->hbondWeight != 0 )
    {
     identifyHbondAcceptorsAndDonors( p->numCentersA, p->xkAOrig, p->ykAOrig, p->zkAOrig, p->radiiA,
	                              p->atNamesA, p->resNumsA, &p->hbondTypeA, p->typeA, true );
     identifyHbondAcceptorsAndDonors( p->numCentersB, p->xkBOrig, p->ykBOrig, p->zkBOrig, p->radiiB,
	                              p->atNamesB, p->resNumsB, &p->hbondTypeB, p->typeB, false );
    }

  if ( ( p->hydrophobicityWeight != 0 ) || ( p->hydroPhobicPhobicWeight != 0 ) || ( p->hydroPhilicPhilicWeight != 0 ) || ( p->hydroPhobicPhilicWeight != 0 )
    || p->applyPseudoGsolFilter )
    {
     int *index;
     int nIndex;
     double hydrophobicityCutoff;

     preprocessElementInformationTable( &index, &nIndex, &hydrophobicityCutoff, p->numTopHydrophobicResidues );

     assignHydrophobicity( p->numCentersA, p->xkAOrig, p->ykAOrig, p->zkAOrig,
	                   p->atNamesA, p->resTypesA, p->typeA, p->chargesA, true,
	                   6.0, index, nIndex, p->useInterfacePropensity, p->perResidueHydrophobicity,
	                   hydrophobicityCutoff, &( p->hydrophobicityA ), &( p->hydrophobicity2A ) );

     assignHydrophobicity( p->numCentersB, p->xkBOrig, p->ykBOrig, p->zkBOrig,
	                   p->atNamesB, p->resTypesB, p->typeB, p->chargesB, false,
	                   6.0, index, nIndex, p->useInterfacePropensity, p->perResidueHydrophobicity,
	                   hydrophobicityCutoff, &( p->hydrophobicityB ), &( p->hydrophobicity2B ) );
    }

  if ( p->randomRotate )
    {
      generateRandomRotationMatrix( p->initRot );
    }

  if ( !decodeSpectrum( p ) ) return false;
  
  return true;
}



bool getComplexType( PARAMS_IN *p, char *paramFile )
{
  char s[ 2000 ];
  char key[ 500 ], val[ 500 ];
  FILE *fp;

  p->paramFile = strdup( paramFile );
  p->complexType = 'U';
  p->antibody = p->enzyme = 0;

  fp = fopen( paramFile, "r" );

  if ( fp == NULL )
    {
      printError( "Failed to open parameter file %s!", paramFile );
      return false;
    }

  while ( fgets( s, 1999, fp ) != NULL )
    {
      if ( sscanf( s, "%s %s", key, val ) != 2 ) continue;

      if ( !strcasecmp( key, "complexType" ) )
         {
           if ( ( val[ 0 ] == 'A' ) || ( val[ 0 ] == 'E' ) || ( val[ 0 ] == 'G' ) || ( val[ 0 ] == 'U' ) ) p->complexType = val[ 0 ];
         }
      else if ( !strcasecmp( key, "staticMoleculePQR" ) ) p->staticMoleculePQR = strdup( val );
      else if ( !strcasecmp( key, "movingMoleculePQR" ) ) p->movingMoleculePQR = strdup( val );
    }

  fclose( fp );
 
  if ( ( p->complexType == 'U' ) || ( p->complexType == 'A' ) )
    {
      if ( isAntibody( p->staticMoleculePQR ) )
         {
           p->complexType == 'A';
           p->antibody = 1;
         }
      else if ( isAntibody( p->movingMoleculePQR ) )
              {
                p->complexType == 'A';
                p->antibody = -1;
              }    
    }
    
  if ( ( p->complexType == 'U' ) || ( p->complexType == 'E' ) )
    {
      int c, t;
    
      if ( countResidues( p->staticMoleculePQR, GLY, &c, &t ) && ( t > 200 ) && ( c > 0.08 * t ) )
         {
           p->complexType == 'E';
           p->enzyme = 1;
         }
      else if ( countResidues( p->movingMoleculePQR, GLY, &c, &t ) && ( t > 200 ) && ( c > 0.08 * t ) )
              {
                p->complexType == 'E';
                p->enzyme = -1;
              }
    }
    
  if ( p->complexType == 'U' ) p->complexType = 'G';

  return true;
}



int main( int argc, char* argv[] )
{
  char *fixedMolFileName, *movingMolFileName, paramFileName[256];

  if (argc<2 || argc>3 ) {
    printf("Usage: F2Dock -score|saveGrid|vdw|effGridFile parameterFile\n");
    return(1);
  }

  if ( argc == 2 ) {
    strcpy(paramFileName, argv[1]);
  } else if (argc == 3) {
    strcpy(paramFileName, argv[2]);
  }

  srand( time( NULL ) );

  // data structure used to pass parameters to DockingMain
  PARAMS_IN pr;
  
  getComplexType( &pr, paramFileName );  
  
// computevdw and vdwSmoothWidth were commented out out since the vdw scores computed were not believeable and did not compare to a python implementation that was using the same libraries.  Actually, computevdw was never defined in the struc may need to do that if vdw score is to be computed in the future.
  // initialize the structure with default parameters
  pr.id = NULL;
  pr.performDocking = 1;
  pr.numThreads = 4;
  pr.breakDownScores = 0;
  pr.numberOfPositions = 20000;
  pr.gridSize = 256;
  pr.gridSizeSpecified = false;
  pr.gridSpacing = 1.2;
  pr.enforceExactGridSpacing = false;
  pr.gridSpacingSpecified = false;
  pr.interpFuncExtentInAngstroms = 0;
  pr.numFreq = 64;
  pr.numFreqSpecified = false;
  pr.smoothSkin = false;
  pr.singleLayerLigandSkin = false;  // if set to 'true', only the surface atoms of the ligand are considered as skin atoms
  pr.pseudoAtomRadius = 1.1;  // if set to a non-negative value, this value overrides the receptor pseudo atoms radii read from the F2D file
  pr.pseudoAtomDistance = 1.7; // distnace of the pseudo atom centers from the vdW surface of the static molecule
  pr.rotateVolume = true;    // by default, rotate the volume instead of the atoms of molecule B
  pr.dockVolume = false;     // by default, dock molecules with explicitly specified atoms instead of volumes

  // initial rotation matrix ( default is < 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 > unless either explicitly set or randomRotate is 'true' )
  pr.initRot.reset( );

  pr.useSparseFFT = true;    // speed up FFT using the sparsity of input and/or output matrices
  pr.narrowBand = true;     // if 'true', consider only the solutions within a narrow band in the output grid

//  pr.gridFile = NULL;
  pr.numEfficientGridSizes = 0;
  pr.efficientGridSizes = NULL;
  pr.minEffGridSize = 0;
  pr.maxEffGridSize = 0;

//  pr.interpFuncExtent = 3;
  pr.numCentersA = 0;
  pr.numCentersB = 0;
  pr.numberOfRotations = 1000000;

  pr.distanceCutoff = 5.0;
  pr.alpha = 1;
  pr.blobbiness = -2.3;
  pr.blobbinessSpecified = false;
  pr.skinSkinWeight = 0.57;
  pr.coreCoreWeight = 5.0;
  pr.skinCoreWeight = -0.23;
  pr.realSCWeight = 1.0;
  pr.imaginarySCWeight = 0.0;

  pr.elecKernelVoidRad = 0.0;
  pr.elecKernelDistLow = 6.0;
  pr.elecKernelDistHigh = 8.0;
  pr.elecKernelValLow = 4.0;
  pr.elecKernelValHigh = 80.0;

  pr.elecScale = 0.72;
  pr.elecRadiusInGrids = 2.9;  // the radius of the sphere inside which a charge is diffused using a Gaussian

  pr.hbondWeight = 0.0;
  pr.hbondDistanceCutoff = 2.0;  // the distance cutoff ( in angstroms ) for hydrogen bonds

  pr.hydrophobicityWeight = 8.5;
  pr.hydrophobicityProductWeight = 0.001;
  pr.hydroRatioTolerance = 10.0;
  pr.hydroMinRatio = 0.5;  
  pr.hydroRatioNumeratorLow = 2.0;
  pr.hydroRatioNumeratorHigh = 100.0;  
  pr.hydroRatioDenominatorLow = 0.2;
  pr.hydroRatioDenominatorHigh = 6.0;  
  pr.twoWayHydrophobicity = true;
  pr.hydroPhobicPhobicWeight = 0;
  pr.hydroPhilicPhilicWeight = 0;
  pr.hydroPhobicPhilicWeight = 0;
  pr.hydroRadExt = 1.5;  // in Angstroms
  pr.useInterfacePropensity = true;    // instead of just hydrophobic residues:
                                       // S. Jones & J. M. Thornton, Analysis of Protein-Protein Interaction Sites using Surface Patches,
	                               // JMB 272, pp. 121-132, 1997
  pr.perResidueHydrophobicity = true;
  pr.numTopHydrophobicResidues = 100;  // only the residues with the topmost `numTopHydrophobicResidues' hydrophobicity
                                       // values will be considered
  pr.staticMolHydroDistCutoff = 4.0;   // distance cutoff from the skin atom centers for the atoms of the static
                                       // molecule which contribute to hydrophobicity computation
                                       // (this is a center-to-surface distance cutoff)

  pr.simpleShapeWeight = 0.000001;    // simplified shape complementarity
  pr.simpleChargeWeight = 2.0;        // simplified charge complementarity
  pr.simpleRadExt = 1.5;              // in Angstroms

  pr.scoreScaleUpFactor = 10000;

  pr.bandwidth = 2;
  pr.gradFactor = 1.1;

  pr.curvatureWeightedStaticMol = true;   // if 'true', construct curvature-weighted receptor skin
  pr.curvatureWeightedMovingMol = false;   // if 'true', construct curvature-weighted ligand skin
  pr.curvatureWeightingRadius = 4.5;       // radius (in angstroms) of the influence zone sphere for curvature weighting

  pr.spreadReceptorSkin = false;
  pr.randomRotate = false;

  pr.staticMoleculePdb = NULL;
  pr.movingMoleculePdb = NULL;
  pr.staticMoleculeF2d = NULL;
  pr.movingMoleculeF2d = NULL;

  pr.staticMoleculePQR = NULL;  
  pr.movingMoleculePQR = NULL;

  pr.staticMoleculeSCReRaw = NULL;
  pr.staticMoleculeSCImRaw = NULL;
  pr.staticMoleculeElecReRaw = NULL;

  pr.movingMoleculeSCReRaw = NULL;
  pr.movingMoleculeSCImRaw = NULL;
  pr.movingMoleculeElecReRaw = NULL;

  pr.outputFilename = NULL;
  pr.numCentersA = 0;
  pr.numCentersB = 0;
  pr.typeA = NULL;
  pr.typeB = NULL;
  pr.atNamesA = NULL;
  pr.atNamesB = NULL;
  pr.resTypesA = NULL;
  pr.resTypesB = NULL;

  pr.resNumsA = NULL;
  pr.resNumsB = NULL;

  pr.chargesA = NULL;
  pr.hydrophobicityA = NULL;
  pr.radiiA = NULL;
  pr.chargesB = NULL;
  pr.hydrophobicityB = NULL;
  pr.radiiB = NULL;

  pr.hbondTypeA = NULL;
  pr.hbondTypeB = NULL;

  pr.rotations = NULL;
  pr.rotGraph = NULL;
  pr.pruneAngle = 0;

  // static molecule
  pr.xkAOrig = NULL;
  pr.ykAOrig = NULL;
  pr.zkAOrig = NULL;

  // moving molecule
  pr.xkBOrig = NULL;
  pr.ykBOrig = NULL;
  pr.zkBOrig = NULL;

  // RMSD calculations
  pr.nbRMSDAtoms = 0; // # of atoms used to compute RMSD
  pr.atNums = NULL;   // indices of atoms in moving molecules
  pr.xRef = NULL;  // X-coord of atoms in reference position
  pr.yRef = NULL;  // Y-coord of atoms in reference position
  pr.zRef = NULL;  // Z-coord of atoms in reference position

  pr.transformationFilename = NULL; // contains 4 x 4 transformation matrices for the moving molecule for postprocessing (e.g., vdW calculation)
  pr.vdWGridSize = 512;  // grid size for vdW potential computation
  pr.compQuadVdW = false;      // if set to true, the vdW potential is also computed using the quadratic time algorithm
  pr.surfaceBasedVdW = false;  // if set to true, only the surface atoms of the moving molecule are used for vdW computation

  pr.applyVdWFilter = true;   // if set to true, on-the-fly filtering based on vdW potential is performed
  pr.vdWCutoffLow = 0;         // when applyVdWFilter is true and #clashes < clashTolerance / 2, poses with vdW potential > vdWCutoffLow are penalized
  pr.vdWCutoffHigh = 5;        // when applyVdWFilter is true and #clashes >= clashTolerance / 2, poses with vdW potential > vdWCutoffHigh are penalized  
  pr.vdWEqmRadScale = 0.3;     // all r_eqm values are multiplied by this factor

  pr.vdWWellWidth = 0;  //  smooth vdW energy well width

  pr.vdWTolerance = 20;           // when pruneAngle is positive and applyVdWFilter is set,
                                  // a node is considered good it has a neighbor with
                                  // vdWScore <= vdWCutoffLow + vdWTolerance

  pr.applyClashFilter = true;    // if set to true, on-the-fly filtering based on number of atomic clashes is performed
  pr.eqmDistFrac = 0.5;           // two atoms clash if distance between atom centers < eqmDistFrac * r_eqm_XY
  pr.clashTolerance = 10;        // maximum number of clashes tolerated
  pr.forbiddenVolClashTolerance = 0; 
  pr.applyForbiddenVolumeFilter = false;
  pr.forbiddenVolumeFileType = -1; 
  pr.clashWeight = -0.5;         // weight given to each clash when added to total score

  pr.applyMiscFilter = false;     // when set to 'true' various minor filters are applied

  pr.applyPseudoGsolFilter = true;    // if set to true, on-the-fly filtering based on pseudo solvation energy is performed
  pr.pseudoGsolCutoff = 0;            // when applyPseudoGsolFilter is set to true, all solutions with pseudo solvation energy above pseudoGsolCutoff are discarded
  pr.pseudoGsolWeight = 0.0;          // after cutoff surviving docking poses have their scores increased (additive) by weighted pseudoGsol value
  pr.pseudoGsolFilterLowestRank = 1500;

  pr.applyDispersionFilter = false;   // if set to true, on-the-fly filtering based on solute-solvent dispersion energy is performed
  pr.dispersionCutoff = 0.0;            // when applyDispersionFilter is set to true, all solutions with dispersion energy change below dispersionCutoff are discarded
  pr.dispersionWeight = 0.0;          // after cutoff surviving docking poses have their scores increased (additive) by weighted disperions energy change
  pr.dispersionEnergyLimit = 1000.0;  // dispersion energy will be trancated to remain within [ dispersionEnergyLimit, dispersionEnergyLimit ]
  pr.dispersionMinAtomRadius = 0.1;   // smallest atom radius for dispersion energy calculation

  pr.filterScaleDownFactor = 0.2;     // if a solution is filtered scale down its score by this factor instead of discarding it

  pr.filterDepth = 4;           // for a given rotation, all configurations that are ranked beyond filterDepth
                                 // when sorted by decreasing order by score are automatically filtered
                                 // ( value -1 means no such automatic filtering will be done )

  pr.peaksPerRotation = 4;  // at most how many solutions per rotation should be retained
                            // initially set to 0; will be set to (pr.numFreq)^3 if the user does not set it
  pr.clusterTransRad = 1.2;   // radius (in angstroms) of a cluster in translational space;
                            // if a solution is kept, no solution within clusterTransRad of it can be retained
  pr.clusterTransSize = 1;  // maximum number of solutions retained within a cluster in translational space
  pr.clusterRotRad = 0;     // radius (in degrees) of a cluster in rotational space;
                            // if a solution is kept, no solution within clusterTransRad and clusterRotRad of it can be retained

  if ( pr.complexType == 'A' )
    {
      pr.skinSkinWeight = 0.73;
      pr.skinCoreWeight = -0.31;
      pr.coreCoreWeight = 31.0;

      pr.simpleChargeWeight = 0.1;
      
      pr.clashTolerance = 2;
      pr.clashWeight = -30; 
      
      pr.hydroMinRatio = 1.5;  
      pr.hydroRatioTolerance = 8.0;      
      pr.hydroRatioNumeratorLow = 1.25;
      pr.hydroRatioDenominatorLow = 0.45;
      pr.hydroRatioDenominatorHigh = 2.5;              

      pr.vdWCutoffHigh = 0;
      
      pr.filterDepth = 3;
      pr.peaksPerRotation = 3;      
    }
  else if ( pr.complexType == 'E' )
    {
      pr.curvatureWeightingRadius = 6.0;    

      pr.skinSkinWeight = 0.78;
      pr.skinCoreWeight = -0.08;
      pr.coreCoreWeight = 5.0;

      pr.elecScale = 0.15;
      pr.elecKernelVoidRad = 3.0;
      pr.elecKernelDistLow = 6.0;
      pr.elecKernelDistHigh = 8.0;
      pr.elecKernelValLow = 1.0;
      pr.elecKernelValHigh = 80.0;

      pr.simpleChargeWeight = 5.5;

      pr.clashTolerance = 9;

      pr.hydrophobicityWeight = 9.0;
      pr.hydroMinRatio = 1.22;  
      pr.hydroRatioTolerance = 8.0;      
      pr.hydroRatioNumeratorLow = 2.0;
      pr.hydroRatioDenominatorLow = 1.0;
      pr.hydroRatioDenominatorHigh = 7.0; 
      
      pr.vdWCutoffHigh = 20;
      
      pr.filterDepth = 2;
      pr.peaksPerRotation = 2;                         
    }

  pr.rerank = false;
  pr.applyAntibodyFilter = true;     
  pr.applyEnzymeFilter = true;
  pr.applyResidueContactFilter = true;
  
#ifdef LIBMOL_FOUND
  pr.applyHbondFilter = false;
#endif

  pr.spectrum = NULL;
  pr.numBands = 0;
  pr.bands = NULL;

  pr.numRerank = 2000;      // number of top positions to rerank
  pr.rerankerPseudoGsolWeight = 1.0;
  pr.rerankerDispersionWeightLow = -1.1;  
  pr.rerankerDispersionWeightHigh = -21.0;
  pr.rerankerF2DockScoreWeight = 100.0;    


  pr.control = 9876;

  //printf("AAAAAAAAAAApr = %p %d\n", &pr, pr.numThreads);

  // overwrite defaults in pr with values from the parameter file
  if (! setParamFromFile( &pr, paramFileName) ) return 1;

//  pr.coreCoreWeight *= ( pr.numCentersB / 3000.0 );

#ifdef WITH_ALL_RMSD
  pr.numThreads = 1;
#endif
  //printf("BNUMTHREAD = %p, %d\n", &pr, pr.numThreads);

  // do the docking
  if ( argc == 2 ) dock( &pr );
  else
    {
     if ( strcasecmp( argv[ 1 ], "-score" ) == 0 ) scoreUntransformed( &pr );
          else if ( strcasecmp( argv[ 1 ], "-savegrid" ) == 0 ) saveGrid( &pr );
               else if ( strcasecmp( argv[ 1 ], "-effGridFile" ) == 0 ) computeEffGrid( pr.minEffGridSize, pr.maxEffGridSize );               
    }
}
