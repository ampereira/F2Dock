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

#include "miscIdent.h"

//
// arand, 8-23,2011
//        modifying this file to handle "more standard" pqr format
//        since pqr files with chains are causing this to break
//
//        If this breaks other example, we need to work harder to 
//        decide how to parse the input.
//

void getString1(char* line, char* str, int pos, int len)
{
	char buffer[80];
	//	assert(len < 80);
	strcpy(str, "");
	strncpy(buffer, &line[pos], len);
	buffer[len] = '\0';
	int i;
	for(i=0; i<len; i++)
	{
		str[i] = buffer[i];
	}
	str[len] = '\0';
}

void getChar1(char* line, char* str, int pos)
{
	char buffer[80];
	*str = ' ';
	strncpy(buffer, &line[pos], 1);
	buffer[1] = '\0';
	sscanf(buffer, "%c", str);
}

void getInt1(char* line, int* i, int pos, int len)
{
	char buffer[80];
	*i = -1;
	strncpy(buffer, &line[pos], len);
	buffer[len] = '\0';
	sscanf(buffer, "%d", i);
}

void getDouble1(char* line, double* d, int pos, int len)
{
	char buffer[80];
	*d = -1;
	strncpy(buffer, &line[pos], len);
	buffer[len] = '\0';
	sscanf(buffer, "%lf", d);
}

void getFloat1(char* line, float* d, int pos, int len)
{
	char buffer[80];
	*d = -1;
	strncpy(buffer, &line[pos], len);
	buffer[len] = '\0';
	sscanf(buffer, "%f", d);
}




int getResidueID( char *resName )
{
  int r = NONE;

  switch ( resName[ 0 ] )
    {
      case 'A' : switch ( resName[ 1 ] )
                   {
                     case 'L' : r = ALA;
                                break;

                     case 'R' : r = ARG;
                                break;

                     case 'S' : if ( resName[ 2 ] == 'N' ) r = ASN;
                                else r = ASP;
                                break;
                   }

                 break;

      case 'C' : r = CYS;
                 break;

      case 'G' : if ( resName[ 1 ] == 'L' )
                   {
                     switch ( resName[ 2 ] )
                       {
                         case 'N' : r = GLN;
                                    break;

                         case 'U' : r = GLU;
                                    break;

                         case 'Y' : r = GLY;
                                    break;
                       }
                   }

                 break;

      case 'H' : r = HIS;
                 break;

      case 'I' : r = ILE;
                 break;

      case 'L' : if ( resName[ 1 ] == 'E' ) r = LEU;
                 else r = LYS;
                 break;

      case 'M' : r = MET;
                 break;

      case 'P' : if ( resName[ 1 ] == 'H' ) r = PHE;
                 else r = PRO;
                 break;

      case 'S' : r = SER;
                 break;

      case 'T' : switch ( resName[ 1 ] )
                   {
                     case 'H' : r = THR;
                                break;

                     case 'R' : r = TRP;
                                break;

                     case 'Y' : r = TYR;
                                break;
                   }
                   
                 break;

      case 'V' : r = VAL;
                 break;
    }
              
   return r;              
}


bool countResidues( char *pqrFile, int res, int *count, int *total )
{
   FILE *fp;

   fp = fopen( pqrFile, "rt" );

   if ( fp == NULL )
     {
        printError( (char *)"Failed to open PQR file (%s)!", pqrFile );
        return false;
     }

   *count = *total = 0;

   char line[ 101 ];
   int l = -1;

   while ( fgets( line, 100, fp ) != NULL )
     {
        int i = 0, j;
        char tmp[ 100 ];

        getString1( line, tmp, 0, 6 );   // get 'ATOM'/'HETATM', and ignore

        if ( tmp[ 0 ] != 'A' ) continue;

        //getInt1( line, &j,6, 5 );            // get atom number, and ignore
        //getString1( line, tmp ,12, 4);        // get atom name, and ignore

        getString1( line, tmp, 17, 3 );   // get residue name

        getInt1( line, &j, 22, 4 );            // get residue number
        
        if ( j != l ) 
          {
            if ( getResidueID( tmp ) == res ) ( *count )++;
            ( *total )++;
          }  
        
        l = j;
     }

   fclose( fp );

   return true;
}


bool isGXY( int nRes, RESIDUE *res, int i )
{
   if ( i + 2 >= nRes ) return false;
   
   int R1 = res[ i ].resID;
   int R2 = res[ i + 1 ].resID;
   int R3 = res[ i + 2 ].resID;        

   if ( R1 != GLY ) return false;
   
   if ( ( R2 != VAL ) && ( R2 != LEU ) && ( R2 != ILE ) && ( R2 != ALA ) && ( R2 != SER ) && ( R2 != THR ) && ( R2 != ASP ) && ( R2 != PRO ) ) return false;
   
   if ( ( R3 != VAL ) && ( R3 != LEU ) && ( R3 != ILE ) && ( R3 != ALA ) && ( R3 != SER ) && ( R3 != THR ) && ( R3 != ASP ) && ( R3 != PRO ) ) return false;   
   
   return true;
}


bool isYXG( int nRes, RESIDUE *res, int i )
{
   if ( i + 2 >= nRes ) return false;
   
   int R3 = res[ i ].resID;
   int R2 = res[ i + 1 ].resID;
   int R1 = res[ i + 2 ].resID;      
 
   if ( R1 != GLY ) return false;
   
   if ( ( R2 != VAL ) && ( R2 != LEU ) && ( R2 != ILE ) && ( R2 != ALA ) && ( R2 != SER ) && ( R2 != THR ) && ( R2 != ASP ) && ( R2 != PRO ) ) return false;
   
   if ( ( R3 != VAL ) && ( R3 != LEU ) && ( R3 != ILE ) && ( R3 != ALA ) && ( R3 != SER ) && ( R3 != THR ) && ( R3 != ASP ) && ( R3 != PRO ) ) return false;   
   
   return true;
}


bool readResidues( char *pqrFile, int *nRes, RESIDUE **res, int *nChn, int **chn )
{
   FILE *fp;

   fp = fopen( pqrFile, (char *)"rt" );

   if ( fp == NULL )
     {
        printError( (char *)"Failed to open PQR file (%s)!", pqrFile );
        return false;
     }

   int numRes = 0, numChains = 1;

   char line[ 101 ];
   int l = -1;

   while ( fgets( line, 100, fp ) != NULL )
     {
        int i = 0, j;
        char tmp[ 100 ];
        char tmp1[ 100 ];

        getString1( line, tmp, 0, 6 );   // get 'ATOM'/'HETATM', and ignore

        if ( tmp[ 0 ] != 'A' ) continue;

        //i = getInt( line, i, &j );            // get atom number, and ignore
        //i = getString( line, i, tmp );        // get atom name, and ignore

        getString1( line, tmp, 17, 3 );   // get residue name, and ignore	

        getInt1( line, &j, 22, 4 );            // get residue number

        if ( j != l ) numRes++;
        
        if ( j < l ) numChains++;

        l = j;
     }

   fclose( fp );

   *nRes = numRes;
   *nChn = numChains;
   ( *res ) = ( RESIDUE * ) malloc( ( *nRes ) * sizeof( RESIDUE ) );
   ( *chn ) = ( int * ) malloc( ( *nChn + 1 ) * sizeof( int ) );   

   if ( ( *res == NULL ) || ( *chn == NULL ) )
     {
        printError( (char *)"Failed to allocate memory!" );
        return false;
     }

   fp = fopen( pqrFile, (char *)"rt" );

   if ( fp == NULL )
     {
         printError( (char *)"Failed to open PQR file (%s)!", pqrFile );
         return false;
     }

   int k = 0, c = 1;
   l = -1;
   
   ( *chn )[ 0 ] = 0;

   while ( fgets( line, 100, fp ) != NULL )
     {
        int i = 0, j;
        char tmp[ 100 ];
        char tmp1[ 100 ];

        i = getAlphaString( line, i, tmp );   // get 'ATOM'/'HETATM', and ignore

        if ( tmp[ 0 ] != 'A' ) continue;

	// old
        //i = getInt( line, i, &j );            // get atom number, and ignore
        //i = getString( line, i, tmp );        // get atom name, and ignore

        getString1( line, tmp, 17, 3 );   // get residue name

        getInt1( line, &j, 22, 4 );            // get residue number

        if ( j != l )
          {
            if ( j < l ) ( *chn )[ c++ ] = k;

            ( *res )[ k ].resNum = j;
            ( *res )[ k ].chainID = c;
            ( *res )[ k ].resID = getResidueID( tmp );

            k++;
            
            l = j;
          }
     }
     
   ( *chn )[ c ] = k;  

   fclose( fp );

   return true;
}



bool readAtomsAndResidues( char *pqrFile, int *nAtm, double **atm, int *nRes, RESIDUE **res, int *nChn, int **chn )
{
   FILE *fp;

   fp = fopen( pqrFile, (char *)"rt" );

   if ( fp == NULL )
     {
        printError( (char *)"Failed to open PQR file (%s)!", pqrFile );
        return false;
     }

   int numAtoms = 0, numRes = 0, numChains = 1;

   char line[ 101 ];
   int l = -1;

   while ( fgets( line, 100, fp ) != NULL )
     {
        int i = 0, j;
        char tmp[ 100 ];
        char tmp1[ 100 ];

        i = getAlphaString( line, i, tmp );   // get 'ATOM'/'HETATM', and ignore

        if ( tmp[ 0 ] != 'A' ) continue;

        //i = getInt( line, i, &j );            // get atom number, and ignore
        //i = getString( line, i, tmp );        // get atom name, and ignore
        //i = getAlphaString( line, i, tmp );   // get residue name, and ignore
        //i = getAlphaString( line, i, tmp1 );   // get residue name, and ignore

        getInt1( line, &j, 22, 4 );            // get residue number

        numAtoms++;

        if ( j != l ) numRes++;
        
        if ( j < l ) numChains++;

        l = j;
     }

   fclose( fp );

   *nAtm = numAtoms;
   *nRes = numRes;
   *nChn = numChains;
   ( *atm ) = ( double * ) malloc( 7 * ( *nAtm ) * sizeof( double ) );
   ( *res ) = ( RESIDUE * ) malloc( ( *nRes ) * sizeof( RESIDUE ) );
   ( *chn ) = ( int * ) malloc( ( *nChn + 1 ) * sizeof( int ) );   

   if ( ( *atm == NULL ) || ( *res == NULL ) || ( *chn == NULL ) )
     {
        printError( (char *)"Failed to allocate memory!" );
        return false;
     }

   fp = fopen( pqrFile, (char *)"rt" );

   if ( fp == NULL )
     {
         printError( (char *)"Failed to open PQR file (%s)!", pqrFile );
         return false;
     }

   int k = 0, n = 0, c = 1;
   l = -1;
   
   ( *chn )[ 0 ] = 0;

   while ( fgets( line, 100, fp ) != NULL )
     {
        int i = 0, j;
        char tmp[ 100 ];
        char tmp1[ 100 ];

        i = getAlphaString( line, i, tmp );   // get 'ATOM'/'HETATM', and ignore

        if ( tmp[ 0 ] != 'A' ) continue;

        //i = getInt( line, i, &j );            // get atom number, and ignore
        //i = getString( line, i, tmp );        // get atom name, and ignore

        getString1( line, tmp, 17, 3 );   // get residue name

        getInt1( line, &j, 22, 4 );            // get residue number

        if ( j != l )
          {
            if ( j < l ) ( *chn )[ c++ ] = k;

            ( *res )[ k ].resNum = j;
            ( *res )[ k ].chainID = c;
            ( *res )[ k ].resID = getResidueID( tmp );

            k++;
            
            l = j;
          }
          
        double v;  
          
        getDouble1( line, &v, 30, 8 );         // get X coordinate
        ( *atm )[ 7 * n + 0 ] = v;    

        getDouble1( line, &v, 38, 8 );         // get Y coordinate
        ( *atm )[ 7 * n + 1 ] = v;    

        getDouble1( line, &v, 46, 8 );         // get Z coordinate
        ( *atm )[ 7 * n + 2 ] = v;    
    
        getDouble1( line, &v, 54, 9 );         // get charge
        ( *atm )[ 7 * n + 3 ] = v;    

        getDouble1( line, &v, 63, 7 );         // get radius
        ( *atm )[ 7 * n + 4 ] = v;   

        ( *atm )[ 7 * n + 5 ] = j;            // residue number  
        
        ( *atm )[ 7 * n + 6 ] = c;            // chain number  
                          
        n++;  
     }
     
   ( *chn )[ c ] = k;  

   fclose( fp );

   return true;
}



bool readAtomsOnly( char *pqrFile, int *nAtm, double **atm )
{
   FILE *fp;

   fp = fopen( pqrFile, (char *)"rt" );

   if ( fp == NULL )
     {
        printError( (char *)"Failed to open PQR file (%s)!", pqrFile );
        return false;
     }

   int numAtoms = 0;

   char line[ 101 ];
   
   while ( fgets( line, 100, fp ) != NULL )
     {
        int i = 0, j;
        char tmp[ 100 ];

        i = getAlphaString( line, i, tmp );   // get 'ATOM'/'HETATM', and ignore

        if ( tmp[ 0 ] != 'A' ) continue;

        numAtoms++;
     }

   fclose( fp );

   *nAtm = numAtoms;
   ( *atm ) = ( double * ) malloc( 5 * ( *nAtm ) * sizeof( double ) );

   if ( *atm == NULL )
     {
        printError( (char *)"Failed to allocate memory!" );
        return false;
     }

   fp = fopen( pqrFile, (char *)"rt" );

   if ( fp == NULL )
     {
         printError( (char *)"Failed to open PQR file (%s)!", pqrFile );
         return false;
     }

   int n = 0;

   while ( fgets( line, 100, fp ) != NULL )
     {
        int i = 0, j;
        char tmp[ 100 ];

        i = getAlphaString( line, i, tmp );   // get 'ATOM'/'HETATM', and ignore

        if ( tmp[ 0 ] != 'A' ) continue;

        //i = getInt( line, i, &j );            // get atom number, and ignore
        //i = getString( line, i, tmp );        // get atom name, and ignore

        getString1( line, tmp, 17, 3 );   // get residue name

        getInt1( line, &j, 22, 4 );            // get residue number
          
        double v;  
          
        getDouble1( line, &v, 30, 8);         // get X coordinate
        ( *atm )[ 5 * n + 0 ] = v;    

        getDouble1( line, &v, 38, 8 );         // get Y coordinate
        ( *atm )[ 5 * n + 1 ] = v;    

        getDouble1( line, &v, 46, 8 );         // get Z coordinate
        ( *atm )[ 5 * n + 2 ] = v;    
    
        getDouble1( line, &v, 54, 9 );         // get charge
        ( *atm )[ 5 * n + 3 ] = v;    

        getDouble1( line, &v, 63, 7 );         // get radius
        ( *atm )[ 5 * n + 4 ] = v;   
        
        n++;  
     }
     
   fclose( fp );

   return true;
}


bool readAtomsWithResidueInfo( char *pqrFile, int *nAtm, double **atm )
{
   FILE *fp;

   fp = fopen( pqrFile, (char *)"rt" );

   if ( fp == NULL )
     {
        printError( (char *)"Failed to open PQR file (%s)!", pqrFile );
        return false;
     }

   int numRes = 0;

   char line[ 101 ];
   int l = -1, r = -1;

   while ( fgets( line, 100, fp ) != NULL )
     {
        int i = 0, j;
        char tmp[ 100 ];

        i = getAlphaString( line, i, tmp );   // get 'ATOM'/'HETATM', and ignore

        if ( tmp[ 0 ] != 'A' ) continue;

        //i = getInt( line, i, &j );            // get atom number, and ignore
        //i = getString( line, i, tmp );        // get atom name, and ignore

        getString1( line, tmp, 17, 3 );   // get residue name
        
        int resID = getResidueID( tmp );

        getInt1( line, &j, 22, 4 );            // get residue number

        if ( ( j != l ) || ( resID != r ) ) numRes++;

        l = j;
        r = resID;
     }

   fclose( fp );
   
   printf( (char *)"\nnumRes = %d\n", numRes );
   fflush( stdout );

   *nAtm = numRes;
   ( *atm ) = ( double * ) malloc( 5 * ( *nAtm ) * sizeof( double ) );

   if ( *atm == NULL )
     {
        printError( (char *)"Failed to allocate memory!" );
        return false;
     }

   fp = fopen( pqrFile, (char *)"rt" );

   if ( fp == NULL )
     {
         printError( (char *)"Failed to open PQR file (%s)!", pqrFile );
         return false;
     }

   int n = 0;
   
   l = -1; r = -1;
   
   bool done = false;

   while ( fgets( line, 100, fp ) != NULL )
     {
        int i = 0, j;
        char tmp[ 100 ], atn[ 100 ];

        i = getAlphaString( line, i, tmp );   // get 'ATOM'/'HETATM', and ignore

        if ( tmp[ 0 ] != 'A' ) continue;

        //i = getInt( line, i, &j );            // get atom number, and ignore

        getString1( line, atn, 12, 4 );        // get atom name

        getString1( line, tmp, 17, 3 );   // get residue name
        
        int resID = getResidueID( tmp );
       
        getInt1( line, &j, 22, 4 );            // get residue number       
	
        if ( ( j != l ) || ( resID != r ) ) done = false;

        l = j;
        r = resID;        
                
	// arand: changed "CA" to " CA " and "CB" to " CB " to match the new method of reading data...
        if ( !done && ( resID != NONE ) && ( ( ( resID == GLY ) && !strcmp( atn, (char *)" CA " ) ) || ( ( resID != GLY ) && !strcmp( atn, (char *)" CB " ) ) ) )
          {              
            double v;  
              
            getDouble1( line, &v, 30, 8 );         // get X coordinate
            ( *atm )[ 5 * n + 0 ] = v;    
    
	    getDouble1( line, &v, 38, 8 );         // get Y coordinate
            ( *atm )[ 5 * n + 1 ] = v;    
    
            getDouble1( line, &v, 46, 8 );         // get Z coordinate
            ( *atm )[ 5 * n + 2 ] = v;    
        
            ( *atm )[ 5 * n + 3 ] = resID - 1;    
            
            ( *atm )[ 5 * n + 4 ] = j;   
            
            n++;
            
            done = true;
          }  
     }
   
   *nAtm = n;
     
   fclose( fp );

   return true;
}


bool readGlycines( char *pqrFile, int *nAtm, double **atm )
{

   int nRes = 0, nChn = 0;
   RESIDUE *res = NULL;
   int *chn = NULL;
   
   if ( !readResidues( pqrFile, &nRes, &res, &nChn, &chn ) ) 
     {
       freeMem( res );
       freeMem( chn);
       return false;
     }  

   FILE *fp;

   fp = fopen( pqrFile, (char *)"rt" );

   if ( fp == NULL )
     {
        printError( (char *)"Failed to open PQR file (%s)!", pqrFile );
        return false;
     }

   int numAtoms = 0;

   char line[ 101 ];
   int k = -1, l = -1;

   while ( fgets( line, 100, fp ) != NULL )
     {
        int i = 0, j;
        char tmp[ 100 ];

        i = getAlphaString( line, i, tmp );   // get 'ATOM'/'HETATM', and ignore

        if ( tmp[ 0 ] != 'A' ) continue;

        numAtoms++;
     }

   fclose( fp );

   *nAtm = numAtoms;
   ( *atm ) = ( double * ) malloc( 5 * ( *nAtm ) * sizeof( double ) );

   if ( *atm == NULL )
     {
        printError( (char *)"Failed to allocate memory!" );
        return false;
     }

   fp = fopen( pqrFile, (char *)"rt" );

   if ( fp == NULL )
     {
         printError( (char *)"Failed to open PQR file (%s)!", pqrFile );
         return false;
     }

   int n = 0;
   k = -1; l = -1;   

   while ( fgets( line, 100, fp ) != NULL )
     {
        int i = 0, j;
        char tmp[ 100 ];

        i = getAlphaString( line, i, tmp );   // get 'ATOM'/'HETATM', and ignore

        if ( tmp[ 0 ] != 'A' ) continue;

        //i = getInt( line, i, &j );            // get atom number, and ignore
        //i = getString( line, i, tmp );        // get atom name, and ignore

        getString1( line, tmp, 17, 3 );   // get residue name

        getInt1( line, &j, 22, 4 );            // get residue number
        
        if ( j != l )
          {
            k++;
            l = j;
          }        
        
        if ( isGXY( nRes, res, k ) || isYXG( nRes, res, k ) 
          || ( ( k >= 1 ) && ( isGXY( nRes, res, k - 1 ) || isYXG( nRes, res, k - 1 ) ) ) 
          || ( ( k >= 2 ) && ( isGXY( nRes, res, k - 2 ) || isYXG( nRes, res, k - 2 ) ) ) )
          {
              double v;  
                
              getDouble1( line, &v, 30, 8 );         // get X coordinate
              ( *atm )[ 5 * n + 0 ] = v;    
        
              getDouble1( line, &v, 38, 8);         // get Y coordinate
              ( *atm )[ 5 * n + 1 ] = v;    
        
              getDouble1( line, &v, 46, 8 );         // get Z coordinate
              ( *atm )[ 5 * n + 2 ] = v;    
          
              getDouble1( line, &v, 54, 9 );         // get charge
              ( *atm )[ 5 * n + 3 ] = v;    
        
              getDouble1( line, &v, 63, 7 );         // get radius
              ( *atm )[ 5 * n + 4 ] = v;   
              
              n++;            
          }          
     }
     
   *nAtm = n;  
     
   fclose( fp );

   freeMem( res );
   freeMem( chn);

   return true;
}


bool getTotalCharge( char *pqrFile, double *tCharge )
{
   FILE *fp;

   fp = fopen( pqrFile, (char *)"rt" );

   if ( fp == NULL )
     {
        printError( (char *)"Failed to open PQR file (%s)!", pqrFile );
        return false;
     }
     
   *tCharge = 0;  
   char line[ 101 ];

   while ( fgets( line, 100, fp ) != NULL )
     {
        int i = 0, j;
        char tmp[ 101 ];

        i = getAlphaString( line, i, tmp );   // get 'ATOM'/'HETATM', and ignore

        if ( tmp[ 0 ] != 'A' ) continue;

        //i = getInt( line, i, &j );            // get atom number, and ignore
        //i = getString( line, i, tmp );        // get atom name, and ignore

	// don't use these either?
        getString1( line, tmp, 17, 3 );   // get residue name
	getInt1( line, &j, 22, 4 );            // get residue number
          
        double v;  
          
        //i = getDouble( line, i, &v );         // get X coordinate
        //i = getDouble( line, i, &v );         // get Y coordinate
        //i = getDouble( line, i, &v );         // get Z coordinate
    
        getDouble1( line, &v, 54, 9 );         // get charge
        
        *tCharge += v;

        //i = getDouble( line, i, &v );         // get radius
     }
     
   fclose( fp );

   return true;
}


bool CDR_L1_start( RESIDUE *res, int i )
{
   if ( res[ i - 1 ].resID == CYS ) return true;
       
   return false;    
}


bool CDR_L1_end( RESIDUE *res, int i, int chainEnd )
{
   if ( i + 3 > chainEnd ) return false;
   
   if ( ( res[ i + 1 ].resID == TRP ) 
     && ( ( ( res[ i + 2 ].resID == TYR ) && ( res[ i + 3 ].resID == GLN ) ) 
       || ( ( res[ i + 2 ].resID == LEU ) && ( res[ i + 3 ].resID == GLN ) )
       || ( ( res[ i + 2 ].resID == PHE ) && ( res[ i + 3 ].resID == GLN ) ) 
       || ( ( res[ i + 2 ].resID == TYR ) && ( res[ i + 3 ].resID == LEU ) ) ) ) return true;
       
   return false;    
}


bool CDR_L2_start( RESIDUE *res, int i )
{
   if ( ( ( res[ i - 2 ].resID == LEU ) && ( res[ i - 1 ].resID == TYR ) )   // LEU-TYR: I added; not in the list
     || ( ( res[ i - 2 ].resID == ILE ) && ( res[ i - 1 ].resID == TYR ) ) 
     || ( ( res[ i - 2 ].resID == VAL ) && ( res[ i - 1 ].resID == TYR ) )
     || ( ( res[ i - 2 ].resID == ILE ) && ( res[ i - 1 ].resID == LYS ) ) 
     || ( ( res[ i - 2 ].resID == ILE ) && ( res[ i - 1 ].resID == PHE ) ) ) return true;
       
   return false;    
}


bool CDR_L3_start( RESIDUE *res, int i )
{
   if ( res[ i - 1 ].resID == CYS ) return true;
       
   return false;    
}


bool CDR_L3_end( RESIDUE *res, int i, int chainEnd )
{
   if ( i + 4 > chainEnd ) return false;
   
   if ( ( res[ i + 1 ].resID == PHE ) && ( res[ i + 2 ].resID == GLY ) && ( res[ i + 4 ].resID == GLY ) ) return true;
       
   return false;    
}


bool CDR_H1_start( RESIDUE *res, int i )
{
   if ( res[ i - 4 ].resID == CYS ) return true;
       
   return false;    
}


bool CDR_H1_end( RESIDUE *res, int i, int chainEnd )
{
   if ( i + 2 > chainEnd ) return false;
   
   if ( ( res[ i + 1 ].resID == TRP ) 
     && ( ( res[ i + 2 ].resID == VAL ) || ( res[ i + 2 ].resID == ILE ) || ( res[ i + 2 ].resID == ALA ) ) ) return true;
       
   return false;    
}

bool CDR_H2_end( RESIDUE *res, int i, int chainEnd )
{
   if ( i + 3 > chainEnd ) return false;
   
   if ( ( ( res[ i + 1 ].resID == LYS ) || ( res[ i + 1 ].resID == ARG ) || ( res[ i + 1 ].resID == SER ) )  // SER: I added; not in the list
     && ( ( res[ i + 2 ].resID == LEU ) ||  ( res[ i + 2 ].resID == ILE ) 
       || ( res[ i + 2 ].resID == VAL ) ||  ( res[ i + 2 ].resID == PHE ) 
       || ( res[ i + 2 ].resID == THR ) ||  ( res[ i + 2 ].resID == ALA ) )
     && ( ( res[ i + 3 ].resID == THR ) ||  ( res[ i + 3 ].resID == SER ) 
       || ( res[ i + 3 ].resID == ILE ) ||  ( res[ i + 3 ].resID == ALA ) ) ) return true;
       
   return false;    
}


bool CDR_H3_start( RESIDUE *res, int i )
{
   if ( res[ i - 3 ].resID == CYS ) return true;
       
   return false;    
}


bool CDR_H3_end( RESIDUE *res, int i, int chainEnd )
{
   if ( i + 4 > chainEnd ) return false;
   
   if ( ( res[ i + 1 ].resID == TRP ) && ( res[ i + 2 ].resID == GLY ) && ( res[ i + 4 ].resID == GLY ) ) return true;
       
   return false;    
}


bool CDR_L3_identified( RESIDUE *res, int c, int l, int s, int *id )
{
   int i = c;
   
   while ( i < c + l )
     {
       if ( res[ i ].resNum == s ) break;
       i++;
     }
   
   if ( ( i == c + l ) || !CDR_L3_start( res, i ) ) return false;

   printf( (char *)"\nL3 = < %d, ", s );
   
   id[ 4 ] = s;

   for ( int e = s + 6; e <= s + 10; e++ )
     {
       int j = i + 1;
       
       while ( j < c + l )
         {
           if ( res[ j ].resNum == e ) break;
           j++;
         }
       
       if ( ( j < c + l ) && CDR_L3_end( res, j, c + l - 1 ) ) 
         {
           printf( (char *)"%d >\n", e );
           id[ 5 ] = e;
           return true;
         }  
     }

   return false;  
}


bool CDR_L2_L3_identified( RESIDUE *res, int c, int l, int s, int *id )
{
   int i = c;
   
   while ( i < c + l )
     {
       if ( res[ i ].resNum == s ) break;
       i++;
     }
   
   if ( ( i == c + l ) || !CDR_L2_start( res, i ) ) return false;

   id[ 2 ] = s;

   int e = s + 6;

   int j = i + 1;
   
   while ( j < c + l )
     {
       if ( res[ j ].resNum == e ) break;
       j++;
     }
   
   if ( j == c + l ) return false;
   
   printf( (char *)"\nL2 = < %d, %d >\n", s, e );
   
   id[ 3 ] = e;
   
   for ( int k = 30; k <= 36; k++ )
     if ( CDR_L3_identified( res, j, l - ( j - c ), e + k, id ) ) return true;

   return false;  
}


bool CDR_L_identified( RESIDUE *res, int c, int l, int *id )
{
   int L1Start[ ] = { 24, 23, 25, 22, 26, 21, 27 };
   
   for ( int s = 0; s < sizeof( L1Start ) / sizeof( L1Start[ 0 ] ); s++ )
     {
       int i = c;
       
       while ( i < c + l )
         {
           if ( res[ i ].resNum == L1Start[ s ] ) break;
           i++;
         }  
       
       if ( ( i == c + l ) || !CDR_L1_start( res, i ) ) continue;
       
       id[ 0 ] = L1Start[ s ];
       
       for ( int e = L1Start[ s ] + 9; e <= L1Start[ s ] + 16; e++ )
         {
           int j = i + 1;
           
           while ( j < c + l )
             {
               if ( res[ j ].resNum == e ) break;
               j++;
             }
           
           if ( ( j == c + l ) || !CDR_L1_end( res, j, c + l - 1 ) ) continue;
           
           id[ 1 ] = e;
           
           printf( (char *)"\nL1 = < %d, %d >\n", L1Start[ s ], e );
           
           if ( CDR_L2_L3_identified( res, j, l - ( j - c ), e + 16, id ) ) return true;    
         }
     }
     
   return false;  
}




bool CDR_H3_identified( RESIDUE *res, int c, int l, int s, int *id )
{
   int i = c;
   
   while ( i < c + l )
     {
       if ( res[ i ].resNum == s ) break;
       i++;
     }
   
   if ( ( i == c + l ) || !CDR_H3_start( res, i ) ) return false;

   id[ 4 ] = s;

   for ( int e = s + 2; e <= s + 24; e++ )
     {
       int j = i + 1;
       
       while ( j < c + l )
         {
           if ( res[ j ].resNum == e ) break;
           j++;
         }
       
       if ( ( j < c + l ) && CDR_H3_end( res, j, c + l - 1 ) ) 
         {
           id[ 5 ] = e;
           
           printf( (char *)"\nH3 = < %d, %d >\n", s, e );
           
           return true;
         }  
     }

   return false;  
}


bool CDR_H2_H3_identified( RESIDUE *res, int c, int l, int s, int *id )
{
   int i = c;
   
   while ( i < c + l )
     {
       if ( res[ i ].resNum == s ) break;
       i++;
     }
   
   if ( i == c + l ) return false;

   id[ 2 ] = s;

   for ( int e = s + 15; e <= s + 18; e++ )
     {
       int j = i + 1;
       
       while ( j < c + l )
         {
           if ( res[ j ].resNum == e ) break;
           j++;
         }
       
       if ( ( j == c + l ) || !CDR_H2_end( res, j, c + l - 1 ) ) continue;
       
       id[ 3 ] = e;
       
       printf( (char *)"\nH2 = < %d, %d >\n", s, e );

       for ( int k = 30; k <= 36; k++ )       
         if ( CDR_H3_identified( res, j, l - ( j - c ), e + k, id ) ) return true;
     }

   return false;  
}


bool CDR_H_identified( RESIDUE *res, int c, int l, int *id )
{
   int H1Start[ ] = { 26, 25, 27, 24, 28, 23, 29 };
   
   for ( int s = 0; s < sizeof( H1Start ) / sizeof( H1Start[ 0 ] ); s++ )
     {
       int i = c;
       
       while ( i < c + l )
         {
           if ( res[ i ].resNum == H1Start[ s ] ) break;
           i++;
         }  
       
       if ( ( i == c + l ) || !CDR_H1_start( res, i ) ) continue;
       
       id[ 0 ] = H1Start[ s ];
       
       for ( int e = H1Start[ s ] + 9; e <= H1Start[ s ] + 11; e++ )
         {
           int j = i + 1;
           
           while ( j < c + l )
             {
               if ( res[ j ].resNum == e ) break;
               j++;
             }
           
           if ( ( j == c + l ) || !CDR_H1_end( res, j, c + l - 1 ) ) continue;
           
           id[ 1 ] = e;
           
           printf( (char *)"\nH1 = < %d, %d >\n", H1Start[ s ], e );
           
           if ( CDR_H2_H3_identified( res, j, l - ( j - c ), e + 15, id ) ) return true;    
         }
     }
     
   return false;  
}


bool isAntibody( char *pqrFile )
{
   bool abody = false;
   
   RESIDUE *res = NULL;
   int *chn = NULL;
   int nRes = 0, nChn = 0;
   int L[ 6 ], H[ 6 ];
   
   if ( readResidues( pqrFile, &nRes, &res, &nChn, &chn ) )
     {
       printf( (char *)"\nnRes = %d, nChn = %d\n", nRes, nChn );
       
       for ( int i = 0; i < nChn; i++ )
        {
         printf( (char *)"\nchn[ i ] = %d, chn[ i + 1 ] - chn[ i ] = %d\n", chn[ i ], chn[ i + 1 ] - chn[ i ] );
          
         if ( CDR_L_identified( res, chn[ i ], chn[ i + 1 ] - chn[ i ], L ) || CDR_H_identified( res, chn[ i ], chn[ i + 1 ] - chn[ i ], H ) )
            {
//              int j = 0;
//              
//              while ( j < nChn )
//                {
//                  if ( ( j != i ) && CDR_H_identified( res, chn[ j ], chn[ j + 1 ] - chn[ j ] ) )
//                    {
                      abody = true;
                      break;
//                    }
//                    
//                  j++;  
//                }    
//                  
//              if ( j < nChn ) break;    
            }
        }    
     }
   
   freeMem( res );
   freeMem( chn );
   
   return abody;
}


bool getAntibodyBindingSite( char *pqrFile, int *numAtoms, double **atm, bool low, bool high, int l1, int l2 )
{
   bool gotSites = true;
   
   *atm = NULL;
   RESIDUE *res = NULL;
   int *chn = NULL;
   int nAtm = 0, nRes = 0, nChn = 0;
   int I[ 6 ];
   
   int k = 0;
   
   if ( readAtomsAndResidues( pqrFile, &nAtm, atm, &nRes, &res, &nChn, &chn ) )
     {
       int n = 0;
       
       for ( int i = 0; i < nChn; i++ )
        {
          if ( ( low && CDR_L_identified( res, chn[ i ], chn[ i + 1 ] - chn[ i ], I ) )
            || ( high && CDR_H_identified( res, chn[ i ], chn[ i + 1 ] - chn[ i ], I ) ) )
            {
              for ( int l = l1 - 1; l < l2; l++ )
                {
                  while ( ( *atm )[ 7 * n + 6 ] != i + 1 ) n++;
                  while ( ( *atm )[ 7 * n + 5 ] != I[ 2 * l ] ) n++;
                  
                  while ( ( n < nAtm ) && ( ( *atm )[ 7 * n + 6 ] == i + 1 ) && ( ( *atm )[ 7 * n + 5 ] <= I[ 2 * l + 1 ] ) )
                    {
                      for ( int j = 0; j < 5; j++ )
                         ( *atm )[ 5 * k + j ] = ( *atm )[ 7 * n + j ];
                         
                      k++;   
                      n++;
                    }
                }
            }
        }
        
       if ( low && high ) printf( (char *)"\nL3H3: %d atoms\n", n );     
       else if ( low ) printf( (char *)"\nL3: %d atoms\n", n );     
            else if ( high ) printf( (char *)"\nH3: %d atoms\n", n );     
     }
   else 
     {
       gotSites = false; 
       freeMem( atm );
     }   
   
   *numAtoms = k;
   
   freeMem( res );
   freeMem( chn );
   
   return gotSites;
}


bool countGXYandYXG( char *pqrFile, int *count, int *total, int *nGly )
{
   bool done = false;
   
   RESIDUE *res = NULL;
   int *chn = NULL;
   int nRes = 0, nChn = 0;
   
   *count = *total = *nGly = 0;
   
   if ( readResidues( pqrFile, &nRes, &res, &nChn, &chn ) )
     {
       for ( int i = 0; i < nRes; i++ )
           if ( isGXY( nRes, res, i ) || isYXG( nRes, res, i ) ) 
             {
               i += 2;
               ( *count )++;       
             }  

       for ( int i = 0; i < nRes; i++ )
           if ( res[ i ].resID == GLY ) ( *nGly )++;
         
       done = true;  
     }
     
   *total = nRes;  
   
   freeMem( res );
   freeMem( chn );
   
   return done;
}



bool initAntibodyClashFilter( char *staticPQR, char *movingPQR, bool low, bool high, int l1, int l2, clashFilter **cFilter )
{
  int numStaticAtoms = 0, numMovingAtoms = 0;

  double *staticAtoms = NULL;
  double *movingAtoms = NULL;

  bool staticAntibody = isAntibody( staticPQR );
  bool movingAntibody = isAntibody( movingPQR ); 
  
  if ( staticAntibody )
    {
      printf( (char *)"\n# static: antibody ( %s )\n", staticPQR );
      
      if ( !getAntibodyBindingSite( staticPQR, &numStaticAtoms, &staticAtoms, low, high, l1, l2 ) ) return false;
    }  
  else
    {
      printf( (char *)"\n# static: not antibody ( %s )\n", staticPQR ); 
         
      if ( !readAtomsOnly( staticPQR, &numStaticAtoms, &staticAtoms ) ) return false;
    }  
      
  if ( movingAntibody )
    {
      printf( (char *)"\n# moving: antibody ( %s )\n", movingPQR );
    
      if ( !getAntibodyBindingSite( movingPQR, &numMovingAtoms, &movingAtoms, low, high, l1, l2 ) ) 
        {
          freeMem( staticAtoms );
          return false;
        }  
    }  
  else
    {
      printf( (char *)"\n# moving: not antibody ( %s )\n", movingPQR );
    
      if ( !readAtomsOnly( movingPQR, &numMovingAtoms, &movingAtoms ) ) 
        {
          freeMem( staticAtoms );        
          return false;
        }  
        
      printf( (char *)"\n# movingAtoms = %d\n", numMovingAtoms );
      fflush( stdout );        
    }  

  if ( ( ( staticAntibody && !movingAntibody ) || ( !staticAntibody && movingAntibody ) ) && ( numStaticAtoms > 0 ) && ( numMovingAtoms > 0 ) )
    {
      ( *cFilter ) = new clashFilter( numStaticAtoms, staticAtoms, numMovingAtoms, movingAtoms, false );

      ( *cFilter )->setProximityFactors( 2.0, 1.0, 1.0 );

      printf( (char *)"\n# antibody: staticAtoms = %d, movingAtoms = %d ( filter initialized )\n", numStaticAtoms, numMovingAtoms );
    }
  else printf( (char *)"\n# antibody: staticAtoms = %d, movingAtoms = %d ( filter not initialized )\n", numStaticAtoms, numMovingAtoms );

  fflush( stdout );

  freeMem( staticAtoms );
  freeMem( movingAtoms );

  return ( ( bool ) ( ( ( staticAntibody && !movingAntibody ) || ( !staticAntibody && movingAntibody ) ) && ( numStaticAtoms > 0 ) && ( numMovingAtoms > 0 ) ) );
}

bool initEnzymeClashFilter( char *staticPQR, char *movingPQR, clashFilter **cFilter )
{
  int numStaticAtoms = 0, numMovingAtoms = 0;

  double *staticAtoms = NULL;
  double *movingAtoms = NULL;
  
  if ( !readGlycines( staticPQR, &numStaticAtoms, &staticAtoms ) ) return false;
      
  if ( !readAtomsOnly( movingPQR, &numMovingAtoms, &movingAtoms ) ) 
    {
      freeMem( staticAtoms );        
      return false;
    }  

  ( *cFilter ) = new clashFilter( numStaticAtoms, staticAtoms, numMovingAtoms, movingAtoms, false );

  ( *cFilter )->setProximityFactors( 2.0, 1.0, 1.0 );
  
  freeMem( staticAtoms );
  freeMem( movingAtoms );

  return true;
}


bool initResContFilter( char *staticPQR, char *movingPQR, char *resContFile, resContFilter **cFilter )
{
  int numStaticAtoms = 0, numMovingAtoms = 0;

  double *staticAtoms = NULL;
  double *movingAtoms = NULL;
  
  if ( !readAtomsWithResidueInfo( staticPQR, &numStaticAtoms, &staticAtoms ) ) return false;
      
  if ( !readAtomsWithResidueInfo( movingPQR, &numMovingAtoms, &movingAtoms ) ) 
    {
      freeMem( staticAtoms );        
      return false;
    }  

  printf( (char *)"\n# resCont: staticAtoms = %d, movingAtoms = %d\n", numStaticAtoms, numMovingAtoms );
  fflush( stdout );

  ( *cFilter ) = new resContFilter( numStaticAtoms, staticAtoms, numMovingAtoms, movingAtoms, resContFile, false );

  freeMem( staticAtoms );
  freeMem( movingAtoms );

  return true;
}

