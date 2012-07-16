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

#include "utils.h"


void printError( char *format, ... )
{
   char eMsg[ 500 ];
   va_list args;
   
   va_start( args, format );
   
   vsprintf( eMsg, format, args );
   
   va_end( args );
   
   printf( "\nError: %s\n\n", eMsg );   
   
   fflush( stdout );   
}


void flushPrint( char *format, ... )
{
   va_list args;
   
   va_start( args, format );   
   vprintf( format, args );
   va_end( args );

   fflush( stdout );
}


void f_printf( FILE *fp, char *format, ... )
{
   va_list args, args2;
   
   va_start( args, format );
   va_copy( args2, args );
   
   vprintf( format, args );
   vfprintf( fp, format, args2 );
   
   va_end( args );
}


double getTime( void )
{
#ifdef _WIN32
   time_t ltime;
   _timeb tstruct;
   time( &ltime );
   _ftime( &tstruct );
   return ( double ) ( ltime + 1e-3 * ( tstruct.millitm ) );
#else
   struct timeval t;
   gettimeofday( &t, NULL );
   return ( double )( t.tv_sec + 1e-6 * t.tv_usec );
#endif
}


int skipWhiteSpaces( char *buf, int i )
{
   int j = i;
   
   while ( buf[ j ] && isspace( buf[ j ] ) ) j++;
   
   return j;
}


int skipInitial( char *s1, char *s2, char *p )
{
   int i, j, l1, l2;
   
   l1 = 0;
   while ( s1[ l1 ] ) l1++;

   l2 = 0;
   while ( s2[ l2 ] ) l2++;

   if ( l2 > l1 ) return 0;
   
   for ( i = 0; i < l2; i++ )
     if ( s1[ i ] != s2[ i ] ) return 0;
     
   i = 0;
   j = l2;
   
   while ( s1[ j ] ) p[ i++ ] = s1[ j++ ];
   
   p[ i ] = 0;  
    
   return 1;
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


bool getDoublesInRange( char *buf, int i1, int i2, double *v )
{ 
   int i = 1, j = 0, k = 0;
   double u;
   
   while ( i < i1 )
     {
       j = getDouble( buf, j, &u );
       if ( j == k ) return false;
       i++;
       k = j;
     } 
     
   for ( int l = 0; l < i2 - i1 + 1; l++ )   
     {
       j = getDouble( buf, j, v + l );
       if ( j == k ) return false;
       k = j;     
     }
   
   return true;
}


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
