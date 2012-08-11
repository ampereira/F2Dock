#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <libmol/libmol.h>

using namespace LibMol;

#ifndef DEFAULT_PRM
  #define DEFAULT_PRM "params/parm.prm"
#endif

#ifndef DEFAULT_RTF
  #define DEFAULT_RTF "params/pdbamino.rtf"
#endif

#ifndef DEFAULT_APRM
  #define DEFAULT_APRM "params/atoms.0.0.6.prm.ms.3cap+0.5ace.Hr0rec"
#endif

#ifndef MAX_ITER
  #define MAX_ITER 100
#endif

#ifndef DEFAULT_VDW_DIST_CUTOFF
  #define DEFAULT_VDW_DIST_CUTOFF 12.0
#endif

#ifndef DEFAULT_ELEC_DIST_CUTOFF
  #define DEFAULT_ELEC_DIST_CUTOFF 15.0
#endif

#ifndef DEFAULT_HBOND_DIST_CUTOFF
  #define DEFAULT_HBOND_DIST_CUTOFF ( MAX_R * 1.0 )
#endif


struct my_par
{
   int nbupdate;
   struct atomgrp *ag;
   struct agsetup *ags;
   double efacv;
   double eface;
   int useVdw, useElec, useHbond;
   double dcVdw, dcElec, dcHbond;
   double aprxdcVdw, aprxdcElec;
   OCTREE_PARAMS *octpar;
};




void print_usage( char *prog )
{
  printf( "Usage: %s [ options ]\n\n", prog );

  printf( "Options:\n" );

  printf( "\t--pdb filename      : input PDB file ( mandatory )\n" );
  printf( "\t--fpdb filename     : fixed part of the input PDB ( optional, default = filename ( -fixed.pdb ) derived from input PDB )\n" );
  printf( "\t--nofixed           : input PDB has no fixed part ( optional, default = has fixed part ( see --fpdb ) )\n" );
  printf( "\t--psf filename      : CHARMM PSF file for the input PDB ( optional, default = filename ( .psf ) derived from input PDB )\n" );
  printf( "\t--mol2 filename     : MOL2 file for the input PDB ( optional, for --hbond only, default = filename ( .mol2 ) derived from input PDB )\n\n" );

  printf( "\t--outn filename     : output PDB file from nblists minimization ( optional, default = filename ( -outn.pdb ) derived from input PDB )\n" );
  printf( "\t--outo filename     : output PDB file from octree minimization ( optional, default = filename ( -outo.pdb ) derived from input PDB )\n\n" );

  printf( "\t--prm filename      : PRM file ( optional, default = %s )\n", DEFAULT_PRM );
  printf( "\t--rtf filename      : RTF file ( optional, default = %s )\n", DEFAULT_RTF );
  printf( "\t--aprm filename     : APRM file ( optional, default = %s )\n\n", DEFAULT_APRM );

  printf( "\t--iter number       : maximum number of iterations of the minimizer ( optional, default = %d )\n\n", MAX_ITER );

  printf( "\t--nonblist          : no nblists based minimization ( optional, default = perform nblists based minimization )\n" );
  printf( "\t--nooct             : no octree based minimization ( optional, default = perform octree based minimization )\n\n" );

  printf( "\t--vdw               : use van der Waals energy ( optional, default = do not use vdw )\n" );
  printf( "\t--elec              : use electrostatics ( optional, default = do not use electrostatics )\n" );
  printf( "\t--hbond             : use hbond energy ( optional, default = do not use hbond )\n\n" );

  printf( "\t--vdwcut value(s)   : distance cutoff(s) ( in angstroms ) for van der Waals energy ( optional, default = %.2lf )\n", DEFAULT_VDW_DIST_CUTOFF );
  printf( "\t--eleccut value(s)  : distance cutoff(s) ( in angstroms ) for electrostatics ( optional, default = %.2lf )\n", DEFAULT_ELEC_DIST_CUTOFF );
  printf( "\t--hbondcut value(s) : distance cutoff(s) ( in angstroms ) for hbond energy ( optional, default = %.2lf )\n\n", DEFAULT_HBOND_DIST_CUTOFF );

  printf( "\t--vdwaprx value     : distance cutoff approximation factor ( in ( 0, 1 ] ) for van der Waals energy ( optional, default = 1.0 )\n" );
  printf( "\t--elecaprx value    : distance cutoff approximation factor ( in ( 0, 1 ] ) for electrostatics ( optional, default = 1.0 )\n\n" );

  printf( "\t--help              : print this help screen\n\n" );
}


int file_exists( char *fname, const char *msg )
{
  FILE *fp = fopen( fname, "rt" );
  
  if ( fp == NULL )
    {
      print_error( "Failed to open %s ( please specify %s filename )!", fname, msg );
      return 0;
    }
  else 
    {
      fclose( fp );
      return 1;
    }  
}


int read_command_line( int argc, char *argv[ ], char **pdbFile, char **pdbFixedFile, char **mol2File, char **psfFile, char **outnFile, char **outoFile,
                       char **prmFile, char **rtfFile, char **aprmFile, int *maxIter,
                       int *noOct, int *noNblist,
                       int *useVdw, int *useElec, int *useHbond,
                       double **dcVdw, double **dcElec, double **dcHbond, int *ndc,
                       double *aprxVdw, double *aprxElec )
{
  if ( argc <= 1 )
    {
      print_usage( argv[ 0 ] );
      exit( 0 );
    }

  ( *pdbFile ) = ( *pdbFixedFile ) = ( *mol2File ) = ( *psfFile ) = ( *outnFile ) = ( *outoFile ) = ( *prmFile ) = ( *rtfFile ) = ( *aprmFile ) = NULL;
  *maxIter = MAX_ITER;
  *noOct = *noNblist = 0;
  *useVdw = *useElec = *useHbond = 0;
  int ndcVdw = 0, ndcElec = 0, ndcHbond = 0;
  *ndc = 0;
  *aprxVdw = *aprxElec = 1.0;

  int i, j, nofixed = 0;

  for ( i = 1; i < argc; )
    {
     j = i;

     if ( !strcasecmp( argv[ i ], "--nofixed" ) )
       {
        nofixed = 1;

        i++;

        if ( i >= argc ) break;
       }

     if ( !strcasecmp( argv[ i ], "--nooct" ) )
       {
        *noOct = 1;

        i++;

        if ( i >= argc ) break;
       }

     if ( !strcasecmp( argv[ i ], "--nonblist" ) )
       {
        *noNblist = 1;

        i++;

        if ( i >= argc ) break;
       }

     if ( !strcasecmp( argv[ i ], "--vdw" ) )
       {
        *useVdw = 1;

        i++;

        if ( i >= argc ) break;
       }

     if ( !strcasecmp( argv[ i ], "--elec" ) )
       {
        *useElec = 1;

        i++;

        if ( i >= argc ) break;
       }

     if ( !strcasecmp( argv[ i ], "--hbond" ) )
       {
        *useHbond = 1;

        i++;

        if ( i >= argc ) break;
       }

     if ( !strcasecmp( argv[ i ], "--vdwcut" ) )
       {
        if ( ( i + 1 >= argc ) || ( argv[ i + 1 ][ 0 ] == '-' ) )
          {
           print_error( "Missing distance cutoff for vdw ( specify --vdwcut value(s) )!" );
           return 0;
          }
          
        i++;
        
        while ( ( i < argc ) && ( argv[ i ][ 0 ] != '-' ) )
          {
            if ( ndcVdw ) ( *dcVdw ) = ( double * ) _mol_realloc( *dcVdw, ( ndcVdw + 1 ) * sizeof( double ) );
            else ( *dcVdw ) = ( double * ) _mol_malloc( ( ndcVdw + 1 ) * sizeof( double ) );
            
            ( *dcVdw )[ ndcVdw ] = atof( argv[ i++ ] );
            
            if ( ( *dcVdw )[ ndcVdw ] <= 0 )
              {
               print_error( "Invalid distance cutoff for vdw ( %lf )!", ( *dcVdw )[ ndcVdw ] );
               return 0;
              }            
              
            ndcVdw++;  
          }  
          
        if ( i >= argc ) break;
       }


     if ( !strcasecmp( argv[ i ], "--eleccut" ) )
       {
        if ( ( i + 1 >= argc ) || ( argv[ i + 1 ][ 0 ] == '-' ) )
          {
           print_error( "Missing distance cutoff for electrostatics ( specify --eleccut value(s) )!" );
           return 0;
          }
          
        i++;
        
        while ( ( i < argc ) && ( argv[ i ][ 0 ] != '-' ) )
          {
            if ( ndcElec ) ( *dcElec ) = ( double * ) _mol_realloc( *dcElec, ( ndcElec + 1 ) * sizeof( double ) );
            else ( *dcElec ) = ( double * ) _mol_malloc( ( ndcElec + 1 ) * sizeof( double ) );
            
            ( *dcElec )[ ndcElec ] = atof( argv[ i++ ] );
            
            if ( ( *dcElec )[ ndcElec ] <= 0 )
              {
               print_error( "Invalid distance cutoff for electrostatics ( %lf )!", ( *dcElec )[ ndcElec ] );
               return 0;
              }            
              
            ndcElec++;  
          }  
          
        if ( i >= argc ) break;
       }


     if ( !strcasecmp( argv[ i ], "--hbondcut" ) )
       {
        if ( ( i + 1 >= argc ) || ( argv[ i + 1 ][ 0 ] == '-' ) )
          {
           print_error( "Missing distance cutoff for hbond ( specify --hbondcut value(s) )!" );
           return 0;
          }
          
        i++;
        
        while ( ( i < argc ) && ( argv[ i ][ 0 ] != '-' ) )
          {
            if ( ndcHbond ) ( *dcHbond ) = ( double * ) _mol_realloc( *dcHbond, ( ndcHbond + 1 ) * sizeof( double ) );
            else ( *dcHbond ) = ( double * ) _mol_malloc( ( ndcHbond + 1 ) * sizeof( double ) );
            
            ( *dcHbond )[ ndcHbond ] = atof( argv[ i++ ] );
            
            if ( ( *dcHbond )[ ndcHbond ] <= 0 )
              {
               print_error( "Invalid distance cutoff for vdw ( %lf )!", ( *dcHbond )[ ndcHbond ] );
               return 0;
              }            
              
            ndcHbond++;  
          }  
          
        if ( i >= argc ) break;
       }


     if ( !strcasecmp( argv[ i ], "--vdwaprx" ) )
       {
        if ( ( i + 1 >= argc ) || ( argv[ i + 1 ][ 0 ] == '-' ) )
          {
           print_error( "Missing distance cutoff approximation factor for vdw ( specify --vdwaprx value )!" );
           return 0;
          }

        *aprxVdw = atof( argv[ i + 1 ] );

        if ( ( *aprxVdw <= 0 ) || ( *aprxVdw > 1 ) )
          {
           print_error( "Invalid distance cutoff approximation factor for vdw ( %lf )!", *aprxVdw );
           return 0;
          }

        i += 2;

        if ( i >= argc ) break;
       }

     if ( !strcasecmp( argv[ i ], "--elecaprx" ) )
       {
        if ( ( i + 1 >= argc ) || ( argv[ i + 1 ][ 0 ] == '-' ) )
          {
           print_error( "Missing distance cutoff approximation factor for electrostatics ( specify --vdwaprx value )!" );
           return 0;
          }

        *aprxElec = atof( argv[ i + 1 ] );

        if ( ( *aprxElec <= 0 ) || ( *aprxElec > 1 ) )
          {
           print_error( "Invalid distance cutoff approximation factor for electrostatics ( %lf )!", *aprxElec );
           return 0;
          }

        i += 2;

        if ( i >= argc ) break;
       }

     if ( !strcasecmp( argv[ i ], "--iter" ) )
       {
        if ( ( i + 1 >= argc ) || ( argv[ i + 1 ][ 0 ] == '-' ) )
          {
           print_error( "Missing number of iterations ( specify --iter number )!" );
           return 0;
          }

        *maxIter = atoi( argv[ i + 1 ] );

        if ( *maxIter < 0 )
          {
           print_error( "Invalid number of iterations ( %d )!", *maxIter );
           return 0;
          }

        i += 2;

        if ( i >= argc ) break;
       }


     if ( !strcasecmp( argv[ i ], "--pdb" ) )
       {
        if ( ( i + 1 >= argc ) || ( argv[ i + 1 ][ 0 ] == '-' ) )
          {
           print_error( "Missing PDB file name ( specify --pdb pdb-file-name )!" );
           return 0;
          }

        ( *pdbFile ) = strdup( argv[ i + 1 ] );

        i += 2;

        if ( i >= argc ) break;
       }


     if ( !strcasecmp( argv[ i ], "--fpdb" ) )
       {
        if ( ( i + 1 >= argc ) || ( argv[ i + 1 ][ 0 ] == '-' ) )
          {
           print_error( "Missing fixed PDB file name ( specify --fpdb pdb-file-name )!" );
           return 0;
          }

        ( *pdbFixedFile ) = strdup( argv[ i + 1 ] );

        i += 2;

        if ( i >= argc ) break;
       }


     if ( !strcasecmp( argv[ i ], "--psf" ) )
       {
        if ( ( i + 1 >= argc ) || ( argv[ i + 1 ][ 0 ] == '-' ) )
          {
           print_error( "Missing PSF file name ( specify --psf psf-file-name )!" );
           return 0;
          }

        ( *psfFile ) = strdup( argv[ i + 1 ] );

        i += 2;

        if ( i >= argc ) break;
       }


     if ( !strcasecmp( argv[ i ], "--outn" ) )
       {
        if ( ( i + 1 >= argc ) || ( argv[ i + 1 ][ 0 ] == '-' ) )
          {
           print_error( "Missing output ( from nblists minimization ) PDB file name ( specify --outn pdb-file-name )!" );
           return 0;
          }

        ( *outnFile ) = strdup( argv[ i + 1 ] );

        i += 2;

        if ( i >= argc ) break;
       }


     if ( !strcasecmp( argv[ i ], "--outo" ) )
       {
        if ( ( i + 1 >= argc ) || ( argv[ i + 1 ][ 0 ] == '-' ) )
          {
           print_error( "Missing output ( from octree minimization ) PDB file name ( specify --outo pdb-file-name )!" );
           return 0;
          }

        ( *outoFile ) = strdup( argv[ i + 1 ] );

        i += 2;

        if ( i >= argc ) break;
       }


     if ( !strcasecmp( argv[ i ], "--prm" ) )
       {
        if ( ( i + 1 >= argc ) || ( argv[ i + 1 ][ 0 ] == '-' ) )
          {
           print_error( "Missing PRM file name ( specify --prm prm-file-name )!" );
           return 0;
          }

        ( *prmFile ) = strdup( argv[ i + 1 ] );

        i += 2;

        if ( i >= argc ) break;
       }


     if ( !strcasecmp( argv[ i ], "--rtf" ) )
       {
        if ( ( i + 1 >= argc ) || ( argv[ i + 1 ][ 0 ] == '-' ) )
          {
           print_error( "Missing output RTF file name ( specify --rtf rtf-file-name )!" );
           return 0;
          }

        ( *rtfFile ) = strdup( argv[ i + 1 ] );

        i += 2;

        if ( i >= argc ) break;
       }


     if ( !strcasecmp( argv[ i ], "--aprm" ) )
       {
        if ( ( i + 1 >= argc ) || ( argv[ i + 1 ][ 0 ] == '-' ) )
          {
           print_error( "Missing output APRM file name ( specify --aprm aprm-file-name )!" );
           return 0;
          }

        ( *aprmFile ) = strdup( argv[ i + 1 ] );

        i += 2;

        if ( i >= argc ) break;
       }


     if ( !strcasecmp( argv[ i ], "-h" ) || !strcasecmp( argv[ i ], "-help" ) || !strcasecmp( argv[ i ], "--help" ) )
       {
        print_usage( argv[ 0 ] );
        exit( 0 );
       }


     if ( i == j )
       {
        print_error( "Unknown option ( %s )!", argv[ i ] );
        return 0;
       }
    }

   int l;

   if ( ( *pdbFile ) == NULL )
     {
      print_error( "Missing input PDB file!" );
      return 0;
     }
   else
     {
      l = strlen( *pdbFile );

      if ( ( ( *pdbFile )[ l - 4 ] != '.' )
        || ( ( ( *pdbFile )[ l - 3 ] != 'p' ) && ( ( *pdbFile )[ l - 3 ] != 'P' ) )
        || ( ( ( *pdbFile )[ l - 2 ] != 'd' ) && ( ( *pdbFile )[ l - 2 ] != 'D' ) )
        || ( ( ( *pdbFile )[ l - 1 ] != 'b' ) && ( ( *pdbFile )[ l - 1 ] != 'B' ) ) )
         {
          print_error( "PDB file extension missing ( %s )!", ( *pdbFile ) );
          return 0;
         }
     }

   ( *pdbFile )[ l - 4 ] = 0;

   if ( ( *pdbFixedFile ) == NULL )
     {
       char s[ 1000 ];

       sprintf( s, "%s-fixed.pdb", ( *pdbFile ) );

       ( *pdbFixedFile ) = strdup( s );
     }

   if ( nofixed )
     {
       if ( ( *pdbFixedFile ) != NULL ) free( *pdbFixedFile );
       ( *pdbFixedFile ) = NULL;
     }

   if ( ( *psfFile ) == NULL )
     {
       char s[ 1000 ];

       sprintf( s, "%s.psf", ( *pdbFile ) );

       ( *psfFile ) = strdup( s );
     }

   if ( ( *mol2File ) == NULL )
     {
       char s[ 1000 ];

       sprintf( s, "%s.mol2", ( *pdbFile ) );

       ( *mol2File ) = strdup( s );
     }

   if ( ( *outnFile ) == NULL )
     {
       char s[ 1000 ];

       sprintf( s, "%s-outn.pdb", ( *pdbFile ) );

       ( *outnFile ) = strdup( s );
     }

   if ( ( *outoFile ) == NULL )
     {
       char s[ 1000 ];

       sprintf( s, "%s-outo.pdb", ( *pdbFile ) );

       ( *outoFile ) = strdup( s );
     }

   ( *pdbFile )[ l - 4 ] = '.';

   if ( ( *prmFile ) == NULL ) ( *prmFile ) = strdup( DEFAULT_PRM );

   if ( ( *rtfFile ) == NULL ) ( *rtfFile ) = strdup( DEFAULT_RTF );

   if ( ( *aprmFile ) == NULL ) ( *aprmFile ) = strdup( DEFAULT_APRM );

   *ndc = ndcVdw;
   if ( ndcElec > *ndc ) *ndc = ndcElec;
   if ( ndcHbond > *ndc ) *ndc = ndcHbond;
   
   if ( *ndc < 1 ) *ndc = 1;
   
   if ( ndcVdw < *ndc )
     {
       if ( ndcVdw ) ( *dcVdw ) = ( double * ) _mol_realloc( *dcVdw, ( *ndc ) * sizeof( double ) );
       else ( *dcVdw ) = ( double * ) _mol_malloc( ( *ndc ) * sizeof( double ) );
       
       for ( i = ndcVdw; i < *ndc; i++ )
         ( *dcVdw )[ i ] = DEFAULT_VDW_DIST_CUTOFF;
     }

   if ( ndcElec < *ndc )
     {
       if ( ndcElec ) ( *dcElec ) = ( double * ) _mol_realloc( *dcElec, ( *ndc ) * sizeof( double ) );
       else ( *dcElec ) = ( double * ) _mol_malloc( ( *ndc ) * sizeof( double ) );
       
       for ( i = ndcElec; i < *ndc; i++ )
         ( *dcElec )[ i ] = DEFAULT_ELEC_DIST_CUTOFF;
     }

   if ( ndcHbond < *ndc )
     {
       if ( ndcHbond ) ( *dcHbond ) = ( double * ) _mol_realloc( *dcHbond, ( *ndc ) * sizeof( double ) );
       else ( *dcHbond ) = ( double * ) _mol_malloc( ( *ndc ) * sizeof( double ) );
       
       for ( i = ndcHbond; i < *ndc; i++ )
         ( *dcHbond )[ i ] = DEFAULT_HBOND_DIST_CUTOFF;
     }

   printf( "\nParemeters:\n" );
   printf( "    --pdb         %s\n", ( *pdbFile ) );
   if ( nofixed ) printf( "    --nofixed\n" );
   else printf( "    --fpdb        %s\n", ( *pdbFixedFile ) );
   printf( "    --psf         %s\n", ( *psfFile ) );
   if ( *useHbond ) printf( "    --mol2        %s\n", ( *mol2File ) );
   if ( !( *noNblist ) ) printf( "    --outn        %s\n", ( *outnFile ) );
   if ( !( *noOct ) )printf( "    --outo        %s\n", ( *outoFile ) );
   printf( "    --prm         %s\n", ( *prmFile ) );
   printf( "    --rtf         %s\n", ( *rtfFile ) );
   printf( "    --aprm        %s\n", ( *aprmFile ) );
   printf( "    --iter        %d\n", *maxIter );
   if ( *noNblist ) printf( "    --nonblist\n" );
   if ( *noOct ) printf( "    --nooct\n" );
   if ( *useVdw ) 
      {
        printf( "    --vdw\n");
        printf( "    --vdwcut      " );
        for ( i = 0; i < *ndc; i++ )
          printf( "%.2lf A%s", ( *dcVdw )[ i ], ( i == ( *ndc ) - 1 ) ? "\n" : "  " );
        printf( "    --vdwaprx     %.3lf\n", *aprxVdw );
      }        
   if ( *useElec ) 
      {
        printf( "    --elec\n" );
        printf( "    --eleccut     " );
        for ( i = 0; i < *ndc; i++ )
          printf( "%.2lf A%s", ( *dcElec )[ i ], ( i == ( *ndc ) - 1 ) ? "\n" : "  " );
        printf( "    --elecaprx    %.3lf\n", *aprxElec );
      }  
   if ( *useHbond ) 
      {
        printf( "    --hbond\n" );
        printf( "    --hbondcut    " );
        for ( i = 0; i < *ndc; i++ )
          printf( "%.2lf A%s", ( *dcHbond )[ i ], ( i == ( *ndc ) - 1 ) ? "\n" : "  " );        
      }  

   printf( "\n" );
   
   if ( !file_exists( *pdbFile, "--pdb" ) || ( !nofixed && !file_exists( *pdbFixedFile, "--fpdb" ) ) 
     || !file_exists( *psfFile, "--psf" ) || ( *useHbond && !file_exists( *mol2File, "--mol2" ) ) 
     || !file_exists( *prmFile, "--prm" ) || !file_exists( *rtfFile, "--rtf" ) 
     || !file_exists( *aprmFile, "--aprm" ) ) return 0;

   return 1;
}




int size_of_nblist( struct atomgrp *ag, struct agsetup *ags )
{
   struct nblist *nblst = ags->nblst;

   int s = sizeof( struct nblist );

   s += ag->natoms * ( 3 * sizeof( float ) + 2 * sizeof( int ) + sizeof( int * ) ) + nblst->npairs * sizeof( int );

   return s;
}



void read_fix( char *ffile, int *nfix, int **fix )
{
   int linesz = 91;
   char *buffer = (char*) mymalloc( sizeof( char ) * linesz );

   *nfix = 0;

   if ( ffile != NULL )
     {
       FILE* fp = fopen( ffile, "r" );

       if ( fp == NULL )
         {
           print_error( "Failed to open fixed PDB %s ( use the --nofixed option if the molecule has no fixed part )!", ffile );
           exit( 0 );
         }

       while ( fgets( buffer, linesz - 1, fp ) != NULL )
         {
           if ( !strncmp( buffer, "ATOM", 4 ) ) ( *nfix )++;
         }
       fclose( fp );
     }

   *fix = (int *) mymalloc( ( *nfix ) * sizeof( int ) );

   if ( ffile != NULL )
     {
       FILE *fp = fopen( ffile, "r" );
       int na = 0;
       while ( fgets( buffer, linesz - 1, fp ) != NULL )
         {
           if ( !strncmp( buffer, "ATOM", 4 ) ) ( *fix )[ na++ ] = atoi( buffer + 4 ) - 1;
         }
       fclose( fp );
     }

   free( buffer );
}


void my_en_grad( int n, double *inp, void *prms, double *en, double *grad )
{
    static int count = 0;
    struct my_par *prm = ( struct my_par * ) prms;
    struct atomgrp *mag = prm->ag;

    if ( inp != NULL)
      {
        if ( n != mag->nactives * 3 )
          {
            print_error( "Mismatch in vector length ( my_en_grad )!");
            exit(0);
          }
        array2ag(inp, mag);
      }

    int mqcheck = ( ( struct my_par * ) prms )->nbupdate;
    struct agsetup *mags = ( ( struct my_par * ) prms )->ags;

    OCTREE_PARAMS *octpar = ( ( struct my_par * ) prms )->octpar;

    if ( octpar != NULL ) reorganize_octree( octpar->octree_static, 1 );
    else
      {
        if ( mqcheck ) check_clusterupdate( mag, mags );
      }

    // insert energy terms here
    *en = 0;
    zero_grads( mag );

    if ( octpar == NULL )
      {
        if ( prm->useVdw ) vdweng( mag,  en, mags->nblst );

        if ( prm->useElec ) eleng( mag, prm->eface, en, mags->nblst );

        if ( prm->useHbond )
          {
            double en1 = 0;
            hbondeng( mag, &en1, mags->nblst );
            *en += en1;
          }
      }
    else
      {
        if ( prm->useVdw )
           *en += octree_accumulation_excluding_far( octpar->octree_static, octpar->octree_moving, prm->dcVdw, prm->aprxdcVdw, octpar->fixed_cull, octpar->trans,
                                                     octpar->proc_func_params, vdweng_octree_single_mol );
        if ( prm->useElec )
           *en += octree_accumulation_excluding_far( octpar->octree_static, octpar->octree_moving, prm->dcElec, prm->aprxdcElec, octpar->fixed_cull, octpar->trans,
                                                     octpar->proc_func_params, eleng_octree_single_mol );

        if ( prm->useHbond )
           *en += octree_accumulation_excluding_far( octpar->octree_static, octpar->octree_moving, prm->dcHbond, prm->dcHbond, octpar->fixed_cull, octpar->trans,
                                                     octpar->proc_func_params, hbondeng_octree_single_mol );
      }

    beng(mag,  en);
    aeng(mag,  en);
    ieng(mag,  en);
//    teng(mag,  en);

    count++;

    int i;
    if ( grad != NULL )
      {
        for ( i = 0; i < n / 3 ; i++ )
          {
            grad[ 3 * i ] = -1 * mag->atoms[ mag->activelist[ i ] ].GX;
            grad[ 3 * i + 1 ] = -1 * mag->atoms[ mag->activelist[ i ] ].GY;
            grad[ 3 * i + 2 ] = -1 * mag->atoms[ mag->activelist[ i ] ].GZ;
          }
      }
}




int main( int argc, char *argv[ ] )
{
   char *pdbFile, *pdbFixedFile, *psfFile, *mol2File, *outnFile, *outoFile, *prmFile, *rtfFile, *aprmFile;
   int maxIter;
   int noOct, noNblist;
   int useVdw, useElec, useHbond;
   double *dcVdw, *dcElec, *dcHbond;
   int ndc;
   double aprxVdw, aprxElec;

   if ( !read_command_line( argc, argv, &pdbFile, &pdbFixedFile, &mol2File, &psfFile, &outnFile, &outoFile, &prmFile, &rtfFile, &aprmFile, &maxIter,
                            &noOct, &noNblist, &useVdw, &useElec, &useHbond, &dcVdw, &dcElec, &dcHbond, &ndc, &aprxVdw, &aprxElec ) ) return 1;

   if ( noOct && noNblist )
     {
       print_error( "Nothing to do ( both --nooct and --nonblist are set )!" );
       return 1;
     }

   if ( !useVdw && !useElec && !useHbond )
     {
       print_error( "Nothing to do ( at least one of --vdw, --elec and --hbond must be set )!" );
       return 1;
     }

   printf( "Processing %s... ", aprmFile ); fflush( stdout );
   struct prm *aprm = read_prm( aprmFile, "0.0.6" );   // read in atomic parameters
   printf( "done\n" );

   printf( "Processing %s... ", pdbFile ); fflush( stdout );
   struct atomgrp *ag = read_pdb( pdbFile, aprm );     // read in input file and it's forcefield parameters
   printf( "done\n" );

   double start, end, charmmTime;

   printf( "Processing %s, %s and %s... ", psfFile, prmFile, rtfFile ); fflush( stdout );
   start = clock( );
   read_ff_charmm( psfFile, prmFile, rtfFile, ag );    // read in CHARMM forcefield parameters
   end = clock( );   
   charmmTime = ( ( double ) ( end - start ) ) / CLOCKS_PER_SEC;
   printf( "%.3lf sec\n", charmmTime );
   
   if ( useHbond )
     {
       printf( "Processing %s... ", mol2File ); fflush( stdout );
       if ( !read_hybridization_states_from_mol2( mol2File, pdbFile, ag ) ) return 1;  // read in hybridization states
       printf( "done\n" );
      
       printf( "Fixing acceptor bases... " ); fflush( stdout );
       fix_acceptor_bases( ag, aprm );  // find base1 and base2 of each acceptor
       printf( "done\n" );
     }

   int nfix, *fix = ( int * ) mymalloc( ag->natoms * sizeof( int ) );

   if ( pdbFixedFile != NULL )
     {
       printf( "Processing %s... ", pdbFixedFile ); fflush( stdout );
     }
   read_fix( pdbFixedFile, &nfix, &fix );  // read in fixed part of the molecule
   fixed_init( ag );
   fixed_update( ag, nfix, fix );
   if ( pdbFixedFile != NULL ) printf( "done\n" );


   zero_grads( ag );  	// finish structure initialization
   fill_ingrp( ag );

   assign_combined_residue_sequence( ag );
   // storing current active list of atoms in array
   int ndim = ag->nactives * 3;
   double *startag = (double*) mymalloc( ndim * sizeof( double ) );
   ag2array( startag, ag );


   double nblist_cons = 0, octree_cons = 0;
   int snblist = 0, soctree = 0;

   struct agsetup *ags = (struct agsetup *) mymalloc( sizeof( struct agsetup ) );
   init_nblst( ag, ags );
   
   int k;

   for ( k = 0; k < ndc; k++ )
     {
       printf( "\n### START: RUN %d ###\n\n", k + 1 );
     
       ags->nblst->nbcof = 0;
    
       if ( useVdw && ( dcVdw[ k ] > ags->nblst->nbcof ) ) ags->nblst->nbcof = dcVdw[ k ];
       if ( useElec && ( dcElec[ k ] > ags->nblst->nbcof ) ) ags->nblst->nbcof = dcElec[ k ];
       if ( useHbond && ( dcHbond[ k ] > ags->nblst->nbcof ) ) ags->nblst->nbcof = dcHbond[ k ];
    
       ags->nblst->nbcut = ags->nblst->nbcof + 1;
    
       if ( !noNblist )
         {
           printf( "Building nblists... " ); fflush( stdout );
           // setup nblist information
           start = clock( );
           update_nblst( ag, ags ); // update nblst
           end = clock( );
           nblist_cons = ( ( double ) ( end - start ) ) / CLOCKS_PER_SEC;
           printf( "%.3lf sec\n", nblist_cons );
    
           snblist = size_of_nblist( ag, ags );
         }

       OCTREE octree;
      
       if ( !noOct )
         {
           printf( "Building octree... " ); fflush( stdout );
           // build octree
           start = clock( );
           if ( !build_octree( &octree, 60, 6, 1.0, ag ) )
             {
               print_error( "Failed to build octree!" );
               return 1;
             }
           end = clock( );
           octree_cons = ( ( double ) ( end - start ) ) / CLOCKS_PER_SEC;
           printf( "%.3lf sec\n", octree_cons );
      
           soctree = get_octree_size( &octree );
         }
      
       printf( "\n" );
    
       OCTREE_PARAMS eng_params;
    
       eng_params.ags = ags;
       eng_params.eps = 1.0;
       eng_params.octree_static = &octree;
       eng_params.octree_moving = &octree;
       eng_params.hdist_cutoff = MAX_R;
       eng_params.fixed_cull = 1;
       eng_params.trans = NULL;
       eng_params.engcat = NULL;
       eng_params.proc_func_params = &eng_params;
    
       struct my_par mupa;
    
       mupa.nbupdate = 1;
       mupa.ag = ag;
       mupa.ags = ags;
       mupa.efacv = 1.0;
       mupa.eface = 1.0;
       mupa.useVdw = useVdw;
       mupa.useElec = useElec;
       mupa.useHbond = useHbond;
       mupa.dcVdw = dcVdw[ k ];
       mupa.dcElec = dcElec[ k ];
       mupa.dcHbond = dcHbond[ k ];
       mupa.aprxdcVdw = aprxVdw * dcVdw[ k ];
       mupa.aprxdcElec = aprxElec * dcElec[ k ];
    
       double nblist_init_E = 0, octree_init_E = 0, nblist_E = 0, octree_E = 0, nblist_min = 0, octree_min = 0;
       const char *min_method[ ] = { "LBFGS" };
    
       if ( !noNblist )
         {
           mupa.octpar = NULL;
           my_en_grad( 0, NULL, ( void * )( &mupa ), &nblist_init_E, NULL );
         
           printf( "\nApplying %s with NBLISTS ( initial energy = %lf kcal/mol )...\n\n", min_method[ 0 ], nblist_init_E ); fflush( stdout );
           start = clock( );
           minimize_ag( MOL_LBFGS, maxIter, 1E-5, ag, ( void * )( &mupa ), my_en_grad );
           end = clock( );
           nblist_min = ( ( double ) ( end - start ) ) / CLOCKS_PER_SEC;
           my_en_grad( 0, NULL, ( void * )( &mupa ), &nblist_E, NULL );
           printf("\ndone ( time = %.3f sec, final energy = %lf kcal/mol )\n\n", nblist_min, nblist_E );
    
           printf( "Writing %s... ", outnFile ); fflush( stdout );
           write_pdb_nopar( ag, pdbFile, outnFile );
           printf( "done\n" );
    
           array2ag( startag, ag );
           update_nblst( ag, ags );
           zero_grads( ag );                      
         }
    
       if ( !noOct )
         {
           mupa.octpar = &eng_params;
           my_en_grad( 0, NULL, ( void * )( &mupa ), &octree_init_E, NULL );
    
           printf( "\nApplying %s with OCTREES ( initial energy = %lf kcal/mol )...\n\n", min_method[ 0 ], octree_init_E ); fflush( stdout );
           start = clock( );
           minimize_ag( MOL_LBFGS, maxIter, 1E-5, ag, ( void * )( &mupa ), my_en_grad );
           end = clock( );
           octree_min = ( ( double ) ( end - start ) ) / CLOCKS_PER_SEC;
           my_en_grad( 0, NULL, ( void * )( &mupa ), &octree_E, NULL );
           printf("\ndone ( time = %.3f sec, final energy = %f kcal/mol )\n\n", octree_min, octree_E );
    
           printf( "Writing %s... ", outoFile ); fflush( stdout );
           write_pdb_nopar( ag, pdbFile, outoFile );
           printf( "done\n" );
           
           array2ag( startag, ag );    
           zero_grads( ag );                       
           
           destroy_octree( &octree );                            
         }
        
       printf( "\n\nMolecule info:\n" );                 
       printf( "\tnatoms  = %d\n", ag->natoms );
       printf( "\tnfixed  = %d ( %.2lf \% )\n", nfix, ( nfix * 100.0 ) / ag->natoms );
       printf( "\tnmoving = %d  ( %.2lf \% )\n\n", ag->natoms - nfix, ( ( ag->natoms - nfix ) * 100.0 ) / ag->natoms );

       printf( "Distance cutoffs:\n" );         
       if ( useVdw )    printf( "\tvdw    = %.2lf A\n", dcVdw[ k ] );         
       if ( useElec )   printf( "\telec   = %.2lf A\n", dcElec[ k ] );         
       if ( useHbond )  printf( "\thbond  = %.2lf A\n", dcHbond[ k ] );                       
       if ( !noNblist ) printf( "\tnblist = %.2lf A\n", ags->nblst->nbcof );
       printf( "\n" );
    
       if ( !noNblist ) 
         {  
           printf( "NBLIST results:\n" );
           printf( "\tconstruction time = %.3lf sec\n", nblist_cons );
           printf( "\trunning time = %.3lf sec\n", nblist_min );       
           printf( "\tsize = %d bytes ( %.2lf KB, %.2lf MB ), npairs = %d\n",
                    snblist, snblist / 1024.0, snblist / ( 1024.0 * 1024.0 ), ags->nblst->npairs );
           printf( "\tinitial energy = %lf kcal/mol\n", nblist_init_E );                         
           printf( "\tfinal energy = %lf kcal/mol\n", nblist_E );         
           printf( "\n" );         
         }           
    
       if ( !noOct ) 
         {
           printf( "OCTREE results:\n" );     
           printf( "\tconstruction time = %.3lf sec\n", octree_cons );
           printf( "\trunning time = %.3lf sec\n", octree_min );                   
           printf( "\tsize = %d bytes ( %.2lf KB, %.2lf MB )\n",
                    soctree, soctree / 1024.0, soctree / ( 1024.0 * 1024.0 ) );
           printf( "\tinitial energy = %lf kcal/mol\n", octree_init_E );                                         
           printf( "\tfinal energy = %lf kcal/mol\n", octree_E );                         
           printf( "\n" );                                    
         }           
    
       if ( !noNblist && !noOct ) 
         {
           printf( "NBLIST vs. OCTREE:\n" );     
           
           if ( octree_cons > 0 ) printf( "\tconstruction-time( nblist ) / construction-time( octree ) = %.2lf\n", nblist_cons / octree_cons );
           else printf( "\tconstruction-time( nblist ) / construction-time( octree ) = inf\n" );
         
           if ( octree_min > 0 ) printf("\trunning-time( nblist ) / running-time( octree ) = %.2lf\n", nblist_min / octree_min );
           else printf("\trunning-time( nblist ) / running-time( octree ) = inf\n" );
           
           if ( soctree > 0 ) printf( "\tsize( nblist ) / size( octree ) = %.2lf\n", ( 1.0 * snblist ) / soctree );
           else printf( "\tsize( nblist ) / size( octree ) = inf\n" );       
           
           printf( "\n" );                                
         }  
         
       printf( "### END: RUN %d ###\n\n", k + 1 );         
     }

   free_agsetup( ags );

   return 0;
}

