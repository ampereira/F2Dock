/**-----------------------------------------------------------------------------                                                                                                  
**                                                                                                                                                                               
**  Copyright (C) : Structural Bioinformatics Laboratory, Boston University.                                                                                                                        
**                                                                                                                                                                               
**  This software was developed at the Boston University 2006-2011, by                                                                                                      
**  Structural Bioinformatics Laboratory, as part of NIH funded research.                                                                                                                      
**                                                                                                                                                                               
**  Explicit permission is hereby granted to US Universities and US                                                                                                     
**  Government supported Research Institutions to copy and modify this                                                                                                           
**  software for educational and research purposes, provided copies include                                                                                                      
**  this notice. This software (or modified copies thereof) may not be                                                                                                           
**  distributed to any other institution without express permission from the                                                                                                     
**  Structural Bioinformatics Laboratory and  Boston University. Requests to use this software (or                                                                                 **  modified copies therof) in any other way should be sent to Dima Kozakov,                                                                                                     
**  Department of Biomedical Engineering: "midas@bu.edu".                                                                                                                  
**                                                                                                                                                                               
**---------------------------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>

#include <libmol/libmol.h>

namespace LibMol{

#ifndef _DEBUG_
void* mymalloc (size_t size)
{
//	_PRINT_DEPRECATED_
//	fprintf (stderr, "\t(please write your own malloc wrapper)\n");

	void* v = (void*) malloc (size);
	if (v == NULL)
	{
		exit (EXIT_FAILURE);
	}
	return v;
}

void* mycalloc (size_t nmemb, size_t size)
{
//	_PRINT_DEPRECATED_
//	fprintf (stderr, "\t(please write your own calloc wrapper)\n");

	void* v = (void*) calloc (nmemb, size);
	if (v == NULL)
	{
		perror ("calloc"), exit (EXIT_FAILURE);
	}
	return v;
}

void* myrealloc (void* ptr, size_t size)
{
//	_PRINT_DEPRECATED_
//	fprintf (stderr, "\t(please write your own realloc wrapper)\n");

	if (size < 1)
	{
		fprintf (stderr, "warning: _mol_realloc called with size 0, no realloc will occur\n");
		return ptr;
		
	}
	void* v = (void*) realloc (ptr, size);
	if (v == NULL && ptr != NULL)
	{
		perror ("realloc"), exit (EXIT_FAILURE);
	}
	return v;
}
#endif /* _DEBUG_ */


#ifdef _DARWIN_

int
getline (char **lineptr, size_t *n, FILE *stream)
{
 char* read=(char*)malloc(10000*sizeof(char));
 char* result= fgets(read, 10000, stream);
 *lineptr=read;
 *n=strlen(read);
 if (result==NULL) return -1;
 return 0;
}

#endif /* _DARWIN_ */


void myexit (int status)
{
	exit (status);
}

FILE* myfopen (const char* path, const char* mode)
{
    if ( strncmp(path, "-", 2) == 0 ) {
        if (strncmp(mode, "r", 1) == 0 ) {
            return stdin;
        } else if (strncmp(mode, "w", 1) == 0) {
            return stdout;
        } else {
            fprintf(stderr, "case in myfopen not handled\n");
            exit( EXIT_FAILURE );
        }
    }
	FILE* fp = fopen (path, mode);
	if (fp == NULL) // file could not be opened
	{
		fprintf (stderr, "fopen error on %s:\n", path);
		perror ("fopen"), exit (EXIT_FAILURE);
	}

	return fp;
}

void myfclose (FILE* fp)
{
    if (fp == stdin || fp == stdout) {
        rewind(fp); //Ryan thinks this is a hack, David thinks it's beautiful
        return;     //if we are using a pdb from stdin, we often need to reuse
    }               //it, so we just rewind to the beginning

	int retval = fclose (fp);

	if (retval != 0) // file could not be closed
	{
		fprintf (stderr, "fclose error:\n");
		perror ("fclose"), exit (EXIT_FAILURE);
	}
}

int iswhiteline (char* line)
{
	int i;
	for (i = 0; line[i] != '\0'; i++)
	{
		int ch = line[i];
		if (ch != ' ' && ch != '\t' && ch != '\n')
		{
			return 0;
		}
	}
	return 1;
}

char* rstrip (char* string)
{
    size_t length = strlen(string);

    char* end = string + length - 1;

    while (end >= string && isspace(*end))
        end--;
    *(end+1) = '\0';

    return string;
}

/**
	Prints a formatted error message.
	[ Rezaul Chowdhury, Nov. 17, 2010 ]	
*/
void print_error( char *format, ... )
{
   char eMsg[ 500 ];
   va_list args;
   
   va_start( args, format );
   
   vsprintf( eMsg, format, args );
   
   va_end( args );
   
   printf( "\nError: %s\n\n", eMsg );   
   
   fflush( stdout );   
}

};
