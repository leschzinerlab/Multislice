/*
		*** memory.h ***

------------------------------------------------------------------------
Copyright 1998 Earl J. Kirkland
The computer code and or data in this file is provided for demonstration
purposes only with no guarantee or warranty of any kind that it is correct
or produces correct results.  By using the code and or data in this file
the user agrees to accept all risks and liabilities associated with the code
and or data. The computer code and or data in this file may be copied (and used)
for non-commercial academic or research purposes only, provided that this
notice is included. This file or any portion of it may not be resold, rented
or distributed without the written permission of the author.
------------------------------------------------------------------------

	common memory allocation routines

These basically call malloc() and do the error checking and 
print an error message.  This saves doing the checking
EVERY time and are mainly for convenience.

	char1D()	= allocate 1D char array
	int1D()		= allocate 1D int array
	float1D()	= allocate 1D float array
	double1D()	= allocate 1D double array

	char2D()	= allocate 2D char array
	short2D()   = allocate 2D short array
	long2D()    = allocate 2D long array
	float2D()	= allocate 2D float array
	double2D()  = allocate 2D double array

started 24-Aug-1995 Earl J. Kirkland
added short2D() 21-feb-1996 ejk
added char2D() 20-jun-96 ejk
added double2D() 16-jan-1997 ejk

*/

#include <stdlib.h>
#include <stdio.h>

/*---------------------------- char1D() -------------------------------*/
/*
	1D array allocator for type char
	make space for m[0...(n-1)]
	printf error message and exit if not successful
	
	this save checking for a NULL return etc every time 
	
*/
char* char1D( int n, const char *message );

/*---------------------------- int1D() -------------------------------*/
/*
	1D array allocator for type int
	make space for m[0...(n-1)]
	printf error message and exit if not successful
	
	this save checking for a NULL return etc every time 
	
*/
int* int1D( int n, const char *message );

/*---------------------------- float1D() -------------------------------*/
/*
	1D array allocator for type float
	make space for m[0...(n-1)]
	printf error message and exit if not successful
	
	this save checking for a NULL return etc every time 
	
*/
float *float1D( int n, const char *message );

/*---------------------------- double1D() -------------------------------*/
/*
	1D array allocator for type double
	make space for m[0...(n-1)]
	printf error message and exit if not successful
	
	this save checking for a NULL return etc every time 
	
*/
double* double1D( int n, const char *message );

/*---------------------------- char2D() -------------------------------*/
/*
	2D array allocator for type char
	make space for m[0...(nx-1)][0..(ny-1)]
	
*/
char **char2D( int nx, int ny, const char *message );

/*---------------------------- short2D() -------------------------------*/
/*
	2D array allocator for type short
	make space for m[0...(nx-1)][0..(ny-1)]
	
	message = char[] with error message
	
*/
short **short2D( int nx, int ny, const char *message );

/*---------------------------- long2D() -------------------------------*/
/*
	2D array allocator for type long
	make space for m[0...(nx-1)][0..(ny-1)]
	
	message = char[] with error message
	
*/
long **long2D( int nx, int ny, const char *message );

/*---------------------------- float2D() -------------------------------*/
/*
	2D array allocator for type float
	make space for m[0...(nx-1)][0..(ny-1)]

*/
float **float2D( int nx, int ny, const char *message );

/*---------------------------- double2D() -------------------------------*/
/*
	2D array allocator for type doubel
	make space for m[0...(nx-1)][0..(ny-1)]

*/
double **double2D( int nx, int ny, const char *message );
