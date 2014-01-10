/*
		*** memory.c ***

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
char* char1D( int n, const char *message )
{
	char *m;
	
	m = (char*) malloc( n * sizeof( char ) );
	if( m == NULL ) {
		printf("char1D() cannot allocate memory size=%d: %s\n",
		       n, message);
		exit( 0 );
	}
	return( m );

}  /* end char1D() */

/*---------------------------- int1D() -------------------------------*/
/*
	1D array allocator for type int
	make space for m[0...(n-1)]
	printf error message and exit if not successful
	
	this save checking for a NULL return etc every time 
	
*/
int* int1D( int n, const char *message )
{
	int *m;
	
	m = (int*) malloc( n * sizeof( int ) );
	if( m == NULL ) {
		printf("int1D() cannot allocate memory size=%d: %s\n",
		       n, message);
		exit( 0 );
	}
	return( m );

}  /* end int1D() */

/*---------------------------- float1D() -------------------------------*/
/*
	1D array allocator for type float
	make space for m[0...(n-1)]
	printf error message and exit if not successful
	
	this save checking for a NULL return etc every time 
	
*/
float *float1D( int n, const char *message )
{
	float *m;
	
	m = (float*) malloc( n * sizeof( float) );
	if( m == NULL ) {
		printf("float1D() cannot allocate memory size=%d: %s\n",
		       n, message);
		exit( 0 );
	}
	return( m );

}  /* end float1D() */

/*---------------------------- double1D() -------------------------------*/
/*
	1D array allocator for type double
	make space for m[0...(n-1)]
	printf error message and exit if not successful
	
	this save checking for a NULL return etc every time 
	
*/
double* double1D( int n, const char *message )
{
	double *m;
	
	m = (double*) malloc( n * sizeof( double ) );
	if( m == NULL ) {
		printf("double1D() cannot allocate memory size=%d: %s\n",
		       n, message);
		exit( 0 );
	}
	return( m );

} /* end double1D() */

/*---------------------------- char2D() -------------------------------*/
/*
	2D array allocator for type char
	make space for m[0...(nx-1)][0..(ny-1)]
	
*/
char **char2D( int nx, int ny, const char *message )
{	char **m;
	int i;

	m = (char**) malloc( nx * sizeof( char* ) ); 
	if( m == NULL ) {
		printf("char2D cannot allocate pointers, size=%d: %s\n",
		       nx, message );
		exit(0);
	}

	for (i=0; i<nx; i++){
		m[i] = (char*) malloc( ny * sizeof( char ) );
		if( m[i] == NULL ){
			printf("char2D cannot allocate arrays, size=%d: %s\n",
			       ny, message );
			exit(0);
		}
	}

	return m;
}  /* end char2D() */

/*---------------------------- short2D() -------------------------------*/
/*
	2D array allocator for type short
	make space for m[0...(nx-1)][0..(ny-1)]
	
	message = char[] with error message
	
*/
short **short2D( int nx, int ny, const char *message )
{	short **m;
	int i;

	m = (short**) malloc( nx * sizeof( short* ) ); 
	if( m == NULL ) {
		printf("short2D cannot allocate pointers, size=%d : %s\n",
		       nx, message );
		exit(0);
	}

	for (i=0; i<nx; i++){
		m[i] = (short *) malloc( ny * sizeof( short ) );
		if( m[i] == NULL ){
			printf("short2D cannot allocate arrays, size=%d : %s\n",
			       ny, message );
			exit(0);
		}
	}

	return m;
}  /* end short2d() */

/*---------------------------- long2D() -------------------------------*/
/*
	2D array allocator for type long
	make space for m[0...(nx-1)][0..(ny-1)]
	
	message = char[] with error message
	
*/
long **long2D( int nx, int ny, const char *message )
{	long **m;
	int i;

	m = (long**) malloc( nx * sizeof( long* ) ); 
	if( m == NULL ) {
		printf("long2D cannot allocate pointers, size=%d : %s\n",
		       nx, message );
		exit(0);
	}

	for (i=0; i<nx; i++){
		m[i] = (long *) malloc( ny * sizeof( long ) );
		if( m[i] == NULL ){
			printf("long2D cannot allocate arrays, size=%d : %s\n",
			       ny, message );
			exit(0);
		}
	}

	return m;
}  /* end long2d() */

/*---------------------------- float2D() -------------------------------*/
/*
	2D array allocator for type float
	make space for m[0...(nx-1)][0..(ny-1)]

*/
float **float2D( int nx, int ny, const char *message )
{	float **m;
	int i;

	m = (float**) malloc( nx * sizeof( float* ) ); 
	if( m == NULL ) {
		printf("float2D cannot allocate pointers, size=%d: %s\n",
		       nx, message );
		exit(0);
	}

	for (i=0; i<nx; i++){
		m[i] = (float *) malloc( ny * sizeof( float ) );
		if( m[i] == NULL ){
			printf("float2D cannot allocate arrays, size=%d: %s\n",
			       ny, message );
			exit(0);
		}
	}

	return m;

}  /* end float2D() */

/*---------------------------- double2D() -------------------------------*/
/*
	2D array allocator for type doubel
	make space for m[0...(nx-1)][0..(ny-1)]

*/
double **double2D( int nx, int ny, const char *message )
{	double **m;
	int i;

	m = (double**) malloc( nx * sizeof( double* ) ); 
	if( m == NULL ) {
		printf("double2D cannot allocate pointers, size=%d: %s\n",
		       nx, message );
		exit(0);
	}

	for (i=0; i<nx; i++){
		m[i] = (double *) malloc( ny * sizeof( double ) );
		if( m[i] == NULL ){
			printf("double2D cannot allocate arrays, size=%d: %s\n",
			       ny, message );
			exit(0);
		}
	}

	return m;

}  /* end double2D() */

