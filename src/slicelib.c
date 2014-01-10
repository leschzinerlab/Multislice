/*		*** slicelib.c ***

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

	function library for use with multislice programs

	askYN()       : ask a yes/no question and return true/false
	bessi0()      : modified Bessel function I0(x)
	bessk0()      : modified Bessel function K0(x)
	cputim()      : return current CPU time in sec
	featom()      : return scattering factor for a given atom
	freqn()       : calculate spatial frequencies
	parlay()      : parse the layer structure
	propagate()   : propagate a 2D wavefunction
	ranflat()     : return a random number with uniform distribution
	rangauss()    : return a random number with Gaussian distribution
	ReadfeTable() : read fe scattering factor table
	ReadLine()    : read a line of char from file
	ReadXYZcoord(): read a set of (x,y,z) coordinates from a file
	seval()       : Interpolate from cubic spline coefficients
	sigma()	      : return the interaction parameter
	splinh()      : fit a quasi-Hermite  cubic spline
	transmit()    : transmit a 2D wavefunction
	vatom()       : return real space atomic potential (NOT projected)
	vzatom()      : return real space projected atomic potential
	wavelength()  : return electron wavelength in A for given keV

	this file is formatted for a tab size of 8 characters

	started may 1996 E. Kirkland
	added askYN() 23-jun-1996 ejk
	add scattering factor routines featom, ReadfeTable, ReadLine
		seval, spinh	26-july-1996 ejk
	move random number generator here 3-aug-1996 ejk
	move freqn(), transmitt() and propagate() to here 4-aug-1996 ejk
	converted to new 12 parameter fe(k) from mcdf 16-jan-1997 ejk
	add bessk0(), bessi0(), vzatom() and wavelength() 18-feb-1997 ejk
	update random number generators 20-may-1997 ejk
	added more directories to look for "fparams.dat" 5-july-1997 ejk
	added vatom() and changed constants in vxatom() in about
			5th sig. fig.  24-nov-1997 ejk
	added sigma() 30-nov-1997 ejk
	change min Z to 1 for hydrogen 15-dec-1997 ejk
	fixed possible log(0) problem in rangaus() 10-jan-1998 ejk
	moved seval(), splinh() and ReadXYZcoord()
			into this library 11-jan-1998 ejk
*/

#include <stdio.h>	/* ANSI C libraries */
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#include "memory.h"	/* memory allocation routines */

/* for the atomic scattering factor tables */

#define NPMAX	12	/* number of parameters for each Z */
#define NZMIN	1	/* min Z in featom.tab */
#define NZMAX	103	/* max Z in featom.tab */
#define	NCMAX	132	/* characters per line to read */
int feTableRead=0;	/* flag to remember if the param file has been read */
double **fparams;	/* to get fe(k) parameters */
const int nl=3, ng=3;	/* number of Lorenzians and Gaussians */

/*--------------------- askYN() -----------------------------------*/
/*
	ask a yes/no question and return 1 or 0 for TRUE/FALSE

	message[] = question to ask
*/
int askYN( const char message[] )
{
	char cline[132];

	printf("%s (y/n) :\n", message);
	scanf("%s", cline );
	if( (strcmp(cline,"Y")==0) || (strcmp(cline,"y")==0) ||
		(strcmp(cline,"1")==0) ) return( 1 );  else return( 0 );

} /* end askYN() */

/*-------------------- bessi0() ---------------*/
/*
    modified Bessel function I0(x)
    see Abramowitz and Stegun page 379

    x = (double) real arguments

    12-feb-1997 E. Kirkland
 */
 double bessi0( double x )
 {
 	int i;
 	double ax, sum, t;
 	
 	double i0a[] = { 1.0, 3.5156229, 3.0899424, 1.2067492,
		0.2659732, 0.0360768, 0.0045813 };

 	double i0b[] = { 0.39894228, 0.01328592, 0.00225319,
 		-0.00157565, 0.00916281, -0.02057706, 0.02635537,
 		-0.01647633, 0.00392377};

	ax = fabs( x );
	if( ax <= 3.75 ) {
		t = x / 3.75;
		t = t * t;
		sum = i0a[6];
		for( i=5; i>=0; i--) sum = sum*t + i0a[i]; 
	} else {
		t = 3.75 / ax;
		sum = i0b[8];
		for( i=7; i>=0; i--) sum = sum*t + i0b[i];
		sum = exp( ax ) * sum / sqrt( ax );
	}
	return( sum );

}  /* end bessi0() */

/*-------------------- bessk0() ---------------*/
/*
    modified Bessel function K0(x)
    see Abramowitz and Stegun page 380
    
    Note: K0(0) is not define and this function
        returns 1E20
 
    x = (double) real arguments
    
    this routine calls bessi0() = Bessel function I0(x)
    
    12-feb-1997 E. Kirkland
 */
 double bessk0( double x )
 {
 	double bessi0(double);
 
 	int i;
 	double ax, x2, sum;
 	double k0a[] = { -0.57721566, 0.42278420, 0.23069756,
 		 0.03488590, 0.00262698, 0.00010750, 0.00000740};
 	        
 	double k0b[] = { 1.25331414, -0.07832358, 0.02189568,
 		 -0.01062446, 0.00587872, -0.00251540, 0.00053208};

	ax = fabs( x );
	if( (ax > 0.0)  && ( ax <=  2.0 ) ){
		x2 = ax/2.0;
		x2 = x2 * x2;
		sum = k0a[6];
		for( i=5; i>=0; i--) sum = sum*x2 + k0a[i];
		sum = -log(ax/2.0) * bessi0(x) + sum;
	} else if( ax > 2.0 ) {
		x2 = 2.0/ax;
		sum = k0b[6];
		for( i=5; i>=0; i--) sum = sum*x2 + k0b[i];
		sum = exp( -ax ) * sum / sqrt( ax );
	} else sum = 1.0e20;
	return ( sum );

}  /* end bessk0() */

/*--------------------- cputim() -----------------------------------*/
/*
   retrieve current CPU time in seconds
*/

double cputim()
{
	return ( ( (double)clock() ) / ( (double)CLOCKS_PER_SEC) );

}  /* end cputim() */

/*--------------------- featom() -----------------------------------*/
/*
	return the electron scattering factor for atomic
	number Z at scattering angle k

	Z = atomic number 1 <= Z <= 103
	k2  = k*k where k =1/d = scattering angle (in 1/A)

  assumed global vars:

#define NZMIN	1	 min Z in featom.tab 
#define NZMAX	103	 max Z in featom.tab 

int feTableRead=0; = flag to remember if the param file has been read 
int nl=3, ng=3; = number of Lorenzians and Gaussians 
double fparams[][] = fe parameters

*/

double featom( int Z, double k2 )
{
   int i, nfe;
   double sum;
   int ReadfeTable( );

   if( (Z<NZMIN) || (Z>NZMAX) ) return( 0.0 );

	/* read in the table from a file if this is the
		first time this is called */
	if( feTableRead == 0 ) nfe = ReadfeTable();

	sum = 0.0;

	/* Lorenztians */
	for( i=0; i<2*nl; i+=2 )
	    sum += fparams[Z][i]/( k2 + fparams[Z][i+1] );

	/* Gaussians */
	for( i=2*nl; i<2*(nl+ng); i+=2 )
	    sum += fparams[Z][i]*exp( - k2 * fparams[Z][i+1] );

	return( sum );

}  /* end featom() */

/*------------------------ freqn() ------------------------*/
/*
	Calculate spatial frequencies for use with fft's
	NOTE: zero freg is in the bottom left corner and
		expands into all other corners - not in the center
		this is required for fft - don't waste time rearanging

	This routine must be called once for each direction

	ko[n]  = real array to get spatial frequencies
	ko2[n] = real array to get k[i]*k[i]
	xo[n]  = real array to get positions 
	nk     = integer number of pixels
	ak     = real full scale size of image in pixels
*/
void freqn( float *ko, float *ko2, float *xo, int nk, double ak )
{
	int i, imid;

	imid = nk/2;

	for( i=0; i<nk; i++) {
		xo[i] = ((float) (i * ak) ) / ((float)(nk-1));
		if ( i > imid ) {
			ko[i]  = ((float)(i-nk)) / ((float)ak);
		} else {
			ko[i]  = ((float)i) / ((float)ak);
		}
		ko2[i] = ko[i] * ko[i];
	}

}  /*  end freqn() */

/*--------------------- parlay() -----------------------------------*/
/*
  subroutine to parse the atomic layer stacking sequence 
  for use with multislice programs.

  This converts layer structure definition of the form:
      2(abc)d
  into a sequence of numerical indices where a=1, b=2, c=3, ...
  The main attraction is the repeat operator n(...) where
  the structure inside the parenthesis (...) is repeated n times.
  for instance the above structure is translated into:
      1 2 3 1 2 3 4
  The parenthesis may be nested up to 100 levels (determined by
  nlmax). For instance  5( 2(abc) 3(cba) ) is also a valid structure
  definition. This is a compact way of specifying the layer structure
  for a multislice calculation. tab's may be present in the structure
  anywhere (they are ignored).

  This is done by pushing the position of each '(' and its repeat
  count onto a stack. Each ')' pops the last entry from the stack
  and invokes a duplication process.

  Layers refer to the distinquishable subset of different
  types of layers, whereas slices refer to the way these
  layers are sequentially stacked.

  fortran version started 22-dec-1989 earl j. kirkland
  added nested parenthesis by stacking entry points
     31-jan-1990 ejk
  added tab handling (in ascii and ebcdic) 1-feb-1990 ejk
  minor changes to error messages to make rs6000's happy
      7-july-1992 ejk
  converted to ANSI-C 26-jun-1995 ejk
  added include of stdlib.h 6-july-1995 ejk

   c             = input character string
   islice[nsmax] = integer array to get stacking sequence indicies
   nsmax         = (integer) size of layer
   lmax          = (integer) maximum allowed layer index
   nslice        = (integer) number of layers
   returned value= (integer) success/error code
                       0 : success
                      -1 : layer out of range
                      -2 : missing left parenthesis
                      -3 : parenthesis nested too deep
                      -4 : bad repeat code
                      -5 : unmatched right parenthesis
                      -6 : invalid character
                      -7 : too many layers
                      -8 : incomplete stacking sequence

   fperr       = (int) if this is not a NULL then write
                   error messages 

  NOTE: islice and nslice may be modified by this routine

  cname determines the mapping of characters into numbers.

*/

#define NSTKMAX 100	/* maximum stack depth (i.e. maximum level
				 of nesting the parenthesis) */
#define NCHARS 52	/* maximum number of character symbols */

int parlay( const char c[], int islice[], int nsmax, int lmax,
	   int *nslice, int fperr )
{
	/* define our own symbol sequence */
	const char cname[] =
	   "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";

	int ic, lenc, name, i, j, istack, n,
		ipoint[NSTKMAX], repeat[NSTKMAX];

	/*  initialize misc constants  */

        *nslice = 0;
        istack = 0;
        lenc = (int) strlen( c );  /* length of c */

	for( ic=0; ic<lenc; ic++ ) {	/* loop over all characters in c */

	   /*  skip over embedded spaces, tabs and wierd characters */
	   while( isspace( c[ic] )  ) ic++;

	   /*  if its a character do this  */		
	   if( isalpha( c[ic] ) ) {
		for(name=0; name<NCHARS; name++) 
		   if( c[ic] == cname[name] ) break;
		if( name > lmax ) {
		   if( fperr != 0 ) 
		      printf("Layer /%c/ out of range in the stacking sequence"
              		" in subroutine parlay.\n", c[ic] );
			return( -1 );
		} else {
		   if( *nslice >= nsmax ) {
			if( fperr != 0 ) 
              		printf("Too many layers generated in the stacking"
              			" sequence in parlay.\n");
        			return( -7 );
        	   }
			islice[ (*nslice)++ ] = name;
		}
        	
	/*  if its a number then extract repeat count up to the '('
    	and save it until a ')' appears (i.e. push on the stack)
	*/
	   } else if ( isdigit( c[ic] ) ) {
		if( istack >= NSTKMAX ) {
			if( fperr != 0 ) 
			printf("Parenthesis nested too deep in "
			   "the stacking sequence in subroutine parlay.\n");
			return( -3 );
		}
        	repeat[istack] = atoi( &c[ic] ) - 1; 
        	while( isdigit(c[ic]) || isspace(c[ic]) ) ic++;
		if( c[ic] != '(' ) {
			if( fperr != 0 ) 
			   printf("Missing left parenthesis in "
			   "the stacking sequence in subroutine parlay.'\n");
			   for(i=0; i<=ic; i++) printf("%c", c[i]);
			   printf("?\n");
        		   return( -2 );
		}
		ipoint[istack++] = *nslice;

		/*  if its a ')' then repeat back to the last '('
  				(i.e. pop the stack) */

  		} else if( c[ic] == ')' ) {
		   if( istack <= 0 ) {
			if( fperr != 0 ) 
			printf("Unmatched right parenthesis in "
				"the stacking sequence in subroutine parlay.\n");
 			return( -5 );
		   }
		   n = *nslice;
		   istack--;
		   for( j=0; j<repeat[istack]; j++)
			for(i=ipoint[istack]; i<n; i++){
			   if( *nslice >= nsmax ) {
                		if( fperr != 0 ) 
              			printf("Too many layers generated in the stacking"
              				" sequence in parlay.\n");
        				return( -7 );
			   }
				islice[ (*nslice)++ ] = islice[ i ];
			   }
		} else {
		   if( fperr != 0 ) 
              	   printf("Invalid character /%c/ encountered in "
			  "the stacking sequence in subroutine parlay.\n",
			  c[ic]);
 		   return( -6 );
 		}

	} /* end for( ic... */
	
	if( istack != 0 ) {
		if( fperr != 0 ) 
     		printf("incomplete stacking sequence in parlay.\n");
		return( -8 ); 
	} else return( 0 );

#undef NSTKMAX 
#undef NCHARS

}  /* end parlay() */

/*------------------------ propagate() ------------------------*/
/*
	propagate the wavefunction thru one layer

	waver,i[ix][iy]  = real and imaginary parts of wavefunction
	propxr,i[ix]     = real and imag. parts of x comp of propagator
	propyr,i[iy]     = real and imag. parts of y comp of propagator

	kx2[], ky2[]     = spatial frequency components
	k2max            = square of maximum k value 
	nx, ny           = size of array
	
	on entrance waver,i and 
		 propxr,i/propyr,i are in reciprocal space
	
	only waver,i will be changed by this routine
*/
void propagate( float** waver, float** wavei, 
	float* propxr, float* propxi, float* propyr, float* propyi,
	float* kx2, float* ky2, float k2max, int nx, int ny )
{
	int ix, iy;
	float pxr, pxi, pyr, pyi, wr, wi, tr, ti;

	/*  multiplied by the propagator function */

	for( ix=0; ix<nx; ix++) {
	   if( kx2[ix] < k2max ) {
		pxr = propxr[ix];
		pxi = propxi[ix];
		for( iy=0; iy<ny; iy++) {
			if( (kx2[ix] + ky2[iy]) < k2max ) {
				pyr = propyr[iy];
				pyi = propyi[iy];
				wr = waver[ix][iy];
				wi = wavei[ix][iy];
				tr = wr*pyr - wi*pyi;
				ti = wr*pyi + wi*pyr;
				waver[ix][iy] = tr*pxr - ti*pxi;
				wavei[ix][iy] = tr*pxi + ti*pxr;
			} else
				waver[ix][iy] = wavei[ix][iy] = 0.0F;
		} /* end for(iy..) */

	   } else for( iy=0; iy<ny; iy++)
		waver[ix][iy] = wavei[ix][iy] = 0.0F;
			
	} /* end for(ix..) */

} /* end propagate() */

/*---------------------------- ranflat -------------------------------*/
/*
	return a random number in the range 0.0->1.0
	with uniform distribution

	the 'Magic Numbers' are from 
		Numerical Recipes 2nd edit pg. 285
*/
double ranflat( unsigned long *iseed )
{
	static unsigned long a=1366, c=150889L, m=714025L;
	
	*iseed = ( a * (*iseed) + c ) % m;
	
	return( ((double) *iseed)/m);

}  /* end ranflat() */

/*-------------------- rangauss() -------------------------- */
/*
	Return a normally distributed random number with 
	zero mean and unit variance using Box-Muller method
	
	ranflat() is the source of uniform deviates

    ref.  Numerical Recipes, 2nd edit. page 289

    added log(0) test 10-jan-1998 E. Kirkland
*/
double rangauss( unsigned long *iseed )
{
	double ranflat( unsigned long* );
	double x1, x2, y;
	static double tpi=2.0*3.141592654;

	/* be careful to avoid taking log(0) */
	do{
		x1 = ranflat( iseed );
		x2 = ranflat( iseed );

	} while ( (x1 < 1.0e-30) || (x2 < 1.0e-30) );

	y = sqrt( - 2.0 * log(x1) ) * cos( tpi * x2 );

	return( y );

}  /* end rangauss() */

/*--------------------- ReadfeTable() -----------------------*/
/*
   read electron scattering factors parameters 
	   from file fparam.dat
  
  the constants that must be defined above are

#define NPMAX	12	= number of parameters for each Z 
#define NZMIN	1	= min Z in featom.tab 
#define NZMAX	103	= max Z in featom.tab 
#define	NCMAX	132	= characters per line to read 

  assumed global vars:

int feTableRead=0; = flag to remember if the param file has been read 
int nl=3, ng=3; = number of Lorenzians and Gaussians 
double fparams[][] = fe parameters

*/
int ReadfeTable( )
{
   char cline[NCMAX], *cstatus;
   int i, j, n, zi, z, na;
   double w;
   
   FILE *fp;

   /* if the file has been read already then just return */
   if( feTableRead == 1 ) return(0);

   /* set next lines to point to the possible places for the
		scattering factor table "fparams.dat"
	try the current directory first and then the others
      Remember that \ is a special character - you must type \\ to
		get a \ in the path in Windows etc. */

   if( (fp = fopen( "fparams.dat", "r") ) != NULL) {}
   else if( (fp = fopen( "C:\\temsim\\fparams.dat", "r") ) != NULL) {}
   else if( (fp = fopen( "/home/kirkland/public/fparams.dat", "r") ) != NULL) {}
   else {
	perror("ReadfeTable() can't open file fparams.dat");
	exit( 0 );
   }

	na = 2*( nl + ng );	/* number of params to fit */
	fparams = double2D( NZMAX+1, na, "fparams" );
	n = 0;

	for( zi=NZMIN; zi<=NZMAX; zi++) {

	    /* find Z delimiter */
	    do { 
		cstatus = fgets( cline, NCMAX, fp );
		if( cstatus == NULL ) break;
	    } while ( strncmp( cline, "Z=", 2 ) != 0 );
	    if( cstatus == NULL ) break;

	    n += 1;
	    sscanf( cline, "Z=%d,  chisq=%lf\n", &z, &w);
	    for( j=0; j<na; j+=4 ) {
		fgets( cline, NCMAX, fp );
		for( i=0; i<4; i++) {
		    sscanf(&cline[i*17],"%le", &fparams[z][i+j] );
		}
	    }
	    if( z != zi ) {  /* test integrity of the file */
		printf( "Warning, Z= %d read when expecting"
			" Z= %d, in ReadfeTable().\n",
			z, zi);
	    }
	}  /* end for(zi=2.. */

	if( n != (NZMAX-NZMIN + 1) ) {
		printf("Warning, only %d elements read in "
			"in feTableRead() (too small).\n", n );
	}
	fclose( fp );

	feTableRead = 1;	/* remember that table has been read */
	return( n );

} /* end ReadfeTable() */

/*--------------------- ReadLine() -----------------------*/
/*
	read a full line from a file and 
	return length of line read
	
	to bad this looks like Pascal but its the easiest
	way to read just whole line because fscanf() ignores
	end of line characters
	
	fpread = pointer to file
	cMax = length of data buffer cRead
	cRead = char[] buffer to read into
	mesg = error message to print if not successful
*/
int ReadLine( FILE* fpRead, char* cRead, int cMax, const char *mesg )
{
	if( fgets( cRead, cMax, fpRead) == NULL ) {
		printf("error reading input file: %s\n", mesg);
		exit( 0 );
	}
	return( strlen( cRead ) );

}  /* end ReadLine */

/*--------------------- ReadXYZcoord() -----------------------*/
/*
	read a set of (x,y,z) coordinates from a file
	and return number of coord read
	
	infile = name of input file to read from
	ncellx,y,z = number of unit cells in x,y,z to replicate
			(must be >= 1 )
	ax, by, cz = will get unit cell size from file
	x,y,z  = pointer to a pointer to get array
			of (x,y,z) coord
	occ    = pointer to a pointer to get array
			of occupancy
	wobble = rms thermal displacement (in Angstroms)
			at Temperature = 300 degrees K
	Znum  = pointer to a pointer to get array
			atomic numbers Z
	line1 = char array to get 1st line of file
		with description
	nline1 = number of char in line1[]

   NOTE: x,y,z,occ and Znum array are allocated by this routine
	because it is not known ahead of time home many points
	will be read in (i.e. this routine figures it out)

	only infile and nline1 are unchanged by this routine

*/
int ReadXYZcoord( const char* infile, const int ncellx, const int ncelly,
	const int ncellz, float *ax, float *by, float *cz, int** Znum, 
	float** x, float** y, float** z, float** occ, float **wobble,
	char* line1, int nline1 )
{
   int i, j, iz, ncoord, na, ntotal;
   char buf[NCMAX];
   FILE *fp;

   /*  abort if ncellx,y,z not valid */
   
   	if( (ncellx<1) || (ncelly<1) || (ncellz<1) ) {
   		printf("invalid ncellx,y,z = %d, %d, %d in ReadXYZcoord()\n",
   			ncellx, ncelly, ncellz );
   		return( 0 );
   	}

   /*  make first pass to count how many coordinates there are */

	fp = fopen( infile, "r" );
	if( fp == NULL ) {
		printf("ReadXYZcoord() cannot open file %s (1)\n",
		       infile );
		exit( 0 );
	}
	/* skip first two line in first pass */
	ReadLine( fp, buf, NCMAX, "in ReadXYZcoord" );
	ReadLine( fp, buf, NCMAX, "in ReadXYZcoord" );

	ncoord = -1;
	do {	ReadLine( fp, buf, NCMAX, "in ReadXYZcoord" );
		sscanf( buf, "%d", &iz );
		ncoord += 1;
	} while( iz > 0 );
	fclose( fp );

  /* now that we know how many coordinates there are
	allocate the arrays */

	ntotal  = ncoord * ncellx * ncelly * ncellz;
	*x      = float1D( ntotal, "x" );
	*y      = float1D( ntotal, "y" );
	*z      = float1D( ntotal, "z" );
	*occ    = float1D( ntotal, "occ" );
	*wobble = float1D( ntotal, "wobble" );
	*Znum   =   int1D( ntotal, "Znum" );

	fp = fopen( infile, "r" );
	if( fp == NULL ) {
		printf("ReadXYZcoord() cannot open file %s (2)\n",
		       infile );
		exit( 0 );
	}

	/* read file description */
	ReadLine( fp, line1, NCMAX, "in ReadXYZcoord" );

	/* read unit cell size */
	*cz = 0.0F;  /* initi to 0 because cz may be blank in file */ 
	ReadLine( fp, buf, NCMAX, "in ReadXYZcoord" );
	sscanf( buf, "%g %g %g", ax, by, cz );

	for( i=0; i<ncoord; i++) {
		ReadLine( fp, buf, NCMAX, "in ReadXYZcoord()" );
		*(*wobble+i) = 0.0F;
		sscanf( buf, "%d %g %g %g %g %g",
		       *Znum+i, *x+i, *y+i, *z+i, *occ+i, *wobble+i );
		if( (*(*Znum+i) < 1 ) || (*(*Znum+i) > NZMAX) )
			printf("Warning bad atomic number %d in ReadXYZcoord()\n",
				*(*Znum+i) );
	}
	fclose( fp );
	
	na = ncoord;

        if( (ncellx > 1) || (ncelly > 1) || (ncellz > 1 ) ) {

	    if( ncellx > 1 ) {
	    	for( i=1; i<ncellx; i++)
	    	for( j=0; j<na; j++) {
	    		*(*x + j + na*i)  = *(*x + j) + i*(*ax);
	    		*(*y + j + na*i)  = *(*y + j);
	    		*(*z + j + na*i)  = *(*z + j);
	    		*(*occ + j + na*i)     = *(*occ + j);
	    		*(*wobble + j + na*i)  = *(*wobble + j);
	    		*(*Znum + j + na*i)    = *(*Znum + j);
	    	}
	    	na = na * ncellx;
	        *ax = (*ax) * ncellx;
	    }

	    if( ncelly > 1 ) {
	    	for( i=1; i<ncelly; i++)
	    	for( j=0; j<na; j++) {
	    		*(*x + j + na*i)  = *(*x + j);
	    		*(*y + j + na*i)  = *(*y + j) + i*(*by);
	    		*(*z + j + na*i)  = *(*z + j);
	    		*(*occ + j + na*i)     = *(*occ + j);
	    		*(*wobble + j + na*i)  = *(*wobble + j);
	    		*(*Znum + j + na*i)    = *(*Znum + j);
	    	}
	    	na = na * ncelly;
	        *by = (*by) * ncelly;
	    }

	    if( ncellz > 1 ) {
	    	for( i=1; i<ncellz; i++)
	    	for( j=0; j<na; j++) {
	    		*(*x + j + na*i)  = *(*x + j);
	    		*(*y + j + na*i)  = *(*y + j);
	    		*(*z + j + na*i)  = *(*z + j) + i*(*cz);
	    		*(*occ + j + na*i)     = *(*occ + j);
	    		*(*wobble + j + na*i)  = *(*wobble + j);
	    		*(*Znum + j + na*i)    = *(*Znum + j);
	    	}
	    	na = na * ncellz;
	        *cz = (*cz) * ncellz;
	    }

	}  /* end if( ncellx,y,z > 1 ) */

	return( na );

}  /* end ReadXYZcoord() */

/*----------------------- seval() ----------------------*/
/*
	Interpolate from cubic spline coefficients

	E. Kirkland 4-JUL-85
	modified to do a binary search for efficiency 13-Oct-1994 ejk
	converted to C 26-jun-1995 ejk
	fixed problem on end-of-range 16-July-1995 ejk

	The inputs are:
		x[n] = array of x values in ascending order, each x[i] must
			be unique
		y[n] = array of y values corresponding to x[n]
		b[n] = array of spline coeficients for (x-x[i])
		c[n] = array of spline coeficients for (x-x[i])**2
		d[n] = array of spline coeficients for (x-x[i])**3
		n  = number of data points
		x0  = the x value to interpolate at
		(x[i] <= x <= x[i+1]) and all inputs remain unchanged

	The value returned is the interpolated y value.

	The coeficients b[i], c[i], d[i] refer to the x[i] to x[i+1]
	interval. NOTE that the last set of coefficients,
	b[n-1], c[n-1], d[n-1] are meaningless.
*/
double seval( double *x, double *y, double *b, double *c,
	     double *d, int n, double x0 )
{
	int i, j, k;
	double z, seval1;

	/*  exit if x0 is outside the spline range */
	if( x0 <= x[0] ) i = 0;
	else if( x0 >= x[n-2] ) i = n-2;
	else { 
		i = 0;
		j = n;
		do{ k = ( i + j ) / 2 ;
			if( x0 < x[k] )  j = k;
			else if( x0 >= x[k] ) i = k;
		} while ( (j-i) > 1 );
	}
        
	z = x0 - x[i];
	seval1 = y[i] + ( b[i] + ( c[i] + d[i] *z ) *z) *z;

	return( seval1 );

} /* end seval() */

/*--------------------- sigma() -----------------------------------*/
/*
	return the interaction parameter sigma in radians/(kv-Angstroms)
	keep this is one place so I don't have to keep typing in these
	constants (that I can never remember anyhow)

	ref: Physics Vade Mecum, 2nd edit, edit. H. L. Anderson
		(The American Institute of Physics, New York) 1989
		page 4.

	kev = electron energy in keV

*/

double sigma( double kev )
{
	double s, pi, wavl, x;
	const double emass=510.99906; /* electron rest mass in keV */
	double wavelength( double kev );  /*  get electron wavelength */

	x = ( emass + kev ) / ( 2.0*emass + kev);
	wavl = wavelength( kev );
	pi = 4.0 * atan( 1.0 );
	
	s = 2.0 * pi * x / (wavl*kev);
	
	return( s );

}  /* end sigma() */


/*------------------ splinh() -----------------------------*/
/*
	fit a quasi-Hermite  cubic spline
	
	[1] Spline fit as in H.Akima, J. ACM 17(1970)p.589-602
		'A New Method of Interpolation and Smooth
		Curve Fitting Based on Local Procedures'

	[2] H.Akima, Comm. ACM, 15(1972)p.914-918

	E. Kirkland 4-JUL-85
	changed zero test to be a small nonzero number 8-jul-85 ejk
	converted to C 24-jun-1995 ejk

	The inputs are:
		x[n] = array of x values in ascending order, each X(I) must
			be unique
		y[n] = array of y values corresponding to X(N)
		n  = number of data points must be 2 or greater

	The outputs are (with z=x-x(i)):
		b[n] = array of spline coeficients for (x-x[i])
		c[n] = array of spline coeficients for (x-x[i])**2
		d[n] = array of spline coeficients for (x-x[i])**3
		( x[i] <= x <= x[i+1] )
	To interpolate y(x) = yi + bi*z + c*z*z + d*z*z*z

	The coeficients b[i], c[i], d[i] refer to the x[i] to x[i+1]
	interval. NOTE that the last set of coefficients,
	b[n-1], c[n-1], d[n-1] are meaningless.
*/
void splinh( double x[], double y[],
	     double b[], double c[], double d[], int n)
{
#define SMALL 1.0e-25

	int i, nm1, nm4;
	double m1, m2, m3, m4, m5, t1, t2, m54, m43, m32, m21, x43;

	if( n < 4) return;

	/* Do the first end point (special case),
	   and get starting values */

	m5 = ( y[3] - y[2] ) / ( x[3] - x[2] );	/* mx = slope at pt x */
	m4 = ( y[2] - y[1] ) / ( x[2] - x[1] );
	m3 = ( y[1] - y[0] ) / ( x[1] - x[0] );

	m2 = m3 + m3 - m4;	/* eq. (9) of reference [1] */
	m1 = m2 + m2 - m3;

	m54 = fabs( m5 - m4);
	m43 = fabs( m4 - m3);
	m32 = fabs( m3 - m2);
	m21 = fabs( m2 - m1);

	if ( (m43+m21) > SMALL )
		t1 = ( m43*m2 + m21*m3 ) / ( m43 + m21 );
	else
		t1 = 0.5 * ( m2 + m3 );

	/*  Do everything up to the last end points */

	nm1 = n-1;
	nm4 = n-4;

	for( i=0; i<nm1; i++) {

		if( (m54+m32) > SMALL )
			t2= (m54*m3 + m32*m4) / (m54 + m32);
		else 
			t2 = 0.5* ( m3 + m4 );
      
		x43 = x[i+1] - x[i];
		b[i] = t1;
		c[i] = ( 3.0*m3 - t1 - t1 - t2 ) /x43;
		d[i] = ( t1 + t2 - m3 - m3 ) / ( x43*x43 );

		m1 = m2;
		m2 = m3;
		m3 = m4;
		m4 = m5;
		if( i < nm4 ) {
			m5 = ( y[i+4] - y[i+3] ) / ( x[i+4] - x[i+3] );
		} else {
			m5 = m4 + m4 - m3;
		}

		m21 = m32;
		m32 = m43;
		m43 = m54;
		m54 = fabs( m5 - m4 );
		t1 = t2;
	}

	return;

} /* end splinh() */

/*------------------------ transmit() ------------------------*/
/*
	transmit the wavefunction thru one layer

	waver,i[ix][iy]  = real and imaginary parts of wavefunction
	transr,i[ix][iy] = real and imag parts of transmission functions

	nx, ny = size of array
	
	on entrance waver,i and transr,i are in real space
	
	only waver,i will be changed by this routine
*/
void transmit( float** waver, float** wavei,
			   float** transr, float** transi,
				int nx, int ny )
{
	int ix, iy;
	float wr, wi, tr, ti;
	
	   for( ix=0; ix<nx; ix++) {
	        for( iy=0; iy<ny; iy++) {
		   wr = waver[ix][iy];
		   wi = wavei[ix][iy];
		   tr = transr[ix][iy];
		   ti = transi[ix][iy];
		   waver[ix][iy] = wr*tr - wi*ti;
		   wavei[ix][iy] = wr*ti + wi*tr;
		} /* end for(iy...) */
	   }  /* end for(ix...) */

} /* end transmit() */

/*--------------------- vatom() -----------------------------------*/
/*
	return the real space atomic potential (NOT projected)
	in volts for atomic number Z at radius r

	Z = atomic number 1 <= Z <= 103
	radius  = radius in Angstroms (MUST be > 0)

  assumed global vars:

#define NZMIN	1	= min Z in featom.tab 
#define NZMAX	103	= max Z in featom.tab 

int feTableRead=0; = flag to remember if the param file has been read 
int nl=3, ng=3; = number of Lorenzians and Gaussians 
double fparams[][] = fe parameters

  al and ag calculated using physical constants from:
	H. L. Anderson, editor "A Physicist's Desk Reference",
		2nd edition, Amer. Instit. Physics, 1989

  started from vzatom() 24-nov-1997 ejk
*/

double vatom( int Z, double radius )
{
   int i, nfe;
   double suml, sumg, x,t, r;
   int ReadfeTable( );

   /* Lorenzian, Gaussian constants */
   const double al=150.4121417, ag=266.5985798;
   const double pi=3.141592654;

   if( (Z<NZMIN) || (Z>NZMAX) ) return( 0.0 );

   /* read in the table from a file if this is the
	first time this is called */
   if( feTableRead == 0 ) nfe = ReadfeTable();

   r = fabs( radius );
   if( r < 1.0e-10 ) r = 1.0e-10;  /* avoid singularity at r=0 */
   suml = sumg = 0.0;

   /* Lorenztians */
   x = 2.0*pi*r;
   for( i=0; i<2*nl; i+=2 )
	suml += fparams[Z][i]* exp( -x*sqrt(fparams[Z][i+1]) );

   /* Gaussians */
   x = pi*r;
   x = x*x;
   for( i=2*nl; i<2*(nl+ng); i+=2 ) {
	t = sqrt( fparams[Z][i+1] );
	t = t*t*t;
	sumg += fparams[Z][i]*exp(-x/fparams[Z][i+1]) / t;
   }

   return( al*suml/r + ag*sumg );

}  /* end vatom() */

/*--------------------- vzatom() -----------------------------------*/
/*
	return the real space projected atomic potential
	in volt-Angstroms for atomic number Z at radius r

	Z = atomic number 1 <= Z <= 103
	radius  = radius in Angstroms (MUST be > 0)

  assumed global vars:

#define NZMIN	1	= min Z in featom.tab 
#define NZMAX	103	= max Z in featom.tab

int feTableRead=0; = flag to remember if the param file has been read 
int nl=3, ng=3; = number of Lorenzians and Gaussians 
double fparams[][] = fe parameters

  al and ag calculated using physical constants from:
	H. L. Anderson, editor "A Physicist's Desk Reference",
		2nd edition, Amer. Instit. Physics, 1989
*/

double vzatom( int Z, double radius )
{
   int i, nfe;
   double suml, sumg, x, r;
   int ReadfeTable( );

   /* Lorenzian, Gaussian constants */
   const double al=300.8242834, ag=150.4121417;
   const double pi=3.141592654;

   if( (Z<NZMIN) || (Z>NZMAX) ) return( 0.0 );

   /* read in the table from a file if this is the
	first time this is called */
   if( feTableRead == 0 ) nfe = ReadfeTable();

   r = fabs( radius );
   if( r < 1.0e-10 ) r = 1.0e-10;  /* avoid singularity at r=0 */
   suml = sumg = 0.0;

   /* Lorenztians */
   x = 2.0*pi*r;
   for( i=0; i<2*nl; i+=2 )
	suml += fparams[Z][i]* bessk0( x*sqrt(fparams[Z][i+1]) );

   /* Gaussians */
   x = pi*r;
   x = x*x;
   for( i=2*nl; i<2*(nl+ng); i+=2 )
	sumg += fparams[Z][i]*exp(-x/fparams[Z][i+1]) / fparams[Z][i+1];

   return( al*suml + ag*sumg );

}  /* end vzatom() */

/*--------------------- wavelength() -----------------------------------*/
/*
	return the electron wavelength (in Angstroms)
	keep this is one place so I don't have to keep typing in these
	constants (that I can never remember anyhow)

	ref: Physics Vade Mecum, 2nd edit, edit. H. L. Anderson
		(The American Institute of Physics, New York) 1989
		page 4.

	kev = electron energy in keV

*/

double wavelength( double kev )
{
	double w;
	const double emass=510.99906; /* electron rest mass in keV */
	const double hc=12.3984244; /* Planck's const x speed of light*/

	/* electron wavelength in Angstroms */
	w = hc/sqrt( kev * ( 2*emass + kev ) );

	return( w );

}  /* end wavelength() */