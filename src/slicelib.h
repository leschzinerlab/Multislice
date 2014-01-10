/*		*** slicelib.h ***

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

/* define symbols for the parameter offsets */

#define pNPIX		0   /* number of pix 1 for real and 2 for complex */
#define	pRMAX		1	/* maximum value of the real part of the image */
#define	pIMAX		2	/* maximum value of the imaginary part of the image */
#define pRMIN		3	/* minimum value of the real part of the image */
#define	pIMIN		4	/* minimum value of the imaginary part of the image */
#define pXBTILT		5	/* x beam tilt in rad */
#define pYBTILT		6	/* y beam tilt in rad */
#define pC			7	/* c unit cell dimension in Angstroms */
#define pRES		8	/* real space resolution in atompot */
#define pXCTILT		9	/* x crystal tilt in rad */
#define pYCTILT		10	/* y crystal tilt in rad */
#define pDEFOCUS	11	/* defocus in Angstroms */
#define pASTIG		12	/* astigmatism in Angstroms */
#define pTHETA		13	/* angle of astigmatism in radians */
#define pDX			14	/* dimension of pixel in x direction in Angstroms */
#define pDY			15	/* dimension of pixel in y direction in Angstroms */
#define pENERGY		16	/* beam energy in keV */
#define pOAPERT		17	/* objective aperture semi-angle in radians */
#define pCS			18	/* spherical aberration in Angstroms */
#define pWAVEL		19	/* electron wavelength in Angstroms */
#define pCAPERT		21	/* condenser (CTEM) illumination angle in radians */
#define pDDF		22	/* defocus spread in Angstroms */
#define pNSLICES	29	/* number of slices */
#define pMINDET		31	/* minimum detector angle (STEM) in radians */
#define pMAXDET		32	/* maximum detector angle (STEM) in radians */

/*--------------------- askYN() -----------------------------------*/
/*
	ask a yes/no question and return 1 or 0 for TRUE/FALSE

	message[] = question to ask
*/
int askYN( const char message[] );

/*-------------------- bessi0() ---------------*/
/*
    modified Bessel function I0(x)
    see Abramowitz and Stegun page 379

    x = (double) real arguments

    12-feb-1997 E. Kirkland
 */
 double bessi0( double x );

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
 double bessk0( double x );

/*--------------------- cputim() -----------------------------------*/
/*
   retrieve current CPU time in seconds
*/

double cputim();

/*--------------------- featom() -----------------------------------*/
/*
	return the electron scattering factor for atomic
	number Z at scattering angle k

	Z = atomic number 2 <= Z <= 103
	k2  = k*k where k =1/d = scattering angle (in 1/A)

  assumed global vars:

int feTableRead=0; = flag to remember if the param file has been read 
int nl=3, ng=3; = number of Lorenzians and Gaussians 
double fparams[][] = fe parameters

*/

double featom( int iz, double s );

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
void freqn( float *ko, float *ko2, float *xo, int nk, double ak );

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

   c             = input character string
*  islice(nsmax) = integer array to get stacking sequence indicies
   nsmax         = (integer) size of layer
   lmax          = (integer) maximum allowed layer index
*  nslice        = (integer) number of layers
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

   a * means that these variables may be modified by this routinec

  cname determines the mapping of characters into numbers.

*/

int parlay( const char c[], int islice[], int nsmax, int lmax, int *nslice,
				 int fperr );
				 
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
	float* kx2, float* ky2, float k2max, int nx, int ny );


/*---------------------------- ranflat -------------------------------*/
/*
	return a random number in the range 0.0->1.0
	with uniform distribution

	the 'Magic Numbers' are from 
		Numerical Recipes 2nd edit pg. 285
*/
double ranflat( unsigned long *iseed );

/*-------------------- rangauss() -------------------------- */
/*
	Return a normally distributed random number with 
	zero mean and unit variance using Box-Muller method
	
	ranflat() is the source of uniform deviates

    ref.  Numerical Recipes, 2nd edit. page 289
*/
double rangauss( unsigned long *iseed );

/*--------------------- ReadfeTable() -----------------------*/
/*
   read electron scattering factors from file 
	featom.tab from the following:

  [1] Table 2.4.6A of The Internation Tables for X-Ray 
       Crystallography IV (1974) for 0<s<2.0, all Z
  [2] Doyle and Turner Acta Cryst A24 (1968) p390-397,
       for 2.5<s<6.0, Z= 2-38, 42, 47-51, 53-56, 63, 
       79,80, 82,83, 86, 92 
  [3] Fox, O'Keefe and Tabbernor, Acta Cryst A45(1989)p.786-793
       for Z= 39-41, 43-46, 52, 57-62, 64-78, 81, 84,85, 87-91,
       93-98 and s=2.5-6.0 using Mott formula applied to Table 1
       (generated using file fot.tab and program fotest.f)

  -- remember that there is a /'+++S/Z'/ in the comments so look
        for '+++S/Z ' instead (i.e. add a space at the end
  -- 1<Z<96 is in same size blocks - read Z=97,98 separately

  s[][]   = float** to get scattering angles NSTMAX x NSSMAX
  fet[][] = float** to get electron scattering factors NSTMAX x NZMAX
  
  the constants that must be defined above are

#define NSTMAX	62	 max s lines in featom.tab 
#define NSSMAX	13	 max groups in featom.tab 
#define NZMAX	98	 max Z in featom.tab 
#define	NCMAX	132	 characters per line to read

*/
int ReadfeTable( );

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
int ReadLine( FILE* fpRead, char* cRead, int cMax, const char *mesg );
 
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
	char* line1, int nline1 );

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
	     double *d, int n, double x0 );

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

double sigma( double kev );

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
	     double b[], double c[], double d[], int n);

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
				int nx, int ny );

/*--------------------- vatom() -----------------------------------*/
/*
	return the real space atomic potential (NOT projected)
	in volts for atomic number Z at radius r

	Z = atomic number 2 <= Z <= 103
	radius  = radius in Angstroms (MUST be > 0)

  assumed global vars:

int feTableRead=0; = flag to remember if the param file has been read 
int nl=3, ng=3; = number of Lorenzians and Gaussians 
double fparams[][] = fe parameters

  al and ag calculated using physical constants from:
	H. L. Anderson, editor "A Physicist's Desk Reference",
		2nd edition, Amer. Instit. Physics, 1989

  started from vzatom() 24-nov-1997 ejk
*/

double vatom( int Z, double radius );

/*--------------------- vzatom() -----------------------------------*/
/*
	return the real space projected atomic potential
	in volt-Angstroms for atomic number Z at radius r

	Z = atomic number 2 <= Z <= 103
	radius  = radius in Angstroms (MUST be > 0)

  assumed global vars:

int feTableRead=0; = flag to remember if the param file has been read 
int nl=3, ng=3; = number of Lorenzians and Gaussians 
double fparams[][] = fe parameters

  al and ag calculated using physical constants from:
	H. L. Anderson, editor "A Physicist's Desk Reference",
		2nd edition, Amer. Instit. Physics, 1989
*/

double vzatom( int Z, double radius );

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

double wavelength( double kev );

