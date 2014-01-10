/*          *** add_poisson.c ***


Adapted from the probe.c program:
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

	ANSI-C version

	Calculate a focused probe wavefunction in real space

	this file is formatted for  a tab size of 8 characters

	rewritten in C 6-dec-1995 ejk
	fixed sign error in aberration function 1-mar-1997 ejk
	removed commas from keyboard input 3-oct-1997 ejk
	
	Calculate a poisson distribution for the electron wave front
	rjhall 2008.
*/

#include <stdio.h>	/*  ANSI-C libraries */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "fft2dc.h"	/*  FFT's */
#include "memory.h"	/* memory allocation routines */
#include "slicelib.h"	/* define parameter offsets */
#include "tiffsubs.h"	/* file I/O libraries */

#define	NCMAX	132	/* characters per line to read featom.tab */
#define NPARAM	64	/* number of parameters */

void main()
{
	char fileout[NCMAX], filein[NCMAX];
	char datetime[20];

	int ix, iy, nx, ny, ixmid, iymid, i, ismoth,npix,nbits[3],samples;
	float rmin, rmax, aimin, aimax, wi, wr, phi, chi;
	float *param, **pixr, **pixi;
	double  epix, ax, by, rx2, ry2, pi, dx, dy, scale, pixel,
		Cs, df, chi1, chi2, sum,  cputime, ex;
	long nxl, nyl;

	srand( (unsigned int) time( NULL ));

/*  Echo version date  */
	printf( "c-poisson version dated 2008 rjh\n");
/*  Get desired image size, etc. */
	printf("Name of file to add poisson noise:\n");
	scanf("%s", filein );

/* Check input file */
	 if( topenFloat( filein ) != 1 ) {
                printf( "Cannot open input file %s\n", filein );
                exit( 0 );
        }
        tsize( &nxl, &nyl, nbits, &samples );
        nx = (int) nxl;
        nx = nx/2;              /* should be complex */
        ny = (int) nyl;
        ixmid = nx/2;
        iymid = ny/2;

        if( (nx != powerof2(nx) ) || ( ny != powerof2(ny) )) {
                printf("Nx, Ny = %d, %d\n", nx,ny );
                printf("   not a power of 2, try again.\n");
                tclose();
                exit( 0 );
        }
///////////////////////////////

	printf("Name of output file:\n");
	scanf("%s", fileout );

	printf("Exposure, electrons per square angstrom:\n");
	scanf("%lf", &ex );

/*  Calculate misc constants  */

	cputime = cputim( );
        pi = 4.0 * atan( 1.0 );
// Electrons per pixel	
	param = float1D( NPARAM, "param" );
        for( ix=0; ix<NPARAM; ix++) param[ix] = 0.0F;
        pixr = float2D( 2*nx, ny, "pixr" );
        pixi = pixr + nx;

        if( treadFloatPix( pixr, nxl, nyl, &npix, datetime, param )
                 != 1 ) {
                printf("Cannot read input file %s.\n", filein);
                exit(0);
        }
        tclose();
        if( npix != 2 ) {
                printf("Input file %s must be complex, can't continue.\n",
                        filein );
                exit( 0 );
        }

        ax = param[pDX] * ((float) nx);
//      printf("ax %g\n",ax);
        by = param[pDY] * ((float) ny);
//
	

	printf("ax = %f by = %f\n",ax,by);


	double pix = nx * ny;
	double ang = ax * by;
 	epix = (ang/pix) * ex;

	printf("ny = %d: nx = %d: electron per pixel = %g\n",nx, ny, epix);	
        
	rx2 = 1.0/ax;
        rx2 = rx2 * rx2;
        ry2 = 1.0/by;
        ry2 = ry2 * ry2;
//        
        ixmid = nx/2;
        iymid = ny/2;
       
//   Calculate noise 

	pixel = ( rx2 + ry2 );

	for( iy=0; iy<ny; iy++) {
	   for( ix=0; ix<nx; ix++) {
			 wi = pixi[ix][iy];
                         wr = pixr[ix][iy];
                         phi = atan2(wi,wr);
                         chi = (wi*wi)+(wr*wr);
                         chi = sqrt(chi);
			 chi = chi *  PoissonRandomNumber(epix);
                         wr = (float) cos(phi);
			 wi = (float) sin(phi);
 			 pixi[ix][iy] = wi * chi;
                         pixr[ix][iy] = wr * chi;


		}
	}

//	fft2d( pixr, pixi, nx, ny, -1);

//:  Output results and find min and max to echo
//    remember that complex pix are stored in the file in FORTRAN
//		order for compatability


        rmin = pixr[0][0];
        rmax = rmin;
        aimin = pixi[0][0];
        aimax = aimin;
        for( iy=0; iy<ny; iy++) {
          for( ix=0; ix<nx; ix++) {
		if( pixr[ix][iy] < rmin ) rmin = pixr[ix][iy];
		if( pixr[ix][iy] > rmax ) rmax = pixr[ix][iy];
		if( pixi[ix][iy] < aimin ) aimin = pixi[ix][iy];
		if( pixi[ix][iy] > aimax ) aimax = pixi[ix][iy];
	   }
	}

	param[pRMAX] = rmax;
       	param[pIMAX] = aimax;
       	param[pRMIN] = rmin;
        param[pIMIN] = aimin;

	if( tcreateFloatPixFile( fileout, pixr, (long) (2*nx), (long) ny,
			 2, param ) != 1 )
		printf( "cannot write an output file.\n");

	printf( "Pix range %15.7g to %15.7g real,\n"
		"      and %15.7g to %15.7g imaginary\n",
		rmin, rmax, aimin, aimax );
	cputime = cputim() - cputime;
	printf("\nCPU time = %f sec\n", cputime );

}   

int PoissonRandomNumber(double lambda)
{
	
  	double p  = ((double)rand()/((double)(RAND_MAX)+(double)(1))); //uniform random number
  	int k=0;                          //Counter
  	int max_k = 1000;                 //k upper limit
  	double P = exp(-lambda);          //probability
  	double sum=P;                     //cumulant
	if (sum>=p) return 0;             //done allready
	for (k=1; k<max_k; ++k) {         //Loop over all k:s
    		P*=lambda/(double)k;            //Calc next prob
    		sum+=P;                         //Increase cumulant
    		if (sum>=p) break;              //Leave loop
  	}
  	return k;                         //return random number
}

