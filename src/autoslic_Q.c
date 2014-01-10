/*
          *** autoslice.c ***

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

  ANSI C and TIFF version

  Transmit an electron wave through a specimen using the
  multislce method with automatic slicing.  Read in the (x,y,z)
  coordinates of the whole specimen and break into slices
  on-the-fly.

  started 24-july-1996 E. Kirkland
  working 19feb-1997 ejk
  last revised 19-feb-1997 ejk
  added look-up-table vzatomLUT() for 3X-4X increase 
		in speed 23-may-1997 ejk
  put bandwith limit inside trlayer() 1-oct-1997 ejk
  added Gaussian thermal displacements 1-oct-1997 ejk
  removed /sqrt(3) from Thermal rms displacements 
	to be consistent with Int'l X-ray tables 22-dec-1997 ejk
  corrected zmin/max error with thermal displac. 24-dec-1997 ejk
  fixed small aliasing problem 5-jan-1998 ejk
  added unit cell replication option and moved ReadXYZcoord()
  	into slicelib.c  11-jan-1998 ejk
  added astigmatism and modify to use different set of
    random offsets on each illum. angle with partial coherence
             5-feb-1998 ejk
  fix typo in z range message with partial coherence and
	thermal vibrations 9-july-1998 ejk
  
  Many fixes and updates 2008 rjh

  ax,by,cz  = unit cell size in x,y
  BW     = Antialiasing bandwidth limit factor
  acmin  = minimum illumination angle
  acmax  = maximum illumination angle
  Cs     = spherical aberration coeficient
  df0    = defocus (mean value)
  sgmaf = defocus spread (standard deviation)
  dfdelt = sampling interval for defocus integration
  
  this file is formatted for a TAB size of 8 characters 
  
*/

#include <stdio.h>	/* ANSI C libraries */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "fft2dc.h"	/* FFT routines */
#include "memory.h"	/* memory allocation routines */
#include "slicelib.h"	/* misc. routines for multislice */
#include "tiffsubs.h"	/* file I/O routines in TIFF format */

#define BW (2.0F/3.0F)	/* bandwidth limit */
#define ABERR 1.0e-4	/* max error for a,b */

#define NSMAX 1000	/* max number of slices */
#define NPARAM	64	/* number of parameters */

#define NCMAX 256	/* max characters in file names */

#define NZMAX	103	/* max atomic number Z */
#define NRMAX	100	/* number of in look-up-table in vzatomLUT */
#define RMIN	0.01	/* r (in Ang) range of LUT for vzatomLUT() */
#define RMAX	5.0

/* spline interpolation coeff. */
int splineInit=0, *nspline;
double  *splinx, **spliny, **splinb, **splinc, **splind;


/* define functions at end of this file (i.e. so main can be 1st) */

void trlayer( const float x[], const float y[], const float occ[],
	const int Znum[], const int natom, 
	const float ax, const float by, const float kev,
	float **transr, float **transi, const long nx, const long ny,
	double *phirms, long *nbeams, const float k2max );
double vzatomLUT( int Z, double r );
void sortByZ( float x[], float y[], float z[], float occ[],
	int Znum[], int natom );
void QsortByZ( float x[], float y[], float z[], float occ[],
        int Znum[], int beg, int end );
/* global data for trlayer() */
float *kx2, *ky2;

void main()
{
	char filein[NCMAX], fileout[NCMAX], filestart[NCMAX],
		filebeam[NCMAX], datetime[20], description[NCMAX];

	int lstart=0, lpartl=0, lbeams=0, lwobble=0;
	int ix, iy, nx, ny, ixmid, iymid, i, nslic0, islice,
		nacx,nacy, iqx, iqy, npix,
		ndf, idf, nbout, ib, nbits[3], samples,
		ncellx, ncelly, ncellz;
	int *hbeam, *kbeam;
	int natom, *Znum, *Znum2, istart, na;
	long nxl, nyl, nbeams, nillum;
	long  ltime;
	unsigned long iseed;

	float *x, *y, *z, *occ, *wobble;
 	float *x2, *y2, *z2, *occ2;
   
	float wmin, wmax, xmin,xmax, ymin, ymax, zmin, zmax;

	float *kx, *ky, *xpos, *ypos, *param, *sparam;
	float k2, k2max, scale, v0, mm0, wavlen, rx, ry,
		ax, by, cz, pi, rmin, rmax, aimin, aimax,
		rx2,ry2, ctiltx, ctilty, tctx, tcty,
		acmin, acmax, Cs, df, df0, sigmaf, dfdelt, aobj,
		qx, qy, qy2, q2, q2min, q2max, sumdf, pdf, k2maxo,
		temperature, dfa2, dfa2phi, dfa3, dfa3phi;

	float tr, ti, wr, wi;

	float **waver, **wavei, **transr, **transi,
		*propxr, *propxi, *propyr, *propyi,
		**tempr, **tempi, **pix;

	double sum, timer, xdf, chi, chi1, chi2, phi, t, zslice,
		deltaz, phirms;
		
	FILE *fp1;

/*  echo version date and get input file name */

	printf("autoslice_Q version dated 2008 rjh\n\n");
	pi = (float) (4.0 * atan( 1.0 ));
	param = float1D( NPARAM, "param" );
	sparam = float1D( NPARAM, "sparam" );

	printf("Name of file with input atomic "
	       "potential in x,y,z format:\n");
	scanf("%s", filein );

/*  get simulation options */

        printf("Replicate unit cell by NCELLX,NCELLY,NCELLZ :\n");
	scanf("%d %d %d", &ncellx, &ncelly, &ncellz);
	if( ncellx < 1 ) ncellx = 1;
	if( ncelly < 1 ) ncelly = 1;
	if( ncellz < 1 ) ncellz = 1;

	printf("Name of file to get binary output of multislice result:\n");
	scanf("%s", fileout );

	lpartl = askYN("Do you want to include partial coherence");

	if( lpartl == 1 ) {
		printf("Illumination angle min, max in mrad:\n");
		scanf("%f %f", &acmin, &acmax);
		acmin  = acmin  * 0.001F;
		acmax  = acmax  * 0.001F;
		printf("Spherical aberration (in mm.):\n");
		scanf("%g", &Cs);
		Cs = Cs * 1.0e7F;
		printf("Defocus, mean, standard deviation, and"
		       " sampling size (in Angstroms) =\n");
		scanf("%f %f %f", &df0, &sigmaf, &dfdelt);
		printf("Objective aperture (in mrad) =\n");
		scanf("%f", &aobj);
		aobj = aobj * 0.001F;
		printf( "Magnitude and angle of 2-fold astigmatism"
			" (in Ang. and degrees):\n");
		scanf( "%f %f", &dfa2, &dfa2phi);
			dfa2phi = dfa2phi * pi /180.0F;
		printf( "Magnitude and angle of 3-fold astigmatism"
			" (in Ang. and degrees):\n");
		scanf( "%f %f", &dfa3, &dfa3phi);
		dfa3phi = dfa3phi * pi /180.0F;
		lstart = 0;
	} else {
		printf("NOTE, the program image must also be run.\n");
		lstart = askYN("Do you want to start from previous result");
	}

	if ( lstart == 1 ) {
		printf("Name of file to start from:\n");
		scanf("%s", filestart);
	} else {
		printf("Incident beam energy in kev:\n");
		scanf("%g", &v0);
		printf("Wavefunction size in pixels, Nx,Ny:\n");
		scanf("%d %d", &nx, &ny );
	}

	printf("Crystal tilt x,y in mrad.:\n");
	scanf("%f %f", &ctiltx, &ctilty);
	ctiltx = ctiltx /1000;
	ctilty = ctilty /1000;

	/*  remember that the slice thickness must be > atom size
		to use projected atomic potential */
	printf("Slice thickness (in Angstroms):\n");
	scanf("%lf", &deltaz );
	if( deltaz < 1.0 ) {
		printf("WARNING: this slice thickness is probably too thin"
			" for autoslice to work properly.\n");
	}

	if( lpartl == 0 ) {
		lbeams = askYN("Do you want to record the (real,imag) value\n"
			" of selected beams vs. thickness");
		if( lbeams == 1 ) {
			printf("Name of file for beams info:\n");
			scanf("%s", filebeam );
			printf("Number of beams:\n");
			scanf("%d", &nbout);
			if( nbout<1 ) nbout = 1;
			hbeam = int1D( nbout, "hbeam" );
			kbeam = int1D( nbout, "kbeam" );
			for( ib=0; ib<nbout; ib++) {
				printf("Beam %d, h,k=\n", ib+1);
				scanf("%d %d", &hbeam[ib], &kbeam[ib] );
			}
		}
	}

	lwobble = askYN("Do you want to include thermal vibrations");
	if( lwobble == 1 ) {
	    printf( "Type the temperature in degrees K:\n");
	    scanf( "%g", &temperature );
	    /* get random number seed from time if available 
			otherwise ask for a seed */
	    if( lwobble == 1 ) {
		ltime = (long) time( NULL );
		iseed = (unsigned) ltime;
		if( ltime == -1 ) {
			printf("Type initial seed for random number generator:\n");
			scanf("%ld", &iseed);
		} else
			printf( "Random number seed initialized to %ld\n", iseed );
		}
	} else temperature = 0.0F;

/* start timing the actual computation just for fun */

	timer = cputim();

/*  get starting value of transmitted wavefunction if required
   (this can only be used in coherent mode)
    remember to save params for final output pix  */

	if ( lstart == 1 ) {
		if( topenFloat( filestart ) != 1 ) {
			printf("Cannot open input file: %s .\n", filestart ); 
			exit( 0 );
	 	}
	 	tsize( &nxl, &nyl, nbits, &samples );
	 	nx = (int) nxl;
	 	ny = (int) nyl;

		waver = float2D( nx, ny, "waver" );
		nx = nx/2;
		wavei = waver + nx;
		if( treadFloatPix( waver, nxl, nyl, &npix, datetime, sparam )
			 != 1 ) {
			printf("Cannot read input file %s.\n", filestart);
			exit(0);
		}
		tclose();
		if( npix != 2 ) {
			printf("Input file %s must be complex, can't continue.\n",
				filestart );
			exit( 0 );
		}
           	if( (nx != powerof2(nx)) || (ny != powerof2(ny)) ) {
			printf("Nx=%d, Ny=%d must be a power of 2\n"
			   "problem in starting image, try again.\n", nx, ny);
			exit( 0 );
		}

		ax = sparam[pDX] * nx;
		by = sparam[pDY] * ny;
		v0     = sparam[pENERGY];
		nslic0 = (int) sparam[pNSLICES];
		printf("Starting pix range %g to %g real\n"
		       "                   %g to %g imag\n",
		       sparam[pRMIN], sparam[pRMAX], sparam[pIMIN], sparam[pIMAX] );
		printf("Beam voltage = %g kV\n", v0);
		printf("Old crystal tilt x,y = %g, %g mrad\n", 1000.*sparam[pXCTILT],
			1000.*sparam[pYCTILT]);

	} else {

	  nslic0 = 0;

	}  /* end if( lstart...) */

/*  calculate relativistic factor and electron wavelength */

	mm0 = 1.0F + v0/511.0F;
	wavlen = (float) wavelength( v0 );
	printf("electron wavelength = %g Angstroms\n", wavlen);

/*  read in specimen coordinates and scattering factors */

	natom = ReadXYZcoord( filein, ncellx, ncelly, ncellz,
		&ax, &by, &cz, &Znum, &x, &y, &z, &occ, &wobble,
		description, NCMAX );

	printf("%d atomic coordinates read in\n", natom );
	printf("%s", description );

	printf("Size in pixels Nx, Ny= %d x %d = %d beams\n",
	       nx,ny, nx*ny);
	printf("Lattice constant a,b = %12.4f, %12.4f\n", ax,by);

	/*  calculate the total specimen volume and echo */
	xmin = xmax = x[0];
	ymin = ymax = y[0];
	zmin = zmax = z[0];
	wmin = wmax = wobble[0];

	for( i=0; i<natom; i++) {
		if( x[i] < xmin ) xmin = x[i];
		if( x[i] > xmax ) xmax = x[i];
		if( y[i] < ymin ) ymin = y[i];
		if( y[i] > ymax ) ymax = y[i];
		if( z[i] < zmin ) zmin = z[i];
		if( z[i] > zmax ) zmax = z[i];
		if( wobble[i] < wmin ) wmin = wobble[i];
		if( wobble[i] > wmax ) wmax = wobble[i];
	}
	printf("Total specimen range is\n %g to %g in x\n"
	       " %g to %g in y\n %g to %g in z\n", xmin, xmax,
	       ymin, ymax, zmin, zmax );
	if( lwobble == 1 )
		printf("Range of thermal rms displacements (300K) = %g to %g\n",
			wmin, wmax );

/*  calculate spatial frequencies and positions for future use */

	rx = 1.0F/ax;
	rx2= rx*rx;
	ry = 1.0F/by;
	ry2= ry*ry;
	ixmid = nx/2;
	iymid = ny/2;
	nxl = nx;
	nyl = ny;

	kx   = float1D( nx, "kx" );
	kx2  = float1D( nx, "kx2" );
	xpos = float1D( nx, "xpos" );
	freqn( kx, kx2, xpos, nx, ax );

	ky   = float1D( ny, "ky" );
	ky2  = float1D( ny, "ky2" );
	ypos = float1D( ny, "ypos" );
	freqn( ky, ky2, ypos, ny, by );

/*  allocate some more arrays and initialize wavefunction */

	transr = float2D( nx, ny, "transr" );
	transi = float2D( nx, ny, "transi" );

	if( lstart == 0 ) {
		waver = float2D( 2*nx, ny, "waver" );
		wavei = waver + nx;
		for( ix=0; ix<nx; ix++)
		for( iy=0; iy<ny; iy++) {
			waver[ix][iy] = 1.0F;
			wavei[ix][iy] = 0.0F;
		}
	}

 /*  calculate propagator function  */
 
	k2max = nx/(2.0F*ax);
	tctx = ny/(2.0F*by);
	if( tctx < k2max ) k2max = tctx;
	k2max = BW * k2max;
	printf("Bandwidth limited to a real space resolution of %f Angstroms\n",
					 1.0F/k2max);
	printf("   (= %.2f mrad)  for symmetrical anti-aliasing.\n",
		 wavlen*k2max*1000.0F);
	k2max = k2max*k2max;

	tctx = (float) (2.0 * tan(ctiltx));
	tcty = (float) (2.0 * tan(ctilty));
	
	propxr = float1D( nx, "propxr" );
	propxi = float1D( nx, "propxi" );
	propyr = float1D( ny, "propyr" );
	propyi = float1D( ny, "propyi" );

	scale = pi * ((float)deltaz);

	for( ix=0; ix<nx; ix++) {
		t = scale * ( kx2[ix]*wavlen - kx[ix]*tctx );
		propxr[ix] = (float)  cos(t);
		propxi[ix] = (float) -sin(t);
	}
	for( iy=0; iy<ny; iy++) {
		t = scale * ( ky2[iy]*wavlen - ky[iy]*tcty );
		propyr[iy] = (float)  cos(t);
		propyi[iy] = (float) -sin(t);
	}


/*  iterate the multislice algorithm proper

   NOTE: zero freg is in the bottom left corner and
     expandes into all other corners - not in the center
     this is required for the FFT - don't waste time rearanging

  partial coherence method
   force the integrals to include the origin and to be symmetric
   about the origin and to have the same periodic boundary
   conditions as the sampling grid
*/
	if( lpartl == 1 ) {

	    printf("Illumination angle sampling (in mrad) = %f, %f\n\n",
			1000.*rx*wavlen, 1000.*ry*wavlen);

	    pix = float2D( nx, ny, "pix" );
	    for( ix=0; ix<nx; ix++)
		for( iy=0; iy<ny; iy++) pix[ix][iy] = 0.0F;
		
	    tempr = float2D( nx, ny, "tempr" );
	    tempi = float2D( nx, ny, "tempi" );

	    ndf = (int) ( ( 2.5F * sigmaf ) / dfdelt );

	    nacx = (int) ( ( acmax / ( wavlen * rx ) ) + 1.5F );
	    nacy = (int) ( ( acmax / ( wavlen * ry ) ) + 1.5F );

	    q2max = acmax / wavlen;
	    q2max = q2max*q2max;

	    q2min = acmin / wavlen;
	    q2min = q2min*q2min;

	    k2maxo = aobj / wavlen;
	    k2maxo = k2maxo*k2maxo;

	    chi1 = pi * wavlen;
	    chi2 = 0.5 * Cs * wavlen *wavlen;
	    nillum = 0;

	    x2      = float1D( natom, "x2" );  /* for Monte Carlo stuff */
	    y2      = float1D( natom, "y2" );
	    z2      = float1D( natom, "z2" );
	    occ2    = float1D( natom, "occ2" );
	    Znum2   =   int1D( natom, "Znum2" );

	    if( lwobble == 0 )
		 //   sortByZ( x, y, z, occ, Znum, natom );
		QsortByZ (x, y, z, occ, Znum, 0, natom); 
/*  integrate over the illumination angles */

	    for( iqy= -nacy; iqy<=nacy; iqy++) {
		qy = iqy * ry;
		qy2 = qy * qy;

		for( iqx= -nacx; iqx<=nacx; iqx++) {
		    qx = iqx * rx;
		    q2 = qx*qx + qy2;

		    if( (q2 <= q2max) && (q2 >= q2min) ) {
			nillum += 1;
			for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++) {
			    t = 2.0*pi*( qx*xpos[ix] + qy*ypos[iy] );
			    waver[ix][iy] = (float) cos(t);
			    wavei[ix][iy] = (float) sin(t);
			}
			/*  add random thermal displacements scaled by temperature
					if requested 
				remember that initial wobble is at 300K for each direction */
			if( lwobble == 1 ){
			    scale = (float) sqrt(temperature/300.0) ;
			    for( i=0; i<natom; i++) {
				x2[i] = x[i] + (float)(wobble[i]*rangauss(&iseed)*scale);
				y2[i] = y[i] + (float)(wobble[i]*rangauss(&iseed)*scale);
				z2[i] = z[i] + (float)(wobble[i]*rangauss(&iseed)*scale);
				occ2[i] = occ[i];
				Znum2[i] = Znum[i];
			    }
			    printf( "Sorting atoms by depth...\n");
			//    sortByZ( x2, y2, z2, occ2, Znum2, natom );
			QsortByZ (x2, y2, z2, occ2, Znum2, 0, natom); 
			    zmin = z2[0];		/* reset zmin/max after wobble */
			    zmax = z2[natom-1];
			    printf("Thickness range with thermal displacements"
				" is %g to %g (in z)\n", zmin, zmax );
			} else for( i=0; i<natom; i++) {
				x2[i] = x[i];
				y2[i] = y[i];
				z2[i] = z[i];
				occ2[i] = occ[i];
				Znum2[i] = Znum[i];
			}

			zslice = zmin + deltaz;
			istart = 0;

			while( istart < natom ) {

			    /* find range of atoms for current slice */
			    na = 0;
			    for(i=istart; i<natom; i++) 
				if( z2[i] < zslice ) na++; else break;

			    /* calculate transmission function */
			    trlayer( &x2[istart], &y2[istart], &occ2[istart],
				&Znum2[istart],na, ax, by, v0, 
				transr, transi, nxl, nyl, &phirms, &nbeams, k2max );

			    transmit( waver, wavei, transr, transi, nx, ny );

			    /* remember: prop needed here to get anti-aliasing
						right */
			    fft2d( waver, wavei, nx, ny, +1);
			    propagate( waver,  wavei, propxr, propxi,
				    propyr, propyi,
				    kx2,  ky2,  k2max, nx, ny );
			    fft2d( waver, wavei, nx, ny, -1);

			    zslice += deltaz;
			    istart += na;

			} /* end while(zslice<=..) */

			scale = 1.0F / ( ((float)nx) * ((float)ny) );
			sum = 0.0;
			for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++)
				sum += waver[ix][iy]*waver[ix][iy]
					+ wavei[ix][iy]*wavei[ix][iy];
			sum = sum * scale;

			printf("Illumination angle = %7.3f, %7.3f mrad",
				1000.*qx*wavlen, 1000.*qy*wavlen);
			printf(", integrated intensity= %f\n", sum );

/*  integrate over +/- 2.5 sigma of defocus */

			fft2d( waver, wavei, nx, ny, +1);
			sumdf = 0.0F;

			for( idf= -ndf; idf<=ndf; idf++) {
			    df = df0 + idf*dfdelt;

			    for( ix=0; ix<nx; ix++)
			    for( iy=0; iy<ny; iy++) {
				k2 = kx2[ix] + ky2[iy];
				if( k2 <= k2maxo ) {
					phi = atan2( ky[iy], kx[ix] );
					chi = chi1*k2* ( chi2*k2 - df 
					   + dfa2*sin( 2.0*(phi-dfa2phi) ) 
					   + 2.0F*dfa3*wavlen*sqrt(k2)*
					   sin( 3.0*(phi-dfa3phi) )/3.0 );
					tr = (float)  cos(chi);
					ti = (float) -sin(chi);
					wr = waver[ix][iy];
					wi = wavei[ix][iy];
					tempr[ix][iy] = wr*tr - wi*ti;
					tempi[ix][iy] = wr*ti + wi*tr;
				} else {
					tempr[ix][iy] = 0.0F;
					tempi[ix][iy] = 0.0F;
				}
			    }

			    fft2d( tempr, tempi, nx, ny, -1);

			    xdf = (double) ( (df - df0) /sigmaf );
			    pdf = (float) exp( -0.5F * xdf*xdf );
			    sumdf += pdf;

			    for( ix=0; ix<nx; ix++)
			    for( iy=0; iy<ny; iy++) {
		  		wr = tempr[ix][iy];
		  		wi = tempi[ix][iy];
				pix[ix][iy] += pdf* ( wr*wr + wi*wi );
			    }

			}/* end for(idf..) */
		    }/* end if( q2...) */

  		} /* end for( iqx..) */
  	    } /* end for( iqy..) */

	    printf("Total number of illumination angle = %ld\n",
				nillum);
	    printf("Total number of defocus values = %d\n", 2*ndf+1);
	    scale = 1.0F / ( ((float)nillum) * sumdf );
	    rmin  = pix[0][0] * scale;
	    rmax  = rmin;
	    aimin = 0.0F;
	    aimax = 0.0F;

	    for( ix=0; ix<nx; ix++)
	    for( iy=0; iy<ny; iy++) {
	    	pix[ix][iy] = pix[ix][iy] * scale;
	    	if( pix[ix][iy] < rmin ) rmin = pix[ix][iy];
	   	if( pix[ix][iy] > rmax ) rmax = pix[ix][iy];
	    }

/* ---- start coherent method below ----------------
	    (remember that waver,i[][] was initialize above) */

	} else {

		if( lbeams ==1 ) {
			fp1 = fopen( filebeam, "w" );
			if( fp1==NULL) {
				printf("can't open file %s\n", filebeam);
				exit(0);
			}
			fprintf( fp1, " (h,k) = ", nbout);
			for(ib=0; ib<nbout; ib++)
				fprintf(fp1," (%d,%d)", hbeam[ib], kbeam[ib]);
			fprintf( fp1, "\n" );
			fprintf( fp1, "nslice, (real,imag) (real,imag) ...\n\n");
			for( ib=0; ib<nbout; ib++) {
				if( hbeam[ib] < 0 ) hbeam[ib] = nx + hbeam[ib];
				if( kbeam[ib] < 0 ) kbeam[ib] = ny + kbeam[ib];
				if( hbeam[ib] < 0 ) hbeam[ib] = 0;
				if( kbeam[ib] < 0 ) kbeam[ib] = 0;
				if( hbeam[ib] > nx-1 ) hbeam[ib] = nx-1;
				if( kbeam[ib] > ny-1 ) kbeam[ib] = ny-1;
			}
		}

		/*  add random thermal displacements scaled by temperature if requested 
			remember that initial wobble is at 300K for each direction */
		if( lwobble == 1 ){
			scale = (float) sqrt(temperature/300.0) ;
			for( i=0; i<natom; i++) {
				x[i] += (float) (wobble[i] * rangauss( &iseed ) * scale);
				y[i] += (float) (wobble[i] * rangauss( &iseed ) * scale);
				z[i] += (float) (wobble[i] * rangauss( &iseed ) * scale);
			}
		}

		printf( "Sorting atoms by depth...\n");
	//Richard edit
	//	for( i=0; i<natom; i++) {
	//	printf("%d: x %f, y %f, z %f, Znum %d\n", i, x[i], y[i], z[i], Znum[i]);
	//	}

	
	
			//sortByZ( x, y, z, occ, Znum, natom );
			QsortByZ (x, y, z, occ, Znum, 0, natom); 
			
			zmin = z[0];		/* reset zmin/max after wobble */
			zmax = z[natom-1];
			printf("Thickness range with thermal displacements"
				" is %g to %g (in z)\n", zmin, zmax );
	//	printf("\n");
	//	for( i=0; i<natom; i++) {
	//	printf("%d: x %f, y %f, z %f, Znum %d\n", i, x[i], y[i], z[i], Znum[i]);
	//	}
		scale = 1.0F / ( ((float)nx) * ((float)ny) );

		zslice = zmin + deltaz/2;  /* ??? */
/*???		zslice = zmin + deltaz;  */
		istart = 0;
		islice = 1;

		while( (istart < natom) && ( zslice < (zmax+deltaz) ) ) {

		    /* find range of atoms for current slice */
		    na = 0;
		    for(i=istart; i<natom; i++) 
			if( z[i] < zslice ) na++; else break;

		    /* calculate transmission function and bandwidth limit */
		    trlayer( &x[istart], &y[istart], &occ[istart],
			&Znum[istart], na, ax, by, v0, transr, transi,
			nxl, nyl, &phirms, &nbeams, k2max );

/*??? printf("average atompot comparison = %g\n", phirms/(wavlen*mm0) ); */

		    transmit( waver, wavei, transr, transi, nx, ny );

		    fft2d( waver, wavei, nxl, nyl, +1);
		    if( lbeams== 1 )  {
			fprintf( fp1, "%5d", islice);
			for( ib=0; ib<nbout; ib++) 
				fprintf(fp1, "%10.6f %10.6f",
					scale*waver[hbeam[ib]][kbeam[ib]],
					scale*wavei[hbeam[ib]][kbeam[ib]]);
				fprintf( fp1, "\n");
		    }
		    /* remember: prop needed here to get anti-aliasing right */
		    propagate( waver,  wavei, propxr, propxi,
				propyr, propyi,	kx2,  ky2,  k2max, nx, ny );
		    fft2d( waver, wavei, nxl, nyl, -1);

		    sum = 0.0;
		    for( ix=0; ix<nx; ix++)
		    for( iy=0; iy<ny; iy++)
		   	sum += waver[ix][iy]*waver[ix][iy] +
		   		wavei[ix][iy]*wavei[ix][iy];
		    sum = sum * scale;

		    printf("z= %f A, %ld beams, %d coord., \n"
		       "     aver. phase= %f, total intensity = %f\n",
		       zslice, nbeams, na, phirms, sum );

		    zslice += deltaz;
		    istart += na;
		    islice++;

		} /* end while(istart<natom..) */
	
		rmin  = waver[0][0];
		rmax  = rmin;
		aimin = wavei[0][0];
		aimax = aimin;

		for( ix=0; ix<nx; ix++)
		for( iy=0; iy<ny; iy++) {
			wr = waver[ix][iy];
			wi = wavei[ix][iy];
			if( wr < rmin ) rmin = wr;
			if( wr > rmax ) rmax = wr;
			if( wi < aimin ) aimin = wi;
			if( wi > aimax ) aimax = wi;
		}

	} /* end else .. coherent section */

/*  output results and find min and max to echo
    remember that complex pix are stored in the file in FORTRAN
		order for compatability */
	if( lstart == 1 )
		for( ix=0; ix<NPARAM; ix++ ) param[ix] = sparam[ix];
	else
		for( ix=0; ix<NPARAM; ix++ ) param[ix] = 0.0F;
	param[pRMAX]  = rmax;
	param[pIMAX]  = aimax;
	param[pRMIN]  = rmin;
	param[pIMIN]  = aimin;
	param[pXCTILT]  = ctiltx;
	param[pYCTILT] = ctilty;
	param[pENERGY] = v0;
	param[pDX] = (float) ( ax/((float)nx) );
	param[pDY] = (float) ( by/((float)ny) );
	param[pWAVEL] = wavlen;
	param[pNSLICES] = 0.0F;  /* ??? */
	if ( lpartl == 1 ) {
	  param[pDEFOCUS] = df0;
	  param[pOAPERT] = aobj;
	  param[pCS] = Cs;
	  param[pCAPERT] = acmax;
	  param[pDDF] = sigmaf;
	}

	if ( lpartl == 1 )
		i = tcreateFloatPixFile( fileout, pix, (long) nx,
			 (long) ny, 1, param );
	else 
		i = tcreateFloatPixFile( fileout, waver, (long) (2*nx),
			 (long) ny, 2, param );

	if( i != 1 ) printf( "autoslice cannot write TIF file %s\n",
			fileout );
	printf( "pix range %g to %g real,\n"
	        "          %g to %g imag\n",  rmin,rmax,aimin,aimax );

	printf("Total CPU time = %f sec.\n", cputim()-timer );

} /* end main() */

/*--------------------- sortByZ() -----------------------*/
/*  sort atom info by coordinate z (Shell sort)

  x[], y[], z[]	= atom coordinates 
  occ[]		= occupancy of each atom
  Znum[]	= atomic number of each atom
  natom		= number of atoms
*/
void sortByZ( float x[], float y[], float z[], float occ[],
	int Znum[], int natom )
{
	int i, j, k, m;
	void iswap( int *i, int *j );
	void fswap( float *a, float *b );

	m = 1;
	j = 2;
	do{ m+=1; j*=2;} while( j <= natom );
	m = j / 2;

	do{		// ??? not most efficient stride but it works 
		k = natom - m;
		for( j=0; j<k; j++)
		for( i=j; i>=0; i-=m) {
			if( z[i+m] < z[i] ) {
				fswap( &x[i], &x[i+m] );
				fswap( &y[i], &y[i+m] );
				fswap( &z[i], &z[i+m] );
				fswap( &occ[i], &occ[i+m] );
				iswap( &Znum[i], &Znum[i+m] );
			}
		}
        	m = m/2;
	} while( m > 0 );
} 
/*------------------------ iswap() ------------------------*/
/*
  Swap 2 int's (for sorting)
  Can also be used for reals (via equivalence)
*/
void iswap( int *i, int *j )
{
        int it;
        it = *i;
        *i = *j;
        *j = it;
        return;
}  /* end iswap() */

/*------------------------ fswap() ------------------------*/
/*
  Swap 2 int's (for sorting)
  Can also be used for reals (via equivalence)
*/
void fswap( float *a, float *b )
{
        float t;
         t = *a;
        *a = *b;
        *b = t;
        return;
} /* end fswap() */


void QsortByZ( float x[], float y[], float z[], float occ[],
        int Znum[], int beg, int end )
{
	int i, j, k;
	void iswap( int *i, int *j );
	void fswap( float *a, float *b );
	if (end > beg +1)
	{
		float piv = z[beg];
		int l = beg + 1;
		int r = end;

		while (l <r)
		{
			if(z[l] <= piv){ 
				l++;
			}
			else{
				k = --r;
				fswap( &x[l], &x[k] );
                                fswap( &y[l], &y[k] );
                                fswap( &z[l], &z[k] );
                                fswap( &occ[l], &occ[k] );
                                iswap( &Znum[l], &Znum[k] );
			}
			
		
		}
		k = --l;	
		fswap( &x[k], &x[beg] );
		fswap( &y[k], &y[beg] );
		fswap( &z[k], &z[beg] );
		fswap( &occ[k], &occ[beg] );
		iswap( &Znum[k], &Znum[beg] );
		QsortByZ( x, y, z, occ, Znum, beg, l);
		QsortByZ( x, y, z, occ, Znum, r, end);
//
	}
		
}  /* end sortByZ()*/ 

/*--------------------- trlayer() -----------------------*/
/*
  Calculate complex specimen transmission function
  for one layer using real space projected atomic potentials

  x[],y[] = real array of atomic coordinates
  occ[]   = real array of occupancies
  Znum[]  = array of atomic numbers
  natom   = number of atoms
  ax, by  = size of transmission function in Angstroms
  kev     = beam energy in keV
  transr  = 2D array to get real part of specimen
		transmission function
  transi  = 2D array to get imag part of specimen
		transmission function
  nx, ny  = dimensions of transmission functions
  *phirms = average phase shift of projected atomic potential
  *nbeams = will get number of Fourier coefficients
  k2max   = square of max k = bandwidth limit

*/
void trlayer( const float x[], const float y[], const float occ[],
	const int Znum[], const int natom, 
	const float ax, const float by, const float kev,
	float **transr, float **transi, const long nx, const long ny,
	double *phirms, long *nbeams, const float k2max )
{
   int idx, idy, i, ixo, iyo, ix, iy, ixw, iyw, nx1, nx2, ny1, ny2;
   float k2;
   double r, rx2, vz, rmin, sum, scale;

   const double rmax=3.0; /* max atomic radius in Angstroms */

	scale = sigma( kev ) / 1000.0;  /* in 1/(volt-Angstroms) */

	/* min radius to avoid  singularity */
	rmin = ax/((double)nx);
	r = by/((double)ny);
	rmin =  0.25 * sqrt( 0.5*(rmin*rmin + r*r) );

	idx = (int) ( nx*rmax/ax ) + 1;
	idy = (int) ( ny*rmax/by ) + 1;

	for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++)
	    transr[ix][iy] = 0.0F;

	for( i=0; i<natom; i++) {
	   ixo = (int) ( nx*x[i]/ax );
	   iyo = (int) ( ny*y[i]/by );
	   nx1 = ixo - idx;
	   nx2 = ixo + idx;
	   ny1 = iyo - idy;
	   ny2 = iyo + idy;

	/* add proj. atomic potential at a local region near its center
	   taking advantage of small range of atomic potential */

	   for( ix=nx1; ix<=nx2; ix++) {
		rx2 = x[i] - ((double)ix)*ax/nx;
		rx2 = rx2 * rx2;
		ixw = ix;
		while( ixw < 0 ) ixw = ixw + nx;
		ixw = ixw % nx;
		for( iy=ny1; iy<=ny2; iy++) {
		   r = y[i] - ((double)iy)*by/ny;
		   r = sqrt( rx2 + r*r );
		   if( r <= rmax ) {
			iyw = iy;
			while( iyw < 0 ) iyw = iyw + ny;
			iyw = iyw % ny;
			if( r < rmin ) r = rmin;
			/* vz = occ[i] * scale * vzatom( Znum[i], r ); slow */
			vz = occ[i] * scale * vzatomLUT( Znum[i], r );
			transr[ixw][iyw] += (float) vz;
		   }
		} /* end for(iy... */
	   }  /* end for(ix... */

  	} /* end for(i=0... */

	/* convert phase to a complex transmission function */
	sum = 0;
	for( ix=0; ix<nx; ix++) 
	for( iy=0; iy<ny; iy++) {
		vz = transr[ix][iy];
		sum += vz;
		transr[ix][iy] = (float) cos( vz );
		transi[ix][iy] = (float) sin( vz );
	}

	*phirms = sum / ( ((double)nx)*((double)ny) );

	/* bandwidth limit the transmission function */
	*nbeams = 0;
	fft2d( transr, transi, nx, ny, +1);
	for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++) {
		k2 = ky2[iy] + kx2[ix];
		if (k2 < k2max) *nbeams += 1;
		else transr[ix][iy] = transi[ix][iy] = 0.0F;
	}
	fft2d( transr, transi, nx, ny, -1);
	
	return;

 }  /* end trlayer() */
 
/*--------------------- vzatomLUT() -----------------------------------*/
/*
	return the (real space) projected atomic potential for atomic
	number Z at radius r (in Angstroms)

	this mimics vzatom() in slicelib.c but uses a look-up-table
	with cubic spline interpolatoin to make it run about 2X-4X faster

	started 23-may-1997 E. Kirkland
	fix Z range to allow Hydrogen 1-jan-1998 ejk

	Z = atomic number 1 <= Z <= 98
	r = radius in Angstroms
*/

double vzatomLUT( int Z, double r )
{
   int i, iz;
   double dlnr, vz;

	if( splineInit == 0 ) {
	   splinx =  (double*) malloc( NRMAX * sizeof( double ) );
	   spliny = (double**) malloc( NZMAX * sizeof( double* ) );
	   splinb = (double**) malloc( NZMAX * sizeof( double* ) );
	   splinc = (double**) malloc( NZMAX * sizeof( double* ) );
	   splind = (double**) malloc( NZMAX * sizeof( double* ) );
	   if( (splinx==NULL) || (spliny==NULL) || (splinb==NULL) 
		   || (splinc==NULL) || (splind==NULL) ) {
		printf( "Cannot allocate spline pointer"
		  "in vzatomLUT()\n");
		exit( 0 );
	    }

	    /*  generate a set of logarithmic r values */
	    dlnr = log(RMAX/RMIN)/(NRMAX-1);
	    for( i=0; i<NRMAX; i++)
			splinx[i] = RMIN * exp( i * dlnr );
	    printf( "fit from r= %g to r= %g\n", splinx[0], splinx[NRMAX-1] );

	    nspline = int1D( NZMAX, "nspline" );
	    for( i=0; i<NZMAX; i++) nspline[i] = 0;
	    splineInit = 1;	/* remember that this has been done */
	}

	iz = Z - 1;  /* convert atomic number to array index */
	if( (Z < 1) || ( Z > NZMAX) ) {
		printf("Bad atomic number %d in vzatomLUT()\n", Z);
		exit( 0 );
	}

	/* if this atomic number has not been called before
		generate the spline coefficients */
	if( nspline[iz] == 0 ) {
	   spliny[iz] = double1D( NRMAX, "spliny" );
	   splinb[iz] = double1D( NRMAX, "splinb" );
	   splinc[iz] = double1D( NRMAX, "splinc" );
	   splind[iz] = double1D( NRMAX, "splind" );

	   for( i=0; i<NRMAX; i++)
		spliny[iz][i] = vzatom( Z, splinx[i] );
	   nspline[iz] = NRMAX;
	   splinh( splinx, spliny[iz], splinb[iz],
		  splinc[iz], splind[iz], NRMAX);
	}

	/* now that everything is set up find the
		scattering factor by interpolation in the table 
	*/

	vz = seval( splinx, spliny[iz], splinb[iz],
			   splinc[iz], splind[iz], nspline[iz], r );

	return( vz );

}  /* end vzatomLUT() */

