/*	*** fft2dc.c ***

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

	1D and 2D Fast Fourier Transforms ( FFT's )
	including a fast radix4,2 complex to complex FFT in C
	
	started 8-jan-1995 E. Kirkland
	added real to complex 1D and 2D FFT's 2-jun-1995 ejk
	reorganized rfft2d() 13-jun-1995 ejk
	added fft42t() and put in fft2d() 16-jul-1995 ejk
	added conditional compilation of fft42t() vs.
		fft42() because some CPUS are faster one
		way and some are faster the other 28-jul-1995 ejk
	add invert2D() and powerof2() 7-sep-1995 ejk
	removed possible free() of unallocated memory
			in fft42t() 21-sep-1997 ejk
	removed some non-essential routines 2-mar-1998 ejk

fft2d()   : 2D complex to complex FFT
rfft2d()  : 2D real to complex FFT
fft42()   : 1D complex to complex radix 4,2 FFT
fft42t()  : 1D complex to complex radix 4,2 FFT
		similar to fft42() but with look-up-tables
		faster for multiple calls of same length
invert2D(): move FFT center from corners to the center
powerof2(): return nearest power of 2 >= argument

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*  the symbol FFT1D must be define as one of:
	fft42	: 1D radix 4,2 complex-to-complex FFT routine
	fft42t	: same as fft42 but with Sin/Cos look-up-table
			- this is faster on some CPU's and slower on others
*/
#define FFT1D	fft42t
void FFT1D( float *, float *, long );

long powerof2( long );

/*-----------------  fft2d() ---------------- */
/*
	2D complex to complex FFT

	pixr[ix][iy] = real part of 2D pix to Fourier transform
	pixi[ix][iy] = imag part of 2D pix to Fourier transform
	nx,ny = (long int) size of array
		ix=0,1,2...(nx-1)
		iy=0,1,2...(ny-1)
	inverse = if > 0 then forward transform
	          if < 0 then forward transform

	On exit pixr and pixi will have the Fourier transform 
	and the orignal data will be lost.	
*/

void fft2d( float **pixr, float **pixi, long nx, long ny, int inverse )
{
	long ix, iy;
	float *tempr, *tempi, scale;

	/*  check that nx and ny are powers of two here */
	if( (powerof2(nx)!=nx) || (powerof2(ny)!=ny) ) {
		printf("nx=%ld, ny=%ld not apower of 2 in fft2D", nx, ny);
		exit(0);
	}

	tempr = (float*) malloc( (size_t) nx * sizeof( float ) );
	tempi = (float*) malloc( (size_t) nx * sizeof( float ) );
	if( (tempi == 0) || (tempr == 0) ) {
		printf("Cannot alloc temp storage in fft2D");
		exit(0);
	}
		
	/* complex conjugate before and after (and scale) to get inverse */
	if( inverse < 0 ) {
		for( ix=0; ix<nx; ix++) for(iy=0; iy<ny; iy++) 
			pixi[ix][iy] = -pixi[ix][iy];
	}

	for ( ix=0; ix<nx; ix++) {	/* transform the rows */
		FFT1D( pixr[ix], pixi[ix], ny );
	}
	
	for ( iy=0; iy<ny; iy++) {	/* transform the cols */
		for( ix=0; ix<nx; ix++) {
			tempr[ix] = pixr[ix][iy];
			tempi[ix] = pixi[ix][iy];
		}
		FFT1D( tempr, tempi, ny );
		for( ix=0; ix<nx; ix++) {
			pixr[ix][iy] = tempr[ix];
			pixi[ix][iy] = tempi[ix];
		}
	}

	/*  complex conjugate and scale if inverse */
	if( inverse < 0 ) {
		scale = (float) (1.0 / ( ((double) nx) * ((double)ny) ));
		for( ix=0; ix<nx; ix++) for(iy=0; iy<ny; iy++) {
			pixr[ix][iy] =  pixr[ix][iy]*scale;
			pixi[ix][iy] = -pixi[ix][iy]*scale;
		}
	}

	free( tempr );
	free( tempi );

}  /* end fft2d()  */

/*--------------------  rfft2d() -------------------- */
/*
	2D real to complex FFT

	The fwd FFT produces the positive half plane with kx=0
	and the the ky frequencies are stored in the normal
	wrap around manner iy=0 is DC, iy=1 is the lowest positive
	frequency, iy=(ny-1) is the lowest neg. frequency.
	The complex values are stored at adjacent ix values.
	for example (pixr[0][iy],pixr[1][iy]) = (real,imag)

	pixr[ix][iy] = real 2D pix to Fourier transform
	   NOTE: array size must be (nx+2) x ny
	nx,ny = (long int) size of original data (before FFT)
		ix=0,1,2...(nx-1) - see note below
		iy=0,1,2...(ny-1)
	inverse = if > 0 then forward transform
	          if < 0 then forward transform

	On exit pixr[][] will have the Fourier transform 
	and the orignal data will be lost.

	Although it is possible to pack the fwd FFT data into the
	same sized output array it is aranged in an awkward
	order so add 2 to the array size in the x direction for
	the high frequency component to make it easy to index
	the frequencies.
*/
void rfft2d( float **pixr, long nx, long ny, int inverse )
{
	long ix, iy, ix2;
	float *tempr, *tempi, scale;

	/*  check that nx and ny are powers of two here */
	if( (powerof2(nx)!=nx) || (powerof2(ny)!=ny) ) {
		printf("nx=%ld, ny=%ld not a power of 2 in fft2D", nx, ny);
		exit(0);
	}

	tempr = (float*) malloc( (size_t) nx * sizeof( float ) );
	tempi = (float*) malloc( (size_t) nx * sizeof( float ) );
	if( (tempr == 0) || (tempi ==0)  ) {
		printf("Cannot alloc temp storage in fft2D");
		exit(0);
	}

	if( inverse > 0 ) {	/* forward transform */

	   /* transform the rows as reals, two at a time  */
	   for ( iy=0; iy<ny; iy+=2) {
	   	for ( ix=0; ix<nx; ix++) {
			tempr[ix] = pixr[ix][iy];
			tempi[ix] = pixr[ix][iy+1];
		}
		FFT1D( tempr, tempi, nx );
		for( ix=1,ix2=(nx-1); ix<=(nx/2); ix++, ix2-- ) {
			pixr[2*ix  ][iy] = 0.5F*( tempr[ix] + tempr[ix2] );
			pixr[2*ix+1][iy] = 0.5F*( tempi[ix] - tempi[ix2] );

			pixr[2*ix  ][iy+1] = 0.5F*( tempi[ix] + tempi[ix2] );
			pixr[2*ix+1][iy+1] = 0.5F*(-tempr[ix] + tempr[ix2] );
		}
		pixr[0][iy] = tempr[0];
		pixr[0][iy+1] = tempi[0];
		pixr[1][iy] = pixr[1][iy+1] = 0.0F;  /* imag part is zero */
	   }  /* end for( iy... */

	   /* NOTE: although the cols at 0 and nx are real
		its easier to transform them as complex and this
		has a negligible effect on CPU time */

	   /* transform the cols as complex: even-ix=real, odd-ix=imag */
	   for ( ix=0; ix<=nx; ix+=2) {
		FFT1D( pixr[ix], pixr[ix+1], ny );
	   }
	   
	} else {	/* inverse transform */

	   /* transform the cols as complex: even-ix=real, odd-ix=imag
		and complex conjugate and scale to get inverse */
	   scale = (float) (1.0/ (((double)ny) *((double)nx)));
	   for ( ix=0; ix<=nx; ix+=2) {
		for( iy=0; iy<ny; iy++) {
			pixr[ix+1][iy] = -pixr[ix+1][iy];
		}
		FFT1D( pixr[ix], pixr[ix+1], ny );
		for( iy=0; iy<ny; iy++) {
			pixr[ix  ][iy] =  pixr[ix  ][iy] * scale;
			pixr[ix+1][iy] = -pixr[ix+1][iy] * scale;
		}
	   }

	   /* transform the rows as real, two at a time, 
		and complex conjugate and scale to get inverse */
	   for ( iy=0; iy<ny; iy+=2) {
		for( ix=1,ix2=(nx-1); ix<=(nx/2); ix++, ix2-- ) {
			tempr[ix] =  pixr[2*ix][iy] - pixr[2*ix+1][iy+1];
			tempi[ix] = -pixr[2*ix+1][iy] - pixr[2*ix][iy+1];

			tempr[ix2] = pixr[2*ix][iy] + pixr[2*ix+1][iy+1];
			tempi[ix2] = pixr[2*ix+1][iy] - pixr[2*ix][iy+1];
		}
		tempr[0] =   pixr[0][iy] - pixr[1][iy+1];
		tempi[0] = -(pixr[1][iy] + pixr[0][iy+1]);
		FFT1D( tempr, tempi, nx );
	   	for ( ix=0; ix<nx; ix++) {
			 pixr[ix][iy]   =  tempr[ix];
			 pixr[ix][iy+1] = -tempi[ix];
		}
	   }  /* end for( iy... */

	}  /* end inverse transform */

	free( tempr );
	free( tempi );

} /* end rfft2d() */

/*------------------------ fft42 --------------------------

	fft42( fr[], fi[], n )       radix-4,2 FFT in C

	fr[], fi[] = (float) real and imaginary array with input data
	n          = (long) size of array

  calculate the complex fast Fourier transform of (fr,fi)   
  input array fr,fi (real,imaginary) indexed from 0 to (n-1)
  on output fr,fi contains the transform 

  started 8-jan-1995 E. Kirkland
  
*/

void fft42 ( float *fr, float *fi, long n  )
{
#define	TWOPI	6.283185307

	long i, j, nv2, nm1, k, k0, k1, k2, k3, kinc, kinc2;
	float qr, qi, rr, ri, sr, si, tr, ti, ur, ui;
	double x1, w0r, w0i, w1r, w1i, w2r, w2i, w3r, w3i;
		
	kinc = n;

	while( kinc >= 4 ) {	/* start radix-4 section */
	
		kinc2 = kinc;
		kinc = kinc / 4;
		
		for( k0=0; k0<n; k0+=kinc2) {
			k1 = k0 + kinc;
			k2 = k1 + kinc;
			k3 = k2 + kinc;
			
			rr =  fr[k0] + fr[k2];    ri = fi[k0] + fi[k2];
			sr =  fr[k0] - fr[k2];    si = fi[k0] - fi[k2];
			tr =  fr[k1] + fr[k3];    ti = fi[k1] + fi[k3];
			ur = -fi[k1] + fi[k3];    ui = fr[k1] - fr[k3];
			
			fr[k0] = rr + tr;    fi[k0] = ri + ti;
			fr[k2] = sr + ur;    fi[k2] = si + ui;
			fr[k1] = rr - tr;    fi[k1] = ri - ti;
			fr[k3] = sr - ur;    fi[k3] = si - ui;
		}

		x1 = TWOPI/( (double) kinc2 );
		w0r = cos( x1 );   w0i = sin( x1 );
		w1r = 1.0;   w1i = 0.0;
		
		for( i=1; i<kinc; i++) {
			 x1 = w0r*w1r - w0i*w1i;    w1i = w0r*w1i + w0i*w1r;
			w1r = x1;
			w2r = w1r*w1r - w1i*w1i;    w2i = w1r*w1i + w1i*w1r;
			w3r = w2r*w1r - w2i*w1i;    w3i = w2r*w1i + w2i*w1r;

			for( k0=i; k0<n; k0+=kinc2) {
				k1 = k0 + kinc;
				k2 = k1 + kinc;
				k3 = k2 + kinc;
			
				rr =  fr[k0] + fr[k2];    ri = fi[k0] + fi[k2];
				sr =  fr[k0] - fr[k2];    si = fi[k0] - fi[k2];
				tr =  fr[k1] + fr[k3];    ti = fi[k1] + fi[k3];
				ur = -fi[k1] + fi[k3];    ui = fr[k1] - fr[k3];
			
				fr[k0] = rr + tr;    fi[k0] = ri + ti;

				qr = sr + ur;    qi = si + ui;
				fr[k2] = (float) (qr*w1r - qi*w1i);
				fi[k2] = (float) (qr*w1i + qi*w1r);

				qr = rr - tr;    qi = ri - ti;
				fr[k1] = (float) (qr*w2r - qi*w2i);
				fi[k1] = (float) (qr*w2i + qi*w2r);

				qr = sr - ur;    qi = si - ui;
				fr[k3] = (float) (qr*w3r - qi*w3i);
				fi[k3] = (float) (qr*w3i + qi*w3r);
			}
		}
		
	}  /*  end radix-4 section */

	while( kinc >= 2 ) {	/* start radix-2 section */
	
		kinc2 = kinc;
		kinc = kinc /2 ;
		
		x1 = TWOPI/( (double) kinc2 );
		w0r = cos( x1 );   w0i = sin( x1 );
		w1r = 1.0;   w1i = 0.0;
		
		for( k0=0; k0<n; k0+=kinc2 ){
			k1 = k0 + kinc;
			tr = fr[k0] - fr[k1];        ti = fi[k0] - fi[k1];
			fr[k0] = fr[k0] + fr[k1];    fi[k0] = fi[k0] + fi[k1];
			fr[k1] = tr;                 fi[k1] = ti;
		}
		
		for( i=1; i<kinc; i++) {
			 x1 = w0r*w1r - w0i*w1i;  w1i = w0r*w1i + w0i*w1r;
			w1r = x1;
			for( k0=i; k0<n; k0+=kinc2 ){
				k1 = k0 + kinc;
				tr = fr[k0] - fr[k1];        ti = fi[k0] - fi[k1];
				fr[k0] = fr[k0] + fr[k1];    fi[k0] = fi[k0] + fi[k1];
				fr[k1] = (float) (tr*w1r - ti*w1i);
				fi[k1] = (float) (tr*w1i + ti*w1r);
			}
		}
			
	}  /* end radix-2 section */

	nv2 = n / 2;
	nm1 = n - 1;
	j = 0;

	for (i=0; i< nm1; i++) {  /* reorder in bit reversed order */
		if( i < j ){
			tr = fr[j];     ti = fi[j];
			fr[j] = fr[i];  fi[j] = fi[i];
			fr[i] = tr;     fi[i] = ti; }
		k = nv2;
		while ( k <= j ) { j -=  k;  k = k>>1; }
		/* while ( k <= j ) { j = j - k;  k = k /2; }  is slower */
		j += k;
	}

#undef TWOPI
}  /* end fft42() */

/*------------------------ fft42t --------------------------

	fft42t( fr[], fi[], n )       radix-4,2 FFT in C

	fr[], fi[] = (float) real and imaginary array with input data
	n          = (long) size of array

  calculate the complex fast Fourier transform of (fr,fi)   
  input array fr,fi (real,imaginary) indexed from 0 to (n-1)
  on output fr,fi contains the transform 

	this is similar to fft42() but uses a look-up-table that
	speed it up a little if there are many calls of the same
	length (i.e. in multidimensional FFT's)

  started from fft42() 16-jul-1995 E. Kirkland
  
*/

void fft42t ( float *fr, float *fi, long n  )
{
#define	TWOPI	6.283185307
#define SWAP(x,y)	{tr=x; x=y; y=tr;}

	long i, j, nv2, nm1, k, k0, k1, k2, k3, kinc, kinc2;
	float qr, qi, rr, ri, sr, si, tr, ti, ur, ui;
	float  w1r, w1i, w2r, w2i, w3r, w3i;
	double x1, wr, wi, wrr, wii, w2;

	static float  *sint, *cost;	/* sin-cos look-up-table */
	static long ntable = -10;	/* length of table */
	
	if( ntable != n ) {	/* make a look-up-table */
		if( ntable < n ) {
		    if( ntable > 0 ) {
			free( cost );
			free( sint );
		    }
		    cost = (float*) malloc( ((size_t)n) * sizeof( float ) );
		    sint = (float*) malloc( ((size_t)n) * sizeof( float ) );
		}
		wr = 1.0;   wi = 0.0;
		x1 = TWOPI/( (double) n );
		wrr = cos( x1 );   wii = sin( x1 );
		for( i=0; i<n; i++) {
			cost[i] = (float) wr;
			sint[i] = (float) wi;
			w2 = wr*wrr - wi*wii;	/* recursive evaluation */
			wi = wr*wii + wi*wrr;
			wr = w2;
		}
		ntable = n;
	} /* end making look-up-table */

	kinc = n;

	while( kinc >= 4 ) {	/* start radix-4 section */
	
		kinc2 = kinc;
		kinc = kinc / 4;
		j = k = n/kinc2;
		
		for( k0=0; k0<n; k0+=kinc2) {
			k1 = k0 + kinc;
			k2 = k1 + kinc;
			k3 = k2 + kinc;
			
			rr =  fr[k0] + fr[k2];    ri = fi[k0] + fi[k2];
			sr =  fr[k0] - fr[k2];    si = fi[k0] - fi[k2];
			tr =  fr[k1] + fr[k3];    ti = fi[k1] + fi[k3];
			ur = -fi[k1] + fi[k3];    ui = fr[k1] - fr[k3];
			
			fr[k0] = rr + tr;    fi[k0] = ri + ti;
			fr[k2] = sr + ur;    fi[k2] = si + ui;
			fr[k1] = rr - tr;    fi[k1] = ri - ti;
			fr[k3] = sr - ur;    fi[k3] = si - ui;
		}

		for( i=1; i<kinc; i++) {
			w1r = cost[j];    w1i = sint[j];
			w2r = cost[j+j];    w2i = sint[j+j];
			w3r = cost[j+j+j];    w3i = sint[j+j+j];
			j += k;

			for( k0=i; k0<n; k0+=kinc2) {
				k1 = k0 + kinc;
				k2 = k1 + kinc;
				k3 = k2 + kinc;
			
				rr =  fr[k0] + fr[k2];    ri = fi[k0] + fi[k2];
				sr =  fr[k0] - fr[k2];    si = fi[k0] - fi[k2];
				tr =  fr[k1] + fr[k3];    ti = fi[k1] + fi[k3];
				ur = -fi[k1] + fi[k3];    ui = fr[k1] - fr[k3];
			
				fr[k0] = rr + tr;    fi[k0] = ri + ti;

				qr = sr + ur;    qi = si + ui;
				fr[k2] = qr*w1r - qi*w1i;
				fi[k2] = qr*w1i + qi*w1r;

				qr = rr - tr;    qi = ri - ti;
				fr[k1] = qr*w2r - qi*w2i;
				fi[k1] = qr*w2i + qi*w2r;

				qr = sr - ur;    qi = si - ui;
				fr[k3] = qr*w3r - qi*w3i;
				fi[k3] = qr*w3i + qi*w3r;
			}
		}
		
	}  /*  end radix-4 section */

	while( kinc >= 2 ) {	/* start radix-2 section */
	
		kinc2 = kinc;
		kinc = kinc /2 ;
		
		j = k = n/kinc2;
		
		for( k0=0; k0<n; k0+=kinc2 ){
			k1 = k0 + kinc;
			tr = fr[k0] - fr[k1];        ti = fi[k0] - fi[k1];
			fr[k0] = fr[k0] + fr[k1];    fi[k0] = fi[k0] + fi[k1];
			fr[k1] = tr;                 fi[k1] = ti;
		}
		
		for( i=1; i<kinc; i++) {
			w1r = cost[j];
			w1i = sint[j];
			j += k;
			for( k0=i; k0<n; k0+=kinc2 ){
				k1 = k0 + kinc;
				tr = fr[k0] - fr[k1];        ti = fi[k0] - fi[k1];
				fr[k0] = fr[k0] + fr[k1];    fi[k0] = fi[k0] + fi[k1];
				fr[k1] = tr*w1r - ti*w1i;    fi[k1] = tr*w1i + ti*w1r;
			}
		}
			
	}  /* end radix-2 section */

	/* remember a table of swap index values was actually
			 slower than this 21-sep-1997 ejk  */

	nv2 = n / 2;
	nm1 = n - 1;
	j = 0;

	for (i=0; i< nm1; i++) {  /* reorder in bit reversed order */
		if( i < j ){
			SWAP( fr[j], fr[i] );
			SWAP( fi[j], fi[i] );
		}
		k = nv2;
		while ( k <= j ) { j -=  k;  k = k>>1; }
		/* while ( k <= j ) { j = j - k;  k = k /2; }  is slower */
		j += k;
	}

#undef TWOPI
#undef SWAP

}  /* end fft42t() */

/*------------------------- powerof2() ------------------------*/
/*
	return the nearest power of 2 greater than or equal to 
	the argument
*/
long powerof2( long n )
{
	int ln;
	long n2;

	ln = 1;
	n2 = 2;
	while( (n2 < n) && (ln<31) ) { n2*=2; ln++; }

	return( n2 );
}

/*------------------------- invert2D() ----------------------*/
/*
	rearrange pix with corners moved to center (a la FFT's)

	 pix[ix][iy] = real array with image
	 nx,ny = range of image 0<ix<(nx-1) and 0<iy<(ny-1)

*/
void invert2D( float** pix, long nx, long ny )
{
#define SWAP(a,b)	{t=a; a=b; b=t;}

	long ix, iy, ixmid, iymid;
	float t;

	ixmid = nx/2;
	iymid = ny/2;

	for( ix=0; ix<nx; ix++) 
	for( iy=0; iy<iymid; iy++)
		SWAP( pix[ix][iy], pix[ix][iy+iymid] );

	for( ix=0; ix<ixmid; ix++) 
	for( iy=0; iy<ny; iy++)
		SWAP( pix[ix][iy], pix[ix+ixmid][iy] );

#undef SWAP
}
