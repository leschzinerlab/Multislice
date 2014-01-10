/*	*** fft2dc.h ***

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
	add IBM/ESSL scft() 8-sep-1995 ejk
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

void fft2d( float **pixr, float **pixi, long nx, long ny, int inverse );


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
		ix=0,1,2...(nx-1)
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
void rfft2d( float **pixr, long nx, long ny, int inverse );


/*------------------------ fft42 --------------------------

	fft42( fr[], fi[], n )       radix-4,2 FFT in C

	fr[], fi[] = (float) real and imaginary array with input data
	n          = (long) size of array

  calculate the complex fast Fourier transform of (fr,fi)   
  input array fr,fi (real,imaginary) indexed from 0 to (n-1)
  on output fr,fi contains the transform 

  started 8-jan-1995 E. Kirkland
  
*/

void fft42 ( float *fr, float *fi, long n );

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

void fft42t ( float *fr, float *fi, long n  );

/*------------------------- invert2D() ----------------------*/
/*
	rearrange pix with corners moved to center (a la FFT's)

	 pix[ix][iy] = real array with image
	 nx,ny = range of image 0<ix<(nx-1) and 0<iy<(ny-1)

*/
void invert2D( float** pix, long nx, long ny );

/*---------------------------- powerof2() ---------------------------*/
/*
	return the nearest power of 2 greater than or equal to 
	the argument
*/
long powerof2( long n );
