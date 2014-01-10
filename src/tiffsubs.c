/*  		*** tiffsubs.c ***

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

   Functions to work on TIFF format image files (in ANSI-C)
   The computer is assumed to have either big-endian or little-endian
   byte ordering.

   The extended TIFF routines tcreateFloatPix() and treadFloatPix()
   assume that the computer hardware uses IEEE floating point (these
   should still work if this is not true but the image files will not
   be transporable).

   These routines are thought to conform to TIFF version 6.0 
   (dated June 3, 1992) from:

	Aldus Developers Desk
	Aldus Corp.
	411 First Ave. South
	Seattle, WA  98104-2871  

   (phone (206) 628-6593 for a copy of the TIFF specifications):

   NOTE-1: This software is thought to be correct but absolutely no guarantee
	whether written or implied is given.  This software is supplied on a
	user-beware "as is" basis.

   NOTE-2: These routines read in either byte order (Intel=little endian or
	Motorola=big endian) but write only the byte order of the computer
	they are running on.  It is assumed that the computer supports
        byte addressing and ASCII.
 
   NOTE-3:  Various symbol definitions (at top of file) can be changed
       for different computers and modes (the symbol CompByteOrder
       will be automatically determined at run time):

CPU		OS	Compiler	WindowsGUI	printfOK
--------------------------------------------------------------------
IBM/PC		WIN 3.1	MSVC++1.0	defined		not defined
IBM/PC		WIN32	MSVC++4.2	not defined	defined
IBM/PC		Linux	gcc		not defined	defined
IBM RS6000	AIX 3.2	xlc		not defined	defined
Apple Mac	Sys 7	Think C++6.0	not defined	defined

The source code is formatted for a tab size of 8.

----------------------------------------------------------

The routines in this file are:

tclose:  close currently open tiff file and
             deallocate internal storage

tcreateFloatPixFile: create a floating point image file with an
	8 bit image at the beginning for viewing (can be read by
	treadFloatPix)

tcreatePixFile:  create an 8 or 16 bit greyscale image file
	in TIFF format

tFloatTest: check if a file has a valid TIFF header
          call before topenFloat() to avoid problems

tifferr:  common error handler - print messages
        (internal use only)

tifftest: check if a file has a valid TIFF header
          call before topen() to avoid problems

tinter:  interpolate in a float pix array 
        (used by tcreateFloatPixFile() )

tlist:	this routine list a few key parameters on the screen
	about the currently open file

topen: open a TIFF file for reading

topenFloat: open an extended TIFF floating point image
		 for reading

tread: read bytes from file with byte reversal if necessary
        (internal use only)

treadFloatPix: read a floating point image file created by
	tcreateFloatPixFile into memory

treadIFD:  read the IFD or directory of currently open file
        (internal use only)

treadPix: read a 'standard' TIFF image (8 or 16 bit integer)
	file into memory

tsetByteOrder:  determine the byte order of the computer that
       this is running on

tsize:  get image size of currently open TIFF file

----------------------------------------------------------
To write and output file requires only one function call
		tcreadtePixFile(...)  for integer images
  or		tcreateFloatPixFile(...) for floating point images
----------------------------------------------------------
To read in an integer image requires four function call:
	topen(...)
	tsize( ... &nx, &ny...)
	< allocate long** array of size nx by ny >
	treadPix(...)
	tclose()
----------------------------------------------------------
To read in a floating point image requires four function call:
	topenFloat(...)
	tsize( ... &nx, &ny...)
	< allocate float** array of size nx by ny >
	treadFloatPix(...)
	tclose()
----------------------------------------------------------

The floating point image routines tcreatFloatPixFile()
and treadFloatPix() use an extended TIFF format (i.e. they
conform to the official TIFF standard but are not a common
permutation).  A floating point image (real or 
complex=real+imaginary) and a parameter array are stored in
a single file.  The first image in the file is a standard
8 bit greyscale image which most TIFF readers should be able
to handle.  This 8 bit image may be expanded in one direction
using bilinear interpolation to get square pixels.  The second
image in the file is stored as 32 bit IEEE floating point (for
simulation data) and the third image is one line with a 64
element floating point parameter array.

The order of the data in the file is:

<TIFF header>
unsigned 8 bit image data with square pixels
<IFD-1>
32 bit floating point image data (pixels may be rectangular)
<IFD-2>
32 bit parameter data
<IFD-3>
----------------------------------------------------------

   started 10-Aug-1992 E. Kirkland
   all in working form  10-Oct-1992 ejk
   fixed bug in tcreate8bit (did not write out all
        of strip offset with large pix   5-Mar-1993 ejk
   fixed bug in writing images smaller than 8k
      (i.e. one strip)  20-aug-1993 ejk
   switch to common error message handler so that printf's etc
       can be turned off easily for GUI's like Windows 20-Nov-1993 ejk
   small changes for MS-Windows
   	(consistent types and #ifndef on printf's)  29-dec-1993 ejk
   changed tcreate8bit to tcreatePixFile and added ix0,iy0,nbits
   	and compress arguments  24-may-1994 ejk
   added 16 bit mode in tcreatePixFile  30-may-1994 ejk
   added 16 bit mode in treadPix  31-May-1994 ejk
   fixed byte ordering problem on II byte order 2-June-1994 ejk
   fixed bug in scaling of 16 bit images in tcreatePixFile
   	 16-july-1994 ejk
   fixed bug in reading 64x64x16bit images (i.e. nbytes was
       not calculated correctly)  4-Oct-1994 ejk
   added tsetByteOrder() to automatically determine byte order
       of computer that this is running on  18-mar-1996 ejk
   fixed aspect ratio in tcreatPix()  18-mar-1996 ejk
   added default to 8 BitsPerSample if a strange value
         is read   18-mar-1996 ejk
   added tcreateFloatPixFile()  22-mar-1996 ejk
   added treadFloatPix()  25-mar-1996 ejk
   added tinter() and topenFloat() 14-apr-1996 ejk
   remove 2nd free() in tcreateFloatPix() 19-feb-1997 ejk
   add tFloatTest() 20-apr-1997 ejk
   changed ulong and ushort to uslong and usshort beacause
      of name conflicts under AIX  25-apr-1997 ejk

*/                                       	

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <time.h>

/* define the byte order types  */
#define tBigEndian	0x4d4d		/* = 'MM' = Motorola */
#define tLittleEndian	0x4949		/*= 'II' = Intel */

/* CompByteOrder will get the byte order type for this computer
    - found dynamically by tsetByteOrder() -
  a value of 0 indicates that it has not been found yet */
short CompByteOrder=0;

/*#define WindowsGUI	   set for MS-Windows graphical user interface */
#define printfOK /*  set if we want to see error messages via printf's */ 

#ifdef WindowsGUI
#include "stdafx.h"   /* requires MSVC++ v1.0, 13-jun-1994 ejk */
#endif				

typedef unsigned short usshort;
typedef unsigned long uslong;
typedef struct { unsigned short tag, type;
		unsigned long length, value; } IFD;
	
FILE* tfp=NULL;		/* pointer to TIFF file */
short FileByteOrder;	/* byte order of the TIFF file */
short Version;		/* TIFF 'version', must be 42 */
uslong NextIFD;		/* pointer to next Image File Directory */
uslong CurrentIFD;	/* pointer to current Image File Directory */
short nIFD=0;		/* number of entries in IFD */
int nstrips=0;		/* number of strips in file */
int nstripcounts=0;	/* number of strip byte counts in file */
IFD* ifd1=NULL;		/* IFD data */

usshort* stempo=NULL; /* temporary short arrays for IFDread */
usshort* stempc=NULL;

const int sizes[] = { 0, 1, 1, 2, 4, 8 };  /* size in bytes of each */
					/* TIFF data type */

/*  TIFF parameters */

static usshort BitsPerSample[3], Compression, FillOrder,
	NewSubFileType,
	MaxSampleValue, MinSampleValue,
	Orientation, PhotometricInterpretation,
	PlanarConfiguration, Predictor, ResolutionUnit,
	SamplesPerPixel, StripByteCount,
	SubfileType, SampleFormat;

static uslong ImageLength, ImageWidth, NewSubfileType,
	RowsPerStrip, XResolution[2], YResolution[2],
	*StripByteCounts, *StripOffsets,
	XPosition[2], YPosition[2];

char DateTime[32]="No-Date-Given";

/*--------------------- tclose -------------------*/
/*
	close a TIFF file opened by topen
	and free malloc'ed data storage

	return 1 for success and </=0 for failure
*/

int tclose( )
{
   void tifferr( const char* );

	if( fclose( tfp ) != 0 ) {
		tifferr( "tclose cannot close the file.");
		return( -1 ); }

	if( ifd1 != NULL ) {
		free( ifd1 );  nIFD = 0;   ifd1 = NULL;
	}

	if( StripByteCounts != NULL) {  nstrips = 0;
		free( StripByteCounts); StripByteCounts = NULL; }

	if( StripOffsets != NULL) { nstripcounts = 0;
		free( StripOffsets );  StripOffsets = NULL; }

	if( stempo != NULL ) { free( stempo ); stempo = NULL; }
	if( stempc != NULL ) { free( stempc ); stempc = NULL; }

	return(1);

}  /* end tclose */

/*--------------------- tcreatePixFile ----------------*/
/*
	create an 8 or 16 bit greyscale image file
	in TIFF format

	Warning: this routine assumes that the CPU has
	enough memory to hold an entire image

	file[]   = (char) name of file to create	
	pix[][]  = long** array holding image indexed as 
			pix[iy][ix], iy=0..(ny-1), ix=0..(nx-1)
	nx,ny    = (int) size of image
	ixo, iy0 = (int) offset of pix in array (i.e. so more than one
			pix may be stored in one array
	nbits    = (int) 8 or 16  = number of bits per sample
	compress = (NOT YET IMPLEMENTED) 
		if 0 then do not compress else do an LZW compression
	psub     = long constant to be subtracted from each pixel
			 before scaling the image as 
			scale * (pix[iy][ix]-offset)
	scale    = double constant to scale the image before writing
	aspect   = aspect ratio = dy/dx
*/

int tcreatePixFile( const char *file, long** pix, long nx, long ny,
		int ixo, int iyo, int nbits, int compress,
		long psub, double scale, double aspect )
{
#define  ntags 13

   void tifferr( const char* );
   void tsetByteOrder();
   unsigned char uctemp, *StripBuf;
   unsigned short *SStripBuf, *stemp, st;
   short version;
   int i, strip, nbytes;
   long offset, ltemp, StripsPerImage, xy[2], nwrite;
   long *Offsets, *ByteCounts, ix, iy, nrows, count;
   IFD ifd2[ntags];
   FILE *tfp2;

	tsetByteOrder();
	if( nbits == 8 ) nbytes = 1;	/* allow 8 or 16 bits per pixel */
	else if( nbits == 16 ) nbytes = 2;
	else {
		tifferr( "Bad nbits in tcreatePixFile()." );
		return( -100 );
	}

	/* rows per strip, make strips about 8k */	
	nrows = 8192 / ( nx * nbytes );
	if( nrows <= 0 ) nrows = 1;
	if( nrows > ny ) nrows = ny;  /* if whole image < 1 strip */

	for( i=0; i<ntags; i++) ifd2[i].length = 1;  /* init tags */

	/* set up easy tags */

	ifd2[0].tag = 254;	/* NewSubfileType */
	ifd2[0].type = 4;
	ifd2[0].value = 0;

	ifd2[1].tag = 256;	/* ImageWidth */
	ifd2[1].type = 4;
	ifd2[1].value = (long) nx;

	ifd2[2].tag = 257;	/* ImageLength */
	ifd2[2].type = 4;
	ifd2[2].value = (long) ny;

	ifd2[3].tag = 258;	/* BitsPerSample */
	ifd2[3].type = 3;
	ifd2[3].value = 8*nbytes;

	ifd2[4].tag = 259;	/* Compression */
	ifd2[4].type = 3;
	ifd2[4].value = 1;

	ifd2[5].tag = 262;	/* PhotmetricInterpretation */
	ifd2[5].type = 3;
	ifd2[5].value = 1;

/*	ifd2[6].tag = 273;	 - StripOffsets - see below */

	ifd2[7].tag = 277;	/* SamplesPerPixel */
	ifd2[7].type = 3;
	ifd2[7].value = 1;

	ifd2[8].tag = 278;	/* RowsPerStrip */
	ifd2[8].type = 4;
	ifd2[8].value = nrows;

/*	ifd2[9].tag = 279;	- StripByteCounts - see below */

/*	ifd2[10].tag = 282;	- XResolution - see below */

/*	ifd2[11].tag = 283;	- YResolution - see below */

	ifd2[12].tag = 296;	/* ResolutionUnit */
	ifd2[12].type = 3;
	ifd2[12].value = 2;	/* inches */

	if( CompByteOrder == tBigEndian ) { /* word swap short values in IFD*/
		for( i=0; i<ntags; i++) 
			if( (ifd2[i].type == 3) && (ifd2[i].length <= 1 ) ) {
				stemp = (unsigned short*) &ifd2[i].value;
				st = stemp[0];  stemp[0] = stemp[1];
				stemp[1] = st;
			}
	}

	StripsPerImage = ( ny + nrows - 1 ) / nrows;

	Offsets = (long*) malloc( (size_t) StripsPerImage * sizeof(long) );
	if( Offsets == NULL ) {
		tifferr( "tcreatePixFile cannot allocate memory.");
		return( -1 ); }
	ByteCounts = (long*) malloc( (size_t) StripsPerImage * sizeof(long) );
	if( ByteCounts == NULL ) {
		free( Offsets );
		tifferr( "tcreatePixFile cannot allocate memory.");
		return( -2 ); }

	count = nrows * nx;			/* pixel count for StripBuf */
	SStripBuf = (unsigned short*) calloc( (size_t) count, nbytes );
	StripBuf  = (unsigned char*) SStripBuf;
	if( (StripBuf == NULL) || (SStripBuf == NULL) ) {
		free( Offsets );  free( ByteCounts );
		tifferr( "tcreatePixFile cannot allocate memory.");
		return( -3 ); }

	tfp2 = fopen( file, "w+b" );
	if( tfp2 == NULL) {
		free( Offsets);  free( ByteCounts ); free( SStripBuf );
		tifferr( "tcreatePixFile cannot open file.");
		return( -4 ); }

	offset = 8;  /* start writing after the header */
			/* i.e. fill in header later */
	if( fseek( tfp2, offset, SEEK_SET ) != 0 ) {
		free( Offsets);  free( ByteCounts ); free( StripBuf );
		tifferr( "Bad file seek in tcreatePixFile().");
		return( -5 ); }

	strip = 0;		/*  write image raster data */
	count = 0;
	for( iy=0; iy<ny; iy++){
	   for( ix=0; ix<nx; ix++) {
		ltemp = (long) ( scale * ( 
		   (double) ( pix[iy+iyo][ix+ixo] - psub ) ) + 0.5 );
		if( ltemp < 0 ) ltemp = 0;
		if( nbytes == 1 ) {
			if( ltemp > 255 ) ltemp = 255;
			StripBuf[count+ix] = (unsigned char) (ltemp & 0x00FF); 
		} else if( nbytes == 2 ) {
			if( ltemp > 65535L ) ltemp = 65535L;
			SStripBuf[count+ix] = (unsigned short) 
						(ltemp & 0x0FFFF);
		}
	   }
	   count += nx;

	   if( ( nrows*((iy+1)/nrows) == (iy+1) ) || ( (iy+1)==ny ) ) { 
	     if( (long) fwrite( StripBuf, sizeof(char)*nbytes,
			       (size_t) count,tfp2) != count) {
		   free( Offsets);  free( ByteCounts ); free( StripBuf );
		   tifferr("tcreatePixFile cannot write image data");
		   return( -6 ); }
		Offsets[ strip ] = offset;
		offset += count*nbytes;
		ByteCounts[strip] = count*nbytes;
		count = 0;  strip += 1;
	   } /* end if nrows */
	}  /* end for iy loop */

	nwrite = 0;
	ifd2[6].tag = 273;	/* StripOffsets */
	ifd2[6].type = 4;
	ifd2[6].length = strip;
	if( strip > 1 ) {
		ifd2[6].value = offset;
		nwrite += fwrite( Offsets, sizeof(long), strip, tfp2);
		offset += strip*sizeof(long);
	} else {
		ifd2[6].value = Offsets[0];
		nwrite += 1;
	}

	ifd2[9].tag = 279;	/* StripByteCounts */
	ifd2[9].type = 4;
	ifd2[9].length = strip;
	if( strip > 1 ) {
		ifd2[9].value = offset;
		nwrite += fwrite( ByteCounts, sizeof(long), strip, tfp2);
		offset += strip*sizeof(long);
	} else {
		ifd2[9].value = ByteCounts[0];
		nwrite += 1;
	}

	xy[0] = (long) ( (double) nx / 0.04 );  /* 4 inch wide image */
	xy[1] = 100;
	nwrite += fwrite( xy, sizeof(long), 2, tfp2);
	ifd2[10].tag = 282;	/* XResolution */
	ifd2[10].type = 5;
	ifd2[10].value = offset;
	offset += 2*sizeof(long);

	if( aspect > 0.0 ) xy[0] =   /* 4 inch wide with aspect ratio */
		(long) ( aspect* ((double)nx) / 0.04 );
	xy[1] = 100;
	nwrite += fwrite( xy, sizeof(long), 2, tfp2);
	ifd2[11].tag = 283;	/* YResolution */
	ifd2[11].type = 5;
	ifd2[11].value = offset;
	offset += 2*sizeof(long);

	/* remember to store IFD on a word boundary */
	if( (offset/2)*2 != offset ) {
		uctemp = 0;
		fwrite( &uctemp, 1, 1, tfp2 );
		offset += 1; }

	version = ntags;		/*  write the IFD */
	nwrite += fwrite( &version, sizeof(short), 1, tfp2);
	nwrite += fwrite( ifd2, sizeof(IFD), ntags, tfp2);
	ltemp = 0;
	nwrite += fwrite( &ltemp, sizeof(long), 1, tfp2);
	if( nwrite != ( 2*strip + ntags + 6 ) )
	   {  	tifferr("tcreatePixFile cannot write IFD.");
		free( Offsets);  free( ByteCounts ); free( StripBuf );
		return( -7 ); }

	fseek( tfp2, 0, SEEK_SET );
	nwrite = 0;
	nwrite += fwrite( &CompByteOrder, sizeof(short), 1, tfp2);
	version = 42;
	nwrite += fwrite( &version, sizeof(short), 1, tfp2);
	nwrite += fwrite( &offset, sizeof(long), 1, tfp2);
	if( nwrite !=  3 )
	   {    tifferr("tcreatePixFile cannot write header.");
		free( Offsets);  free( ByteCounts ); free( StripBuf );
		return( -8 ); }

	if( fclose( tfp2 ) != 0 ) {
		tifferr( "tcreatePixFile cannot close file.");
		free( Offsets);  free( ByteCounts ); free( StripBuf );
		return( -9 ); }

	free( Offsets);  free( ByteCounts ); free( StripBuf );
	return( 1 );

#undef ntags

}  /* end of tcreatePixFile */

/*--------------------- tcreateFloatPixFile ----------------*/
/*
	create a 32 bit floating point image file with an
	8 bit greyscale image in the front in TIFF format
	(the first image is a standard 8 bit greyscale image
	for a standard TIFF reader and the 2nd image is a
	32 bit floating point image for number crunching
	 in extended TIFF format - both represent the same pix )

	WARNING: This conforms to the TIFF standard but
	is an uncommon usage.  Also it assumes that the
	computer it is running on uses IEE floating point.

	file[]   = (char) name of file to create	
	pix[][]  = float** array holding image indexed as 
			pix[ix][iy], ix=0..(nx-1), iy=0..(ny-1)
		   if the image is complex then the real part
		   is contained in 0<ix<(nx/2-1) and the imaginary
		   part is within nx<ix<(nx-1) (see npix below)
		(NOTE: the order of ix,iy is opposite of the
			integer only version tcreatePixFile() )
	nx,ny    = (long) size of image
			 (nx includes BOTH real and imag parts)
	npix     = (int) real or complex (i.e. number of images)
			1 for real valued image
			2 for complex valued image

	param[64]= floating array of parameters

   NOTE: scaling etc. are controlled via elements of param[]

	param[0] is reserved for npix

	param[1] thru param[4] are used to scale the 8 bit image:
		param[1] = max real part
		param[2] = max imag part
		param[3] = min real part
		param[4] = min imag part

	param[14]=dx and param[15]=dy determine the aspect ratio

   started 18-mar-1996 ejk
   made 8-bit image pixels be square via interpolation/expansion
		14-apr-1996 ejk
   remove extra free()'s of non-existant memory 19-feb-1997 ejk
*/

int tcreateFloatPixFile( const char *file, float** pix,
	long nx, long ny, int npix, float param[] )
{
#define  ntags 14

   void tifferr( const char* );
   void tsetByteOrder();
   float tinter( float **pix, long nx, long ny, double x, double y );

   unsigned char uctemp, *StripBuf;
   unsigned short *stemp, st;
   short version;
   int i, strip;
   unsigned long offset, offsetIFD1, offsetIFD2, offsetIFD1to2,
	offsetIFD3, offsetIFD2to3;
   long ltemp, StripsPerImage, xy[2], nwrite, nx2, ny2;
   long *Offsets, *ByteCounts, ix, iy, nrows, count;
   IFD ifd2[ntags];
   FILE *tfp2;
   float scaler, rmin, scalei, imin, dx, dy, *FStripBuf;
   double x, y, z;

   time_t caltime;
   struct tm *mytime;

	tsetByteOrder();

   /* test for valid arguments */

	scaler = param[1] - param[3];
	if( scaler < 0.0F ) scaler = -scaler;
	if( scaler  < 1.0E-20 ) return( -1 );

	if( (npix<1) && (npix>2) ) return( -2 );

   /* first write the 8 bit image */
   /* expand to get square pixels so simple viewers will work
   	nx2, ny2 = new output size (in pixels)  */

	dx = param[14];
	if( dx <= 0.0F ) dx = 1.0F;
	dy = param[15];
	if( dy <= 0.0F ) dy = dx;
	if( dx > dy ) {
		nx2 = (long) ( nx * dx/dy + 0.5F);
		ny2 = ny;
	} else if( dx < dy ) {
		nx2 = nx;
		ny2 = (long) ( ny * dy/dx + 0.5F);;
	} else {
		nx2 = nx;
		ny2 = ny;
	}
	
	/* rows per strip, make strips about 8k */	
	nrows = 8192 / (nx2) ;
	if( nrows <= 0 ) nrows = 1;
	if( nrows > ny2 ) nrows = ny2;  /* if whole image < 1 strip */

	for( i=0; i<ntags; i++) ifd2[i].length = 1;  /* init tags */

	/* set up easy tags */

	ifd2[0].tag = 254;	/* NewSubfileType */
	ifd2[0].type = 4;
	ifd2[0].value = 0;

	ifd2[1].tag = 256;	/* ImageWidth */
	ifd2[1].type = 4;
	ifd2[1].value = (long) nx2;

	ifd2[2].tag = 257;	/* ImageLength */
	ifd2[2].type = 4;
	ifd2[2].value = (long) ny2;

	ifd2[3].tag = 258;	/* BitsPerSample */
	ifd2[3].type = 3;
	ifd2[3].value = 8;

	ifd2[4].tag = 259;	/* Compression */
	ifd2[4].type = 3;
	ifd2[4].value = 1;	/* no compression */

	ifd2[5].tag = 262;	/* PhotmetricInterpretation */
	ifd2[5].type = 3;
	ifd2[5].value = 1;	/* black is zero */

/*	ifd2[6].tag = 273;	 - StripOffsets - see below */

	ifd2[7].tag = 277;	/* SamplesPerPixel */
	ifd2[7].type = 3;
	ifd2[7].value = 1;

	ifd2[8].tag = 278;	/* RowsPerStrip */
	ifd2[8].type = 4;
	ifd2[8].value = nrows;

/*	ifd2[9].tag = 279;	- StripByteCounts - see below */

/*	ifd2[10].tag = 282;	- XResolution - see below */

/*	ifd2[11].tag = 283;	- YResolution - see below */

	ifd2[12].tag = 296;	/* ResolutionUnit */
	ifd2[12].type = 3;
	ifd2[12].value = 2;	/* inches */

	ifd2[13].tag = 339;	/* SampleFormat */
	ifd2[13].type = 3;
	ifd2[13].value = 1;	/* unsigned integer data */

	{ /* word swap short values in IFD*/
	if( CompByteOrder == tBigEndian )
		for( i=0; i<ntags; i++) 
		   if( (ifd2[i].type == 3) && (ifd2[i].length <= 1 ) ) {
			stemp = (unsigned short*) &ifd2[i].value;
			st = stemp[0];  stemp[0] = stemp[1];
			stemp[1] = st;
		}
	}

	StripsPerImage = ( ny2 + nrows - 1 ) / nrows;

	Offsets = (long*) malloc( (size_t) StripsPerImage * sizeof(long) );
	if( Offsets == NULL ) {
		tifferr( "tcreateFloatPixFile cannot allocate memory.");
		return( -3 ); }
	ByteCounts = (long*) malloc( (size_t) StripsPerImage * sizeof(long) );
	if( ByteCounts == NULL ) {
		free( Offsets );
		tifferr( "tcreateFloatPixFile cannot allocate memory.");
		return( -4 ); }

	count = nrows * nx2;	/* pixel count for StripBuf */
	StripBuf = (unsigned char*) malloc( (size_t) count );
	if( StripBuf == NULL ) {
		free( Offsets );  free( ByteCounts );
		tifferr( "tcreateFloatPixFile cannot allocate memory.");
		return( -5 ); }

	tfp2 = fopen( file, "w+b" );
	if( tfp2 == NULL) {
		free( Offsets);  free( ByteCounts ); free( StripBuf );
		tifferr( "tcreateFloatPixFile cannot open file.");
		return( -6 ); }

	offset = 8;	/* start writing after the header */
			/* i.e. fill in header later */
	if( fseek( tfp2, offset, SEEK_SET ) != 0 ) {
		free( Offsets);  free( ByteCounts ); free( StripBuf );
		tifferr( "Bad file seek in tcreateFloatPixFile().");
		return( -6 ); }

	/*  write image raster data --
		note this will interpolate inbetween real and image
			images in a strange way */
	strip = 0;
	count = 0;
	rmin = param[3];
	imin = param[4];
	scaler = 255.0F/(param[1] - rmin);
	if( npix > 1 ) scalei = 255.0F/(param[2] - imin);
	for( iy=0; iy<ny2; iy++){
	   y = (ny-1) * ((double)iy)/((double)(ny2-1)) ;
	   for( ix=0; ix<nx2; ix++) {
	   	x = (nx-1) * ((double)ix)/((double)(nx2-1)) ;
	   	z = tinter( pix, nx, ny, x, y);
		if( ix < nx2/npix ) ltemp = (long) ( scaler*(z-rmin) + 0.5 );
			       else ltemp = (long) ( scalei*(z-imin) + 0.5 );
		if( ltemp < 0 ) ltemp = 0;
		if( ltemp > 255 ) ltemp = 255;
		StripBuf[count+ix] = (unsigned char) (ltemp & 0x00FF); 
	   }
	   count += nx2;

	   if( ( nrows*((iy+1)/nrows) == (iy+1) ) || ( (iy+1)==ny2 ) ) { 
	     if( (long) fwrite( StripBuf, sizeof(char),
			       (size_t) count,tfp2) != count) {
		   free( Offsets);  free( ByteCounts ); free( StripBuf );
		   tifferr("tcreateFloatPixFile cannot write image data");
		   return( -8 ); }
		Offsets[ strip ] = offset;
		offset += count;
		ByteCounts[strip] = count;
		count = 0;  strip += 1;
	   } /* end if nrows */
	}  /* end for iy loop */

	nwrite = 0;
	ifd2[6].tag = 273;	/* StripOffsets */
	ifd2[6].type = 4;
	ifd2[6].length = strip;
	if( strip > 1 ) {
		ifd2[6].value = offset;
		nwrite += fwrite( Offsets, sizeof(long), strip, tfp2);
		offset += strip*sizeof(long);
	} else {
		ifd2[6].value = Offsets[0];
		nwrite += 1;
	}

	ifd2[9].tag = 279;	/* StripByteCounts */
	ifd2[9].type = 4;
	ifd2[9].length = strip;
	if( strip > 1 ) {
		ifd2[9].value = offset;
		nwrite += fwrite( ByteCounts, sizeof(long), strip, tfp2);
		offset += strip*sizeof(long);
	} else {
		ifd2[9].value = ByteCounts[0];
		nwrite += 1;
	}

	xy[0] = (long) ( ((double) nx2) / 0.04 );  /* 4 inch wide image */
	xy[1] = 100;
	nwrite += fwrite( xy, sizeof(long), 2, tfp2);
	ifd2[10].tag = 282;	/* XResolution */
	ifd2[10].type = 5;
	ifd2[10].value = offset;
	offset += 2*sizeof(long);

	/* remember that the aspect ratio is 1 for first 8 bit pix */
	nwrite += fwrite( xy, sizeof(long), 2, tfp2);
	ifd2[11].tag = 283;	/* YResolution */
	ifd2[11].type = 5;
	ifd2[11].value = offset;
	offset += 2*sizeof(long);

	/* remember to store IFD on a word boundary */
	if( (offset/2)*2 != offset ) {
		uctemp = 0;
		fwrite( &uctemp, 1, 1, tfp2 );
		offset += 1; }

	offsetIFD1 = offset;	/* pointer to first IFD */

	version = ntags;	/*  write the first IFD */
	nwrite += fwrite( &version, sizeof(short), 1, tfp2);
	nwrite += fwrite( ifd2, sizeof(IFD), ntags, tfp2);

	/* save location of pointer to 2nd IFD, and write 0 for now */
	offset += sizeof(short) + ntags*sizeof(IFD);
	offsetIFD1to2 = offset;
	ltemp = 0;
	nwrite += fwrite( &ltemp, sizeof(long), 1, tfp2);
	if( nwrite != ( 2*strip + ntags + 6 ) )
	   {  	tifferr("tcreateFloatPixFile cannot write IFD.");
		free( Offsets);  free( ByteCounts ); free( StripBuf );
		return( -9 ); }
	offset += sizeof(long);

	free( StripBuf );
 
   /* next write the 32 bit floating point image */

	/* rows per strip, make strips about 8k */	
	nrows = 8192 / ( nx * sizeof(float) );
	if( nrows <= 0 ) nrows = 1;
	if( nrows > ny ) nrows = ny;  /* if whole image < 1 strip */

	/* change appropriate tags - most are the same */

	ifd2[1].tag = 256;	/* ImageWidth */
	ifd2[1].type = 4;
	ifd2[1].value = (long) nx;

	ifd2[2].tag = 257;	/* ImageLength */
	ifd2[2].type = 4;
	ifd2[2].value = (long) ny;

	ifd2[3].tag = 258;	/* BitsPerSample */
	ifd2[3].type = 3;
	ifd2[3].value = 8 * sizeof(float);

	ifd2[4].tag = 259;	/* Compression */
	ifd2[4].type = 3;
	ifd2[4].value = 1;	/* no compression */

	ifd2[5].tag = 262;	/* PhotmetricInterpretation */
	ifd2[5].type = 3;
	ifd2[5].value = 1;	/* black is zero */

/*	ifd2[6].tag = 273;	 - StripOffsets - see below */

	ifd2[7].tag = 277;	/* SamplesPerPixel */
	ifd2[7].type = 3;
	ifd2[7].value = 1;

	ifd2[8].tag = 278;	/* RowsPerStrip */
	ifd2[8].type = 4;
	ifd2[8].value = nrows;

/*	ifd2[9].tag = 279;	- StripByteCounts - see below */

/*	ifd2[10].tag = 282;	- XResolution - see below */

/*	ifd2[11].tag = 283;	- YResolution - see below */

	ifd2[12].tag = 296;	/* ResolutionUnit */
	ifd2[12].type = 3;
	ifd2[12].value = 2;	/* inches */

	ifd2[13].tag = 339;	/* SampleFormat */
	ifd2[13].type = 3;
	ifd2[13].value = 3;	/* IEEE floating point */

	if( CompByteOrder == tBigEndian ) { /* word swap short values in IFD*/
	   for( i=3; i<ntags; i++) {
	      if( (ifd2[i].type == 3) && (ifd2[i].length <= 1 ) ) {
		stemp = (unsigned short*) &ifd2[i].value;
		st = stemp[0];  stemp[0] = stemp[1];
		stemp[1] = st;
	      }
	   }
	}

	StripsPerImage = ( ny + nrows - 1 ) / nrows;

	Offsets = (long*) realloc( Offsets,
			(size_t) StripsPerImage * sizeof(long) );
	if( Offsets == NULL ) {
		tifferr( "tcreateFloatPixFile cannot allocate memory.");
		return( -10 ); }
	ByteCounts = (long*) realloc( ByteCounts,
			(size_t) StripsPerImage * sizeof(long) );
	if( ByteCounts == NULL ) {
		free( Offsets );
		tifferr( "tcreateFloatPixFile cannot allocate memory.");
		return( -11 ); }

	count = nrows * nx ;	/* pixel count for StripBuf */
	FStripBuf = (float*) calloc( (size_t) count, sizeof(float) );
	if( FStripBuf == NULL ) {
		free( Offsets );  free( ByteCounts );
		tifferr( "tcreateFloatPixFile cannot allocate memory.");
		return( -12 ); }

	/* continue writing where we left off - after the first IFD */
	/* it should be positioned OK from 8 bit write above */

	strip = 0;		/*  write image raster data */
	count = 0;
	for( iy=0; iy<ny; iy++){
	   for( ix=0; ix<nx; ix++) {
		FStripBuf[count+ix] = pix[ix][iy];
	   }
	   count += nx;

	   if( ( nrows*((iy+1)/nrows) == (iy+1) ) || ( (iy+1)==ny ) ) { 
	     if( (long) fwrite( FStripBuf, sizeof(float),
			       (size_t) count,tfp2) != count) {
		   free( Offsets);  free( ByteCounts ); free( FStripBuf );
		   tifferr("tcreateFloatPixFile cannot write image data");
		   return( -13 ); }
		Offsets[ strip ] = offset;
		offset += count*sizeof(float);
		ByteCounts[strip] = count*sizeof(float);
		count = 0;  strip += 1;
	   } /* end if nrows */
	}  /* end for iy loop */

	free( FStripBuf );
	nwrite = 0;
	ifd2[6].tag = 273;	/* StripOffsets */
	ifd2[6].type = 4;
	ifd2[6].length = strip;
	if( strip > 1 ) {
		ifd2[6].value = offset;
		nwrite += fwrite( Offsets, sizeof(long), strip, tfp2);
		offset += strip*sizeof(long);
	} else {
		ifd2[6].value = Offsets[0];
		nwrite += 1;
	}

	ifd2[9].tag = 279;	/* StripByteCounts */
	ifd2[9].type = 4;
	ifd2[9].length = strip;
	if( strip > 1 ) {
		ifd2[9].value = offset;
		nwrite += fwrite( ByteCounts, sizeof(long), strip, tfp2);
		offset += strip*sizeof(long);
	} else {
		ifd2[9].value = ByteCounts[0];
		nwrite += 1;
	}

	xy[0] = (long) ( (double) nx / 0.04 );  /* 4 inch wide image */
	xy[1] = 100;
	nwrite += fwrite( xy, sizeof(long), 2, tfp2);
	ifd2[10].tag = 282;	/* XResolution */
	ifd2[10].type = 5;
	ifd2[10].value = offset;
	offset += 2*sizeof(long);

	dx = param[14];  /* determine the aspect ratio */
	dy = param[15];
	if( (dx > 0.0F) && (dy > 0.0F) )
		xy[0] =	(long) ( (dy/dx)* ((double)nx) / 0.04 );
	xy[1] = 100;
	nwrite += fwrite( xy, sizeof(long), 2, tfp2);
	ifd2[11].tag = 283;	/* YResolution */
	ifd2[11].type = 5;
	ifd2[11].value = offset;
	offset += 2*sizeof(long);

	/* remember to store IFD on a word boundary */
	if( (offset/2)*2 != offset ) {
		uctemp = 0;
		fwrite( &uctemp, 1, 1, tfp2 );
		offset += 1; }

	offsetIFD2 = offset;
	version = ntags;	/*  write the second IFD */
	nwrite += fwrite( &version, sizeof(short), 1, tfp2);
	nwrite += fwrite( ifd2, sizeof(IFD), ntags, tfp2);
	offset += sizeof(short) + ntags*sizeof(IFD);
	offsetIFD2to3 = offset;
	ltemp = 0;
	nwrite += fwrite( &ltemp, sizeof(long), 1, tfp2);
	if( nwrite != ( 2*strip + ntags + 6 ) )
	   {  	tifferr("tcreateFloatPixFile cannot write IFD2.");
		free( Offsets);  free( ByteCounts );
		return( -14 ); }
	offset += sizeof(long);

  /* next write the 32 bit floating point parameters as a 
	one line image */

	/* one strip and one row per strip */	
	nrows = 1;

	/* set up easy tags */

	ifd2[1].tag = 256;	/* ImageWidth */
	ifd2[1].type = 4;
	ifd2[1].value = 64L;

	ifd2[2].tag = 257;	/* ImageLength */
	ifd2[2].type = 4;
	ifd2[2].value = 1L;

	ifd2[3].tag = 258;	/* BitsPerSample */
	ifd2[3].type = 3;
	ifd2[3].value = 8 * sizeof(float);

	ifd2[8].tag = 278;	/* RowsPerStrip */
	ifd2[8].type = 4;
	ifd2[8].value = nrows;

	if( CompByteOrder == tBigEndian ) { /* word swap short values in IFD*/
	   i = 3;
	   if( (ifd2[i].type == 3) && (ifd2[i].length <= 1 ) ) {
		stemp = (unsigned short*) &ifd2[i].value;
		st = stemp[0];  stemp[0] = stemp[1];
		stemp[1] = st;
	   }
	}

	StripsPerImage = 1;

	/* continue writing where we left off - after the second IFD */
	/* it should be positioned OK from 8 bit write above */

	count = 64;
	strip = 0;
	param[0] = (float) npix;	/* make sure this is right */
	if( (long) fwrite( param, sizeof(float),
		       (size_t) count,tfp2) != count) {
	   free( Offsets);  free( ByteCounts );
	   tifferr("tcreateFloatPixFile cannot write image data");
	   return( -15 ); }
	Offsets[ 0 ] = offset;
	offset += count*sizeof(float);
	ByteCounts[ 0 ] = count*sizeof(float);

	ifd2[6].tag = 273;	/* StripOffsets */
	ifd2[6].type = 4;
	ifd2[6].length = 1;
	ifd2[6].value = Offsets[0];

	ifd2[9].tag = 279;	/* StripByteCounts */
	ifd2[9].type = 4;
	ifd2[9].length = 1;
	ifd2[9].value = ByteCounts[0];

	nwrite = 0;
	xy[0] = (long) ( (double) nx / 0.04 );  /* 4 inch wide image */
	xy[1] = 100;
	nwrite += fwrite( xy, sizeof(long), 2, tfp2);
	ifd2[10].tag = 282;	/* XResolution */
	ifd2[10].type = 5;
	ifd2[10].value = offset;
	offset += 2*sizeof(long);

	dx = param[14];  /* determine the aspect ratio */
	dy = param[15];
	if( (dx > 0.0F) && (dy > 0.0F) )
		xy[0] =	(long) ( (dy/dx)* ((double)nx) / 0.04 );
	xy[1] = 100;
	nwrite += fwrite( xy, sizeof(long), 2, tfp2);
	ifd2[11].tag = 283;	/* YResolution */
	ifd2[11].type = 5;
	ifd2[11].value = offset;
	offset += 2*sizeof(long);

	/* use the ResolutionUnit tag for the date/time */
	caltime = time( NULL );
	mytime = localtime( &caltime );
	strftime( DateTime, 20, "%Y:%m:%d %H:%M:%S", mytime );
	ifd2[12].tag = 306;	/* DateTime */
	ifd2[12].type = 2;	/* ASCII */
	ifd2[12].value = offset;
	nwrite += fwrite( DateTime, sizeof(char), 20, tfp2);
	offset += 20*sizeof(char);

	/* remember to store IFD on a word boundary */
	if( (offset/2)*2 != offset ) {
		uctemp = 0;
		fwrite( &uctemp, 1, 1, tfp2 );
		offset += 1; }

	offsetIFD3 = offset;
	version = ntags;	/*  write the second IFD */
	nwrite += fwrite( &version, sizeof(short), 1, tfp2);
	nwrite += fwrite( ifd2, sizeof(IFD), ntags, tfp2);
	ltemp = 0;
	nwrite += fwrite( &ltemp, sizeof(long), 1, tfp2);
	if( nwrite != ( ntags + 26 ) )
	   {  	tifferr("tcreateFloatPixFile cannot write IFD3.");
		free( Offsets);  free( ByteCounts );
		return( -16 ); }

  /*  all done write the header */

	fseek( tfp2, 0, SEEK_SET );
	nwrite = 0;
	nwrite += fwrite( &CompByteOrder, sizeof(short), 1, tfp2);
	version = 42;
	nwrite += fwrite( &version, sizeof(short), 1, tfp2);
	nwrite += fwrite( &offsetIFD1, sizeof(long), 1, tfp2);
	if( nwrite !=  3 )
	   {    tifferr("tcreateFloatPixFile cannot write header.");
		free( Offsets);  free( ByteCounts );
		return( -17 ); }

   /*  fix up the IFD pointers */

	fseek( tfp2, offsetIFD1to2, SEEK_SET );
	nwrite = fwrite( &offsetIFD2, sizeof(long), 1, tfp2);
	if( nwrite !=  1 )
	   {    tifferr("tcreateFloatPixFile cannot write header.");
		free( Offsets);  free( ByteCounts );
		return( -18 ); }

	fseek( tfp2, offsetIFD2to3, SEEK_SET );
	nwrite = fwrite( &offsetIFD3, sizeof(long), 1, tfp2);
	if( nwrite !=  1 )
	   {    tifferr("tcreateFloatPixFile cannot write header.");
		free( Offsets);  free( ByteCounts );
		return( -19 ); }

   /* close and exit */
	if( fclose( tfp2 ) != 0 ) {
		tifferr( "tcreateFloatPixFile cannot close file.");
		free( Offsets);  free( ByteCounts );
		return( -20 ); }

	free( Offsets);  free( ByteCounts );
	return( 1 );

#undef ntags

}  /* end of tcreateFloatPixFile() */


/*--------------------- tFloatTest -------------------*/
/*
   Test that a file is an extended TIFF file and
   has a floating point image (not an exhaustive test)

   filename  = pointer to string with filename

   return 1 for success and </=0 for failure
*/

int tFloatTest( const char *filename )
{
   void tifferr( const char* );
   void tsetByteOrder();

   int iread;
   int tifftest( const char *filename );
   int tread( void *bufptr, int size, int n, long offset);
   int treadIFD( unsigned long IFDoffset );

	tsetByteOrder();
	if( tifftest( filename ) != 1 ) return( -1 );

	tfp = fopen( filename, "rb" );
	if( tfp == NULL ) {
		tifferr( "tFloatTest() bad filename." );
		return( -2 ); }

	iread  = fread( &FileByteOrder, sizeof(short), 1, tfp);

	if( (iread != 1) || 
	   !( (FileByteOrder==tLittleEndian) || 
                      (FileByteOrder==tBigEndian) ) ) {
		fclose( tfp );;
		return( -3 ); }

	tread( &Version, 2, 1, 2L );
	if( Version != 42 ) {    	/* the TIFF magic number */
		fclose( tfp );;
		return( -3 ); }

	/* get location of IFD  and read it */

	tread( &CurrentIFD, 4, 1, 4L );
	iread = treadIFD( CurrentIFD );

	if( iread != 1 ) {
		tclose();
		return( -4 ); }

	/* find the first floating point image   */
	while( ( SampleFormat != 3 ) &&
	       ( BitsPerSample[0] != 8*sizeof(float) ) &&
	       ( NextIFD > 0 ) ) {
		if( treadIFD( NextIFD ) < 1 ) {
			tclose( );;
			return( -5 );
		}
	}

	tclose();
	if( ( SampleFormat != 3 ) ||
	       ( BitsPerSample[0] != 8*sizeof(float) ) ) return( -6 );

	return( 1 );

}  /* end tFloatTest() */


/*--------------------- tifferr ------------------*/
/*
    common error handler - print messages
    only for internal use 
    DO NOT CALL FROM MAIN PROGRAM
*/

void tifferr( const char* error_text )
{
#ifndef WindowsGUI

	if( errno != 0 )
		fprintf( stderr, "System Error: %s\n", strerror( errno ) );
	fprintf( stderr, "TIFF Error: %s\n", error_text );
	
#else
	char stemp[500];
	
	if( errno != 0 ) {
	   sprintf( stemp, "System Error: %s\n", strerror( errno ) );
 		AfxMessageBox( stemp, MB_OK, 0);
 	}
	sprintf( stemp, "TIFF Error: %s\n", error_text );
 	AfxMessageBox( stemp, MB_OK, 0);
 
#endif

	return;

}  /* end tifferr() */


/*--------------------- tifftest -------------------*/
/*
   check if a file has a valid TIFF header

   return a 1 if it is a valid TIFF header, </=0 otherwise
   filename  = pointer to string with filename
*/
int tifftest( const char *filename )
{
   void tifferr( const char* );
   void tsetByteOrder();
   FILE *tfp3;
   short byte01, byte23;
   int iread;

	tsetByteOrder();

	tfp3 = fopen( filename, "rb" );
	if( tfp3 == NULL ) {
	   tifferr( "Bad filename in tifftest");
	   return( -1 ); }

	iread   = fread( &byte01, sizeof(short), 1, tfp3);
	iread  += fread( &byte23, sizeof(short), 1, tfp3);

	if( fclose( tfp3 ) != 0) {
	   tifferr( "tifftest cannot close file.");
	   return( -2 ); }

	if( iread != 2 ) {
	   tifferr( "tifftest cannot read header.");
	   return(-3); }

   /*   tLittleEndian = "II" = least significant byte first 
	tBigEndian    = "MM" = most significant byte first
	42 = 0x2A = TIFF magic number  */

	if( (byte01!=tLittleEndian) && 
	    (byte01!=tBigEndian) ) return(-4);
	
	if( (byte01==CompByteOrder) && (byte23!=0x002A) ) return(-5);
	if( (byte01!=CompByteOrder) && (byte23!=0x2A00) ) return(-5);

	return( +1 );

}  /* end tifftest */

/*--------------------- tinter() ----------------------------*/
/*
  Bilinear interpolation from pix array

  pix[][]  = real input array with the image
  nx,ny    = integer dimension of the image
  x,y      = real coord of image point to interpolate at
  		0->(nx-1) and 0->(ny-1)
*/
float tinter( float **pix, long nx, long ny, double x, double y )
{
	long ix, iy, ix2, iy2;
	float  a, b, c, d, x1, y1, ans;

 	ix = (int) x;
	iy = (int) y;
	if( ix > (nx-2) ) ix = (nx-2);
	if( ix <    0 )   ix = 0;
	if( iy > (ny-2) ) iy = (ny-2);
	if( iy <    0 )   iy = 0;
	ix2 = ix + 1;
	iy2 = iy + 1;
	x1 = (float) ix;
	y1 = (float) iy;

	d = pix[ix][iy]  - pix[ix][iy2] - pix[ix2][iy]
				 + pix[ix2][iy2];
	c = pix[ix][iy2] - pix[ix][iy]  - d*x1;
	b = pix[ix2][iy] - pix[ix][iy]  - d*y1;
	a = pix[ix][iy]  - b*x1 - c*y1 - d*x1*y1;

	ans = (float) ( a + b*x + c*y + d*x*y );
	
	return( ans );

}  /* end tinter() */

/*--------------------- tlist -------------------*/
/*
	this routine list a few key parameters on the screen
	about the currently open file
*/

void tlist( )
{
#ifdef	printfOK	/*  remember printf's don't work inside a GUI !! */
   int i;
   char *c;
   
   	c = (char*) &FileByteOrder;
	printf( "File ByteOrder = %x (hex) = %c%c (ascii)\n", FileByteOrder,
	       c[0], c[1] );
   	c = (char*) &CompByteOrder;
	printf( "Computer ByteOrder =  %x (hex) = %c%c (ascii)\n",
	       CompByteOrder, c[0], c[1] );
	printf( "TIFF Version # = %d (decimal)\n", Version );
	printf( "Number of entries in IFD = %d (decimal)\n", nIFD );
	printf( "Current IFD at = %ld (decimal)\n", CurrentIFD );
	printf( "Next IFD at    = %ld (decimal)\n", NextIFD );
	printf( "nstrips = %d, nstripcounts = %d\n", nstrips,
		nstripcounts);

	if( nIFD > 0 ) {
	   printf("\nTIFF Image File Directory parameters,\n");
	   printf("    (some may be defaults):\n");
	   printf("BitsPerSample(258) = ");
		for(i=0; i<((int)SamplesPerPixel); i++)
			printf("%u", BitsPerSample[i] );
		printf("\n");
	   printf("Compression(259) = %u\n", Compression );
	   printf("DateTime(306) = %s\n", DateTime);
	   printf("ImageLength(257) = %ld\n", ImageLength );
	   printf("ImageWidth(256) = %ld\n", ImageWidth );
	   printf("MaxSampleValue(281) = %u\n", MaxSampleValue );
	   printf("MinSampleValue(280) = %u\n", MinSampleValue );
	   printf("NewSubfileType(254) = %lu\n", NewSubfileType );
	   printf("Orientation(274) = %u\n", Orientation );
	   printf("PhotometricInterpretation(262) = %u\n",
	  		PhotometricInterpretation);
	   printf("PlanarConfiguration(284) = %u\n", PlanarConfiguration );
	   printf("Predictor(317) = %u\n", Predictor );
	   printf("ResolutionUnit(296) = %u\n", ResolutionUnit );
	   printf("RowPerStrip(278) = %lu\n", RowsPerStrip );
	   printf("SampleFormat(339) = %d\n", SampleFormat );
	   printf("SamplesPerPixel(277) = %u\n", SamplesPerPixel );
	   if( nstripcounts > 0 ) {
		printf("StripByteCounts length = %d\n", nstripcounts);
		printf("StripByteCounts(279) =\n   ");
		for( i=0; i<nstripcounts; i++) {
			printf("%12lu ", StripByteCounts[i]);
			if( ((i+1)/5)*5 == (i+1) ) printf("\n   "); }
		printf("\n");
	   } else printf( "StripByteCounts missing (BAD!)\n");

	   printf("StripOffsets length = %d\n", nstrips);
	   printf("StripOffsets(273) =\n   ");
	   for( i=0; i<nstrips; i++) {
		printf("%12lu ", StripOffsets[i]);
		if( ((i+1)/5)*5 == (i+1) ) printf("\n   "); }
	   printf("\n");
	   printf("SubfileType(255) = %u\n", SubfileType );
	   printf("XResolution(282) = %lu, %lu\n", 
			XResolution[0], XResolution[1] );
	   printf("YResolution(283) = %lu, %lu\n", 
			YResolution[0], YResolution[1] );
	   printf("XPosition(286) = %lu, %lu\n", 
			XPosition[0], XPosition[1] );
	   printf("YPosition(287) = %lu, %lu\n", 
			YPosition[0], YPosition[1] );

	   printf("\n Current IFD (type=3, length=1 values word swapped): \n");
	   printf("     #    tag   type       length      value(dec, hex) \n");
	   for( i=0; i<nIFD; i++)
		printf("%6d %6d %6d %12lu %10lu %10lX\n",
		      i, ifd1[i].tag, ifd1[i].type, ifd1[i].length,
			ifd1[i].value, ifd1[i].value );  }
#endif
		return;

}  /* end tlist */


/*--------------------- topen -------------------*/
/*
   open a TIFF file for reading and check that the
   header is a valid TIFF header

   filename  = pointer to string with filename

   return 1 for success and </=0 for failure
*/

int topen( const char *filename )
{
   void tifferr( const char* );
   void tsetByteOrder();

   int iread;
   int tread( void *bufptr, int size, int n, long offset);
   int treadIFD( unsigned long IFDoffset );

	tsetByteOrder();
	tfp = fopen( filename, "rb" );
	if( tfp == NULL ) {
		tifferr( "topen() bad filename." );
		return( -1 ); }

	iread  = fread( &FileByteOrder, sizeof(short), 1, tfp);

	if( (iread != 1) || 
	   !( (FileByteOrder==tLittleEndian) || 
                      (FileByteOrder==tBigEndian) ) ) {
		tifferr( "Not a TIFF file in topen().");
		return( -2 ); }

	tread( &Version, 2, 1, 2L );
	if( Version != 42 ) {    		/* the TIFF magic number */
		tifferr( "Not a TIFF file in topen().");
		return( -3 ); }

	/* get location of IFD  and read it */

	tread( &CurrentIFD, 4, 1, 4L );
	iread = treadIFD( CurrentIFD );

	if( iread != 1 ) {
		tifferr( "topen() can't read IFD.");
		return( iread ); }

	return( 1 );

}  /* end topen() */


/*--------------------- topenFlaot -------------------*/
/*
   open a extended TIFF file for reading as a floating
   point image and check that the header is a valid TIFF header

   filename  = pointer to string with filename

   return 1 for success and </=0 for failure
*/

int topenFloat( const char *filename )
{
   void tifferr( const char* );
   void tsetByteOrder();

   int iread;
   int tread( void *bufptr, int size, int n, long offset);
   int treadIFD( unsigned long IFDoffset );

	tsetByteOrder();
	tfp = fopen( filename, "rb" );
	if( tfp == NULL ) {
		tifferr( "topenFloat() bad filename." );
		return( -1 ); }

	iread  = fread( &FileByteOrder, sizeof(short), 1, tfp);

	if( (iread != 1) || 
	   !( (FileByteOrder==tLittleEndian) || 
                      (FileByteOrder==tBigEndian) ) ) {
		tifferr( "Not a TIFF file in topenFloat().");
		return( -2 ); }

	tread( &Version, 2, 1, 2L );
	if( Version != 42 ) {    		/* the TIFF magic number */
		tifferr( "Not a TIFF file in topenFloat().");
		return( -3 ); }

	/* get location of IFD  and read it */

	tread( &CurrentIFD, 4, 1, 4L );
	iread = treadIFD( CurrentIFD );

	if( iread != 1 ) {
		tifferr( "topenFloat() can't read first IFD.");
		return( iread ); }

	/* find the first floating point image   */
	while( ( SampleFormat != 3 ) &&
	       ( BitsPerSample[0] != 8*sizeof(float) ) &&
	       ( NextIFD > 0 ) ) {
		if( treadIFD( NextIFD ) < 1 ) {
			tifferr("Cannot read IFD in treadFloatPix.");
			return( -4 );
		}
	}

	if( ( SampleFormat != 3 ) ||
	       ( BitsPerSample[0] != 8*sizeof(float) ) ) return( -5 );

	return( 1 );

}  /* end topenFloat() */


/*--------------------- tread -------------------*/
/*
	this routine mimics the ANSI fread function
	however it fixes the byte order if required

	read n elements of length size into the buffer
	pointed to by bufptr from the current TIFF file
	at offset position (in bytes)

	if offset < 0 then do sequential reads from
	current position

	return 1 for success and </=0 for failure
	
	for internal use only
	DO NOT CALL FROM MAIN PROGRAM
*/

int tread( void *bufptr, int size, int n, long offset )
{
   void tifferr( const char* );
   int iread, i, j;
   unsigned char *temp, *ctmp;

	if( offset > 0 ) if( fseek( tfp, offset, SEEK_SET ) != 0 ) {
		tifferr("Bad file seek in tread().");
		return( -1 ); }
	
	iread = fread( bufptr, size, n, tfp);
	if( iread != n) {
		tifferr("Bad file read in tread().");
		return( -2 ); }

	if ( ( FileByteOrder != CompByteOrder ) && ( size > 1 ) )
	  {	ctmp = (unsigned char*) malloc( size );
	  	temp = (unsigned char*) bufptr;
		for( i=0; i<n; i++) {
			for( j=0; j<size; j++) ctmp[j] = temp[i*size + j];
			for( j=0; j<size; j++) temp[i*size + j] = 
					ctmp[size-1-j];
		}
		free(ctmp);
	  }

	return( 1 );
		
}  /* end tread */

/*--------------------- treadFloatPix -------------------*/
/*
	read in the current image file opend by topen()
	(i.e. topen() must be called before this routine)
	--- assumed to be a floating point image made by
		tcreateFloatPixFile() 

	this routine returns +1 for success 
		and </= 0 for failure

	NOTE: this routine makes some non-standard assumptions
	about the image file (i.e. floating point pixel data
	in a specific order - as created by tcreatFloatPix)

	pix[][]      = float** array that will get the image
			indexed as pix[ix][iy]
			must be allocated before calling this routine
	nx, ny       = actual size of pix[][]
	*npix        = 1 for real and 2 for complex (real+imag)
	param[64]    = float array to get parameter string
			param[0] = 1 for real and 2 for complex
	datetime[20] = char array to get date and time

	return +1 for success, 0 or negative for fatal errors
		and >+1 for non-fatal errors

   started 25-mar-1996 ejk

*/
int treadFloatPix( float** pix, long nx, long ny, int *npix,
	char datetime[], float param[] )
{
   void tifferr( const char* );
   
   float *FStripBuf;
   uslong ix, iy, iyi, iyf, nbytes;
   int istrip, BytesPerPixel, status;
   int treadIFD( unsigned long IFDoffset );

	/* check that this is a floating point image */

	if( (BitsPerSample[0] != 8*sizeof(float)) || (SampleFormat!=3) ) {
	   tifferr("treadFloatPix cannot find floating point image.");
	   return( -2 );
	}

	if( nstrips != nstripcounts ) {
	   tifferr( "nstrips and nstripcounts do not agree in treadFloatPix.");
	   /* return( -3 );  leave return for rigorous TIFF */
	}

	if( SamplesPerPixel != 1 ) {
	   tifferr("treadFloatPix cannot handle more than 1 sample/pixel.");
	   return( -4 ); }

	if( (nx < (signed) ImageWidth ) || ( ny < (signed) ImageLength ) ) {
	   tifferr("Pix array size too small in treadFloatPix.");
	   return( -5 ); }
	   
	BytesPerPixel = sizeof(float);

	/* if the image is properly broken into strips then
	    use them, otherwise read it in one row at a time
	    to minimize the size of the intermediate storage
	    buffer so that this will work on DOS type machines too */

	if( nstrips > 1 )
	   nbytes = RowsPerStrip * ImageWidth;
	else
	   nbytes = ImageWidth;

	FStripBuf = (float*) calloc( (size_t) nbytes, BytesPerPixel );
	if( FStripBuf == NULL ) {
	   tifferr( "treadFloatPix unable to allocate pix temporary memory.");
	   return( -6 ); }
	nbytes = nbytes * BytesPerPixel;

	if( nstrips > 1 ) {    /* use stipes if present  */
	   iyi = 0;
	   for( istrip=0; istrip<nstrips; istrip++ ) {
		if( tread( FStripBuf, BytesPerPixel, 
			(int) (StripByteCounts[istrip]/BytesPerPixel),
			(long) StripOffsets[istrip]  ) != 1) {
 		 		free( FStripBuf );
		 		tifferr("Bad file read in treadFloatPix.");
		 		return( -7 ); }
	      iyf = iyi + RowsPerStrip;
	      if( iyf > ImageLength ) iyf = ImageLength;
	      for( iy=iyi; iy<iyf; iy++) {
		for( ix=0; ix<ImageWidth; ix++)
		   pix[ix][iy] = FStripBuf[ix + (iy-iyi)*ImageWidth ];
	      }
	      iyi = iyf;
	   } /* end for istrip loop */

	} else {		/* otherwise break it up into rows */
	   if( fseek( tfp, StripOffsets[0], SEEK_SET ) != 0 ) {
	      	free( FStripBuf );
	      	tifferr("Bad file seek in treadFloatPix.");
	      	return( -8 ); }
	   for( iy=0; iy<ImageLength; iy++ ) {
	      if( tread( FStripBuf, BytesPerPixel,   /* -1 = don't seek */
			(int)(nbytes/BytesPerPixel), -1L ) != 1 ) {
		      free( FStripBuf );
		      tifferr("Bad file read in treadFloatPix");
		      return( -9 ); }
	     for( ix=0; ix<ImageWidth; ix++) pix[ix][iy] = FStripBuf[ix];
	   }
	} /* end if nstrip */

	free( FStripBuf );

   /*  read the parameter array from the 3rd IFD */

	if( NextIFD > 0 ) {
	   if( status=treadIFD( NextIFD ) < 1 ) {
		tifferr("Cannot read parameter IFD in treadFloatPix.");
		return( +2 );  /* wrong but not fatal */
	   }
	} else {
		tifferr("No parameters available in treadFloatPix.");
		return( +3 ); /* wrong but non-fatal */
	}

	if( (BitsPerSample[0] != 8*sizeof(float)) || (SampleFormat!=3) ) {
	   tifferr("treadFloatPix cannot read parameters.");
	   return( +4 ); /* wrong but not fatal */
	}

	strncpy( datetime, DateTime, 20 );

	if( StripByteCounts[0]/BytesPerPixel <= 64 ) {
	   if( tread( param, BytesPerPixel, 
		(int) (StripByteCounts[0]/BytesPerPixel),
		(long) StripOffsets[0]  ) != 1) {
	 		tifferr("Bad file read in treadFloatPix.");
	 		return( -10 ); }
	    *npix = (int) ( param[0] + 0.5 );
	} else return( +5 );  /* wrong but not fatal */

   /*  all done */

   return( 1 );  /* success */

}  /* end treadFloatPix */

/*--------------------- treadIFD -------------------*/
/*
	read the tiff directory structure (IFD)
	this is called by topen()

	return 1 for success and </=0 for failure

	for internal use only
	DO NOT CALL FROM MAIN PROGRAM
*/
int treadIFD( unsigned long IFDoffset )
{
   void tifferr( const char* );

   char s[100];		/* s = temp storage for error text */
   char s2[100];	/* s2 = another temp storage for error text */
   unsigned short n1, *stemp;
   int i, j, n, ierr;
   unsigned long ltemp;
   int tread( void *bufptr, int size, int n, long offset );

	tread( &n1, 2, 1, IFDoffset);
	n = (int) n1;

	if( nIFD <= 0 )
		ifd1  = (IFD*) malloc( (size_t) n * sizeof(IFD));
	else if( n > nIFD )
		ifd1  = (IFD*) realloc( ifd1, n * sizeof(IFD));

	if( ifd1 == NULL ){
		tifferr( "Unable to allocate memory in treadIFD(1).");
		return( -1 ); }

	nIFD = n;

	for( i=0; i<nIFD; i++ ){
		tread( &ifd1[i].tag   , 2, 1, -1L );
		tread( &ifd1[i].type  , 2, 1, -1L );
		tread( &ifd1[i].length, 4, 1, -1L );
		/* to get word order right */
		if( (ifd1[i].type == 3) && (ifd1[i].length <= 2 ) )
			tread( &ifd1[i].value, 2, 2, -1L );
		else 
			tread( &ifd1[i].value, 4, 1, -1L ); 
	}

	stemp = (unsigned short*) &ltemp;  /* for data type 3 separations */

	tread( &NextIFD, 4, 1, -1L);  /* point to next IFD */

	/* set defaults */

	BitsPerSample[0] = 1;
	BitsPerSample[1] = 0;
	BitsPerSample[2] = 0;
	Compression = 1;    /* no compression */
	strcpy( DateTime, "No-Date-Given");
	FillOrder = 1;
	MaxSampleValue = 1; /* tag no longer recommended */
	MinSampleValue = 0; /* tag no longer recommended */
	NewSubFileType = 0;
	nstrips = 0;
	nstripcounts = 0;
	Orientation = 1; /* tags not recommended */
	PhotometricInterpretation = 1;  /* 0=black, 1+=white */
	PlanarConfiguration = 1;
	Predictor = 1;
	ResolutionUnit = 2;
	RowsPerStrip = 0; /* default is really to infinite- see below */
	SamplesPerPixel = 1;
	SampleFormat = 1;  /* default to unsigned integer */
	SubfileType = 0; /* tag no longer recommended */
	XPosition[0] = 0;
	XPosition[1] = 0;
	YPosition[0] = 0;
	YPosition[1] = 0;
	XResolution[0] = 0;
	XResolution[1] = 0;
	YResolution[0] = 0;
	YResolution[1] = 0;

	for( i=0; i<nIFD; i++ ) { ierr = 0;
	   switch ( ifd1[i].tag ) {
	   case 254:                       /* NewSubfileType tag */
		if( (ifd1[i].type != 4) || 
		    (ifd1[i].length != 1 ) ) {
			ierr = -254; break; }
		NewSubfileType = ifd1[i].value;
		break;
	   case 255:                       /* SubfileType tag */
		if( (ifd1[i].type != 3) || (ifd1[i].length != 1 ) )
			{ ierr = -255; break; }
		ltemp = ifd1[i].value;
		SubfileType =  stemp[0];
		break;
	   case 256:                       /* ImageWidth tag */
		if((ifd1[i].length != 1 ) ) {
			ierr = -256; break; }
		if( ifd1[i].type == 3 ) {  /* type = short */
		   ltemp = ifd1[i].value;
		   ImageWidth =  stemp[0]; }
		else if( ifd1[i].type == 4 )   /* type = long */
		   ImageWidth = ifd1[i].value;
		else { ierr = -256; break; }
		break;
	   case 257:                       /* ImageLength tag */
		if((ifd1[i].length != 1 ) ) {
			ierr = -257; break; }
		if( ifd1[i].type == 3 ) {  /* type = short */
		   ltemp = ifd1[i].value;
		   ImageLength =  stemp[0]; }
		else if( ifd1[i].type == 4 )   /* type = long */
		   ImageLength = ifd1[i].value;
		else { ierr = -257; break; }
		break;
	   case 258:                       /* BitsPerSample tag */
		if( (ifd1[i].type != 3) || 
		    (ifd1[i].length < 1 ) || (ifd1[i].length > 3) ) {
			ierr = -258; 
			sprintf(s, "BitsPerSample type=%d, length=%d ??\n",
			       ifd1[i].type, ifd1[i].length); break; }
		if( ifd1[i].length == 1 ) {
		   ltemp = ifd1[i].value;
		   BitsPerSample[0] = stemp[0];
		   if( (BitsPerSample[0] < 1)
		     || (BitsPerSample[0] > 32 ) )  {
		      sprintf(s2, "Bad BitsPerSample = %d, will try 8.",
			      BitsPerSample[0] );
		      tifferr( s2 );
		      BitsPerSample[0] = 8;  }
		   }
		else if( ifd1[i].length == 2 ) {
		   ltemp = ifd1[i].value;
		   BitsPerSample[0] = stemp[0];
		   BitsPerSample[1] = stemp[1];  }
		else if( ifd1[i].length == 3 )
		   tread( BitsPerSample, 2, 3, ifd1[i].value );
		break;
	   case 259:                       /* Compression tag */
		if( (ifd1[i].type != 3) || 
		    (ifd1[i].length != 1 ) ) {
			ierr = -259;
			sprintf(s,"Compression type=%d, length=%d ??\n",
			      ifd1[i].type, ifd1[i].length); break; }
		ltemp = ifd1[i].value;
		Compression = stemp[0];
		break;
	   case 262:             /* PhotometricInterpretation tag */
		if( (ifd1[i].type != 3) || 
		    (ifd1[i].length != 1 ) ) {
			ierr = -262; break; }
		ltemp = ifd1[i].value;
		PhotometricInterpretation = stemp[0];
		break;
	   case 266:			/* FillOrder tag */
		if( (ifd1[i].type != 3) || 
		    (ifd1[i].length != 1 ) ) {
			ierr = -266; break; }
		ltemp = ifd1[i].value;
		FillOrder = stemp[0];
		break;
	   case 273:                       /* StripOffsets tag */
		n = (int) nstrips;
		nstrips = (int) ifd1[i].length;
		if( StripOffsets == NULL ) {

		   StripOffsets = (uslong*) malloc( 
				(size_t) nstrips * sizeof(long) );
		   stempo = (usshort*) malloc(
				(size_t) nstrips * sizeof(short));
		} else if( nstrips > n ) {
		   StripOffsets = (uslong*) realloc( StripOffsets,
			    nstrips * sizeof(long));
		   stempo = (usshort*) realloc( stempo,
			    nstrips * sizeof(short)); }
		if( (StripOffsets == NULL) || ( stempo == NULL) ) {
		   tifferr( "treadIFD is unable to allocate Offsets"
			   " memory(2).");
		   return( -2 ); }

		/*  check if IFDvalues is a value or an offset */
		if( (ifd1[i].length * sizes[ ifd1[i].type ]) > 4 ) {
		   if( ifd1[i].type == 3 ) {  /* type = short */
		      tread( stempo, 2, nstrips, ifd1[i].value );
		      for( j=nstrips-1; j>=0; j--)
			   StripOffsets[j] = (uslong) stempo[j]; }
		   else if( ifd1[i].type == 4 )   /* type = long */
		      tread( StripOffsets, 4, nstrips, ifd1[i].value );
		   else { ierr = -273; break; } }
		else {
		   if( ifd1[i].type == 3 ) {  /* type = short */
			ltemp = ifd1[i].value;
			StripOffsets[0] = stemp[0]; }
		   else if( ifd1[i].type == 4 )   /* type = long */
			StripOffsets[0] = ifd1[i].value;
		   else { ierr = -273; break; } }
		break;
	   case 274:                       /* Orientation tag */
		if( (ifd1[i].type != 3) || (ifd1[i].length != 1 ) )
			{ ierr = -274; break; }
		ltemp = ifd1[i].value;
		Orientation =  stemp[0];
		break;
	   case 277:                       /* SamplesPerPixel tag */
		if( (ifd1[i].type != 3) || 
		    (ifd1[i].length != 1 ) ) {
			ierr = -277; break; }
		ltemp = ifd1[i].value;
		SamplesPerPixel = stemp[0];
		break;
	   case 278:                       /* RowsPerStrip tag */
		if( ifd1[i].length != 1 ) {
			ierr = -278;
			sprintf(s,"RowsPerStrip type=%d, length=%d ??\n",
			      ifd1[i].type, ifd1[i].length); break; }
		if( ifd1[i].type == 3 ) {  /* type = short */
		   ltemp = ifd1[i].value;
		   RowsPerStrip =  stemp[0]; }
		else if( ifd1[i].type == 4 )   /* type = long */
		   RowsPerStrip = ifd1[i].value;
		else { ierr = -278; break; }
		break;
	   case 279:                       /* StripByteCounts tag */
		n = nstripcounts;
		nstripcounts = (int) ifd1[i].length;
		if( StripByteCounts == NULL ) {
		   StripByteCounts = (uslong*) malloc( 
				(size_t) nstrips * sizeof(long));
		   stempc = (usshort*) malloc( 
				(size_t) nstrips * sizeof(short)); }
		else if( nstripcounts > n ) {
		   StripByteCounts = 
			   (uslong*) realloc( StripByteCounts,
				nstripcounts * sizeof(long));
		   stempc = (usshort*) realloc( stempc, 
				nstripcounts * sizeof(short)); }
		if( (StripByteCounts == NULL) || ( stempc == NULL) ) {
		   tifferr( "treadIFD is unable to allocate memory(3).");
		   return( -3 ); }

		/*  check if IFDvalues is a value or an offset */
		if( (ifd1[i].length * sizes[ ifd1[i].type ]) > 4 ) {
		   if( ifd1[i].type == 3 ) {  /* type = short */
		      tread( stempc, 2, nstrips, ifd1[i].value );
		      for( j=nstrips-1; j>=0; j--)
			   StripByteCounts[j] = (uslong) stempc[j]; }
		   else if( ifd1[i].type == 4 )   /* type = long */
		      tread( StripByteCounts, 4, nstrips, ifd1[i].value );
		   else { ierr = -279; break; } }
		else {
		   if( ifd1[i].type == 3 ) {  /* type = short */
			ltemp = ifd1[i].value;
			StripByteCounts[0] = stemp[0]; }
		   else if( ifd1[i].type == 4 )   /* type = long */
			StripByteCounts[0] = ifd1[i].value;
		   else { ierr = -279; break; } }
		break;
	   case 280:                       /* MinSampleValue tag */
		if( (ifd1[i].type != 3) || (ifd1[i].length != 1 ) )
			{ ierr = -280; break; }
		ltemp = ifd1[i].value;
		MinSampleValue =  stemp[0];
		break;
	   case 281:                       /* MaxSampleValue tag */
		if( (ifd1[i].type != 3) || (ifd1[i].length != 1 ) )
			{ ierr = -281; break; }
		ltemp = ifd1[i].value;
		MaxSampleValue =  stemp[0];
		break;
	   case 282:                       /* XResolution tag */
		if( (ifd1[i].length != 1 ) ||
		    (ifd1[i].type != 5 )   ) {
			ierr = -282; break; }
		tread( &XResolution[0], 4, 1, ifd1[i].value );
		tread( &XResolution[1], 4, 1, -1L );
		break;
	   case 283:                       /* YResolution tag */
		if( (ifd1[i].length != 1 ) ||
		    (ifd1[i].type != 5 )   ) {
			ierr = -283; break; }
		tread( &YResolution[0], 4, 1, ifd1[i].value );
		tread( &YResolution[1], 4, 1, -1L );
		break;
	   case 284:                       /* PlanarConfiguration tag */
		if( (ifd1[i].type != 3) || 
		    (ifd1[i].length != 1 ) ) {
			ierr = -284; 
			sprintf(s,"PlanarConfiguration type=%d, "
				"length=%d ??\n", ifd1[i].type,
				ifd1[i].length); break; }
		ltemp = ifd1[i].value;
		PlanarConfiguration = stemp[0];
		break;
	   case 286:                       /* XPosition tag */
		if( (ifd1[i].length != 1 ) ||
		    (ifd1[i].type != 5 )   ) {
			ierr = -286; break; }
		tread( &XPosition[0], 4, 1, ifd1[i].value );
		tread( &XPosition[1], 4, 1, -1L );
		break;
	   case 287:                       /* YPosition tag */
		if( (ifd1[i].length != 1 ) ||
		    (ifd1[i].type != 5 )   ) {
			ierr = -287; break; }
		tread( &YPosition[0], 4, 1, ifd1[i].value );
		tread( &YPosition[1], 4, 1, -1L );
		break;
	   case 296:             /* ResolutionUnit tag */
		if( (ifd1[i].type != 3) || 
		    (ifd1[i].length != 1 ) ) {
			ierr = -296; break; }
		ltemp = ifd1[i].value;
		ResolutionUnit = stemp[0];
		break;
	   case 306:             /* (ASCII) DateTime tag */
		if( ifd1[i].type != 2 ) {
			ierr = -296; break; }
		tread( DateTime, 1, 20, ifd1[i].value);
		break;
	   case 317:                       /* Predictor tag */
		if( (ifd1[i].type != 3) || 
		    (ifd1[i].length != 1 ) ) {
			ierr = -317; break; }
		ltemp = ifd1[i].value;
		Predictor = stemp[0];
		break;
	   case 339:             /* SampleFormat tag */
		if( (ifd1[i].type != 3) || 
		    (ifd1[i].length != 1 ) ) {
			ierr = -339; break; }
		ltemp = ifd1[i].value;
		SampleFormat = stemp[0];
		break;
	   default:
		ierr = -1000;
		sprintf(s, "warning, treadIFD unknown TIFF tag # %d",
				ifd1[i].tag);
		break;
	   }
	  if ( ierr != 0 ) {
		tifferr( s );
		if( ierr > -999 ) return( ierr ); }

	}  /* end for(i=0; i<nIFD; i++) */

	/* if RowsPerStrip is missing make a guess */
	if( RowsPerStrip == 0 ) RowsPerStrip = ImageLength;

	return( 1 );

}  /* end treadIFD */


/*--------------------- treadPix -------------------*/
/*
	read in the current image file opend by topen
	(i.e. topen() must be called before this routine)

	this routine returns +1 for success 
		and </= 0 for failure

	WARNING: this routine assumes that the image
	storage array has already been allocated

	NOTE: this routine ignores many TIFF options!
		it thinks everything is 8 bit greyscale

	pix[][] = long** array that will get the image
		indexed as pix[iy][ix]

*/
int treadPix( long** pix )
{
   void tifferr( const char* );
   
   unsigned char *StripBuf;
   unsigned short *SStripBuf;
   uslong ix, iy, iyi, iyf, nbytes;
   int istrip, BytesPerPixel;

	if( nstrips != nstripcounts ) {
		tifferr( "nstrips and nstripcounts do not agree in treadPix.");
		/* return( -1 );  leave return for rigorous TIFF */
	}

	if( SamplesPerPixel != 1 ) {
	   tifferr("treadPix cannot handle more than 1 sample/pixel.");
	   return( -2 ); }
	   
	if( BitsPerSample[0] == 8 ) BytesPerPixel = 1;
	else  if( BitsPerSample[0] == 16 ) BytesPerPixel = 2;
	else {
	   tifferr("treadPix cannot handle anything except 1"
		   " or 2 bytes/pixel.");
	   return( -3 );
	}

	/* if the image is properly broken into strips then
	    use them, otherwise read it in one row at a time
	    to minimize the size of the intermediate storage
	    buffer so that this will work on DOS type machines too */

	if( nstrips > 1 )
	   nbytes = RowsPerStrip * ImageWidth;
	else
	   nbytes = ImageWidth;

	SStripBuf = (unsigned short*) calloc( (size_t) nbytes, BytesPerPixel );
	StripBuf = (unsigned char*) SStripBuf;
	if( (StripBuf == NULL) || (SStripBuf == NULL) ) {
	   tifferr( "treadPix unable to allocate pix temporary memory.");
	   return( -4 ); }
	nbytes = nbytes * BytesPerPixel;

	if( nstrips > 1 ) {    /* use stipes if present  */
	   iyi = 0;
	   for( istrip=0; istrip<nstrips; istrip++ ) {
		if( tread( StripBuf, BytesPerPixel, 
			(int) (StripByteCounts[istrip]/BytesPerPixel),
			(long) StripOffsets[istrip]  ) != 1) {
 		 		free( StripBuf );
		 		tifferr("Bad file read in treadPix.");
		 		return( -5 ); }
	      iyf = iyi + RowsPerStrip;
	      if( iyf > ImageLength ) iyf = ImageLength;
	      for( iy=iyi; iy<iyf; iy++) {
	      	  if( BytesPerPixel == 1 ) {
			for( ix=0; ix<ImageWidth; ix++)
			pix[iy][ix] = (long) StripBuf[ix 
					+ (iy-iyi)*ImageWidth ];
		  } else if( BytesPerPixel == 2 ) {
			for( ix=0; ix<ImageWidth; ix++)
			pix[iy][ix] = (long) SStripBuf[ix 
					+ (iy-iyi)*ImageWidth ];
		  }
	      }
	      iyi = iyf;
	   } /* end for istrip loop */

	} else {		/* otherwise break it up into rows */
	   if( fseek( tfp, StripOffsets[0], SEEK_SET ) != 0 ) {
	      	free( StripBuf );
	      	tifferr("Bad file seek in treadPix.");
	      	return( -6 ); }
	   for( iy=0; iy<ImageLength; iy++ ) {
	      if( tread( StripBuf, BytesPerPixel,   /* -1 = don't seek */
			(int)(nbytes/BytesPerPixel), -1L ) != 1 ) {
		      free( StripBuf );
		      tifferr("Bad file read in treadPix");
		      return( -7 ); }
	      	  if( BytesPerPixel == 1 ) {
			for( ix=0; ix<ImageWidth; ix++)
			pix[iy][ix] = (long) StripBuf[ix];
		  } else if( BytesPerPixel == 2 ) {
			for( ix=0; ix<ImageWidth; ix++)
			pix[iy][ix] = (long) SStripBuf[ix];

		  }
	   }
	} /* end if nstrip */

	free( StripBuf );

	if( ( PhotometricInterpretation == 0 ) &&
	    ( BitsPerSample[0] == 8 ) ) {
	   for( iy=0; iy<ImageLength; iy++)
		for( ix=0; ix<ImageWidth; ix++)
		       pix[iy][ix] = 255 - pix[iy][ix];
	}  /* end if photometric */

   return( 1 );

}  /* end treadPix */

/*--------------------- tsetByteOrder -------------------*/
/*
   determine the byte order (big-endian or little-endian)
   of the computer this is running on
*/
void tsetByteOrder()
{
   unsigned short us;
   unsigned char *uc;

   uc = (unsigned char*) &us;
   uc[0] = 1;
   uc[1] = 0;
   if( us > 1 ) CompByteOrder = tBigEndian;
          else  CompByteOrder = tLittleEndian;

}  /* end tsetByteOrder() */

/*--------------------- tsize -------------------*/
/*
   nx, ny will get the size of the image in pixels
      of the currently open file (see topen)
      (i.e. topen must be called before this routine)
   nbits[3] will get the number of bits/sample
   samples will get the number of samples/pixel
*/
void tsize( long *nx, long *ny, int *nbits, int *samples )
{
   int i;
	*nx = ImageWidth;   *ny = ImageLength;
	for( i=0; i<3; i++) nbits[i] = (int) BitsPerSample[i];
	*samples = (int) SamplesPerPixel;

}  /* end tsize() */




	
