/*  		*** tiffsubs.h ***

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

   Header file to go with tiffsubs.c

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

*/                                       	

/*--------------------- tclose -------------------*/
/*
	close a TIFF file opened by topen
	and free malloc'ed data storage

	return 1 for success and </=0 for failure
*/

int tclose( );

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
		long psub, double scale, double aspect );

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
*/

int tcreateFloatPixFile( const char *file, float** pix,
	long nx, long ny, int npix, float param[] );

/*--------------------- tFloatTest -------------------*/
/*
   Test that a file is an extended TIFF file and
   has a floating point image (not an exhaustive test)

   filename  = pointer to string with filename

   return 1 for success and </=0 for failure
*/

int tFloatTest( const char *filename );

/*--------------------- tifferr ------------------*/
/*
    common error handler - print messages
    only for internal use 
    DO NOT CALL FROM MAIN PROGRAM
*/

void tifferr( const char* error_text );


/*--------------------- tifftest -------------------*/
/*
   check if a file has a valid TIFF header

   return a 1 if it is a valid TIFF header, </=0 otherwise
   filename  = pointer to string with filename
*/
int tifftest( const char *filename );

/*--------------------- tinter() ----------------------------*/
/*
  Bilinear interpolation from pix array

  pix[][]  = real input array with the image
  nx,ny    = integer dimension of the image
  x,y      = real coord of image point to interpolate at
*/
float tinter( float **pix, long nx, long ny, double x, double y );

/*--------------------- tlist -------------------*/
/*
	this routine list a few key parameters on the screen
	about the currently open file
*/

void tlist( );

/*--------------------- topen -------------------*/
/*
   open a TIFF file for reading and check that the
   header is a valid TIFF header

   filename  = pointer to string with filename

   return 1 for success and </=0 for failure
*/

int topen( const char *filename );

/*--------------------- topenFlaot -------------------*/
/*
   open a extended TIFF file for reading as a floating
   point image and check that the header is a valid TIFF header

   filename  = pointer to string with filename

   return 1 for success and </=0 for failure
*/

int topenFloat( const char *filename );

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

int tread( void *bufptr, int size, int n, long offset );

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
	char datetime[], float param[] );

/*--------------------- treadIFD -------------------*/
/*
	read the tiff directory structure (IFD)
	this is called by topen()

	return 1 for success and </=0 for failure

	for internal use only
	DO NOT CALL FROM MAIN PROGRAM
*/
int treadIFD( unsigned long IFDoffset );


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
int treadPix( long** pix );

/*--------------------- tsetByteOrder -------------------*/
/*
   determine the byte order (big-endian or little-endian)
   of the computer this is running on
*/
void tsetByteOrder();

/*--------------------- tsize -------------------*/
/*
   nx, ny will get the size of the image in pixels
      of the currently open file (see topen)
      (i.e. topen must be called before this routine)
   nbits[3] will get the number of bits/sample
   samples will get the number of samples/pixel
*/
void tsize( long *nx, long *ny, int *nbits, int *samples );



	
