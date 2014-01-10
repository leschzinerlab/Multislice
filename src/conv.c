/*          *** conv.c ***


*/
#include <stdio.h>	/* ANSI C libraries */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "fft2dc.h"     /* FFT routines */
#include "memory.h"     /* memory allocation routines */
#include "slicelib.h"   /* misc. routines for multislice */
#include "tiffsubs.h"   /* file I/O routines in TIFF format */

#define NCMAX	132	/* max number of characters per line */
#define NPARAM	64	/* number of parameters */

#define MCPIX	0	/* coherent image mode */
#define MPCPIX	1	/* partial coherent mode */
#define MDIFFP	2	/* diffraction mode */
#define MPPPIX  3       /* phase plate mode */

void main()
{
	char filein[NCMAX], fileout[NCMAX], fileout_hed[NCMAX], fileout_img[NCMAX];
	const char version[] = "May 2008 (rjh)";
        char datetime[20];
	char h15[4];
	int ix, iy, nx, ny, ixmid, iymid, npix, ic,
	                     itens, mode, nsum, nbits[3], samples;
	int lcenter=0, lapert=0, count;
	long nxl, nyl;

	float kx,ky, ky2, k2, k2max, v0, wavlen, scale, pixc,
		ax, by, cz, rx, ry, pi, rmin, rmax, aimin, aimax,
		Cs,df,  alpha0, ddf, objlx, objly,
		tr, ti, wr, wi, dfa2, dfa2phi, dfa3, dfa3phii, pprad, ishift, oshift;
	float *param;
        float  **pixr, **pixi, **pix2r, **pix2i;
 	double sum, time, chi, chi1, chi2, phi, clog;
	FILE *hed, *img;
	size_t fret;
	int h1,h2,h4,h5,h6,h7,h8,h9,h10,h12,h13,h14,h61,h69, hnull;	
/*  echo version date */

	printf( "image version dated %s\n", version );
	
/*  get input file name */

	printf("Name of input file:\n");
        scanf("%s", filein );

	printf("Name of output file, no extension:\n");
	scanf("%s", fileout);
	strcpy(fileout_hed, fileout);
	strcpy(fileout_img, fileout);
	strcat ( fileout_hed, ".hed");
	strcat ( fileout_img, ".img");
	printf("%s\n", fileout_hed);
	printf("%s\n", fileout_img);
	
/*  open input file and check sizes etc. */

        if( topenFloat( filein ) != 1 ) {
        	printf( "Cannot open input file %s\n", filein );
        	exit( 0 );
        }
        tsize( &nxl, &nyl, nbits, &samples );
        nx = (int) nxl;
        ny = (int) nyl;
						
/*  read in specimen parameters and multislice result  */
        param = float1D( NPARAM, "param" );
	for( ix=0; ix<NPARAM; ix++) param[ix] = 0.0F;
        pixr = float2D( nx, ny, "pixr" );

	if( treadFloatPix( pixr, nxl, nyl, &npix, datetime, param )!= 1 ) {
		printf("Cannot read input file %s.\n", filein);
		exit(0);
	}
	tclose();
//
	rmax  = param[pRMAX];
	aimax = param[pIMAX];
	rmin  = param[pRMIN];
	aimin = param[pIMIN];
//
//
	hed = fopen(fileout_hed,"w");
	img = fopen(fileout_img,"w");
	
	count = 1;	
	      
	if (hed == NULL || img == NULL){
	       printf("Error opening output file \n");
		exit(0);	
	}
	
	if (npix != 2 ) {
//write header
		printf("Input file not complex.\n" );
		strncpy(h15, "REAL",4);
		hnull = 0;
		h1 = 1;		
		h2 = 0;
		h4 = 1;
		h5 = 7;
		h6 = 27;
		h7 = 2008;
		h8 = 12;
		h9 = 30;
		h10 = 0;
		h12 = (int) ( nxl*nyl);
		h13 = (int) nyl;
		h14 = (int) nxl;
		h61 = 1;
		h69 = 33686018;

//write image
		for (iy=0; iy <ny ; iy++){
			for (ix=0; ix < nx; ix++){
				fret = fwrite(&pixr[ix][iy], sizeof(float), (size_t) count, img);
			}
		}
	}
	else {
		printf("Input file complex.\n");
		strncpy (h15, "COMP",4);
		hnull = 0;
		h1 = 1;		
		h2 = 0;
		h4 = 1;
		h5 = 7;
		h6 = 27;
		h7 = 2008;
		h8 = 12;
		h9 = 30;
		h10 = 0;
		h12 = (int) ( (nxl/2)*nyl);
		h13 = (int) nyl;
		h14 = (int) nxl/2;
		h61 = 1;
		h69 = 33686018;
		nx = nx/2;
		for (iy=0; iy <ny ; iy++){
			for (ix=0; ix < nx; ix++){
				fret = fwrite(&pixr[ix][iy], sizeof(float), (size_t) count, img);
				fret = fwrite(&pixr[ix+nx][iy], sizeof(float), (size_t) count, img);
			}
		}
		printf("y %d x %d\n",h13,h14);

	}

		fret = fwrite(&h1,sizeof(float), (size_t) 1, hed);
		fret = fwrite(&h2,sizeof(float), (size_t) 1, hed);
		fret = fwrite(&hnull, sizeof(float), (size_t) 1, hed);
		fret = fwrite(&h4,sizeof(float), (size_t) 1, hed);
		fret = fwrite(&h5,sizeof(float), (size_t) 1, hed);
		fret = fwrite(&h6,sizeof(float), (size_t) 1, hed);
		fret = fwrite(&h7,sizeof(float), (size_t) 1, hed);
		fret = fwrite(&h8,sizeof(float), (size_t) 1, hed);
		fret = fwrite(&h9,sizeof(float), (size_t) 1, hed);
		fret = fwrite(&h10,sizeof(float), (size_t) 1, hed);
		fret = fwrite(&hnull,sizeof(float), (size_t) 1, hed);
		fret = fwrite(&h12,sizeof(float), (size_t) 1, hed);
		fret = fwrite(&h13,sizeof(float), (size_t) 1, hed);
		fret = fwrite(&h14,sizeof(float), (size_t) 1, hed);
		fret = fwrite(&h15,sizeof(float), (size_t) 4, hed);
		for (ic = 1; ic <=46 ;ic++){
			fret = fwrite(&hnull,sizeof(float), (size_t) 1, hed);
		}
		fret = fwrite(&h61,sizeof(float), (size_t) 1, hed);
		for (ic = 1; ic <=7 ;ic++){
			fret = fwrite(&hnull,sizeof(float), (size_t) 1, hed);
		}
		fret = fwrite(&h69,sizeof(float), (size_t) 1, hed);
		for (ic = 1; ic <=187 ;ic++){
			fret = fwrite(&hnull,sizeof(float), (size_t) 1, hed);
		}

	
	fclose(hed);
	fclose(img);
} /* end main() */
