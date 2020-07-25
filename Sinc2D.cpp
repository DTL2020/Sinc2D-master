// Sinc2D.cpp : Defines the entry point for the console application.
//


#include "stdafx.h"
#include "windows.h"
#include "commctrl.h"

#include "tiffio.h"

#include "math.h"
#include "stdlib.h"


BYTE *g_pImageBuffer = 0, *g_pFilteredImageBuffer = 0;

float *g_pfImageBuffer = 0, *g_pfFilteredImageBuffer = 0;

float *g_pfKernel = 0;

int iKernelSize;
int iDiv = 4;
int iTaps = 11;

int iWidth=2000;
int iHeight=2000;

//-----------------------------------------------------------------------------
// Msg: Display an error if needed
//-----------------------------------------------------------------------------
void Msg(TCHAR *szFormat, ...)
{
    printf(szFormat);
}

BOOL LoadTIFF8(void)
{
	TIFF *in_tiff;
	unsigned short bps, spp, cmpr;
	tstrip_t strip;
	uint32* bc;
	uint32 stripsize, rowspstrip;
	uint32 uiWidth, uiHeight;
	uint32 row, col;

	unsigned char *ucSrc_data;
	unsigned char *ucImageBuffer = (unsigned char*)g_pImageBuffer;

	// Open the src TIFF image
	if((in_tiff = TIFFOpen("in.tif", "r")) == NULL)
	{
		Msg(TEXT("Could not open src TIFF image"));
		return FALSE;
	}
	
	// Check that it is of a type that we support
	if((TIFFGetField(in_tiff, TIFFTAG_COMPRESSION, &cmpr) ==0 ) || (cmpr != COMPRESSION_NONE))
	{
		Msg(TEXT("Compressed TIFFs not supported"));
		TIFFClose(in_tiff);
		return FALSE;
	}

	if((TIFFGetField(in_tiff, TIFFTAG_BITSPERSAMPLE, &bps) == 0) || (bps != 8))
	{
		Msg(TEXT("Either undefined or unsupported number of bits per sample - must be 8"));
		TIFFClose(in_tiff);
		return FALSE;
	}

	if((TIFFGetField(in_tiff, TIFFTAG_SAMPLESPERPIXEL, &spp) == 0) || (spp != 3))
	{
		Msg(TEXT("Either undefined or unsupported number of samples per pixel - must be 3"));
		TIFFClose(in_tiff);
		return FALSE;
	}

	TIFFGetField(in_tiff, TIFFTAG_IMAGEWIDTH, &uiWidth);
	iWidth = uiWidth;

	TIFFGetField(in_tiff, TIFFTAG_IMAGELENGTH, &uiHeight);
	iHeight = uiHeight;

	ucSrc_data = (unsigned char*)malloc(uiWidth*uiHeight*3);

	TIFFGetField(in_tiff, TIFFTAG_STRIPBYTECOUNTS, &bc);
	stripsize = bc[0];
	// test if all strips are equal
	for (strip = 0; strip < TIFFNumberOfStrips(in_tiff); strip++) 
	{
		if (bc[strip] != stripsize) 
		{
			Msg(TEXT("Strip sizes unequal"));
//			return FALSE;
		}
	}

	TIFFGetField(in_tiff, TIFFTAG_ROWSPERSTRIP, &rowspstrip);

	// load from tiff to temp buffer
	for(row = 0; row < uiHeight; row+=rowspstrip)
	{
		TIFFReadRawStrip(in_tiff, row/rowspstrip, &ucSrc_data[row*uiWidth*3], stripsize);
	}

	// separate planes
	for(row=0; row < uiHeight; row++)
	{
		for(col=0; col < uiWidth; col++)
		{
			g_pImageBuffer[(row*uiWidth + col)] = ucSrc_data[(row*uiWidth+col)*3+0]; // R
			g_pImageBuffer[(uiWidth*uiHeight + row*uiWidth + col)] = ucSrc_data[(row*uiWidth+col)*3+1]; // G
			g_pImageBuffer[(uiWidth*uiHeight*2+ row*uiWidth + col)] = ucSrc_data[(row*uiWidth+col)*3+2]; // B
		}
	}

	free(ucSrc_data);
	TIFFClose(in_tiff);

	return TRUE;
}

int SaveTIFF8()
{
	TIFF *out_tiff;
	// Open the dst TIFF image
	if((out_tiff = TIFFOpen("out.tif", "w")) == NULL)
	{
		fprintf(stderr, "Could not open dst TIFF image\n");
		return(1);
	}


	TIFFSetField(out_tiff, TIFFTAG_IMAGEWIDTH, iWidth);
	TIFFSetField(out_tiff, TIFFTAG_BITSPERSAMPLE, 8);
	TIFFSetField(out_tiff, TIFFTAG_SAMPLESPERPIXEL, 3);
	TIFFSetField(out_tiff, TIFFTAG_ROWSPERSTRIP, 1);
	TIFFSetField(out_tiff, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
	TIFFSetField(out_tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
	TIFFSetField(out_tiff, TIFFTAG_FILLORDER, FILLORDER_MSB2LSB);
	TIFFSetField(out_tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);

	TIFFSetField(out_tiff, TIFFTAG_XRESOLUTION, 300);
	TIFFSetField(out_tiff, TIFFTAG_YRESOLUTION, 300);
	TIFFSetField(out_tiff, TIFFTAG_RESOLUTIONUNIT, RESUNIT_INCH);

	TIFFSetField(out_tiff, TIFFTAG_IMAGELENGTH, iHeight);

	unsigned char *ucDst_data  = (unsigned char*)malloc(iWidth*3);

	for (int row=0; row < iHeight; row++) // lines counter
		{

			for (int i = 0; i < iWidth; i++)
			{
				ucDst_data[i*3] = g_pFilteredImageBuffer[(row*iWidth + i)]; // R channel
				ucDst_data[i*3+1] = g_pFilteredImageBuffer[(iWidth*iHeight + row*iWidth + i)]; // G channel
				ucDst_data[i*3+2] = g_pFilteredImageBuffer[(iWidth*iHeight*2+ row*iWidth + i)]; // B channel
			}
		
			TIFFWriteRawStrip(out_tiff, row, ucDst_data, /*g_stripsize*/iWidth*3);
	}

	TIFFClose(out_tiff);

	free (ucDst_data);

	return 0;

}


void fill2DKernel(void)
{
	float fPi = 3.14159265358979;
	int i,j;

	for (i = 0; i < iKernelSize; i++)
	{
		for (j = 0; j < iKernelSize; j++)
		{

			float fDist= sqrtf((float(iKernelSize/2)-j)*(float(iKernelSize/2)-j) + (float(iKernelSize/2)-i)*(float(iKernelSize/2)-i));
			if (fDist <= (iKernelSize/2)*0.95) // FLOAT KRN TRUNCATION - significally changes edge effects !
			{
				float fArg = fPi*fDist/iDiv;

				if (fArg != 0)
				{
				  g_pfKernel[i*iKernelSize+j] = sinf(fArg)/fArg;
				}
				else
					g_pfKernel[i*iKernelSize+j] = 1.0f;
			}
			else
				g_pfKernel[i*iKernelSize+j] = 0.0f;
		}
	}

	// normalize to 1
	float fSum = 0.0f;
	for (i = 0; i < iKernelSize; i++)
	{
		for (j = 0; j < iKernelSize; j++)
		{
			fSum += g_pfKernel[i*iKernelSize+j];
		}
	}

	for (i = 0; i < iKernelSize; i++)
	{
		for (j = 0; j < iKernelSize; j++)
		{
			g_pfKernel[i*iKernelSize+j] /= fSum;
		}
	}


}

void KernelProc(void)
{
	int row,col, k_row,k_col;
	// R-channel
	for (row=0; row < iHeight; row++) // lines counter
	{
		for (col = 0; col < iWidth; col++)
		{
			g_pfImageBuffer[(row*iWidth + col)] = (float)g_pImageBuffer[(row*iWidth + col)]; // R channel
		}		
	}

	// set zero out
	ZeroMemory(g_pfFilteredImageBuffer,iWidth * iHeight * 3 * sizeof(float));

	// 2d convolution pass
	for (row=0+iKernelSize/2; row < iHeight-iKernelSize/2; row++) // lines counter
	{
		for (col = 0+iKernelSize/2; col < iWidth-iKernelSize/2; col++)
		{
			float fInpSample=g_pfImageBuffer[(row*iWidth + col)];
/*			if (fInpSample != 16.0f)
			{
				fInpSample *= 10;
			} */
			for (k_row = 0; k_row < iKernelSize; k_row++)
			{
				for (k_col = 0; k_col < iKernelSize; k_col++)
				{
					float fDist_sq=(k_row-(iKernelSize/2))*(k_row-(iKernelSize/2)) + (k_col-(iKernelSize/2))*(k_col-(iKernelSize/2));
//					if (fDist_sq > ((iKernelSize/2)^2)*0.9f) continue;

					float fKrn=g_pfKernel[(k_row)*iKernelSize + (k_col)];
					// out
					g_pfFilteredImageBuffer[(row-(iKernelSize/2)+k_row)*iWidth + (col-(iKernelSize/2)+k_col)] += fInpSample*fKrn;
				}
			}
		}
		printf("Just %d row of %d done. Be patient...\r",row,iHeight);
	}

    // R-channel
	for (row=0; row < iHeight; row++) // lines counter
	{
		for (col = 0; col < iWidth; col++)
		{
			unsigned char ucVal;
			float fVal = g_pfFilteredImageBuffer[(row*iWidth + col)]; // R channel
			//fVal*=1.2f;
			if (fVal > 255.0f) 
			{
				fVal= 255.0f;
			}
			if (fVal < 0.0f) 
			{
				fVal = 0.0f;
			}
			ucVal = (unsigned char) fVal;
			g_pFilteredImageBuffer[(row*iWidth + col)] = ucVal; // R channel
			g_pFilteredImageBuffer[(iWidth*iHeight + row*iWidth + col)] = ucVal; //temp to G 
			g_pFilteredImageBuffer[(iWidth*iHeight*2 + row*iWidth + col)] = ucVal; //temp to B 
		}		
	}
}

int main(int argc, char* argv[])
{
	printf("Hello World!\n");

	g_pImageBuffer = (BYTE*)malloc (iWidth * iHeight * 3);
	g_pFilteredImageBuffer = (BYTE*)malloc (iWidth * iHeight * 3);

	ZeroMemory(g_pFilteredImageBuffer,(iWidth * iHeight * 3));

	g_pfImageBuffer = (float*)malloc (iWidth * iHeight * 3 * sizeof(float));
	g_pfFilteredImageBuffer = (float*)malloc (iWidth * iHeight * 3 * sizeof(float));

	iKernelSize = iDiv*iTaps*2;

	g_pfKernel = (float*)malloc (iKernelSize * iKernelSize * sizeof(float));

	ZeroMemory(g_pfKernel, iKernelSize * iKernelSize * sizeof(float));

	fill2DKernel();

	if (LoadTIFF8() == FALSE)
		goto end;

	KernelProc();

	SaveTIFF8();
end:
	free (g_pImageBuffer);
    free (g_pFilteredImageBuffer);

	free (g_pfImageBuffer);
	free (g_pfFilteredImageBuffer);

	free (g_pfKernel);

	return 0;
}

