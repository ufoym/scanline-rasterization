#pragma once

//*****************************************************************************//
// Written by Josiah Manson, September 2012
//*****************************************************************************//

#include "color_space.h"
#include <stdio.h>

template <int N>
void flip_img(Array2D< vect<N, double> > &img)
{
	for (int j = 0; j < img.s[1]/2; j++)
	{
		for (int i = 0; i < img.s[0]; i++)
		{
			swap(img(i,j), img(i,img.s[1]-1-j));
		}
	}
}

template <int N>
void save_pgm(Array2D< vect<N, double> > &img, string fn, int FILTER_W)
{
	//float color_padding = 0.1;

	FILE *f = fopen(fn.c_str(), "wb");
	if (f == 0)
		return;
	
	int border = FILTER_W >> 1;
	int extra = border << 1;

	fprintf(f, "P5 %d %d 255\n", img.s[0] - extra, img.s[1] - extra);

	for (int j = border; j < img.s[1] - border; j++)
	{
		for (int i = border; i < img.s[0] - border; i++)
		{
			for (int k = 0; k < N; k++)
			{
				//int val = (int)floor((img(i,j)[k]*(1-2*color_padding) + color_padding)*255 + .5);
				int val = RGB32F_to_sRGB8(img(i,j)[k]);
				if (val < 0)
					val = 0;
				else if (val > 255)
					val = 255;
				fprintf(f, "%c", val);
			}
		}
	}

	fclose(f);
}

template <int N>
void save_ppm(Array2D< vect<N, double> > &img, string fn, int FILTER_W)
{
	//float color_padding = 0.1;

	FILE *f = fopen(fn.c_str(), "wb");
	if (f == 0)
		return;
	
	int border = FILTER_W >> 1;
	int extra = border << 1;

	fprintf(f, "P6 %d %d 255\n", img.s[0] - extra, img.s[1] - extra);

	for (int j = border; j < img.s[1] - border; j++)
	{
		for (int i = border; i < img.s[0] - border; i++)
		{
			for (int k = 0; k < N; k++)
			{
				//int val = (int)floor((img(i,j)[k]*(1-2*color_padding) + color_padding)*255 + .5);
				int val = RGB32F_to_sRGB8(img(i,j)[k]);
				if (val < 0)
					val = 0;
				else if (val > 255)
					val = 255;
				fprintf(f, "%c", val);
			}
		}
	}

	fclose(f);
}
