#pragma once

//*****************************************************************************//
// Written by Josiah Manson, September 2012
//*****************************************************************************//

#include "rasterizer_shared.h"


//=============================================================================//
//================== Data Shared Between Rasterizations =======================//
//=============================================================================//

template <int W, int COLORCOMP, int COLORORDER>
struct RasterizerInstance
{
	RasterizerInstance()
	{
		is_created = false;
	}

	RasterizerInstance(Array2D<vect<COLORCOMP,double> > &img, bool curve22, bool curve32, bool curve42, bool curve33)
	{
		create(img, curve22, curve32, curve42, curve33);
	}
	
	~RasterizerInstance()
	{
		destroy();
	}

	void create(Array2D<vect<COLORCOMP,double> > &img_in, bool curve22, bool curve32, bool curve42, bool curve33)
	{
		is_created = true;
		img = &img_in;

		// allocate raster data for image and curve types
		if (curve22)
			raster_data.push_back(new RasterData<2,2,W,COLORCOMP,COLORORDER>(img->s));
		if (curve32) 
			raster_data.push_back(new RasterData<3,2,W,COLORCOMP,COLORORDER>(img->s));
		if (curve42) 
			raster_data.push_back(new RasterData<4,2,W,COLORCOMP,COLORORDER>(img->s));
		if (curve33) 
			raster_data.push_back(new RasterData<3,3,W,COLORCOMP,COLORORDER>(img->s));
	}

	void destroy()
	{
		if (!is_created)
			return;
		is_created = false;

		// free raster data
		for (int iras = 0; iras < raster_data.s; iras++)
		{
			delete raster_data[iras];
		}
		raster_data.s = 0;
	}

	int is_created;
	Array2D<vect<COLORCOMP,double> > *img;
	Array<IRasterData*> raster_data;
};

//=============================================================================//
//====================== Main Rasterization Function ==========================//
//=============================================================================//

// I could use nicer syntax for the filter function pointers if VS11 supported C++11 fully. The structure workaround doesn't work, because typenames need to be concrete with template aliases.
//template <int PTS, int DIM> using FilterPiece = void (*)(double *pixel, double *accum, Curve<PTS,DIM> *curve); // C++11 feature that isn't supported by the compiler

template <int W, int COLORCOMP, int COLORORDER>
void rasterize(RasterizerInstance<W, COLORCOMP, COLORORDER> &inst,
			   ColorFunction<COLORCOMP, COLORORDER> &color, 
			   Array<Curve<2,2> > *curve22, void (**filter22)(double *pixel, double *accum, Curve<2,2> *curve), 
			   Array<Curve<3,2> > *curve32, void (**filter32)(double *pixel, double *accum, Curve<3,2> *curve), 
			   Array<Curve<4,2> > *curve42, void (**filter42)(double *pixel, double *accum, Curve<4,2> *curve), 
			   Array<Curve<3,3> > *curve33, void (**filter33)(double *pixel, double *accum, Curve<3,3> *curve))
{
	// allocate raster data for image and curve types
	// this builds the active edge table
	for (int iras = 0; iras < inst.raster_data.s; iras++)
	{
		if (inst.raster_data[iras]->curve_type == 22 && curve22 && curve22->s)
			((RasterData<2,2,W,COLORCOMP,COLORORDER>*)inst.raster_data[iras])->create_scanline_table(curve22);
		else if (inst.raster_data[iras]->curve_type == 32 && curve32 && curve32->s) 
			((RasterData<3,2,W,COLORCOMP,COLORORDER>*)inst.raster_data[iras])->create_scanline_table(curve32);
		else if (inst.raster_data[iras]->curve_type == 42 && curve42 && curve42->s) 
			((RasterData<4,2,W,COLORCOMP,COLORORDER>*)inst.raster_data[iras])->create_scanline_table(curve42);
		else if (inst.raster_data[iras]->curve_type == 33 && curve33 && curve33->s) 
			((RasterData<3,3,W,COLORCOMP,COLORORDER>*)inst.raster_data[iras])->create_scanline_table(curve33);
	}

	// find scanline range occupied by curves
	int scan_start = 2000000000;
	int scan_end = -2000000000;
	
	for (int iras = 0; iras < inst.raster_data.s; iras++)
	{
		if (inst.raster_data[iras]->scan_start < scan_start)
			scan_start = inst.raster_data[iras]->scan_start;
		
		if (inst.raster_data[iras]->scan_end > scan_end)
			scan_end = inst.raster_data[iras]->scan_end;
	}

	// iterate over the scanlines
	for (int iscan = scan_start; iscan <= scan_end; iscan++)
	{
		int pixel_end = -2000000000;

		// cut curves to the next scanline / pixel buckets and find common extent
		for (int iras = 0; iras < inst.raster_data.s; iras++)
		{
			inst.raster_data[iras]->cut_next_scanline(iscan);

			if (inst.raster_data[iras]->pixel_end > pixel_end)
				pixel_end = inst.raster_data[iras]->pixel_end;
		}

		// evaluate pixel values. all curve types need raster up to the same ending point
		for (int iras = 0; iras < inst.raster_data.s; iras++)
		{
			inst.raster_data[iras]->pixel_end = pixel_end;
			
			if (inst.raster_data[iras]->curve_type == 22)
				((RasterData<2,2,W,COLORCOMP,COLORORDER>*)inst.raster_data[iras])->eval_pixels(*inst.img, color, filter22, iscan);
			else if (inst.raster_data[iras]->curve_type == 32)
				((RasterData<3,2,W,COLORCOMP,COLORORDER>*)inst.raster_data[iras])->eval_pixels(*inst.img, color, filter32, iscan);
			else if (inst.raster_data[iras]->curve_type == 42)
				((RasterData<4,2,W,COLORCOMP,COLORORDER>*)inst.raster_data[iras])->eval_pixels(*inst.img, color, filter42, iscan);
			else if (inst.raster_data[iras]->curve_type == 33)
				((RasterData<3,3,W,COLORCOMP,COLORORDER>*)inst.raster_data[iras])->eval_pixels(*inst.img, color, filter33, iscan);
		}
	}
}

