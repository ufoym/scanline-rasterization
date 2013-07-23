#pragma once

//*****************************************************************************//
// Written by Josiah Manson, September 2012
//*****************************************************************************//

#include <omp.h>
#include <windows.h>
#include "rasterizer_shared.h"

//=============================================================================//
//================== Data Shared Between Rasterizations =======================//
//=============================================================================//

template <int W, int COLORCOMP, int COLORORDER>
struct PRasterizerInstance
{
	PRasterizerInstance()
	{
		is_created = false;
	}

	PRasterizerInstance(Array2D<vect<COLORCOMP,double> > &img, bool curve22, bool curve32, bool curve42, bool curve33)
	{
		create(img, curve22, curve32, curve42, curve33);
	}
	
	~PRasterizerInstance()
	{
		destroy();
	}

	void create(Array2D<vect<COLORCOMP,double> > &img_in, bool curve22, bool curve32, bool curve42, bool curve33)
	{
		is_created = true;
		img = &img_in;

		//omp_set_num_threads(1);
		thread_num = omp_get_max_threads();
		
		// allocate raster data for image and curve types
		raster_data.resize(thread_num);
		last_scanline.resize(thread_num);

		if (curve22)
			init_data.push_back(new RasterData<2,2,W,COLORCOMP,COLORORDER>(img->s));
		if (curve32) 
			init_data.push_back(new RasterData<3,2,W,COLORCOMP,COLORORDER>(img->s));
		if (curve42) 
			init_data.push_back(new RasterData<4,2,W,COLORCOMP,COLORORDER>(img->s));
		if (curve33) 
			init_data.push_back(new RasterData<3,3,W,COLORCOMP,COLORORDER>(img->s));
		
		for (int thread_id = 0; thread_id < thread_num; thread_id++)
		{
			raster_data[thread_id].resize(init_data.s);
			
			for (int iras = 0; iras < init_data.s; iras++)
			{
				if (init_data[iras]->curve_type == 22)
					raster_data[thread_id][iras] = new RasterData<2,2,W,COLORCOMP,COLORORDER>(img->s);
				else if (init_data[iras]->curve_type == 32)
					raster_data[thread_id][iras] = new RasterData<3,2,W,COLORCOMP,COLORORDER>(img->s);
				else if (init_data[iras]->curve_type == 42)
					raster_data[thread_id][iras] = new RasterData<4,2,W,COLORCOMP,COLORORDER>(img->s);
				else if (init_data[iras]->curve_type == 33)
					raster_data[thread_id][iras] = new RasterData<3,3,W,COLORCOMP,COLORORDER>(img->s);
			}
		}
		
#pragma omp parallel
		{
			SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_TIME_CRITICAL);
		} 
	}

	void destroy()
	{
		if (!is_created)
			return;
		is_created = false;

		// free raster data
		for (int thread_id = 0; thread_id < raster_data.s; thread_id++)
		{
			for (int iras = 0; iras < raster_data[thread_id].s; iras++)
			{
				delete raster_data[thread_id][iras];
			}
			raster_data[thread_id].s = 0;
		}
		raster_data.s = 0;

		for (int iras = 0; iras < init_data.s; iras++)
		{
			delete init_data[iras];
		}
		init_data.s = 0;
	}

	int is_created;
	Array2D<vect<COLORCOMP,double> > *img;
	int thread_num;
	Array<IRasterData*> init_data;
	Array<Array<IRasterData*> > raster_data;
	Array<int> last_scanline;
};

//=============================================================================//
//====================== Main Rasterization Function ==========================//
//=============================================================================//

template <int W, int COLORCOMP, int COLORORDER>
void rasterize_parallel(PRasterizerInstance<W, COLORCOMP, COLORORDER> &inst,
						ColorFunction<COLORCOMP, COLORORDER> &color_in, 
						Array<Curve<2,2> > *curve22, void (**filter22)(double *pixel, double *accum, Curve<2,2> *curve), 
						Array<Curve<3,2> > *curve32, void (**filter32)(double *pixel, double *accum, Curve<3,2> *curve), 
						Array<Curve<4,2> > *curve42, void (**filter42)(double *pixel, double *accum, Curve<4,2> *curve), 
						Array<Curve<3,3> > *curve33, void (**filter33)(double *pixel, double *accum, Curve<3,3> *curve))
{
	//-- serial section to set up data for parallel execution

	// allocate raster data for image and curve types
	// this builds the active edge table
	for (int iras = 0; iras < inst.init_data.s; iras++)
	{
		if (inst.init_data[iras]->curve_type == 22 && curve22 && curve22->s)
			((RasterData<2,2,W,COLORCOMP,COLORORDER>*)inst.init_data[iras])->create_scanline_table(curve22);
		else if (inst.init_data[iras]->curve_type == 32 && curve32 && curve32->s) 
			((RasterData<3,2,W,COLORCOMP,COLORORDER>*)inst.init_data[iras])->create_scanline_table(curve32);
		else if (inst.init_data[iras]->curve_type == 42 && curve42 && curve42->s) 
			((RasterData<4,2,W,COLORCOMP,COLORORDER>*)inst.init_data[iras])->create_scanline_table(curve42);
		else if (inst.init_data[iras]->curve_type == 33 && curve33 && curve33->s) 
			((RasterData<3,3,W,COLORCOMP,COLORORDER>*)inst.init_data[iras])->create_scanline_table(curve33);
	}

	// find scanline range occupied by curves
	int scan_start = 2000000000;
	int scan_end = -2000000000;
	
	for (int iras = 0; iras < inst.init_data.s; iras++)
	{
		if (inst.init_data[iras]->scan_start < scan_start)
			scan_start = inst.init_data[iras]->scan_start;
		
		if (inst.init_data[iras]->scan_end > scan_end)
			scan_end = inst.init_data[iras]->scan_end;
	}

	//-- parallel section
	// threads render the image in multiple passes so that they don't interfere with eachother.
	// box filter is a special case, and can be done in one pass.
	
	// determine how work should be broken up
	int bar_width = W - 1; // minimal width for 2 passes is W - 1. minimal width for W passes is 1
	int num_passes = 2;

	if (W == 1)
	{
		// box filter doesn't overlap, so only one pass is needed 
		bar_width = max(bar_width, 1);
		num_passes = 1;
	}

	int stride = bar_width * num_passes;

	// run the rasterizer in parallel
	for (int pass = 0; pass < num_passes; pass++)
	{
		int bar_offset = pass * bar_width;

		// copy initialized raster data for each thread
#pragma omp parallel
		{
			int thread_id = omp_get_thread_num();
			inst.last_scanline[thread_id] = scan_start;

			for (int iras = 0; iras < inst.init_data.s; iras++)
			{
				if (inst.init_data[iras]->curve_type == 22)
					((RasterData<2,2,W,COLORCOMP,COLORORDER>*)inst.raster_data[thread_id][iras])->copy((RasterData<2,2,W,COLORCOMP,COLORORDER>*)inst.init_data[iras]);
				else if (inst.init_data[iras]->curve_type == 32)
					((RasterData<3,2,W,COLORCOMP,COLORORDER>*)inst.raster_data[thread_id][iras])->copy((RasterData<3,2,W,COLORCOMP,COLORORDER>*)inst.init_data[iras]);
				else if (inst.init_data[iras]->curve_type == 42)
					((RasterData<4,2,W,COLORCOMP,COLORORDER>*)inst.raster_data[thread_id][iras])->copy((RasterData<4,2,W,COLORCOMP,COLORORDER>*)inst.init_data[iras]);
				else if (inst.init_data[iras]->curve_type == 33)
					((RasterData<3,3,W,COLORCOMP,COLORORDER>*)inst.raster_data[thread_id][iras])->copy((RasterData<3,3,W,COLORCOMP,COLORORDER>*)inst.init_data[iras]);
			}
		}

#pragma omp parallel for schedule (dynamic)
		for (int bar = scan_start + bar_offset; bar <= scan_end; bar += stride)
		{
			int thread_id = omp_get_thread_num();

			int range_start = bar;
			int range_end = min(scan_end, bar + bar_width-1);

			ColorFunction<COLORCOMP, COLORORDER> color = color_in; // thread local color, because the color is modified

			// skip to beginning
			if (range_start > inst.last_scanline[thread_id])
			{
				for (int iscan = inst.last_scanline[thread_id]; iscan < range_start; iscan++)
				{
					for (int iras = 0; iras < inst.raster_data[thread_id].s; iras++)
					{
						inst.raster_data[thread_id][iras]->skip_to_next_scanline(iscan);
					}
				}

				// set end point of $y$ curve segment and its end parameter value
				for (int iras = 0; iras < inst.raster_data[thread_id].s; iras++)
				{
					inst.raster_data[thread_id][iras]->init_end_points(range_start);
				}
			}

			// iterate over the scanlines
			for (int iscan = range_start; iscan <= range_end; iscan++)
			{
				int pixel_end = -2000000000;

				// cut curves to the next scanline / pixel buckets and find common extent
				for (int iras = 0; iras < inst.raster_data[thread_id].s; iras++)
				{
					inst.raster_data[thread_id][iras]->cut_next_scanline(iscan);

					if (inst.raster_data[thread_id][iras]->pixel_end > pixel_end)
						pixel_end = inst.raster_data[thread_id][iras]->pixel_end;
				}

				// evaluate pixel values. all curve types need raster up to the same ending point
				for (int iras = 0; iras < inst.raster_data[thread_id].s; iras++)
				{
					inst.raster_data[thread_id][iras]->pixel_end = pixel_end;

					if (inst.raster_data[thread_id][iras]->curve_type == 22)
						((RasterData<2,2,W,COLORCOMP,COLORORDER>*)inst.raster_data[thread_id][iras])->eval_pixels(*inst.img, color, filter22, iscan);
					else if (inst.raster_data[thread_id][iras]->curve_type == 32)
						((RasterData<3,2,W,COLORCOMP,COLORORDER>*)inst.raster_data[thread_id][iras])->eval_pixels(*inst.img, color, filter32, iscan);
					else if (inst.raster_data[thread_id][iras]->curve_type == 42)
						((RasterData<4,2,W,COLORCOMP,COLORORDER>*)inst.raster_data[thread_id][iras])->eval_pixels(*inst.img, color, filter42, iscan);
					else if (inst.raster_data[thread_id][iras]->curve_type == 33)
						((RasterData<3,3,W,COLORCOMP,COLORORDER>*)inst.raster_data[thread_id][iras])->eval_pixels(*inst.img, color, filter33, iscan);
				}
			}

			inst.last_scanline[thread_id] = range_end + 1;
		}
	}
}
