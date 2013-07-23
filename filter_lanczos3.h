#pragma once

//*****************************************************************************//
// Written by Josiah Manson, September 2012
//*****************************************************************************//

#include "rasterizer_common.h"

template <int COLORORDER>
struct FilterLanczos3
{
	static const int W = 6;

	void (*filter22[W*W])(double *pixel, double *accum, Curve<2,2> *curve);
	void (*filter32[W*W])(double *pixel, double *accum, Curve<3,2> *curve);
	void (*filter42[W*W])(double *pixel, double *accum, Curve<4,2> *curve);
	void (*filter33[W*W])(double *pixel, double *accum, Curve<3,3> *curve);

	FilterLanczos3();
};
