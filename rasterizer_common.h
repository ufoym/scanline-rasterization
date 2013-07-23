#pragma once

//*****************************************************************************//
// Written by Josiah Manson, September 2012
//*****************************************************************************//

#include "vect.h"

//=============================================================================//
//=============================== Base Curve ==================================//
//=============================================================================//

// A curve with PTS control points, and DIM dimensions. 2 dimensions for bezier, and 3 for rational bezier.
template <int PTS, int DIM>
struct Curve
{
	vect<DIM, double> p[PTS];
};

//=============================================================================//
//============================ Filter Indexing ================================//
//=============================================================================//

template <int W>
inline int filter_elem(int x, int y)
{
	return y*W + x;
}

template <int W>
inline int filter_elem_flip(int x, int y)
{
	return filter_elem<W>(x, W-1-y);
}
