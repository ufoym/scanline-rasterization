#pragma once

//*****************************************************************************//
// Written by Josiah Manson, September 2012
//*****************************************************************************//


#include <stdio.h>
#include <math.h>
#include <algorithm>

#include "rasterizer_common.h"
#include "boundedarray.h"
#include "array.h"
#include "array2d.h"

using namespace std;

//=============================================================================//
//================================= Utility ===================================//
//=============================================================================//

template <int PTS, int DIM>
void print_curve(Curve<PTS, DIM> &c)
{
	printf("{");
	for (int ip = 0; ip < PTS; ip++)
	{
		if (ip > 0)
			printf(",");

		if (DIM == 2)
			printf("{%f,%f}", c.p[ip][0], c.p[ip][1]);
		else if (DIM == 3)
			printf("{%f,%f,%f}", c.p[ip][0], c.p[ip][1], c.p[ip][2]);
	}
	printf("},");
}

//=============================================================================//
//================================= Colors ====================================//
//=============================================================================//

// coefficients for the color function. a constant function will have one variable, and a linear will have three
template <int COLORORDER>
struct ColorComponent
{
	static const int NUM = COLORORDER*(COLORORDER+1)/2;
	double c[NUM];
};

// color function that can have N components, monochrom, rgb, rgba
template <int COLORCOMP, int COLORORDER>
struct ColorFunction
{
	ColorComponent<COLORORDER> comp[COLORCOMP]; // coefficients for each component
};


//=============================================================================//
//================================== Solve ====================================//
//=============================================================================//

// template for bezier solving
template <int dir, int PTS, int DIM>
double solve_bez(Curve<PTS, DIM> *bez, double h)
{
	return 0;
}

//-- linear specialization
template <int dir>
inline double solve_bez(Curve<2, 2> *bez, double h)
{
	double &p0 = bez->p[0][dir];
	double &p1 = bez->p[1][dir];

	return (h - p0) / (p1 - p0);
}

//-- quadratic specialization

// Solving a quadratic bezier reduced to its two fundamental degrees of freedom, normalized for y0 = 0, y2 = 1.
// The issue with y1~=.5 is that both the numerator and the denominator approach zero, and precision is lost.
// In this case, I use a Taylor expansion of the solution, which is quite accurate.
inline double solve_bez_quad_normalized(const double y1, const double h)
{
	const double e = y1 - .5;
	if (-.0125 < e && e < .0125)
	{
#if 1
		// taylor expansions of increasing order 
		// (I have put tolerances required for accuracy to 6 significant figures)
		//double t = h; // |e| < .000002
		//double t = h + 2*e*(-1 + h)*h; // |e| < .002
		double t = h*(1 + 2*e*(-1 + h)*(1 + e*(-2 + 4*h))); // |e| < .0125
		//double t = h*(1 + 2*e*(-1 + h)*(1 + 2*e*(-1 + 2*h + 2*e*(1 + 5*(-1 + h)*h)))); // |e| < .04
#else
		// newton iterations are a bit slow
		double t = h;
		for (int i = 0; i < 2; i++)
			t -= (h + t*(-t + 2*(-1 + t)*y1))/(2*(t*y2 - y1));
#endif

		return t;
	}
	else
	{
		// standard solution
		return (y1 - sqrt(h - 2*h*y1 + y1*y1))/(2*e);
	}
}

template <int dir>
inline double solve_bez(Curve<3, 2> *bez, double h)
{
	double &p0 = bez->p[0][dir];
	double &p1 = bez->p[1][dir];
	double &p2 = bez->p[2][dir];

	double Y1 = (p1 - p0) / (p2 - p0);
	double H = (h - p0) / (p2 - p0);

	return solve_bez_quad_normalized(Y1, H);
}

//-- cubic

// Reduce degrees of freedom by normalizing y0=0, y3=1. The strategy I use in this this function is to 
// initialize from a quadratic approximation and to then perform Newton iterations until convergence.
inline double solve_bez_cube_normalized(const double y1, const double y2, const double h)
{
	// Match the ends. This works best when there is an inflection.
	//double y1appx = y1*(1-h) + y2*h;
	
	// Best fit quadratic. This works best when there is no inflection.
	double y1appx = (-1 + 3*y1 + 3*y2)*.25; 
	
	double t = solve_bez_quad_normalized(y1appx, h);

#if 1
	// use newton iterations. this should converge very quickly with a good guess
	double dt = 1;
	while (dt > 1e-6 || dt < -1e-6)
	{
		dt = (-h + t*(3*y1 + t*(t - 6*y1 + 3*t*y1 + 3*y2 - 3*t*y2)))/(3.*(y1 + t*(t - 4*y1 + 3*t*y1 + 2*y2 - 3*t*y2)));
		t -= dt;
	}
#else
	// Use bounded newton iterations. This should converge quickly.
	// Steps are bounded, with the bound decreasing each step
	// so that in the worst case there is linear convergence.
	// This may be useful when cutting at inflections is
	// disabled, because inflections can cause very bad steps that
	// this prevents.
	double b = .25;
	double dt = 1;
	while (dt > 1e-6 || dt < -1e-6)
	{
		dt = (-h + t*(3*y1 + t*(t - 6*y1 + 3*t*y1 + 3*y2 - 3*t*y2)))/(3.*(y1 + t*(t - 4*y1 + 3*t*y1 + 2*y2 - 3*t*y2)));
		
		double sign = dt < 0 ? -1 : 1;
		if (sign*dt > b)
			t -= sign*b;
		else
			t -= dt;

		b *= .75;
	}
#endif

	return t;
}

template <int dir>
inline double solve_bez(Curve<4, 2> *bez, double h)
{
	double &p0 = bez->p[0][dir];
	double &p1 = bez->p[1][dir];
	double &p2 = bez->p[2][dir];
	double &p3 = bez->p[3][dir];

	double Y1 = (p1 - p0) / (p3 - p0);
	double Y2 = (p2 - p0) / (p3 - p0);
	double H = (h - p0) / (p3 - p0);

	return solve_bez_cube_normalized(Y1, Y2, H);
}

//-- rational quadratic
template <int dir>
inline double solve_bez(Curve<3, 3> *bez, double h)
{
	double &p0 = bez->p[0][dir];
	double &p1 = bez->p[1][dir];
	double &p2 = bez->p[2][dir];
	
	double &w0 = bez->p[0][2];
	double &w1 = bez->p[1][2];
	double &w2 = bez->p[2][2];
		
	double c0 = (p0 - h) * w0;
	double c1 = (p1 - h) * w1;
	double c2 = (p2 - h) * w2;

	double Y1 = (c1 - c0) / (c2 - c0);
	double H = (-c0) / (c2 - c0);

	return solve_bez_quad_normalized(Y1, H);
}

//=============================================================================//
//================================ Evaluate ===================================//
//=============================================================================//

// template for cutting a curve to the parameter range $[s,t]$. the first point that depends only on $s$ is already known, so only solve for the remainder of control points that depend on $t$.
template <int PTS, int DIM>
void eval_bez(vect<DIM, double> &p, Curve<PTS, DIM> *original, double &t)
{
}

// linear specialization
template <>
inline void eval_bez(vect<2, double> &p, Curve<2, 2> *original, double &t)
{
	Curve<2, 2> &b = *original;
	p = b.p[0] * (1-t) + b.p[1] * t;
}

// quadratic specialization
template <>
inline void eval_bez(vect<2, double> &p, Curve<3, 2> *original, double &t)
{
	Curve<3, 2> &b = *original;

	Curve<2, 2> bt;
	bt.p[0] = b.p[0] * (1-t) + b.p[1] * t;
	bt.p[1] = b.p[1] * (1-t) + b.p[2] * t;
	
	p = bt.p[0] * (1-t) + bt.p[1] * t;
}

// cubic specialization
template <>
inline void eval_bez(vect<2, double> &p, Curve<4, 2> *original, double &t)
{
	Curve<4, 2> &b = *original;

	Curve<3, 2> bt;
	bt.p[0] = b.p[0] * (1-t) + b.p[1] * t;
	bt.p[1] = b.p[1] * (1-t) + b.p[2] * t;
	bt.p[2] = b.p[2] * (1-t) + b.p[3] * t;
	
	Curve<2, 2> btt;
	btt.p[0] = bt.p[0] * (1-t) + bt.p[1] * t;
	btt.p[1] = bt.p[1] * (1-t) + bt.p[2] * t;
	
	p = btt.p[0] * (1-t) + btt.p[1] * t;
}

// rational quadratic specialization
template <>
inline void eval_bez(vect<3, double> &p, Curve<3, 3> *original, double &t)
{
	// weights are not premultiplied into the numerator
	Curve<3, 3> b;
	for (int i = 0; i < 3; i++)
	{
		b.p[i][0] = original->p[i][0] * original->p[i][2];
		b.p[i][1] = original->p[i][1] * original->p[i][2];
		b.p[i][2] = original->p[i][2];
	}

	// standard bezier cutting, now that weights are corrected
	Curve<2, 3> bt;
	bt.p[0] = b.p[0] * (1-t) + b.p[1] * t;
	bt.p[1] = b.p[1] * (1-t) + b.p[2] * t;
	
	p = bt.p[0] * (1-t) + bt.p[1] * t;
	
	// we need to project the numerator so that it has geometric meaning in 2D
	p[0] /= p[2];
	p[1] /= p[2];
}


//=============================================================================//
//================================== Cut ======================================//
//=============================================================================//

// template for cutting a curve to the parameter range $[s,t]$. the first point that depends only on $s$ is already known, so only solve for the remainder of control points that depend on $t$.
template <int PTS, int DIM>
void cut_bez_span(Curve<PTS, DIM> &curve, Curve<PTS, DIM> *original, double &s, double &t)
{
}

// linear specialization
template <>
inline void cut_bez_span(Curve<2, 2> &curve, Curve<2, 2> *original, double &s, double &t)
{
	Curve<2, 2> &b = *original;
	curve.p[1] = b.p[0] * (1-t) + b.p[1] * t;
}

// quadratic specialization
template <>
inline void cut_bez_span(Curve<3, 2> &curve, Curve<3, 2> *original, double &s, double &t)
{
	Curve<3, 2> &b = *original;

	Curve<2, 2> bt;
	bt.p[0] = b.p[0] * (1-t) + b.p[1] * t;
	bt.p[1] = b.p[1] * (1-t) + b.p[2] * t;
	
	curve.p[1] = bt.p[0] * (1-s) + bt.p[1] * s;
	curve.p[2] = bt.p[0] * (1-t) + bt.p[1] * t;
}

// cubic specialization
template <>
inline void cut_bez_span(Curve<4, 2> &curve, Curve<4, 2> *original, double &s, double &t)
{
	Curve<4, 2> &b = *original;

	Curve<3, 2> bt;
	bt.p[0] = b.p[0] * (1-t) + b.p[1] * t;
	bt.p[1] = b.p[1] * (1-t) + b.p[2] * t;
	bt.p[2] = b.p[2] * (1-t) + b.p[3] * t;
	
	Curve<2, 2> btt;
	btt.p[0] = bt.p[0] * (1-t) + bt.p[1] * t;
	btt.p[1] = bt.p[1] * (1-t) + bt.p[2] * t;
	
	Curve<2, 2> bst;
	bst.p[0] = bt.p[0] * (1-s) + bt.p[1] * s;
	bst.p[1] = bt.p[1] * (1-s) + bt.p[2] * s;

	curve.p[1] = bst.p[0] * (1-s) + bst.p[1] * s;
	curve.p[2] = btt.p[0] * (1-s) + btt.p[1] * s;
	curve.p[3] = btt.p[0] * (1-t) + btt.p[1] * t;
}

// rational quadratic specialization
template <>
inline void cut_bez_span(Curve<3, 3> &curve, Curve<3, 3> *original, double &s, double &t)
{
	// weights are not premultiplied into the numerator
	Curve<3, 3> b;
	for (int i = 0; i < 3; i++)
	{
		b.p[i][0] = original->p[i][0] * original->p[i][2];
		b.p[i][1] = original->p[i][1] * original->p[i][2];
		b.p[i][2] = original->p[i][2];
	}

	// standard bezier cutting, now that weights are corrected
	Curve<2, 3> bt;
	bt.p[0] = b.p[0] * (1-t) + b.p[1] * t;
	bt.p[1] = b.p[1] * (1-t) + b.p[2] * t;
	
	curve.p[1] = bt.p[0] * (1-s) + bt.p[1] * s;
	curve.p[2] = bt.p[0] * (1-t) + bt.p[1] * t;
	
	// we need to project the numerator so that it has geometric meaning in 2D
	for (int i = 1; i < 3; i++)
	{
		curve.p[i][0] /= curve.p[i][2];
		curve.p[i][1] /= curve.p[i][2];
	}
}

// cut into N+1 pieces for N t_values
template <int PTS, int DIM>
void cut_bez_pieces(Array<Curve<PTS,DIM> > *dest, Curve<PTS,DIM> &orig, Array<double> &t_vals)
{
	// if no cuts, just copy and quit
	if (t_vals.s == 0)
	{
		dest->push_back(orig);
		return;
	}

	// sort the parameter values, and add the last point t=1. the first point t=0 is implied
	sort(t_vals.data, t_vals.data + t_vals.s);
	t_vals.push_back(1);

	// cut into pieces
	double s = 0;
	Curve<PTS, DIM> seg;
	seg.p[0] = orig.p[0];

	for (int it = 0; it < t_vals.s; it++)
	{
		if (t_vals[it] - s < 1e-6)
			continue;
		cut_bez_span(seg, &orig, s, t_vals[it]);
		dest->push_back(seg);
		seg.p[0] = seg.p[PTS-1];
		s = t_vals[it];
	}
}

//=============================================================================//
//=========================== Monotonic Cutting ===============================//
//=============================================================================//

// this nonspecialized function is a placeholder (it also works for lines, because lines are monotonic)
template <int PTS, int DIM>
void find_monotonic_params(Array<double> &t_vals, Curve<PTS,DIM> &orig, int dir)
{
}

template <>
inline void find_monotonic_params(Array<double> &t_vals, Curve<3,2> &orig, int dir)
{
	const double y0 = orig.p[0][dir];
	const double y1 = orig.p[1][dir];
	const double y2 = orig.p[2][dir];
	
	double t = (y0 - y1)/(y0 - 2*y1 + y2);

	if (t > 0+1e-6 && t < 1-1e-6)
		t_vals.push_back(t);
}

template <>
inline void find_monotonic_params(Array<double> &t_vals, Curve<4,2> &orig, int dir)
{
	const double y0 = orig.p[0][dir];
	const double y1 = orig.p[1][dir];
	const double y2 = orig.p[2][dir];
	const double y3 = orig.p[3][dir];

	double a = y0 - 2*y1 + y2;
	double b = y0 - 3*y1 + 3*y2 - y3;
	double c = y1*y1 - y0*y2 - y1*y2 + y2*y2 + y0*y3 - y1*y3;

	if (b != 0)
	{
		// inflection
		double t = a/b;
		if (t > 0+1e-6 && t < 1-1e-6)
			t_vals.push_back(t);

		// direction changes
		if (c > 0)
		{
			double d = sqrt(c);

			t = (a-d)/b;
			if (t > 0+1e-6 && t < 1-1e-6)
				t_vals.push_back(t);

			t = (a+d)/b;
			if (t > 0+1e-6 && t < 1-1e-6)
				t_vals.push_back(t);
		}
	}
}

template <>
inline void find_monotonic_params(Array<double> &t_vals, Curve<3,3> &orig, int dir)
{
	const double w0 = orig.p[0][2];
	const double w1 = orig.p[1][2];
	const double w2 = orig.p[2][2];
	const double p0 = w0*orig.p[0][dir];
	const double p1 = w1*orig.p[1][dir];
	const double p2 = w2*orig.p[2][dir];
	
	double a = -2*p1*w0 + p2*w0 + 2*p0*w1 - p0*w2;
	double d = sqrt(p2*p2*w0*w0 + w2*(4*p1*p1*w0 - 4*p0*p1*w1 + p0*p0*w2) - 2*p2*(2*p1*w0*w1 - 2*p0*w1*w1 + p0*w0*w2));
	double b = 2*(p2*(w0 - w1) + p0*(w1 - w2) + p1*(-w0 + w2));
	
	double t;

	t = (a - d) / b;
	if (t > 0+1e-6 && t < 1-1e-6)
		t_vals.push_back(t);

	t = (a + d) / b;
	if (t > 0+1e-6 && t < 1-1e-6)
		t_vals.push_back(t);
}


// add cut pieces of the original curve to the destination list
// the curves will be cut to be monotonic segments in both $x$ and $y$, as well as inflection free (for nicer solutions to cubics)
template <int PTS, int DIM>
void cut_monotonic(Array<Curve<PTS,DIM> > *dest, Curve<PTS,DIM> &orig)
{
	static Array<double> t_vals;
	t_vals.s = 0;

	//-- get parameters to cut at
	find_monotonic_params(t_vals, orig, 0);
	find_monotonic_params(t_vals, orig, 1);

	//-- cut the curve into pieces
	cut_bez_pieces(dest, orig, t_vals);
}


//=============================================================================//
//=========================== Specialized Curves ==============================//
//=============================================================================//

// a linked list of curves that are cut to pixels
template <int PTS, int DIM>
struct BucketCurve
{
	BucketCurve<PTS, DIM> *next; // next element in the linked list
	Curve<PTS, DIM> curve; // the curve contained in the pixel. this is stored so that whatever flipping was done earlier is undone
};

// a segment of a curve that is incrementally cut
template <int PTS, int DIM>
struct IncrementalCurveSegment
{
	Curve<PTS, DIM> curve;
	double t_seg_start, t_seg_end; // parameter value of the start / end of the old segment
	double t_curve_end; // parameter value of the end of the curve (this is always 1 for cutting in y, but is the end value of the y segment when cutting in x)
	int end_cell; // cell index to stop cutting at
	int flipped; // 0 / 1 if curve is not / is flipped compared to the original curve
};

// a curve that is incrementally cut. this structure is used for cutting in both the X and Y directions, and stores a segment for each
template <int PTS, int DIM>
struct IncrementalCurve
{
	IncrementalCurve<PTS, DIM> *next; // next element in the linked list
	Curve<PTS, DIM> *original; // the original curve that is being cut against
	IncrementalCurveSegment<PTS, DIM> seg[2]; // segments that are cut in the x=0 and y=1 directions

	template <int dir>
	void increment(int curr_cell)
	{
		//-- find the starting parameter $s$ and ending parameter $t$ of the next segment
		double s = seg[dir].t_seg_end;
		double t;
		if (seg[dir].end_cell <= curr_cell)
			t = seg[dir].t_curve_end;
		else
			t = solve_bez<dir>(original, curr_cell + 1);

		//-- find the new control points. 
		// we already know that the starting point of the new segment is the ending point of the old segment
		seg[dir].curve.p[0] = seg[dir].curve.p[PTS-1];

		// calculate the control points that depend on $t$
		cut_bez_span(seg[dir].curve, original, s, t);

		//-- store the parameter values of the segment
		if (dir)
			seg[dir].t_seg_start = s;
		seg[dir].t_seg_end = t;
	}

	void initial_cut(int curr_cell) // TODO: seems to not work for higher order curves
	{
		//-- find the ending parameter $t$ of the segment, which will be incrementally updated.
		double t;
		
		if (seg[1].end_cell <= curr_cell - 1)
			t = seg[1].t_curve_end;
		else
			t = solve_bez<1>(original, curr_cell);

		//-- find the new end control point. 
		eval_bez(seg[1].curve.p[PTS-1], original, t);

		//-- store the parameter values of the segment
		seg[1].t_seg_end = t;
	}
};

// stores curves in each pixel
template <int PTS, int DIM>
struct PixelTable
{
	Array<BucketCurve<PTS, DIM>*> buckets; // stores curves in each pixel (pointers to heap)
	BoundedArray<BucketCurve<PTS, DIM>, 100000> heap; // stores the curves in the buckets THIS COULD OVERFLOW
};

// stores the curves that start on each scanline
template <int PTS, int DIM>
struct ScanlineTable
{
	Array<IncrementalCurve<PTS, DIM>*> buckets; // stores the curves that start on each scanline (pointers to heap)
	Array<IncrementalCurve<PTS, DIM> > heap; // stores the curves in the buckets
};


//=============================================================================//
//======================= Curve Specific Rasterizer ===========================//
//=============================================================================//

// stores curve data that is specific to a particular curve type
struct IRasterData 
{
	virtual ~IRasterData() {}

	int curve_type; // what sort of curve it is. 2 digits, first is number of points and second is dimension.
	int pixel_start, pixel_end; // starting and ending pixels in the x direction
	int scan_start, scan_end; // starting and ending scanlines in the y direction
	
	virtual void init_end_points(int curr_y) = 0;
	virtual void skip_to_next_scanline(int curr_y) = 0;
	virtual void cut_next_scanline(int curr_y) = 0;
};


template <int PTS, int DIM, int W, int COLORCOMP, int COLORORDER>
struct RasterData : public IRasterData
{
	PixelTable<PTS, DIM> ptab; // buckets of cut curves, with one bucket per pixel
	ScanlineTable<PTS, DIM> stab; // table of curves stored on the scanline they start at
	IncrementalCurve<PTS, DIM> *active; // list of active curves for cutting in the x and y directions
	Array<Curve<PTS,DIM> > *curves; // the list of curves that are potentially cut. this is not used for lines, and should only be freed for higher order curves

	void copy(RasterData<PTS,DIM,W,COLORCOMP,COLORORDER> *in)
	{
		scan_start = in->scan_start;
		scan_end = in->scan_end;

		curves = 0;
		active = 0;

		ptab.heap.s = 0;
		ptab.buckets.resize(in->ptab.buckets.s);
		ptab.buckets = 0;

		stab.heap.reserve(in->stab.heap.s);
		stab.heap.s = 0;
		stab.buckets.resize(in->stab.buckets.s);
		stab.buckets = 0;

		// copy stab (active edge table, because the incremental curves must be thread local)
		// this is slightly trickier because we must make new incremental curves and have the linked lists be the same
		for (int iscan = scan_start; iscan <= scan_end; iscan++)
		{
			IncrementalCurve<PTS, DIM> *iter = in->stab.buckets[iscan];
			while (iter)
			{
				IncrementalCurve<PTS, DIM> &curve = stab.heap.push_back();

				curve.original = iter->original;
				curve.seg[1] = iter->seg[1];
				curve.next = stab.buckets[iscan];
				assert(iscan >=0 && iscan < stab.buckets.s);
				stab.buckets[iscan] = &curve;

				//-- update iterator
				iter = iter->next;
			}
		}
	}

	void create_scanline_table(Array<Curve<PTS,DIM> > *curves_in)
	{
		scan_start = 2000000000;
		scan_end = -2000000000;
		active = 0;

		// cut monotonic
		if (PTS == 2 && W == 1)
		{
			curves = curves_in;
		}
		else
		{
			if (curves == 0)
				curves = new Array<Curve<PTS,DIM> >();
			curves->s = 0;
			curves->reserve(curves_in->s * 2);
			for (int ic = 0; ic < curves_in->s; ic++)
				cut_monotonic(curves, curves_in->get(ic));
		}

		// insert curves into table
		stab.heap.resize(curves->s);
		stab.buckets = 0;

		for (int icur = 0; icur < curves->s; icur++)
		{
			Curve<PTS, DIM> &orig = curves->get(icur);
			IncrementalCurve<PTS, DIM> &inc = stab.heap[icur];

			// offset curve either 0 or 1/2 pixels
			for (int ip = 0; ip < PTS; ip++)
				for (int id = 0; id < 2; id++)
					orig.p[ip][id] += (W-1)*.5 - (int)((W-1)*.5);

			// initialize data common between x and y segments
			inc.original = &orig;

			// initialize segment in y-direction so that cutting it once will make it the first segment
			int start_y;

			if (orig.p[PTS-1][1] < orig.p[0][1])
			{
				// the curve is flipped wrt y
				inc.seg[1].flipped = 1;
				start_y = (int)orig.p[PTS-1][1];
				inc.seg[1].end_cell = (int)orig.p[0][1];
				if (orig.p[0][1] == floor(orig.p[0][1]))
					inc.seg[1].end_cell--;
				
				inc.seg[1].curve.p[PTS-1] = orig.p[PTS-1];
				inc.seg[1].t_seg_end = 1;
				inc.seg[1].t_curve_end = 0;
			}
			else
			{
				// the curve is not flipped wrt y
				inc.seg[1].flipped = 0;
				start_y = (int)orig.p[0][1];
				inc.seg[1].end_cell = (int)orig.p[PTS-1][1];
				if (orig.p[PTS-1][1] == floor(orig.p[PTS-1][1]))
					inc.seg[1].end_cell--;
				
				inc.seg[1].curve.p[PTS-1] = orig.p[0];
				inc.seg[1].t_seg_end = 0;
				inc.seg[1].t_curve_end = 1;
			}
			
			if (start_y < 0)
				start_y = 0;
			if (start_y > stab.buckets.s - W)
				start_y = stab.buckets.s - W;
			if (inc.seg[1].end_cell < start_y)
				inc.seg[1].end_cell = start_y;
			if (inc.seg[1].end_cell >= stab.buckets.s)
				inc.seg[1].end_cell = stab.buckets.s - 1;
			
			if (start_y < scan_start)
				scan_start = start_y;
			if (inc.seg[1].end_cell > scan_end)
				scan_end = inc.seg[1].end_cell;

			// stick in front of head
			assert(start_y >=0 && start_y < stab.buckets.s);
			inc.next = stab.buckets[start_y];
			stab.buckets[start_y] = &inc;
		}
	}

	RasterData(RasterData<PTS,DIM,W,COLORCOMP,COLORORDER> *in)
	{
		curve_type = in->curve_type;
		curves = 0;

		copy(in);
	}

	RasterData(int img_s[2])
	{
		curve_type = PTS * 10 + DIM;
		curves = 0;

		// init some vars
		ptab.buckets.resize(img_s[0]);
		ptab.buckets = 0;

		stab.buckets.resize(img_s[1]);
		stab.buckets = 0;

		scan_start = 2000000000;
		scan_end = -2000000000;
		
		active = 0;
	}

	RasterData(Array<Curve<PTS,DIM> > *curves_in, int img_s[2])
	{
		curve_type = PTS * 10 + DIM;
		curves = 0;

		// image dependent vars
		ptab.buckets.resize(img_s[0]);
		ptab.buckets = 0;
		
		stab.buckets.resize(img_s[1]);

		// create scanline table
		create_scanline_table(curves_in);
	}

	virtual ~RasterData()
	{
		if (PTS == 2 && W == 1)
		{
		}
		else
		{
			// delete monotonic cut curves
			delete curves;
		}
	}
	
	void update_active_edge_list(int curr_y)
	{
		// remove old edges
		IncrementalCurve<PTS, DIM> *iter = active, *prev = 0;
		while (iter)
		{
			if (curr_y > iter->seg[1].end_cell)
			{
				if (prev == 0 && iter->next == 0)
				{
					active = 0;
					break;
				}

				if (prev)
					prev->next = iter->next;
				else
					active = iter->next;

				iter = iter->next;
				continue;
			}
			
			// update iterator
			prev = iter;
			iter = iter->next;
		}
		
		// add curves starting on this scanline to the active edge list
		if (prev)
			prev->next = stab.buckets[curr_y];
		else
			active = stab.buckets[curr_y];
	}

	virtual void skip_to_next_scanline(int curr_y)
	{
		update_active_edge_list(curr_y);
	}
	
	virtual void init_end_points(int curr_y)
	{
		IncrementalCurve<PTS, DIM> *iter = active;
		while (iter)
		{
			//-- init curve
			iter->initial_cut(curr_y);
			//iter->increment<1>(curr_y-1);
			
			//-- update iterator
			iter = iter->next;
		}
	}

	virtual void cut_next_scanline(int curr_y)
	{
		update_active_edge_list(curr_y);

		// increment and cut all of the curves in the active curve list
		ptab.heap.clear();
		pixel_start = 2000000000;
		pixel_end = -2000000000;

		IncrementalCurve<PTS, DIM> *iter = active;
		while (iter)
		{
			//-- increment curve
			iter->increment<1>(curr_y);
			Curve<PTS, DIM> &curve = iter->seg[1].curve;
			Curve<PTS, DIM> &orig = *iter->original;

			//-- find pixel extents
			int start_x;

			if (orig.p[PTS-1][0] < orig.p[0][0])
			{
				// the curve is flipped wrt x
				iter->seg[0].flipped = 1;
				
				if (iter->seg[1].flipped)
				{
					start_x = (int)iter->seg[1].curve.p[0][0];
					iter->seg[0].end_cell = (int)iter->seg[1].curve.p[PTS-1][0];
					if (iter->seg[1].curve.p[PTS-1][0] == floor(iter->seg[1].curve.p[PTS-1][0]))
						iter->seg[0].end_cell--;
					iter->seg[0].curve.p[PTS-1] = iter->seg[1].curve.p[0];
					iter->seg[0].t_seg_end = iter->seg[1].t_seg_start;
					iter->seg[0].t_curve_end = iter->seg[1].t_seg_end;
				}
				else
				{
					start_x = (int)iter->seg[1].curve.p[PTS-1][0];
					iter->seg[0].end_cell = (int)iter->seg[1].curve.p[0][0];
					if (iter->seg[1].curve.p[0][0] == floor(iter->seg[1].curve.p[0][0]))
						iter->seg[0].end_cell--;
					iter->seg[0].curve.p[PTS-1] = iter->seg[1].curve.p[PTS-1];
					iter->seg[0].t_seg_end = iter->seg[1].t_seg_end;
					iter->seg[0].t_curve_end = iter->seg[1].t_seg_start;
				}
			}
			else
			{
				// the curve is not flipped wrt x
				iter->seg[0].flipped = 0;

				if (iter->seg[1].flipped)
				{
					start_x = (int)iter->seg[1].curve.p[PTS-1][0];
					iter->seg[0].end_cell = (int)iter->seg[1].curve.p[0][0];
					if (iter->seg[1].curve.p[0][0] == floor(iter->seg[1].curve.p[0][0]))
						iter->seg[0].end_cell--;
					iter->seg[0].curve.p[PTS-1] = iter->seg[1].curve.p[PTS-1];
					iter->seg[0].t_seg_end = iter->seg[1].t_seg_end;
					iter->seg[0].t_curve_end = iter->seg[1].t_seg_start;
				}
				else
				{
					start_x = (int)iter->seg[1].curve.p[0][0];
					iter->seg[0].end_cell = (int)iter->seg[1].curve.p[PTS-1][0];
					if (iter->seg[1].curve.p[PTS-1][0] == floor(iter->seg[1].curve.p[PTS-1][0]))
						iter->seg[0].end_cell--;
					iter->seg[0].curve.p[PTS-1] = iter->seg[1].curve.p[0];
					iter->seg[0].t_seg_end = iter->seg[1].t_seg_start;
					iter->seg[0].t_curve_end = iter->seg[1].t_seg_end;
				}
			}
			
			if (start_x < 0)
				start_x = 0;
			if (start_x > ptab.buckets.s - W)
				start_x = ptab.buckets.s - W;
			if (iter->seg[0].end_cell < start_x)
				iter->seg[0].end_cell = start_x;
			if (iter->seg[0].end_cell >= ptab.buckets.s)
				iter->seg[0].end_cell = ptab.buckets.s - 1;
			
			if (start_x < pixel_start)
				pixel_start = start_x;
			if (iter->seg[0].end_cell > pixel_end)
				pixel_end = iter->seg[0].end_cell;

			//-- cut to pixels
			cut_to_pixels(iter, start_x, curr_y);

			//-- update iterator
			iter = iter->next;
		}
	}

	void cut_to_pixels(IncrementalCurve<PTS, DIM> *inc, int curr_x, int curr_y)
	{
		while (curr_x <= inc->seg[0].end_cell)
		{
			//-- increment the curve
			inc->increment<0>(curr_x);
			Curve<PTS, DIM> &curve = inc->seg[0].curve;
			
			//-- store the curve
			BucketCurve<PTS, DIM> &bc = ptab.heap.push_back();

			// copy data over and put into local coordinates. x coordinates and order are reversed to correct for integration direction
			// if statements is more verbose than combining, but clearer and faster.
			if (inc->seg[0].flipped)
			{
				for (int ip = 0; ip < PTS; ip++)
				{
					bc.curve.p[ip].set(1 - (curve.p[ip][0] - curr_x), curve.p[ip][1] - curr_y);
					if (DIM == 3)
						bc.curve.p[ip][2] = curve.p[ip][2];
				}
			}
			else
			{
				for (int ip = 0; ip < PTS; ip++)
				{
					bc.curve.p[PTS-1-ip].set(1 - (curve.p[ip][0] - curr_x), curve.p[ip][1] - curr_y);
					if (DIM == 3)
						bc.curve.p[PTS-1-ip][2] = curve.p[ip][2];
				}
			}

			// stick at head of bucket list
			assert(curr_x >=0 && curr_x < ptab.buckets.s);
			bc.next = ptab.buckets[curr_x];
			ptab.buckets[curr_x] = &bc;

			//-- update iterator
			curr_x++;
		}
	}

	void eval_pixels(Array2D<vect<COLORCOMP,double> > &img, ColorFunction<COLORCOMP, COLORORDER> &color, void (**filter)(double *pixel, double *accum, Curve<PTS,DIM> *curve), int curr_y)
	{
		// init
		int curr_x = pixel_start;
		
		const int COLORCOEFFS = ColorComponent<COLORORDER>::NUM;
		double accum[W][W][COLORCOMP][COLORORDER];
		vect<COLORCOMP,double> *pixel[W];
		
		for (int i = 0; i < W; i++)
			for (int j = 0; j < W; j++)
				for (int k = 0; k < COLORCOMP; k++)
					for (int l = 0; l < COLORORDER; l++)
						accum[i][j][k][l] = 0;
		
		for (int y = 0; y < W; y++)
			pixel[y] = &img(curr_x, curr_y + y);
		
		// adjust color
		double color_const_save[COLORCOMP];
		if (COLORORDER == 2)
		{
			for (int k = 0; k < COLORCOMP; k++)
			{
				color_const_save[k] = color.comp[k].c[0];
				color.comp[k].c[0] += curr_x*color.comp[k].c[1] + curr_y*color.comp[k].c[2];
			}
		}
		
		// process each pixel in the scanline
		while (curr_x <= pixel_end)
		{
			assert(curr_x >= 0 && curr_x <= img.s[0] - W);
			assert(curr_y >= 0 && curr_y <= img.s[1] - W);

			// clear pixel buffer and add accum (need to add accum here to avoid writing past the end of the scanline)
			for (int y = 0; y < W; y++)
			{
				for (int x = 0; x < W; x++)
				{
					for (int c = 0; c < COLORCOMP; c++)
					{
						// add linear accum term
						if (COLORORDER == 2)
							accum[x][y][c][0] += accum[x][y][c][1];

						// add constant term
						pixel[y][x][c] += accum[x][y][c][0];
					}
				}
			}

			// evaluate curve segments
			BucketCurve<PTS, DIM> *bc = ptab.buckets[curr_x];
			while (bc)
			{
				for (int y = 0; y < W; y++)
				{
					for (int x = 0; x < W; x++)
					{
						double pbuf[COLORCOEFFS];
						double abuf[COLORCOEFFS*COLORORDER];

						// Box filter specializations because box is particularly simple and doesn't take much space.
						if (COLORORDER == 1 && W == 1 && PTS == 2)
						{	
							const double &p00 = bc->curve.p[0].v[0];const double &p01 = bc->curve.p[0].v[1];const double &p10 = bc->curve.p[1].v[0];const double &p11 = bc->curve.p[1].v[1];
							*pbuf = ((1.0/2.0)*(p00+p10)*((-1.0*p01)+p11));
							*abuf = ((-1.0*p01)+p11);
						}
						else if (COLORORDER == 1 && W == 1 && PTS == 3 && DIM == 2)
						{	
							const double &p00 = bc->curve.p[0].v[0];const double &p01 = bc->curve.p[0].v[1];const double &p10 = bc->curve.p[1].v[0];const double &p11 = bc->curve.p[1].v[1];const double &p20 = bc->curve.p[2].v[0];const double &p21 = bc->curve.p[2].v[1];
							*pbuf = ((1.0/6.0)*((-2.0*p11*p20)+(-1.0*p01*((2.0*p10)+p20))+(2.0*p10*p21)+(3.0*p20*p21)+(p00*((-3.0*p01)+(2.0*p11)+p21))));
							*abuf = ((-1.0*p01)+p21);
						}
						else if (COLORORDER == 1 && W == 1 && PTS == 4)
						{	
							const double &p00 = bc->curve.p[0].v[0];const double &p01 = bc->curve.p[0].v[1];const double &p10 = bc->curve.p[1].v[0];const double &p11 = bc->curve.p[1].v[1];const double &p20 = bc->curve.p[2].v[0];const double &p21 = bc->curve.p[2].v[1];const double &p30 = bc->curve.p[3].v[0];const double &p31 = bc->curve.p[3].v[1];
							*pbuf = ((1.0/20.0)*((-3.0*p11*p20)+(-3.0*p11*p30)+(-6.0*p21*p30)+(-1.0*p01*((6.0*p10)+(3.0*p20)+p30))+(6.0*p20*p31)+(10.0*p30*p31)+(3.0*p10*(p21+p31))+(p00*((-10.0*p01)+(6.0*p11)+(3.0*p21)+p31))));
							*abuf = (p31 - p01);
						}
						// general filters via function pointer tables.
						else
						{
							const int idx = filter_elem<W>(x,y);
							filter[idx](pbuf, abuf, &bc->curve);
						}

						// add buffers into pixels
						for (int c = 0; c < COLORCOMP; c++)
						{
							for (int k = 0; k < COLORCOEFFS; k++)
							{
								pixel[y][x][c] += pbuf[k] * color.comp[c].c[k];
								for (int l = 0; l < COLORORDER; l++)
									accum[x][y][c][l] += abuf[k*COLORORDER+l] * color.comp[c].c[k];
							}
						}
					}
				}
				bc = bc->next;
			}
			
			// progress
			if (COLORORDER > 1)
			{
				for (int c = 0; c < COLORCOMP; c++)
					color.comp[c].c[0] += color.comp[c].c[1];
			}

			ptab.buckets[curr_x] = 0;
			curr_x++;
			for (int y = 0; y < W; y++)
				pixel[y]++;
		}

		// restore color
		if (COLORORDER == 2)
		{
			for (int k = 0; k < COLORCOMP; k++)
				color.comp[k].c[0] = color_const_save[k];
		}
	}

};
