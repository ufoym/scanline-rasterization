#include <stdio.h>
#include <string>
#include "rasterizer.h"
#include "filter_box.h"


#define DLLEXPORT extern "C" __declspec(dllexport)

DLLEXPORT double * rasterize_bezigon(const double * data, const int n, 
									 const int w, const int h)
{
	//-- filter
	typedef FilterBox<1> FilterType;
	FilterType filter;

	//-- curves
	Array<Curve<2,2> > lines;
	Array<Curve<3,2> > quads;
	Array<Curve<4,2> > cubics;
	Array<Curve<3,3> > rquads;

	for (int i = 0; i < n; i += 6) {
		Curve<4,2> cubic;
		cubic.p[0].set(data[i],			data[i+1]);
		cubic.p[1].set(data[i+2],		data[i+3]);
		cubic.p[2].set(data[i+4],		data[i+5]);
		cubic.p[3].set(data[(i+6)%n],	data[(i+7)%n]);
		cubics.push_back(cubic);
	}

	//-- image
	Array2D<vect1d> img;
	int extra = (filter.W >> 1) << 1;
	img.resize(w + extra, h + extra);
	vect1d zero; zero = 0;
	img = zero;

	//-- color
	ColorFunction<1,1> color;
	color.comp[0].c[0] = 1;

	//-- rasterize
	RasterizerInstance<FilterType::W, 1, 1> inst(
		img, lines.s, quads.s, cubics.s, rquads.s);
	rasterize(inst, color, &lines, &filter.filter22[0], 
		&quads, &filter.filter32[0], &cubics, &filter.filter42[0], 
		&rquads, &filter.filter33[0]);

	//-- output
	double * alpha_arr = new double [h * w];
	for (int i = 0; i < h; i++) {
		for (int j = 0; j < w; j++) {
			int idx = i*w + j;
			alpha_arr[idx] = img.data[idx].getitem(0);
		}
	}
	return alpha_arr;
}

int main(int argc, char **argv)
{
	const int n = 12, h = 8, w = 8;
	double data[] = {1,1, 7,1, 7,7, 1,7, 1,3, 1,2};
	double * alpha_arr = rasterize_bezigon(data, n, w, h);
	for (int i = 0; i < h * w; i++)
		printf("%f ", alpha_arr[i]);
	return 0;
}
