#include <stdio.h>
#include <string>
#include "rasterizer.h"
#include "filter_box.h"
#include "save_ppm.h"

void create_bezigon(const double * data, const int n, Array<Curve<4,2> > &cubics)
{
	for (int i = 0; i < n; i += 6) {
		Curve<4,2> cubic;
		cubic.p[0].set(data[i],			data[i+1]);
		cubic.p[1].set(data[i+2],		data[i+3]);
		cubic.p[2].set(data[i+4],		data[i+5]);
		cubic.p[3].set(data[(i+6)%n],	data[(i+7)%n]);
		cubics.push_back(cubic);
	}
}

int main(int argc, char **argv)
{
	//-- filter
	typedef FilterBox<1> FilterType;
	FilterType filter;

	int img_w = 8;
	int img_h = img_w;

	//-- curves
	Array<Curve<2,2> > lines;
	Array<Curve<3,2> > quads;
	Array<Curve<4,2> > cubics;
	Array<Curve<3,3> > rquads;

	const int n = 12;
	double data[] = {1,1, 7,1, 7,7, 1,7, 1,3, 1,2};
	create_bezigon(data, n, cubics);
		
	//-- image
	Array2D<vect1d> img;
	int extra = (filter.W >> 1) << 1;
	img.resize(img_w + extra, img_h + extra);
	vect1d zero; zero = 0;
	img = zero;

	//-- color
	ColorFunction<1,1> color;
	color.comp[0].c[0] = 1;
		
	//-- rasterize
	RasterizerInstance<FilterType::W, 1, 1> inst(img, lines.s, quads.s, cubics.s, rquads.s);
	rasterize(inst, color, &lines, &filter.filter22[0], &quads, &filter.filter32[0], &cubics, &filter.filter42[0], &rquads, &filter.filter33[0]);
	
	//-- display
	for (int i = 0; i < img_h; i++) {
		for (int j = 0; j < img_w; j++)
			printf("%f ", img.data[i*img_w + j]);
		printf("\n");
	}

	//-- save
	save_pgm(img, "output.pgm", FilterType::W);

	return 0;
}
