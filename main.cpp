#include <windows.h>
#include <stdio.h>
#include <string>
#include <fstream>
#include "rasterizer.h"
#include "rasterizer_parallel.h"
#include "filter_box.h"
#include "filter_tent.h"
#include "filter_bspline2.h"
#include "filter_mitchell.h"
#include "filter_lanczos2.h"
#include "filter_lanczos3.h"
#include "filter_radial2.h"
#include "filter_radial3.h"
#include "timer.h"
#include "save_ppm.h"

using namespace std;

void create_svg(Array<Curve<2,2> > &lines, Array<Curve<3,2> > &quads, Array<Curve<4,2> > &cubics, Array<Curve<3,3> > &rquads);

template <int PTS, int DIM>
void scale_curves(Array<Curve<PTS, DIM> > &curves, double scale)
{
	for (int i = 0; i < curves.s; i++)
	{
		Curve<PTS, DIM> &c = curves[i];
		for (int j = 0; j < PTS; j++)
		{
			c.p[j][0] = c.p[j][0] * scale;
			c.p[j][1] = c.p[j][1] * scale;
		}
	}
}

template <int COMP, int ORDER>
void scale_color(ColorFunction<COMP,ORDER> &color, double scale)
{
	for (int i = 0; i < COMP; i++)
	{
		if (ORDER > 1)
		{
			color.comp[i].c[1] /= scale;
			color.comp[i].c[2] /= scale;
		}
	}
}

void add_circle(Array<Curve<3,3> > &lines, vect2d center, double radius)
{
	Curve<3,3> l;
	l.p[0][2] = 1;
	l.p[1][2] = sqrt(.5);
	l.p[2][2] = 1;

	if (radius >= 0)
	{
		l.p[0].set(center[0]+radius,center[1]);
		l.p[1].set(center[0]+radius,center[1]+radius);
		l.p[2].set(center[0],center[1]+radius);
		lines.push_back(l);

		l.p[0].set(center[0],center[1]+radius);
		l.p[1].set(center[0]-radius,center[1]+radius);
		l.p[2].set(center[0]-radius,center[1]);
		lines.push_back(l);

		l.p[0].set(center[0]-radius,center[1]);
		l.p[1].set(center[0]-radius,center[1]-radius);
		l.p[2].set(center[0],center[1]-radius);
		lines.push_back(l);

		l.p[0].set(center[0],center[1]-radius);
		l.p[1].set(center[0]+radius,center[1]-radius);
		l.p[2].set(center[0]+radius,center[1]);
		lines.push_back(l);
	}
	else
	{
		l.p[2].set(center[0]+radius,center[1]);
		l.p[1].set(center[0]+radius,center[1]+radius);
		l.p[0].set(center[0],center[1]+radius);
		lines.push_back(l);

		l.p[2].set(center[0],center[1]+radius);
		l.p[1].set(center[0]-radius,center[1]+radius);
		l.p[0].set(center[0]-radius,center[1]);
		lines.push_back(l);

		l.p[2].set(center[0]-radius,center[1]);
		l.p[1].set(center[0]-radius,center[1]-radius);
		l.p[0].set(center[0],center[1]-radius);
		lines.push_back(l);

		l.p[2].set(center[0],center[1]-radius);
		l.p[1].set(center[0]+radius,center[1]-radius);
		l.p[0].set(center[0]+radius,center[1]);
		lines.push_back(l);
	}
}

void create_apollonian(Array<Curve<2,2> > &lines, Array<Curve<3,2> > &quads, Array<Curve<4,2> > &cubics, Array<Curve<3,3> > &rquads, double img_size)
{
	ifstream f("data/apollonian.txt");
	for (int i = 0; ; i++)
	{
		int n1,d1,n2,d2,r;
		f >> n1 >> d1 >> n2 >> d2 >> r;

		if (f.eof())
			break;

		double radius = 1.0 / r;
		vect2d center;
		center[0] = (double)n1 / d1;
		center[1] = (double)n2 / d2;

		if (i > 0) add_circle(rquads, center, .875*radius);
	}
	scale_curves(rquads, img_size);
	f.close();
}

void create_circle(Array<Curve<2,2> > &lines, Array<Curve<3,2> > &quads, Array<Curve<4,2> > &cubics, Array<Curve<3,3> > &rquads, double img_size)
{
	double radius = 2;
	vect2d center;
	center[0] = 2.25;
	center[1] = 2.5;
	add_circle(rquads, center, radius);

	scale_curves(rquads, img_size / 5);
}


void add_diamond(Array<Curve<2,2> > &lines, vect2d center, double radius)
{
	Curve<2,2> line;

	line.p[0].set(center[0],-radius+center[1]);
	line.p[1].set(radius+center[0],center[1]);
	lines.push_back(line);

	line.p[0].set(radius+center[0],center[1]);
	line.p[1].set(center[0],radius+center[1]);
	lines.push_back(line);

	line.p[0].set(center[0],radius+center[1]);
	line.p[1].set(-radius+center[0],center[1]);
	lines.push_back(line);

	line.p[0].set(-radius+center[0],center[1]);
	line.p[1].set(center[0],-radius+center[1]);
	lines.push_back(line);
}

void create_diamond(Array<Curve<2,2> > &lines, Array<Curve<3,2> > &quads, Array<Curve<4,2> > &cubics, Array<Curve<3,3> > &rquads, double img_size)
{
	vect2d center;
	center.set(.25, .25);
	add_diamond(lines, center, .2);
	center.set(.75, .25);
	add_diamond(lines, center, .2);
	center.set(.25, .75);
	add_diamond(lines, center, .2);
	center.set(.75, .75);
	add_diamond(lines, center, .2);

	scale_curves(lines, img_size);
}


void create_radial(Array<Curve<2,2> > &lines, Array<Curve<3,2> > &quads, Array<Curve<4,2> > &cubics, Array<Curve<3,3> > &rquads, double img_size)
{
	int n = 128;
	const double pi = 3.1415926535;
	for (int i = 0; i < n; i++)
	{
		double t1 = i * 2.0 * pi / n;
		double t2 = (i + .5) * 2.0 * pi / n;

		vect3d p[4];
		p[0].set(.5, .5);
		p[1].set(.5 + cos(t1)*.49, .5 + sin(t1)*.49);
		p[2].set(.5 + cos(t2)*.49, .5 + sin(t2)*.49);
		p[3].set(.5, .5);

		for (int ip = 0; ip < 3; ip++)
		{
			Curve<2,2> line;
			line.p[0] = p[ip];
			line.p[1] = p[ip+1];
			lines.push_back(line);
		}
	}

	scale_curves(lines, img_size);
}

void create_radial_paper(Array<Curve<2,2> > &lines, Array<Curve<3,2> > &quads, Array<Curve<4,2> > &cubics, Array<Curve<3,3> > &rquads, double img_size)
{
	int n = 128;
	const double pi = 3.1415926535;
	for (int i = 0; i < n/4; i++)
	{
		double t1 = i * 2.0 * pi / n;
		double t2 = (i + .5) * 2.0 * pi / n;

		vect3d p[4];
		p[0].set(0, 0);
		p[1].set(cos(t1)*.98, sin(t1)*.98);
		p[2].set(cos(t2)*.98, sin(t2)*.98);
		p[3].set(0, 0);

		for (int ip = 0; ip < 3; ip++)
		{
			Curve<2,2> line;
			line.p[0] = p[ip];
			line.p[1] = p[ip+1];
			lines.push_back(line);
		}
	}

	scale_curves(lines, img_size);
}

void create_quadratic_arc(Array<Curve<2,2> > &lines, Array<Curve<3,2> > &quads, Array<Curve<4,2> > &cubics, Array<Curve<3,3> > &rquads, double img_size)
{
	Curve<2,2> a;
	Curve<3,2> b;

	a.p[0].set(0,1);
	a.p[1].set(0,0);
	lines.push_back(a);

	b.p[0].set(1,0);
	b.p[1].set(1,1);
	b.p[2].set(0,1);
	quads.push_back(b);

	scale_curves(lines, img_size);
	scale_curves(quads, img_size);
}

void create_cubic_arc(Array<Curve<2,2> > &lines, Array<Curve<3,2> > &quads, Array<Curve<4,2> > &cubics, Array<Curve<3,3> > &rquads, double img_size)
{
	// cubic approximation of quadrant of circle
	Curve<2,2> a;
	a.p[0].set(0,1);
	a.p[1].set(0,0);
	lines.push_back(a);

	Curve<4,2> cubic;
	double cw = 2.0 * 0.551784; // least squares approximation of a circle quadrant
	cubic.p[0].set(2.25,0.5);
	cubic.p[1].set(2.25 + cw,0.5);
	cubic.p[2].set(4.25,2.5 - cw);
	cubic.p[3].set(4.25,2.5);
	cubics.push_back(cubic);

	scale_curves(lines, img_size);
	scale_curves(cubics, img_size);
}

void create_zone_paper(Array<Curve<2,2> > &lines, Array<Curve<3,2> > &quads, Array<Curve<4,2> > &cubics, Array<Curve<3,3> > &rquads, double img_size)
{
	Array<Curve<3,2> > proto;
	Curve<3,2> a;
	Curve<2,2> b;

	//b.p[0].set(1,0);
	//b.p[1].set(1,1);
	//b.p[2].set(0,1);
	//proto.push_back(b);

	const double s2 = sqrt(.5);
	const double S2 = sqrt(2.0);

	a.p[0].set(1,0);
	a.p[1].set(1,S2-1);
	a.p[2].set(s2,s2);
	proto.push_back(a);

	a.p[0].set(s2,s2);
	a.p[1].set(S2-1,1);
	a.p[2].set(0,1);
	proto.push_back(a);

	for (int i = 1; i < 64; i+=2)
	{
		double w = .003*4;
		double factor1 = sqrt(i*w + w*w*i*i*.25);
		double factor2 = sqrt((i+1)*w + w*w*(i+1)*(i+1)*.25);
		if (sqrt((i+2)*w + w*w*(i+2)*(i+2)*.25) > 1)
			break;

		// outside
		a = proto[0];
		for (int k = 0; k < 3; k++)
			a.p[k] = a.p[k] * factor2;
		//swap(a.p[0], a.p[2]);
		quads.push_back(a);

		a = proto[1];
		for (int k = 0; k < 3; k++)
			a.p[k] = a.p[k] * factor2;
		//swap(a.p[0], a.p[2]);
		quads.push_back(a);

		b.p[0] = proto[1].p[2] * factor2;
		b.p[1] = proto[1].p[2] * factor1;
		lines.push_back(b);


		// inside
		a = proto[1];
		for (int k = 0; k < 3; k++)
			a.p[k] = a.p[k] * factor1;
		swap(a.p[0], a.p[2]);
		quads.push_back(a);

		a = proto[0];
		for (int k = 0; k < 3; k++)
			a.p[k] = a.p[k] * factor1;
		swap(a.p[0], a.p[2]);
		quads.push_back(a);

		b.p[0] = proto[0].p[0] * factor1;
		b.p[1] = proto[0].p[0] * factor2;
		lines.push_back(b);
	}

	scale_curves(lines, img_size);
	scale_curves(quads, img_size);
}

void poly_to_bound(Array<Curve<2,2> > &lines, Array<vect2d> &poly)
{
	for (int i = 0; i < poly.s; i++)
	{
		Curve<2,2> c;
		c.p[0] = poly[i];
		c.p[1] = poly[(i+1)%poly.s];
		lines.push_back(c);
	}
}

void create_perspective_paper(Array<Curve<2,2> > &lines, Array<Curve<3,2> > &quads, Array<Curve<4,2> > &cubics, Array<Curve<3,3> > &rquads, double img_size)
{
	Array<vect2d> poly;
	vect2d pt;

	int fake_img_size = 128;
	double px = 1.0 / fake_img_size;
	int n = fake_img_size * 2;
	double dx = 1.0/n;
	double dx2 = dx * 8;
	for (int i = 0; i < n; i++)
	{
		poly.s = 0;

		if ((i+.5)*dx2 <= 1)
		{
			pt.set(i*dx, 0); poly.push_back(pt);
			pt.set((i+.5)*dx, 0); poly.push_back(pt);
			pt.set((i+.5)*dx2, 1); poly.push_back(pt);
			pt.set(i*dx2, 1); poly.push_back(pt);
			poly_to_bound(lines, poly);
		}
		else
		{
			pt.set(i*dx, 0); poly.push_back(pt);
			pt.set((i+.5)*dx, 0); poly.push_back(pt);

			{
				double x1 = (i+.5)*dx;
				double x2 = (i+.5)*dx2;
				double y1 = 0;
				double y2 = 1;
				double t = (1 - x1) / (x2 - x1);
				if (y1 + t*(y2 - y1) <= 0)
					continue;
				pt.set(1, y1 + t*(y2 - y1)); poly.push_back(pt);
			}
			{
				double x1 = i*dx;
				double x2 = i*dx2;
				double y1 = 0;
				double y2 = 1;
				double t = (1 - x1) / (x2 - x1);
				if (y1 + t*(y2 - y1) <= 0)
					continue;
				pt.set(1, y1 + t*(y2 - y1)); poly.push_back(pt);
			}

			poly_to_bound(lines, poly);
		}
	}

	scale_curves(lines, img_size);
}

void set_color(ColorFunction<1,1> &color)
{
	color.comp[0].c[0] = 1;
}

void set_color(ColorFunction<1,2> &color)
{
	color.comp[0].c[0] = 0;
	color.comp[0].c[1] = 1;
	color.comp[0].c[2] = 0;
}

void set_color(ColorFunction<3,1> &color)
{
	color.comp[0].c[0] = 0;
	color.comp[1].c[0] = 0;
	color.comp[2].c[0] = 1;
}

void set_color(ColorFunction<3,2> &color)
{
	color.comp[0].c[0] = 0;
	color.comp[0].c[1] = 1;
	color.comp[0].c[2] = 0;
	color.comp[1].c[0] = 0;
	color.comp[1].c[1] = 0;
	color.comp[1].c[2] = 1;
	color.comp[2].c[0] = 1;
	color.comp[2].c[1] = 0;
	color.comp[2].c[2] = 0;
}

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
	int timing_iters = 5;

	if (argc > 1)
		timing_iters = atoi(argv[1]);

	//freopen("stdout.txt", "w", stdout);
	//SetPriorityClass(GetCurrentProcess(), REALTIME_PRIORITY_CLASS);

	//-- filter
	typedef FilterBox<1> FilterType;
	FilterType filter;


	//for (int n_test = 0; n_test < 3; n_test++)
	//int n_test = 1;
	{
		int img_w = 1280;
		int img_h = img_w;

		//-- curves
		Array<Curve<2,2> > lines;
		Array<Curve<3,2> > quads;
		Array<Curve<4,2> > cubics;
		Array<Curve<3,3> > rquads;

		const int n = 12;
		double data[] = {100,100, 700,100, 700,700, 100,700, 100,300, 100,200};
		create_bezigon(data, n, cubics);
		
		//-- image
		Array2D<vect1d> img;
		int extra = (filter.W >> 1) << 1;
		img.resize(img_w + extra, img_h + extra);
		vect1d zero; zero = 0;
		img = zero;

		//-- color
		ColorFunction<1,1> color; set_color(color);
		scale_color(color, img_w);
		
		//-- rasterize
		RasterizerInstance<FilterType::W, 1, 1> inst(img, lines.s, quads.s, cubics.s, rquads.s);
		rasterize(inst, color, &lines, &filter.filter22[0], &quads, &filter.filter32[0], &cubics, &filter.filter42[0], &rquads, &filter.filter33[0]);
		//rasterize<FilterType::W>(img, color, &lines, &filter.filter22[0], &quads, &filter.filter32[0], &cubics, &filter.filter42[0], &rquads, &filter.filter33[0]);

		//-- save
		flip_img(img);
		if (img.data[0].DIMENSION == 1)
			save_pgm(img, "output.pgm", FilterType::W);
		else if (img.data[0].DIMENSION == 3)
			save_ppm(img, "output.ppm", FilterType::W);
	}

	return 0;
}
