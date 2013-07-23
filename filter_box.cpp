#include "filter_box.h"

static void line_00(double *pixel, double *accum, Curve<2,2> *curve)
{
	const double &p00 = curve->p[0].v[0];
	const double &p01 = curve->p[0].v[1];
	const double &p10 = curve->p[1].v[0];
	const double &p11 = curve->p[1].v[1];
	*pixel = ((1.0/2.0)*(p00+p10)*((-1.0*p01)+p11));
	*accum = ((-1.0*p01)+p11);
}

static void quad_00(double *pixel, double *accum, Curve<3,2> *curve)
{
	const double &p00 = curve->p[0].v[0];
	const double &p01 = curve->p[0].v[1];
	const double &p10 = curve->p[1].v[0];
	const double &p11 = curve->p[1].v[1];
	const double &p20 = curve->p[2].v[0];
	const double &p21 = curve->p[2].v[1];
	*pixel = ((1.0/6.0)*((-2.0*p11*p20)+(-1.0*p01*((2.0*p10)+p20))+(2.0*p10*p21)+(3.0*p20*p21)+(p00*((-3.0*p01)+(2.0*p11)+p21))));
	*accum = ((-1.0*p01)+p21);
}

static void cube_00(double *pixel, double *accum, Curve<4,2> *curve)
{
	const double &p00 = curve->p[0].v[0];
	const double &p01 = curve->p[0].v[1];
	const double &p10 = curve->p[1].v[0];
	const double &p11 = curve->p[1].v[1];
	const double &p20 = curve->p[2].v[0];
	const double &p21 = curve->p[2].v[1];
	const double &p30 = curve->p[3].v[0];
	const double &p31 = curve->p[3].v[1];
	*pixel = ((1.0/20.0)*((-3.0*p11*p20)+(-3.0*p11*p30)+(-6.0*p21*p30)+(-1.0*p01*((6.0*p10)+(3.0*p20)+p30))+(6.0*p20*p31)+(10.0*p30*p31)+(3.0*p10*(p21+p31))+(p00*((-10.0*p01)+(6.0*p11)+(3.0*p21)+p31))));
	*accum = (p31 - p01);
}

static void rquad_00(double *pixel, double *accum, Curve<3,3> *curve)
{
	const double &p00 = curve->p[0].v[0];
	const double &p01 = curve->p[0].v[1];
	const double &w0 = curve->p[0].v[2];
	const double &p10 = curve->p[1].v[0];
	const double &p11 = curve->p[1].v[1];
	const double &w1 = curve->p[1].v[2];
	const double &p20 = curve->p[2].v[0];
	const double &p21 = curve->p[2].v[1];
	const double &w2 = curve->p[2].v[2];
	
	double aaa = w0*w2 - w1*w1;

	if (aaa > 0.001)
	{
		double sqrt_aaa = sqrt(aaa);
		*pixel = (sqrt_aaa*((p01*p10 + p00*(p01 - p11) + p11*p20 - (p10 + p20)*p21)*w1*w1 - (p00 + p20)*(p01 - p21)*w0*w2) + 
				(p01*(p10 - p20) + p11*p20 - p10*p21 + p00*(-p11 + p21))*w0*w1*w2*(atan((-w0 + w1)/sqrt_aaa) + atan((w1 - w2)/sqrt_aaa)))/
					(2.*sqrt_aaa*aaa);
	}
	else
	{
		// If the weights are too close to constant the rational formula breaks, so revert to quadratic bezier.
		*pixel = ((1.0/6.0)*((-2.0*p11*p20)+(-1.0*p01*((2.0*p10)+p20))+(2.0*p10*p21)+(3.0*p20*p21)+(p00*((-3.0*p01)+(2.0*p11)+p21))));
	}

	*accum = (p21 - p01);
}

template<>
FilterBox<1>::FilterBox()
{
	filter22[0] = line_00;
	filter32[0] = quad_00;
	filter42[0] = cube_00;
	filter33[0] = rquad_00;
}
