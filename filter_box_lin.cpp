#include "filter_box.h"

static void line_00(double *pixel, double *accum, Curve<2,2> *curve)
{
	const double &p00 = curve->p[0].v[0];
	const double &p01 = curve->p[0].v[1];
	const double &p10 = curve->p[1].v[0];
	const double &p11 = curve->p[1].v[1];
	const double p002 = p00*p00;
	const double p012 = p01*p01;
	const double p102 = p10*p10;
	const double p112 = p11*p11;
	pixel[0] = ((1.0/2.0)*(p00+p10)*((-1.0*p01)+p11));
	accum[0] = ((-1.0*p01)+p11);
	accum[1] = 0.0;
	pixel[1] = ((1.0/6.0)*(p002+(p00*p10)+p102)*((-1.0*p01)+p11));
	accum[2] = ((1.0/2.0)*(p01+(-1.0*p11)));
	accum[3] = ((-1.0*p01)+p11);
	pixel[2] = ((1.0/6.0)*((-1.0*p01)+p11)*((p00*((2.0*p01)+p11))+(p10*(p01+(2.0*p11)))));
	accum[4] = ((-1.0/2.0)*(p01+(-1.0*p11))*(p01+p11));
	accum[5] = 0.0;
}

static void quad_00(double *pixel, double *accum, Curve<3,2> *curve)
{
	const double &p00 = curve->p[0].v[0];
	const double &p01 = curve->p[0].v[1];
	const double &p10 = curve->p[1].v[0];
	const double &p11 = curve->p[1].v[1];
	const double &p20 = curve->p[2].v[0];
	const double &p21 = curve->p[2].v[1];
	const double p002 = p00*p00;
	const double p012 = p01*p01;
	const double p102 = p10*p10;
	const double p112 = p11*p11;
	const double p202 = p20*p20;
	const double p212 = p21*p21;
	pixel[0] = ((1.0/6.0)*((-2.0*p11*p20)+(-1.0*p01*((2.0*p10)+p20))+(2.0*p10*p21)+(3.0*p20*p21)+(p00*((-3.0*p01)+(2.0*p11)+p21))));
	accum[0] = ((-1.0*p01)+p21);
	accum[1] = 0.0;
	pixel[1] = ((1.0/30.0)*((-2.0*p10*p11*p20)+(-4.0*p11*p202)+(-1.0*p01*((2.0*p102)+(2.0*p10*p20)+p202))+(2.0*p102*p21)+(4.0*p10*p20*p21)+(5.0*p202*p21)+(p002*((-5.0*p01)+(4.0*p11)+p21))+(p00*((-1.0*p01*((4.0*p10)+p20))+(p20*p21)+(2.0*p10*(p11+p21))))));
	accum[2] = ((1.0/2.0)*(p01+(-1.0*p21)));
	accum[3] = ((-1.0*p01)+p21);
	pixel[2] = ((1.0/30.0)*((-1.0*p012*((4.0*p10)+p20))+(2.0*((-1.0*p112*p20)+(p11*(p10+(-2.0*p20))*p21)+(((2.0*p10)+(5.0*p20))*p212)))+(-1.0*p01*((2.0*p10*p11)+(p20*((2.0*p11)+p21))))+(p00*((-10.0*p012)+(2.0*p112)+(2.0*p11*p21)+p212+(p01*((4.0*p11)+p21))))));
	accum[4] = ((1.0/2.0)*((-1.0*p012)+p212));
	accum[5] = 0.0;
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
	const double p002 = p00*p00;
	const double p012 = p01*p01;
	const double p102 = p10*p10;
	const double p112 = p11*p11;
	const double p202 = p20*p20;
	const double p212 = p21*p21;
	const double p302 = p30*p30;
	const double p312 = p31*p31;
	pixel[0] = ((1.0/20.0)*((-3.0*p11*p20)+(3.0*p10*p21)+(-3.0*p11*p30)+(-6.0*p21*p30)+(-1.0*p01*((6.0*p10)+(3.0*p20)+p30))+(3.0*p10*p31)+(6.0*p20*p31)+(10.0*p30*p31)+(p00*((-10.0*p01)+(6.0*p11)+(3.0*p21)+p31))));
	accum[0] = ((-1.0*p01)+p31);
	accum[1] = 0.0;
	pixel[1] = ((1.0/840.0)*((-27.0*p10*p11*p20)+(-27.0*p11*p202)+(27.0*p102*p21)+(27.0*p10*p20*p21)+(-18.0*p10*p11*p30)+(-45.0*p11*p20*p30)+(-45.0*p20*p21*p30)+(-30.0*p11*p302)+(-105.0*p21*p302)+(-1.0*p01*((45.0*p102)+(45.0*p10*p20)+(18.0*p202)+(12.0*p10*p30)+(15.0*p20*p30)+(5.0*p302)))+(18.0*p102*p31)+(45.0*p10*p20*p31)+(45.0*p202*p31)+(30.0*p10*p30*p31)+(105.0*p20*p30*p31)+(140.0*p302*p31)+(5.0*p002*((-28.0*p01)+(21.0*p11)+(6.0*p21)+p31))+(p00*((18.0*p20*p21)+(-3.0*p11*p30)+(3.0*p21*p30)+(-5.0*p01*((21.0*p10)+(6.0*p20)+p30))+(12.0*p20*p31)+(5.0*p30*p31)+(15.0*p10*((3.0*p11)+(3.0*p21)+p31))))));
	accum[2] = ((1.0/2.0)*(p01+(-1.0*p31)));
	accum[3] = ((-1.0*p01)+p31);
	pixel[2] = ((1.0/840.0)*((-27.0*p112*p20)+(27.0*p10*p11*p21)+(-27.0*p11*p20*p21)+(27.0*p10*p212)+(-18.0*p112*p30)+(-45.0*p11*p21*p30)+(-45.0*p212*p30)+(-5.0*p012*((21.0*p10)+(6.0*p20)+p30))+(18.0*p10*p11*p31)+(45.0*p10*p21*p31)+(45.0*p20*p21*p31)+(-30.0*p11*p30*p31)+(-105.0*p21*p30*p31)+(30.0*p10*p312)+(105.0*p20*p312)+(280.0*p30*p312)+(-1.0*p01*((45.0*p10*p11)+(45.0*p11*p20)+(18.0*p20*p21)+(15.0*p11*p30)+(12.0*p21*p30)+(-3.0*p10*p31)+(3.0*p20*p31)+(5.0*p30*p31)))+(p00*((-280.0*p012)+(45.0*p112)+(45.0*p11*p21)+(18.0*p212)+(12.0*p11*p31)+(15.0*p21*p31)+(5.0*p312)+(5.0*p01*((21.0*p11)+(6.0*p21)+p31))))));
	accum[4] = ((1.0/2.0)*((-1.0*p012)+p312));
	accum[5] = 0.0;
}

static void rquad(double *pixel, double *accum, Curve<3,3> *curve)
{
}

template<>
FilterBox<2>::FilterBox()
{
	filter22[0] = line_00;
	filter32[0] = quad_00;
	filter42[0] = cube_00;
	filter33[0] = rquad;
}
