#include "filter_tent.h"

//=============================================================================//
//================================= Lines =====================================//
//=============================================================================//

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
	*pixel = ((1.0/24.0)*((-1.0*p01)+p11)*((2.0*p00*p10*(p01+p11))+(p002*((3.0*p01)+p11))+(p102*(p01+(3.0*p11)))));
	*accum = ((-1.0/4.0)*(p01+(-1.0*p11))*(p01+p11));
}
static void line_10(double *pixel, double *accum, Curve<2,2> *curve)
{
	const double &p00 = curve->p[0].v[0];
	const double &p01 = curve->p[0].v[1];
	const double &p10 = curve->p[1].v[0];
	const double &p11 = curve->p[1].v[1];
	const double p002 = p00*p00;
	const double p012 = p01*p01;
	const double p102 = p10*p10;
	const double p112 = p11*p11;
	*pixel = ((1.0/24.0)*(p01+(-1.0*p11))*((p002*((3.0*p01)+p11))+(2.0*p00*((p01*(-4.0+p10))+((-2.0+p10)*p11)))+(p10*((p01*(-4.0+p10))+((-8.0+(3.0*p10))*p11)))));
	*accum = ((-1.0/4.0)*(p01+(-1.0*p11))*(p01+p11));
}
static void line_01(double *pixel, double *accum, Curve<2,2> *curve)
{
	const double &p00 = curve->p[0].v[0];
	const double &p01 = curve->p[0].v[1];
	const double &p10 = curve->p[1].v[0];
	const double &p11 = curve->p[1].v[1];
	const double p002 = p00*p00;
	const double p012 = p01*p01;
	const double p102 = p10*p10;
	const double p112 = p11*p11;
	*pixel = ((1.0/24.0)*(p01+(-1.0*p11))*((2.0*p00*p10*(-2.0+p01+p11))+(p002*(-4.0+(3.0*p01)+p11))+(p102*(-4.0+p01+(3.0*p11)))));
	*accum = ((1.0/4.0)*(p01+(-1.0*p11))*(-2.0+p01+p11));
}
static void line_11(double *pixel, double *accum, Curve<2,2> *curve)
{
	const double &p00 = curve->p[0].v[0];
	const double &p01 = curve->p[0].v[1];
	const double &p10 = curve->p[1].v[0];
	const double &p11 = curve->p[1].v[1];
	const double p002 = p00*p00;
	const double p012 = p01*p01;
	const double p102 = p10*p10;
	const double p112 = p11*p11;
	*pixel = ((-1.0/24.0)*(p01+(-1.0*p11))*((2.0*p00*(6.0+(p01*(-4.0+p10))+(p10*(-2.0+p11))+(-2.0*p11)))+(p002*(-4.0+(3.0*p01)+p11))+(p10*(12.0+(p01*(-4.0+p10))+(-8.0*p11)+(p10*(-4.0+(3.0*p11)))))));
	*accum = ((1.0/4.0)*(p01+(-1.0*p11))*(-2.0+p01+p11));
}


//=============================================================================//
//=============================== Quadratics ==================================//
//=============================================================================//

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
	*pixel = ((1.0/840.0)*((-16.0*p10*p112*p20)+(-20.0*p112*p202)+(-1.0*p012*((20.0*p102)+(12.0*p10*p20)+(3.0*p202)))+(16.0*p102*p11*p21)+(-60.0*p11*p202*p21)+(20.0*p102*p212)+(60.0*p10*p20*p212)+(105.0*p202*p212)+(p002*((-105.0*p012)+(20.0*p112)+(12.0*p11*p21)+(3.0*p212)+(10.0*p01*((6.0*p11)+p21))))+(-2.0*p01*((8.0*p102*p11)+(4.0*p10*p20*((3.0*p11)+p21))+(p202*((6.0*p11)+(5.0*p21)))))+(2.0*p00*((-5.0*p012*((6.0*p10)+p20))+(p20*p21*((4.0*p11)+(5.0*p21)))+(p01*((-4.0*p11*p20)+(4.0*p10*p21)))+(2.0*p10*((4.0*p112)+(6.0*p11*p21)+(3.0*p212)))))));
	*accum = ((1.0/4.0)*((-1.0*p012)+p212));
}
static void quad_10(double *pixel, double *accum, Curve<3,2> *curve)
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
	*pixel = ((1.0/840.0)*((-56.0*p112*p20)+(16.0*p10*p112*p20)+(20.0*p112*p202)+(p012*((20.0*p102)+(4.0*p10*(-28.0+(3.0*p20)))+(p20*(-28.0+(3.0*p20)))))+(56.0*p10*p11*p21)+(-16.0*p102*p11*p21)+(-112.0*p11*p20*p21)+(60.0*p11*p202*p21)+(112.0*p10*p212)+(-20.0*p102*p212)+(280.0*p20*p212)+(-60.0*p10*p20*p212)+(-105.0*p202*p212)+(p002*((105.0*p012)+(-20.0*p112)+(-12.0*p11*p21)+(-3.0*p212)+(-10.0*p01*((6.0*p11)+p21))))+(2.0*p00*(((28.0+(-8.0*p10))*p112)+(5.0*p012*(-28.0+(6.0*p10)+p20))+(-4.0*p11*(-7.0+(3.0*p10)+p20)*p21)+((14.0+(-6.0*p10)+(-5.0*p20))*p212)+(2.0*p01*((2.0*p11*(14.0+p20))+((7.0+(-2.0*p10))*p21)))))+(2.0*p01*((8.0*p102*p11)+(4.0*p10*((p11*(-7.0+(3.0*p20)))+(p20*p21)))+(p20*((p11*(-28.0+(6.0*p20)))+((-14.0+(5.0*p20))*p21)))))));
	*accum = ((1.0/4.0)*((-1.0*p012)+p212));
}
static void quad_01(double *pixel, double *accum, Curve<3,2> *curve)
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
	*pixel = ((1.0/840.0)*((-56.0*p10*p11*p20)+(16.0*p10*p112*p20)+(-112.0*p11*p202)+(20.0*p112*p202)+(p012*((20.0*p102)+(12.0*p10*p20)+(3.0*p202)))+(56.0*p102*p21)+(-16.0*p102*p11*p21)+(112.0*p10*p20*p21)+(140.0*p202*p21)+(60.0*p11*p202*p21)+(-20.0*p102*p212)+(-60.0*p10*p20*p212)+(-105.0*p202*p212)+(p002*((105.0*p012)+(-20.0*p112)+((28.0+(-3.0*p21))*p21)+(-10.0*p01*(14.0+(6.0*p11)+p21))+(-4.0*p11*(-28.0+(3.0*p21)))))+(2.0*p01*((4.0*p102*(-7.0+(2.0*p11)))+(4.0*p10*p20*(-7.0+(3.0*p11)+p21))+(p202*(-14.0+(6.0*p11)+(5.0*p21)))))+(2.0*p00*((5.0*p012*((6.0*p10)+p20))+(p20*(14.0+(-4.0*p11)+(-5.0*p21))*p21)+(-2.0*p01*(((7.0+(-2.0*p11))*p20)+(2.0*p10*(14.0+p21))))+(-2.0*p10*((4.0*p112)+(p21*(-14.0+(3.0*p21)))+(2.0*p11*(-7.0+(3.0*p21)))))))));
	*accum = ((1.0/4.0)*((-2.0*p01)+p012+(-1.0*(-2.0+p21)*p21)));
}
static void quad_11(double *pixel, double *accum, Curve<3,2> *curve)
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
	*pixel = ((1.0/840.0)*((-280.0*p11*p20)+(56.0*p10*p11*p20)+(56.0*p112*p20)+(-16.0*p10*p112*p20)+(112.0*p11*p202)+(-20.0*p112*p202)+(p012*((-20.0*p102)+((28.0+(-3.0*p20))*p20)+(-4.0*p10*(-28.0+(3.0*p20)))))+(280.0*p10*p21)+(-56.0*p102*p21)+(-56.0*p10*p11*p21)+(16.0*p102*p11*p21)+(420.0*p20*p21)+(-112.0*p10*p20*p21)+(112.0*p11*p20*p21)+(-140.0*p202*p21)+(-60.0*p11*p202*p21)+(-112.0*p10*p212)+(20.0*p102*p212)+(-280.0*p20*p212)+(60.0*p10*p20*p212)+(105.0*p202*p212)+(p002*((-105.0*p012)+(20.0*p112)+(10.0*p01*(14.0+(6.0*p11)+p21))+(4.0*p11*(-28.0+(3.0*p21)))+(p21*(-28.0+(3.0*p21)))))+(-2.0*p01*((4.0*p102*(-7.0+(2.0*p11)))+(4.0*p10*(35.0+(p11*(-7.0+(3.0*p20)))+(p20*(-7.0+p21))))+(p20*(70.0+(-14.0*p20)+(p11*(-28.0+(6.0*p20)))+(-14.0*p21)+(5.0*p20*p21)))))+(-2.0*p00*(((28.0+(-8.0*p10))*p112)+(5.0*p012*(-28.0+(6.0*p10)+p20))+(((p10*(28.0+(-6.0*p21)))+(p20*(14.0+(-5.0*p21)))+(14.0*(-5.0+p21)))*p21)+(2.0*p01*((2.0*p11*(14.0+p20))+(-2.0*p10*(14.0+p21))+(7.0*(15.0+(-1.0*p20)+p21))))+(-4.0*p11*(35.0+((-7.0+p20)*p21)+(p10*(-7.0+(3.0*p21)))))))));
	*accum = ((1.0/4.0)*((-2.0*p01)+p012+(-1.0*(-2.0+p21)*p21)));
}

//=============================================================================//
//================================= Cubics ====================================//
//=============================================================================//

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
	*pixel = ((1.0/18480.0)*((-1.0*p012*((378.0*p102)+(63.0*p202)+(30.0*p20*p30)+(5.0*p302)+(42.0*p10*((6.0*p20)+p30))))+(p002*((-2310.0*p012)+(378.0*p112)+(63.0*p212)+(30.0*p21*p31)+(5.0*p312)+(42.0*p11*((6.0*p21)+p31))+(28.0*p01*((45.0*p11)+(9.0*p21)+p31))))+(-2.0*p00*((-126.0*p10*p112)+(-189.0*p10*p11*p21)+(-54.0*p11*p20*p21)+(-81.0*p10*p212)+(-45.0*p20*p212)+(9.0*p112*p30)+(-9.0*p212*p30)+(14.0*p012*((45.0*p10)+(9.0*p20)+p30))+(-54.0*p10*p11*p31)+(-30.0*p11*p20*p31)+(-60.0*p10*p21*p31)+(-54.0*p20*p21*p31)+(-6.0*p11*p30*p31)+(-21.0*p21*p30*p31)+(-15.0*p10*p312)+(-21.0*p20*p312)+(-14.0*p30*p312)+(3.0*p01*((2.0*p21*p30)+(7.0*p11*((4.0*p20)+p30))+(-2.0*p20*p31)+(-7.0*p10*((4.0*p21)+p31))))))+(-2.0*p01*((45.0*p202*p21)+(54.0*p20*p21*p30)+(21.0*p21*p302)+(3.0*p11*((27.0*p202)+(20.0*p20*p30)+(5.0*p302)))+(9.0*p102*((14.0*p11)+(-1.0*p31)))+(9.0*p202*p31)+(21.0*p20*p30*p31)+(14.0*p302*p31)+(3.0*p10*((9.0*p11*((7.0*p20)+(2.0*p30)))+(2.0*((9.0*p20*p21)+(5.0*p21*p30)+(p30*p31)))))))+(3.0*((-3.0*p112*((15.0*p202)+(18.0*p20*p30)+(7.0*p302)))+(-2.0*p11*((27.0*p202*p21)+(42.0*p302*(p21+p31))+(7.0*p20*p30*((9.0*p21)+(4.0*p31)))))+(3.0*p102*((15.0*p212)+(18.0*p21*p31)+(7.0*p312)+(2.0*p11*((9.0*p21)+(5.0*p31)))))+(p10*((-6.0*p112*((9.0*p20)+(5.0*p30)))+(28.0*p30*p31*((2.0*p21)+(3.0*p31)))+(-36.0*p11*((p21*p30)+(-1.0*p20*p31)))+(6.0*p20*((9.0*p212)+(21.0*p21*p31)+(14.0*p312)))))+(14.0*((3.0*p202*p31*((2.0*p21)+(3.0*p31)))+(-6.0*p20*p30*(p212+(-5.0*p312)))+(p302*((-9.0*p212)+(-30.0*p21*p31)+(55.0*p312)))))))));
	*accum = ((1.0/4.0)*((-1.0*p012)+p312));
}
static void cube_10(double *pixel, double *accum, Curve<4,2> *curve)
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
	*pixel = ((1.0/18480.0)*((p012*((378.0*p102)+(63.0*p202)+(30.0*p20*p30)+(5.0*p302)+(42.0*p10*((6.0*p20)+p30))))+(-1.0*p002*((-2310.0*p012)+(378.0*p112)+(63.0*p212)+(30.0*p21*p31)+(5.0*p312)+(42.0*p11*((6.0*p21)+p31))+(28.0*p01*((45.0*p11)+(9.0*p21)+p31))))+(2.0*p00*((-126.0*p10*p112)+(-189.0*p10*p11*p21)+(-54.0*p11*p20*p21)+(-81.0*p10*p212)+(-45.0*p20*p212)+(9.0*p112*p30)+(-9.0*p212*p30)+(14.0*p012*((45.0*p10)+(9.0*p20)+p30))+(-54.0*p10*p11*p31)+(-30.0*p11*p20*p31)+(-60.0*p10*p21*p31)+(-54.0*p20*p21*p31)+(-6.0*p11*p30*p31)+(-21.0*p21*p30*p31)+(-15.0*p10*p312)+(-21.0*p20*p312)+(-14.0*p30*p312)+(3.0*p01*((2.0*p21*p30)+(7.0*p11*((4.0*p20)+p30))+(-2.0*p20*p31)+(-7.0*p10*((4.0*p21)+p31))))))+(-22.0*((27.0*p112*p20)+(-27.0*p10*p11*p21)+(27.0*p11*p20*p21)+(-27.0*p10*p212)+(18.0*p112*p30)+(45.0*p11*p21*p30)+(45.0*p212*p30)+(5.0*p012*((21.0*p10)+(6.0*p20)+p30))+(-18.0*p10*p11*p31)+(-45.0*p10*p21*p31)+(-45.0*p20*p21*p31)+(30.0*p11*p30*p31)+(105.0*p21*p30*p31)+(-30.0*p10*p312)+(-105.0*p20*p312)+(-280.0*p30*p312)+(p01*((45.0*p10*p11)+(45.0*p11*p20)+(18.0*p20*p21)+(15.0*p11*p30)+(12.0*p21*p30)+(-3.0*p10*p31)+(3.0*p20*p31)+(5.0*p30*p31)))+(p00*((280.0*p012)+(-45.0*p112)+(-45.0*p11*p21)+(-18.0*p212)+(-12.0*p11*p31)+(-15.0*p21*p31)+(-5.0*p312)+(-5.0*p01*((21.0*p11)+(6.0*p21)+p31))))))+(2.0*p01*((45.0*p202*p21)+(54.0*p20*p21*p30)+(21.0*p21*p302)+(3.0*p11*((27.0*p202)+(20.0*p20*p30)+(5.0*p302)))+(9.0*p102*((14.0*p11)+(-1.0*p31)))+(9.0*p202*p31)+(21.0*p20*p30*p31)+(14.0*p302*p31)+(3.0*p10*((9.0*p11*((7.0*p20)+(2.0*p30)))+(2.0*((9.0*p20*p21)+(5.0*p21*p30)+(p30*p31)))))))+(-3.0*((-3.0*p112*((15.0*p202)+(18.0*p20*p30)+(7.0*p302)))+(-2.0*p11*((27.0*p202*p21)+(42.0*p302*(p21+p31))+(7.0*p20*p30*((9.0*p21)+(4.0*p31)))))+(3.0*p102*((15.0*p212)+(18.0*p21*p31)+(7.0*p312)+(2.0*p11*((9.0*p21)+(5.0*p31)))))+(p10*((-6.0*p112*((9.0*p20)+(5.0*p30)))+(28.0*p30*p31*((2.0*p21)+(3.0*p31)))+(-36.0*p11*((p21*p30)+(-1.0*p20*p31)))+(6.0*p20*((9.0*p212)+(21.0*p21*p31)+(14.0*p312)))))+(14.0*((3.0*p202*p31*((2.0*p21)+(3.0*p31)))+(-6.0*p20*p30*(p212+(-5.0*p312)))+(p302*((-9.0*p212)+(-30.0*p21*p31)+(55.0*p312)))))))));
	*accum = ((1.0/4.0)*((-1.0*p012)+p312));
}
static void cube_01(double *pixel, double *accum, Curve<4,2> *curve)
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
	*pixel = ((1.0/18480.0)*((p012*((378.0*p102)+(63.0*p202)+(30.0*p20*p30)+(5.0*p302)+(42.0*p10*((6.0*p20)+p30))))+(-1.0*p002*((-2310.0*p012)+(378.0*p112)+(63.0*p212)+(30.0*p21*p31)+(5.0*p312)+(42.0*p11*((6.0*p21)+p31))+(28.0*p01*((45.0*p11)+(9.0*p21)+p31))))+(-22.0*((45.0*p01*p102)+(45.0*p01*p10*p20)+(27.0*p10*p11*p20)+(18.0*p01*p202)+(27.0*p11*p202)+(-27.0*p102*p21)+(-27.0*p10*p20*p21)+(12.0*p01*p10*p30)+(18.0*p10*p11*p30)+(15.0*p01*p20*p30)+(45.0*p11*p20*p30)+(45.0*p20*p21*p30)+(5.0*p01*p302)+(30.0*p11*p302)+(105.0*p21*p302)+(5.0*p002*((28.0*p01)+(-21.0*p11)+(-6.0*p21)+(-1.0*p31)))+(-18.0*p102*p31)+(-45.0*p10*p20*p31)+(-45.0*p202*p31)+(-30.0*p10*p30*p31)+(-105.0*p20*p30*p31)+(-140.0*p302*p31)+(p00*((-18.0*p20*p21)+(3.0*p11*p30)+(-3.0*p21*p30)+(5.0*p01*((21.0*p10)+(6.0*p20)+p30))+(-12.0*p20*p31)+(-5.0*p30*p31)+(-15.0*p10*((3.0*p11)+(3.0*p21)+p31))))))+(2.0*p00*((-126.0*p10*p112)+(-189.0*p10*p11*p21)+(-54.0*p11*p20*p21)+(-81.0*p10*p212)+(-45.0*p20*p212)+(9.0*p112*p30)+(-9.0*p212*p30)+(14.0*p012*((45.0*p10)+(9.0*p20)+p30))+(-54.0*p10*p11*p31)+(-30.0*p11*p20*p31)+(-60.0*p10*p21*p31)+(-54.0*p20*p21*p31)+(-6.0*p11*p30*p31)+(-21.0*p21*p30*p31)+(-15.0*p10*p312)+(-21.0*p20*p312)+(-14.0*p30*p312)+(3.0*p01*((2.0*p21*p30)+(7.0*p11*((4.0*p20)+p30))+(-2.0*p20*p31)+(-7.0*p10*((4.0*p21)+p31))))))+(2.0*p01*((45.0*p202*p21)+(54.0*p20*p21*p30)+(21.0*p21*p302)+(3.0*p11*((27.0*p202)+(20.0*p20*p30)+(5.0*p302)))+(9.0*p102*((14.0*p11)+(-1.0*p31)))+(9.0*p202*p31)+(21.0*p20*p30*p31)+(14.0*p302*p31)+(3.0*p10*((9.0*p11*((7.0*p20)+(2.0*p30)))+(2.0*((9.0*p20*p21)+(5.0*p21*p30)+(p30*p31)))))))+(-3.0*((-3.0*p112*((15.0*p202)+(18.0*p20*p30)+(7.0*p302)))+(-2.0*p11*((27.0*p202*p21)+(42.0*p302*(p21+p31))+(7.0*p20*p30*((9.0*p21)+(4.0*p31)))))+(3.0*p102*((15.0*p212)+(18.0*p21*p31)+(7.0*p312)+(2.0*p11*((9.0*p21)+(5.0*p31)))))+(p10*((-6.0*p112*((9.0*p20)+(5.0*p30)))+(28.0*p30*p31*((2.0*p21)+(3.0*p31)))+(-36.0*p11*((p21*p30)+(-1.0*p20*p31)))+(6.0*p20*((9.0*p212)+(21.0*p21*p31)+(14.0*p312)))))+(14.0*((3.0*p202*p31*((2.0*p21)+(3.0*p31)))+(-6.0*p20*p30*(p212+(-5.0*p312)))+(p302*((-9.0*p212)+(-30.0*p21*p31)+(55.0*p312)))))))));
	*accum = ((1.0/4.0)*((-2.0*p01)+p012+(-1.0*(-2.0+p31)*p31)));
}
static void cube_11(double *pixel, double *accum, Curve<4,2> *curve)
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
	*pixel = ((1.0/18480.0)*((-1.0*p012*((378.0*p102)+(63.0*p202)+(30.0*p20*p30)+(5.0*p302)+(42.0*p10*((6.0*p20)+p30))))+(p002*((-2310.0*p012)+(378.0*p112)+(63.0*p212)+(30.0*p21*p31)+(5.0*p312)+(42.0*p11*((6.0*p21)+p31))+(28.0*p01*((45.0*p11)+(9.0*p21)+p31))))+(22.0*((45.0*p01*p102)+(45.0*p01*p10*p20)+(27.0*p10*p11*p20)+(18.0*p01*p202)+(27.0*p11*p202)+(-27.0*p102*p21)+(-27.0*p10*p20*p21)+(12.0*p01*p10*p30)+(18.0*p10*p11*p30)+(15.0*p01*p20*p30)+(45.0*p11*p20*p30)+(45.0*p20*p21*p30)+(5.0*p01*p302)+(30.0*p11*p302)+(105.0*p21*p302)+(5.0*p002*((28.0*p01)+(-21.0*p11)+(-6.0*p21)+(-1.0*p31)))+(-18.0*p102*p31)+(-45.0*p10*p20*p31)+(-45.0*p202*p31)+(-30.0*p10*p30*p31)+(-105.0*p20*p30*p31)+(-140.0*p302*p31)+(p00*((-18.0*p20*p21)+(3.0*p11*p30)+(-3.0*p21*p30)+(5.0*p01*((21.0*p10)+(6.0*p20)+p30))+(-12.0*p20*p31)+(-5.0*p30*p31)+(-15.0*p10*((3.0*p11)+(3.0*p21)+p31))))))+(-2.0*p00*((-126.0*p10*p112)+(-189.0*p10*p11*p21)+(-54.0*p11*p20*p21)+(-81.0*p10*p212)+(-45.0*p20*p212)+(9.0*p112*p30)+(-9.0*p212*p30)+(14.0*p012*((45.0*p10)+(9.0*p20)+p30))+(-54.0*p10*p11*p31)+(-30.0*p11*p20*p31)+(-60.0*p10*p21*p31)+(-54.0*p20*p21*p31)+(-6.0*p11*p30*p31)+(-21.0*p21*p30*p31)+(-15.0*p10*p312)+(-21.0*p20*p312)+(-14.0*p30*p312)+(3.0*p01*((2.0*p21*p30)+(7.0*p11*((4.0*p20)+p30))+(-2.0*p20*p31)+(-7.0*p10*((4.0*p21)+p31))))))+(-22.0*((-27.0*p112*p20)+(27.0*p10*p11*p21)+(-27.0*p11*p20*p21)+(27.0*p10*p212)+(-18.0*p112*p30)+(-45.0*p11*p21*p30)+(-45.0*p212*p30)+(-5.0*p012*((21.0*p10)+(6.0*p20)+p30))+(18.0*p10*p11*p31)+(45.0*p10*p21*p31)+(45.0*p20*p21*p31)+(-30.0*p11*p30*p31)+(-105.0*p21*p30*p31)+(30.0*p10*p312)+(105.0*p20*p312)+(280.0*p30*p312)+(42.0*((3.0*p11*p20)+(-3.0*p10*p21)+(3.0*p11*p30)+(6.0*p21*p30)+(p01*((6.0*p10)+(3.0*p20)+p30))+(p00*((10.0*p01)+(-6.0*p11)+(-3.0*p21)+(-1.0*p31)))+(-3.0*p10*p31)+(-6.0*p20*p31)+(-10.0*p30*p31)))+(-1.0*p01*((18.0*p20*p21)+(12.0*p21*p30)+(15.0*p11*((3.0*p20)+p30))+(p10*((45.0*p11)+(-3.0*p31)))+(3.0*p20*p31)+(5.0*p30*p31)))+(-1.0*p00*((280.0*p012)+(-45.0*p112)+(-45.0*p11*p21)+(-18.0*p212)+(-12.0*p11*p31)+(-15.0*p21*p31)+(-5.0*p312)+(-5.0*p01*((21.0*p11)+(6.0*p21)+p31))))))+(-2.0*p01*((45.0*p202*p21)+(54.0*p20*p21*p30)+(21.0*p21*p302)+(3.0*p11*((27.0*p202)+(20.0*p20*p30)+(5.0*p302)))+(9.0*p102*((14.0*p11)+(-1.0*p31)))+(9.0*p202*p31)+(21.0*p20*p30*p31)+(14.0*p302*p31)+(3.0*p10*((9.0*p11*((7.0*p20)+(2.0*p30)))+(2.0*((9.0*p20*p21)+(5.0*p21*p30)+(p30*p31)))))))+(3.0*((-3.0*p112*((15.0*p202)+(18.0*p20*p30)+(7.0*p302)))+(-2.0*p11*((27.0*p202*p21)+(42.0*p302*(p21+p31))+(7.0*p20*p30*((9.0*p21)+(4.0*p31)))))+(3.0*p102*((15.0*p212)+(18.0*p21*p31)+(7.0*p312)+(2.0*p11*((9.0*p21)+(5.0*p31)))))+(p10*((-6.0*p112*((9.0*p20)+(5.0*p30)))+(28.0*p30*p31*((2.0*p21)+(3.0*p31)))+(-36.0*p11*((p21*p30)+(-1.0*p20*p31)))+(6.0*p20*((9.0*p212)+(21.0*p21*p31)+(14.0*p312)))))+(14.0*((3.0*p202*p31*((2.0*p21)+(3.0*p31)))+(-6.0*p20*p30*(p212+(-5.0*p312)))+(p302*((-9.0*p212)+(-30.0*p21*p31)+(55.0*p312)))))))));
	*accum = ((1.0/4.0)*((-2.0*p01)+p012+(-1.0*(-2.0+p31)*p31)));
}


//=============================================================================//
//=========================== Rational Quadratics =============================//
//=============================================================================//

static void rquad(double *pixel, double *accum, Curve<3,3> *curve)
{
}

//=============================================================================//
//============================== Constructors =================================//
//=============================================================================//

template<>
FilterTent<1>::FilterTent()
{
	filter22[filter_elem_flip<W>(0,0)] = line_00;
	filter22[filter_elem_flip<W>(1,0)] = line_10;
	filter22[filter_elem_flip<W>(0,1)] = line_01;
	filter22[filter_elem_flip<W>(1,1)] = line_11;

	filter32[filter_elem_flip<W>(0,0)] = quad_00;
	filter32[filter_elem_flip<W>(1,0)] = quad_10;
	filter32[filter_elem_flip<W>(0,1)] = quad_01;
	filter32[filter_elem_flip<W>(1,1)] = quad_11;

	filter42[filter_elem_flip<W>(0,0)] = cube_00;
	filter42[filter_elem_flip<W>(1,0)] = cube_10;
	filter42[filter_elem_flip<W>(0,1)] = cube_01;
	filter42[filter_elem_flip<W>(1,1)] = cube_11;

	filter33[filter_elem_flip<W>(0,0)] = rquad;
	filter33[filter_elem_flip<W>(1,0)] = rquad;
	filter33[filter_elem_flip<W>(0,1)] = rquad;
	filter33[filter_elem_flip<W>(1,1)] = rquad;
}
