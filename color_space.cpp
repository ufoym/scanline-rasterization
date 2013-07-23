#include <math.h>
#include "color_space.h"

float RGB32F_to_sRGB32F(float c)
{
	if (c < 0)
		return 0;
	else if (c < 0.0031308)
		return 12.92 * c;
	else if (c < 1)
		return 1.055 * pow(c, 0.41666f) - 0.055;
	else
		return 1;
}
float sRGB32F_to_RGB32F(float ic)
{
	if (ic < .04045)
		return ic / (12.92);
	else
		return pow((ic + 0.055)/1.055, 2.4);
}


int RGB32F_to_sRGB8(float c)
{
	if (c < 0)
		return 0;
	else if (c < 0.0031308)
		return 12.92 * 255.0 * c;
	else if (c < 1)
		return 1.055 * 255.0 * pow(c, 0.41666f) - 0.055 * 255.0;
	else
		return 255;
}
float sRGB8_to_RGB32F(int ic)
{
	if (ic < (unsigned char)(.04045*255.0))
		return ic / (12.92*255.0);
	else
		return pow((ic*(1.0/255.0) + 0.055)/1.055, 2.4);
}
