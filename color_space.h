
//** sRGB <-> linear RGB
// all components are treated the same, so only scalar functions are needed
// alpha in a texture should not have this conversion applied.

float RGB32F_to_sRGB32F(float c);
float sRGB32F_to_RGB32F(float ic);

int RGB32F_to_sRGB8(float c);
float sRGB8_to_RGB32F(int ic);
