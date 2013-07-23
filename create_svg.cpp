#include "rasterizer.h"
#include <fstream>

void create_svg(Array<Curve<2,2> > &lines, Array<Curve<3,2> > &quads, Array<Curve<4,2> > &cubics, Array<Curve<3,3> > &rquads)
{
	ifstream f("data/floral2.txt");
	if (f.fail())
	{
		printf("failed to open input file\n");
		exit(0);
	}

	for (;;)
	{
		int n;
		f >> n;

		if (f.eof())
			break;

		if (n == 2)
		{
			Curve<2,2> c;
			for (int i = 0; i < n; i++)
				f >> c.p[i][0] >> c.p[i][1];

			if (f.eof())
				break;

			lines.push_back(c);
		}
		else if (n == 3)
		{
			Curve<3,2> c;
			for (int i = 0; i < n; i++)
				f >> c.p[i][0] >> c.p[i][1];

			if (f.eof())
				break;

			quads.push_back(c);
		}
		else if (n == 4)
		{
			Curve<4,2> c;
			for (int i = 0; i < n; i++)
				f >> c.p[i][0] >> c.p[i][1];

			if (f.eof())
				break;

			cubics.push_back(c);
		}
	}

	f.close();
}
