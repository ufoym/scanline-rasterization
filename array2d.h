#pragma once
#include <assert.h>

template <class T>
struct Array2D
{
	Array2D(){data = 0; clear();}
	Array2D(const Array2D<T> &a){data = 0; copy(a);}
	~Array2D(){clear();}

	T *data;
	int data_size;
	int s[2]; // size

	void operator=(const Array2D<T> &a) {copy(a);}
	void copy(const Array2D<T> &a)
	{
		resize(a.s);
		for (int i = 0; i < data_size; i++)
			data[i] = a.data[i];
	}

	void resize(int w, int h)
	{
		clear();

		data_size = w*h;
		s[0] = w;
		s[1] = h;

		data = new T[data_size];
	}
	template <class S> void resize(S p) {resize(p[0], p[1]);}

	T &operator()(int x, int y)
	{
		return data[y*s[0] + x];
	}

	T &get(int x, int y)
	{
		return data[y*s[0] + x];
	}

	T &get_clamp(int x, int y)
	{
		if (x < 0)
			x = 0;
		else if (x > s[0]-1)
			x = s[0]-1;

		if (y < 0)
			y = 0;
		else if (y > s[1]-1)
			y = s[1]-1;

		return data[y*s[0] + x];
	}

	T &get_wrap(int x, int y)
	{
		x = x % s[0];
		y = y % s[1];
		
		if (x < 0)
			x += s[0];
		if (y < 0)
			y += s[1];

		return data[y*s[0] + x];
	}

	template <class S>
	T &operator[](S p)
	{
		return data[p[1]*s[0] + p[0]];
	}

	void clear()
	{
		delete[] data;
		data = 0;
		s[0] = s[1] = 0;
		data_size = 0;
	}

	// math operators
	void operator+=(const Array2D<T> &a)
	{
		assert(a.s[0] == s[0] && a.s[1] == s[1]);
		for (int i = 0; i < data_size; i++)
			data[i] += a.data[i];
	}

	void operator-=(const Array2D<T> &a)
	{
		assert(a.s[0] == s[0] && a.s[1] == s[1]);
		for (int i = 0; i < data_size; i++)
			data[i] -= a.data[i];
	}

	void operator=(const T &a)
	{
		for (int i = 0; i < data_size; i++)
			data[i] = a;
	}
};
