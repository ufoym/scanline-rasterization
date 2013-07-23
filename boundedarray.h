#pragma once

#include "assert.h"

template <class T, int capacity = 9>
struct BoundedArray
{
	BoundedArray(){init();}

	int s; // size
	T data[capacity];

	void init(){s = 0;}
	void push_back(T &e){data[s++] = e;	assert(s <= capacity);}
	T& push_back(){assert(s+1 <= capacity); return data[s++];}
	void pop_back(){s--; assert(s >= 0);}
	T &operator[](int i){return data[i]; assert(i >= 0 && i < s);}
	int size(){return s;}
	void resize(int s){this->s = s; assert(s <= capacity);}
	void clear(){s = 0;}
	void reserve(int s){}
	
	T* begin(){return data;}
	T* end(){return data + s;}

	T &first()
	{
		return data[0];
	}

	T &last()
	{
		return data[s-1];
	}

	// math operators
	void operator+=(const BoundedArray<T, capacity> &a)
	{
		assert(a.s == s);
		for (int i = 0; i < s; i++)
			data[i] += a.data[i];
	}

	void operator-=(const BoundedArray<T, capacity> &a)
	{
		assert(a.s == s);
		for (int i = 0; i < s; i++)
			data[i] -= a.data[i];
	}

	void operator=(const T &a)
	{
		for (int i = 0; i < s; i++)
			data[i] = a;
	}

};
