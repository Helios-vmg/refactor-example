#pragma once

#include <type_traits>
#include <cmath>

class complex{
public:
	double real;
	double imag;

	double abs() const{
		return sqrt(real * real + imag * imag);
	}
};

static_assert(std::is_pod<complex>::value && sizeof(complex) == sizeof(double) * 2, "oops!");

class point3d{
public:
	double x;
	double y;
	double z;

	double norm() const{
		return sqrt(x * x + y * y + z * z);
	}
};
