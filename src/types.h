#pragma once

#include <type_traits>
#include <cmath>

class complex2d;
class point2d;

class complex{
public:
	double real;
	double imag;

	double abs() const{
		return sqrt(real * real + imag * imag);
	}
	complex &operator*=(double x){
		this->real *= x;
		this->imag *= x;
		return *this;
	}
	complex operator*(double x) const{
		auto ret = *this;
		ret *= x;
		return ret;
	}
	complex2d operator*(const point2d &x) const;
};

static_assert(std::is_pod<complex>::value && sizeof(complex) == sizeof(double) * 2, "oops!");

class complex2d{
public:
	complex x;
	complex y;
};

class point2d{
public:
	double x;
	double y;

	double sqnorm() const{
		return x * x + y * y;
	}
	double norm() const{
		return sqrt(this->sqnorm());
	}
};

class point3d{
public:
	double x;
	double y;
	double z;

	double norm() const{
		return sqrt(x * x + y * y + z * z);
	}
};

complex2d complex::operator*(const point2d &x) const{
	return {
		*this * x.x,
		*this * x.y,
	};
}
