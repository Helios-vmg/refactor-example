#pragma once

#include <type_traits>
#include <cmath>

class complex2d;
class point2d;

class complex{
public:
	double real;
	double imag;

	double sqabs() const{
		return real * real + imag * imag;
	}
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
	complex &operator/=(double x){
		this->real /= x;
		this->imag /= x;
		return *this;
	}
	complex operator/(double x) const{
		auto ret = *this;
		ret /= x;
		return ret;
	}
	complex &operator+=(const complex &c){
		this->real += c.real;
		this->imag += c.imag;
		return *this;
	}
	complex operator+(const complex &c) const{
		auto ret = *this;
		ret += c;
		return ret;
	}
	complex &operator-=(const complex &c){
		this->real -= c.real;
		this->imag -= c.imag;
		return *this;
	}
	complex operator-(const complex &c) const{
		auto ret = *this;
		ret -= c;
		return ret;
	}
	complex operator-() const{
		return { -real, -imag };
	}
	complex2d operator*(const point2d &x) const;
};

static_assert(std::is_pod<complex>::value && sizeof(complex) == sizeof(double) * 2, "oops!");

class complex2d{
public:
	complex x;
	complex y;
	complex2d &operator+=(const complex2d &c){
		this->x += c.x;
		this->y += c.y;
		return *this;
	}
	complex2d operator+(const complex2d &c) const{
		auto ret = *this;
		ret += c;
		return ret;
	}
	complex2d &operator-=(const complex2d &c){
		this->x -= c.x;
		this->y -= c.y;
		return *this;
	}
	complex2d operator-(const complex2d &c) const{
		auto ret = *this;
		ret -= c;
		return ret;
	}
	complex2d operator-() const{
		return { -x, -y };
	}
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

template <typename T>
class Freq{
public:
	//ion-ion
	T ii;
	//ion-electron
	T ie;
	//electron-ion
	T ei;
	//electron-neutral (did you mean "neutron"?)
	T en;
	//ion-neutral (did you mean "neutron"?)
	T in;
	//electron-electron
	T ee;
};
