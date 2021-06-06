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
	complex2d &operator*=(const point2d &x);
	complex2d operator*(const point2d &x) const{
		auto ret = *this;
		ret *= x;
		return ret;
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
	point2d &operator*=(double x){
		this->x *= x;
		this->y *= x;
		return *this;
	}
	point2d &operator*=(const point2d &other){
		this->x *= other.x;
		this->y *= other.y;
		return *this;
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

inline complex2d complex::operator*(const point2d &x) const{
	return {
		*this * x.x,
		*this * x.y,
	};
}

inline complex2d &complex2d::operator*=(const point2d &x){
	this->x *= x.x;
	this->y *= x.y;
	return *this;
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

static_assert(std::is_pod<complex>::value && sizeof(complex) == sizeof(double) * 2, "oops!");
static_assert(std::is_pod<complex2d>::value && sizeof(complex2d) == sizeof(double) * 4, "oops!");
static_assert(std::is_pod<point2d>::value && sizeof(point2d) == sizeof(double) * 2, "oops!");
static_assert(std::is_pod<Freq<double>>::value && sizeof(Freq<double>) == sizeof(double) * 6, "oops!");
static_assert(std::is_pod<Freq<complex>>::value && sizeof(Freq<complex>) == sizeof(double) * 12, "oops!");

inline bool nan_check(double x){
	return std::isnan(x);
}

inline bool nan_check(const complex &x){
	return nan_check(x.real) || nan_check(x.imag);
}

inline bool nan_check(const point2d &x){
	return nan_check(x.x) || nan_check(x.y);
}

inline bool nan_check(const point3d &x){
	return nan_check(x.x) || nan_check(x.y) || nan_check(x.z);
}

inline bool nan_check(const complex2d &x){
	return nan_check(x.x) || nan_check(x.y);
}

template <typename T>
bool nan_check(const Freq<T> &x){
	return nan_check(x.ii) || nan_check(x.ie) || nan_check(x.ei) || nan_check(x.en) || nan_check(x.in) || nan_check(x.ee);
}
