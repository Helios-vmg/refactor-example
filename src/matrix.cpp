#include "matrix.h"

std::ostream &operator<<(std::ostream &stream, const complex &c){
	return stream << "(" << c.real << ", " << c.imag << ")";
}

std::ostream &operator<<(std::ostream &stream, const point2d &p){
	return stream << "(" << p.x << ", " << p.y << ")";
}

std::ostream &operator<<(std::ostream &stream, const complex2d &p){
	return stream << "(" << p.x << ", " << p.y << ")";
}

Matrix<complex> to_fourier(const Matrix<double> &real){
	real.nan_check();
	Matrix<complex> ret(real.cols(), real.rows() / 2 + 1);
	auto copy = real;
	auto r2c = fftw_plan_dft_r2c_2d(real.cols(), real.rows(), &copy.get(0, 0), (fftw_complex *)&ret.get(0, 0), FFTW_ESTIMATE);
	fftw_execute(r2c);
	fftw_destroy_plan(r2c);
	fftw_cleanup();
	return ret;
}

Matrix<complex2d> to_fourier(const Matrix<point2d> &real){
	real.nan_check();
	Matrix<complex2d> ret(real.cols(), real.rows() / 2 + 1);
	int dims[2] = {real.cols(), real.rows()};
	auto r2c = fftw_plan_many_dft_r2c(2, dims, 2, (double *)&real.get(0, 0), nullptr, 2, 1, (fftw_complex *)&ret.get(0, 0), nullptr, 2, 1, FFTW_ESTIMATE);
	fftw_execute(r2c);
	fftw_destroy_plan(r2c);
	fftw_cleanup();
	return ret;
}

Matrix<Freq<complex>> to_fourier(const Matrix<Freq<double>> &real){
	real.nan_check();
	Matrix<Freq<complex>> ret(real.cols(), real.rows() / 2 + 1);
	int dims[2] = {real.cols(), real.rows()};
	auto r2c = fftw_plan_many_dft_r2c(2, dims, 6, (double *)&real.get(0, 0), nullptr, 6, 1, (fftw_complex *)&ret.get(0, 0), nullptr, 6, 1, FFTW_ESTIMATE);
	fftw_execute(r2c);
	fftw_destroy_plan(r2c);
	fftw_cleanup();
	return ret;
}

Matrix<double> from_fourier(const Matrix<complex> &c){
	c.nan_check();
	auto copy = c;
	Matrix<double> ret(c.cols(), (c.rows() - 1) * 2);
	auto c2r = fftw_plan_dft_c2r_2d(copy.cols(), copy.rows(), (fftw_complex *)&copy.get(0, 0), &ret.get(0, 0), FFTW_ESTIMATE);
	fftw_execute(c2r);
	fftw_destroy_plan(c2r);
	fftw_cleanup();
	ret *= 1.0 / (ret.cols() * ret.rows());
	return ret;
}

Matrix<point2d> from_fourier(const Matrix<complex2d> &c){
	c.nan_check();
	Matrix<point2d> ret(c.cols(), (c.rows() - 1) * 2);
	int dims[2] = {c.cols(), c.rows()};
	auto c2r = fftw_plan_many_dft_c2r(2, dims, 2, (fftw_complex *)&c.get(0, 0), dims, 2, 1, &ret.get(0, 0).x, dims, 2, 1, FFTW_ESTIMATE);
	fftw_execute(c2r);
	fftw_destroy_plan(c2r);
	fftw_cleanup();
	ret *= 1.0 / (ret.cols() * ret.rows());
	return ret;
}

Matrix<complex> operator*(const Matrix<complex> &a, const Matrix<double> &b){
	auto ret = a;
	ret.for_each([&b](auto j, auto i, auto &p){
		p *= b.get(j, i);
	});
	return ret;
}

Matrix<complex2d> operator*(const Matrix<complex> &a, const Matrix<point2d> &b){
	Matrix<complex2d> ret(a.cols(), a.rows());
	ret.for_each([&a, &b](auto j, auto i, complex2d &p){
		p = a.get(j, i) * b.get(j, i);
	});
	return ret;
}

Matrix<complex2d> operator*(const Matrix<complex2d> &a, const Matrix<point2d> &b){
	Matrix<complex2d> ret = a;
	ret.for_each([&b](auto j, auto i, complex2d &p){
		p *= b.get(j, i);
	});
	return ret;
}

void square(complex &c){
	auto temp = c.real;
	c.real = -c.imag;
	c.imag = temp;
}

void square(Matrix<complex> &m){
	m.for_each([](auto j, auto i, auto &p){
		square(p);
	});
}

void square(Matrix<complex2d> &m){
	m.for_each([](auto j, auto i, auto &p){
		square(p.x);
		square(p.y);
	});
}

template <typename T, typename T2>
auto basic_derivk(const Matrix<T> &a, const T2 &b){
	a.nan_check();
	b.nan_check();
	auto ret = a * b;
	square(ret);
	return ret;
}

Matrix<complex> derivk(const Matrix<complex> &a, const Matrix<double> &b){
	return basic_derivk(a, b);
}

Matrix<complex2d> derivk(const Matrix<complex> &a, const Matrix<point2d> &b){
	return basic_derivk(a, b);
}

Matrix<complex2d> derivk(const Matrix<complex2d> &a, const Matrix<point2d> &b){
	return basic_derivk(a, b);
}

Matrix<complex> laplaciank(const Matrix<complex> &a, const Matrix<double> &b){
	a.nan_check();
	b.nan_check();
	auto ret = a * b;
	ret.for_each_unordered([](auto &p){
		p *= -1;
	});
	return ret;
}

Matrix<complex2d> fourier_division(const Matrix<complex2d> &a, const Matrix<complex> &b){
	a.nan_check();
	b.nan_check();
	auto a2 = from_fourier(a);
	auto b2 = from_fourier(b);
	a2.for_each_with(b2, [](auto j, auto i, point2d &p, double x){
		p *= 1.0 / x;
	});
	return to_fourier(a2);
}
