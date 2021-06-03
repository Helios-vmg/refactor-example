#include "matrix.h"

Matrix<complex> to_fourier(const Matrix<double> &real){
	Matrix<complex> ret(real.cols(), real.rows() / 2 + 1);
	auto r2c = fftw_plan_dft_r2c_2d(real.cols(), real.rows(), (double *)&real.get(0, 0), (fftw_complex *)&ret.get(0, 0), FFTW_ESTIMATE);
	fftw_execute(r2c);
	fftw_destroy_plan(r2c);
	fftw_cleanup();
	return ret;
}

Matrix<double> from_fourier(const Matrix<complex> &c){
	auto copy = c;
	Matrix<double> ret(c.cols(), (c.rows() - 1) * 2);
	auto c2r = fftw_plan_dft_c2r_2d(copy.cols(), copy.rows(), (fftw_complex *)&copy.get(0, 0), &ret.get(0, 0), FFTW_ESTIMATE);
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

template <typename T>
auto basic_derivk(const Matrix<complex> &a, const T &b){
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

Matrix<complex> laplaciank(const Matrix<complex> &a, const Matrix<double> &b){
	auto ret = a * b;
	ret.for_each_unordered([](auto &p){
		p *= -1;
	});
	return ret;
}

Matrix<complex> convolve2d(const Matrix<complex> &a, const Matrix<complex> &b){
	auto temp = from_fourier(a);
	temp.elementwise_multiplication(from_fourier(b));
	return to_fourier(temp);
}
