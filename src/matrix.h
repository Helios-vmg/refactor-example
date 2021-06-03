#pragma once

#include "fftw3.h"
#include "types.h"
#include <type_traits>
#include <memory>
#include <algorithm>
#include <stdexcept>

#define CHECK_ACCESSES

template <typename T>
class Matrix{
	static_assert(std::is_pod<T>::value, "T must be a POD");

	struct Releaser{
		void operator()(T *p){
			if (p)
				fftw_free(p);
		}
	};

	std::unique_ptr<T, Releaser> matrix;
	size_t w = 0;
	size_t h = 0;

	static std::unique_ptr<T, Releaser> allocate(size_t w, size_t h){
		std::unique_ptr<T, Releaser> ret((T *)fftw_malloc(w * h * sizeof(T)));
		if (!ret)
			throw std::bad_alloc();
		return ret;
	}
public:
	Matrix() = default;
	Matrix(size_t w, size_t h, const T &init = {}): w(w), h(h){
		this->matrix = allocate(w, h);
		auto p = this->matrix.get();
		std::fill(p, p + w * h, init);
	}
	Matrix(const Matrix &other){
		if (!other)
			return;
		this->w = other.w;
		this->h = other.h;
		this->matrix = allocate(this->w, this->h);
		memcpy(this->matrix.get(), other.matrix.get(), this->w * this->h * sizeof(T));
	}
	Matrix &operator=(const Matrix &other){
		if (&other == this)
			return *this;
		if (!other){
			this->matrix.reset();
			this->w = this->h = 0;
			return *this;
		}
		if (other.w != this->w || other.h != this->h){
			this->w = other.w;
			this->h = other.h;
			this->matrix = allocate(this->w, this->h);
		}
		memcpy(this->matrix.get(), other.matrix.get(), this->w * this->h * sizeof(T));
		return *this;
	}
	Matrix(Matrix &&other){
		*this = std::move(other);
	}
	Matrix &operator=(Matrix &&other){
		this->matrix = std::move(other.matrix);
		this->w = other.w;
		this->h = other.h;
		other.w = other.h = 0;
		return *this;
	}
	bool operator!() const{
		return !this->matrix;
	}
	T &get(size_t row, size_t col){
#ifdef CHECK_ACCESSES
		if (col >= this->w || row >= this->h)
			throw std::runtime_error("bad access");
#endif
		return this->matrix.get()[col * this->h + row];
	}
	const T &get(size_t row, size_t col) const{
#ifdef CHECK_ACCESSES
		if (col >= this->w || row >= this->h)
			throw std::runtime_error("bad access");
#endif
		return this->matrix.get()[col * this->h + row];
	}
	size_t cols() const{
		return this->w;
	}
	size_t rows() const{
		return this->h;
	}
	template <typename F>
	void for_each(const F &f){
		for (size_t i = 0; i < w; i++)
			for (size_t j = 0; j < h; j++)
				f(j, i, this->get(j, i));
	}
	template <typename F>
	void for_each_unordered(const F &f){
		auto p = this->matrix.get();
		for (size_t i = w * h; i--;)
			f(p[i]);
	}
	template <typename T2>
	Matrix<T> &operator*=(const T2 &x){
		this->for_each_unordered([x](auto &p){
			p *= x;
		});
		return *this;
	}
	template <typename T2>
	Matrix<T> operator*(const T2 &x) const{
		auto ret = *this;
		ret *= x;
		return ret;
	}
	void elementwise_multiplication(const Matrix &other){
		this->for_each([other](auto j, auto i, auto &p){
			p *= other.get(j, i);
		});
	}
};

Matrix<complex> to_fourier(const Matrix<double> &);
Matrix<Freq<complex>> to_fourier(const Matrix<Freq<double>> &);
Matrix<double> from_fourier(const Matrix<complex> &);
Matrix<complex> derivk(const Matrix<complex> &, const Matrix<double> &);
Matrix<complex2d> derivk(const Matrix<complex> &, const Matrix<point2d> &);
Matrix<complex> laplaciank(const Matrix<complex> &, const Matrix<double> &);
Matrix<complex> convolve2d(const Matrix<complex> &, const Matrix<complex> &);
