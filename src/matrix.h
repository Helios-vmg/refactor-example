#pragma once

#include "fftw3.h"
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
		std::fill(p[0], p[w * h], init);
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
	T &get(size_t x, size_t y){
#ifdef CHECK_ACCESSES
		if (x >= this->w || y >= this->h)
			throw std::runtime_error("bad access");
#endif
		return this->matrix.get()[x + y * this->w];
	}
	const T &get(size_t x, size_t y) const{
#ifdef CHECK_ACCESSES
		if (x >= this->w || y >= this->h)
			throw std::runtime_error("bad access");
#endif
		return this->matrix.get()[x + y * this->w];
	}
};
