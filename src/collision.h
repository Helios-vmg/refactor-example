#pragma once

#include "matrix.h"

class FrequencyResults{
public:
	Matrix<Freq<complex>> nuk;
	Matrix<complex> isigPk;
	Matrix<complex> invnk;

	FrequencyResults() = default;
	FrequencyResults(const FrequencyResults &) = default;
	FrequencyResults &operator=(const FrequencyResults &) = default;
	FrequencyResults(FrequencyResults &&other){
		*this = std::move(other);
	}
	FrequencyResults &operator=(FrequencyResults &&other){
		this->nuk = std::move(other.nuk);
		this->isigPk = std::move(other.isigPk);
		this->invnk = std::move(other.invnk);
		return *this;
	}
};

class CollisionFreqCalculator{
public:
	FrequencyResults calculate_collision_frequencies(const Matrix<complex> &nek, const Matrix<complex> &Tik, const Matrix<complex> &Tek) const;
};
