#pragma once

#include "matrix.h"

class FrequencyResults;

class FourierMesh2d{
public:
	Matrix<point2d> k;
	Matrix<double> ksqu;
	Matrix<double> ninvksqu;

	FourierMesh2d() = default;
	FourierMesh2d(double Lx, double Ly, size_t w, size_t h, int FourierMeshType);
	FourierMesh2d(const FourierMesh2d &) = delete;
	FourierMesh2d &operator=(const FourierMesh2d &) = delete;
	FourierMesh2d(FourierMesh2d &&other){
		*this = std::move(other);
	}
	FourierMesh2d &operator=(FourierMesh2d &&other){
		this->k = std::move(other.k);
		this->ksqu = std::move(other.ksqu);
		this->ninvksqu = std::move(other.ninvksqu);
		return *this;
	}
};

class PotentialSourceCalculator{
public:
	double Oci;
	double Oce;
	point3d B;
	point3d u;
	const Matrix<complex2d> *dndk;
	
	Matrix<complex> calculate_potential_source(
		const FrequencyResults &,
		const Matrix<double> &ksqu,
		const Matrix<complex> &Pik,
		const Matrix<complex> &Pek) const;
	int potentialk(const FrequencyResults &, const FourierMesh2d &, Matrix<complex> &phik, const Matrix<complex> &potential_source, double err_max, int max_iter);
};
