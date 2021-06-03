#pragma once

#include "matrix.h"

class FrequencyResults;

class PotentialSourceCalculator{
public:
	double Oci;
	double Oce;
	point3d B;
	point3d u;
	
	Matrix<complex> calculate_potential_source(
		const FrequencyResults &,
		const Matrix<double> &ksqu,
		const Matrix<complex> &Pik,
		const Matrix<complex> &Pek,
		const Matrix<complex2d> &dndk) const;
};
