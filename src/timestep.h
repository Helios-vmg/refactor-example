#pragma once

#include "matrix.h"

class TimeStepCalculator{
public:
	Matrix<point2d> *vexb = nullptr;
	Matrix<point2d> *vdmi = nullptr;
	Matrix<point2d> *vdme = nullptr;

	double calculate_time_step() const;
};
