#include "constants.h"
#include "timestep.h"

double max(const Matrix<double> &m){
	auto ret = m.get(0, 0);
	m.for_each_unordered([&ret](auto x){
		ret = std::max(ret, x);
	});
	return ret;
}

double max(const point2d &p){
	return std::max(p.x, p.y);
}

double max(const point3d &p){
	return std::max(std::max(p.x, p.y), p.z);
}

double max(const Matrix<point2d> &m){
	auto ret = max(m.get(0, 0));
	m.for_each_unordered([&ret](auto p){
		ret = std::max(ret, max(p));
	});
	return ret;
}

double TimeStepCalculator::calculate_time_step() const{
	auto m = max(*this->vexb);
	m = std::max(m, max(*this->vdmi));
	m = std::max(m, max(*this->vdme));
	m = std::max(m, max(u));
	return std::min(CFL / (m * kmax), dt_max);
}
