#include "constants.h"
#include "matrix.h"
#include "types.h"
#include "collision.h"
#include "potential.h"
#include "timestep.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cstdio>

// Make x coordinates. We want them to change as we change the zeroth index and be constant as we change the first index.
Matrix<point2d> make_spatial_mesh(size_t cols, size_t rows){
	Matrix<point2d> ret(cols, rows);
	ret.for_each([](auto j, auto i, auto &p){
		p.x = i * dx;
		p.y = j * dy;
	});
	return ret;
}

void print_matrix(const char *path, const Matrix<double> &m){
	std::ofstream stream(path);
	if (!stream)
		throw std::runtime_error("failed to open file");
	auto cs = m.cols();
	auto rs = m.rows();
	for (size_t i = 0; i < cs; i++){
		for (size_t j = 0; j < rs; j++){
			char temp[128];
			sprintf(temp, "%+4.2le  ", m.get(j, i));
			stream << temp;
		}
		stream << std::endl;
	}
}

void print_matrix(const char *path_x, const char *path_y, const Matrix<point2d> &m){
	std::ofstream stream_x(path_x);
	std::ofstream stream_y(path_y);
	if (!stream_x || !stream_y)
		throw std::runtime_error("failed to open file");
	auto cs = m.cols();
	auto rs = m.rows();
	for (size_t i = 0; i < cs; i++){
		for (size_t j = 0; j < rs; j++){
			auto &p = m.get(j, i);
			char temp[128];
			sprintf(temp, "%+4.2le  ", p.x);
			stream_x << temp;
			sprintf(temp, "%+4.2le  ", p.y);
			stream_y << temp;
		}
		stream_x << std::endl;
		stream_y << std::endl;
	}
}

Matrix<complex2d> calcV_ExBk(const Matrix<complex2d> &dphidk){
	auto ret = dphidk;
	ret.for_each_unordered([k = B.z / B2](auto &p){
		auto copy = p;
		p.x = copy.y * -k;
		p.y = copy.x * k;
	});
	return ret;
}

const Matrix<complex2d> calc_diamag(const Matrix<complex2d> &dpdk, double qa, const Matrix<complex> &nek){
	auto prediv = dpdk;
	prediv.for_each_unordered([k = B.z / (B2 * qa)](auto &p){
		auto copy = p;
		p.x = copy.y * -k;
		p.y = copy.x * k;
	});
	return fourier_division(prediv, nek);
}

Matrix<complex> calc_residualn(const Matrix<complex2d> &vexbk, const Matrix<point2d> &k, const Matrix<complex> &nek){
	auto dnekk = derivk(nek, k);
	auto mult = convolve2d(vexbk, dnekk);
	Matrix<complex> ret(mult.geom());
	ret.for_each([&mult](auto j, auto i, auto &p){
		auto s = mult.get(j, i);
		p = s.x + s.y;
	});
	return ret;
}

Matrix<complex> calc_residualt(const Matrix<complex2d> &vok, const Matrix<complex> &tempink, const Matrix<point2d> &k){
	auto dtempink = derivk(tempink, k);
	auto dvok = derivk(vok, k);
	Matrix<complex> divvo(dvok.cols(), dvok.rows());
	divvo.for_each([&dvok](auto j, auto i, auto &p){
		auto s = dvok.get(j, i);
		p = s.x + s.y;
	});
	auto mult = convolve2d(dtempink, vok);
	auto ret = convolve2d(tempink, divvo) * (2.0 / 3.0);
	ret.for_each([&mult](auto j, auto i, auto &p){
		auto s = mult.get(j, i);
		p += s.x + s.y;
	});
	return ret;
}

Matrix<complex> calc_sourcen(const Matrix<double> &ksqu, const Matrix<complex> &nek, double d){
	return laplaciank(nek, ksqu) * d;
}

Matrix<complex> RK4_1(double dt, const Matrix<complex> &residual, const Matrix<complex> &source, int stage){
	auto alpha = dt / (4 - stage);
	auto ret = residual;
	ret.for_each([alpha, &residual, &source](auto j, auto i, auto &p){
		p = -(residual.get(j, i) - source.get(j, i)) * alpha;
	});
	return ret;
}

Matrix<complex> RK4(const Matrix<complex> &f, double dt, const Matrix<complex> &residual, const Matrix<complex> &source, int stage){
	auto alpha = dt / (4 - stage);
	auto ret = f;
	ret.for_each([alpha, &residual, &source](auto j, auto i, auto &p){
		p -= (residual.get(j, i) - source.get(j, i)) * alpha;
	});
	return ret;
}

std::string get_filename(const char *s, int n, int w = 8){
	std::stringstream stream;
	stream << s << std::setw(w) << std::setfill('0') << n << ".gkyl";
	return stream.str();
}

void print_binary_matrix(const std::string &filename, const Matrix<double> &m){
	std::ofstream file(filename, std::ios::binary);
	if (!file)
		throw std::runtime_error("failed to open file");
	
	const std::uint64_t real_type = 2;
    const std::uint64_t dimension_count = 2;
    std::uint64_t dimensions[] = { m.cols(), m.rows() };
	const double lower[] = { 0, 0 };
	const double upper[] = { 1, 1 };
	std::uint64_t double_size = sizeof(double);
	std::uint64_t cell_count = dimensions[0] * dimensions[1];

	file.write((const char *)&real_type, sizeof(real_type));
	file.write((const char *)&dimension_count, sizeof(dimension_count));
	file.write((const char *)&dimensions, sizeof(dimensions));
	file.write((const char *)&lower, sizeof(lower));
	file.write((const char *)&upper, sizeof(upper));
	file.write((const char *)&double_size, sizeof(double_size));
	file.write((const char *)&cell_count, sizeof(cell_count));
	file.write((const char *)&m.get(0, 0), double_size * cell_count);
}

void print_binary_matrix(const std::string &filename, const Matrix<complex> &m){
	std::ofstream file(filename, std::ios::binary);
	if (!file)
		throw std::runtime_error("failed to open file");
	
	const std::uint64_t real_type = 2;
    const std::uint64_t dimension_count = 2;
    std::uint64_t dimensions[] = { m.cols(), m.rows() };
	const double lower[] = { 0, 0 };
	const double upper[] = { 1, 1 };
	std::uint64_t double_size = sizeof(complex);
	std::uint64_t cell_count = dimensions[0] * dimensions[1];

	file.write((const char *)&real_type, sizeof(real_type));
	file.write((const char *)&dimension_count, sizeof(dimension_count));
	file.write((const char *)&dimensions, sizeof(dimensions));
	file.write((const char *)&lower, sizeof(lower));
	file.write((const char *)&upper, sizeof(upper));
	file.write((const char *)&double_size, sizeof(double_size));
	file.write((const char *)&cell_count, sizeof(cell_count));
	file.write((const char *)&m.get(0, 0), double_size * cell_count);
}

void print_binary_matrix(const std::string &filename, const Matrix<complex2d> &m){
	std::ofstream file(filename, std::ios::binary);
	if (!file)
		throw std::runtime_error("failed to open file");
	
	const std::uint64_t real_type = 2;
    const std::uint64_t dimension_count = 2;
    std::uint64_t dimensions[] = { m.cols(), m.rows() };
	const double lower[] = { 0, 0 };
	const double upper[] = { 1, 1 };
	std::uint64_t double_size = sizeof(complex2d);
	std::uint64_t cell_count = dimensions[0] * dimensions[1];

	file.write((const char *)&real_type, sizeof(real_type));
	file.write((const char *)&dimension_count, sizeof(dimension_count));
	file.write((const char *)&dimensions, sizeof(dimensions));
	file.write((const char *)&lower, sizeof(lower));
	file.write((const char *)&upper, sizeof(upper));
	file.write((const char *)&double_size, sizeof(double_size));
	file.write((const char *)&cell_count, sizeof(cell_count));
	file.write((const char *)&m.get(0, 0), double_size * cell_count);
}

void print_binary_matrix(const std::string &filename, const Matrix<point2d> &m){
	std::ofstream file(filename, std::ios::binary);
	if (!file)
		throw std::runtime_error("failed to open file");
	
	const std::uint64_t real_type = 2;
    const std::uint64_t dimension_count = 2;
    std::uint64_t dimensions[] = { m.cols(), m.rows() };
	const double lower[] = { 0, 0 };
	const double upper[] = { 1, 1 };
	std::uint64_t double_size = sizeof(point2d);
	std::uint64_t cell_count = dimensions[0] * dimensions[1];

	file.write((const char *)&real_type, sizeof(real_type));
	file.write((const char *)&dimension_count, sizeof(dimension_count));
	file.write((const char *)&dimensions, sizeof(dimensions));
	file.write((const char *)&lower, sizeof(lower));
	file.write((const char *)&upper, sizeof(upper));
	file.write((const char *)&double_size, sizeof(double_size));
	file.write((const char *)&cell_count, sizeof(cell_count));
	file.write((const char *)&m.get(0, 0), double_size * cell_count);
}

int main(){
	Matrix<double> Ti(nx, ny, 1000);
	Matrix<double> Te(nx, ny, 1000);
	Matrix<double> ne(nx, ny);
	auto Pi = ne;
	auto Pe = ne;
	auto phi = ne;
	auto tem = ne;
	auto tem2 = ne;
	Matrix<complex> base(nx, nyk);

	FourierMesh2d fourier_mesh(Lx, Ly, nx, nyk, 1);

	{
		static const char * const XXgrid = "X.txt";
		static const char * const YYgrid = "Y.txt";
		auto mesh2d = make_spatial_mesh(nx, ny);
		mesh2d.for_each([&ne, &Ti, &Te, &Pi, &Pe](auto j, auto i, auto p){
			auto &nep = ne.get(j, i);

			auto t = tanh(b * (p.x + c));
			auto bg1 = -bg * (p.x - xg) * (p.x - xg);
			auto bg2 = -bg * (p.x - Lx + xg) * (p.x - Lx + xg);
			auto expsum = exp(bg1) + exp(bg2);

			nep = a * t;
			nep += d;
			nep += a2 * t;
			nep += d2;
			nep += .02 * cos(2 * tau * p.y / Ly) * expsum;
			nep *= 1.E11;

			Pe.get(j, i) = Pi.get(j, i) = nep * (1000 * 1.38E-23);
		});
		print_matrix(XXgrid, YYgrid, mesh2d);
	}
	
	auto nek = to_fourier(ne);
	auto Tik = to_fourier(Ti);
	auto Tek = to_fourier(Te);
	auto phik = to_fourier(phi);
	auto Pik = to_fourier(Pi);
	auto Pek = to_fourier(Pe);

	auto dndk = derivk(nek, fourier_mesh.k);
	auto dphidk = derivk(phik, fourier_mesh.k);
	auto dpedk = derivk(Pek, fourier_mesh.k);
	auto dpidk = derivk(Pik, fourier_mesh.k);

	// Calculate ion and electron gyrofrequencies - qb/m
	
	CollisionFreqCalculator cfc;
	
	auto collision_frequencies = cfc.calculate_collision_frequencies(nek, Tik, Tek);

	PotentialSourceCalculator psc;
	psc.dndk = &dndk;

	auto potSourcek = psc.calculate_potential_source(collision_frequencies, fourier_mesh.ksqu, Pik, Pek);

	print_binary_matrix(get_filename("potSourcek", 0, 5), potSourcek);
	
	int phi_iter = psc.potentialk(collision_frequencies, fourier_mesh, phik, potSourcek, err_max, phi_iter_max);
	
	// Initialize a time vector of length iter_max+1.
	std::vector<double> time;
	time.reserve(iter_max);
	time.push_back(0);

	ne = from_fourier(nek);
	Te = from_fourier(Tek);
	Ti = from_fourier(Tik);
	phi = from_fourier(phik);
	
	static const char * const neinitial = "ne_initial.txt";
	static const char * const Teinitial = "Te_initial.txt";
	static const char * const Tiinitial = "Ti_initial.txt";
	static const char * const phiinitial = "phi_initial.txt";

	print_matrix(neinitial, ne);
	print_matrix(Teinitial, Te);
	print_matrix(Tiinitial, Ti);
	print_matrix(phiinitial, phi);

	print_binary_matrix(get_filename("phik", 0, 5), phik);
	print_binary_matrix(get_filename("dphidk", 0, 5), dphidk);
	
	auto vexbk = calcV_ExBk(dphidk);
	auto vexb = from_fourier(vexbk);
	auto vdmek = calc_diamag(dpedk, -e, nek);
	auto vdmik = calc_diamag(dpidk, e, nek);

	auto vdme = from_fourier(vdmek);
	auto vdmi = from_fourier(vdmik);

	TimeStepCalculator tsc;
	tsc.vexb = &vexb;
	tsc.vdmi = &vdmi;
	tsc.vdme = &vdme;
	
	int saveNum = 1;
	for (int iter = 0; iter < iter_max; iter++){
		// Calculate pressures
		Pek = convolve2d(nek, Tek) * kb;
		Pik = convolve2d(nek, Tik) * kb;

		dndk = derivk(nek, fourier_mesh.k);
		dphidk = derivk(phik, fourier_mesh.k);
		dpedk = derivk(Pek, fourier_mesh.k);
		dpidk = derivk(Pik, fourier_mesh.k);

		print_binary_matrix(get_filename("dphidk", iter + 1, 5), dphidk);
		
		collision_frequencies = cfc.calculate_collision_frequencies(nek, Tik, Tek);
		potSourcek = psc.calculate_potential_source(collision_frequencies, fourier_mesh.ksqu, Pik, Pek);
		
		phi_iter = psc.potentialk(collision_frequencies, fourier_mesh, phik, potSourcek, err_max, phi_iter_max);
		// Calculate all  velocities
		vexbk = calcV_ExBk(dphidk);
		vexb = from_fourier(vexbk);
		vdmek = calc_diamag(dpedk, -e, nek);
		vdmik = calc_diamag(dpidk, e, nek);
		vdme = from_fourier(vdmek);
		vdmi = from_fourier(vdmik);

		// Calculate time step
		double dt = tsc.calculate_time_step();
		// Update time
		time.push_back(time.back() + dt);

		auto nek_old = nek;
		auto Tek_old = Tek;
		auto Tik_old = Tik;
		
		// Begin RK method
	
		for (int stage = 0; stage < 4; stage++){
			// Calculate diamagnetic drifts
			Pek = convolve2d(nek, Tek) * kb;
			Pik = convolve2d(nek, Tik) * kb;

			dpedk = derivk(Pek, fourier_mesh.k);
			dpidk = derivk(Pik, fourier_mesh.k);

			vdmek = calc_diamag(dpedk, -e, nek);
			vdmik = calc_diamag(dpidk, e, nek);

			// Get total velocity
			auto veok = vdmek + vexbk;
			auto viok = vdmik + vexbk;
			
			{
				print_binary_matrix(get_filename("ne0_", saveNum, 5), from_fourier(nek));
				print_binary_matrix(get_filename("vexbk", saveNum, 5), vexbk);
				print_binary_matrix(get_filename("k", saveNum, 5), fourier_mesh.k);
			}
			// Get all residuals
			auto residualnk = calc_residualn(vexbk, fourier_mesh.k, nek);
			auto residualtek = calc_residualt(veok, Tek, fourier_mesh.k);
			auto residualtik = calc_residualt(viok, Tik, fourier_mesh.k);
			
			// Get all source terms (only density for now is fine)
			auto sourcenk = calc_sourcen(fourier_mesh.ksqu, nek, Dart);
			Matrix<complex> sourcetk(nx, nyk, {0, 0});
					
			// Update variables using RK method
			{
				print_binary_matrix(get_filename("residualnk", saveNum, 5), residualnk);
				print_binary_matrix(get_filename("sourcenk", saveNum, 5), sourcenk);
			}
			nek = RK4(nek_old, dt, residualnk, sourcenk, stage);
			Tik = RK4(Tik_old, dt, residualtik, sourcetk, stage);
			{
				print_binary_matrix(get_filename("ne", saveNum, 5), from_fourier(nek));
			}
			saveNum++;
		}

		std::cout << "Iteration = " << iter << "    t = " << time.back() << "   phi_iter = " << phi_iter << std::endl;

		// Save data every saveFrequency time steps
		if (iter / saveFrequency - saveNum == 0){
			ne = from_fourier(nek);
			Te = from_fourier(Tek);
			Ti = from_fourier(Tik);
			phi = from_fourier(phik);

			print_binary_matrix(get_filename("Te", saveNum, 5), Te);
			print_binary_matrix(get_filename("Ti", saveNum, 5), Ti);
			print_binary_matrix(get_filename("phi", saveNum, 5), phi);
			print_binary_matrix(get_filename("ne", saveNum, 5), ne);
			
			saveNum++;
		}
	}
	return 0;
}
