#include "constants.h"
#include "potential.h"
#include "collision.h"

Matrix<complex> PotentialSourceCalculator::calculate_potential_source(
		const FrequencyResults &fr,
		const Matrix<double> &ksqu,
		const Matrix<complex> &Pik,
		const Matrix<complex> &Pek) const{

	auto &params = *this;
	
	auto d2Pik = laplaciank(Pik, ksqu);
	auto d2Pek = laplaciank(Pek, ksqu);

	auto dummy = from_fourier(d2Pik);

	Matrix<complex> pikTerm(fr.nuk.cols(), fr.nuk.rows());
	auto pekTerm = pikTerm;

	pikTerm.for_each([&params, &pekTerm, &fr](auto j, auto i, auto &pikTermp){
		auto &pekTermp = pekTerm.get(j, i);
		auto &nukp = fr.nuk.get(j, i);
		pikTermp = (nukp.ie - nukp.in) / params.Oci + nukp.ei / params.Oce;
		pekTermp = nukp.ie / params.Oci + (nukp.ei  + nukp.en) / params.Oce;
	});

	pikTerm = convolve2d(fr.isigPk, convolve2d(d2Pik, pikTerm));
	pekTerm = convolve2d(fr.isigPk, convolve2d(d2Pek, pekTerm));

	this->dndk->nan_check();
	nan_check(params.u);
	nan_check(params.B);
	
	Matrix<complex> ret(pikTerm.cols(), pikTerm.rows());
	ret.for_each([&](auto j, auto i, auto &retp){
		auto &dndkp = this->dndk->get(j, i);
		
		auto x = dndkp.y * params.u.z * params.B.x;
		auto y = -dndkp.x * params.u.z * params.B.y;
		auto z = (dndkp.x * params.u.y - dndkp.y * params.u.x) * params.B.z;

		retp = x + y + z + pikTerm.get(j, i) + pekTerm.get(j, i);
	});

	return ret;
}

FourierMesh2d::FourierMesh2d(double Lx, double Ly, size_t w, size_t h, int FourierMeshType){
	this->k = Matrix<point2d>(w, h);

	for (size_t i = 0; i < w / 2; i++)
		for (size_t j = 0; j < h; j++)
			this->k.get(j, i).x = tau * i / Lx;
	for (size_t i = w / 2; i < w; i++){
		auto k_counter = i - w;
		for (size_t j = 0; j < h; j++)
			this->k.get(j, i).x = tau * k_counter / Lx;
	}
		
	// Make ky. Because we are using special FFTs for purely real valued functions, the last dimension (y in this case) only need n/2 + 1 points due to symmetry about the imaginary axis.
	this->k.for_each([Ly](auto j, auto i, auto &p){
		p.y = tau * j / Ly;
	});

	double dx = Lx / w;
	double dy = Ly / ((h - 1) * 2);
	double dx2 = dx / 2;
	double dy2 = dy / 2;
	this->ksqu = Matrix<double>(w, h);
	this->ninvksqu = ksqu;
	if (FourierMeshType == 0){
		this->k.for_each([&](auto j, auto i, const auto &p){
			auto sqnorm = p.sqnorm();
			this->ksqu.get(j, i) = sqnorm;
			this->ninvksqu.get(j, i) = -1.0 / sqnorm;
		});
	}else if (FourierMeshType == 1){
		this->k.for_each([&](auto j, auto i, auto &p){
			//Use approximations for Kx, Ky, K^2
			auto x = sin(p.x * dx2) / dx2;
			auto y = sin(p.y * dx2) / dx2;
			p = { x, y };
			auto sqnorm = p.sqnorm();
			this->ksqu.get(j, i) = sqnorm;
			this->ninvksqu.get(j, i) = -1.0 / sqnorm;
		});
	}

	//Account for the [0][0] point since there is a divide by 0 issue.
	this->ninvksqu.get(0, 0) = -1;
}

double max_absComp(const Matrix<complex> &arr3D){
	double max = arr3D.get(0, 0).sqabs();
	arr3D.for_each_unordered([&max](auto &p){
		max = std::max(max, p.sqabs());
	});
	return sqrt(max);
}

int PotentialSourceCalculator::potentialk(const FrequencyResults &fr, const FourierMesh2d &fm, Matrix<complex> &phik, const Matrix<complex> &potential_source, double err_max, int max_iter){
	double it_error;
	int iterations = 0;
	potential_source.nan_check();
	auto RHS = potential_source;
	RHS.nan_check();
	do{
		auto grad_n_grad_phi = convolve2d(*this->dndk, derivk(phik, fm.k));
		grad_n_grad_phi.nan_check();

		RHS.for_each([&potential_source, &grad_n_grad_phi](auto j, auto i, auto &p){
			complex2d g = grad_n_grad_phi.get(j, i);
			p = potential_source.get(j, i) - g.x - g.y;
		});

		RHS.nan_check();
		
		RHS = convolve2d(RHS, fr.invnk);
		
		auto phik_max_old = max_absComp(phik);
		phik = std::move(RHS);
		phik.elementwise_multiplication(fm.ninvksqu);
		auto phik_max = max_absComp(phik);

		it_error = fabs((phik_max - phik_max_old) / phik_max);
		iterations++;
	}while (it_error > err_max && iterations <= max_iter);
	return iterations;
}
