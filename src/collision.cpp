#include "constants.h"
#include "collision.h"

const auto e2 = e * e;
const auto e4 = e2 * e2;
const auto kbeps0 = eps0 * kb;
const auto meeps0 = eps0 * me;
const auto meeps02 = meeps0 * meeps0;
const auto me_over_mi = me / mi;
const auto sqrt_me_over_mi = sqrt(me_over_mi);

static int fileidx = 0;

std::string get_filename(const char *s, int n, int w = 8);
void print_binary_matrix(const std::string &filename, const Matrix<double> &m);

FrequencyResults CollisionFreqCalculator::calculate_collision_frequencies(const Matrix<complex> &nek, const Matrix<complex> &Tik, const Matrix<complex> &Tek) const{
	//The reason we have to convert all of this to real space is because we can't take square roots or reciprocals in Fourier space
	auto ne = from_fourier(nek);
	auto Ti = from_fourier(Tik);
	auto Te = from_fourier(Tek);

	Matrix<Freq<double>> nu(ne.cols(), ne.rows());
	Matrix<double> isigP(ne.cols(), ne.rows());
	Matrix<double> invn(ne.cols(), ne.rows());

	Ti.nan_check();
	Te.nan_check();
	ne.nan_check();

	print_binary_matrix(get_filename("calculate_collision_frequencies_ne", fileidx++, 5), ne);

	ne.for_each([&](auto j, auto i, auto nep){
		auto Tip = Ti.get(j, i);
		auto Tep = Te.get(j, i);
		auto &freq = nu.get(j, i);
		
		// Calculate thermal velocities
		auto Vthi = sqrt(2 * kb * Tip / mi);
		if (nan_check(Vthi)){
			std::stringstream stream;
			stream << "2 * kb * Tip / mi = " << (2 * kb * Tip / mi) << ", Tip = " << Tip;
			throw std::runtime_error(stream.str());
		}
		auto Vthe = sqrt(2 * kb * Tep / me);
		if (nan_check(Vthe)){
			std::stringstream stream;
			stream << "2 * kb * Tep / mi = " << (2 * kb * Tep / mi) << ", Tep = " << Tep;
			throw std::runtime_error(stream.str());
		}

		// Calculate ion-neutral and electron-neutral collision frequencies
		freq.in = ri + rn;
		freq.in *= freq.in;
		freq.in *= nn * Vthi * pi;
		freq.en = nn * Vthi * pi * rn * rn;
		
		// Calculate Debye length
		auto lambdaD = sqrt(kbeps0 * Tep / (nep * e2));
		if (nan_check(lambdaD)){
			std::stringstream stream;
			stream << "kbeps0 * Tep / (nep * e2) = " << (kbeps0 * Tep / (nep * e2)) << ", Tep = " << Tep << ", nep = " << nep;
			throw std::runtime_error(stream.str());
		}
		
		// Calculate plasma parameter
		auto Lambda = 12 * pi * nep * (lambdaD * lambdaD * lambdaD);

		if (nan_check(Lambda) || Lambda < 0){
			std::stringstream stream;
			stream << "12 * pi * nep * (lambdaD * lambdaD * lambdaD) = " << Lambda << " < 0, nep = " << nep << ", lambdaD = " << lambdaD;
			throw std::runtime_error(stream.str());
		}

		// Calculate electron-electron collision frequency
		freq.ee = nep * e4 * log(Lambda / 3) / (tau * meeps02 * (Vthe * Vthe * Vthe));

		if (nan_check(freq.ee)){
			std::stringstream stream;
			stream << "nep * e4 * log(Lambda / 3) / (tau * meeps02 * (Vthe * Vthe * Vthe)) = " << freq.ee << std::endl;
			stream << "nep = " << nep << std::endl;
			stream << "Lambda = " << Lambda << std::endl;
			stream << "tau = " << tau << std::endl;
			stream << "meeps02 = " << meeps02 << std::endl;
			stream << "Vthe = " << Vthe << std::endl;
			throw std::runtime_error(stream.str());
		}

		// Calculate ion-ion collision frequency
		freq.ii = freq.ee * sqrt_me_over_mi;
		
		// Calculate ion-electron collision frequency
		freq.ie = freq.ee * 0.5 * me_over_mi;
		
		// Calculate electron-ion collision frequency
		freq.ei = freq.ee;
		
		// Calculate "inverse Pedersen conductivity"
		isigP.get(j, i) = 1.0 / (e * (freq.in / Oci + freq.en / Oce));

		// Calculate the inverse of the density
		// inverse of ne in Fourier space (which is needed for several terms in the temperature equation )	
		invn.get(j, i) = 1.0 / nep;
	});

	FrequencyResults ret;
	ret.nuk = to_fourier(nu);
	ret.isigPk = to_fourier(isigP);
	ret.invnk = to_fourier(invn);

	return ret;
}
