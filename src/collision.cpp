#include "constants.h"
#include "collision.h"

FrequencyResults CollisionFreqCalculator::calculate_collision_frequencies(const Matrix<complex> &nek, const Matrix<complex> &Tik, const Matrix<complex> &Tek) const{
	//The reason we have to convert all of this to real space is because we can't take square roots or reciprocals in Fourier space
	auto ne = from_fourier(nek);
	auto Ti = from_fourier(Tik);
	auto Te = from_fourier(Tek);

	std::cout << "Ti = " << Ti;

	Matrix<Freq<double>> nu(ne.cols(), ne.rows());
	Matrix<double> isigP(ne.cols(), ne.rows());
	Matrix<double> invn(ne.cols(), ne.rows());

	auto e2 = e * e;
	auto e4 = e2 * e2;
	auto kbeps0 = eps0 * kb;
	auto meeps0 = eps0 * me;
	auto meeps02 = meeps0 * meeps0;
	auto me_over_mi = me / mi;
	auto sqrt_me_over_mi = sqrt(me_over_mi);

	Ti.nan_check();
	Te.nan_check();
	ne.nan_check();

	ne.for_each([&, e2, e4, meeps02, me_over_mi, sqrt_me_over_mi](auto j, auto i, auto nep){
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
		if (nan_check(Vthe)){
			std::stringstream stream;
			stream << "kbeps0 * Tep / (nep * e2) = " << (kbeps0 * Tep / (nep * e2)) << ", Tep = " << Tep << ", nep = " << nep;
			throw std::runtime_error(stream.str());
		}
		
		// Calculate plasma parameter
		auto Lambda = 12 * pi * nep * (lambdaD * lambdaD * lambdaD);

		// Calculate electron-electron collision frequency
		freq.ee = nep * e4 * log(Lambda / 3) / (tau * meeps02 * (Vthe * Vthe * Vthe));

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
