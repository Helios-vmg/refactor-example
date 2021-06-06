#include "constants.h"
#include "collision.h"

FrequencyResults CollisionFreqCalculator::calculate_collision_frequencies(const Matrix<complex> &nek, const Matrix<complex> &Tik, const Matrix<complex> &Tek) const{
	//The reason we have to convert all of this to real space is because we can't take square roots or reciprocals in Fourier space
	auto ne = from_fourier(nek);
	auto Ti = from_fourier(Tik);
	auto Te = from_fourier(Tek);

	Matrix<Freq<double>> nu(ne.cols(), ne.rows());
	Matrix<double> isigP(ne.cols(), ne.rows());
	Matrix<double> invn(ne.cols(), ne.rows());

	auto &params = *this;
	
	auto e2 = params.e * params.e;
	auto e4 = e2 * e2;
	auto kbeps0 = params.eps0 * params.kb;
	auto meeps0 = params.eps0 * params.me;
	auto meeps02 = meeps0 * meeps0;
	auto me_over_mi = params.me / params.mi;
	auto sqrt_me_over_mi = sqrt(me_over_mi);

	ne.for_each([&, e2, e4, meeps02, me_over_mi, sqrt_me_over_mi](auto j, auto i, auto nep){
		auto Tip = Ti.get(j, i);
		auto Tep = Te.get(j, i);
		auto &freq = nu.get(j, i);
		
		// Calculate thermal velocities
		auto Vthi = sqrt(2 * params.kb * Tip / params.mi);
		auto Vthe = sqrt(2 * params.kb * Tep / params.me);

		// Calculate ion-neutral and electron-neutral collision frequencies
		freq.in = params.ri + params.rn;
		freq.in *= freq.in;
		freq.in *= params.nn * Vthi * pi;
		freq.en = params.nn * Vthi * pi * params.rn * params.rn;
		
		// Calculate Debye length
		auto lambdaD = sqrt(kbeps0 * Tep / (nep * e2));
		
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
		isigP.get(j, i) = 1.0 / (params.e * (freq.in / params.Oci + freq.en / params.Oce));

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
