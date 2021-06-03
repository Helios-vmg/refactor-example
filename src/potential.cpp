#include "potential.h"
#include "collision.h"

Matrix<complex> PotentialSourceCalculator::calculate_potential_source(
		const FrequencyResults &fr,
		const Matrix<double> &ksqu,
		const Matrix<complex> &Pik,
		const Matrix<complex> &Pek,
		const Matrix<complex> &isigPk,
		const Matrix<complex2d> &dndk) const{

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

	pikTerm = convolve2d(isigPk, convolve2d(d2Pik, pikTerm));
	pekTerm = convolve2d(isigPk, convolve2d(d2Pek, pekTerm));

	Matrix<complex> ret(pikTerm.cols(), pikTerm.rows());
	ret.for_each([&](auto j, auto i, auto &retp){
		auto &dndkp = dndk.get(j, i);
		
		auto x = dndkp.y * params.u.z * params.B.x;
		auto y = -dndkp.x * params.u.z * params.B.y;
		auto z = (dndkp.x * params.u.y - dndkp.y * params.u.x) * params.B.z;

		retp = x + y + z + pikTerm.get(j, i) + pekTerm.get(j, i);
	});

	return ret;
}
