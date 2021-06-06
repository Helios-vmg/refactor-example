#include <cstdlib>
#include <fftw3.h>

extern "C"{

void *fftw_malloc(size_t n){
	return malloc(n);
}

void fftw_free(void *p){
	free(p);
}

fftw_plan fftw_plan_dft_r2c_2d(int, int, double *, fftw_complex *, unsigned){
	return {};
}

fftw_plan fftw_plan_dft_c2r_2d(int, int, fftw_complex *, double *, unsigned){
	return {};
}

fftw_plan fftw_plan_many_dft_r2c(int, const int *, int, double *, const int *, int, int, fftw_complex *, const int *, int, int, unsigned){
	return {};
}

fftw_plan fftw_plan_many_dft_c2r(int, const int *, int, fftw_complex *, const int *, int, int, double *, const int *, int, int, unsigned){
	return {};
}

void fftw_execute(fftw_plan){}
void fftw_destroy_plan(fftw_plan){}
void fftw_cleanup(){}
	
}
