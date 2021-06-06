#pragma once

#define _USE_MATH_DEFINES // for C++

#include "types.h"
#include <cmath>
#include <algorithm>

static const int nx = 256;
static const int ny = 256;
static const int nyk = ny/2 + 1;
static const int ncomp = 2;
static const double pi = M_PI;
static const double tau = 2 * pi;
static const double saveFrequency = 1.0;    // Save data every this many time steps
static const double dt_max = 0.1;        // Set max allowable time step
static const double tend = 1000.;     //1000.     // Set the end time
static const double err_max = 1.e-8;      // Set max allowable error for potential solver
static const double CFL = 3.;             // Set CFL number. We can just leave this at 3.
static const double Dart = 0;      // Set artifical diffusion constants
// change Dart to 1e3 and 7e3
//Try 7e5
//modify D artif 1e3: different than GDI TGI, func off grid size, tot box size
static const int phi_iter_max = 500;      // Max allowable iterations for potential solver
// Calculated parameters	
static const int iter_max = tend / dt_max; 
//int fftw_threads_init(void); // should only be called once and (in your main()function), performs any one-time initialization required to use threads on your system. It returns zeros if succesful and a non-zero if not (error)
// hyperbolic tan paramters 
static const double L = tau;	//useless 
static const double Lx = 200000., Ly = 80000.; //different
//double Lx = L, Ly = L;
static const double dx = Lx / nx;
static const double dy = Ly / ny;
//I.Cs:  Parameters for initializing hyperbolic tangent IC
static const double xg = Lx * (19. / 24); //xg center of hyperbolic tangent func
static const double m = 0.5;
static const double lg = 12000.; //lg is width of hyperbolic tangent func
static const double o = 1.;
static const double a = (o - m) / 2.; //difference from background?
static const double c = -xg;
static const double d = (o + m) / 2.; //size of the cloud ??
static const double a2 = -a;
static const double c2 = -Lx - c;
static const double d2 = d - m;
static const double b = abs((tanh(log(m / o) / 4) * a + d) * cosh(log(m/o)/4) * cosh(log(m / o) / 4) / (a * lg));
//b = abs(((((tanh(log(m/o)/4)) * a) + d) * pow(cosh(log(m/o)/4), 2))/(a * lg));
static const double bg = 2.0 / lg / lg;

// Set physical constants
static const double e = 1.602E-19;
static const double kb = 1.38E-23;
static const double me = 9.109E-31;
static const double eps0 = 8.854E-12;
static const double nn = 1.E14;   // Neutral particle number density
static const double mi = 16. * 1.672E-27;   // Mass O+
static const double mn = mi;               // Mass O
static const double ri = 152.E-12;          // Effective collision radius
static const double rn = ri;
// Set magnetic field. Also calculate different properties of it
static const point3d B = {0,0,5E-5}; 
static const auto Bmag = B.norm();
static const double B2 = Bmag * Bmag;
// Set neutral wind. Generally, the last element will be 0.
static const point3d u = {-500, 0, 0};
static const double kmax = 1.0 / (2.0 * std::min(dx, dy));
static const double Oci = e * Bmag / mi;
static const double Oce = e * Bmag / me;
