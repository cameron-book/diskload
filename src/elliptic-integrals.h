#ifndef __ELLIPTIC_INTEGRALS__
#define __ELLIPTIC_INTEGRALS__

#include <gsl/gsl_mode.h>
#include <gsl/gsl_complex.h>

gsl_complex gsl_sf_ellint_RJ_z (double x, double y, double z, double p, gsl_mode_t mode);

double gsl_sf_ellint_Kcomp_extended (double k, gsl_mode_t mode);
double gsl_sf_ellint_Pcomp_extended (double k, double n, gsl_mode_t mode);
double gsl_sf_ellint_F_extended (double phi, double k, gsl_mode_t mode);

gsl_complex gsl_sf_ellint_Fz (gsl_complex u, double k, gsl_mode_t mode);

#endif
