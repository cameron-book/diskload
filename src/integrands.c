#include <stdio.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_ellint.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_sum.h>
#include <time.h>

double hyperg_z_GT1 (double a, double b, double c, double z);

double f_hypergeometric_ty (double t, double y)
{
  double z = (1-t*t)*(1-y*y)/((1-t*y)*(1-t*y));

  // BADBAD: this is broken
  if (fabs(z - 1.0) < 1e-10)
    z = 1.0 - 1e-10;

  
  if ((z < -1) || (z > 1))
    return hyperg_z_GT1(0.25, 0.75, 1.0, z)/sqrt(2-2*t*y);
  else
    return gsl_sf_hyperg_2F1(0.25, 0.75, 1.0, z)/sqrt(2-2*t*y);
}

double f_elliptic_original (double t, double y) {
  double z = 0.5*(1 - fabs((t*y-1)/(t-y)));
  double factor = 2.0;
  
  if (z < 0) {
    factor *= sqrt(-1.0/(z-1.0));
    z = z/(z-1);
  }

  if (z > 1) {
    factor /= sqrt(z);
    z = 1.0/z;
  }

  return factor*gsl_sf_ellint_Kcomp( sqrt(z), GSL_PREC_DOUBLE ) / M_PI / sqrt(2.0*fabs(y-t));
}


double f_elliptic_ty (double t, double y) {
  double z = (fabs((t*y - 1)/(t - y)) - 1)/(fabs((t*y - 1)/(t - y)) + 1);

  double factor = 2.0;
  
  if (z < 0) {
    factor *= sqrt(-1.0/(z-1.0));
    z = z/(z-1);
  }

  if (z > 1) {
    factor /= sqrt(z);
    z = 1.0/z;
  }

  return factor*gsl_sf_ellint_Kcomp( sqrt(z), GSL_PREC_DOUBLE ) / M_PI / sqrt(fabs(t-y) + fabs(t*y-1));
}

int main( void ) {
  int i;
  double range = 0.001;
  
  for(i = 0; i<1000; i++ ) {
    double t = cos((double)rand()/(double)(RAND_MAX/range));
    double y = cos((double)rand()/(double)(RAND_MAX/range));

    double e = f_elliptic_original(t,y);
    double g = f_hypergeometric_ty(t,y);

    printf( "%e %e %e %e %e\n", t, y, e, g, (e-g)/g );
  }
}
