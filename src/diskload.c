#include <stdio.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_monte_vegas.h>
#include <string.h>
#include <gsl/gsl_complex_math.h>
#include <math.h>
#include <gsl/gsl_sum.h>
#include <time.h>
#include "elliptic-integrals.h"
#include "diskload.h"
#include "tanhsinh.h"

double diskload_core_H_slow(double x, double y);

const double Degrees = M_PI / 180.0;

const EarthModel DefaultEarthModel = { 6371.0, 9.8046961 };

gsl_integration_workspace *integration_workspace = NULL;
gsl_integration_workspace *double_integration_workspace = NULL;

#define LEVIN_U_TERMS 500
gsl_sum_levin_u_workspace *levin_u_workspace = NULL;

/**
 * @brief Initialize workspaces for diskload calculations
 *
 * @usage This function must be called before other `diskload_`
 * functions.
 *
 **/
void diskload_initialize(void) {
  integration_workspace = gsl_integration_workspace_alloc (10000);
  levin_u_workspace = gsl_sum_levin_u_alloc (LEVIN_U_TERMS);

  double_integration_workspace = gsl_integration_workspace_alloc (10000);
  
  return;
}


/*----------------------------------------------------------------*/
/** @defgroup errors Error handling
 *  Provide both human-readable strings and error codes
 *  @{
 */

struct {
    DiskLoadError code;
    char *message;
} error_descriptions[] = {
    { E_SUCCESS, "No error" },
    { E_LOVE_NUMBER_VECTOR_TOO_SHORT, "The computation would access entries beyond the end of the Love Number vectors" },
    { E_LOVE_NUMBER_L0_AND_K0_ARE_NONZERO, "The computation assumes l_0 and k_0 are zero" },
    { 0, NULL },
};

DiskLoadError diskload_errno = E_SUCCESS;
#define RETURN_ERROR(e) return (diskload_errno = e);

/**
 * @brief print most recent error
 *
 * @usage Interprets the most recent error as an error message, and
 * prints it to stderr, optionally preceding it with the custom
 * message specified in `str` if `str` is non-`NULL`.
 *
 * @param str the optional string to preceed the output
 *
 **/
void diskload_perror( const char* str ) {
  // Search for the appropriate error message
  int i = 0;
  while ( error_descriptions[i].message != NULL ) {
    if (error_descriptions[i].code == diskload_errno)
      break;
    else
      i++;
  }

  // Print the error message to stderr
  char* e = error_descriptions[i].message;
  if (e == NULL)
    e = "Unknown error";
    
  if (str != NULL)
    fprintf( stderr, "%s: %s\n", str, e );
  else
    fprintf( stderr, "%s\n", e );

  return;
}

/** @} */ // end of errors

/*----------------------------------------------------------------*/
/** @defgroup love Love numbers
 *  Love numbers can be loaded from external files and extrapolated beyond available data
 *  @{
 */

/**
 * @brief Replace occurances of one character with another
 *
 * @usage In `str` replace all occurances of the character `find` with the character `replace`.
 * 
 * @param str the string to be modified
 * @param find the character to search for
 * @param replace the character to replace `find` by
 *
 * @return the string `str`
 *
 **/
char* replace_char(char* str, char find, char replace) {
    char *current_pos = strchr(str,find);
    while (current_pos){
        *current_pos = replace;
        current_pos = strchr(current_pos,find);
    }
    return str;
}

/**
 * @brief Load Love numbers from an external ASCII file.
 *
 * @usage Read Love numbers from `filename` assuming that each line
 * includes four numbers (an integer index, followed by the three Love
 * numbers h, l, k).  Allocate storage space for a LoveNumber.
 *
 * This is the format used in
 * `REF_6371_loading_love_numbers_0_40000.txt` so you can load this
 * file by simply using `LoveNumbers* ln =
 * read_love_numbers("REF_6371_loading_love_numbers_0_40000.txt");`
 * 
 * @param filename the text filename to read
 *
 * @return A newly allocated LoveNumbers struct, or NULL on error. It
 * is the responsibility of the caller to `free` the arrays
 * `LoveNumbers->h` and `LoveNumbers->l` and `LoveNumbers->k` as well
 * as the returned `LoveNumbers` struct itself.
 *
 **/
LoveNumbers* diskload_read_love_numbers(const char* filename) {
  FILE* f = fopen(filename, "r");

  if (ferror(f)) {
    return NULL;
  }
  
  LoveNumbers* result = malloc(sizeof(LoveNumbers));
  if (result == NULL)
    return NULL;
  
  result->degrees = 0;
  result->h = NULL;
  result->l = NULL;
  result->k = NULL; 
  int allocated_size = 0;
  
  const size_t line_size = 255;
  char* line = alloca(line_size);

  double x, y, z;
  int n;

  // Loop over each line in the REF text file
  while (fgets(line, line_size, f) != NULL)  {
    replace_char( line, 'D', 'E' );
    if (sscanf( line, "%d %lf %lf %lf\n", &n, &x, &y, &z ) == 4) {

      // Grow our arrays to a larger size if needed
      if (n >= allocated_size) {
        allocated_size = allocated_size * 2 + 1;
        result->h = realloc( result->h, allocated_size * sizeof(double) );
        result->l = realloc( result->l, allocated_size * sizeof(double) );
        result->k = realloc( result->k, allocated_size * sizeof(double) );

        if ((result->h == NULL) || (result->l == NULL) || (result->k == NULL)) {
          // BADBAD: we leak some memory in this error
          return NULL;
        }
      }

      // Store this line in the struct
      result->h[n] = x;
      result->l[n] = y;
      result->k[n] = z;
      if (n + 1 > result->degrees)
        result->degrees = n + 1;
    }
  }

  return result;
}

/**
 * @brief extrapolate the load Love numbers via asymptotic approximation
 *
 * @usage
 * This function assumes the the elastic load Love numbers `h,k,l` complete 
 * through degree nmax were computed directly, and that the asymptotic 
 * limit has been reached, i.e. for any n > nmax, h(n)=h(nmax), 
 * n k(n)= nmax k(nmax), and n l(n)= nmax l(nmax). The function then
 * computes h,k and l thru degree nmax.
 *
 * @param love the elastic loading Love numbers
 * @param nmax the desired (and larger) count of Love numbers
 *
 **/
void diskload_extrapolate_love_numbers( LoveNumbers* love, int nmax ) {
  int oo = love->degrees-1;
  double h_oo = love->h[oo];
  double l_oo = love->l[oo];
  double k_oo = love->k[oo];  

  int allocated_size = nmax + 1;
  love->h = realloc( love->h, allocated_size * sizeof(double) );
  love->l = realloc( love->l, allocated_size * sizeof(double) );
  love->k = realloc( love->k, allocated_size * sizeof(double) );  
  
  int n;
  for( n=love->degrees; n <= nmax; n++ ) {
    love->h[n] = h_oo;
    love->k[n] = oo * k_oo / n;
    love->l[n] = oo * l_oo / n;
  }

  love->degrees = nmax + 1;

  return;
}

/** @} */ // end of love numbers

/*----------------------------------------------------------------*/
/** @defgroup computations Computations
 *  There are various algorithms for computing the geoelastic response to a uniform spherical disk load
 *  @{
 */

/**
 * @brief compute disk load via a truncated sum
 *
 * @usage  
 * This function computes the response to a uniform surface pressure load 
 * imposed in a disc of angular radius `alpha`. The elastic response is found 
 * at one or more stations located on the surface of the earth at angular 
 * distance(s) theta from the center of the disc load.
 *
 * The elastic response is computed using user-supplied elastic loading 
 * Love numbers (h,k,l) generated using a specific elastic structure model
 * for the earth. If three output arguments are invoked, this function
 * also computes the change in the height of the geoid at each station.
 * The pressure load imposed within the disk is expressed in terms of the 
 * equivalent depth height, thickness) of liquid water (density=1000 kg/m3).
 *
 * This function is a modified verison of function diskload.m
 * associated with the publication \cite originalDiskload
 *
 * @note input w can be positive or negative allowing the user to
 * model the incremental response to incremental pressure changes
 *
 * @note It is easy to switch from the state to the rate problem. If
 * input `h` is actually the rate of change of the load (in m/yr w.e.),
 * then outputs u,v and g will become velocities (in mm/yr).
 *
 * @warning All elements of nmax must be <= n, the maximum order
 * provided for the elastic loading Love numbers
 *
 * @param[in] alpha angular radius of disk stated in degrees
 * @param[in] icomp compute using either a compensated (1) or uncompensated (0) disk load
 * @param[in] theta angular distances of stations from disc center in degrees
 * @param[in] w pressure imposed on the surface within the spherical disk expressed as the height or depth of equivalent water load (in meters)
 * @param[in] nmax maximum harmonic degrees of the expansion to be used. nmax may be a scalar or a vector with multiple truncation points.  If nmax=[], then nmax will be set to 5*(360/alpha) if the LNs are complete thru that degree,so as to satisfy Bevis et al. (2016)s Rule of Thump (ROT) for a good approximation. Otherwise nmax will be set to the highest value of n in LN.
 * @param[in] LN the eastic loading Love numbers (LN), a structure with fields LN.h  (n+1)-vector containing Love number h for degrees 0:n LN.k  (n+1)-vector containing Love number k for degrees 0:n LN.k  (n+1)-vector containing Love number l for degrees 0:n
 * @param[in] EM spherical Earth model parameters; if NULL, defaults to DefaultEarthModel
 * @param[out] u output radial or 'vertical' elastic displacement (mm)
 * @param[out] v output tangential or 'horizontal' elastic displacement (mm)
 * @param[out] g output geoid change (mm) 
 *
 * @return a `DiskLoadError` which is E_SUCCESS if successful
 ****************************************************************************************/
DiskLoadError diskload_truncated(double alpha,DiskLoadType icomp,double theta,double w,int nmax,LoveNumbers* LN,const EarthModel* earth, double *u, double *v, double *g) {

  if (earth == NULL) earth = &DefaultEarthModel;

  // The 2014 CODATA-recommended value of the gravitational constant in SI units
  double newtonsG =6.6740831e-11;

  if (nmax > LN->degrees) {
    RETURN_ERROR(E_LOVE_NUMBER_VECTOR_TOO_SHORT);
  }
  
  if ((LN->l[0] != 0) || (LN->k[0] != 0)) {
    RETURN_ERROR(E_LOVE_NUMBER_L0_AND_K0_ARE_NONZERO);
  }
  
  double x = cos(theta * Degrees);
  double sintheta = sin(theta * Degrees);  
  double cosalpha = cos(alpha * Degrees);

  double p0 = 1.0;
  double p1 = x;
 
  double pp0 = 1.0;
  double pp1 = cosalpha;
  double pp2 = 0.0;
  
  double sigma;
  if (icomp == Compensated) sigma = 0.0;
  if (icomp == Uncompensated) sigma = (1-cosalpha)/2.0;

#ifndef DISKLOAD_NO_U  
  *u = LN->h[0] * sigma * p0;
#endif
#ifndef DISKLOAD_NO_V  
  *v = 0.0;
#endif
#ifndef DISKLOAD_NO_G  
  *g = sigma;
#endif
  
  // keeping track of low-order bits
  double Uc = 0.0, Vc = 0.0, Gc = 0.0;
  
  double dp1, q1, p2, input, t, z;

  for( int n = 1; n <= nmax; n++ ) {
    // Recursively compute the loading factor
    pp2 = ((2*n+1)*cosalpha*pp1 - n*pp0)/(n+1);

    if (icomp == Compensated) sigma = (pp0 - pp2)/(2*n+1)/(1+cosalpha);
    if (icomp == Uncompensated) sigma = (pp0 - pp2)/(2*n+1)/2.0;

#ifndef DISKLOAD_NO_U      
    // Kahan summation for U
    z = (LN->h[n] * sigma * p1) - Uc;
    t = *u + z;
    Uc = (t - *u) - z;
    *u = t;
#endif

    // Compute the derivative of P_n(cos(x)) with respect to x
    dp1=(-n/sintheta)*(p0 - x*p1);

#ifndef DISKLOAD_NO_V    
    // Kahan summation for V
    z = (LN->l[n] * sigma * dp1) - Vc;
    t = *v + z;
    Vc = (t - *v) - z;
    *v = t;
#endif

#ifndef DISKLOAD_NO_G    
    // Kahan summation for the geoid term
    z = ((1.0 + LN->k[n]) * sigma * p1) - Gc;
    t = *g + z;
    Vc = (t - *g) - z;
    *g = t;
#endif
    
    pp0 = pp1;
    pp1 = pp2;
    
    // Recursively compute the Legendre polynomial values
    p2 = ((2*n+1)*x*p1 - n*p0)/(n+1);
    p0 = p1;
    p1 = p2;
  }

  double radiusMeters = earth->radius * 1000;
  double rhoWater = 1000; // density of pure water(kg/m^3)
  double rhoEarth = 3.0*earth->gravity/4.0/newtonsG/M_PI/radiusMeters; // Average Earth density in (kg/m^3)
  double metersToMillimeters = 1000;
  double outputScale = (3*rhoWater/rhoEarth) * w * metersToMillimeters;

#ifndef DISKLOAD_NO_U  
  *u = *u * outputScale;
#endif
#ifndef DISKLOAD_NO_V  
  *v = *v * outputScale;
#endif
#ifndef DISKLOAD_NO_G  
  *g = *g * outputScale;
#endif
  
  return 0;
}

/**
 * @brief compute disk load via a computational "core"
 *
 * @usage This function is exactly like diskload_truncated except for
 * the additional parameters `coreG` and `coreH` which point to
 * functions providing the "core" contribution to the infinite sum
 * computing the disk load.
 *
 * @return a `DiskLoadError` which is E_SUCCESS if successful
 ****************************************************************************************/
DiskLoadError diskload_core(double alpha,DiskLoadType icomp,double theta,double w,ComputationalCore coreG, ComputationalCore coreH,int nmax,LoveNumbers* LN,const EarthModel* earth, double *u, double *v, double *g) {

  if (earth == NULL) earth = &DefaultEarthModel;

  // The 2014 CODATA-recommended value of the gravitational constant in SI units
  double newtonsG =6.6740831e-11;

  if (nmax > LN->degrees) {
    RETURN_ERROR(E_LOVE_NUMBER_VECTOR_TOO_SHORT);
  }
  
  if ((LN->l[0] != 0) || (LN->k[0] != 0)) {
    RETURN_ERROR(E_LOVE_NUMBER_L0_AND_K0_ARE_NONZERO);
  }
  
  double x = cos(theta * Degrees);
  double sintheta = sin(theta * Degrees);  
  double cosalpha = cos(alpha * Degrees);

  double p0 = 1.0;
  double p1 = x;
 
  double pp0 = 1.0;
  double pp1 = cosalpha;
  double pp2 = 0.0;

#if !defined(DISKLOAD_NO_U) || !defined(DISKLOAD_NO_G)
  double coreValue = (*coreG)(alpha * Degrees,theta * Degrees);
#endif
  
#ifndef DISKLOAD_NO_V    
  double coreDerivative = (*coreH)(alpha * Degrees,theta * Degrees);
  printf( "coreDerivative = %g\n", coreDerivative );
#endif
  
  double sigma;
  if (icomp == Compensated) sigma = 0.0;
  if (icomp == Uncompensated) sigma = (1.0-cosalpha)/2.0;

  double coreFactor;
  if (icomp == Compensated) coreFactor = 1.0+cosalpha;
  if (icomp == Uncompensated) coreFactor = 2.0;

#ifndef DISKLOAD_NO_U
  *u = LN->h[nmax] * coreValue / coreFactor + (LN->h[0] - LN->h[nmax]) * sigma * p0;
#endif

#ifndef DISKLOAD_NO_V  
  double l_oo = (nmax) * LN->l[nmax];
  *v = - sintheta * l_oo * coreDerivative / coreFactor;
  printf( "initial v = %g\n", *v );
#endif

#ifndef DISKLOAD_NO_G
  *g = coreValue / coreFactor;
#endif
  
  if (icomp == Compensated) {
    // BADBAD: this differs from the paper -- we should be subtracting the first term?
    *u -= LN->h[nmax] * (1 - cosalpha) / coreFactor;
    *g -= (1 - cosalpha) / coreFactor;
  }
  
  // keeping track of low-order bits
  double Uc = 0.0, Vc = 0.0, Gc = 0.0;
  
  double dp1, q1, p2, input, t, z;

  for( int n = 1; n <= nmax; n++ ) {
    // Recursively compute the loading factor
    pp2 = ((2*n+1)*cosalpha*pp1 - n*pp0)/(n+1);

    sigma = (pp0 - pp2)/(2*n+1)/coreFactor;

#ifndef DISKLOAD_NO_U
    // Kahan summation for U
    z = ((LN->h[n] - LN->h[nmax]) * sigma * p1) - Uc;
    t = *u + z;
    Uc = (t - *u) - z;
    *u = t;
#endif

    // Compute the derivative of P_n(cos(x)) with respect to x
    dp1=(-n/sintheta)*(p0 - x*p1);

#ifndef DISKLOAD_NO_V
    // Kahan summation for V
    // this is /n and not /(n+1)
    z = ((LN->l[n] - l_oo/n) * sigma * dp1) - Vc;
    t = *v + z;
    Vc = (t - *v) - z;
    *v = t;
#endif

#ifndef DISKLOAD_NO_G    
    // Kahan summation for the geoid term
    z = (LN->k[n] * sigma * p1) - Gc;
    t = *g + z;
    Vc = (t - *g) - z;
    *g = t;
#endif
    
    pp0 = pp1;
    pp1 = pp2;
    
    // Recursively compute the Legendre polynomial values
    p2 = ((2*n+1)*x*p1 - n*p0)/(n+1);
    p0 = p1;
    p1 = p2;
  }

  double radiusMeters = earth->radius * 1000;
  double rhoWater = 1000; // density of pure water(kg/m^3)
  double rhoEarth = 3.0*earth->gravity/4.0/newtonsG/M_PI/radiusMeters; // Average Earth density in (kg/m^3)
  double metersToMillimeters = 1000;
  double outputScale = (3*rhoWater/rhoEarth) * w * metersToMillimeters;

#ifndef DISKLOAD_NO_U  
  *u = *u * outputScale;
#endif
#ifndef DISKLOAD_NO_V
  *v = *v * outputScale;
#endif
#ifndef DISKLOAD_NO_G
  *g = *g * outputScale;
#endif
  
  return 0;
}

double hyperg_z_GT1 (double a, double b, double c, double z) {
  double coef1,coef2;
  
  coef1=gsl_sf_gamma(c)*gsl_sf_gamma(b-a)*pow(1-z,-a)/(gsl_sf_gamma(b)*gsl_sf_gamma(c-a));
  coef2=gsl_sf_gamma(c)*gsl_sf_gamma(a-b)*pow(1-z,-b)/(gsl_sf_gamma(a)*gsl_sf_gamma(c-b));
  double result = coef1*gsl_sf_hyperg_2F1(a,c-b,a-b+1,1.0/(1-z))+coef2*gsl_sf_hyperg_2F1(b,c-a,b-a+1,1.0/(1-z));

  return result;
}

double f_hypergeometric (double t, void *params)
{
  double y = *(double *) params;
  double z = (1-t*t)*(1-y*y)/((1-t*y)*(1-t*y));

  // BADBAD: this is broken
  /*
  if (fabs(z - 1.0) < 1e-10)
    z = 1.0 - 1e-10;
  */

  if ((z >= 1.0) && (z < 1.0 + 1e-10))
    z = 1.0 + 1e-10;

  if ((z < 1.0) && (z > 1.0 - 1e-10))
    z = 1.0 - 1e-10;  
    
  
  if ((z < -1) || (z > 1))
    return hyperg_z_GT1(0.25, 0.75, 1.0, z)/sqrt(2-2*t*y);
  else
    return gsl_sf_hyperg_2F1(0.25, 0.75, 1.0, z)/sqrt(2-2*t*y);
}

double f_elliptic (double t, void *params) {
  double y = *(double *) params;    
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

double diskload_core_terms[LEVIN_U_TERMS];

double diskload_core_H_bad(double x,double y) {
  return diskload_core_H_truncated(x,y) - 0.000001149195;
}

double diskload_core_H_truncated(double x,double y) {
  x = cos(x);
  y = cos(y);
  
  int nmax = 400000;
  double result = 0.0;
  int n;

  double p0 = 1.0;
  double p1 = x;
  double p2 = (3.0*x*x - 1.0)/2;

  double pp0 = 1.0;
  double pp1 = y;
  double pp2 = 0.0;

  double q1;

  result = 0.0;

  // This looks like an off-by-one error, but the first term is
  // actually 0, and we don't include it in the accelerated series
  for( n=1; n <= nmax; n++ ) {
    p2 = ((2*n+1)*x*p1 - n*p0)/(n+1);    
    q1 = (p0 - p2)/(2*n+1);

    double dp1 = (y*pp1 - pp0) * n / (y*y - 1);

    result += q1 * dp1 / (n);

    p0 = p1;
    p1 = p2;
    
    pp2 = ((2*n+1)*y*pp1 - n*pp0)/(n+1);
    pp0 = pp1;
    pp1 = pp2;
  }
  
  return result;
}


double diskload_core_H_truncated_poorly(double x,double y) {
  //int nmax = LEVIN_U_TERMS;
  int nmax = 80000;
  double result = 0.0;
  int n;

  double p0 = 1.0;
  double p1 = x;
  double p2 = (3.0*x*x - 1.0)/2;

  double pp0 = 1.0;
  double pp1 = y;
  double pp2 = 0.0;

  double q1;

  result = 0.0;

  double largest = -1e300;
  double smallest = 1e300;  
  double cesaro = 0.0;
  
  // This looks like an off-by-one error, but the first term is
  // actually 0, and we don't include it in the accelerated series
  for( n=1; n <= nmax; n++ ) {
    p2 = ((2*n+1)*x*p1 - n*p0)/(n+1);    
    q1 = (p0 - p2)/(2*n+1);

    double dp1 = (y*pp1 - pp0) * n / (y*y - 1);

    //diskload_core_terms[n-1] = q1 * (y*pp1 - pp0) / (y*y - 1) * ((double)n / (n+1));
    
    //cesaro += ((double)nmax/(n+1) - 1) * q1 * dp1;
    result += q1 * dp1 / (n+1);
    cesaro += result;
    if (n > 1000) {
      if (result > largest)
        largest = result;
      if (result < smallest)
        smallest = result;
    }

    p0 = p1;
    p1 = p2;
    
    pp2 = ((2*n+1)*y*pp1 - n*pp0)/(n+1);
    pp0 = pp1;
    pp1 = pp2;
  }

  double sum_accel, err;
  //  gsl_sum_levin_u_accel (diskload_core_terms, nmax, levin_u_workspace,
  //&sum_accel, &err );

  printf( "result = %e\n", result );
  printf( "cesaro = %e\n", cesaro / nmax );  
  printf( "large  = %e\n", largest );
  printf( "small  = %e\n", smallest );
  printf( "average  = %e\n", (largest + smallest)/2.0 );
  printf( "geom  = %e\n", sqrt(largest * smallest) );
  //printf( "other  = %e\n", levin_u_workspace->sum_plain );
  //printf( "sumacc = %e\n", sum_accel );
  //printf( "error  = %e\n", err );
  printf( "truth  = %e\n", diskload_core_H_truncated(x,y) );
  
  return result / nmax;
}

double diskload_hypergeometric_core_G(double x,double y) {
  x = cos(x);
  y = cos(y);
  
  double lower_limit = x;	/* lower limit a */
  double upper_limit = 1;	/* upper limit b */
  double abs_error = 1.0e-6;	/* to avoid round-off problems */
  double rel_error = 1.0e-6;	/* the result will usually be much better */
  double result;		/* the result from the integration */
  double error;			/* the estimated error from the integration */
  
  gsl_function My_function;
  void *params_ptr = &y;
  My_function.function = &f_elliptic;
  My_function.params = params_ptr;
  
  double pts[3];
  int npts;

  if ((y > x) && (y < 1)) {
    pts[0] = lower_limit;
    pts[1] = y;
    pts[2] = upper_limit;
    npts = 3;
  } else {
    pts[0] = lower_limit;
    pts[1] = upper_limit;
    npts = 2;    
  }

  gsl_integration_qagp (&My_function, pts, npts,
			abs_error, rel_error, 1000,
                        integration_workspace, &result,
			&error);

  return result;
}

double diskload_core_G(double x, double y) {
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  
  x = x/2;
  y = y/2;  
  
  long double factor = (cosl(x)/sinl(y) + sinl(y)/cosl(x) - (cosl(y)/tanl(y)/cosl(x) + sinl(x)*tanl(x)/sinl(y)))/M_PI;
  long double k = tanl(x)/tanl(y);
  gsl_complex u = gsl_complex_arcsin_real(tanl(y)/tanl(x));
  long double n = (sinl(x)/sinl(y))*(sinl(x)/sinl(y));

  long double complete = gsl_sf_ellint_Kcomp_extended(k, mode) - 
                         gsl_sf_ellint_Pcomp_extended(k, -n, mode);
  
  printf( "mpmath.ellipk(%.12Le) = %.12e\n", k*k, gsl_sf_ellint_Kcomp_extended(k,mode) );
  printf( "mpmath.ellippi(%.12e,%.12e) = %.12e\n", (double)(n), (double)(k*k), gsl_sf_ellint_Pcomp_extended(k,-n,mode) );

  long double incomplete = GSL_REAL(gsl_sf_ellint_Fz(u,k,mode)) -
                           gsl_sf_ellint_Pcomp_extended(1.0/k,-n/k/k,mode) / k;

  printf( "mpmath.ellipf(%.12e + %.12e * I,%.12Le) = %.12e\n", GSL_REAL(u), GSL_IMAG(u), k*k, GSL_REAL(gsl_sf_ellint_Fz(u,k,mode)) );

  printf( "mpmath.ellippi(%.12e,%.12e)/%.12e = %.12e\n",
          (double)(n/k/k), (double)(1.0/k/k), (double)(k),
          (double)(gsl_sf_ellint_Pcomp_extended(1.0/k,-n/k/k,mode) / k) );
  printf( "gsl_sf_ellint_Pcomp_extended(%.12e,%.12e,mode)\n", 
          (double)(1.0/k),(double)(-n/k/k) );

  printf( "1/k = %.12Le\n", 1.0/k );
  printf( "n/k/k = %.12Le\n", n/k/k );
  
  printf( "incomplete = %.12Le\n", incomplete );
  
  return 1.0 - factor*(incomplete + complete);
}


/**
 * @brief compute disk load via integrals of hypergeometric functions
 *
 * @usage This function is numerically akin to diskload_truncated, but
 * uses a rather different algorithm internally.
 *
 * @return a `DiskLoadError` which is E_SUCCESS if successful
 ****************************************************************************************/
DiskLoadError diskload_hypergeometric(double alpha,DiskLoadType icomp,double theta,double w,int nmax,LoveNumbers* LN,const EarthModel* earth, double *u, double *v, double *g) {
  return diskload_core(alpha,icomp,theta,w,diskload_core_G, diskload_core_H, nmax,LN,earth,u,v,g);
}

/**
 * @brief compute point load 
 *
 * @usage This function is numerically akin to diskload_truncated, but
 * uses a rather different algorithm internally.
 *
 * @return a `DiskLoadError` which is E_SUCCESS if successful
 ****************************************************************************************/
DiskLoadError diskload_point(double theta,double w,int nmax,LoveNumbers* LN,const EarthModel* earth, double *u, double *v, double *g) {
  if (earth == NULL) earth = &DefaultEarthModel;

  // The 2014 CODATA-recommended value of the gravitational constant in SI units
  double newtonsG =6.6740831e-11;

  if (nmax > LN->degrees) {
    RETURN_ERROR(E_LOVE_NUMBER_VECTOR_TOO_SHORT);
  }
  
  if ((LN->l[0] != 0) || (LN->k[0] != 0)) {
    RETURN_ERROR(E_LOVE_NUMBER_L0_AND_K0_ARE_NONZERO);
  }
  
  double alpha = 0.0;
  double x = cos(theta * Degrees);
  double sintheta = sin(theta * Degrees);  
  double cosalpha = 1.0;

  double p0 = 1.0;
  double p1 = x;
 
  double pp0 = 1.0;
  double pp1 = cosalpha;
  double pp2 = 0.0;

  double coreValue = 1.0 / (2*sin(theta*Degrees/2.0));
  //double coreDerivative = (*coreH)(alpha * Degrees,theta * Degrees);
  double coreDerivative = 17;
  
  double sigma = 0.0;
  double coreFactor = 2.0;

  *u = (LN->h[nmax] / (2.0*sin(theta*Degrees/2.0)) + (LN->h[0] - LN->h[nmax]) * p0)/coreFactor;
  double l_oo = (nmax + 1) * LN->l[nmax];
  *v = - sintheta * l_oo * coreDerivative / coreFactor;
  *g = coreValue / coreFactor;
   
  // keeping track of low-order bits
  double Uc = 0.0, Vc = 0.0, Gc = 0.0;
  
  double dp1, q1, p2, input, t, z;

  for( int n = 1; n <= nmax; n++ ) {
    // Recursively compute the loading factor
    pp2 = ((2*n+1)*cosalpha*pp1 - n*pp0)/(n+1);
    //printf( "%g\n", pp2 );

    sigma = (pp0 - pp2)/(2*n+1)/coreFactor;
    
    // Kahan summation for U
    z = ((LN->h[n] - LN->h[nmax]) * p1 / coreFactor) - Uc;
    t = *u + z;
    Uc = (t - *u) - z;
    *u = t;

    // Compute the derivative of P_n(cos(x)) with respect to x
    dp1=(-n/sintheta)*(p0 - x*p1);

    // Kahan summation for V
    z = ((LN->l[n] - l_oo/(n+1)) * sigma * dp1) - Vc;
    t = *v + z;
    Vc = (t - *v) - z;
    *v = t;
        
    // Kahan summation for the geoid term
    z = (LN->k[n] * sigma * p1) - Gc;
    t = *g + z;
    Vc = (t - *g) - z;
    *g = t;

    pp0 = pp1;
    pp1 = pp2;
    
    // Recursively compute the Legendre polynomial values
    p2 = ((2*n+1)*x*p1 - n*p0)/(n+1);
    p0 = p1;
    p1 = p2;
  }

  double radiusMeters = earth->radius * 1000;
  double rhoWater = 1000; // density of pure water(kg/m^3)
  double rhoEarth = 3.0*earth->gravity/4.0/newtonsG/M_PI/radiusMeters; // Average Earth density in (kg/m^3)
  double metersToMillimeters = 1000;
  double outputScale = (3*rhoWater/rhoEarth) * w * metersToMillimeters;

  *u = *u * outputScale;
  *v = *v * outputScale;
  *g = *g * outputScale;
  
  return 0;
}  
 


/** @} */ // end of computations

double F(double x,double y,double z) {
  double d = 1 - 2*x*y*z + z*z;
  double w = 4.0*(1-x*x)*(1-y*y)*z*z / (d*d);

  if ((w >= 1.0) && (w < 1.0 + 1e-10))
    w = 1.0 + 1e-10;

  if ((w < 1.0) && (w > 1.0 - 1e-10))
    w = 1.0 - 1e-10;    
  
  if ((w < -1) || (w > 1))
    return hyperg_z_GT1(0.25, 0.75, 1.0, w)/sqrt(d);
  else
    return gsl_sf_hyperg_2F1(0.25, 0.75, 1.0, w)/sqrt(d);
}

double FF(double t, void *params) {
  double x = ((double *) params)[0];
  double y = ((double *) params)[1];  
  double z = t;
  double result = F(x,y,z);
  //printf( "(%f,%f,%f) = %f\n", x, y, z, result );
  return result;
}

double int_f_dz (double t, void *params) {
  double y = *(double *) params;    

  double lower_limit = 0;	/* lower limit a */
  double upper_limit = 1;	/* upper limit b */
  double abs_error = 1.0e-4;	/* to avoid round-off problems */
  double rel_error = 1.0e-4;	/* the result will usually be much better */
  double result;		/* the result from the integration */
  double error;			/* the estimated error from the integration */

  gsl_function My_function;
  double params_array[2] = { t, y };
  void *params_ptr = params_array;
  My_function.function = &FF;
  My_function.params = params_ptr;
  
  gsl_integration_qag (&My_function, lower_limit, upper_limit,
                       abs_error, rel_error, 1000, GSL_INTEG_GAUSS61,
                       double_integration_workspace, &result,
                       &error);
  return result; 
}

double H(double y, void *params) {
  double x = ((double *) params)[0];
  
  double lower_limit = x;	/* lower limit a */
  double upper_limit = 1;	/* upper limit b */
  double abs_error = 1.0e-4;	/* to avoid round-off problems */
  double rel_error = 1.0e-4;	/* the result will usually be much better */
  double result;		/* the result from the integration */
  double error;			/* the estimated error from the integration */

  gsl_function My_function;
  void *params_ptr = &y;
  My_function.function = &int_f_dz;
  My_function.params = params_ptr;
  
  double pts[3];
  int npts;

  if ((y > x) && (y < 1)) {
    pts[0] = lower_limit;
    pts[1] = y;
    pts[2] = upper_limit;
    npts = 3;
  } else {
    pts[0] = lower_limit;
    pts[1] = upper_limit;
    npts = 2;    
  }

  gsl_integration_qagp (&My_function, pts, npts,
			abs_error, rel_error, 1000,
                        integration_workspace, &result,
			&error);

  return result; 
}

double diskload_core_H_deriv(double x, double y) {
  x = cos(x);
  y = cos(y);
  
  gsl_function F;
  
  double result, abserr;
  F.function = &H;
  void *params_ptr = &x;  
  F.params = params_ptr;

  gsl_deriv_backward (&F, y, 1e-8, &result, &abserr);
  
  return result;

}

double integrand (double t, const void *params) {
  double x = ((double*)params)[0];
  double y = ((double*)params)[1];
  
  /*
    Effectively this is:
    
  if ((x < t) && (y < t)) 
    return t*sqrt((-x + t)/(-y + t))/sqrt(1-t*t);
  
  if ((x > t) && (y > t)) 
    return -t*sqrt((-x + t)/(-y + t))/sqrt(1-t*t);

  return -sqrt(-(-x + t)/(-y + t));
  */

  double sign = 1.0;
  if (signbit(t)) sign = -1.0;
  if ((x > t) && (y > t))
    sign = -sign;
  
  if (((x < t) && (y < t)) || ((x > t) && (y > t)))
    return sign*sqrt(t*t*(-x + t)/(-y + t)/(1-t*t));

  return -sqrt(-(-x + t)/(-y + t));

}

// https://github.com/Rufflewind/tanhsinh

double diskload_core_M_slow(double x, double y) {
  x = cos(x);
  y = cos(y);
  
  double abs_error = 1.0e-6;	/* to avoid round-off problems */  
  double params[3] = {x,y};

  double result = 0.0;
  result += tanhsinh_quad(integrand, params,
                          -1 + DBL_EPSILON, fmin(x,y) - DBL_EPSILON, abs_error,
                          NULL, NULL) / M_PI;
  
  result += -fabs(x-y)/2.0;
  
  result += tanhsinh_quad(integrand, params,
                          fmax(x,y) + DBL_EPSILON, 1.0 - DBL_EPSILON, abs_error,
                          NULL, NULL) / M_PI;
  return result;
}

// Provided x or y is negative, this is the integral of sqrt((t*t*(-x + t)/(-y + t)/(1-t*t))) from t = -1 to t = min(x,y); note that side(-x,-y) is the integral of the same from t = max(x,y) to 1
double one_side(double x,double y) {
  const gsl_mode_t mode = GSL_PREC_DOUBLE;

  double m = ((-1 + x)*(1 + y))/((1 + x)*(-1 + y));
  double k = sqrt(m);
  gsl_complex u = gsl_complex_arcsin_real(sqrt(((1 + x)*(-1 + y))/((-1 + x)*(1 + y))));

  // EllipticE[ArcCsc[Sqrt[m]], m] // FunctionExpand  
  //double EE = mpmath.ellipe(phi, m);
  //double EE = k*(gsl_sf_ellint_Ecomp_extended(1/k,mode) + (-1 + 1/m)*gsl_sf_ellint_Kcomp_extended(1/k,mode));
  double EE = gsl_sf_ellint_Ecomp_extended(k,mode);

  double EF = GSL_REAL(gsl_sf_ellint_Fz(u,k,mode));

  double n = (-1 + x)/(-1 + y);
  //double EPI = mpmath.ellippi(n, phi, m);
  // (EllipticPi[n/m, m^(-1)]/Sqrt[m]) == EllipticPi[n, ArcCsc[Sqrt[m]], m]
  double EPI = gsl_sf_ellint_Pcomp_extended(1.0/k,-n/k/k,mode) / k;
  
  double result = -(EE*(1 + x)*(-1 + y) + (x - y)*(EF + EPI*(x-y) + EF*y))/sqrt(-((1 + x)*(-1 + y)));

  return result;
}

double diskload_core_M(double x, double y) {
  if ((x == 0) && (y == 0)) return 2;
  
  x = cos(x);
  y = cos(y);  
  
  double middle = -M_PI * fabs(x-y)/2.0;

  return (middle + one_side(x,y) + one_side(-x,-y)) / M_PI;
}

double diskload_core_H(double x, double y) {
  long double G = diskload_core_G(x,y);
  printf( "gdiff = %.12Lf\n", G - 0.00025156762368938846588587262536634 );
  long double M = diskload_core_M(x,y);
  printf( "mdiff = %.12Lf\n", M - 0.000251567593076294593039689259839147 );
  long double factor = -sinl(y)*tanl(y);
  long double sine2 = sinl(y)*sinl(y);
  
  //x = cos(x);
  //y = cos(y);

  printf("LHS=%.12Lf\n", (G - (1-cosl(x)))/factor );
  printf("RHS=%.12Lf\n", M/sine2 );

  printf("sin*Gblagh=%.12Lf\n", sine2*(G - (1-cosl(x)))/factor );
  printf("M=%.12Lf\n", M );

  //return (- (G - (1-cosl(x)))/tanl(y) + M/sinl(y))/sinl(y);
  return (M - G*cosl(y) + cosl(y) - cosl(x)*cosl(y))/sine2;
}

double diskload_core_M_together(double x, double y) {
  if ((x == 0) && (y == 0)) return 2;

  long double result = -fabsl(cosl(x)-cosl(y))/2.0;
  
  const gsl_mode_t mode = GSL_PREC_DOUBLE;

  x = x/2;
  y = y/2;  
  
  long double m = tanl(x)*tanl(x)/tanl(y)/tanl(y);
  long double k = fabsl(tanl(x)/tanl(y));

  gsl_complex u = gsl_complex_arcsin_real(tanl(y)/tanl(x));
  long double EE = gsl_sf_ellint_Ecomp_extended(k,mode);
  long double EF = GSL_REAL(gsl_sf_ellint_Fz(u,k,mode));
  long double n = sinl(x)*sinl(x)/sinl(y)/sinl(y);
  long double EPI = gsl_sf_ellint_Pcomp_extended(1.0/k,-n/k/k,mode) / k;
  
  // same thing with -x and -y
  long double EE2 = gsl_sf_ellint_Ecomp_extended(1.0/k,mode);

  u = gsl_complex_arcsin_real(tanl(x)/tanl(y));  
  long double EF2 = GSL_REAL(gsl_sf_ellint_Fz(u,1.0/k,mode));

  n = cosl(x)*cosl(x)/cosl(y)/cosl(y);
  long double EPI2 = gsl_sf_ellint_Pcomp_extended(k,-n*k*k,mode) * k;

  return result + (fabsl(cosl(x))*fabsl(sinl(y))*
                   (EE2 - EPI2*cosl(2*x)*cosl(2*x) + 
          cosl(2*x)*(-EE2 + EF2 - 
             (EE2 + EF2 - 2*EPI2)*cosl(2*y)) + 
          cosl(2*y)*(EE2 - EF2 + (EF2 - EPI2)*cosl(2*y))) + 
       fabsl(cosl(y))*fabsl(sinl(x))*
        (-((cosl(2*x) - cosl(2*y))*
             (EF + EPI*cosl(2*x) + (EF - EPI)*cosl(2*y))) + 
         4*EE*cosl(x)*cosl(x)*sinl(y)*sinl(y)))/
      (2.*M_PI*fabsl(cosl(x)*cosl(y)*sinl(x)*sinl(y)));
}


////////////////////////////////////////////////////////////////

double diskload_core_H_separated(double x, double y) {
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  long double factor = -sinl(y)*tanl(y);
  long double sine2 = sinl(y)*sinl(y);

  long double G = diskload_core_G(x,y);
  long double M = diskload_core_M(x,y);

  return (G - (1-cosl(x)))/factor + M/sine2;
}

double diskload_core_H_fast(double x, double y) {
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  double factor2 = -sin(y)*tan(y);
  double sine2 = sin(y)*sin(y);
  
  x = x/2;
  y = y/2;  

  double k = fabs(tan(x)/tan(y));
  //printf( "k = %.12f\n", k );
  gsl_complex u = gsl_complex_arcsin_real(tan(y)/tan(x));
  double n = (sin(x)/sin(y))*(sin(x)/sin(y));
  //printf( "n = %.12f\n", n );
  
  double EK = gsl_sf_ellint_Kcomp_extended(k,mode);
  double EE = gsl_sf_ellint_Ecomp_extended(k,mode);
  
  double EF = GSL_REAL(gsl_sf_ellint_Fz(u,k,mode));
  n = sin(x)/sin(y)*sin(x)/sin(y);
  double EPI = gsl_sf_ellint_Pcomp_extended(1.0/k,-n/k/k,mode) / k;
  double EPI3 = gsl_sf_ellint_Pcomp_extended(k,-n,mode);

  // same thing with -x and -y
  double EE2 = gsl_sf_ellint_Ecomp_extended(1.0/k,mode);
  u = gsl_complex_arcsin_real(tan(x)/tan(y));  
  double EF2 = GSL_REAL(gsl_sf_ellint_Fz(u,1.0/k,mode));
  n = cos(x)*cos(x)/cos(y)/cos(y);
  double EPI2 = gsl_sf_ellint_Pcomp_extended(k,-n*k*k,mode) * k;
  
  double middle = -fabs(cos(2*x)-cos(2*y))/2.0/sin(2*y)/sin(2*y);

  double M = (EF-EPI+EK-EPI3)*(cos(2*x)-cos(2*y))*cos(2*y)/sin(y)/sin(y)/sin(y)/cos(x)/cos(y)/cos(y)/4/M_PI + (fabs(cos(x))*fabs(sin(y))*(EE2 - EPI2*cos(2*x)*cos(2*x) + 
          cos(2*x)*(-EE2 + EF2 - 
             (EE2 + EF2 - 2*EPI2)*cos(2*y)) + 
          cos(2*y)*(EE2 - EF2 + (EF2 - EPI2)*cos(2*y))) + 
       fabs(cos(y))*fabs(sin(x))*
        (-((cos(2*x) - cos(2*y))*
             (EF + EPI*cos(2*x) + (EF - EPI)*cos(2*y))) + 
         4*EE*cos(x)*cos(x)*sin(y)*sin(y)))/
             (2*M_PI*fabs(4*cos(x)*cos(y)*cos(y)*cos(y)*sin(x)*sin(y)*sin(y)*sin(y)));

  //printf("factortest=%.12f\n", (-(cos(2*x)-cos(2*y))*(cos(2*x)-cos(2*y))/sin(y)/sin(y)/2.0) );
  
      M = ((2*(EE2*(2*sin(x)*sin(x)/tan(y)/tan(y)) + 
                  EPI2*(-(cos(2*x)-cos(2*y))*(cos(2*x)-cos(2*y))/sin(y)/sin(y)/2.0) + 
                  (EF2)*(cos(2*x) - cos(2*y)) 
               ))/cos(y)/cos(y)/cos(y)/sin(x) + 
           ((-cos(2*x) + cos(2*y))/cos(y)*
          (2*(-EF - EK + EPI + EPI3)*cos(2*y) + 
           (EF + EF*cos(2*y) + (EPI)*(cos(2*x)-cos(2*y)))) + 
            4*EE*cos(x)*cos(x)*sin(y)*tan(y))/sin(y)/sin(y)/sin(y)/cos(x)/cos(y))/(8*M_PI);

      double lhs = EE2*4*sin(x)/sin(y)/sin(y)/cos(y) -
                   EPI2*(cos(2*x)-cos(2*y))*(cos(2*x)-cos(2*y))/sin(x)/sin(y)/sin(y)/cos(y)/cos(y)/cos(y) +
                   2*(EF2)*(cos(2*x) - cos(2*y))/cos(y)/cos(y)/cos(y)/sin(x);
      double rhs = ((-cos(2*x) + cos(2*y))*
                    (2*(- EK + EPI3)*cos(2*y) + 
                     EF*(1-cos(2*y)) +
                     (EPI)*(cos(2*x)+cos(2*y)))/cos(y) + 
                    4*EE*cos(x)*cos(x)*sin(y)*tan(y))
                   /sin(y)/sin(y)/sin(y)/cos(x)/cos(y);

      double sum = 0.0;
      double diff = -cos(2*x) + cos(2*y);
      
      //printf( "Sum = %.12f\n", sum );
            
      //printf( "Sum = %.12f\n", sum );
      sum += EE2*4*sin(x)/sin(y)  /sin(y)/cos(y);
      //printf( "Sum = %.12f\n", sum );
      
      sum += EE*4*cos(x)/cos(y) /sin(y)/cos(y);
      //printf( "Sum = %.12f\n", sum );

      double diffed = 0.0;
      diffed -= EPI2*diff*diff/sin(x)/cos(y)/cos(y)/cos(y)/sin(y)/sin(y);
      //printf( "Suma = %.12f\n", diffed );
      diffed -= 2*(EF2)*diff/cos(y)/cos(y)/cos(y)/sin(x);
      //printf( "Suma = %.12f\n", diffed );
      diffed += EF*(1-cos(2*y))*diff/cos(y)/sin(y)/sin(y)/sin(y)/cos(x)/cos(y);
      //printf( "Suma = %.12f\n", diffed );

      sum += diffed;
      sum += 2*(- EK + EPI3)*cos(2*y)*diff/cos(y)/sin(y)/sin(y)/sin(y)/cos(x)/cos(y);
      //printf( "Sum = %.12f\n", sum );
      sum += 8*M_PI*middle;
      //printf( "Sum = %.12f\n", sum );

      //printf( "EPI = %.12e\n", EPI );
      sum += (EPI)*(cos(2*x)+cos(2*y))*(diff)/cos(y)/cos(x)/cos(y)/sin(y)/sin(y)/sin(y);
      //printf( "adding %.12e\n", (EPI)*(cos(2*x)+cos(2*y))*(diff)/cos(y)/sin(y)/sin(y)/sin(y)/cos(x)/cos(y) );
      //printf( "Sum = %.12e\n", sum );
      sum -= 8*M_PI*cos(2*x)/sin(2*y)/tan(2*y);            
      //printf( "Sum = %.12f\n", sum );
      
      return sum/8/M_PI;  
      //return (G - (1-x))/factor + M/sine2;
}

