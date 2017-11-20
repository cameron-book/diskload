#include <stdio.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "diskload.h"

const double Degrees = M_PI / 180.0;

const EarthModel DefaultEarthModel = { 6371.0, 9.8046961 };

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
  
  *u = LN->h[0] * sigma * p0;
  *v = 0.0;
  *g = sigma;
  
  // keeping track of low-order bits
  double Uc = 0.0, Vc = 0.0, Gc = 0.0;
  
  double dp1, q1, p2, input, t, z;

  for( int n = 1; n <= nmax; n++ ) {
    // Recursively compute the loading factor
    pp2 = ((2*n+1)*cosalpha*pp1 - n*pp0)/(n+1);

    if (icomp == Compensated) sigma = (pp0 - pp2)/(2*n+1)/(1+cosalpha);
    if (icomp == Uncompensated) sigma = (pp0 - pp2)/(2*n+1)/2.0;
    
    // Kahan summation for U
    z = (LN->h[n] * sigma * p1) - Uc;
    t = *u + z;
    Uc = (t - *u) - z;
    *u = t;

    // Compute the derivative of P_n(cos(x)) with respect to x
    dp1=(-n/sintheta)*(p0 - x*p1);

    // Kahan summation for V
    z = (LN->l[n] * sigma * dp1) - Vc;
    t = *v + z;
    Vc = (t - *v) - z;
    *v = t;
        
    // Kahan summation for the geoid term
    z = ((1.0 + LN->k[n]) * sigma * p1) - Gc;
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

/**
 * @brief compute disk load via a computational "core"
 *
 * @usage This function is exactly like diskload_truncated except for
 * the additional parameter `core` which is a pointer to a function
 * which computes the "core" contribution to the infinite sum
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

  double coreValue = (*coreG)(cos(alpha * Degrees),cos(theta * Degrees));
  double coreDerivative = (*coreH)(cos(alpha * Degrees),cos(theta * Degrees));
  
  double sigma;
  if (icomp == Compensated) sigma = 0.0;
  if (icomp == Uncompensated) sigma = (1.0-cosalpha)/2.0;

  double coreFactor;
  if (icomp == Compensated) coreFactor = 1.0+cosalpha;
  if (icomp == Uncompensated) coreFactor = 2.0;

  *u = LN->h[nmax] * coreValue / coreFactor + (LN->h[0] - LN->h[nmax]) * sigma * p0;
  double l_oo = (nmax + 1) * LN->l[nmax];
  *v = - sintheta * l_oo * coreDerivative / coreFactor;
  *g = coreValue / coreFactor;
 
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
    
    // Kahan summation for U
    z = ((LN->h[n] - LN->h[nmax]) * sigma * p1) - Uc;
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

double hyperg_z_GT1 (double a, double b, double c, double z) {
  double coef1,coef2;
  coef1=gsl_sf_gamma(c)*gsl_sf_gamma(b-a)*pow(1-z,-a)/(gsl_sf_gamma(b)*gsl_sf_gamma(c-a));
  coef2=gsl_sf_gamma(c)*gsl_sf_gamma(a-b)*pow(1-z,-b)/(gsl_sf_gamma(a)*gsl_sf_gamma(c-b));
  double result = coef1*gsl_sf_hyperg_2F1(a,c-b,a-b+1,1/(1-z))+coef2*gsl_sf_hyperg_2F1(b,c-a,b-a+1,1/(1-z));

  return result;
}

double f (double t, void *params)
{
  double y = *(double *) params;
  double z = (1-t*t)*(1-y*y)/((1-t*y)*(1-t*y));

  // BADBAD: this is broken
  if (fabs(z - 1.0) < 1e-10)
    z = 1.0 - 1e-10;
  
  if ((z < -1) || (z > 1))
    return hyperg_z_GT1(0.25, 0.75, 1.0, z)/sqrt(2-2*t*y);
  else
    return gsl_sf_hyperg_2F1(0.25, 0.75, 1.0, z)/sqrt(2-2*t*y);
}

double diskload_core_H(double x,double y) {
  int nmax = 40000;
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
  
  for( n=1; n <= nmax; n++ ) {
    p2 = ((2*n+1)*x*p1 - n*p0)/(n+1);    
    q1 = (p0 - p2)/(2*n+1);

    double dp1 = (y*pp1 - pp0) * n / (y*y - 1);

    result += q1 * dp1 / (n+1);

    p0 = p1;
    p1 = p2;
    
    pp2 = ((2*n+1)*y*pp1 - n*pp0)/(n+1);
    pp0 = pp1;
    pp1 = pp2;
  }
  return result;
}

double diskload_hypergeometric_core_G(double x,double y) {
  gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc (10000);

  double lower_limit = x;	/* lower limit a */
  double upper_limit = 1;	/* upper limit b */
  double abs_error = 1.0e-6;	/* to avoid round-off problems */
  double rel_error = 1.0e-6;	/* the result will usually be much better */
  double result;		/* the result from the integration */
  double error;			/* the estimated error from the integration */

  double expected = -4.0;	// exact answer

  gsl_function My_function;
  void *params_ptr = &y;
  My_function.function = &f;
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
			abs_error, rel_error, 1000, work_ptr, &result,
			&error);

  return result;
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
  return diskload_core(alpha,icomp,theta,w,diskload_hypergeometric_core_G, diskload_core_H, nmax,LN,earth,u,v,g);
}


/** @} */ // end of computations
