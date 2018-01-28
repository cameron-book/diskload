int
main (void)
{
  load_love_numbers();
  
  clock_t begin = clock();

  double degrees = M_PI / 180.0;
  double x = 0.1*degrees, y;

  for( y=x/4.0; y<2*x; y=y+(x/200.0)) {
    double combined = uu_combined(x,y);
    double d = uu_spherical(x,y); 
    //double combined_bad = uu_combined_to_n(x,y,400);    
    //double kahan_big = uu_kahan_to_n(x,y, 4000000);
    //double kahan = uu_kahan(x,y);
    //double kahan = uu_kahan_to_n(x,y,400000);
    //double kahan_bad = uu_kahan_to_n(x,y, 4000);
    //double kahan_maybe = uu_kahan_to_n(x,y, 40000);
    //double kahan_really_bad = uu_kahan_to_n(x,y, 400);        
    //printf ("%e %e %e (err %e)\n", y/x, combined, kahan, fabs(combined-kahan) );
    //int n = uu_terms_needed(x, y, combined, 1e-8);
    //printf ("%e %e %e %e\n", y/x, combined, combined_bad, combined-combined_bad );
    //printf ("%e %e %e %e\n", y/x, kahan, d, (kahan - d)/kahan );
    //printf ("%e %e %e %e\n", y/x, kahan, combined, (kahan - combined)/kahan );
    printf ("%e %e %e %e\n", x, y, combined, d/combined );
  }

  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf( "%e seconds", time_spent );
  
  return 0;
}


/****************************************************************/
/* Naive version of Uu */

double uu_plain(double alpha, double theta) {
  double x = cos(theta);
  double y = cos(alpha);
    
  double p0 = 1;
  double p1 = x;
 
  double pp0 = 1;
  double pp1 = y;
  double q0 = 1-y;
  
  double result = h[0] * q0 * p0;
  
  double pp2, q1, p2;
  
  for( int n = 1; n < 4000000; n++ ) {
    pp2 = ((2*n+1)*y*pp1 - n*pp0)/(n+1);
    q1 = (pp0 - pp2)/(2*n+1);
    pp0 = pp1;
    pp1 = pp2;

    if (n <= 40000)
      result = result + h[n] * q1 * p1;
    else
      result = result + h[40000] * q1 * p1;

    //printf("%d %e\n", n, result);
    
    p2 = ((2*n+1)*x*p1 - n*p0)/(n+1);
    p0 = p1;
    p1 = p2;
  }
  
  return result;
}

double
spherical_integrand (double *k, size_t dim, void *params)
{
  (void)(dim);
  double phi = k[0];
  double theta = k[1];
  double tr = ((double*)params)[0];
  double result = (sin(theta)/(4.0 * M_PI)) * ( sqrt(2.0/(1.0-sin(theta)*cos(phi)*sin(tr) - cos(theta)*cos(tr))) );
  return result;
}

double spherical_core( double alpha, double tr ) {
  gsl_monte_function G = { &spherical_integrand, 2, &tr };
  
  gsl_monte_vegas_state *s = gsl_monte_vegas_alloc ( 2 );

  double xl[2] = { 0, 0 };
  double xu[2] = { 2*M_PI, alpha };
  
  const gsl_rng_type *T;
  gsl_rng *r;

  size_t calls = 10000;

  gsl_rng_env_setup ();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  
  double res, err;
    
  //gsl_monte_vegas_integrate (&G, xl, xu, 2, 10000, r, s, &res, &err);
  gsl_monte_vegas_integrate (&G, xl, xu, 2, calls, r, s, &res, &err);

  /*
  do
  {
    gsl_monte_vegas_integrate (&G, xl, xu, 2, calls/5, r, s, &res, &err);
    printf ("result = % .6f sigma = % .6f "
    "chisq/dof = %.1f\n", res, err, gsl_monte_vegas_chisq (s));
    
    
  }
  while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
  */
  //printf ("vegas final error = %.6f\n", err);
  //printf( "error = %g\n", err );
  
  gsl_monte_vegas_free (s);
  gsl_rng_free (r);

  return res;
}

double uu_spherical(double alpha, double theta) {
  double x = cos(theta);
  double y = cos(alpha);
#define N 40000
  double core = h[N] * spherical_core(alpha,theta);

  double p0 = 1;
  double p1 = x;
 
  double pp0 = 1;
  double pp1 = y;
  double q0 = 1-y;
  
  double result = (h[0] - h[N]) * q0 * p0;
  double c = 0.0; // compensation for low-order bits
  
  double pp2, q1, p2, input, t, z;
  
  for( int n = 1; n < N; n++ ) {
    pp2 = ((2*n+1)*y*pp1 - n*pp0)/(n+1);
    q1 = (pp0 - pp2)/(2*n+1);
    pp0 = pp1;
    pp1 = pp2;

    input = (h[n] - h[N]) * q1 * p1;

    z = input - c;
    t = result + z;
    c = (t - result) - z;
    result = t;

    p2 = ((2*n+1)*x*p1 - n*p0)/(n+1);
    p0 = p1;
    p1 = p2;
  }
  
  return core + result;
}


double uu_kahan(double alpha, double theta) {
  double x = cos(theta);
  double y = cos(alpha);
    
  double p0 = 1;
  double p1 = x;
 
  double pp0 = 1;
  double pp1 = y;
  double q0 = 1-y;
  
  double result = h[0] * q0 * p0;
  double c = 0.0; // compensation for low-order bits
  
  double pp2, q1, p2, input, t, z;
  
  for( int n = 1; n < 4000000; n++ ) {
    pp2 = ((2*n+1)*y*pp1 - n*pp0)/(n+1);
    q1 = (pp0 - pp2)/(2*n+1);
    pp0 = pp1;
    pp1 = pp2;

    input = h[40000] * q1 * p1;
    if (n <= 40000)
      input = h[n] * q1 * p1;

    z = input - c;
    t = result + z;
    c = (t - result) - z;
    result = t;

    p2 = ((2*n+1)*x*p1 - n*p0)/(n+1);
    p0 = p1;
    p1 = p2;
  }
  
  return result;
}


double uu_kahan_to_n(double alpha, double theta, int M) {
  double x = cos(theta);
  double y = cos(alpha);
    
  double p0 = 1;
  double p1 = x;
 
  double pp0 = 1;
  double pp1 = y;
  double q0 = 1-y;
  
  double result = h[0] * q0 * p0;
  double c = 0.0; // compensation for low-order bits
  
  double pp2, q1, p2, input, t, z;
  
  for( int n = 1; n < M; n++ ) {
    pp2 = ((2*n+1)*y*pp1 - n*pp0)/(n+1);
    q1 = (pp0 - pp2)/(2*n+1);
    pp0 = pp1;
    pp1 = pp2;

    input = h[40000] * q1 * p1;
    if (n <= 40000)
      input = h[n] * q1 * p1;

    z = input - c;
    t = result + z;
    c = (t - result) - z;
    result = t;

    p2 = ((2*n+1)*x*p1 - n*p0)/(n+1);
    p0 = p1;
    p1 = p2;
  }
  
  return result;
}



double uu_terms_needed(double alpha, double theta, double goal, double epsilon) {
  double x = cos(theta);
  double y = cos(alpha);
    
  double p0 = 1;
  double p1 = x;
 
  double pp0 = 1;
  double pp1 = y;
  double q0 = 1-y;
  
  double result = h[0] * q0 * p0;
  double c = 0.0; // compensation for low-order bits
  
  double pp2, q1, p2, input, t, z;
  
  for( int n = 1; n < 4000000; n++ ) {
    pp2 = ((2*n+1)*y*pp1 - n*pp0)/(n+1);
    q1 = (pp0 - pp2)/(2*n+1);
    pp0 = pp1;
    pp1 = pp2;

    input = h[40000] * q1 * p1;
    if (n <= 40000)
      input = h[n] * q1 * p1;

    z = input - c;
    t = result + z;
    c = (t - result) - z;
    result = t;

    p2 = ((2*n+1)*x*p1 - n*p0)/(n+1);
    p0 = p1;
    p1 = p2;
    
    if (fabs(result - goal) < epsilon)
      return n;
  }
  
  return -1;
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

double G(double x,double y) {
  gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc (10000);

  double lower_limit = x;	/* lower limit a */
  double upper_limit = 1;	/* upper limit b */
  double abs_error = 1.0e-6;	/* to avoid round-off problems */
  double rel_error = 1.0e-6;	/* the result will usually be much better */
  double result;		/* the result from the integration */
  double error;			/* the estimated error from the integration */

  double alpha = 1.0;		// parameter in integrand
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

double uu_combined(double alpha, double theta) {
  double x = cos(theta);
  double y = cos(alpha);
#define N 40000
  double core = h[N] * G(y,x);

  double p0 = 1;
  double p1 = x;
 
  double pp0 = 1;
  double pp1 = y;
  double q0 = 1-y;
  
  double result = (h[0] - h[N]) * q0 * p0;
  double c = 0.0; // compensation for low-order bits
  
  double pp2, q1, p2, input, t, z;
  
  for( int n = 1; n < N; n++ ) {
    pp2 = ((2*n+1)*y*pp1 - n*pp0)/(n+1);
    q1 = (pp0 - pp2)/(2*n+1);
    pp0 = pp1;
    pp1 = pp2;

    input = (h[n] - h[N]) * q1 * p1;

    z = input - c;
    t = result + z;
    c = (t - result) - z;
    result = t;

    p2 = ((2*n+1)*x*p1 - n*p0)/(n+1);
    p0 = p1;
    p1 = p2;
  }
  
  return core + result;
}

double uu_combined_to_n(double alpha, double theta, int M) {
  double x = cos(theta);
  double y = cos(alpha);
  double core = h[40000] * G(y,x);

  double p0 = 1;
  double p1 = x;
 
  double pp0 = 1;
  double pp1 = y;
  double q0 = 1-y;
  
  double result = (h[0] - h[40000]) * q0 * p0;
  double c = 0.0; // compensation for low-order bits
  
  double pp2, q1, p2, input, t, z;
  
  for( int n = 1; n < M; n++ ) {
    pp2 = ((2*n+1)*y*pp1 - n*pp0)/(n+1);
    q1 = (pp0 - pp2)/(2*n+1);
    pp0 = pp1;
    pp1 = pp2;

    input = (h[n] - h[40000]) * q1 * p1;

    z = input - c;
    t = result + z;
    c = (t - result) - z;
    result = t;

    p2 = ((2*n+1)*x*p1 - n*p0)/(n+1);
    p0 = p1;
    p1 = p2;
  }
  
  return core + result;
}


