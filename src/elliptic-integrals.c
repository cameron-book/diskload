#include <math.h>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "elliptic-integrals.h"
#include "float.h"

double gsl_sf_ellint_Kcomp_extended (double k, gsl_mode_t mode) {
  if (k < 1) return gsl_sf_ellint_Kcomp (k, mode);

  return gsl_sf_ellint_Kcomp (1/k, mode) / k;
}

#define cot(x) (1.0/tan(x))
#define acot(x) (atan(1.0/(x)))
#define csc(x) (1.0/sin(x))
#define square(x) ((x)*(x))

// following https://github.com/moiseevigor/elliptic/blob/master/elliptic123.m
// so this code is likewise licensed under the GPL
gsl_complex gsl_sf_ellint_Fz (gsl_complex u, double k, gsl_mode_t mode) {
  double phi = GSL_REAL(u);
  double psi = GSL_IMAG(u);
  double m = k*k;

  if (k > 1) {
    gsl_complex sine = gsl_complex_mul_real( gsl_complex_sin(u), k );
    return gsl_complex_div_real( gsl_sf_ellint_Fz(gsl_complex_arcsin(sine),1.0/k,mode), k);
  }
  
  // If real...
  if (fabs(GSL_IMAG(u)) < 1e-15) {
    return gsl_complex_rect( gsl_sf_ellint_F (GSL_REAL(u), k, mode), 0 );
  }

  // This avoids the singularity of cot
  if (fabs(phi) < 1e-15)
    phi = 1e-15;

  // finding the roots of the equation X^2 - (cot(phi)^2+m*sinh(psi)^2*csc(phi)^2-1+m)X - (1-m)*cot(phi)^2 = 0
  double b = -(square(cot(phi)) + m*square(sinh(psi)*csc(phi))-1+m);
  double c = -(1-m)*square(cot(phi));

  double x = -b/2 + sqrt(square(b)/4-c);
  if (x<0) x = -b/2 - sqrt(square(b)/4-c);
  
  double lambda = acot( sqrt(x) );
  double mu = atan( sqrt(1/m*(square(tan(phi)*cot(lambda)) - 1)) );

  // change of variables taking into account periodicity ceil to the right
  double sign = (((int)(floor(phi/M_PI*2)) % 2) == 0) ? (1) : (-1);
  lambda = sign*lambda + M_PI*ceil(phi/M_PI-0.5+DBL_EPSILON);
  mu = copysign( mu, psi );
  
  return gsl_complex_rect (gsl_sf_ellint_F( lambda, k, mode ),
                           gsl_sf_ellint_F( mu,  sqrt(1-k*k), mode));
}

double gsl_sf_ellint_F_extended (double phi, double k, gsl_mode_t mode) {
  if (k < 1) return gsl_sf_ellint_F (phi, k, mode);

  double sine = k*sin(phi);
  if (sine <= 1) return gsl_sf_ellint_F(asin(sine),1.0/k,mode)/k;
  return GSL_REAL(gsl_sf_ellint_Fz(gsl_complex_arcsin_real(sine),1.0/k,mode))/k;
}

gsl_complex gsl_sf_ellint_Pcomp_z (double k, double n, gsl_mode_t mode) {
  printf( "computing PI(k=%f,n=%f)\n", k, n );

  if (fabs(n) > 10e100) return gsl_complex_rect( 0, 0 );
  if (fabs(k) > 10e100) return gsl_complex_rect( 0, 0 );  
  
  if ((k == 0) && (n == 0)) {
    return gsl_complex_rect( M_PI/2, 0 );
  }

  if (k < 0) {
    return gsl_sf_ellint_Pcomp_z (-k, n, mode);
  }
  
  if (k < 1) {
    printf( "PI = %f\n", gsl_sf_ellint_Pcomp (k, n, mode) );

    return gsl_complex_rect( gsl_sf_ellint_Pcomp (k, n, mode), 0 );
  }

  gsl_complex lhs = gsl_complex_div_real(gsl_sf_ellint_Pcomp_z( 1/k, n/k/k, mode ), k );
  gsl_complex rhs = gsl_sf_ellint_Pcomp_z(sqrt(1-k*k), -(-n * (1 - k*k))/(-n - k*k), mode);
  gsl_complex factor = gsl_complex_rect( 0, -(k*k/(k*k+n)) );

  gsl_complex result  = gsl_complex_add(lhs, gsl_complex_mul(factor, rhs));

  printf( "PI(k=%f,n=%f) = %f + I %f\n", k, n, GSL_REAL(result), GSL_IMAG(result) );
  printf( "EllipticPi[%f,%f] - (%f + I*%f)\n", n, k*k, GSL_REAL(result), GSL_IMAG(result) );  
  
  return gsl_complex_add(lhs, gsl_complex_mul(factor, rhs));
}

// RC1y and gsl_sf_ellint_RJ_z are from mpmath;
// Copyright (c) 2005-2017 Fredrik Johansson and mpmath contributors
gsl_complex RC1y(gsl_complex y) {
  // This computes RC(1,y)
  gsl_complex v;
  
  if ((GSL_REAL(y) == 0) && (GSL_IMAG(y) == 0)) {
    return gsl_complex_rect( INFINITY, 0 );
  }
  if ((GSL_REAL(y) == 1) && (GSL_IMAG(y) == 0)) {
    return y;
  }

  gsl_complex x = gsl_complex_rect( 1, 0 );
    
  double extraprec = 2*fmax(0,-gsl_complex_abs(gsl_complex_sub(x,y))+gsl_complex_abs(x));
  
  if (GSL_IMAG(y) == 0) {
    double rx = GSL_REAL(x);
    double ry = GSL_REAL(y);
    gsl_complex a = gsl_complex_sqrt_real(rx/ry);
    if (rx < ry) {
      double b = sqrt(ry-rx);
      v = gsl_complex_div_real(gsl_complex_arccos(a),b);
    } else {
      double b = sqrt(rx-ry);
      v = gsl_complex_div_real(gsl_complex_arccosh(a),b);      
    }
  } else {
    gsl_complex sx = gsl_complex_sqrt(x);
    gsl_complex sy = gsl_complex_sqrt(y);    
    v = gsl_complex_arccos(gsl_complex_div(sx,sy));
    gsl_complex one = gsl_complex_rect(1,0);
    v = gsl_complex_div( v, gsl_complex_mul( gsl_complex_sqrt(gsl_complex_sub(one, gsl_complex_div(x,y))), sy ) );
  }
  return v;
} 
             
gsl_complex gsl_sf_ellint_RJ_z (double x, double y, double z, double p, gsl_mode_t mode) {
  // this is a tolerance and it should depend on mode
  double r = 1e-8;
  
  if (!(isnormal(x) && isnormal(y) && isnormal(z) && isnormal(p))) {
    if (isnan(x) || isnan(y) || isnan(z) || isnan(p)) {
      return gsl_complex_rect( x*y*z, 0 );
    }
    if (isinf(x) || isinf(y) || isinf(z) || isinf(p)) {
      return gsl_complex_rect( 0, 0 );
    }
  }

  if (!p) {
    return gsl_complex_rect( INFINITY, 0 );
  }
  
  gsl_complex xm=gsl_complex_rect(x,0),
              ym=gsl_complex_rect(y,0),
              zm=gsl_complex_rect(z,0),
              pm=gsl_complex_rect(p,0);
  
  double A0 = (x + y + z + 2*p)/5;
  gsl_complex Am = gsl_complex_rect( A0, 0 );
  double delta = (p-x)*(p-y)*(p-z);

  double Q = pow(0.25*r, -1.0/6.0) * fmax(fmax(fmax(fabs(A0-x),fabs(A0-y)),fabs(A0-z)),fabs(A0-p));
  int m = 0;
  double g = 0.25;
  double pow4 = 1.0;
  gsl_complex S = gsl_complex_rect( 0, 0 );
  for(;;) {
    gsl_complex sx = gsl_complex_sqrt(xm);
    gsl_complex sy = gsl_complex_sqrt(ym);
    gsl_complex sz = gsl_complex_sqrt(zm);
    gsl_complex sp = gsl_complex_sqrt(pm);
    gsl_complex lm = gsl_complex_add(gsl_complex_add(gsl_complex_mul(sx,sy),
                                                     gsl_complex_mul(sx,sz)),
                                     gsl_complex_mul(sy,sz));
    gsl_complex Am1 = gsl_complex_mul_real(gsl_complex_add(Am,lm),g);

    xm = gsl_complex_mul_real(gsl_complex_add(xm,lm),g);
    ym = gsl_complex_mul_real(gsl_complex_add(ym,lm),g);
    zm = gsl_complex_mul_real(gsl_complex_add(zm,lm),g);
    pm = gsl_complex_mul_real(gsl_complex_add(pm,lm),g);    

    gsl_complex dm = gsl_complex_mul(gsl_complex_mul(gsl_complex_add(sp,sx),
                                                     gsl_complex_add(sp,sy)),
                                     gsl_complex_add(sp,sz));
    
    gsl_complex em = gsl_complex_div( gsl_complex_rect( delta * pow(4, -3*m), 0 ),
                                      gsl_complex_mul(dm,dm));
    
    if (pow4 * Q < gsl_complex_abs(Am))
      break;
    
    gsl_complex T = gsl_complex_div(gsl_complex_mul_real(RC1y(gsl_complex_add_real(em,1)), pow4), dm);
    
    S = gsl_complex_add(S, T);
    pow4 *= g;
    m++;
    Am = Am1;
  }
  
  gsl_complex t = gsl_complex_div( gsl_complex_rect( ldexp(1,-2*m), 0 ), Am );
  gsl_complex X = gsl_complex_mul_real( t, (A0-x) );
  gsl_complex Y = gsl_complex_mul_real( t, (A0-y) );
  gsl_complex Z = gsl_complex_mul_real( t, (A0-z) );
  gsl_complex P = gsl_complex_mul_real(gsl_complex_add(gsl_complex_add(X,Y),Z),-0.5);
  
  gsl_complex E2 = gsl_complex_add( gsl_complex_add( gsl_complex_add( gsl_complex_mul(X,Y) , gsl_complex_mul(X,Z) ), gsl_complex_mul(Y,Z) ), gsl_complex_mul_real(gsl_complex_mul(P,P),-3) );
  
  gsl_complex E3 = gsl_complex_add(gsl_complex_add(gsl_complex_mul(gsl_complex_mul(X,Y),Z),
                                                   gsl_complex_mul_real(gsl_complex_mul(E2,P),2)),
                                   gsl_complex_mul_real(gsl_complex_mul(gsl_complex_mul(P,P),P),4));
  
  gsl_complex E4 = gsl_complex_mul(P,
                                   gsl_complex_add(gsl_complex_add(gsl_complex_mul(gsl_complex_mul(gsl_complex_mul_real(X,2),Y),Z),
                                                                   gsl_complex_mul(E2,P)),
                                                   gsl_complex_mul_real(gsl_complex_mul(gsl_complex_mul(P,P),P),3))
                                   );

  gsl_complex E5 = gsl_complex_mul(gsl_complex_mul(gsl_complex_mul(gsl_complex_mul(X,Y),Z),P),P);
  
  P = gsl_complex_rect(24024,0);
  P = gsl_complex_sub( P, gsl_complex_mul_real(E2,5148) );
  P = gsl_complex_add( P, gsl_complex_mul_real(gsl_complex_mul(E2,E2),2457) );
  P = gsl_complex_add( P, gsl_complex_mul_real(E3,4004) );
  P = gsl_complex_sub( P, gsl_complex_mul_real(gsl_complex_mul(E2,E3),4158) );
  P = gsl_complex_sub( P, gsl_complex_mul_real(E4,3276) );
  P = gsl_complex_sub( P, gsl_complex_mul_real(gsl_complex_mul(E2,E3),4158) );  
  P = gsl_complex_add( P, gsl_complex_mul_real(E5,2772) );
  
  Q = 24024;
  
  gsl_complex v1 = gsl_complex_div_real(
      gsl_complex_mul(gsl_complex_mul_real(gsl_complex_pow_real(Am, -1.5), pow(g,m)),P),
      Q);
  
  gsl_complex v2 = gsl_complex_mul_real(S,6);
  
  return gsl_complex_add(v1,v2);
}
                
double gsl_sf_ellint_Pcomp_extended (double k, double n, gsl_mode_t mode) {
   const double y = 1.0 - k*k;
   n = -n;
   return   gsl_sf_ellint_Kcomp_extended(k,mode)
       + (n/3.0) * GSL_REAL(gsl_sf_ellint_RJ_z(0.0, y, 1.0, 1.0 - n, mode));
}
