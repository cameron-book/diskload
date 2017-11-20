#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "diskload.h"

#define RED     "\x1b[31m"
#define GREEN   "\x1b[32m"
#define YELLOW  "\x1b[33m"
#define BLUE    "\x1b[34m"
#define MAGENTA "\x1b[35m"
#define CYAN    "\x1b[36m"
#define RESET   "\x1b[0m"

#define CHECKMARK "\xE2\x9C\x93"
#define CROSS "x"
#define ASSERT(x) tests++; if(x) {successes++;  printf( " " GREEN CHECKMARK RESET " %s\n", #x ); } else printf( " " RED CROSS RESET " %s\n", #x )
#define SECTION(x) printf( BLUE "\n   %s\n\n" RESET, x );

int main( void ) {
  LoveNumbers* love = diskload_read_love_numbers("../REF_6371_loading_love_numbers_0_40000.txt");
  int tests = 0;
  int successes = 0;
  
  SECTION( "Verify that we loaded Love numbers from 0 to 40000" );
  ASSERT( love != NULL );
  ASSERT( love->degrees == 40001 );

  SECTION( "Check some inputs to be sure that they were loaded" )
  ASSERT( love->h != NULL );  
  ASSERT( love->l != NULL );
  ASSERT( love->k != NULL );  
  
  ASSERT( love->h[0] == -0.216 );
  ASSERT( love->l[0] ==  0.000 );
  ASSERT( love->k[0] ==  0.000 );

  ASSERT( love->h[3901] == -6.2144186090000 );
  ASSERT( love->l[3901] ==  0.0004851105826 );
  ASSERT( love->k[3901] == -0.0007828382211 );

  SECTION( "Evaluate a truncated sum (to 40,000) at a specific input" );
  double alpha = 0.1;
  double w = 1.0;
  double theta = 0.05794132044;
  double u, v, g;
  DiskLoadError e = diskload_truncated( alpha, Uncompensated, theta, w, 40000, love, &DefaultEarthModel,
                     &u, &v, &g );
  ASSERT (e == E_SUCCESS);

  const double epsilon = 10e-10;

  SECTION( "Verify that we're close to the values from the MATLAB code" );
  ASSERT( fabs( ( -2.0449592462e+00 - u ) / u ) < epsilon );
  ASSERT( fabs( ( -2.0687091115e-01 - v ) / v ) < epsilon );
  ASSERT( fabs( ( 4.2902044944e-01 - g ) / g ) < epsilon );    

  SECTION( "Check error handling if we compute without first extrapolating" );
  e = diskload_truncated( alpha, Uncompensated, theta, w, 400000, love, &DefaultEarthModel,
                                        &u, &v, &g );
  // Fail because we haven't extrapolated the Love numbers
  ASSERT( e == E_LOVE_NUMBER_VECTOR_TOO_SHORT );

  SECTION( "Verify the truncated sum with N = 400k" );
  diskload_extrapolate_love_numbers(love, 400000 );
  double uU, vU, gU;
  e = diskload_truncated( alpha, Uncompensated, theta, w, 400000, love, &DefaultEarthModel,
                                        &uU, &vU, &gU );
  ASSERT( e == E_SUCCESS );
  ASSERT( fabs( ( -2.0446039399e+00 - uU ) / uU ) < epsilon );
  ASSERT( fabs( ( -2.0707066586e-01 - vU ) / vU ) < epsilon );
  ASSERT( fabs( ( 4.2896329466e-01  - gU ) / gU ) < epsilon );

  SECTION( "Test the globally compensated case with N = 4000" );
  e = diskload_truncated( alpha, Compensated, theta, w, 4000, love, &DefaultEarthModel,
                                        &u, &v, &g );
  ASSERT( e == E_SUCCESS );
  ASSERT( fabs( ( -2.0997389244e+00 - u ) / u ) < epsilon );
  ASSERT( fabs( ( -2.1688909314e-01 - v ) / v ) < epsilon );
  ASSERT( fabs( (  4.3742883583e-01 - g ) / g ) < epsilon );

  SECTION( "Test the globally compensated case with N = 40k" );  
  e = diskload_truncated( alpha, Compensated, theta, w, 40000, love, &DefaultEarthModel,
                                        &u, &v, &g );
  ASSERT( e == E_SUCCESS );
  ASSERT( fabs( ( -2.0448711591e+00 - u ) / u ) < epsilon );
  ASSERT( fabs( ( -2.0687106866e-01 - v ) / v ) < epsilon );
  ASSERT( fabs( (  4.2860575533e-01 - g ) / g ) < epsilon );  
  
  SECTION( "Uncompensated hypergeometric core (N = 40k) approximates truncated (N = 400k)" );  
  e = diskload_hypergeometric( alpha, Uncompensated, theta, w, 40000, love, &DefaultEarthModel,
                                        &u, &v, &g );
  ASSERT( e == E_SUCCESS );
  ASSERT( fabs( ( uU - u ) / u ) < 10e-7 );
  ASSERT( fabs( ( vU - v ) / v ) < 10e-7 );  
  ASSERT( fabs( ( gU - g ) / g ) < 10e-7 );

  printf( "%e %e %e\n", uU, vU, gU );
  printf( "%e %e %e\n", u, v, g );

  int i;
  for(i = 0; i<5;i++ ) {
    double range = 0.2;
    alpha = (double)rand()/(double)(RAND_MAX/range);
    theta = (double)rand()/(double)(RAND_MAX/range);    
    
    SECTION( "Test the globally compensated case via truncation with N = 400k" );
    double uC, vC, gC;
    
    e = diskload_truncated( alpha, Compensated, theta, w, 400000, love, &DefaultEarthModel,
                            &uC, &vC, &gC );
    ASSERT( e == E_SUCCESS );

    SECTION( "Compensated hypergeometric (N = 40k) near truncated (N = 400k)" );
    printf( BLUE "   (alpha = %e, theta = %e)\n\n" RESET, alpha, theta );
    
    e = diskload_hypergeometric( alpha, Compensated, theta, w, 40000, love, &DefaultEarthModel,
                                 &u, &v, &g );
    
    //printf( "%e %e %e\n", uC, vC, gC );
    //printf( "%e %e %e\n", u, v, g );
    //printf( "%e %e %e\n", fabs( ( uC - u ) / u ), fabs( ( vU - v ) / v ), fabs( ( gC - g ) / g ) );    
    
    ASSERT( e == E_SUCCESS );
    ASSERT( fabs( ( uC - u ) / u ) < 10e-5 );
    ASSERT( fabs( ( vC - v ) / v ) < 10e-5 );
    ASSERT( fabs( ( gC - g ) / g ) < 10e-5 );
  }
  
  /*
  e = diskload_core( alpha, Compensated, theta, w, diskload_hypergeometric_core, 40000, love, &DefaultEarthModel,
                                        &u, &v, &g );
  ASSERT( e == E_SUCCESS );

  */
  /* All done! */
  printf( CYAN "\n   %d tests; %d passed.\n\n" RESET, tests, successes );

  return 0;
}
