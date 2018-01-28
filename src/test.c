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
  diskload_initialize();
  
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

  const double epsilon = 1e-9;

  SECTION( "Verify that we're close to the values from the MATLAB code" );
  ASSERT( fabs( ( -2.0449592462e+00 - u ) / u ) < epsilon );
  ASSERT( fabs( ( -2.0687091115e-01 - v ) / v ) < epsilon );
  ASSERT( fabs( ( 4.2902044944e-01 - g ) / g ) < epsilon );    

  SECTION( "Check error handling if we compute without first extrapolating" );
  e = diskload_truncated( alpha, Uncompensated, theta, w, 400000, love, &DefaultEarthModel,
                                        &u, &v, &g );
  // Fail because we haven't extrapolated the Love numbers
  ASSERT( e == E_LOVE_NUMBER_VECTOR_TOO_SHORT );

  SECTION( "Verify the truncated sum with N = 4M" );
  diskload_extrapolate_love_numbers(love, 4000000 );
  double uU, vU, gU;
  e = diskload_truncated( alpha, Uncompensated, theta, w, 4000000, love, &DefaultEarthModel,
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
  ASSERT( fabs( ( uU - u ) / u ) < 1e-6 );
  ASSERT( fabs( ( vU - v ) / v ) < 1e-3 );  
  ASSERT( fabs( ( gU - g ) / g ) < 1e-6 );

  printf( "%e %e %e\n", uU, vU, gU );
  printf( "%e %e %e\n", u, v, g );

  int i;
  for(i = 0; i<5;i++ ) {
    double range = 0.2;
    alpha = (double)rand()/(double)(RAND_MAX/range);
    theta = (double)rand()/(double)(RAND_MAX/range);    
    alpha += 0.2;
    theta += 0.5;    
    
    SECTION( "Test the globally compensated case via truncation with N = 400k" );
    double uC, vC, gC;
    
    e = diskload_truncated( alpha, Compensated, theta, w, 400000, love, &DefaultEarthModel,
                            &uC, &vC, &gC );
    ASSERT( e == E_SUCCESS );

    SECTION( "Compensated hypergeometric (N = 40k) near truncated (N = 400k)" );
    printf( BLUE "   (alpha = %e, theta = %e)\n\n" RESET, alpha, theta );
    
    e = diskload_core(alpha, Compensated, theta, w, diskload_core_G, diskload_core_H,
                      40000, love, &DefaultEarthModel,
                      &u, &v, &g );
    
    ASSERT( fabs(diskload_core_H(alpha,theta) - diskload_core_H_truncated(alpha,theta)) < 1e-4 );
    const double Degrees = M_PI / 180.0;
    printf( "%.12f,%.12f\n", alpha*Degrees,theta*Degrees );
    printf( "coreH     = %.12f\n", diskload_core_H(alpha*Degrees,theta*Degrees) );
    printf( "coreM     = %.12f\n", diskload_core_M(alpha*Degrees,theta*Degrees) );
    printf( "coreG     = %.12f\n", diskload_core_G(alpha*Degrees,theta*Degrees) );        
    //printf( "slow      = %.12f\n", diskload_core_H_slow(alpha*Degrees,theta*Degrees) );
    printf( "bad       = %.12f\n", diskload_core_H_bad(alpha*Degrees,theta*Degrees) );        
    printf( "truncated = %.12f\n", diskload_core_H_truncated(alpha*Degrees,theta*Degrees) );
    printf( "truncated - bad  = %.12f\n", diskload_core_H_bad(alpha*Degrees,theta*Degrees) - diskload_core_H_truncated(alpha*Degrees,theta*Degrees) );
    printf( "core - truncated = %.12f\n", diskload_core_H(alpha*Degrees,theta*Degrees) - diskload_core_H_truncated(alpha*Degrees,theta*Degrees) );    
    printf( "vC = %.12f\n", vC );
    printf( "v = %.12f\n", v );    

    ASSERT( e == E_SUCCESS );
    ASSERT( fabs( ( uC - u ) / u ) < 1e-4 );
    ASSERT( fabs( ( vC - v ) / v ) < 1e-4 );
    ASSERT( fabs( ( gC - g ) / g ) < 1e-4 );
  }

  SECTION( "Core G(x,y)" );
  
  ASSERT( fabs(diskload_core_G(0.1,0.2) - 0.025898198318164894) < 1e-6 );

  SECTION( "Core H(x,y)" );

  printf( "core h = %.12f\n", diskload_core_H(0.1,0.2) );
  printf( "core h - h' = %.12f\n", diskload_core_H(0.05,0.15) - diskload_core_H_truncated(0.05,0.15) );
  printf( "core h - h' = %.12f\n", diskload_core_H(0.001,0.0001) - diskload_core_H_truncated(0.001,0.0001) );

  {
  double x = 0.2; double y = -0.3;
  printf("M(%f,%f) = %.12f\n", x, y, diskload_core_M(x,y) );
  }


  {
    const double Degrees = M_PI / 180.0;    
    double x = 0.02*Degrees; double y = 0.01*Degrees;
    printf("H(%f,%f) = %.12f\n", x/Degrees, y/Degrees, diskload_core_H(x,y) );
    printf("Ht(%f,%f) = %.12f\n", x/Degrees, y/Degrees, diskload_core_H_truncated(x,y) );  
  }

  {
  double x = 0.2; double y = -0.5;
  printf("M(%f,%f) = %.12f\n", x, y, diskload_core_M(x,y) );
  }


    {
  double x = -0.2; double y = 0.6;
  printf("M(%f,%f) = %.12f\n", x, y, diskload_core_M(x,y) );
  }

  
  ASSERT( fabs(diskload_core_H(0.1,0.2) - 0.13628256838730426) < 1e-6 );
  ASSERT( fabs(diskload_core_H_truncated(0.1,0.2) - 0.13628256838730426) < 1e-6 );    
  
  ASSERT( fabs(diskload_hypergeometric_core_G(0.1,0.2) - 0.025898198318164894) < 1e-6 );

  /*
  e = diskload_core( alpha, Compensated, theta, w, diskload_hypergeometric_core, 40000, love, &DefaultEarthModel,
                                        &u, &v, &g );
  ASSERT( e == E_SUCCESS );

  */
  /* All done! */
  printf( CYAN "\n   %d tests; %d passed.\n\n" RESET, tests, successes );

  return 0;
}
