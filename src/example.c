#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "diskload.h"

#define CHECKMARK "\xE2\x9C\x93"
#define ASSERT(x) tests++; assert(x); printf( " " CHECKMARK " %s\n", #x )

int main( void ) {
  LoveNumbers* love = diskload_read_love_numbers("../REF_6371_loading_love_numbers_0_40000.txt");
  int tests = 0;
  
  // Verify that we loaded Love numbers from 0 to 40000...
  ASSERT( love != NULL );
  ASSERT( love->degrees == 40001 );

  // and check some random places to see if they were loaded
  // correctly.
  ASSERT( love->h != NULL );  
  ASSERT( love->l != NULL );
  ASSERT( love->k != NULL );  
  
  ASSERT( love->h[0] == -0.216 );
  ASSERT( love->l[0] ==  0.000 );
  ASSERT( love->k[0] ==  0.000 );

  ASSERT( love->h[3901] == -6.2144186090000 );
  ASSERT( love->l[3901] ==  0.0004851105826 );
  ASSERT( love->k[3901] == -0.0007828382211 );

  // At some random input...
  double alpha = 0.1;
  double w = 1.0;
  double theta = 0.05794132044;
  double u, v, g;
  
  // Evaluate via a truncated sum to N = 40000
  DiskLoadError e = diskload_truncated( alpha, Uncompensated, theta, w, 40000, love, &DefaultEarthModel,
                     &u, &v, &g );
  assert (e == E_SUCCESS);

  const double epsilon = 10e-11;
  
  // Verify that we're close to the values from the MATLAB code
  ASSERT( fabs( -2.0449592462e+00 - u ) / u < epsilon );
  ASSERT( fabs( -2.0687091115e-01 - v ) / v < epsilon );
  ASSERT( fabs( 4.2902044944e-01 - g ) / g < epsilon );    

  // Attempt to compute to N = 400,000...
  e = diskload_truncated( alpha, Uncompensated, theta, w, 400000, love, &DefaultEarthModel,
                                        &u, &v, &g );
  // and fail because we haven't extrapolated the Love numbers
  ASSERT( e == E_LOVE_NUMBER_VECTOR_TOO_SHORT );

  // Extrapolate and recompute...
  diskload_extrapolate_love_numbers(love, 400000 );
  e = diskload_truncated( alpha, Uncompensated, theta, w, 400000, love, &DefaultEarthModel,
                                        &u, &v, &g );
  
  // and check that we're close
  ASSERT( fabs( -2.0446039399e+00 - u ) / u < epsilon );
  ASSERT( fabs( -2.0707066586e-01 - v ) / v < epsilon );
  ASSERT( fabs( 4.2896329466e-01  - g ) / g < epsilon );

  ASSERT( false );
  
  printf( "\n   %d tests passed.\n\n", tests );
  
  return 0;
}
