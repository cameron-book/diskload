#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include "diskload.h"

#define RED     "\x1b[31m"
#define GREEN   "\x1b[32m"
#define YELLOW  "\x1b[33m"
#define BLUE    "\x1b[34m"
#define MAGENTA "\x1b[35m"
#define CYAN    "\x1b[36m"
#define RESET   "\x1b[0m"

#define MICRO "\xc2\xb5"

int main( void ) {
  diskload_initialize();
  
  LoveNumbers* love = diskload_read_love_numbers("../REF_6371_loading_love_numbers_0_40000.txt");
  int tests = 0;
  int successes = 0;
  
  assert( love != NULL );
  assert( love->degrees == 40001 );

  assert( love->h != NULL );  
  assert( love->l != NULL );
  assert( love->k != NULL );  

  diskload_extrapolate_love_numbers(love, 4000000 );

  double range = 0.2;
  int N = 1000;
  int i;
  double w = 1.0;
  double theta;
  double u, v, g;
  DiskLoadError e;

  double alpha = 0.1;
  //double min_theta = alpha / 2.0;
  //double max_theta = alpha * 4.0;

  double min_theta = alpha * 0.99;
  double max_theta = alpha * 1.01;

  double base = exp( log( max_theta / min_theta ) / (N-1) );

  FILE* f = fopen("data.csv", "w");
  fprintf( f, "alpha theta alpha/theta U4K V4K G4K U40K V40K G40K U V G U400K V400K G400K U4M V4M G4M\n" );
  
  for(i = 0; i<N; i++ ) {
    theta = pow(base, i) * min_theta;

    fprintf(f, "%e %e %e ", alpha, theta, theta/alpha);
    
    e = diskload_truncated( alpha, Uncompensated, theta, w, 4000, love, &DefaultEarthModel,
                            &u, &v, &g );

    fprintf(f, "%e %e %e ", u, v, g);

    e = diskload_truncated( alpha, Uncompensated, theta, w, 40000, love, &DefaultEarthModel,
                            &u, &v, &g );

    fprintf(f, "%e %e %e ", u, v, g);    
    
    e = diskload_hypergeometric( alpha, Uncompensated, theta, w, 40000, love, &DefaultEarthModel,
                            &u, &v, &g );

    fprintf(f, "%e %e %e ", u, v, g);
    
    e = diskload_truncated( alpha, Uncompensated, theta, w, 400000, love, &DefaultEarthModel,
                            &u, &v, &g );

    fprintf(f, "%e %e %e ", u, v, g);

    e = diskload_truncated( alpha, Uncompensated, theta, w, 4000000, love, &DefaultEarthModel,
                            &u, &v, &g );

    fprintf(f, "%e %e %e ", u, v, g);    
    
    fprintf(f, "\n");
  }

  
  return 0;
}
