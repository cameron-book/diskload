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
  LoveNumbers* love = diskload_read_love_numbers("../REF_6371_loading_love_numbers_0_40000.txt");
  int tests = 0;
  int successes = 0;
  
  assert( love != NULL );
  assert( love->degrees == 40001 );

  assert( love->h != NULL );  
  assert( love->l != NULL );
  assert( love->k != NULL );  

  diskload_extrapolate_love_numbers(love, 400000 );

  double range = 0.2;
  int N = 50;
  int i;
  double w = 1.0;
  double alpha, theta;
  double u, v, g;
  DiskLoadError e;

  {
  clock_t begin = clock();

  for(i = 0; i<N; i++ ) {
    alpha = (double)rand()/(double)(RAND_MAX/range);
    theta = (double)rand()/(double)(RAND_MAX/range);    
    
    e = diskload_truncated( alpha, Uncompensated, theta, w, 400000, love, &DefaultEarthModel,
                            &u, &v, &g );
  }

  clock_t end = clock();
  double time_spent = (double)(end - begin) * 10e3 / CLOCKS_PER_SEC / N;
  printf( "%.2f ms/truncated diskload\n", time_spent );
  }

  {
  clock_t begin = clock();

  for(i = 0; i<N; i++ ) {
    alpha = (double)rand()/(double)(RAND_MAX/range);
    theta = (double)rand()/(double)(RAND_MAX/range);    
    
    e = diskload_hypergeometric( alpha, Uncompensated, theta, w, 40000, love, &DefaultEarthModel,
                            &u, &v, &g );
  }

  clock_t end = clock();
  double time_spent = (double)(end - begin) * 10e3 / CLOCKS_PER_SEC / N;
  printf( "%.2f ms/hypergeometric diskload\n", time_spent );
  }
  
  return 0;
}
