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

  diskload_extrapolate_love_numbers(love, 400000 );

  double range = 0.2;
  int N = 50;
  int i;
  double w = 1.0;
  double alpha, theta;
  double u, v, g;
  DiskLoadError e;


#define TRUNCATE(cutoff,name) \
  { \
  clock_t begin = clock(); \
  \
  for(i = 0; i<N; i++ ) { \
    alpha = (double)rand()/(double)(RAND_MAX/range); \
    theta = (double)rand()/(double)(RAND_MAX/range); \
    \
    e = diskload_truncated( alpha, Uncompensated, theta, w, cutoff, love, &DefaultEarthModel, &u, &v, &g ); \
  } \
  clock_t end = clock(); \
  double time_spent = (double)(end - begin) * 10e3 / CLOCKS_PER_SEC / N; \
  printf( "%.2f ms/truncated diskload\n", time_spent ); \
  fprintf( table, "truncated at $N = %s$ & %.2f \\\\\n", name, time_spent ); \
  }

  TRUNCATE(40000,"40\\mathrm{k}")      
  TRUNCATE(200000,"200\\mathrm{k}")      
  TRUNCATE(400000,"400\\mathrm{k}")      
  
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
  fprintf( table, "core method & %.2f  \\\\\n", time_spent );
  }

  fprintf( table, "\\end{tabular}\n" );
fprintf( table, "\\caption{Time per disk load computation for various methods; these are computed by sampling %d random points in the range $0 \\leq \\alpha, \\theta \\leq %.2f$}\n", N, range );
  fprintf( table, "\\label{table:timings}\n" );
  fprintf( table, "\\end{table}\n" );

  return 0;
}
