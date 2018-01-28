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

#ifndef POSTFIX
#define POSTFIX ""
#endif

#ifndef TITLE
#define TITLE ""
#endif

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

  FILE* table = fopen("table-timings" POSTFIX ".tex", "w");
  fprintf( table, "\\begin{table}[ht]\n" );
  fprintf( table, "\\centering\n" );  
  fprintf( table, "\\begin{tabular}{lS[table-format=3.2]}\n" );
  fprintf( table, "method & {time (ms)} \\\\\n" );
  fprintf( table, "\\hline\n" );  
 
  double range = 0.2;
  int N = 700;
  int i;
  double w = 1.0;
  double alpha, theta;
  double u, v, g;
  DiskLoadError e;

  FILE* f = fopen("timings" POSTFIX ".csv", "w");
  fprintf( f, "cutoff milliseconds\n" );

  int cutoff;

  // BADBAD: for( cutoff=10000; cutoff<=400000; cutoff += 10000 ) {
  for( cutoff=10000; cutoff<=80000; cutoff += 10000 ) {
    clock_t begin = clock();
    
    for(i = 0; i<N; i++ ) {
      alpha = (double)rand()/(double)(RAND_MAX/range); 
      theta = (double)rand()/(double)(RAND_MAX/range); 
      
      e = diskload_truncated( alpha, Uncompensated, theta, w, cutoff, love, &DefaultEarthModel, &u, &v, &g ); 
    } 
    clock_t end = clock(); 
    double time_spent = (double)(end - begin) * 10e3 / CLOCKS_PER_SEC / N; 
    fprintf( f, "%d %e\n", cutoff, time_spent );

    char *name = NULL;
    if (cutoff == 40000) name = "40\\mathrm{k}";
    if (cutoff == 100000) name = "100\\mathrm{k}";
    if (cutoff == 200000) name = "200\\mathrm{k}";
    if (cutoff == 400000) name = "400\\mathrm{k}";            
    if (name != NULL) {  
      printf( "%.2f ms/truncated diskload\n", time_spent );             \
      fprintf( table, "truncated at $N = %s$ & %.2f \\\\\n", name, time_spent ); \
    }
  }

  double core_time;
  
  {
    clock_t begin = clock();
    
    for(i = 0; i<N; i++ ) {
      alpha = (double)rand()/(double)(RAND_MAX/range);
      theta = (double)rand()/(double)(RAND_MAX/range);    
      
      e = diskload_hypergeometric( alpha, Uncompensated, theta, w, 40000, love, &DefaultEarthModel,
                                   &u, &v, &g );
    }
    clock_t end = clock();
    core_time = (double)(end - begin) * 10e3 / CLOCKS_PER_SEC / N;
    printf( "%.2f ms/hypergeometric diskload\n", core_time );
    fprintf( table, "core method & %.2f  \\\\\n", core_time );
    
  }

  fprintf( table, "\\end{tabular}\n" );
  fprintf( table, "\\caption{Time per disk load computation" TITLE " for various methods; these are computed by sampling %d random points in the range $0 \\leq \\alpha, \\theta \\leq %.2f$}\n", N, range );
  fprintf( table, "\\label{table:timings" POSTFIX "}\n" );
  fprintf( table, "\\end{table}\n" );
  
  FILE* gp = fopen("timings" POSTFIX ".gp", "w");
  fprintf( gp, "set terminal tikz\n" 
           "set output 'timings" POSTFIX ".tex'\n"
           "set xlabel 'terms in truncated series'\n"
           "set ylabel 'time (milliseconds)'\n"
           "set key off\n"
           "set style line 1 dt 1 lc 'black'\n"
           "set arrow from graph 0, first %f to graph 1, first %f nohead lc rgb 'gray'\n"
           "set style textbox opaque noborder\n"
           "set title 'Speed of core method and series truncation" TITLE "'\n"
           "set label 'core method' at graph 0.7, first %f boxed front\n"
           "plot 'timings" POSTFIX ".csv' using (column('cutoff')):(column('milliseconds')) with lines ls 1 title 'series truncation'\n", core_time, core_time, core_time );
  
  
  return 0;
}
