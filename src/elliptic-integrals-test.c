#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "elliptic-integrals.h"

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

#define EPSILON 1e-8
//#define EQUALS(x,y) (printf("%f\n",x),fabs(x-y) < EPSILON)
#define EQUALS(x,y) (fabs(x-y) < EPSILON)

int main( void ) {
  const gsl_mode_t mode = GSL_PREC_DOUBLE;
  int tests = 0;
  int successes = 0;
  
  SECTION( "Verify that K(k) is defined for small k." );

  ASSERT( EQUALS(gsl_sf_ellint_Kcomp_extended( 0.1, mode ), 1.5747455615173559 ) );
  ASSERT( EQUALS(gsl_sf_ellint_Kcomp_extended( 0.9, mode ), 2.280549138422770204613 ) );

  SECTION( "Verify that K(k) is defined for large k." );
  
  ASSERT( EQUALS(gsl_sf_ellint_Kcomp_extended( 1.1, mode ), 2.1108400334328372 ) );
  ASSERT( EQUALS(gsl_sf_ellint_Kcomp_extended( 17, mode ), 0.09247987048715488822381078398 ) );

  SECTION( "Verify that F(phi,k) is defined for small k." );

  ASSERT( EQUALS(gsl_sf_ellint_F_extended( 0.2, 0.3, mode ), 0.20011923478531490626734858 ) );
  ASSERT( EQUALS(gsl_sf_ellint_F_extended( 0.9, 0.8, mode ), 0.9825430626671948124255033932023 ) );

  SECTION( "Verify that F(phi,k) is defined for not-too-large k." );

  ASSERT( EQUALS(gsl_sf_ellint_F_extended( 0.2, 1.1, mode ), 0.201635965422543380601632865996 ) );
  ASSERT( EQUALS(gsl_sf_ellint_F_extended( 0.1, 1.2, mode ), 0.100241081626030391233176040190719 ) );
  ASSERT( EQUALS(gsl_sf_ellint_F_extended( 0.3, 1.4, mode ), 0.309411047335846597526 ) );
  
  SECTION( "Verify that F(phi,k) is defined for large k." );  

  ASSERT( EQUALS(gsl_sf_ellint_F_extended( 0.8, 2.1, mode ),  0.79684355188395534476670777133 ) );
  ASSERT( EQUALS(gsl_sf_ellint_F_extended( 1.5, 100, mode ),   0.015708355989121522360264 ) );
  ASSERT( EQUALS(gsl_sf_ellint_F_extended( 1.5, 17, mode ), 0.0924798704871548882238107839899 ) );

  SECTION( "Verify that Π(k,n) is defined when k < 1 and n < 1." );  

  // Pcomp's n is the NEGATIVE of Mathematica's n
  ASSERT( EQUALS(gsl_sf_ellint_Pcomp_extended( 0.8, 0.7, mode ),  1.48069123246253303762443580 ) );
  ASSERT( EQUALS(gsl_sf_ellint_Pcomp_extended( 0.3, 0.4, mode ), 1.3563643538969762390660261840 ) );
  ASSERT( EQUALS(gsl_sf_ellint_Pcomp_extended( 0.9, 0.1, mode ), 2.15378685138752863562656029131 ) );
  ASSERT( EQUALS(gsl_sf_ellint_Pcomp_extended( 0.2, 0.6, mode ), 1.25303306759145565040719 ) );

  SECTION( "Verify that Π(k,n) is defined when k > 1 and n < 1." );  

  // Pcomp's n is the NEGATIVE of Mathematica's n
  ASSERT( EQUALS(gsl_sf_ellint_Pcomp_extended( 1.2, 0.6, mode ), 1.411442879361272150997977 ) );
  ASSERT( EQUALS(gsl_sf_ellint_Pcomp_extended( 2.1, 0.2, mode ), 0.7788119950567064340950532578956 ) );

  SECTION( "Verify that Π(k,n) is defined when k < 1 and n > 1." );  

  ASSERT( EQUALS(gsl_sf_ellint_Pcomp_extended( 0.8, 1.6, mode ),  1.167396979278389171140633433 ) );
  ASSERT( EQUALS(gsl_sf_ellint_Pcomp_extended( 0.2, 8.5, mode ),  0.5121783455185635577049018151 ) );

  SECTION( "Verify that Π(k,n) is defined when k > 1 and n > 1." );  

  ASSERT( EQUALS(gsl_sf_ellint_Pcomp_extended( 2.0, 1.00754, mode ), 0.750295441598245225632298674336121119129 ) );

  SECTION( "Verify that Π(k,n) is defined when k > 1 and n < -1." );  

  ASSERT( EQUALS(gsl_sf_ellint_Pcomp_extended( 2.0, -1.00754, mode ), -9.517077071840008254392703704619 ) );

  printf( CYAN "\n   %d tests; %d passed.\n\n" RESET, tests, successes );  
  
  return 0;
}
