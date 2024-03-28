/* besselk_cpp is a wrapper used for
 importing the cyl_bessel_k function from boost
 or std and fixing its template variables.
 WARNING: the cyl_bessel_k function
          accepts only scalar variables */
#include "besselk_cpp.h"
/*Modified bessel function of the second kind
  with fractional order
  inputs:
    nu: (double) bessel function fractional order
    x:  (double) bessel function variable
	outputs:
    besselk_cpp: (double) modified bessel
                 function of the 2nd kind at x
                 and of order nu */
//double besselk_cpp(double const  &nu, double const &x){
double besselk_cpp(double const &nu, double const &x){
  #ifdef USE_BOOST
   // no_sigabrt_policy policy_to_use;
    return boost::math::cyl_bessel_k(nu,x,no_sigabrt_policy());
  #elif USE_STD_BESSELK
    return std::cyl_bessel_k(nu,x);
  #else
    return 0.e0;
  #endif
}
