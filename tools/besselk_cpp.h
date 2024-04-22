/* Header file containing the definitions of 
the besselk_cpp function */
#ifndef BESSELK_CPP_H
#define BESSELK_CPP_H
  #ifdef USE_BOOST
    #include "boost_exception_handlers_cpp.h"
    #include <boost/math/special_functions/bessel.hpp>
    using boost::math::cyl_bessel_k;
    using boost::math::policies::policy;
    using boost::math::policies::promote_float;
    using boost::math::policies::promote_double;
    using boost::math::policies::domain_error;
    using boost::math::policies::pole_error;
    using boost::math::policies::overflow_error;
    using boost::math::policies::underflow_error;
    using boost::math::policies::evaluation_error;
    using boost::math::policies::user_error;
    /* defining a policy for avoiding a SIGABRT in case the
       cyl_bessel_k function of boost does not converges */
    typedef policy <
      promote_float<false>,        //< promote float to double for increasing precision
      promote_double<false>,       //< promote double to long double for increasing precision
      domain_error<user_error>,    //< use user defined domain_error exception handler
      pole_error<user_error>,      //< use user defined pole_error exception handler
      overflow_error<user_error>,  //< use user defined overflow_error exception handler
      underflow_error<user_error>, //< use user defined underflow_error exception handler
      evaluation_error<user_error> //< use user defined evalauation_error exception handler
    > no_sigabrt_policy;
    BOOST_MATH_DECLARE_SPECIAL_FUNCTIONS(no_sigabrt_policy)
  #elif USE_STD_BESSELK
    #include <cmath>
    using std::cyl_bessel_k;
  #endif
  #ifdef __cplusplus
    extern "C" {
  #endif
  
  /* wrapper of the modified bessel function of the 
      2nd kind and with fractional order */
  //double besselk_cpp(double const &nu, double const &x);
  double besselk_cpp(double const &nu, double const &x);
  #ifdef __cplusplus
    }
  #endif
#endif

