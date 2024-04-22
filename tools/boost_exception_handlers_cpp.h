/*Special exception handlers for boost functions */
#ifndef BOOST_EXCEPTION_HANDLERS_CPP_H
#define BOOST_EXCEPTION_HANDLERS_CPP_H
#ifdef USE_BOOST
  #include <limits>
  #include <iostream>
  #include <iomanip>
  #include <boost/math/special_functions.hpp>
  #include <boost/format.hpp>
  using boost::io::group;
  using boost::format;
  namespace boost{ namespace math{ namespace policies{

  // define domain error handler
  template <class T>
  T user_domain_error(const char* function, const char* message, const T& val){
    //int prec = 2 + (std::numeric_limits<T>::digits*30103UL)/100000L;
    int prec = std::numeric_limits<T>::max_digits10;
    std::string msg("Error in function");
    if(function==0) function = "Unknown function with arguments of type %1%";
    if(message==0)  message = "Cause unknown with bad argument %1%";
    msg += (boost::format(function) % typeid(T).name()).str(); msg += ": \n";
    msg += (boost::format(message) % boost::io::group(std::setprecision(prec),val)).str();
    std::cerr << msg << std::endl;
    return std::numeric_limits<T>::quiet_NaN();
  }

  // define pole error handler
  template <class T>
  T user_pole_error(const char* function, const char* message, const T& val){
    return user_domain_error(function,message,val);
  }

  // define overlflow error handler
  template <class T>
   T user_overflow_error(const char* function, const char* message, const T& val){
    std::string msg("Error in function");
    if(function==0) function = "Unknown function with arguments of type %1%";
    if(message==0)  message = "Result of function is too larger to represent: overflow!";
    msg += (boost::format(function) % typeid(T).name()).str(); msg += ": \n";
    msg += message;
    std::cerr << msg << std::endl;
    return val;
  } 

  // define underflow error handler
  template <class T>
   T user_underflow_error(const char* function, const char* message, const T& val){
    std::string msg("Error in function");
    if(function==0) function = "Unknown function with arguments of type %1%";
    if(message==0)  message = "Result of function is too small to represent: underflow!";
    msg += (boost::format(function) % typeid(T).name()).str(); msg += ": \n";
    msg += message;
    std::cerr << msg << std::endl;
    return val;
  } 

  // define denormalised error handler
  template <class T>
   T user_denorm_error(const char* function, const char* message, const T& val){
    std::string msg("Error in function");
    if(function==0) function = "Unknown function with arguments of type %1%";
    if(message==0)  message = "Result of function is denormalised";
    msg += (boost::format(function) % typeid(T).name()).str(); msg += ": \n";
    msg += message;
    std::cerr << msg << std::endl;
    return val;
  } 

  // define evaluation error exception handler
  template <class T>
  T user_evaluation_error(const char* function, const char* message, const T& val){
    //int prec = 2 + (std::numeric_limits<T>::digits*30103UL)/100000L;
    int prec = std::numeric_limits<T>::max_digits10;
    std::string msg("Error in function");
    if(function==0) function = "Unknown function with arguments of type %1%";
    if(message==0)  message = "An internal evaluation error occurred with"
                              "the best value calculated so far of %1%";
    msg += (boost::format(function) % typeid(T).name()).str(); msg += ": \n";
    msg += (boost::format(message) % boost::io::group(std::setprecision(prec),val)).str();
    std::cerr << msg << std::endl;
    return std::numeric_limits<T>::quiet_NaN();
  }
}}}
#endif
#endif


