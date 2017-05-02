/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 4.0.0
 *
 * This file is not intended to be easily readable and contains a number of
 * coding conventions designed to improve portability and efficiency. Do not make
 * changes to this file unless you know what you are doing--modify the SWIG
 * interface file instead.
 * ----------------------------------------------------------------------------- */

#ifdef __cplusplus
/* SwigValueWrapper is described in swig.swg */
template<typename T> class SwigValueWrapper {
  struct SwigMovePointer {
    T *ptr;
    SwigMovePointer(T *p) : ptr(p) { }
    ~SwigMovePointer() { delete ptr; }
    SwigMovePointer& operator=(SwigMovePointer& rhs) { T* oldptr = ptr; ptr = 0; delete oldptr; ptr = rhs.ptr; rhs.ptr = 0; return *this; }
  } pointer;
  SwigValueWrapper& operator=(const SwigValueWrapper<T>& rhs);
  SwigValueWrapper(const SwigValueWrapper<T>& rhs);
public:
  SwigValueWrapper() : pointer(0) { }
  SwigValueWrapper& operator=(const T& t) { SwigMovePointer tmp(new T(t)); pointer = tmp; return *this; }
  operator T&() const { return *pointer.ptr; }
  T *operator&() { return pointer.ptr; }
};

template <typename T> T SwigValueInit() {
  return T();
}
#endif

/* -----------------------------------------------------------------------------
 *  This section contains generic SWIG labels for method/variable
 *  declarations/attributes, and other compiler dependent labels.
 * ----------------------------------------------------------------------------- */

/* template workaround for compilers that cannot correctly implement the C++ standard */
#ifndef SWIGTEMPLATEDISAMBIGUATOR
# if defined(__SUNPRO_CC) && (__SUNPRO_CC <= 0x560)
#  define SWIGTEMPLATEDISAMBIGUATOR template
# elif defined(__HP_aCC)
/* Needed even with `aCC -AA' when `aCC -V' reports HP ANSI C++ B3910B A.03.55 */
/* If we find a maximum version that requires this, the test would be __HP_aCC <= 35500 for A.03.55 */
#  define SWIGTEMPLATEDISAMBIGUATOR template
# else
#  define SWIGTEMPLATEDISAMBIGUATOR
# endif
#endif

/* inline attribute */
#ifndef SWIGINLINE
# if defined(__cplusplus) || (defined(__GNUC__) && !defined(__STRICT_ANSI__))
#   define SWIGINLINE inline
# else
#   define SWIGINLINE
# endif
#endif

/* attribute recognised by some compilers to avoid 'unused' warnings */
#ifndef SWIGUNUSED
# if defined(__GNUC__)
#   if !(defined(__cplusplus)) || (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4))
#     define SWIGUNUSED __attribute__ ((__unused__))
#   else
#     define SWIGUNUSED
#   endif
# elif defined(__ICC)
#   define SWIGUNUSED __attribute__ ((__unused__))
# else
#   define SWIGUNUSED
# endif
#endif

#ifndef SWIG_MSC_UNSUPPRESS_4505
# if defined(_MSC_VER)
#   pragma warning(disable : 4505) /* unreferenced local function has been removed */
# endif
#endif

#ifndef SWIGUNUSEDPARM
# ifdef __cplusplus
#   define SWIGUNUSEDPARM(p)
# else
#   define SWIGUNUSEDPARM(p) p SWIGUNUSED
# endif
#endif

/* internal SWIG method */
#ifndef SWIGINTERN
# define SWIGINTERN static SWIGUNUSED
#endif

/* internal inline SWIG method */
#ifndef SWIGINTERNINLINE
# define SWIGINTERNINLINE SWIGINTERN SWIGINLINE
#endif

/* exporting methods */
#if defined(__GNUC__)
#  if (__GNUC__ >= 4) || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4)
#    ifndef GCC_HASCLASSVISIBILITY
#      define GCC_HASCLASSVISIBILITY
#    endif
#  endif
#endif

#ifndef SWIGEXPORT
# if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
#   if defined(STATIC_LINKED)
#     define SWIGEXPORT
#   else
#     define SWIGEXPORT __declspec(dllexport)
#   endif
# else
#   if defined(__GNUC__) && defined(GCC_HASCLASSVISIBILITY)
#     define SWIGEXPORT __attribute__ ((visibility("default")))
#   else
#     define SWIGEXPORT
#   endif
# endif
#endif

/* calling conventions for Windows */
#ifndef SWIGSTDCALL
# if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
#   define SWIGSTDCALL __stdcall
# else
#   define SWIGSTDCALL
# endif
#endif

/* Deal with Microsoft's attempt at deprecating C standard runtime functions */
#if !defined(SWIG_NO_CRT_SECURE_NO_DEPRECATE) && defined(_MSC_VER) && !defined(_CRT_SECURE_NO_DEPRECATE)
# define _CRT_SECURE_NO_DEPRECATE
#endif

/* Deal with Microsoft's attempt at deprecating methods in the standard C++ library */
#if !defined(SWIG_NO_SCL_SECURE_NO_DEPRECATE) && defined(_MSC_VER) && !defined(_SCL_SECURE_NO_DEPRECATE)
# define _SCL_SECURE_NO_DEPRECATE
#endif

/* Deal with Apple's deprecated 'AssertMacros.h' from Carbon-framework */
#if defined(__APPLE__) && !defined(__ASSERT_MACROS_DEFINE_VERSIONS_WITHOUT_UNDERSCORES)
# define __ASSERT_MACROS_DEFINE_VERSIONS_WITHOUT_UNDERSCORES 0
#endif

/* Intel's compiler complains if a variable which was never initialised is
 * cast to void, which is a common idiom which we use to indicate that we
 * are aware a variable isn't used.  So we just silence that warning.
 * See: https://github.com/swig/swig/issues/192 for more discussion.
 */
#ifdef __INTEL_COMPILER
# pragma warning disable 592
#endif


#include <string>
#include <algorithm>
#include <stdexcept>


#include <algorithm>


#include <stdexcept>


#include <string>



// Fill a Fortran string from a std::string; with whitespace after
void std_string_copyout(const std::string& str, char* s, size_t count)
{
    if (str.size() > count)
        throw std::range_error("string size too small");

    s = std::copy(str.begin(), str.end(), s);
    std::fill_n(s, count - str.size(), ' ');
}



#include "trilinos_handle.hpp"

#ifdef __cplusplus
extern "C" {
#endif
SWIGEXPORT void* swigc_new_TrilinosHandle() {
  void* fresult = 0 ;
  ForTrilinos::TrilinosHandle *result = 0 ;
  
  result = (ForTrilinos::TrilinosHandle *)new ForTrilinos::TrilinosHandle();
  fresult = result; 
  return fresult;
}


SWIGEXPORT void swigc_TrilinosHandle_init__SWIG_0(void* farg1) {
  ForTrilinos::TrilinosHandle *arg1 = (ForTrilinos::TrilinosHandle *) 0 ;
  
  arg1 = (ForTrilinos::TrilinosHandle *)(farg1); 
  (arg1)->init();
}


SWIGEXPORT void swigc_TrilinosHandle_init__SWIG_1(void* farg1, int* farg2) {
  ForTrilinos::TrilinosHandle *arg1 = (ForTrilinos::TrilinosHandle *) 0 ;
  MPI_Comm arg2 ;
  
  arg1 = (ForTrilinos::TrilinosHandle *)(farg1); 
  
  arg2 = (MPI_Comm)(MPI_Comm_f2c(*(MPI_Fint *)(farg2)));
  
  (arg1)->init(arg2);
}


SWIGEXPORT void swigc_TrilinosHandle_setup_matrix(void* farg1, int* farg2, int* farg3, int* farg4, int* farg5, int* farg6, double* farg7) {
  ForTrilinos::TrilinosHandle *arg1 = (ForTrilinos::TrilinosHandle *) 0 ;
  int arg2 ;
  int *arg3 = (int *) 0 ;
  int *arg4 = (int *) 0 ;
  int arg5 ;
  int *arg6 = (int *) 0 ;
  double *arg7 = (double *) 0 ;
  
  arg1 = (ForTrilinos::TrilinosHandle *)(farg1); 
  arg2 = *farg2;
  arg3 = farg3;
  arg4 = farg4;
  arg5 = *farg5;
  arg6 = farg6;
  arg7 = farg7;
  (arg1)->setup_matrix(arg2,(int const *)arg3,(int const *)arg4,arg5,(int const *)arg6,(double const *)arg7);
}


SWIGEXPORT void swigc_TrilinosHandle_setup_operator(void* farg1, int* farg2, int* farg3, void* farg4) {
  ForTrilinos::TrilinosHandle *arg1 = (ForTrilinos::TrilinosHandle *) 0 ;
  int arg2 ;
  int *arg3 = (int *) 0 ;
  ForTrilinos::TrilinosHandle::OperatorCallback arg4 = (ForTrilinos::TrilinosHandle::OperatorCallback) 0 ;
  
  arg1 = (ForTrilinos::TrilinosHandle *)(farg1); 
  arg2 = *farg2;
  arg3 = farg3;
  arg4 = (ForTrilinos::TrilinosHandle::OperatorCallback)(farg4); 
  (arg1)->setup_operator(arg2,(int const *)arg3,arg4);
}


SWIGEXPORT void swigc_TrilinosHandle_setup_solver(void* farg1, void * farg2) {
  ForTrilinos::TrilinosHandle *arg1 = (ForTrilinos::TrilinosHandle *) 0 ;
  Teuchos::RCP< Teuchos::ParameterList > *arg2 = 0 ;
  Teuchos::RCP< Teuchos::ParameterList > tempnull2 ;
  
  arg1 = (ForTrilinos::TrilinosHandle *)(farg1); 
  arg2 = farg2 ? (Teuchos::RCP< Teuchos::ParameterList > *)farg2 : &tempnull2;
  (arg1)->setup_solver((Teuchos::RCP< Teuchos::ParameterList > const &)*arg2);
}


SWIGEXPORT void swigc_TrilinosHandle_solve(void* farg1, int* farg2, double* farg3, double* farg4) {
  ForTrilinos::TrilinosHandle *arg1 = (ForTrilinos::TrilinosHandle *) 0 ;
  int arg2 ;
  double *arg3 = (double *) 0 ;
  double *arg4 = (double *) 0 ;
  
  arg1 = (ForTrilinos::TrilinosHandle *)(farg1); 
  arg2 = *farg2;
  arg3 = farg3;
  arg4 = farg4;
  ((ForTrilinos::TrilinosHandle const *)arg1)->solve(arg2,(double const *)arg3,arg4);
}


SWIGEXPORT void swigc_TrilinosHandle_finalize(void* farg1) {
  ForTrilinos::TrilinosHandle *arg1 = (ForTrilinos::TrilinosHandle *) 0 ;
  
  arg1 = (ForTrilinos::TrilinosHandle *)(farg1); 
  (arg1)->finalize();
}


SWIGEXPORT void swigc_delete_TrilinosHandle(void* farg1) {
  ForTrilinos::TrilinosHandle *arg1 = (ForTrilinos::TrilinosHandle *) 0 ;
  
  arg1 = (ForTrilinos::TrilinosHandle *)(farg1); 
  delete arg1;
}


#ifdef __cplusplus
}
#endif

