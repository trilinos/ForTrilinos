/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 4.0.2+fortran
 *
 * This file is not intended to be easily readable and contains a number of
 * coding conventions designed to improve portability and efficiency. Do not make
 * changes to this file unless you know what you are doing--modify the SWIG
 * interface file instead.
 * ----------------------------------------------------------------------------- */

/*
 * Copyright 2017-2020, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */



#ifndef SWIGFORTRAN
#define SWIGFORTRAN
#endif


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


#ifndef SWIGEXTERN
# ifdef __cplusplus
#   define SWIGEXTERN extern
# else
#   define SWIGEXTERN
# endif
#endif


#define SWIG_exception_impl(DECL, CODE, MSG, RETURNNULL) \
 { throw std::logic_error("In " DECL ": " MSG); }


#ifdef __cplusplus
extern "C" {
#endif
SWIGEXPORT void SWIG_check_unhandled_exception_impl(const char* decl);
SWIGEXPORT void SWIG_store_exception(const char* decl, int errcode, const char *msg);
#ifdef __cplusplus
}
#endif


#undef SWIG_exception_impl
#define SWIG_exception_impl(DECL, CODE, MSG, RETURNNULL) \
    SWIG_store_exception(DECL, CODE, MSG); RETURNNULL;

/*  Errors in SWIG */
#define  SWIG_UnknownError    	   -1
#define  SWIG_IOError        	   -2
#define  SWIG_RuntimeError   	   -3
#define  SWIG_IndexError     	   -4
#define  SWIG_TypeError      	   -5
#define  SWIG_DivisionByZero 	   -6
#define  SWIG_OverflowError  	   -7
#define  SWIG_SyntaxError    	   -8
#define  SWIG_ValueError     	   -9
#define  SWIG_SystemError    	   -10
#define  SWIG_AttributeError 	   -11
#define  SWIG_MemoryError    	   -12
#define  SWIG_NullReferenceError   -13




enum SwigMemFlags {
    SWIG_MEM_OWN = 0x01,
    SWIG_MEM_RVALUE = 0x02,
};


namespace swig {
enum AssignmentType {
  ASSIGNMENT_DEFAULT,
  ASSIGNMENT_NODESTRUCT,
  ASSIGNMENT_SMARTPTR
};
}

#define SWIGPOLICY_Belos_DefaultSolverParameters swig::ASSIGNMENT_DEFAULT

#include <stdexcept>


/* Support for the `contract` feature.
 *
 * Note that RETURNNULL is first because it's inserted via a 'Replaceall' in
 * the fortran.cxx file.
 */
#define SWIG_contract_assert(RETURNNULL, EXPR, MSG) \
 if (!(EXPR)) { SWIG_exception_impl("$decl", SWIG_ValueError, MSG, RETURNNULL); } 


#define SWIGVERSION 0x040002 
#define SWIG_VERSION SWIGVERSION


#define SWIG_as_voidptr(a) const_cast< void * >(static_cast< const void * >(a)) 
#define SWIG_as_voidptrptr(a) ((void)SWIG_as_voidptr(*a),reinterpret_cast< void** >(a)) 


#include <string>


#include <BelosTypes.hpp>


#include <stdlib.h>
#ifdef _MSC_VER
# ifndef strtoull
#  define strtoull _strtoui64
# endif
# ifndef strtoll
#  define strtoll _strtoi64
# endif
#endif


struct SwigArrayWrapper {
    void* data;
    size_t size;
};


SWIGINTERN SwigArrayWrapper SwigArrayWrapper_uninitialized() {
  SwigArrayWrapper result;
  result.data = NULL;
  result.size = 0;
  return result;
}


#include <string.h>


struct SwigClassWrapper {
    void* cptr;
    int cmemflags;
};


SWIGINTERN SwigClassWrapper SwigClassWrapper_uninitialized() {
    SwigClassWrapper result;
    result.cptr = NULL;
    result.cmemflags = 0;
    return result;
}


namespace swig {

template<class T, AssignmentType A>
struct DestructorPolicy {
  static SwigClassWrapper destroy(SwigClassWrapper self) {
    delete static_cast<T*>(self.cptr);
    return SwigClassWrapper_uninitialized();
  }
};
template<class T>
struct DestructorPolicy<T, ASSIGNMENT_NODESTRUCT> {
  static SwigClassWrapper destroy(SwigClassWrapper) {
    SWIG_exception_impl("assignment", SWIG_TypeError, "Invalid assignment: class type has private destructor", return SwigClassWrapper_uninitialized());
  }
};
}


namespace swig {

SWIGINTERN SwigClassWrapper capture(SwigClassWrapper other) {
  other.cmemflags &= ~SWIG_MEM_RVALUE;
  return other;
}

template<class T, AssignmentType A>
struct AssignmentPolicy {
  static SwigClassWrapper destroy(SwigClassWrapper self) {
    return DestructorPolicy<T, A>::destroy(self);
  }
  static SwigClassWrapper alias(SwigClassWrapper other) {
    SwigClassWrapper self = other;
    self.cmemflags &= ~SWIG_MEM_OWN;
    return self;
  }
  static SwigClassWrapper move_alias(SwigClassWrapper self, SwigClassWrapper other) {
    if (self.cmemflags & SWIG_MEM_OWN) {
      destroy(self);
    }
    return capture(other);
  }
  static SwigClassWrapper copy_alias(SwigClassWrapper self, SwigClassWrapper other) {
    if (self.cmemflags & SWIG_MEM_OWN) {
      destroy(self);
    }
    return capture(other);
  }
};

template<class T>
struct AssignmentPolicy<T, ASSIGNMENT_SMARTPTR> {
  static SwigClassWrapper destroy(SwigClassWrapper self) {
    return DestructorPolicy<T, ASSIGNMENT_SMARTPTR>::destroy(self);
  }
  static SwigClassWrapper alias(SwigClassWrapper other) {
    SwigClassWrapper self;
    self.cptr = new T(*static_cast<T*>(other.cptr));
    self.cmemflags = other.cmemflags | SWIG_MEM_OWN;
    return self;
  }
  static SwigClassWrapper move_alias(SwigClassWrapper self, SwigClassWrapper other) {
    self = copy_alias(self, other);
    self.cmemflags = other.cmemflags & ~SWIG_MEM_RVALUE;
    destroy(other);
    return self;
  }
  static SwigClassWrapper copy_alias(SwigClassWrapper self, SwigClassWrapper other) {
    // LHS and RHS should both 'own' their shared pointers
    T *pself = static_cast<T*>(self.cptr);
    T *pother = static_cast<T*>(other.cptr);
    *pself = *pother;
    return self;
  }
};

} // end namespace swig

template<class T, swig::AssignmentType A>
SWIGINTERN void SWIG_assign(SwigClassWrapper* self, SwigClassWrapper other) {
  typedef swig::AssignmentPolicy<T, A> Policy_t;

  if (self->cptr == NULL) {
    /* LHS is unassigned */
    if (other.cmemflags & SWIG_MEM_RVALUE) {
      /* Capture pointer from RHS, clear 'moving' flag */
      *self = swig::capture(other);
    } else {
      /* Aliasing another class; clear ownership or copy smart pointer */
      *self = Policy_t::alias(other);
    }
  } else if (other.cptr == NULL) {
    /* Replace LHS with a null pointer */
    *self = Policy_t::destroy(*self);
  } else if (self->cptr == other.cptr) {
    /* Self-assignment: ignore */
  } else if (other.cmemflags & SWIG_MEM_RVALUE) {
    /* Transferred ownership from a variable that's about to be lost.
     * Move-assign and delete the transient data */
    *self = Policy_t::move_alias(*self, other);
  } else {
    /* RHS shouldn't be deleted, alias to LHS */
    *self = Policy_t::copy_alias(*self, other);
  }
}

template<class T, swig::AssignmentType A>
SWIGINTERN void SWIG_free_rvalue(SwigClassWrapper other) {
  typedef swig::AssignmentPolicy<T, A> Policy_t;
  if (other.cmemflags & SWIG_MEM_RVALUE 
      && other.cmemflags & SWIG_MEM_OWN) {
    /* We own *and* are being passed an expiring value */
    Policy_t::destroy(other);
  }
}


extern "C" {
SWIGEXPORT SwigArrayWrapper _wrap_convertReturnTypeToString(int const *farg1) {
  SwigArrayWrapper fresult ;
  Belos::ReturnType arg1 ;
  std::string result;
  
  arg1 = (Belos::ReturnType)(*farg1);
  result = Belos::convertReturnTypeToString(arg1);
  fresult.size = (&result)->size();
  if (fresult.size > 0) {
    fresult.data = malloc(fresult.size);
    memcpy(fresult.data, (&result)->c_str(), fresult.size);
  } else {
    fresult.data = NULL;
  }
  return fresult;
}


SWIGEXPORT SWIGEXTERN const int _wrap_BelosPassed = (int)(Belos::Passed);

SWIGEXPORT SWIGEXTERN const int _wrap_BelosFailed = (int)(Belos::Failed);

SWIGEXPORT SWIGEXTERN const int _wrap_BelosUndefined = (int)(Belos::Undefined);

SWIGEXPORT SWIGEXTERN const int _wrap_BelosProblem = (int)(Belos::Problem);

SWIGEXPORT SWIGEXTERN const int _wrap_BelosRecycleSubspace = (int)(Belos::RecycleSubspace);

SWIGEXPORT SwigArrayWrapper _wrap_convertStatusTypeToString(int const *farg1) {
  SwigArrayWrapper fresult ;
  Belos::StatusType arg1 ;
  std::string result;
  
  arg1 = (Belos::StatusType)(*farg1);
  result = Belos::convertStatusTypeToString(arg1);
  fresult.size = (&result)->size();
  if (fresult.size > 0) {
    fresult.data = malloc(fresult.size);
    memcpy(fresult.data, (&result)->c_str(), fresult.size);
  } else {
    fresult.data = NULL;
  }
  return fresult;
}


SWIGEXPORT int _wrap_convertStringToStatusType(SwigArrayWrapper *farg1) {
  int fresult ;
  std::string *arg1 = 0 ;
  std::string tempstr1 ;
  Belos::StatusType result;
  
  tempstr1 = std::string(static_cast<char *>(farg1->data), farg1->size);
  arg1 = &tempstr1;
  result = (Belos::StatusType)Belos::convertStringToStatusType((std::string const &)*arg1);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_convertStringToNormType(SwigArrayWrapper *farg1) {
  int fresult ;
  std::string *arg1 = 0 ;
  std::string tempstr1 ;
  Belos::NormType result;
  
  tempstr1 = std::string(static_cast<char *>(farg1->data), farg1->size);
  arg1 = &tempstr1;
  result = (Belos::NormType)Belos::convertStringToNormType((std::string const &)*arg1);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT int _wrap_convertStringToScaleType(SwigArrayWrapper *farg1) {
  int fresult ;
  std::string *arg1 = 0 ;
  std::string tempstr1 ;
  Belos::ScaleType result;
  
  tempstr1 = std::string(static_cast<char *>(farg1->data), farg1->size);
  arg1 = &tempstr1;
  result = (Belos::ScaleType)Belos::convertStringToScaleType((std::string const &)*arg1);
  fresult = (int)(result);
  return fresult;
}


SWIGEXPORT SwigArrayWrapper _wrap_convertScaleTypeToString(int const *farg1) {
  SwigArrayWrapper fresult ;
  Belos::ScaleType arg1 ;
  std::string result;
  
  arg1 = (Belos::ScaleType)(*farg1);
  result = Belos::convertScaleTypeToString(arg1);
  fresult.size = (&result)->size();
  if (fresult.size > 0) {
    fresult.data = malloc(fresult.size);
    memcpy(fresult.data, (&result)->c_str(), fresult.size);
  } else {
    fresult.data = NULL;
  }
  return fresult;
}


SWIGEXPORT SWIGEXTERN const int _wrap_BelosErrors = (int)(Belos::Errors);

SWIGEXPORT SWIGEXTERN const int _wrap_BelosWarnings = (int)(Belos::Warnings);

SWIGEXPORT SWIGEXTERN const int _wrap_BelosIterationDetails = (int)(Belos::IterationDetails);

SWIGEXPORT SWIGEXTERN const int _wrap_BelosOrthoDetails = (int)(Belos::OrthoDetails);

SWIGEXPORT SWIGEXTERN const int _wrap_BelosFinalSummary = (int)(Belos::FinalSummary);

SWIGEXPORT SWIGEXTERN const int _wrap_BelosTimingDetails = (int)(Belos::TimingDetails);

SWIGEXPORT SWIGEXTERN const int _wrap_BelosStatusTestDetails = (int)(Belos::StatusTestDetails);

SWIGEXPORT SWIGEXTERN const int _wrap_BelosDebug = (int)(Belos::Debug);

SWIGEXPORT SwigArrayWrapper _wrap_convertMsgTypeToString(int const *farg1) {
  SwigArrayWrapper fresult ;
  Belos::MsgType arg1 ;
  std::string result;
  
  arg1 = (Belos::MsgType)(*farg1);
  result = Belos::convertMsgTypeToString(arg1);
  fresult.size = (&result)->size();
  if (fresult.size > 0) {
    fresult.data = malloc(fresult.size);
    memcpy(fresult.data, (&result)->c_str(), fresult.size);
  } else {
    fresult.data = NULL;
  }
  return fresult;
}


SWIGEXPORT double _wrap_DefaultSolverParameters_convTol_get() {
  double fresult ;
  double result;
  
  result = (double)(double)Belos::DefaultSolverParameters::convTol;
  fresult = (double)(result);
  return fresult;
}


SWIGEXPORT double _wrap_DefaultSolverParameters_polyTol_get() {
  double fresult ;
  double result;
  
  result = (double)(double)Belos::DefaultSolverParameters::polyTol;
  fresult = (double)(result);
  return fresult;
}


SWIGEXPORT double _wrap_DefaultSolverParameters_orthoKappa_get() {
  double fresult ;
  double result;
  
  result = (double)(double)Belos::DefaultSolverParameters::orthoKappa;
  fresult = (double)(result);
  return fresult;
}


SWIGEXPORT double _wrap_DefaultSolverParameters_resScaleFactor_get() {
  double fresult ;
  double result;
  
  result = (double)(double)Belos::DefaultSolverParameters::resScaleFactor;
  fresult = (double)(result);
  return fresult;
}


SWIGEXPORT double _wrap_DefaultSolverParameters_impTolScale_get() {
  double fresult ;
  double result;
  
  result = (double)(double)Belos::DefaultSolverParameters::impTolScale;
  fresult = (double)(result);
  return fresult;
}


SWIGEXPORT SwigClassWrapper _wrap_new_DefaultSolverParameters() {
  SwigClassWrapper fresult ;
  Belos::DefaultSolverParameters *result = 0 ;
  
  result = (Belos::DefaultSolverParameters *)new Belos::DefaultSolverParameters();
  fresult.cptr = (void*)result;
  fresult.cmemflags = SWIG_MEM_RVALUE | (1 ? SWIG_MEM_OWN : 0);
  return fresult;
}


SWIGEXPORT void _wrap_delete_DefaultSolverParameters(SwigClassWrapper *farg1) {
  Belos::DefaultSolverParameters *arg1 = (Belos::DefaultSolverParameters *) 0 ;
  
  arg1 = (Belos::DefaultSolverParameters *)farg1->cptr;
  delete arg1;
}


SWIGEXPORT void _wrap_DefaultSolverParameters_op_assign__(SwigClassWrapper *farg1, SwigClassWrapper *farg2) {
  Belos::DefaultSolverParameters *arg1 = (Belos::DefaultSolverParameters *) 0 ;
  Belos::DefaultSolverParameters *arg2 = 0 ;
  
  (void)sizeof(arg1);
  (void)sizeof(arg2);
  SWIG_assign<Belos::DefaultSolverParameters, SWIGPOLICY_Belos_DefaultSolverParameters>(farg1, *farg2);
  
}


} // extern

