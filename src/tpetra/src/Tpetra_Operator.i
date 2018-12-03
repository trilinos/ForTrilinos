/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */

// Dependencies
%import <Teuchos_Comm.i>
%import <Tpetra_MultiVector.i>

%{
#include "Tpetra_Operator.hpp"
%}

// =======================================================================
// Instantiate the Operator
// =======================================================================

//%include "Tpetra_Operator.hpp"
namespace Tpetra {
  template <class Scalar = ::Tpetra::Details::DefaultTypes::scalar_type,
            class LocalOrdinal = ::Tpetra::Details::DefaultTypes::local_ordinal_type,
            class GlobalOrdinal = ::Tpetra::Details::DefaultTypes::global_ordinal_type,
            class Node = ::Tpetra::Details::DefaultTypes::node_type>
class Operator {
    Operator(); // Disable constructor
public:
};
}

%teuchos_rcp(Tpetra::Operator<SC,LO,GO,NO>)
%template(TpetraOperator) Tpetra::Operator<SC,LO,GO,NO>;

// =======================================================================
// Add the subclass
// =======================================================================

%teuchos_rcp(ForTpetraOperator)

%insert("header") %{
extern "C" {
/* Fortran BIND(C) function */
void swigd_ForTpetraOperator_apply(
        SwigClassWrapper const *fself,
        SwigClassWrapper const *farg1,
        SwigClassWrapper const *farg2,
        int *farg3,
        double const *farg4,
        double const *farg5
        );
SwigClassWrapper swigd_ForTpetraOperator_getDomainMap(
        SwigClassWrapper const *fself
        );
SwigClassWrapper swigd_ForTpetraOperator_getRangeMap(
        SwigClassWrapper const *fself
        );
}
%}

%fortranprepend ForTpetraOperator::~ForTpetraOperator() %{
  type(C_PTR) :: fself_ptr
  type(ForTpetraOperatorHandle), pointer :: handle
  fself_ptr = swigc_ForTpetraOperator_fhandle(self%swigdata)
  call c_f_pointer(cptr=fself_ptr, fptr=handle)
%}

%fortranappend ForTpetraOperator::~ForTpetraOperator() %{
  ! Release the allocated handle
  deallocate(handle)
%}

%inline %{
  class ForTpetraOperator : public Tpetra::Operator<SC,LO,GO,NO> {
    // Pointer to polymorphic fortran pointer
    void* fhandle_;
   public:
    /* DIRECTOR FUNCTIONS */
    const void* fhandle() const { assert(fhandle_); return this->fhandle_; }
    void init(void* fh) { fhandle_ = fh; }

    /* TPETRA */
    typedef Tpetra::MultiVector<SC,LO,GO,NO> vector_type;
    typedef Tpetra::Map<LO,GO,NO> map_type;
    ForTpetraOperator() : fhandle_(NULL) { /* * */ }


    virtual Teuchos::RCP<const map_type> getDomainMap() const {
      /* construct "this" pointer */
      Teuchos::RCP<ForTpetraOperator> tempthis(
             const_cast<ForTpetraOperator*>(this) SWIG_NO_NULL_DELETER_0);
      SwigClassWrapper self;
      self.ptr = &tempthis;
      self.mem = SWIG_CREF; // since this function is const

      SwigClassWrapper fresult = swigd_ForTpetraOperator_getDomainMap(&self);

      Teuchos::RCP<const map_type>* smartresult = static_cast< Teuchos::RCP<const map_type>* >(fresult.ptr);
      return *smartresult;
    }

    virtual Teuchos::RCP<const map_type> getRangeMap() const {
      /* construct "this" pointer */
      Teuchos::RCP<ForTpetraOperator> tempthis(
             const_cast<ForTpetraOperator*>(this) SWIG_NO_NULL_DELETER_0);
      SwigClassWrapper self;
      self.ptr = &tempthis;
      self.mem = SWIG_CREF; // since this function is const

      SwigClassWrapper fresult = swigd_ForTpetraOperator_getRangeMap(&self);

      Teuchos::RCP<const map_type>* smartresult = static_cast< Teuchos::RCP<const map_type>* >(fresult.ptr);
      return *smartresult;
    }

    virtual void apply(const vector_type &X, vector_type &Y,
                       Teuchos::ETransp mode, SC alpha, SC beta) const
    {
      /* construct "this" pointer */
      Teuchos::RCP<ForTpetraOperator> tempthis(
             const_cast<ForTpetraOperator*>(this) SWIG_NO_NULL_DELETER_0);
      SwigClassWrapper self;
      self.ptr = &tempthis;
      self.mem = SWIG_CREF; // since this function is const

      /* convert X -> class wrapper */
      Teuchos::RCP<const Tpetra::MultiVector<SC,LO,GO,NO> > temprcp1(&X SWIG_NO_NULL_DELETER_0);

      SwigClassWrapper farg1;
      farg1.ptr = &temprcp1;
      farg1.mem = SWIG_CREF; // X is const

      Teuchos::RCP< Tpetra::MultiVector<SC,LO,GO,NO> > temprcp2(&Y SWIG_NO_NULL_DELETER_0);

      SwigClassWrapper farg2;
      farg2.ptr = &temprcp2;
      farg2.mem = SWIG_REF; // Y is mutable

      /* convert scalars to wrappers */
      int farg3 = mode;
      double farg4 = alpha;
      double farg5 = beta;

      swigd_ForTpetraOperator_apply(&self, &farg1, &farg2, &farg3, &farg4, &farg5);
    }
  };
%}

%insert("ftypes") %{
  type :: ForTpetraOperatorHandle
    class(ForTpetraOperator), pointer :: data
  end type
%}

%insert("fpublic") %{
public :: init_ForTpetraOperator
%}

%insert("fwrapper") %{
! Convert a ISO-C class pointer struct into a user Fortran native pointer
subroutine c_f_pointer_ForTpetraOperator(clswrap, fptr)
  type(SwigClassWrapper), intent(in) :: clswrap
  class(ForTpetraOperator), pointer, intent(out) :: fptr
  type(ForTpetraOperatorHandle), pointer :: handle
  type(C_PTR) :: fself_ptr
  ! Convert C handle to fortran pointer
  fself_ptr = swigc_ForTpetraOperator_fhandle(clswrap)
  ! *** NOTE *** : gfortran 5 through 7 falsely claim the next line is not standards compliant. Since 'handle' is a scalar and
  ! not an array it should be OK, but TS29113 explicitly removes the interoperability requirement for fptr.
  ! Error: TS 29113/TS 18508: Noninteroperable array FPTR at (1) to C_F_POINTER: Expression is a noninteroperable derived type
  ! see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=84924
  call c_f_pointer(cptr=fself_ptr, fptr=handle)
  if (.not. associated(handle)) stop 1
  ! Access the pointer inside that
  fptr => handle%data
  if (.not. associated(fptr)) stop 2
end subroutine

! Convert SWIG array wrapper to temporary Fortran string, pass to the fortran
! cb function, convert back to SWIG array wrapper.
! This function must have input/output arguments compatible with ISO C, and it must be marked with "bind(C)"
subroutine swigd_ForTpetraOperator_apply(fself, farg1, farg2, farg3, farg4, farg5) &
    bind(C, name="swigd_ForTpetraOperator_apply")
  use, intrinsic :: ISO_C_BINDING
  implicit none
  type(SwigClassWrapper), intent(in) :: fself
  type(SwigClassWrapper), intent(in) :: farg1
  type(SwigClassWrapper), intent(in) :: farg2
  integer(C_INT), intent(in) :: farg3
  real(C_DOUBLE), intent(in) :: farg4
  real(C_DOUBLE), intent(in) :: farg5

  class(ForTpetraOperator), pointer :: self
  type(TpetraMultiVector) :: x
  type(TpetraMultiVector) :: y
  integer(kind(TeuchosETransp)) :: mode
  real(C_DOUBLE) :: alpha
  real(C_DOUBLE) :: beta

  ! Get pointer to Fortran object from class wrapper
  call c_f_pointer_ForTpetraOperator(fself, self)
  if (.not. associated(self)) stop 3

  ! Convert class references to fortran proxy references
  x%swigdata = farg1
  y%swigdata = farg2

  ! Copy scalars
  mode  = int(farg3, kind(TeuchosETransp))
  alpha = farg4
  beta  = farg5

  ! Call fortran function pointer with native fortran input/output
  call self%apply(x, y, mode, alpha, beta)
end subroutine

function swigd_ForTpetraOperator_getDomainMap(fself) &
    bind(C, name="swigd_ForTpetraOperator_getDomainMap") &
    result(fresult)
  use, intrinsic :: ISO_C_BINDING
  implicit none
  type(SwigClassWrapper), intent(in) :: fself
  type(SwigClassWrapper) :: fresult

  class(ForTpetraOperator), pointer :: self
  type(TpetraMap) :: result

  ! Get pointer to Fortran object from class Handle
  call c_f_pointer_ForTpetraOperator(fself, self)
  if (.not. associated(self)) stop 3

  result = self%getDomainMap()

  fresult = result%swigdata
end function

function swigd_ForTpetraOperator_getRangeMap(fself) &
    bind(C, name="swigd_ForTpetraOperator_getRangeMap") &
    result(fresult)
  use, intrinsic :: ISO_C_BINDING
  implicit none
  type(SwigClassWrapper), intent(in) :: fself
  type(SwigClassWrapper) :: fresult

  class(ForTpetraOperator), pointer :: self
  type(TpetraMap) :: result

  ! Get pointer to Fortran object from class Handle
  call c_f_pointer_ForTpetraOperator(fself, self)
  if (.not. associated(self)) stop 3

  result = self%getRangeMap()

  fresult = result%swigdata
end function

subroutine init_ForTpetraOperator(self)
  ! Note: subclass should call `self = ForTpetraOperator()` in its
  ! initialization code *before* doing this
  class(ForTpetraOperator), target :: self
  type(ForTpetraOperatorHandle), pointer :: handle
  allocate(handle)
  handle%data => self
  call swigc_ForTpetraOperator_init(self%swigdata, c_loc(handle))
end subroutine
%}

