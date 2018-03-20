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

%include "Tpetra_Operator.hpp"

%teuchos_rcp(Tpetra::Operator<SC,LO,GO,NO>)
%template() Tpetra::Operator<SC,LO,GO,NO>;

// =======================================================================
// Add the subclass
// =======================================================================

%teuchos_rcp(TpetraOperator)

%insert("header") %{
extern "C" {
/* Fortran BIND(C) function */
SwigArrayWrapper swigd_TpetraOperator_apply(
        SwigClassWrapper const *fself,
        SwigClassWrapper const *farg1,
        SwigClassWrapper const *farg2,
        int *farg3,
        double const *farg4,
        double const *farg5
        );
}
%}

%fortranprepend TpetraOperator::~TpetraOperator() %{
  type(C_PTR) :: fself_ptr
  type(TpetraOperatorWrapper), pointer :: handle
  fself_ptr = swigc_TpetraOperator_fhandle(self%swigdata)
  call c_f_pointer(cptr=fself_ptr, fptr=handle)
%}

%fortranappend TpetraOperator::~TpetraOperator() %{
  ! Release the allocated handle
  deallocate(handle)
%}

%inline %{
  class TpetraOperator : public Tpetra::Operator<SC,LO,GO,NO> {
    // Pointer to polymorphic fortran pointer
    void* fhandle_;
   public:
    /* DIRECTOR FUNCTIONS */
    const void* fhandle() const { return this->fhandle_; }
    void init(void* fh) { fhandle_ = fh; }

    /* TPETRA */
    typedef Tpetra::MultiVector<SC,LO,GO,NO> MV_type;
    typedef Tpetra::Map<LO,GO,NO> Map_type;
    TpetraOperator() : fhandle_(NULL) { /* * */ }


    virtual Teuchos::RCP<const Map_type> getDomainMap() const { return Teuchos::null; }
    virtual Teuchos::RCP<const Map_type> getRangeMap() const { return Teuchos::null; }

    virtual void apply(const MV_type &X, MV_type &Y,
                       Teuchos::ETransp mode, SC alpha, SC beta) const
    {
      /* construct "this" pointer */
      SwigClassWrapper self;
      self.ptr = const_cast<TpetraOperator*>(this);
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

      swigd_TpetraOperator_apply(&self, &farg1, &farg2, &farg3, &farg4, &farg5);
    }
  };
%}

%insert("ftypes") %{
  type :: TpetraOperatorWrapper
    class(TpetraOperator), pointer :: data
  end type
%}

%insert("fpublic") %{
public :: init_TpetraOperator
%}

%insert("fwrapper") %{
! Convert a ISO-C class pointer struct into a user Fortran native pointer
subroutine c_f_pointer_TpetraOperator(clswrap, fptr)
  type(SwigClassWrapper), intent(in) :: clswrap
  class(TpetraOperator), pointer, intent(out) :: fptr
  type(TpetraOperatorWrapper), pointer :: handle
  type(C_PTR) :: fself_ptr
  ! Convert C handle to fortran pointer
  fself_ptr = swigc_TpetraOperator_fhandle(clswrap)
  ! *** NOTE *** : gfortran 5 through 7 falsely claim the next line is not standards compliant. Since 'handle' is a scalar and
  ! not an array it should be OK, but TS29113 explicitly removes the interoperability requirement for fptr.
  ! Error: TS 29113/TS 18508: Noninteroperable array FPTR at (1) to C_F_POINTER: Expression is a noninteroperable derived type
  ! see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=84924
  call c_f_pointer(cptr=fself_ptr, fptr=handle)
  ! Access the pointer inside that
  fptr => handle%data
end subroutine

! Convert SWIG array wrapper to temporary Fortran string, pass to the fortran
! cb function, convert back to SWIG array wrapper.
! This function must have input/output arguments compatible with ISO C, and it must be marked with "bind(C)"
subroutine swigd_TpetraOperator_apply(fself, farg1, farg2, farg3, farg4, farg5) &
    bind(C, name="swigd_TpetraOperator_apply")
  use, intrinsic :: ISO_C_BINDING
  implicit none
  type(SwigClassWrapper), intent(in) :: fself
  type(SwigClassWrapper), intent(in) :: farg1
  type(SwigClassWrapper), intent(in) :: farg2
  integer(C_INT), intent(in) :: farg3
  real(C_DOUBLE), intent(in) :: farg4
  real(C_DOUBLE), intent(in) :: farg5

  class(TpetraOperator), pointer :: self
  type(TpetraMultiVector) :: x
  type(TpetraMultiVector) :: y
  integer(kind(TeuchosETransp)) :: mode
  real(C_DOUBLE) :: alpha
  real(C_DOUBLE) :: beta

  ! Get pointer to Fortran object from class wrapper
  call c_f_pointer_TpetraOperator(fself, self)
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

subroutine init_TpetraOperator(self)
  class(TpetraOperator), target :: self
  type(TpetraOperatorWrapper), pointer :: handle
  allocate(handle)
  handle%data => self
  self%swigdata = swigc_new_TpetraOperator()
  call swigc_TpetraOperator_init(self%swigdata, c_loc(handle))
end subroutine
%}

