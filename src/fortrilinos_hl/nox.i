
/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
%include <boost_shared_ptr.i>

%{
#include "fortrilinos_hl/model_evaluator.hpp"
%}

%include "model_evaluator.hpp"

%teuchos_rcp(ForTrilinos::ModelEvaluator<SC,LO,GO,NO>)
%template(ForTrilinosModelEvaluator) ForTrilinos::ModelEvaluator<SC,LO,GO,NO>;

// =======================================================================
// Add the subclass
// =======================================================================
%teuchos_rcp(ForModelEvaluator)

%insert("header") %{
extern "C" {
/* Fortran BIND(C) function */
void swigd_ForModelEvaluator_evaluate_residual(
        SwigClassWrapper const *fself,
        SwigClassWrapper const *farg1,
        SwigClassWrapper const *farg2
        );
void swigd_ForModelEvaluator_evaluate_jacobian(
        SwigClassWrapper const *fself,
        SwigClassWrapper const *farg1,
        SwigClassWrapper const *farg2
        );
void swigd_ForModelEvaluator_evaluate_preconditioner(
        SwigClassWrapper const *fself,
        SwigClassWrapper const *farg1,
        SwigClassWrapper const *farg2
        );
SwigClassWrapper swigd_ForModelEvaluator_get_x_map(
        SwigClassWrapper const *fself
        );
SwigClassWrapper swigd_ForModelEvaluator_get_f_map(
        SwigClassWrapper const *fself
        );
SwigClassWrapper swigd_ForModelEvaluator_create_operator(
        SwigClassWrapper const *fself
        );
}
%}

%fortranprepend ForModelEvaluator::~ForModelEvaluator() %{
  type(C_PTR) :: fself_ptr
  type(ForModelEvaluatorHandle), pointer :: handle
  fself_ptr = swigc_ForModelEvaluator_fhandle(self%swigdata)
  call c_f_pointer(cptr=fself_ptr, fptr=handle)
%}

%fortranappend ForModelEvaluator::~ForModelEvaluator() %{
  ! Release the allocated handle
  deallocate(handle)
%}

%inline %{
  // FIXME for some reason SWIG_NO_NULL_DELETER is included *after* this class definition
#define SWIG_NO_NULL_DELETER_0 ,Teuchos::RCP_WEAK_NO_DEALLOC
  class ForModelEvaluator : public ForTrilinos::ModelEvaluator<SC,LO,GO,NO> {
    // Pointer to polymorphic fortran pointer
    void* fhandle_;
   public:
    /* DIRECTOR FUNCTIONS */
    const void* fhandle() const { assert(fhandle_); return this->fhandle_; }
    void init(void* fh) { fhandle_ = fh; }

    /* ModelEvaluator */
    typedef Tpetra::Map<LO,GO,NO> map_type;
    typedef Tpetra::MultiVector<SC,LO,GO,NO> multivector_type;
    typedef Tpetra::Operator<SC,LO,GO,NO> operator_type;

    void setup(Teuchos::RCP<Teuchos::ParameterList>& plist) {
      ForTrilinos::ModelEvaluator<SC,LO,GO,NO>::setup(plist);
    }

    virtual void evaluate_residual(const Teuchos::RCP<const multivector_type>& x,
                                   Teuchos::RCP<multivector_type>& f) const override {
      /* construct "this" pointer */
      Teuchos::RCP<ForModelEvaluator> tempthis(
             const_cast<ForModelEvaluator*>(this) SWIG_NO_NULL_DELETER_0);
      SwigClassWrapper self;
      self.cptr = &tempthis;
      self.cmemflags = SWIG_MEM_CONST; // since this function is const

      /* convert x -> class wrapper */
      SwigClassWrapper farg1;
      farg1.cptr = const_cast<Teuchos::RCP<const multivector_type>*>(&x);
      farg1.cmemflags = SWIG_MEM_CONST; // x is const

      /* convert f -> class wrapper */
      SwigClassWrapper farg2;
      farg2.cptr = &f;
      farg2.cmemflags = 0; // f is mutable

      swigd_ForModelEvaluator_evaluate_residual(&self, &farg1, &farg2);
    }

    virtual void evaluate_jacobian(const Teuchos::RCP<const multivector_type>& x,
                                   Teuchos::RCP<operator_type>& J) const override {
      /* construct "this" pointer */
      Teuchos::RCP<ForModelEvaluator> tempthis(
             const_cast<ForModelEvaluator*>(this) SWIG_NO_NULL_DELETER_0);
      SwigClassWrapper self;
      self.cptr = &tempthis;
      self.cmemflags = SWIG_MEM_CONST; // since this function is const

      /* convert x -> class wrapper */
      SwigClassWrapper farg1;
      farg1.cptr = const_cast<Teuchos::RCP<const multivector_type>*>(&x);
      farg1.cmemflags = SWIG_MEM_CONST; // x is const

      /* convert J -> class wrapper */
      SwigClassWrapper farg2;
      farg2.cptr = &J;
      farg2.cmemflags = 0; // f is mutable

      swigd_ForModelEvaluator_evaluate_jacobian(&self, &farg1, &farg2);
    }

    virtual void evaluate_preconditioner(const Teuchos::RCP<const multivector_type>& x,
                                         Teuchos::RCP<operator_type>& M) const override {
      /* construct "this" pointer */
      Teuchos::RCP<ForModelEvaluator> tempthis(
             const_cast<ForModelEvaluator*>(this) SWIG_NO_NULL_DELETER_0);
      SwigClassWrapper self;
      self.cptr = &tempthis;
      self.cmemflags = SWIG_MEM_CONST; // since this function is const

      /* convert x -> class wrapper */
      SwigClassWrapper farg1;
      farg1.cptr = const_cast<Teuchos::RCP<const multivector_type>*>(&x);
      farg1.cmemflags = SWIG_MEM_CONST; // x is const

      /* convert M -> class wrapper */
      SwigClassWrapper farg2;
      farg2.cptr = &M;
      farg2.cmemflags = 0; // f is mutable

      swigd_ForModelEvaluator_evaluate_preconditioner(&self, &farg1, &farg2);
    }

    virtual Teuchos::RCP<const map_type> get_x_map() const override {
      /* construct "this" pointer */
      Teuchos::RCP<ForModelEvaluator> tempthis(
             const_cast<ForModelEvaluator*>(this) SWIG_NO_NULL_DELETER_0);
      SwigClassWrapper self;
      self.cptr = &tempthis;
      self.cmemflags = SWIG_MEM_CONST; // since this function is const

      SwigClassWrapper fresult = swigd_ForModelEvaluator_get_x_map(&self);

      Teuchos::RCP<const map_type>* smartresult = static_cast< Teuchos::RCP<const map_type>* >(fresult.cptr);
      return *smartresult;
    }

    virtual Teuchos::RCP<const map_type> get_f_map() const override {
      /* construct "this" pointer */
      Teuchos::RCP<ForModelEvaluator> tempthis(
             const_cast<ForModelEvaluator*>(this) SWIG_NO_NULL_DELETER_0);
      SwigClassWrapper self;
      self.cptr = &tempthis;
      self.cmemflags = SWIG_MEM_CONST; // since this function is const

      SwigClassWrapper fresult = swigd_ForModelEvaluator_get_f_map(&self);

      Teuchos::RCP<const map_type>* smartresult = static_cast< Teuchos::RCP<const map_type>* >(fresult.cptr);
      return *smartresult;
    }

    virtual Teuchos::RCP<operator_type> create_operator() const override {
      /* construct "this" pointer */
      Teuchos::RCP<ForModelEvaluator> tempthis(
             const_cast<ForModelEvaluator*>(this) SWIG_NO_NULL_DELETER_0);
      SwigClassWrapper self;
      self.cptr = &tempthis;
      self.cmemflags = SWIG_MEM_CONST; // since this function is const

      SwigClassWrapper fresult = swigd_ForModelEvaluator_create_operator(&self);

      Teuchos::RCP<operator_type>* smartresult = static_cast< Teuchos::RCP<operator_type>* >(fresult.cptr);
      return *smartresult;
    }
  };
#undef SWIG_NO_NULL_DELETER_0
%}

%insert("fdecl") %{
  type :: ForModelEvaluatorHandle
    class(ForModelEvaluator), pointer :: data
  end type
%}

%insert("fdecl") %{
public :: init_ForModelEvaluator
%}

%insert("fsubprograms") %{
! Convert a ISO-C class pointer struct into a user Fortran native pointer
subroutine c_f_pointer_ForModelEvaluator(clswrap, fptr)
  type(SwigClassWrapper), intent(in) :: clswrap
  class(ForModelEvaluator), pointer, intent(out) :: fptr
  type(ForModelEvaluatorHandle), pointer :: handle
  type(C_PTR) :: fself_ptr
  ! Convert C handle to fortran pointer
  fself_ptr = swigc_ForModelEvaluator_fhandle(clswrap)
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

! This function must have input/output arguments compatible with ISO C, and it must be marked with "bind(C)"
subroutine swigd_ForModelEvaluator_evaluate_residual(fself, farg1, farg2) &
    bind(C, name="swigd_ForModelEvaluator_evaluate_residual")
  use, intrinsic :: ISO_C_BINDING
  implicit none
  type(SwigClassWrapper), intent(in) :: fself
  type(SwigClassWrapper), intent(inout) :: farg1
  type(SwigClassWrapper), intent(inout) :: farg2

  class(ForModelEvaluator), pointer :: self
  type(TpetraMultiVector) :: x
  type(TpetraMultiVector) :: f

  ! Get pointer to Fortran object from class wrapper
  call c_f_pointer_ForModelEvaluator(fself, self)
  if (.not. associated(self)) stop 3

  ! Convert class references to fortran proxy references
  x%swigdata = farg1
  f%swigdata = farg2

  ! Call fortran function pointer with native fortran input/output
  call self%evaluate_residual(x, f)
end subroutine

! This function must have input/output arguments compatible with ISO C, and it must be marked with "bind(C)"
subroutine swigd_ForModelEvaluator_evaluate_jacobian(fself, farg1, farg2) &
    bind(C, name="swigd_ForModelEvaluator_evaluate_jacobian")
  use, intrinsic :: ISO_C_BINDING
  implicit none
  type(SwigClassWrapper), intent(in) :: fself
  type(SwigClassWrapper), intent(in) :: farg1
  type(SwigClassWrapper), intent(in) :: farg2

  class(ForModelEvaluator), pointer :: self
  type(TpetraMultiVector) :: x
  type(ForTpetraOperator) :: J

  ! Get pointer to Fortran object from class wrapper
  call c_f_pointer_ForModelEvaluator(fself, self)
  if (.not. associated(self)) stop 3

  ! Convert class references to fortran proxy references
  x%swigdata = farg1
  J%swigdata = farg2

  ! Call fortran function pointer with native fortran input/output
  call self%evaluate_jacobian(x, J)
end subroutine

! This function must have input/output arguments compatible with ISO C, and it must be marked with "bind(C)"
subroutine swigd_ForModelEvaluator_evaluate_preconditioner(fself, farg1, farg2) &
    bind(C, name="swigd_ForModelEvaluator_evaluate_preconditioner")
  use, intrinsic :: ISO_C_BINDING
  implicit none
  type(SwigClassWrapper), intent(in) :: fself
  type(SwigClassWrapper), intent(in) :: farg1
  type(SwigClassWrapper), intent(in) :: farg2

  class(ForModelEvaluator), pointer :: self
  type(TpetraMultiVector) :: x
  type(ForTpetraOperator) :: M

  ! Get pointer to Fortran object from class wrapper
  call c_f_pointer_ForModelEvaluator(fself, self)
  if (.not. associated(self)) stop 3

  ! Convert class references to fortran proxy references
  x%swigdata = farg1
  M%swigdata = farg2

  ! Call fortran function pointer with native fortran input/output
  call self%evaluate_preconditioner(x, M)
end subroutine

function swigd_ForModelEvaluator_get_x_map(fself) &
    bind(C, name="swigd_ForModelEvaluator_get_x_map") &
    result(fresult)
  use, intrinsic :: ISO_C_BINDING
  implicit none
  type(SwigClassWrapper), intent(in) :: fself
  type(SwigClassWrapper) :: fresult

  class(ForModelEvaluator), pointer :: self
  type(TpetraMap) :: result

  ! Get pointer to Fortran object from class Handle
  call c_f_pointer_ForModelEvaluator(fself, self)
  if (.not. associated(self)) stop 3

  result = self%get_x_map()

  fresult = result%swigdata
end function

function swigd_ForModelEvaluator_get_f_map(fself) &
    bind(C, name="swigd_ForModelEvaluator_get_f_map") &
    result(fresult)
  use, intrinsic :: ISO_C_BINDING
  implicit none
  type(SwigClassWrapper), intent(in) :: fself
  type(SwigClassWrapper) :: fresult

  class(ForModelEvaluator), pointer :: self
  type(TpetraMap) :: result

  ! Get pointer to Fortran object from class Handle
  call c_f_pointer_ForModelEvaluator(fself, self)
  if (.not. associated(self)) stop 3

  result = self%get_f_map()

  fresult = result%swigdata
end function

function swigd_ForModelEvaluator_create_operator(fself) &
    bind(C, name="swigd_ForModelEvaluator_create_operator") &
    result(fresult)
  use, intrinsic :: ISO_C_BINDING
  implicit none
  type(SwigClassWrapper), intent(in) :: fself
  type(SwigClassWrapper) :: fresult

  class(ForModelEvaluator), pointer :: self
  type(TpetraOperator) :: result

  ! Get pointer to Fortran object from class Handle
  call c_f_pointer_ForModelEvaluator(fself, self)
  if (.not. associated(self)) stop 3

  result = self%create_operator()

  fresult = result%swigdata
end function

subroutine init_ForModelEvaluator(self)
  ! Note: subclass should call `self = ForModelEvaluator()` in its
  ! initialization code *before* doing this
  class(ForModelEvaluator), target :: self
  type(ForModelEvaluatorHandle), pointer :: handle
  allocate(handle)
  handle%data => self
  call swigc_ForModelEvaluator_init(self%swigdata, c_loc(handle))
end subroutine
%}

// All enums should be prefaced with NOX
%rename("NOX%s", %$isenumitem) "";
%rename("NOX%s", %$isenum)     "";

%{
#include <NOX_StatusTest_Generic.H>
%}

%ignore NOX::Solver::Generic;
%ignore NOX::StatusTest::Generic;
%ignore NOX::StatusTest::operator<<;

%include "NOX_StatusTest_Generic.H"

%{
#include "fortrilinos_hl/nox_solver.hpp"
%}

%include "nox_solver.hpp"

%teuchos_rcp(ForTrilinos::NOXSolver<SC,LO,GO,NO>)
%template(NOXSolver) ForTrilinos::NOXSolver<SC,LO,GO,NO>;
