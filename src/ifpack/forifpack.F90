!*********************************************************************
! ForTrilinos: Object-Oriented Fortran 2003 interface to Trilinos
!                Copyright 2010 Sandia Corporation
!
! Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
! the U.S. Government retains certain rights in this software.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright
!    notice, this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright
!    notice, this list of conditions and the following disclaimer in the
!    documentation and/or other materials provided with the distribution.
!
! 3. Neither the name of the Corporation nor the names of the
!    contributors may be used to endorse or promote products derived from
!    this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
! EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
! PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! Questions? Contact M. Nicole Lemaster (mnlemas@sandia.gov) 
!*********************************************************************

#include "ForTrilinos_config.h"
#ifdef HAVE_FORTRILINOS_IFPACK

module forifpack
  use iso_c_binding ,only : c_int,c_double,c_char,c_bool,c_ptr,c_long,c_float
  use ForTrilinos_enums
  use ForTrilinos_enum_wrappers
  implicit none   ! Prevent implicit typing

  ! This file provides Fortran interface blocks that bind the argument types,
  ! return value types, and procedure names to those in the C function prototypes
  ! in each of the CTrilinos/src/ifpack/CIfpack*.h header files.  The Fortran
  ! 2003 standard guarantees that the types and names used in these bindings
  ! interoperate with a standard-conforming, companion C compiler.

  ! Since this file contains only interface bodies, this interface block ends at
  ! the bottom of the file.

  interface

!> @name Ifpack interface
!! @{

  ! _________________ Ifpack interface bodies _________________


  !> <BR> Original C++ prototype:
  !! Ifpack();
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Ifpack_ID_t Ifpack_Create (  );

  type(FT_Ifpack_ID_t) function Ifpack_Create (  ) bind(C,name='Ifpack_Create')
    import :: FT_Ifpack_ID_t
    
  end function


  !> <BR> Original C++ prototype:
  !! ~Ifpack();
  !> <BR> <BR> CTrilinos prototype:
  !! void Ifpack_Destroy ( CT_Ifpack_ID_t * selfID );

  subroutine Ifpack_Destroy ( selfID ) bind(C,name='Ifpack_Destroy')
    import :: FT_Ifpack_ID_t
    
    type(FT_Ifpack_ID_t)                                          :: selfID
  end subroutine


  !> <BR> Original C++ prototype:
  !! static const char* toString(const EPrecType precType);
  !> <BR> <BR> CTrilinos prototype:
  !! const char * Ifpack_toString ( const CT_EPrecType_E_t precType );

  type(c_ptr) function Ifpack_toString ( precType ) bind(C,name='Ifpack_toString')
    import :: c_ptr ,FT_EPrecType_E_t
    
    integer(FT_EPrecType_E_t)   ,intent(in)   ,value              :: precType
  end function


  !> <BR> Original C++ prototype:
  !! static Ifpack_Preconditioner* Create( EPrecType PrecType, Epetra_RowMatrix* Matrix, const int overlap = 0 );
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Ifpack_Preconditioner_ID_t Ifpack_CreatePreconditioner_UsingType ( CT_EPrecType_E_t PrecType, CT_Epetra_RowMatrix_ID_t MatrixID, const int overlap );

  type(FT_Ifpack_Preconditioner_ID_t) function Ifpack_CreatePreconditioner_UsingType ( &
        PrecType, MatrixID, overlap ) bind(C,name='Ifpack_CreatePreconditioner_UsingType')
    import :: FT_Ifpack_Preconditioner_ID_t ,FT_EPrecType_E_t ,FT_Epetra_RowMatrix_ID_t , &
          c_int
    
    integer(FT_EPrecType_E_t)   ,intent(in)   ,value              :: PrecType
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: MatrixID
    integer(c_int)              ,intent(in)   ,value              :: overlap
  end function


  !> <BR> Original C++ prototype:
  !! Ifpack_Preconditioner* Create(const string PrecType, Epetra_RowMatrix* Matrix, const int overlap = 0);
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Ifpack_Preconditioner_ID_t Ifpack_CreatePreconditioner_UsingName ( CT_Ifpack_ID_t selfID, const char PrecType[], CT_Epetra_RowMatrix_ID_t MatrixID, const int overlap );

  type(FT_Ifpack_Preconditioner_ID_t) function Ifpack_CreatePreconditioner_UsingName ( &
        selfID, PrecType, MatrixID, overlap ) &
        bind(C,name='Ifpack_CreatePreconditioner_UsingName')
    import :: FT_Ifpack_Preconditioner_ID_t ,FT_Ifpack_ID_t ,c_char , &
          FT_Epetra_RowMatrix_ID_t ,c_int
    
    type(FT_Ifpack_ID_t)        ,intent(in)   ,value              :: selfID
    character(kind=c_char)      ,intent(in)         ,dimension(*) :: PrecType
    type(FT_Epetra_RowMatrix_ID_t),intent(in)   ,value              :: MatrixID
    integer(c_int)              ,intent(in)   ,value              :: overlap
  end function


  !> <BR> Original C++ prototype:
  !! int SetParameters(int argc, char* argv[], Teuchos::ParameterList& List, string& PrecType, int& Overlap);
  !> <BR> <BR> CTrilinos prototype:
  !! int Ifpack_SetParameters ( CT_Ifpack_ID_t selfID, int argc, char * argv[], CT_Teuchos_ParameterList_ID_t ListID, char * PrecType[], int * Overlap );

  integer(c_int) function Ifpack_SetParameters ( selfID, argc, argv, ListID, PrecType, &
        Overlap ) bind(C,name='Ifpack_SetParameters')
    import :: c_int ,FT_Ifpack_ID_t ,c_char ,FT_Teuchos_ParameterList_ID_t
    
    type(FT_Ifpack_ID_t)        ,intent(in)   ,value              :: selfID
    integer(c_int)              ,intent(in)   ,value              :: argc
    character(kind=c_char)                          ,dimension(*) :: argv
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: ListID
    character(kind=c_char)                          ,dimension(*) :: PrecType
    integer(c_int)              ,intent(inout)                    :: Overlap
  end function


!> @}


!> @name Ifpack_Preconditioner interface
!! @{

  ! _________________ Ifpack_Preconditioner interface bodies _________________


  !> <BR> CTrilinos prototype:
  !! CT_Ifpack_Preconditioner_ID_t Ifpack_Preconditioner_Degeneralize ( CTrilinos_Universal_ID_t id );

  type(FT_Ifpack_Preconditioner_ID_t) function Ifpack_Preconditioner_Degeneralize ( id ) &
        bind(C,name='Ifpack_Preconditioner_Degeneralize')
    import :: FT_Ifpack_Preconditioner_ID_t ,ForTrilinos_Universal_ID_t
    
    type(ForTrilinos_Universal_ID_t)   ,intent(in)   ,value              :: id
  end function


  !> <BR> CTrilinos prototype:
  !! CTrilinos_Universal_ID_t Ifpack_Preconditioner_Generalize ( CT_Ifpack_Preconditioner_ID_t id );

  type(ForTrilinos_Universal_ID_t) function Ifpack_Preconditioner_Generalize ( id ) &
        bind(C,name='Ifpack_Preconditioner_Generalize')
    import :: ForTrilinos_Universal_ID_t ,FT_Ifpack_Preconditioner_ID_t
    
    type(FT_Ifpack_Preconditioner_ID_t),intent(in)   ,value              :: id
  end function


  !> <BR> Original C++ prototype:
  !! virtual int SetParameters(Teuchos::ParameterList& List) = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Ifpack_Preconditioner_SetParameters ( CT_Ifpack_Preconditioner_ID_t selfID, CT_Teuchos_ParameterList_ID_t ListID );

  integer(c_int) function Ifpack_Preconditioner_SetParameters ( selfID, ListID ) &
        bind(C,name='Ifpack_Preconditioner_SetParameters')
    import :: c_int ,FT_Ifpack_Preconditioner_ID_t ,FT_Teuchos_ParameterList_ID_t
    
    type(FT_Ifpack_Preconditioner_ID_t),intent(in)   ,value              :: selfID
    type(FT_Teuchos_ParameterList_ID_t),intent(in)   ,value              :: ListID
  end function


  !> <BR> Original C++ prototype:
  !! virtual int Initialize() = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Ifpack_Preconditioner_Initialize ( CT_Ifpack_Preconditioner_ID_t selfID );

  integer(c_int) function Ifpack_Preconditioner_Initialize ( selfID ) &
        bind(C,name='Ifpack_Preconditioner_Initialize')
    import :: c_int ,FT_Ifpack_Preconditioner_ID_t
    
    type(FT_Ifpack_Preconditioner_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual bool IsInitialized() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Ifpack_Preconditioner_IsInitialized ( CT_Ifpack_Preconditioner_ID_t selfID );

  integer(FT_boolean_t) function Ifpack_Preconditioner_IsInitialized ( selfID ) &
        bind(C,name='Ifpack_Preconditioner_IsInitialized')
    import :: FT_boolean_t ,FT_Ifpack_Preconditioner_ID_t
    
    type(FT_Ifpack_Preconditioner_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual int Compute() = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Ifpack_Preconditioner_Compute ( CT_Ifpack_Preconditioner_ID_t selfID );

  integer(c_int) function Ifpack_Preconditioner_Compute ( selfID ) &
        bind(C,name='Ifpack_Preconditioner_Compute')
    import :: c_int ,FT_Ifpack_Preconditioner_ID_t
    
    type(FT_Ifpack_Preconditioner_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual bool IsComputed() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! boolean Ifpack_Preconditioner_IsComputed ( CT_Ifpack_Preconditioner_ID_t selfID );

  integer(FT_boolean_t) function Ifpack_Preconditioner_IsComputed ( selfID ) &
        bind(C,name='Ifpack_Preconditioner_IsComputed')
    import :: FT_boolean_t ,FT_Ifpack_Preconditioner_ID_t
    
    type(FT_Ifpack_Preconditioner_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual double Condest() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! double Ifpack_Preconditioner_Condest ( CT_Ifpack_Preconditioner_ID_t selfID );

  real(c_double) function Ifpack_Preconditioner_Condest ( selfID ) &
        bind(C,name='Ifpack_Preconditioner_Condest')
    import :: c_double ,FT_Ifpack_Preconditioner_ID_t
    
    type(FT_Ifpack_Preconditioner_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Ifpack_Preconditioner_ApplyInverse ( CT_Ifpack_Preconditioner_ID_t selfID, CT_Epetra_MultiVector_ID_t XID, CT_Epetra_MultiVector_ID_t YID );

  integer(c_int) function Ifpack_Preconditioner_ApplyInverse ( selfID, XID, YID ) &
        bind(C,name='Ifpack_Preconditioner_ApplyInverse')
    import :: c_int ,FT_Ifpack_Preconditioner_ID_t ,FT_Epetra_MultiVector_ID_t
    
    type(FT_Ifpack_Preconditioner_ID_t),intent(in)   ,value              :: selfID
    type(FT_Epetra_MultiVector_ID_t)   ,intent(in)   ,value              :: XID
    type(FT_Epetra_MultiVector_ID_t)   ,intent(in)   ,value              :: YID
  end function


  !> <BR> Original C++ prototype:
  !! virtual const Epetra_RowMatrix& Matrix() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! CT_Epetra_RowMatrix_ID_t Ifpack_Preconditioner_Matrix ( CT_Ifpack_Preconditioner_ID_t selfID );

  type(FT_Epetra_RowMatrix_ID_t) function Ifpack_Preconditioner_Matrix ( selfID ) &
        bind(C,name='Ifpack_Preconditioner_Matrix')
    import :: FT_Epetra_RowMatrix_ID_t ,FT_Ifpack_Preconditioner_ID_t
    
    type(FT_Ifpack_Preconditioner_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual int NumInitialize() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Ifpack_Preconditioner_NumInitialize ( CT_Ifpack_Preconditioner_ID_t selfID );

  integer(c_int) function Ifpack_Preconditioner_NumInitialize ( selfID ) &
        bind(C,name='Ifpack_Preconditioner_NumInitialize')
    import :: c_int ,FT_Ifpack_Preconditioner_ID_t
    
    type(FT_Ifpack_Preconditioner_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual int NumCompute() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Ifpack_Preconditioner_NumCompute ( CT_Ifpack_Preconditioner_ID_t selfID );

  integer(c_int) function Ifpack_Preconditioner_NumCompute ( selfID ) &
        bind(C,name='Ifpack_Preconditioner_NumCompute')
    import :: c_int ,FT_Ifpack_Preconditioner_ID_t
    
    type(FT_Ifpack_Preconditioner_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual int NumApplyInverse() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! int Ifpack_Preconditioner_NumApplyInverse ( CT_Ifpack_Preconditioner_ID_t selfID );

  integer(c_int) function Ifpack_Preconditioner_NumApplyInverse ( selfID ) &
        bind(C,name='Ifpack_Preconditioner_NumApplyInverse')
    import :: c_int ,FT_Ifpack_Preconditioner_ID_t
    
    type(FT_Ifpack_Preconditioner_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual double InitializeTime() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! double Ifpack_Preconditioner_InitializeTime ( CT_Ifpack_Preconditioner_ID_t selfID );

  real(c_double) function Ifpack_Preconditioner_InitializeTime ( selfID ) &
        bind(C,name='Ifpack_Preconditioner_InitializeTime')
    import :: c_double ,FT_Ifpack_Preconditioner_ID_t
    
    type(FT_Ifpack_Preconditioner_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual double ComputeTime() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! double Ifpack_Preconditioner_ComputeTime ( CT_Ifpack_Preconditioner_ID_t selfID );

  real(c_double) function Ifpack_Preconditioner_ComputeTime ( selfID ) &
        bind(C,name='Ifpack_Preconditioner_ComputeTime')
    import :: c_double ,FT_Ifpack_Preconditioner_ID_t
    
    type(FT_Ifpack_Preconditioner_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual double ApplyInverseTime() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! double Ifpack_Preconditioner_ApplyInverseTime ( CT_Ifpack_Preconditioner_ID_t selfID );

  real(c_double) function Ifpack_Preconditioner_ApplyInverseTime ( selfID ) &
        bind(C,name='Ifpack_Preconditioner_ApplyInverseTime')
    import :: c_double ,FT_Ifpack_Preconditioner_ID_t
    
    type(FT_Ifpack_Preconditioner_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual double InitializeFlops() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! double Ifpack_Preconditioner_InitializeFlops ( CT_Ifpack_Preconditioner_ID_t selfID );

  real(c_double) function Ifpack_Preconditioner_InitializeFlops ( selfID ) &
        bind(C,name='Ifpack_Preconditioner_InitializeFlops')
    import :: c_double ,FT_Ifpack_Preconditioner_ID_t
    
    type(FT_Ifpack_Preconditioner_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual double ComputeFlops() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! double Ifpack_Preconditioner_ComputeFlops ( CT_Ifpack_Preconditioner_ID_t selfID );

  real(c_double) function Ifpack_Preconditioner_ComputeFlops ( selfID ) &
        bind(C,name='Ifpack_Preconditioner_ComputeFlops')
    import :: c_double ,FT_Ifpack_Preconditioner_ID_t
    
    type(FT_Ifpack_Preconditioner_ID_t),intent(in)   ,value              :: selfID
  end function


  !> <BR> Original C++ prototype:
  !! virtual double ApplyInverseFlops() const = 0;
  !> <BR> <BR> CTrilinos prototype:
  !! double Ifpack_Preconditioner_ApplyInverseFlops ( CT_Ifpack_Preconditioner_ID_t selfID );

  real(c_double) function Ifpack_Preconditioner_ApplyInverseFlops ( selfID ) &
        bind(C,name='Ifpack_Preconditioner_ApplyInverseFlops')
    import :: c_double ,FT_Ifpack_Preconditioner_ID_t
    
    type(FT_Ifpack_Preconditioner_ID_t),intent(in)   ,value              :: selfID
  end function


!> @}


  end interface
end module forifpack

#endif /* HAVE_FORTRILINOS_IFPACK */
