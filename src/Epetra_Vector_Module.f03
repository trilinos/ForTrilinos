!module Epetra_Vector_Module
!
!  use ,intrinsic :: iso_c_binding ,only : c_int
!  implicit none
!  private
!
!  public :: EPetra_Vector
!
!  type EPetra_Vector 
!    private
!    integer(c_int) :: id
!  end type EPetra_Vector
!
!contains
!
!  function create(map) result(
!    integer :: map
!    b = FEpetra_Vector_Create(map)
!  end function
!
!  call PutScalar(b, two)
!
!  B = FEpetra_CrsMatrix_Create(...)
!  class PutScalar(B, two)
!
!module FEpetra_Vector_Module

