module Epetra_Vector_Module

  use ,intrinsic :: iso_c_binding     ,only : c_int ,c_double
  use            :: epetra_map_module ,only : Epetra_Map 
  use            :: forepetra
  implicit none
  private

  type ,public :: EPetra_Vector 
    private
    integer(c_int) :: id
  contains
    procedure :: Create
    procedure :: PutScalar
    procedure :: Update
    procedure :: Norm2
    procedure :: Destroy
  end type 

contains

  subroutine Create(vector,map) 
    type(Epetra_Vector) ,intent(out) :: vector
    type(Epetra_Map)    ,intent(in)  :: map
    vector%id = FEpetra_Vector_Create(map%GetMapID())
  end subroutine

  subroutine PutScalar(vector, scalar) 
    type(Epetra_Vector) ,intent(in) :: vector
    real(c_double)      ,intent(in) :: scalar
    call FEpetra_Vector_PutScalar(vector%id,scalar) 
  end subroutine 

  subroutine Update(x, xScalar, b, bScalar) 
    type(Epetra_Vector) ,intent(in) :: x,b
    real(c_double)      ,intent(in) :: xScalar,bScalar
    call FEpetra_Vector_Update(x%id, xScalar, b%id, bScalar) 
  end subroutine 

  function Norm2(vector) 
    type(Epetra_Vector) ,intent(in) :: vector
    real(c_double)                  :: Norm2
    Norm2 = FEpetra_Vector_Norm2(vector%id)
  end function 

  subroutine Destroy(vector) 
    type(Epetra_Vector) :: vector
    call FEpetra_Vector_Destroy(vector%id)
  end subroutine 
end module Epetra_Vector_Module
