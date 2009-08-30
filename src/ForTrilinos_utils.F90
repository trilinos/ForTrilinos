module fortrilinos_utils
implicit none
contains
  
  ! Fortran-style strings take two forms: one with no C counterpart and a second that is
  ! an array of character values.  A C string is interoperable with the second type of 
  ! Fortran string with the primary distinction being that Fortran strings carry their
  ! own string-length information, whereas the length of C strings is indicated by a 
  ! terminating null character. The "count" procedure calculates C string length by 
  ! searching for the null terminator.

  integer(c_int) function count(ptr) bind(C,name="forLinkingOnly") 
    use ,intrinsic :: iso_c_binding ,only: c_char ,c_int  ,c_null_char
    implicit none 
    character(kind=c_char) :: ptr(*) 
    count = 0 
    do 
       if(ptr(count+1) == c_null_char) return 
       count = count+1 
    end do 
  end function 

  ! This is a Fortran implementation of the functionality in the Epetra_*_Abstract procedures
  ! in each CTrilinos/src/CEpetra* file.  It effectively casts any given Epetra derived type
  ! to a general type that can represent any of the Epetra derived types in 
  ! ForTrilinos/src/ForTrilinos_enums.F90.

  type(ForTrilinos_Object_ID_t) function generalize_all(object_id) bind(C,name="for_linking_only") 
    use ForTrilinos_enums ,only : ForTrilinos_Object_ID_t
    use ,intrinsic :: iso_c_binding ,only: c_ptr,c_f_pointer
    type(c_ptr) ,value :: object_id
    type(ForTrilinos_Object_ID_t), pointer :: local_ptr

    call c_f_pointer (object_id, local_ptr)
    generalize_all = local_ptr
  end function
end module
