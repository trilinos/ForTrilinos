module fortrilinos_utils
implicit none
contains

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

  type(ForTrilinos_Object_ID_t) function generalize_all(object_id) bind(C,name="for_linking_only") 
    use ForTrilinos_enums ,only : ForTrilinos_Object_ID_t
    use ,intrinsic :: iso_c_binding ,only: c_ptr,c_f_pointer
    type(c_ptr) ,value :: object_id
    type(ForTrilinos_Object_ID_t), pointer :: local_ptr

    call c_f_pointer (object_id, local_ptr)
    generalize_all = local_ptr
  end function

end module
