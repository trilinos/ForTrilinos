module fortrilinos_utils
implicit none
contains
  integer(c_int) function count(ptr) bind(C,name="forLinkingOnly") 
    use ,intrinsic :: iso_c_binding ,only: c_char ,c_int 
    implicit none 
    character(kind=c_char) :: ptr(*) 
    count = 0 
    do 
       if(ptr(count+1) == achar(0)) return 
       count = count+1 
    end do 
  end function 
end module
