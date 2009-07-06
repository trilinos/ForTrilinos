interface 
  integer(c_int) function string_length(ptr) bind(C,name="forLinkingOnly") 
    use iso_c_binding, only: c_ptr, c_int 
    implicit none 
    type(c_ptr), value :: ptr 
  end function 
end interface 
