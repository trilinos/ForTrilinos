interface 
  function return_string_len(ptr) bind(C,name="forLinkingOnly") 
    use iso_c_binding, only: c_ptr, c_int 
    implicit none 
    integer(c_int) return_string_len 
    type(c_ptr), value :: ptr 
  end function 
end interface 
