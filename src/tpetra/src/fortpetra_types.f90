module fortpetra_types
  use, intrinsic :: iso_c_binding, only : &
    global_ordinal_type => c_int, &  ! c_long_long, &
    global_size_t => c_long, &
    local_ordinal_type => c_int, &
    size_type => c_size_t, &
    bool_type => c_bool, &
    int_type => c_int, &
    scalar_type => c_double, &
    mag_type => c_double, &
    norm_type => c_double
end module fortpetra_types

