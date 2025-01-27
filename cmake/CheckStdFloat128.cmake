function(check_std_float128_t)
    set(CMAKE_EXTRA_INCLUDE_FILES "stdfloat")
    set(CMAKE_REQUIRED_FLAGS "-std=c++23")
    include(CheckTypeSize)
    check_type_size("std::float128_t" SIZEOF_STDFLOAT128 LANGUAGE CXX)
    if(NOT HAVE_SIZEOF_STDFLOAT128)
        message(FATAL_ERROR "DMRG_USE_FLOAT128:${DMRG_USE_FLOAT128} could not be satisfied: std::float128_t not found.")
    endif()
endfunction()

