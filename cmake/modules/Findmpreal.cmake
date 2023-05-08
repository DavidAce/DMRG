
function(find_mpreal)
    find_path(MPREAL_INCLUDE_DIR
                    NAMES mpreal.h
                    PATH_SUFFIXES
                        include
                        include/mpreal
                        include/mpfrc++
                        include/mpfrc
                        include/mpfr
                    )   
endfunction()

find_package(mpfr REQUIRED)
find_mpreal()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(mpreal DEFAULT_MSG MPREAL_INCLUDE_DIR)

if(mpreal_FOUND)
    if(NOT TARGET mpreal::mpreal)
        add_library(mpreal::mpreal INTERFACE IMPORTED)
    endif()

    target_include_directories(mpreal::mpreal SYSTEM INTERFACE ${MPREAL_INCLUDE_DIR})
    target_link_libraries(mpreal::mpreal INTERFACE mpfr::mpfr)
endif()