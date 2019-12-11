

function(CheckLapackeCompiles TAG REQUIRED_FLAGS REQUIRED_DEFINITIONS REQUIRED_LIBRARIES_UNPARSED REQUIRED_INCLUDES REQUIRED_TARGET)
#    message(STATUS "Checking if lapacke headers are working")
    set(REQUIRED_LIBRARIES)
    include(cmake-modules/getExpandedTarget.cmake)
    foreach(elem ${REQUIRED_LIBRARIES_UNPARSED})
        if(TARGET ${elem})
            expand_target_libs(${elem} expanded_libs)
            list(APPEND REQUIRED_LIBRARIES ${expanded_libs})
        else()
            list(APPEND REQUIRED_LIBRARIES ${elem})
        endif()
    endforeach()

    foreach(elem ${REQUIRED_TARGET})
        if(TARGET ${elem})
            expand_target_libs(${elem} expanded_libs)
            expand_target_incs(${elem} expanded_incs)
            expand_target_opts(${elem} expanded_opts)
            list(APPEND REQUIRED_LIBRARIES ${expanded_libs})
            list(APPEND REQUIRED_INCLUDES  ${expanded_incs})
            list(APPEND REQUIRED_FLAGS     ${expanded_opts})
        endif()
    endforeach()
    string(REPLACE ";" " " REQUIRED_FLAGS      "${REQUIRED_FLAGS}") # Needs to be a space-separated list


#   Test features
    include(CheckCXXSourceCompiles)
    set(CMAKE_REQUIRED_FLAGS        ${REQUIRED_FLAGS})
    set(CMAKE_REQUIRED_DEFINITIONS  ${REQUIRED_DEFINITIONS})
    set(CMAKE_REQUIRED_LIBRARIES    ${REQUIRED_LIBRARIES})
    set(CMAKE_REQUIRED_INCLUDES     ${REQUIRED_INCLUDES})
    message(STATUS "LAPACKE TEST COMPILE CMAKE_REQUIRED_FLAGS        ${CMAKE_REQUIRED_FLAGS}")
    message(STATUS "LAPACKE TEST COMPILE CMAKE_REQUIRED_DEFINITIONS  ${CMAKE_REQUIRED_DEFINITIONS}")
    message(STATUS "LAPACKE TEST COMPILE CMAKE_REQUIRED_LIBRARIES    ${CMAKE_REQUIRED_LIBRARIES}")
    message(STATUS "LAPACKE TEST COMPILE CMAKE_REQUIRED_INCLUDES     ${CMAKE_REQUIRED_INCLUDES}")

    check_cxx_source_compiles("
        #ifdef MKL_AVAILABLE
        #include <mkl_lapacke.h>
        #else
        #include <lapacke.h>
        #endif

        int main (int argc, const char * argv[])
        {
           double a[5][3] = {1,1,1,2,3,4,3,5,2,4,2,5,5,4,3};
           double b[5][2] = {-10,-3,12,14,14,12,16,16,18,16};
           lapack_int info,m,n,lda,ldb,nrhs;
           int i,j;
           m = 5;
           n = 3;
           nrhs = 2;
           lda = 3;
           ldb = 2;
           info = LAPACKE_dgels(LAPACK_ROW_MAJOR,'N',m,n,nrhs,*a,lda,*b,ldb);
           return(info);
        }
        " LAPACKE_COMPILES_${TAG})
    if(LAPACKE_COMPILES_${TAG})
        set(LAPACKE_COMPILES_${TAG} TRUE PARENT_SCOPE)
    else()
        set(LAPACKE_COMPILES_${TAG} FALSE PARENT_SCOPE)
    endif()
endfunction()




if (NOT TARGET lapacke)
    # Find from MKL
    if(TARGET mkl)
        # Try finding lapacke in MKL library
        message(STATUS "Searching for Lapacke in Intel MKL.")
        if(MKL_INCLUDE_DIR)
            CheckLapackeCompiles("MKL" ""  "-DMKL_AVAILABLE"  "" "" "mkl")
        endif()

        if(LAPACKE_COMPILES_MKL)
            add_library(lapacke INTERFACE)
            target_link_libraries(lapacke INTERFACE mkl)
            message(STATUS "Searching for Lapacke in Intel MKL - Success")
        else()
            message(STATUS "Searching for Lapacke in Intel MKL - failed")
        endif()
    endif()
endif()




if (NOT TARGET lapacke)
    if (TARGET OpenBLAS)
        message(STATUS "Searching for Lapacke in OpenBLAS")
        CheckLapackeCompiles("OpenBLAS" "" "" "" "" "OpenBLAS")
        if(LAPACKE_COMPILES_OpenBLAS)
            add_library(lapacke INTERFACE)
            target_link_libraries(lapacke INTERFACE OpenBLAS)
            message(STATUS "Searching for Lapacke in OpenBLAS - Success")
        else()
            message(STATUS "Searching for Lapacke in OpenBLAS - failed")
        endif()
    endif()
endif()



if (NOT TARGET lapacke)
    if(TARGET lapack)
        message(STATUS "Searching for Lapacke in system")
        find_path(LAPACKE_INCLUDE_DIR
                NAMES lapacke.h
                HINTS ${DIRECTORY_HINTS}
                PATHS
                    $ENV{EBROOTOPENBLAS}
                    $ENV{EBROOTBLAS}
                    $ENV{BLAS_DIR}    ${BLAS_DIR}
                    $ENV{blas_DIR}    ${blas_DIR}
                    $ENV{LAPACKE_DIR} ${LAPACKE_DIR}
                    $ENV{lapacke_DIR} ${lapacke_DIR}
                PATH_SUFFIXES
                    OpenBLAS openblas openblas/include OpenBLAS/include lapack
                )
        find_library(LAPACKE_LIBRARY
                NAMES liblapacke${CUSTOM_SUFFIX}
                HINTS ${DIRECTORY_HINTS}
                PATHS
                    $ENV{EBROOTOPENBLAS}
                    $ENV{EBROOTBLAS}
                    $ENV{BLAS_DIR}    ${BLAS_DIR}
                    $ENV{blas_DIR}    ${blas_DIR}
                    $ENV{LAPACKE_DIR} ${LAPACKE_DIR}
                    $ENV{lapacke_DIR} ${lapacke_DIR}
                PATH_SUFFIXES
                    OpenBLAS openblas openblas/lib OpenBLAS/lib lapack lapack/lib blas blas/lib
                )
        if(LAPACKE_INCLUDE_DIR AND LAPACKE_LIBRARY)
            CheckLapackeCompiles("SYSTEM" " "   " "
                    "${LAPACKE_LIBRARY}"
                    "${LAPACKE_INCLUDE_DIR}"
                    "lapack"
                    )
        if(NOT LAPACKE_COMPILES_SYSTEM AND LAPACKE_INCLUDE_DIR)
            CheckLapackeCompiles("SYSTEM" " "   " "
                    ""
                    "${LAPACKE_INCLUDE_DIR}"
                    "lapack"
                    )
        if(NOT LAPACKE_COMPILES_SYSTEM AND LAPACKE_LIBRARY)
            CheckLapackeCompiles("SYSTEM" " "   " "
                    "${LAPACKE_LIBRARY}"
                    ""
                    "lapack"
                    )
        endif()
        if(LAPACKE_COMPILES_SYSTEM)
            add_library(lapacke INTERFACE IMPORTED)
            if(LAPACKE_LIBRARY)
                target_link_libraries(lapacke INTERFACE ${LAPACKE_LIBRARY} lapack)
            endif()
            if(LAPACKE_INCLUDE_DIR)
                target_include_directories(lapacke SYSTEM INTERFACE ${LAPACKE_INCLUDE_DIR})
            endif()
            target_link_libraries(lapacke INTERFACE lapack)
            message(STATUS "Searching for Lapacke in system - Success: ${LAPACKE_LIBRARY}")
        else()
            message(STATUS "Searching for Lapacke in system - failed")
        endif()
    endif()
endif()



if(NOT TARGET lapacke)
    message(FATAL_ERROR "Library LAPACKE could not be found.")
endif()