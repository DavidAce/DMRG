function(lapacke_message TYPE MSG)
    if(Lapacke_FIND_VERBOSE)
        message(${TYPE} ${MSG})
    endif()
endfunction()

function(CheckLapackeCompiles TAG REQUIRED_FLAGS REQUIRED_DEFINITIONS REQUIRED_LIBRARIES_UNPARSED REQUIRED_INCLUDES REQUIRED_TARGET)
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
    if(DMRG_PRINT_CHECKS OR Lapacke_FIND_VERBOSE)
        message(STATUS "LAPACKE TEST COMPILE CMAKE_REQUIRED_FLAGS        ${CMAKE_REQUIRED_FLAGS}")
        message(STATUS "LAPACKE TEST COMPILE CMAKE_REQUIRED_DEFINITIONS  ${CMAKE_REQUIRED_DEFINITIONS}")
        message(STATUS "LAPACKE TEST COMPILE CMAKE_REQUIRED_LIBRARIES    ${CMAKE_REQUIRED_LIBRARIES}")
        message(STATUS "LAPACKE TEST COMPILE CMAKE_REQUIRED_INCLUDES     ${CMAKE_REQUIRED_INCLUDES}")
    endif()
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



function(find_Lapacke)

    if (NOT TARGET lapacke::lapacke)
        # Find from MKL
        if(TARGET mkl::mkl)
            # Try finding lapacke in MKL library
#            if(Lapacke_FIND_VERBOSE)
                lapacke_message(STATUS "Searching for Lapacke in Intel MKL")
#            endif()
            if(MKL_INCLUDE_DIR)
                CheckLapackeCompiles("MKL" ""  "-DMKL_AVAILABLE"  "" "" "mkl::mkl")
            endif()

            if(LAPACKE_COMPILES_MKL)
                add_library(lapacke::lapacke INTERFACE IMPORTED)
                get_target_property(LAPACKE_INCLUDE_DIR mkl::mkl INTERFACE_SYSTEM_INCLUDE_DIRECTORIES)
                target_link_libraries(lapacke::lapacke INTERFACE mkl::mkl)
                target_include_directories(lapacke::lapacke SYSTEM INTERFACE ${LAPACKE_INCLUDE_DIR})
                lapacke_message(STATUS "Searching for Lapacke in Intel MKL - Success")
            else()
                lapacke_message(STATUS "Searching for Lapacke in Intel MKL - failed")
            endif()
        endif()
    endif()




    if (NOT TARGET lapacke::lapacke)
        if (TARGET openblas::openblas)
            lapacke_message(STATUS "Searching for Lapacke in OpenBLAS")
            if(LAPACKE_DEBUG)
                include(cmake-modules/PrintTargetProperties.cmake)
                print_target_properties(openblas::openblas)
            endif()
            CheckLapackeCompiles("OpenBLAS" "" "" "" "" "openblas::openblas")
            if(LAPACKE_COMPILES_OpenBLAS)
                add_library(lapacke::lapacke INTERFACE IMPORTED)
                get_target_property(LAPACKE_INCLUDE_DIR openblas::openblas INTERFACE_SYSTEM_INCLUDE_DIRECTORIES)
                target_link_libraries(lapacke::lapacke INTERFACE openblas::openblas)
                target_include_directories(lapacke::lapacke SYSTEM INTERFACE ${LAPACKE_INCLUDE_DIR})
                lapacke_message(STATUS "Searching for Lapacke in OpenBLAS - Success")
            else()
                lapacke_message(STATUS "Searching for Lapacke in OpenBLAS - failed")
            endif()
        endif()
    endif()



    if (NOT TARGET lapacke::lapacke)
        if(TARGET lapack::lapack)
            lapacke_message(STATUS "Searching for Lapacke in system")
            find_path(LAPACKE_INCLUDE_DIR
                    NAMES lapacke.h
                    HINTS ${CMAKE_INSTALL_PREFIX} $ENV{EBROOTOPENBLAS} ${CONDA_HINTS}
                    PATHS $ENV{EBROOTBLAS} $ENV{BLAS_DIR} $ENV{BLAS_ROOT}
                    PATH_SUFFIXES
                        OpenBLAS openblas openblas/include OpenBLAS/include lapack)
            find_library(LAPACKE_LIBRARY
                    NAMES lapacke
                    HINTS ${CMAKE_INSTALL_PREFIX} $ENV{EBROOTOPENBLAS} ${CONDA_HINTS}
                    PATHS $ENV{EBROOTBLAS} $ENV{BLAS_DIR} $ENV{BLAS_ROOT}
                    PATH_SUFFIXES
                        OpenBLAS openblas openblas/lib OpenBLAS/lib lapack lapack/lib blas blas/lib
                    )
            if(LAPACKE_INCLUDE_DIR AND LAPACKE_LIBRARY)
                CheckLapackeCompiles("lib_header" " "   " "
                        "${LAPACKE_LIBRARY}"
                        "${LAPACKE_INCLUDE_DIR}"
                        "lapack::lapack"
                        )
                if(LAPACKE_COMPILES_lib_header)
                    add_library(lapacke::lapacke ${LINK_TYPE} IMPORTED)
                    set_target_properties(lapacke::lapacke PROPERTIES IMPORTED_LOCATION "${LAPACKE_LIBRARY}")
                    target_include_directories(lapacke::lapacke SYSTEM INTERFACE ${LAPACKE_INCLUDE_DIR})
                endif()
            endif()
            if(NOT TARGET lapacke::lapacke AND LAPACKE_INCLUDE_DIR)
                CheckLapackeCompiles("header" " "   " "
                        ""
                        "${LAPACKE_INCLUDE_DIR}"
                        "lapack::lapack"
                        )
                if(LAPACKE_COMPILES_header)
                    add_library(lapacke::lapacke INTERFACE IMPORTED)
                    target_include_directories(lapacke::lapacke SYSTEM INTERFACE ${LAPACKE_INCLUDE_DIR})
                endif()
            endif()
            if(NOT TARGET lapacke::lapacke AND LAPACKE_LIBRARY)
                CheckLapackeCompiles("lib" " "   " "
                        "${LAPACKE_LIBRARY}"
                        ""
                        "lapack::lapack"
                        )
                if(LAPACKE_COMPILES_lib)
                    add_library(lapacke::lapacke ${LINK_TYPE} IMPORTED)
                    set_target_properties(lapacke::lapacke PROPERTIES IMPORTED_LOCATION "${LAPACKE_LIBRARY}")
                endif()
            endif()
            if(TARGET lapacke::lapacke)
                target_link_libraries(lapacke::lapacke INTERFACE blas::blas lapack::lapack gfortran::gfortran)
                lapacke_message(STATUS "Searching for Lapacke in system - Success: ${LAPACKE_LIBRARY}")
            else()
                lapacke_message(STATUS "Searching for Lapacke in system - failed")
                lapacke_message(STATUS "Tried LAPACKE_LIBRARY     : ${LAPACKE_LIBRARY} ")
                lapacke_message(STATUS "Tried LAPACKE_INCLUDE_DIR : ${LAPACKE_INCLUDE_DIR} ")
            endif()
        endif()
    endif()

endfunction()

find_Lapacke()
if(TARGET lapacke::lapacke)
    include(cmake-modules/PrintTargetInfo.cmake)
    print_target_info(lapacke::lapacke)
    get_target_property(LAPACKE_INCLUDE_DIR lapacke::lapacke INTERFACE_SYSTEM_INCLUDE_DIRECTORIES)
    message("LAPACKE_INCLUDE_DIR ${LAPACKE_INCLUDE_DIR}")
    set(Lapacke_TARGET      lapacke::lapacke)
    set(Lapacke_FOUND       TRUE)
else()
    set(Lapacke_FOUND       FALSE)
endif()



include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
        Lapacke
        FOUND_VAR Lapacke_FOUND
        REQUIRED_VARS LAPACKE_INCLUDE_DIR Lapacke_TARGET
        HANDLE_COMPONENTS
        FAIL_MESSAGE "Failed to find Lapacke"
)