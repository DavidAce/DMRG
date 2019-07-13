

function(CheckLapackeCompiles REQUIRED_FLAGS REQUIRED_DEFINITIONS REQUIRED_LIBRARIES REQUIRED_INCLUDES)
    message(STATUS "Checking for working lapacke headers")
    set(LAPACKE_LIBRARY ${LAPACKE_LAPACK_LIBRARY})
    #   Test features
    include(CheckCXXSourceCompiles)
    set(CMAKE_REQUIRED_FLAGS        ${REQUIRED_FLAGS})
    set(CMAKE_REQUIRED_DEFINITIONS  ${REQUIRED_DEFINITIONS})
    set(CMAKE_REQUIRED_LIBRARIES    ${REQUIRED_LIBRARIES})
    set(CMAKE_REQUIRED_INCLUDES     ${REQUIRED_INCLUDES})
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
        " LAPACKE_COMPILES)
    if(LAPACKE_COMPILES)
        message(STATUS "Checking for working lapacke headers - success")
        set(LAPACKE_COMPILES TRUE PARENT_SCOPE)
    else()
        message(STATUS "Lapacke does not work with the following parameters:")
        message(STATUS " -- flags        : ${REQUIRED_FLAGS}")
        message(STATUS " -- definitions  : ${REQUIRED_DEFINITIONS}")
        message(STATUS " -- libraries    : ${REQUIRED_LIBRARIES}")
        message(STATUS " -- includes     : ${REQUIRED_INCLUDES}")
        set(LAPACKE_COMPILES FALSE PARENT_SCOPE)
    endif()
endfunction()








if (NOT TARGET lapacke)
    # Find from MKL
    if(TARGET mkl)
        # Try finding lapacke in MKL library
        message(STATUS "Attempting to find LAPACKE in Intel MKL.")
        get_target_property(MKL_LIBRARIES mkl INTERFACE_LINK_LIBRARIES)
        get_target_property(MKL_INCLUDE_DIR mkl INTERFACE_INCLUDE_DIRECTORIES)
        if(MKL_INCLUDE_DIR)
            CheckLapackeCompiles(" "   "-DMKL_AVAILABLE"  "${MKL_LIBRARIES}" "${MKL_INCLUDE_DIR}")
        endif()
        if(LAPACKE_COMPILES)
            add_library(lapacke INTERFACE)
            add_dependencies(lapacke INTERFACE mkl)
            target_link_libraries(lapacke INTERFACE ${MKL_LIBRARIES})
            target_include_directories(lapacke INTERFACE ${MKL_INCLUDE_DIR})
            return()
        endif()
    endif()
endif()

if (NOT TARGET lapacke)
    if (TARGET lapack)
        message(STATUS "Attempting to find LAPACKE together with lapack (e.g. OpenBLAS)")
        get_target_property(LAPACKE_LAPACK_LIBRARY   lapack INTERFACE_LINK_LIBRARIES)
        get_target_property(LAPACKE_LAPACK_INCLUDE   lapack INTERFACE_INCLUDE_DIRECTORIES)
        find_path(LAPACKE_INCLUDE_DIR
                NAMES lapacke.h
                PATHS ${LAPACKE_LAPACK_INCLUDE}
                NO_DEFAULT_PATH)
        if(LAPACKE_LAPACK_LIBRARY AND LAPACKE_LAPACK_INCLUDE)
            CheckLapackeCompiles(" "   " "  "${LAPACKE_LAPACK_LIBRARY}" "${LAPACKE_INCLUDE_DIR}")
        endif()
        if(LAPACKE_COMPILES)
            add_library(lapacke INTERFACE)
            add_dependencies(lapacke INTERFACE lapack)
            target_link_libraries(lapacke INTERFACE lapack)
            target_include_directories(lapacke INTERFACE ${LAPACKE_INCLUDE_DIR})
            return()
        endif()
    endif()
endif()

if (NOT TARGET lapacke)
    if (TARGET lapack)
        message(STATUS "Attempting to find LAPACKE in system")
        find_path(LAPACKE_INCLUDE_DIR
                NAMES lapacke.h
                PATHS
                ${LAPACKE_LAPACK_INCLUDE}
                $ENV{BLAS_DIR}/include
                $ENV{HOME}/anaconda3/
                $ENV{HOME}/anaconda3/include
                $ENV{HOME}/.conda/include
                /usr/include
                /usr/include/x86_64-linux-gnu
                )
        find_library(LAPACKE_LIBRARY
                NAMES liblapacke${CUSTOM_SUFFIX}
                PATHS
                $ENV{BLAS_DIR}/lib
                $ENV{HOME}/.conda/lib
                $ENV{HOME}/anaconda3/
                $ENV{HOME}/anaconda3/lib
                /usr/lib/x86_64-linux-gnu
                )
        if(LAPACKE_INCLUDE_DIR OR LAPACKE_LIBRARY)
            CheckLapackeCompiles(" "   " "  "${LAPACKE_LIBRARY} " "${LAPACKE_INCLUDE_DIR} ")
        endif()
        if(LAPACKE_COMPILES)
            add_library(lapacke INTERFACE)
            add_dependencies(lapacke INTERFACE lapack)
            target_link_libraries(lapacke INTERFACE ${LAPACKE_LIBRARY})
            target_include_directories(lapacke INTERFACE ${LAPACKE_INCLUDE_DIR})
            return()
        endif()
    endif()
endif()



if(NOT TARGET lapacke)
    message(WARNING "Library LAPACKE could not be found.")
endif()