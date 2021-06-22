function(check_lapacke_compiles TARGETS LIBS INCS OPTS DEFS)
    list(APPEND CMAKE_REQUIRED_LIBRARIES     ${LIBS} ${TARGETS})
    list(APPEND CMAKE_REQUIRED_INCLUDES      ${INCS})
    list(APPEND CMAKE_REQUIRED_FLAGS         ${OPTS})
    list(APPEND CMAKE_REQUIRED_DEFINITIONS   ${DEFS})

    list(TRANSFORM "CMAKE_REQUIRED_DEFINITIONS" PREPEND "-D")  # Definitions should start with "-D"
    string(REPLACE ";" " "  CMAKE_REQUIRED_FLAGS  "${CMAKE_REQUIRED_FLAGS}")        # Needs to be a space-separated list

    if(DMRG_PRINT_CHECKS)
        message(STATUS "LAPACKE COMPILE TEST CMAKE_REQUIRED_LIBRARIES    ${CMAKE_REQUIRED_LIBRARIES}")
        message(STATUS "LAPACKE COMPILE TEST CMAKE_REQUIRED_INCLUDES     ${CMAKE_REQUIRED_INCLUDES}")
        message(STATUS "LAPACKE COMPILE TEST CMAKE_REQUIRED_FLAGS        ${CMAKE_REQUIRED_FLAGS}")
        message(STATUS "LAPACKE COMPILE TEST CMAKE_REQUIRED_DEFINITIONS  ${CMAKE_REQUIRED_DEFINITIONS}")
    endif()
    #   Test features
    include(CheckCXXSourceCompiles)
    check_cxx_source_compiles("
        #include <complex>
        #ifndef lapack_complex_float
            #define lapack_complex_float  std::complex<float>
        #endif
        #ifndef lapack_complex_double
            #define lapack_complex_double std::complex<double>
        #endif
        #if __has_include(<mkl_lapacke.h>)
        #include <mkl_lapacke.h>
        #elif __has_include(<openblas/lapacke.h>)
        #include <openblas/lapacke.h>
        #else
        #include <lapacke.h>
        #endif

        int main (int argc, const char * argv[])
        {
           double a[5][3] = {{1,1,1},{2,3,4},{3,5,2},{4,2,5},{5,4,3}};
           double b[5][2] = {{-10,-3},{12,14},{14,12},{16,16},{18,16}};
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
    set(LAPACKE_COMPILES ${LAPACKE_COMPILES} PARENT_SCOPE)
    if(NOT LAPACKE_COMPILES)
        unset(LAPACKE_COMPILES CACHE)
        unset(LAPACKE_COMPILES PARENT_SCOPE)
        if(DMRG_PRINT_CHECKS AND EXISTS "${CMAKE_BINARY_DIR}/CMakeFiles/CMakeError.log")
            file(READ "${CMAKE_BINARY_DIR}/CMakeFiles/CMakeError.log" ERROR_LOG)
            message(STATUS ${ERROR_LOG})
        endif()
    endif()
endfunction()
