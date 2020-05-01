function(CheckLapackeCompiles TARGETS LIBS INCS OPTS DEFS)
    list(APPEND CMAKE_REQUIRED_LIBRARIES     ${LIBS})
    list(APPEND CMAKE_REQUIRED_INCLUDES      ${INCS})
    list(APPEND CMAKE_REQUIRED_FLAGS         ${OPTS})
    list(APPEND CMAKE_REQUIRED_DEFINITIONS   ${DEFS})
    include(cmake-modules/getExpandedTarget.cmake)
    foreach(elem ${TARGETS})
        if(TARGET ${elem})
            expand_target_libs(${elem} libs)
            expand_target_incs(${elem} incs)
            expand_target_opts(${elem} opts)
            expand_target_defs(${elem} defs)
            list(APPEND CMAKE_REQUIRED_LIBRARIES     ${libs})
            list(APPEND CMAKE_REQUIRED_INCLUDES      ${incs})
            list(APPEND CMAKE_REQUIRED_FLAGS         ${opts})
            list(APPEND CMAKE_REQUIRED_DEFINITIONS   ${defs})
        endif()
    endforeach()


    list(TRANSFORM "CMAKE_REQUIRED_DEFINITIONS" PREPEND "-D")  # Definitions should start with "-D"
    string(REPLACE ";" " "  CMAKE_REQUIRED_FLAGS          "${CMAKE_REQUIRED_FLAGS}")        # Needs to be a space-separated list

    include(CheckCXXSourceCompiles)
    if(DMRG_PRINT_CHECKS OR LAPACKE_FIND_VERBOSE)
        message(STATUS "LAPACKE TEST COMPILE CMAKE_REQUIRED_LIBRARIES    ${CMAKE_REQUIRED_LIBRARIES}")
        message(STATUS "LAPACKE TEST COMPILE CMAKE_REQUIRED_INCLUDES     ${CMAKE_REQUIRED_INCLUDES}")
        message(STATUS "LAPACKE TEST COMPILE CMAKE_REQUIRED_FLAGS        ${CMAKE_REQUIRED_FLAGS}")
        message(STATUS "LAPACKE TEST COMPILE CMAKE_REQUIRED_DEFINITIONS  ${CMAKE_REQUIRED_DEFINITIONS}")
    endif()
    #   Test features
    check_cxx_source_compiles("
        #include <complex>
        #ifdef MKL_AVAILABLE
        #include <mkl_lapacke.h>
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
endfunction()
