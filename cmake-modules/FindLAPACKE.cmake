

function(CheckLapackeCompiles REQUIRED_FLAGS REQUIRED_DEFINITIONS REQUIRED_LIBRARIES_UNPARSED REQUIRED_INCLUDES)
#    message(STATUS "Checking if lapacke headers are working")
    set(REQUIRED_LIBRARIES)
    foreach(elem ${REQUIRED_LIBRARIES_UNPARSED})
        if(TARGET ${elem})
            get_target_property(lib ${elem} INTERFACE_LINK_LIBRARIES)
            list(APPEND REQUIRED_LIBRARIES ${lib})
        else()
            list(APPEND REQUIRED_LIBRARIES ${elem})
        endif()
    endforeach()


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
#        message(STATUS "Lapacke compiles with the following parameters:")
        set(LAPACKE_COMPILES TRUE PARENT_SCOPE)
    else()
#        message(STATUS "Lapacke does not compile with the following parameters:")
        set(LAPACKE_COMPILES FALSE PARENT_SCOPE)
    endif()
#    message(STATUS "    flags        : " ${REQUIRED_FLAGS})
#    message(STATUS "    definitions  : " ${REQUIRED_DEFINITIONS})
#    message(STATUS "    libraries    : " ${REQUIRED_LIBRARIES})
#    message(STATUS "    includes     : " ${REQUIRED_INCLUDES})

endfunction()




if (NOT TARGET lapacke)
    # Find from MKL
    if(TARGET mkl)
        # Try finding lapacke in MKL library
        message(STATUS "Searching for Lapacke in Intel MKL.")
        get_target_property(MKL_INCLUDE_DIR mkl INTERFACE_INCLUDE_DIRECTORIES)
        get_target_property(MKL_LIBRARIES mkl INTERFACE_LINK_LIBRARIES)

        if(MKL_INCLUDE_DIR)
            CheckLapackeCompiles(" "   "-DMKL_AVAILABLE"  "${MKL_LIBRARIES}" "${MKL_INCLUDE_DIR}")
        endif()

        if(LAPACKE_COMPILES)
            add_library(lapacke INTERFACE)
            add_dependencies(lapacke INTERFACE mkl)
            target_link_libraries(lapacke INTERFACE ${MKL_LIBRARIES})
            target_include_directories(lapacke SYSTEM INTERFACE ${MKL_INCLUDE_DIR})
            message(STATUS "Searching for Lapacke in Intel MKL - Success")
            return()
        else()
            message(STATUS "Searching for Lapacke in Intel MKL - failed")
        endif()
    endif()
endif()




if (NOT TARGET lapacke)
    if (TARGET OpenBLAS)
        message(STATUS "Searching for Lapacke in OpenBLAS")
        get_target_property(LAPACKE_OPENBLAS_LIBRARY   OpenBLAS INTERFACE_LINK_LIBRARIES)
        get_target_property(LAPACKE_OPENBLAS_INCLUDE   OpenBLAS INTERFACE_INCLUDE_DIRECTORIES)
        if(LAPACKE_OPENBLAS_LIBRARY AND LAPACKE_OPENBLAS_INCLUDE)
            CheckLapackeCompiles(" "   " "  "${LAPACKE_OPENBLAS_LIBRARY}" "${LAPACKE_OPENBLAS_INCLUDE}")
        endif()
        if(LAPACKE_COMPILES)
            add_library(lapacke INTERFACE)
            add_dependencies(lapacke INTERFACE OpenBLAS)
            target_link_libraries(lapacke INTERFACE OpenBLAS)
            message(STATUS "Searching for Lapacke in OpenBLAS - Success")
            return()
        else()
            message(STATUS "Searching for Lapacke in OpenBLAS - failed")
        endif()
    endif()
endif()



if (NOT TARGET lapacke)
    message(STATUS "Searching for Lapacke in system")
    find_path(LAPACKE_INCLUDE_DIR
            NAMES lapacke.h
            HINTS ${LAPACKE_DIR} $ENV{CONDA_PREFIX}
            PATHS
            ${CMAKE_INSTALL_PREFIX}
            ${CMAKE_INSTALL_PREFIX}/OpenBLAS
            $ENV{BLAS_DIR}/include
            /usr/include
            /usr/include/x86_64-linux-gnu
            )
    find_library(LAPACKE_LIBRARY
            NAMES liblapacke${CUSTOM_SUFFIX}
            HINTS ${LAPACKE_DIR} $ENV{CONDA_PREFIX}
            PATHS
            $ENV{BLAS_DIR}/lib
            ${CMAKE_INSTALL_PREFIX}
            ${CMAKE_INSTALL_PREFIX}/OpenBLAS
            /usr/lib
            /usr/local/lib
            )
    if(LAPACKE_INCLUDE_DIR OR LAPACKE_LIBRARY)
        CheckLapackeCompiles(" "   " "  "${LAPACKE_LIBRARY} " "${LAPACKE_INCLUDE_DIR} ")
    endif()
    if(LAPACKE_COMPILES)
        add_library(lapacke INTERFACE IMPORTED)
        target_link_libraries(lapacke INTERFACE ${LAPACKE_LIBRARY})
        target_include_directories(lapacke SYSTEM INTERFACE ${LAPACKE_INCLUDE_DIR})
        message(STATUS "Searching for Lapacke in system - Success: ${LAPACKE_LIBRARY}")
        return()
    else()
        message(STATUS "Searching for Lapacke in system - failed")
    endif()
endif()



if(NOT TARGET lapacke)
    message(FATAL_ERROR "Library LAPACKE could not be found.")
endif()