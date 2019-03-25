
set(SUFFIX ${CMAKE_STATIC_LIBRARY_SUFFIX})




# Find from MKL
if(MKL_FOUND)
    # Try finding lapacke as module library
    message(STATUS "Attempting to find LAPACKE in loaded environment modules.")
    find_path(LAPACKE_INCLUDE_DIRS
            NAMES mkl_lapacke.h
            PATHS "${MKL_INCLUDE_DIR}"
            NO_DEFAULT_PATH
            )
    if(DEFINED LAPACKE_INCLUDE_DIRS)
        add_library(lapacke INTERFACE)
        set(LAPACKE_FOUND TRUE)
        message(STATUS "Found LAPACKE in Intel MKL")
        target_link_libraries(lapacke INTERFACE mkl::lapacke)
        target_include_directories(lapacke INTERFACE ${LAPACKE_INCLUDE_DIRS})
    endif()
endif()


if(TARGET lapack)
    get_target_property(LAPACKE_LAPACK_LIBRARY   lapack INTERFACE_LINK_LIBRARIES)
    get_target_property(LAPACKE_LAPACK_INCLUDE   lapack INTERFACE_INCLUDE_DIRECTORIES)
endif()



if(NOT LAPACKE_FOUND)

    find_path(LAPACKE_INCLUDE_DIRS
            NAMES lapacke.h
            ${LAPACKE_LAPACK_INCLUDE}
            $ENV{BLAS_DIR}/include
            $ENV{HOME}/.conda/include
            $ENV{HOME}/anaconda3/include
            /usr/include
            /usr/include/x86_64-linux-gnu
            )
    if(LAPACKE_INCLUDE_DIRS)
        message(STATUS "Found lapacke.h headers in: ${LAPACKE_INCLUDE_DIRS}")
    endif()

endif()



add_library(lapacke INTERFACE)
add_dependencies(lapacke INTERFACE lapack)

if(TARGET openblas::lapacke)
    set(LAPACKE_FOUND TRUE)
    target_link_libraries(lapacke INTERFACE openblas::lapacke)
    target_include_directories(lapacke INTERFACE ${LAPACKE_INCLUDE_DIRS} )
endif()






if (NOT LAPACKE_FOUND )
    # Try finding lapacke in OpenBLAS
    message(STATUS "Attempting to find LAPACKE in OpenBLAS")
    if(LAPACKE_LAPACK_LIBRARY AND LAPACKE_INCLUDE_DIRS)
        set(LAPACKE_LIBRARY ${LAPACKE_LAPACK_LIBRARY})
        #   Test features
        include(CheckCXXSourceCompiles)
        set(CMAKE_REQUIRED_LIBRARIES ${LAPACKE_LIBRARY})
        set(CMAKE_REQUIRED_INCLUDES  ${LAPACKE_INCLUDE_DIRS}  ${LAPACKE_LAPACK_INCLUDE})
        check_cxx_source_compiles("
        #include <lapacke.h>
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
        " LAPACKE_COMPILES_OPENBLAS)

        if(NOT LAPACKE_COMPILES_OPENBLAS)
            unset(LAPACKE_FOUND)
        else()
            set(LAPACKE_FOUND TRUE)
        endif()
    else()
        message(STATUS "Could not find LAPACKE in OpenBLAS")
    endif()

endif()


# Find installed in system
if(NOT LAPACKE_FOUND)
    message(STATUS "Attempting to find LAPACKE in system")
    find_library(LAPACKE_LIBRARY_SYSTEM
            NAMES liblapacke${SUFFIX}
            PATHS
            $ENV{BLAS_DIR}/lib
            $ENV{BLAS_DIR}/lib
            $ENV{HOME}/.conda/lib
            $ENV{HOME}/anaconda3/lib
            /usr/lib/x86_64-linux-gnu
            )

    if(LAPACKE_LIBRARY_SYSTEM AND LAPACKE_INCLUDE_DIRS)
        set(LAPACKE_LIBRARY         ${LAPACKE_LIBRARY_SYSTEM} ${LAPACKE_LAPACK_LIBRARY})

        include(CheckCXXSourceCompiles)
        set(CMAKE_REQUIRED_LIBRARIES )
        set(CMAKE_REQUIRED_LIBRARIES ${LAPACKE_LIBRARY})
        set(CMAKE_REQUIRED_INCLUDES  ${LAPACKE_INCLUDE_DIRS} ${LAPACKE_LAPACK_INCLUDE})

        check_cxx_source_compiles("
        #include <lapacke.h>
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
           LAPACKE_dgels(LAPACK_ROW_MAJOR,'N',m,n,nrhs,*a,lda,*b,ldb);
           return 0;
        }
        " LAPACKE_COMPILES_SYSTEM)

        if(NOT LAPACKE_COMPILES_SYSTEM)
            unset(LAPACKE_FOUND)
        else()
            set(LAPACKE_FOUND TRUE)
        endif()
    else()
        message(STATUS "Could not find LAPACKE in system")
    endif()
endif()


if(LAPACKE_FOUND)
    message(STATUS "LAPACKE LIBRARY : ${LAPACKE_LIBRARY}")
    message(STATUS "LAPACKE INCLUDE : ${LAPACKE_INCLUDE_DIRS}")

    target_link_libraries(lapacke INTERFACE ${LAPACKE_LIBRARY})
    target_include_directories(lapacke INTERFACE ${LAPACKE_INCLUDE_DIRS})
else()
    message(WARNING "Could not find working package LAPACKE")
endif()

