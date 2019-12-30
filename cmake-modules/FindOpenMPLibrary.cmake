function(find_package_openmp)
    include(${PROJECT_SOURCE_DIR}/cmake-modules/CheckOMPCompiles.cmake)
    if(NOT TARGET openmp::openmp AND BUILD_SHARED_LIBS)
        find_package(OpenMP)
        if(TARGET OpenMP::OpenMP_CXX)
            include(cmake-modules/PrintTargetProperties.cmake)
            if(OMP_COMPILES)
                set(OpenMP_FOUND TRUE PARENT_SCOPE)
                add_library(openmp::openmp INTERFACE IMPORTED)
                target_link_libraries(openmp::openmp INTERFACE OpenMP::OpenMP_CXX)
                message(STATUS "Found working OpenMP" )
                return()
            endif()
        endif()
    endif()
    if(NOT TARGET openmp::openmp)
        find_library(omp_lib_iomp5 NAMES
                iomp5
                PATHS
                /usr/lib/llvm-9/include
                /usr/lib/llvm-8/include
                /usr/lib/llvm-7/include
                /usr/include /usr/lib /usr/local
                ${MKLROOT}/../intel/lib/intel64
                $ENV{MKLROOT}/../intel/lib/intel64
                $ENV{EBROOTIMKL}/../intel/lib/intel64
                /opt/intel/lib/intel64)

        if(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
            list(APPEND ${omp_lib_iomp5})
        endif()
        list(APPEND omp_lib_candidates gomp omp iomp5 ${omp_lib_iomp5})
        foreach(lib ${omp_lib_candidates})
            check_omp_compiles("-D_OPENMP" "-static;${lib};pthread;rt;dl" "")
            if(OMP_COMPILES)
                set(OpenMP_FOUND TRUE PARENT_SCOPE)
                add_library(openmp::openmp INTERFACE IMPORTED)
                target_link_libraries(openmp::openmp INTERFACE ${lib} pthread rt dl)
                target_compile_definitions(openmp::openmp INTERFACE -D_OPENMP)
                message(STATUS "OpenMP library found: ${lib}" )
                break()
            endif()
        endforeach()
    endif()
    if(NOT OMP_COMPILES)
        message(STATUS "Could not compile simple OpenMP program." )
        set(OpenMP_FOUND FALSE PARENT_SCOPE)
    endif()
endfunction()


