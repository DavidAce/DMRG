function(find_package_openmp)
    include(${PROJECT_SOURCE_DIR}/cmake-modules/CheckOMPCompiles.cmake)
    if(NOT TARGET OpenMP AND BUILD_SHARED_LIBS)
        find_package(OpenMP)
        if(TARGET OpenMP::OpenMP_CXX)
            include(cmake-modules/PrintTargetInfo.cmake)
            print_target_info(OpenMP::OpenMP_CXX)
            check_omp_compiles("" "OpenMP::OpenMP_CXX" "")
            if(OMP_COMPILES)
                set(OpenMP_FOUND TRUE PARENT_SCOPE)
                add_library(OpenMP ALIAS  OpenMP::OpenMP_CXX)
                message(STATUS "Found working OpenMP" )
                return()
            endif()
        endif()
    endif()
    if(NOT TARGET OpenMP)
        find_library(omp_lib_iomp5 NAMES
                libiomp5${CMAKE_STATIC_LIBRARY_SUFFIX}
                PATHS
                /usr/lib/llvm-9/include
                /usr/lib/llvm-8/include
                /usr/lib/llvm-7/include
                /usr/include /usr/lib /usr/local
                ${MKLROOT}/../intel/lib/intel64
                $ENV{MKLROOT}/../intel/lib/intel64
                $ENV{EBROOTIMKL}/../intel/lib/intel64
                /opt/intel/lib/intel64)

        set(omp_lib_candidates -lgomp -lomp liomp5 ${omp_lib_iomp5})
        foreach(lib ${omp_lib_candidates})
            check_omp_compiles("-D_OPENMP" "-static;${lib};-lpthread;-lrt;-ldl" "")
            if(OMP_COMPILES)
                set(OpenMP_FOUND TRUE PARENT_SCOPE)
                add_library(OpenMP INTERFACE IMPORTED)
                target_link_libraries(OpenMP INTERFACE ${lib} -lpthread -lrt -ldl)
                target_compile_definitions(OpenMP INTERFACE -D_OPENMP)
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


