



function(check_omp_compiles REQUIRED_FLAGS REQUIRED_LIBRARIES REQUIRED_INCLUDES)
    include(CheckIncludeFileCXX)
    set(CMAKE_REQUIRED_FLAGS     "${REQUIRED_FLAGS}" )
    set(CMAKE_REQUIRED_LIBRARIES "${REQUIRED_LIBRARIES}")
    set(CMAKE_REQUIRED_INCLUDES  "${REQUIRED_INCLUDES}")
    message(STATUS "REQFLAG: ${REQUIRED_FLAGS}")
    message(STATUS "REQLIBS: ${REQUIRED_LIBRARIES}")
    message(STATUS "REQINCS: ${REQUIRED_INCLUDES}")

    unset(has_omp_h CACHE)
    unset(OMP_COMPILES CACHE)
    unset(OMP_COMPILES)
    check_include_file_cxx(omp.h    has_omp_h)
    if(NOT has_omp_h)
        message(WARNING "Header missing: omp.h")
    endif()

    include(CheckCXXSourceCompiles)
    check_cxx_source_compiles("
            #include <omp.h>
            #include <iostream>
            int main() {
                omp_set_num_threads(4);
                std::cout << \"OMP Threads \" << omp_get_max_threads() << std::endl;
                return 0;
            }
            " OMP_COMPILES
            )
    if(NOT OMP_COMPILES)
        message(STATUS "Setting USE_OpenMP OFF -- Unable to compile a simple OpenMP program")
        set(USE_OpenMP OFF PARENT_SCOPE)
    endif()
endfunction()



if (NOT USE_OpenMP)
    return()
endif()

find_package(OpenMP)
if (OpenMP_FOUND)
    message(STATUS "OpenMP ON" )
#    get_cmake_property(_variableNames VARIABLES)
#    foreach (_variableName ${_variableNames})
#        string(TOLOWER  "${_variableName}" varnamelower)
#        if( "${varnamelower}" MATCHES "openmp")
#            message(STATUS "${_variableName}=${${_variableName}}")
#        endif()
#    endforeach()
#
#    include(cmake-modules/PrintTargetProperties.cmake)
#    print_target_properties(OpenMP::OpenMP_CXX)

    set(OpenMP_LIBRARIES  ${OpenMP_gomp_LIBRARY} ${OpenMP_omp_LIBRARY} ${OpenMP_iomp_LIBRARY})
    set(OpenMP_FLAGS ${OpenMP_CXX_FLAGS})

    if(NOT BUILD_SHARED_LIBS)
        list(INSERT OpenMP_LIBRARIES 0 -static)
    endif()
    check_omp_compiles("${OpenMP_FLAGS}" "${OpenMP_LIBRARIES};Threads::Threads" "${OpenMP_INCLUDE_DIR}")
else()
    message(WARNING "Setting USE_OpenMP OFF -- Could not be found." )
    set(USE_OpenMP OFF)
endif()




#if("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
#    set(TEST_REQUIRED_FLAGS     "" )
#    set(TEST_REQUIRED_LIBRARIES "")
#    set(TEST_REQUIRED_INCLUDES "${OpenMP_INCLUDE_DIR}")
#elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
#    set(TEST_REQUIRED_FLAGS     -static -v -stdlib=libstdc++ -std=c++17 --gcc-toolchain=/usr/lib/gcc/x86_64-linux-gnu/8 -lpthread -ldl)
#    #        set(TEST_REQUIRED_FLAGS     "-fopenmp" )
#    set(TEST_REQUIRED_LIBRARIES -static -v ${OpenMP_LIBRARIES} -lpthread -ldl)
#    set(TEST_REQUIRED_INCLUDES  ${OpenMP_INCLUDE_DIR} /usr/include/c++/8)
#endif()

#
#add_library(OpenMP INTERFACE)
#if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND NOT ${BUILD_SHARED_LIBS})
#    message(STATUS "Clang cannot link \"libomp\" statically. Attempting to link with libiomp")
#    find_library(OpenMP_LIBRARIES
#            NAMES libiomp5.a  libomp.a libgomp.a iomp5 omp gomp
#            PATH_SUFFIXES
#            lib lib64
#            PATHS
#            ${GCC_TOOLCHAIN}
#            $ENV{MKLROOT}/../compilers_and_libraries/linux/lib/intel64/
#            )
#    find_path(
#            OpenMP_INCLUDE_DIR
#            NAMES omp.h
#            PATH_SUFFIXES include
#            PATHS
#            /usr
#            /usr/lib
#            /usr/local
#            /usr/lib/gcc/x86_64-linux-gnu/8
#            ${GCC_TOOLCHAIN} $ENV{MKLROOT}/../compilers_and_libraries/linux/ $ENV{EBROOTGLOG} $ENV{MKL_ROOT} $ENV{MKLROOT}
#    )
#    #        find_package(OpenMP)
#    set(OpenMP_FLAGS  -static -stdlib=libstdc++ -std=c++17)
#    list(APPEND OpenMP_LIBRARIES  -stdlib=libstdc++ -std=c++17 -lpthread -ldl)
#    target_compile_options(OpenMP INTERFACE ${OpenMP_FLAGS})
#    target_link_libraries(OpenMP INTERFACE ${OpenMP_LIBRARIES})
#    target_include_directories(OpenMP INTERFACE ${OpenMP_INCLUDE_DIR})
#
#    if(OpenMP_LIBRARIES)
#        message(STATUS "OpenMP ON" )
#        set(OpenMP_FOUND ON)
#        target_compile_options    (${PROJECT_NAME} PRIVATE -fopenmp=libiomp5)
#        target_link_libraries     (${PROJECT_NAME} PRIVATE ${OpenMP_LIBRARIES})
#    else()
#        message(WARNING "
#            Could not find libiomp.a
#            Please use GNU C++ instead to link \"libgomp.a\" statically, or turn off OpenMP
#            by passing \"-DUSE_OpenMP:BOOL=OFF\" to CMake to avoid this warning.
#            Turning off OpenMP.")
#        message(STATUS "OpenMP OFF")
#        set(USE_OpenMP OFF)
#        unset(OpenMP_LIBRARIES)
#    endif()
#else()
#    find_package(OpenMP)
#    if (OpenMP_FOUND)
#        message(STATUS "OpenMP ON" )
#        set(OpenMP_LIBRARIES  ${OpenMP_gomp_LIBRARY} ${OpenMP_omp_LIBRARY} ${OpenMP_iomp_LIBRARY})
#        target_compile_options    (OpenMP PRIVATE ${OpenMP_CXX_FLAGS})
#        target_link_libraries     (OpenMP PRIVATE ${OpenMP_LIBRARIES})
#    else()
#        message(WARNING "OpenMP OFF -- Could not be found." )
#    endif()
#endif()


