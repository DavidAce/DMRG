function(find_package_openmp_internal omp_paths omp_names)
    include(${PROJECT_SOURCE_DIR}/cmake-modules/CheckOMPCompiles.cmake)
    # Start by trying to find it using the normal methods
    # find_package(OpenMP) Will at least find the correct header?
    find_package(OpenMP)
    if (OpenMP_FOUND)
        set(OpenMP_LIBRARIES  ${OpenMP_gomp_LIBRARY} ${OpenMP_omp_LIBRARY} ${OpenMP_iomp_LIBRARY})
        set(OpenMP_FLAGS ${OpenMP_CXX_FLAGS})

        # Try to find the header omp.h
        get_filename_component(OpenMP_ROOT "${OpenMP_LIBRARIES}" DIRECTORY)
        find_path(OpenMP_INCLUDE_DIR
              NAMES omp.h
              HINTS ${OpenMP_ROOT}
              PATHS ${omp_paths}
              PATH_SUFFIXES
                /../../include /../include include openmp)

        if(OpenMP_LIBRARIES AND OpenMP_INCLUDE_DIR)
            list(APPEND OpenMP_LIBRARIES  -lpthread -ldl)
            message(STATUS "OpenMP libraries found: ${OpenMP_LIBRARIES}" )
            if(BUILD_SHARED_LIBS)
              check_omp_compiles("${OpenMP_FLAGS}" "${OpenMP_LIBRARIES}" "${OpenMP_INCLUDE_DIR}")
            else()
              check_omp_compiles("${OpenMP_FLAGS}" "-static;${OpenMP_LIBRARIES}" "${OpenMP_INCLUDE_DIR}")
            endif()
        endif()
    endif()


    if(NOT OMP_COMPILES)
        message(STATUS "Checking for OpenMP in given paths: ${omp_paths}")

        unset(OpenMP_FOUND)
        unset(OpenMP_FOUND CACHE)
        unset(OpenMP_FLAGS)
        unset(OpenMP_FLAGS CACHE)
        unset(OpenMP_LIBRARIES)
        unset(OpenMP_LIBRARIES CACHE)
        unset(OpenMP_INCLUDE_DIR)
        unset(OpenMP_INCLUDE_DIR CACHE)

        find_library(OpenMP_LIBRARIES NAMES ${omp_names} PATHS ${omp_paths} ${OpenMP_ROOT} )
        find_path(OpenMP_INCLUDE_DIR NAMES omp.h
            HINTS ${OpenMP_ROOT}
            PATHS
                ${omp_paths}
                /usr/lib/llvm-9/include  /usr/lib/llvm-8/include /usr/lib/llvm-7/include  /usr/include
                PATH_SUFFIXES
                openmp /../../include /../include include )

        if(OpenMP_LIBRARIES AND OpenMP_INCLUDE_DIR)
            list(APPEND OpenMP_LIBRARIES -lpthread -ldl)
            set(OpenMP_DEFS "-D_OPENMP")
            if(${BUILD_SHARED_LIBS})
                check_omp_compiles("${OpenMP_DEFS}" "${OpenMP_LIBRARIES}" "${OpenMP_INCLUDE_DIR}")
            else()
                check_omp_compiles("${OpenMP_DEFS}" "-static;${OpenMP_LIBRARIES}" "${OpenMP_INCLUDE_DIR}")
            endif()
        endif()
    endif()

    if(OMP_COMPILES)
        set(OpenMP_FOUND        TRUE                    PARENT_SCOPE)
        add_library(OpenMP INTERFACE IMPORTED)
        target_link_libraries(OpenMP INTERFACE ${OpenMP_LIBRARIES})
        target_compile_definitions(OpenMP INTERFACE ${OpenMP_DEFS})
        target_include_directories(OpenMP INTERFACE ${OpenMP_INCLUDE_DIR})
        message(STATUS "Found working OpenMP" )
    else()
        message(STATUS "Could not compile simle OpenMP program." )
        set(OMP_COMPILES FALSE PARENT_SCOPE)
    endif()
endfunction()


function(find_package_openmp)
    set(OMP_LIBRARY_NAMES)
    set(OMP_LIBRARY_PATHS)

    if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" AND NOT BUILD_SHARED_LIBS)
        # Append possible locations for iomp5
        list(APPEND OMP_LIBRARY_NAMES  libiomp5${CMAKE_STATIC_LIBRARY_SUFFIX} libiomp${CMAKE_STATIC_LIBRARY_SUFFIX})
        foreach(omp_name ${OMP_LIBRARY_NAMES})
            execute_process(COMMAND ${CMAKE_CXX_COMPILER} -print-file-name=${omp_name}
                    OUTPUT_VARIABLE OMP_LIB
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
            get_filename_component(OMP_PATH ${OMP_LIB} DIRECTORY)
            if(OMP_PATH)
                list(APPEND OMP_LIBRARY_PATHS ${OMP_PATH})
            endif()
        endforeach()
        list(APPEND OMP_LIBRARY_PATHS
                ${MKLROOT}/../intel/lib/intel64
                $ENV{MKLROOT}/../intel/lib/intel64
                $ENV{EBROOTIMKL}/../intel/lib/intel64
                /opt/intel/lib/intel64)
    endif()
    find_package_openmp_internal("${OMP_LIBRARY_PATHS}" "${OMP_LIBRARY_NAMES}" "${BUILD_SHARED_LIBS}")
endfunction()

