function(build_dependency dep_name install_dir extra_flags)
    set(build_dir    ${CMAKE_BINARY_DIR}/dmrg-deps-build/${dep_name})
    if (DMRG_DEPS_IN_SUBDIR) # h5pp is run with append libsuffix so we don't need to append it again
        set(install_dir ${install_dir}/${dep_name})
        mark_as_advanced(install_dir)
    endif()
    include(cmake-modules/GetNumThreads.cmake)
    get_num_threads(num_threads)
    execute_process( COMMAND  ${CMAKE_COMMAND} -E remove ${build_dir}/CMakeCache.txt)
    execute_process( COMMAND  ${CMAKE_COMMAND} -E make_directory ${build_dir})
    execute_process(
            COMMAND  ${CMAKE_COMMAND}
            --parallel ${num_threads}
            -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
            -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
            -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
            -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
            -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
            -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=${CMAKE_POSITION_INDEPENDENT_CODE}
            -DCMAKE_INSTALL_PREFIX:PATH=${install_dir}
            -DCMAKE_CXX_STANDARD:STRING=${CMAKE_CXX_STANDARD}
            -DCMAKE_CXX_STANDARD_REQUIRED:BOOL=${CMAKE_CXX_STANDARD_REQUIRED}
            -DCMAKE_CXX_EXTENSIONS:BOOL=${CMAKE_CXX_EXTENSIONS}
            -DCMAKE_CXX_FLAGS_INIT:STRING=${CMAKE_CXX_FLAGS}
            -DCMAKE_CXX_FLAGS_RELEASE_INIT:STRING=${CMAKE_CXX_FLAGS_RELEASE}
            -DCMAKE_CXX_FLAGS_DEBUG_INIT:STRING=${CMAKE_CXX_FLAGS_DEBUG}
            -DCMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT:STRING=${CMAKE_CXX_FLAGS_RELWITHDEBINFO}
            -DCMAKE_CXX_FLAGS_MINSIZEREL_INIT:STRING=${CMAKE_CXX_FLAGS_MINSIZEREL}
            ${extra_flags}
            -G "${CMAKE_GENERATOR}"
            -DCMAKE_GENERATOR_PLATFORM=${CMAKE_GENERATOR_PLATFORM}
            ${PROJECT_SOURCE_DIR}/cmake-modules/external_${dep_name}
            WORKING_DIRECTORY ${build_dir}
            RESULT_VARIABLE config_result
    )
    if(config_result)
        message(STATUS "Got non-zero exit code while configuring ${dep_name}")
        message(STATUS  "build_dir         : ${build_dir}")
        message(STATUS  "install_dir       : ${install_dir}")
        message(STATUS  "extra_flags       : ${extra_flags}")
        message(STATUS  "config_result     : ${config_result}")
        message(FATAL_ERROR "Failed to configure ${dep_name}")
    endif()



    set(ENV{CMAKE_BUILD_PARALLEL_LEVEL} ${num_threads})
    message(STATUS "Starting build with ${num_threads} threads")
    execute_process(COMMAND  ${CMAKE_COMMAND} --build . --parallel ${num_threads}
            WORKING_DIRECTORY "${build_dir}"
            RESULT_VARIABLE build_result
            )

    if(build_result)
        message(STATUS "Got non-zero exit code while building ${dep_name}")
        message(STATUS  "build_dir         : ${build_dir}")
        message(STATUS  "install_dir       : ${install_dir}")
        message(STATUS  "extra_flags       : ${extra_flags}")
        message(STATUS  "build_result      : ${build_result}")
        message(FATAL_ERROR "Failed to build ${dep_name}")
    endif()

    # Copy the install manifest if it exists
    file(GLOB_RECURSE INSTALL_MANIFEST "${build_dir}/*/install_manifest*.txt")
    foreach(manifest ${INSTALL_MANIFEST})
        get_filename_component(manifest_filename ${manifest} NAME_WE)
        message(STATUS "Copying install manifest: ${manifest}")
        configure_file(${manifest} ${CMAKE_CURRENT_BINARY_DIR}/${manifest_filename}_${dep_name}.txt)
    endforeach()

endfunction()