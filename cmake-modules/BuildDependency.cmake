function(build_dependency dep_name install_dir extra_flags)
    set(build_dir    ${CMAKE_BINARY_DIR}/dmrg-deps-build/${dep_name})
    if (DMRG_DEPS_IN_SUBDIR) # h5pp is run with append libsuffix so we don't need to append it again
        set(install_dir ${install_dir}/${dep_name})
        mark_as_advanced(install_dir)
    endif()
    include(cmake-modules/GetNumThreads.cmake)
    get_num_threads(num_threads)

    execute_process( COMMAND  ${CMAKE_COMMAND} -E make_directory ${build_dir})
    execute_process(
            COMMAND  ${CMAKE_COMMAND}
            --parallel ${num_threads}
            -DCMAKE_INSTALL_PREFIX:PATH=${install_dir}
            -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
            -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
            ${extra_flags}
            -G "${CMAKE_GENERATOR}"
            -DCMAKE_GENERATOR_PLATFORM=${CMAKE_GENERATOR_PLATFORM}
            ${PROJECT_SOURCE_DIR}/cmake-modules/external_${dep_name}
            WORKING_DIRECTORY ${build_dir}
            RESULT_VARIABLE config_result
            #ERROR_VARIABLE  config_error
    )
    if(config_result)
        message(STATUS "Got non-zero exit code while configuring ${dep_name}")
        message(STATUS  "build_dir         : ${build_dir}")
        message(STATUS  "install_dir       : ${install_dir}")
        message(STATUS  "extra_flags       : ${extra_flags}")
        message(STATUS  "config_result     : ${config_result}")
        #message(STATUS  "Output saved to ${build_dir}/stdout and ${build_dir}/stderr")
        #file(APPEND ${build_dir}/stdout ${config_result})
        #file(APPEND ${build_dir}/stderr ${config_error})
        #if(DMRG_PRINT_INFO OR CMAKE_VERBOSE_MAKEFILE OR DMRG_PRINT_CHECKS)
        #    message(STATUS "Contents of stdout: \n  ${config_result} \n")
        #    message(STATUS "Contents of stderr: \n  ${config_error}  \n")
        #endif()
        message(FATAL_ERROR "Failed to configure ${dep_name}")
    endif()



    set(ENV{CMAKE_BUILD_PARALLEL_LEVEL} ${num_threads})
    execute_process(COMMAND  ${CMAKE_COMMAND} --build . --parallel ${num_threads}
            WORKING_DIRECTORY "${build_dir}"
            RESULT_VARIABLE build_result
            #ERROR_VARIABLE  build_error
            )

    if(build_result)
        message(STATUS "Got non-zero exit code while building ${dep_name}")
        message(STATUS  "build_dir         : ${build_dir}")
        message(STATUS  "install_dir       : ${install_dir}")
        message(STATUS  "extra_flags       : ${extra_flags}")
        message(STATUS  "build_result      : ${build_result}")
#        message(STATUS  "Output saved to ${build_dir}/stdout and ${build_dir}/stderr")
#        file(APPEND ${build_dir}/stdout ${build_result})
#        file(APPEND ${build_dir}/stderr ${build_error})
#        if(DMRG_PRINT_INFO OR CMAKE_VERBOSE_MAKEFILE OR DMRG_PRINT_CHECKS)
#            message(STATUS "Contents of stdout: \n  ${build_result} \n")
#            message(STATUS "Contents of stderr: \n  ${build_error} \n")
#        endif()
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