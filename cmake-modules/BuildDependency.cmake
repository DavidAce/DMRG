function(build_dependency dep_name install_dir extra_flags)
    set(build_dir    ${CMAKE_BINARY_DIR}/external-deps/${dep_name})
    execute_process( COMMAND  ${CMAKE_COMMAND} -E make_directory ${build_dir})
    execute_process(
            COMMAND  ${CMAKE_COMMAND}
            -DCMAKE_INSTALL_PREFIX:PATH=${install_dir}
            -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
            ${extra_flags}
            -G "CodeBlocks - Unix Makefiles"
            ${PROJECT_SOURCE_DIR}/cmake-modules/external_${dep_name}
            WORKING_DIRECTORY ${build_dir}
            RESULT_VARIABLE config_result
    )
    if(${config_result})
        message(STATUS "Got non-zero exit code ${config_result} while configuring ${dep_name}")
        message(STATUS  "build_dir         : ${build_dir}")
        message(STATUS  "install_dir       : ${install_dir}")
        message(STATUS  "extra_flags       : ${extra_flags}")
        message(STATUS  "config_result     : ${config_result}")
        message(STATUS  "Output saved to ${build_dir}/stdout and ${build_dir}/stderr")
    endif()


    execute_process(COMMAND  ${CMAKE_COMMAND} --build . --target all  --parallel
            WORKING_DIRECTORY "${build_dir}"
            RESULT_VARIABLE build_result
            )

    if(${build_result})
        message(STATUS "Got non-zero exit code ${build_result} while building ${dep_name}")
        message(STATUS  "build_dir         : ${build_dir}")
        message(STATUS  "install_dir       : ${install_dir}")
        message(STATUS  "extra_flags       : ${extra_flags}")
        message(STATUS  "build_result      : ${build_result}")
        message(STATUS  "Output saved to ${build_dir}/stdout and ${build_dir}/stderr")
    endif()

endfunction()