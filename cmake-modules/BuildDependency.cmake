function(build_dependency dep_name extra_flags)
    set(build_dir    ${CMAKE_BINARY_DIR}/external-deps/${dep_name})
    set(install_dir  ${CMAKE_INSTALL_PREFIX})

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

    execute_process(COMMAND  ${CMAKE_COMMAND} --build . --target all
            WORKING_DIRECTORY "${build_dir}"
            RESULT_VARIABLE build_result
    )

    set(config_result ${config_result} PARENT_SCOPE)
    set(build_result ${build_result} PARENT_SCOPE)
endfunction()