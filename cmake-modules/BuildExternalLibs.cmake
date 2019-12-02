function(build_external_libs library_name config_dir build_dir install_dir extra_flags)
    message(STATUS "library_name      : ${library_name}")
    message(STATUS "config_dir        : ${config_dir}")
    message(STATUS "build_dir         : ${build_dir}")
    message(STATUS "install_dir       : ${install_dir}")
    message(STATUS "extra_flags       : ${extra_flags}")
    message(STATUS "CMAKE_SOURCE_DIR  : ${CMAKE_SOURCE_DIR}")
    execute_process( COMMAND  ${CMAKE_COMMAND} -E make_directory ${config_dir}/${library_name})
    execute_process(
            COMMAND  ${CMAKE_COMMAND}
            -DEXTERNAL_BUILD_DIR:PATH=${build_dir}
            -DEXTERNAL_INSTALL_DIR:PATH=${install_dir}
            -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
            -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
            ${extra_flags}
            -G "CodeBlocks - Unix Makefiles"
            ${CMAKE_SOURCE_DIR}/cmake-modules/external_${library_name}
            WORKING_DIRECTORY ${config_dir}/${library_name}
            RESULT_VARIABLE config_result
    )

    execute_process(COMMAND  ${CMAKE_COMMAND}  --build . --target  all --parallel
            WORKING_DIRECTORY "${config_dir}/${library_name}"
            RESULT_VARIABLE build_result
    )

    set(config_result ${config_result} PARENT_SCOPE)
    set(build_result ${build_result} PARENT_SCOPE)
endfunction()