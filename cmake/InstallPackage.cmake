cmake_minimum_required(VERSION 3.15)
include(cmake/CheckCompile.cmake)

# Dumps cached variables to PKG_INIT_CACHE_FILE so that we can propagate
# the current build configuration to dependencies
function(generate_init_cache)
    set(PKG_INIT_CACHE_FILE ${CMAKE_BINARY_DIR}/CMakeFiles/CMakeTmp/init-cache.cmake)
    set(PKG_INIT_CACHE_FILE ${PKG_INIT_CACHE_FILE} PARENT_SCOPE)
    file(WRITE  ${PKG_INIT_CACHE_FILE} "# These variables will initialize the CMake cache for subprocesses.\n")
    get_cmake_property(vars CACHE_VARIABLES)
    foreach(var ${vars})
        if(var MATCHES "CMAKE_CACHE|CMAKE_HOME|CMAKE_EXTRA|CMAKE_PROJECT|MACRO")
            continue()
        endif()
        get_property(help CACHE "${var}" PROPERTY HELPSTRING)
        get_property(type CACHE "${var}" PROPERTY TYPE)
        string(REPLACE "\\" "/" ${var} "${${var}}") # Fix windows backslash paths
        string(REPLACE "\"" "\\\"" help "${help}") #Fix quotes on some cuda-related descriptions
        string(REPLACE "\"" "\\\"" ${var} "${${var}}") #Fix quotes
        file(APPEND ${PKG_INIT_CACHE_FILE} "set(${var} \"${${var}}\" CACHE ${type} \"${help}\" FORCE)\n")
    endforeach()
endfunction()


# This function will configure, build and install a package at configure-time
# by running cmake in a subprocess. The current CMake configuration is transmitted
# by setting the flags manually.
function(install_package pkg_name)
    file(LOCK $ENV{HOME} DIRECTORY GUARD FUNCTION TIMEOUT 600)

    set(options CONFIG MODULE CHECK DEBUG)
    set(oneValueArgs VERSION INSTALL_DIR BUILD_DIR TARGET_NAME PREFIX_PKGNAME)
    set(multiValueArgs HINTS PATHS PATH_SUFFIXES COMPONENTS DEPENDS CMAKE_ARGS)
    cmake_parse_arguments(PARSE_ARGV 1 PKG "${options}" "${oneValueArgs}" "${multiValueArgs}")

    # Set defaults
    if(NOT PKG_BUILD_DIR)
        set(PKG_BUILD_DIR ${CMAKE_BINARY_DIR}/pkg-build)
    endif()
    if(NOT PKG_INSTALL_DIR)
        set(PKG_INSTALL_DIR ${CMAKE_INSTALL_PREFIX})
    endif()

    # Further parsing
    set(build_dir ${PKG_BUILD_DIR}/${pkg_name})
    set(install_dir ${PKG_INSTALL_DIR})
    set(target_name ${pkg_name}::${pkg_name})

    if(PKG_PREFIX_PKGNAME)
        set(install_dir ${install_dir}/${pkg_name})
    endif()
    if(PKG_TARGET_NAME)
        set(target_name ${PKG_TARGET_NAME})
    endif()
    if(PKG_CONFIG)
        set(CONFIG CONFIG)
    endif()
    if(PKG_COMPONENTS)
        set(COMPONENTS COMPONENTS)
    endif()

    foreach (tgt ${PKG_DEPENDS})
        if(NOT TARGET ${tgt})
            list(APPEND PKG_MISSING_TARGET ${tgt})
        endif()
    endforeach()
    if(PKG_MISSING_TARGET)
        message(FATAL_ERROR "Could not install ${pkg_name}: dependencies missing [${PKG_MISSING_TARGET}]")
    endif()

    # We set variables here that allows us to find packages with CMAKE_PREFIX_PATH
    list(APPEND CMAKE_PREFIX_PATH ${install_dir} $ENV{CMAKE_PREFIX_PATH} ${PKG_INSTALL_DIR} ${install_dir} ${CMAKE_INSTALL_PREFIX})
    list(REMOVE_DUPLICATES CMAKE_PREFIX_PATH)
    set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH}" CACHE STRING "" FORCE)
    set(CMAKE_FIND_PACKAGE_PREFER_CONFIG TRUE)

    if(PKG_MODULE)
        find_package(${pkg_name} ${PKG_VERSION} ${COMPONENTS} ${PKG_COMPONENTS})
    else()
        find_package(${pkg_name} ${PKG_VERSION}
                HINTS ${PKG_HINTS}
                PATHS ${PKG_PATHS}
                PATH_SUFFIXES ${PKG_PATH_SUFFIXES}
                ${COMPONENTS} ${PKG_COMPONENTS}
                ${CONFIG}
                # These lets us ignore system packages when pkg manager matches "cmake"
                NO_SYSTEM_ENVIRONMENT_PATH #5
                NO_CMAKE_PACKAGE_REGISTRY #6
                NO_CMAKE_SYSTEM_PATH #7
                NO_CMAKE_SYSTEM_PACKAGE_REGISTRY #8
                )
    endif()


    if(${pkg_name}_FOUND)
        if(PKG_DEPENDS)
            target_link_libraries(${target_name} INTERFACE ${PKG_DEPENDS})
        endif()
        if(PKG_CHECK)
            check_compile(${pkg_name} ${target_name} ${PROJECT_SOURCE_DIR}/cmake/compile/${pkg_name}.cpp)
            if(PKG_DEBUG AND NOT check_compile_${pkg_name} AND EXISTS "${CMAKE_BINARY_DIR}/CMakeFiles/CMakeError.log")
                file(READ "${CMAKE_BINARY_DIR}/CMakeFiles/CMakeError.log" ERROR_LOG)
                message(STATUS "CMakeError.log: \n ${ERROR_LOG}")
            endif()
        endif()
        message(STATUS "Found ${pkg_name}")
        return()
    endif()

    message(STATUS "${pkg_name} will be installed into ${install_dir}")

    set(CMAKE_CXX_STANDARD 17 CACHE STRING "")
    set(CMAKE_CXX_STANDARD_REQUIRED TRUE CACHE BOOL "")
    set(CMAKE_CXX_EXTENSIONS FALSE CACHE BOOL "")

    # Set policies for CMakeLists in packages that require older CMake versions
    set(CMAKE_POLICY_DEFAULT_CMP0074 NEW CACHE STRING "Honor <PackageName>_ROOT")
    set(CMAKE_POLICY_DEFAULT_CMP0091 NEW CACHE STRING "Use MSVC_RUNTIME_LIBRARY") # Fixes spdlog on MSVC

    if(CMAKE_GENERATOR MATCHES "Ninja")
        set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE CACHE BOOL "Otherwise fails with -G Ninja" FORCE)
    endif()

# Generate an init cache to propagate the current configuration
    generate_init_cache()

    # Configure the package
    execute_process( COMMAND  ${CMAKE_COMMAND} -E make_directory ${build_dir})
    execute_process( COMMAND  ${CMAKE_COMMAND} -E remove ${build_dir}/CMakeCache.txt)
    execute_process(
            COMMAND
            ${CMAKE_COMMAND}
            -C ${PKG_INIT_CACHE_FILE}                # For the subproject in external_<libname>
            -DINIT_CACHE_FILE=${PKG_INIT_CACHE_FILE} # For externalproject_add inside the subproject
            -DCMAKE_INSTALL_PREFIX:PATH=${install_dir}
            ${PKG_CMAKE_ARGS}
            ${PROJECT_SOURCE_DIR}/cmake/external_${pkg_name}
            WORKING_DIRECTORY ${build_dir}
            RESULT_VARIABLE config_result
    )
    if(config_result)
        message(STATUS "Got non-zero exit code while configuring ${pkg_name}")
        message(STATUS  "build_dir         : ${build_dir}")
        message(STATUS  "install_dir       : ${install_dir}")
        message(STATUS  "extra_flags       : ${extra_flags}")
        message(STATUS  "config_result     : ${config_result}")
        message(FATAL_ERROR "Failed to configure ${pkg_name}")
    endif()


    # Make sure to do multithreaded builds if possible
    include(cmake/GetNumThreads.cmake)
    get_num_threads(num_threads)

    # Build the package
    if(CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
        # This branch is for multi-config generators such as Visual Studio 16 2019
        foreach(config ${CMAKE_CONFIGURATION_TYPES})
            execute_process(COMMAND  ${CMAKE_COMMAND} --build . --parallel ${num_threads} --config ${config}
                    WORKING_DIRECTORY "${build_dir}"
                    RESULT_VARIABLE build_result
                    )
            if(build_result)
                message(STATUS "Got non-zero exit code while building package: ${pkg_name}")
                message(STATUS  "build_type        : ${config}")
                message(STATUS  "build_dir         : ${build_dir}")
                message(STATUS  "install_dir       : ${install_dir}")
                message(STATUS  "cmake_args        : ${PKG_CMAKE_ARGS}")
                message(STATUS  "build_result      : ${build_result}")
                message(FATAL_ERROR "Failed to build package: ${pkg_name}")
            endif()
        endforeach()
    else()
        # This is for single-config generators such as Unix Makefiles and Ninja
        execute_process(COMMAND  ${CMAKE_COMMAND} --build . --parallel ${num_threads}
                WORKING_DIRECTORY "${build_dir}"
                RESULT_VARIABLE build_result
                )
        if(build_result)
            message(STATUS "Got non-zero exit code while building package: ${pkg_name}")
            message(STATUS  "build_dir         : ${build_dir}")
            message(STATUS  "install_dir       : ${install_dir}")
            message(STATUS  "cmake_args        : ${PKG_CMAKE_ARGS}")
            message(STATUS  "build_result      : ${build_result}")
            message(FATAL_ERROR "Failed to build package: ${pkg_name}")
        endif()
    endif()


    # Copy the install manifest if it exists
    file(GLOB_RECURSE INSTALL_MANIFEST "${build_dir}/*/install_manifest*.txt")
    foreach(manifest ${INSTALL_MANIFEST})
        get_filename_component(manifest_filename ${manifest} NAME_WE)
        message(STATUS "Copying install manifest: ${manifest}")
        configure_file(${manifest} ${CMAKE_CURRENT_BINARY_DIR}/${manifest_filename}_${dep_name}.txt)
    endforeach()


    # Find the package again
    if(PKG_MODULE)
        find_package(${pkg_name} ${PKG_VERSION} ${COMPONENTS} ${PKG_COMPONENTS} REQUIRED)
    else()
        find_package(${pkg_name} ${PKG_VERSION}
                HINTS ${install_dir}
                PATH_SUFFIXES ${PKG_PATH_SUFFIXES}
                ${COMPONENTS} ${PKG_COMPONENTS}
                NO_DEFAULT_PATH REQUIRED)
    endif()

    if(PKG_DEPENDS)
        target_link_libraries(${target_name} INTERFACE ${PKG_DEPENDS})
    endif()
    if(PKG_CHECK)
        check_compile(${pkg_name} ${target_name} ${PROJECT_SOURCE_DIR}/cmake/compile/${pkg_name}.cpp)
        if(PKG_DEBUG AND NOT check_compile_${pkg_name} AND EXISTS "${CMAKE_BINARY_DIR}/CMakeFiles/CMakeError.log")
            file(READ "${CMAKE_BINARY_DIR}/CMakeFiles/CMakeError.log" ERROR_LOG)
            message(STATUS "CMakeError.log: \n ${ERROR_LOG}")
        endif()
    endif()
endfunction()

