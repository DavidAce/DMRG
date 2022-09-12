cmake_minimum_required(VERSION 3.15)

function(pad_string OUT_VARIABLE DESIRED_LENGTH FILL_CHAR VALUE)
    string(LENGTH "${VALUE}" VALUE_LENGTH)
    math(EXPR REQUIRED_PADS "${DESIRED_LENGTH} - ${VALUE_LENGTH}")
    set(PAD ${VALUE})
    if(REQUIRED_PADS GREATER 0)
        math(EXPR REQUIRED_MINUS_ONE "${REQUIRED_PADS} - 1")
        foreach(FOO RANGE ${REQUIRED_MINUS_ONE})
            set(PAD "${PAD}${FILL_CHAR}")
        endforeach()
    endif()
    set(${OUT_VARIABLE} "${PAD}" PARENT_SCOPE)
endfunction()

function(print_padded_option opt)
    pad_string(popt "24" " " "${opt}" )
    get_property(adv CACHE ${opt} PROPERTY ADVANCED SET)
    if(adv)
        message(TRACE   "| ${popt}: ${${opt}}")
    else()
        message(DEBUG "| ${popt}: ${${opt}}")
    endif()
endfunction()

# Print host properties
cmake_host_system_information(RESULT _host_name QUERY HOSTNAME)
cmake_host_system_information(RESULT _proc_type QUERY PROCESSOR_DESCRIPTION)
cmake_host_system_information(RESULT _os_name QUERY OS_NAME)
cmake_host_system_information(RESULT _os_release QUERY OS_RELEASE)
cmake_host_system_information(RESULT _os_version QUERY OS_VERSION)
cmake_host_system_information(RESULT _os_platform QUERY OS_PLATFORM)

message(VERBOSE "| HOST INFO")
message(VERBOSE "|--------------------------")
message(VERBOSE "| CMake Version ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}")
message(VERBOSE "| ${_host_name}")
message(VERBOSE "| ${_os_name} ${_os_platform} ${_os_release}")
message(VERBOSE "| ${_proc_type}")
message(VERBOSE "| ${_os_version}")
message(VERBOSE "|--------------------------")

# Print all CMake options
message(VERBOSE "| OPTIONS")
message(VERBOSE "|--------------------------")
print_padded_option(CMAKE_BUILD_TYPE)
print_padded_option(CMAKE_PREFIX_PATH)
print_padded_option(CMAKE_INSTALL_PREFIX)
print_padded_option(BUILD_SHARED_LIBS)
print_padded_option(CMAKE_VERBOSE_MAKEFILE)
get_cmake_property(proj_opts VARIABLES)
mark_as_advanced(proj_opts)
foreach (opt ${proj_opts})
    if(opt MATCHES "DMRG_|COMPILER_" AND NOT opt MATCHES "VALID|CANDIDATE|_COMPILER")
        print_padded_option(${opt})
    endif()
endforeach()
message(VERBOSE "|--------------------------")

