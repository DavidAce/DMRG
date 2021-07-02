

if(DMRG_PRINT_INFO)

    # Print host properties
    cmake_host_system_information(RESULT _host_name QUERY HOSTNAME)
    cmake_host_system_information(RESULT _proc_name QUERY PROCESSOR_NAME)
    cmake_host_system_information(RESULT _proc_type QUERY PROCESSOR_DESCRIPTION)
    cmake_host_system_information(RESULT _os_name QUERY OS_NAME)
    cmake_host_system_information(RESULT _os_release QUERY OS_RELEASE)
    cmake_host_system_information(RESULT _os_version QUERY OS_VERSION)
    cmake_host_system_information(RESULT _os_platform QUERY OS_PLATFORM)
    message(STATUS "| DMRG BUILD INFO:\n"
            "-- |----------------\n"
            "-- | ${_host_name}\n"
            "-- | ${_os_name} ${_os_platform} ${_os_release}\n"
            "-- | ${_proc_type} | ${_proc_name}\n"
            "-- | ${_os_version}\n"
            "-- | CMake Version ${CMAKE_VERSION}"
            )

endif ()
