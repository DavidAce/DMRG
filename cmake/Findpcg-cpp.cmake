function (find_pcg_cpp)
    find_path(PCG_CPP_INCLUDE_DIR
            include/pcg_random.h
            HINTS ${DMRG_DEPS_INSTALL_DIR}
            PATH_SUFFIXES include pcg-cpp/include
            NO_CMAKE_ENVIRONMENT_PATH
            NO_SYSTEM_ENVIRONMENT_PATH
            NO_CMAKE_SYSTEM_PATH
            )
endfunction()
find_pcg_cpp()
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(pcg-cpp DEFAULT_MSG PCG_CPP_INCLUDE_DIR)


if (pcg_cpp_FOUND)
    add_library(pcg-cpp::pcg-cpp UNKNOWN IMPORTED)
    target_include_directories(pcg-cpp::pcg-cpp SYSTEM INTERFACE ${PCG_CPP_INCLUDE_DIR})
endif()