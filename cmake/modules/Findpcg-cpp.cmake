function (find_pcg_cpp)
    find_path(PCG_CPP_INCLUDE_DIR
            pcg_random.hpp
            HINTS ${DMRG_DEPS_INSTALL_DIR}
            PATH_SUFFIXES pcg-cpp include/pcg-cpp pcg-cpp/include
            NO_CMAKE_ENVIRONMENT_PATH
            NO_SYSTEM_ENVIRONMENT_PATH
            NO_CMAKE_SYSTEM_PATH
            )
endfunction()
find_pcg_cpp()
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(pcg-cpp DEFAULT_MSG PCG_CPP_INCLUDE_DIR)


if (pcg-cpp_FOUND AND NOT TARGET pcg-cpp::pcg-cpp)
    message(DEBUG "Defining target pcg-cpp::pcg-cpp")
    add_library(pcg-cpp::pcg-cpp INTERFACE IMPORTED)
    target_include_directories(pcg-cpp::pcg-cpp SYSTEM INTERFACE ${PCG_CPP_INCLUDE_DIR})
endif()