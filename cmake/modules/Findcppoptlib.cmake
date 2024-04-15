function(find_cppoptlib)
    find_path(CPPOPTLIB_INCLUDE_DIR
            cppoptlib/solver/lbfgs.h
            HINTS ${DMRG_DEPS_INSTALL_DIR}
            PATH_SUFFIXES cppoptlib include/cppoptlib cppoptlib/include
            NO_CMAKE_ENVIRONMENT_PATH
            NO_SYSTEM_ENVIRONMENT_PATH
            NO_CMAKE_SYSTEM_PATH
            )
endfunction()
find_cppoptlib()
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(cppoptlib DEFAULT_MSG CPPOPTLIB_INCLUDE_DIR)


if (cppoptlib_FOUND AND NOT TARGET cppoptlib::cppoptlib)
    message(DEBUG "Defining target cppoptlib::cppoptlib")
    add_library(cppoptlib::cppoptlib INTERFACE IMPORTED)
    target_include_directories(cppoptlib::cppoptlib SYSTEM INTERFACE ${CPPOPTLIB_INCLUDE_DIR})
    find_dependency(Eigen3 REQUIRED)
    target_link_libraries(cppoptlib::cppoptlib INTERFACE Eigen3::Eigen)
endif()