
function(find_optim)
    find_library(OPTIM_LIBRARY
            optim
            HINTS ${DMRG_DEPS_INSTALL_DIR}
            PATH_SUFFIXES lib optim/lib
            NO_CMAKE_ENVIRONMENT_PATH
            NO_SYSTEM_ENVIRONMENT_PATH
            NO_CMAKE_SYSTEM_PATH
            )
    find_path(OPTIM_INCLUDE_DIR
            optim/optim.h
            HINTS ${DMRG_DEPS_INSTALL_DIR}
            PATH_SUFFIXES include optim/include
            NO_CMAKE_ENVIRONMENT_PATH
            NO_SYSTEM_ENVIRONMENT_PATH
            NO_CMAKE_SYSTEM_PATH
            )
endfunction()

find_optim()
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(optim
        DEFAULT_MSG
        OPTIM_LIBRARY OPTIM_INCLUDE_DIR)


if (optim_FOUND)
    add_library(optim::optim UNKNOWN IMPORTED)
    set_target_properties(optim::optim PROPERTIES IMPORTED_LOCATION "${OPTIM_LIBRARY}")
    target_include_directories(optim::optim SYSTEM INTERFACE ${OPTIM_INCLUDE_DIR})
    target_compile_definitions(optim::optim INTERFACE OPTIM_INT_SIZE=32)
endif()