
function(find_h5pp)
    if(DMRG_PACKAGE_MANAGER STREQUAL "find")
        set(REQUIRED REQUIRED)
    endif()
    if(DMRG_PACKAGE_MANAGER STREQUAL "cmake")
        set(NO_DEFAULT_PATH NO_DEFAULT_PATH)
    endif()
    set(HDF5_FIND_VERBOSE ON)
    find_package(h5pp 1.9.0
            HINTS ${h5pp_ROOT} ${DMRG_DEPS_INSTALL_DIR}
            ${NO_DEFAULT_PATH}
            ${REQUIRED})
    if(NOT h5pp_FOUND AND DMRG_PACKAGE_MANAGER MATCHES "cmake")
        include(cmake/InstallPackage.cmake)
        list(APPEND H5PP_CMAKE_OPTIONS  -DEigen3_ROOT:PATH=${DMRG_DEPS_INSTALL_DIR})
        list(APPEND H5PP_CMAKE_OPTIONS  -DH5PP_PACKAGE_MANAGER:STRING=cmake)
        list(APPEND H5PP_CMAKE_OPTIONS  -DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE})
        install_package(h5pp "${DMRG_DEPS_INSTALL_DIR}" "${H5PP_CMAKE_OPTIONS}")
        find_package(h5pp 1.9.0
                HINTS ${h5pp_ROOT} ${DMRG_DEPS_INSTALL_DIR}
                NO_DEFAULT_PATH REQUIRED)
    endif()
endfunction()



find_h5pp()