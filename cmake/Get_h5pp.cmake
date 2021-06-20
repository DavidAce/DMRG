
function(find_h5pp)
    if(DMRG_PACKAGE_MANAGER STREQUAL "find")
        set(H5PP_REQUIRED REQUIRED)
    endif()
    find_package(h5pp 1.8.5 HINTS ${DMRG_DEPS_INSTALL_DIR} ${H5PP_REQUIRED})

    if(NOT h5pp_FOUND AND DMRG_PACKAGE_MANAGER MATCHES "cmake")
        include(cmake/InstallPackage.cmake)
        list(APPEND H5PP_CMAKE_OPTIONS  -DEigen3_ROOT:PATH=${DMRG_DEPS_INSTALL_DIR})
        list(APPEND H5PP_CMAKE_OPTIONS  -DH5PP_PACKAGE_MANAGER:STRING=cmake)
        list(APPEND H5PP_CMAKE_OPTIONS  -DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE})
        install_package(h5pp "${DMRG_DEPS_INSTALL_DIR}" "${H5PP_CMAKE_OPTIONS}")
        find_package(h5pp 1.9.0 HINTS ${DMRG_DEPS_INSTALL_DIR} REQUIRED)
    endif()
endfunction()



find_h5pp()