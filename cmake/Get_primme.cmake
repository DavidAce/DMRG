function(find_primme)
    unset(PRIMME_LIBRARY)
    unset(PRIMME_LIBRARY CACHE)
    unset(PRIMME_INCLUDE_DIR)
    unset(PRIMME_INCLUDE_DIR CACHE)
    find_library(PRIMME_LIBRARY
            primme
            HINTS ${DMRG_DEPS_INSTALL_DIR}
            PATH_SUFFIXES lib primme/lib
            NO_CMAKE_ENVIRONMENT_PATH
            NO_SYSTEM_ENVIRONMENT_PATH
            NO_CMAKE_SYSTEM_PATH
            )
    find_path(PRIMME_INCLUDE_DIR
            primme/primme.h
            HINTS ${DMRG_DEPS_INSTALL_DIR}
            PATH_SUFFIXES include primme/include
            NO_CMAKE_ENVIRONMENT_PATH
            NO_SYSTEM_ENVIRONMENT_PATH
            NO_CMAKE_SYSTEM_PATH
            )

    if (PRIMME_LIBRARY AND PRIMME_INCLUDE_DIR)
        add_library(primme::primme UNKNOWN IMPORTED)
        set_target_properties(primme::primme PROPERTIES IMPORTED_LOCATION "${PRIMME_LIBRARY}")
        target_include_directories(primme::primme SYSTEM INTERFACE ${PRIMME_INCLUDE_DIR})
        target_compile_definitions(primme::primme INTERFACE PRIMME_INT_SIZE=32)
        message(STATUS "Found primme: ${PRIMME_LIBRARY} ${PRIMME_INCLUDE_DIR}")
    else()
        message(STATUS "Could not find primme: ${PRIMME_LIBRARY} ${PRIMME_INCLUDE_DIR}")
    endif()

endfunction()

find_primme()
if(NOT TARGET primme::primme AND DMRG_PACKAGE_MANAGER MATCHES "cmake|conan")
    message(STATUS "primme will be installed into ${CMAKE_BINARY_DIR}/dmrg-deps-install/primme on first build.")
    include(cmake/InstallPackage.cmake)
    install_package(primme "${DMRG_DEPS_INSTALL_DIR}" "")
    find_primme()
    if (TARGET primme::primme)
        message(STATUS "Successfully installed primme")
    else()
        message(FATAL_ERROR "Could not install primme")
    endif()
endif()
