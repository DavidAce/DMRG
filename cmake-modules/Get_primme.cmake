function(find_primme)
    unset(PRIMME_LIBRARY)
    unset(PRIMME_LIBRARY CACHE)
    unset(PRIMME_INCLUDE_DIR)
    unset(PRIMME_INCLUDE_DIR CACHE)
    find_library(PRIMME_LIBRARY
            primme
            HINTS ${CMAKE_INSTALL_PREFIX}/primme
            PATH_SUFFIXES lib
            NO_SYSTEM_PATH
            )
    find_path(PRIMME_INCLUDE_DIR
            primme/primme.h
            HINTS ${CMAKE_INSTALL_PREFIX}/primme
            PATH_SUFFIXES include
            NO_SYSTEM_PATH
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


if(NOT TARGET primme::primme AND DMRG_PACKAGE_MANAGER STREQUAL "find")
    find_primme()
endif()

if(NOT TARGET primme::primme AND DMRG_PACKAGE_MANAGER MATCHES "cmake|conan")
    message(STATUS "primme will be installed into ${CMAKE_BINARY_DIR}/dmrg-deps-install/primme on first build.")
    include(${PROJECT_SOURCE_DIR}/cmake-modules/BuildDependency.cmake)
    build_dependency(primme "${CMAKE_INSTALL_PREFIX}" "")
    find_primme()
    if (TARGET primme::primme)
        message(STATUS "Successfully installed primme")
    else()
        message(FATAL_ERROR "Could not install primme")
    endif()
endif()
