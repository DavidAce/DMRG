


include(cmake-modules/FindArpack++.cmake)
find_Arpackpp()

if (TARGET arpack::arpack++)
    message(STATUS "Arpack++ found")
elseif(NOT ${DMRG_DOWNLOAD_METHOD} MATCHES "none")
    message(STATUS "Arpack++ will be installed into ${CMAKE_BINARY_DIR}/dmrg-deps-install/arpack++ on first build.")

    include(${PROJECT_SOURCE_DIR}/cmake-modules/BuildDependency.cmake)
    build_dependency(arpack++ "${CMAKE_INSTALL_PREFIX}/arpack++" "")
    include(cmake-modules/FindArpack++.cmake)
    find_Arpackpp()
    if (TARGET arpack::arpack++)
        message(STATUS "Successfully installed Arpack++")
    else()
        message(FATAL_ERROR "arpack-ng could not be downloaded.")
    endif()

else()
    message(FATAL_ERROR "Dependency Arpack++ not found and DMRG_DOWNLOAD_METHOD = ${DMRG_DOWNLOAD_METHOD}")
endif()



