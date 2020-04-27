find_package(Lapacke)                                           # Lapacke needed by arpack++
if(NOT TARGET arpack::arpack++ AND DMRG_DOWNLOAD_METHOD MATCHES "find|fetch|native|conan")
    include(cmake-modules/FindArpack++.cmake)
    find_Arpackpp()
    if (TARGET arpack::arpack++)
        message(STATUS "Arpack++ found")
    endif()
endif()

if(NOT TARGET arpack::arpack++ AND DMRG_DOWNLOAD_METHOD MATCHES "fetch|native|conan")
    message(STATUS "Arpack++ will be installed into ${CMAKE_BINARY_DIR}/dmrg-deps-install/arpack++ on first build.")
    include(${PROJECT_SOURCE_DIR}/cmake-modules/BuildDependency.cmake)
    build_dependency(arpack++ "${CMAKE_INSTALL_PREFIX}/arpack++" "")
    include(cmake-modules/FindArpack++.cmake)
    find_Arpackpp()
    if (TARGET arpack::arpack++)
        message(STATUS "Successfully installed Arpack++")
    else()
        message(FATAL_ERROR "arpack-ng could not be installed")
    endif()
endif()



