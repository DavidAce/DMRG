include(cmake-modules/filterTarget.cmake)
find_package(gflags
        PATHS ${CMAKE_INSTALL_PREFIX}/gflags $ENV{EBROOTGFLAGS} $ENV{GFLAGS_DIR} $ENV{gflags_DIR}
        NO_DEFAULT_PATH)

if(TARGET gflags)
    message(STATUS "gflags found")
    #Copy the lib to where it belongs: INTERFACE_LINK_LIBRARIES
    string(TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE)
    get_target_property(gflaglib  gflags IMPORTED_LOCATION_${BUILD_TYPE})
    target_link_libraries(gflags INTERFACE ${gflaglib})
    remove_shared(gflags)
    remove_pthread(gflags)
elseif(DOWNLOAD_MISSING)
    message(STATUS "gflags will be installed into ${CMAKE_INSTALL_PREFIX}")
    include(${PROJECT_SOURCE_DIR}/cmake-modules/BuildDependency.cmake)
    build_dependency(gflags "")
    find_package(gflags HINTS ${CMAKE_INSTALL_PREFIX}/gflags NO_DEFAULT_PATH)
    if(TARGET gflags)
        message(STATUS "gflags installed successfully")
        #Copy the lib to where it belongs: INTERFACE_LINK_LIBRARIES
        get_target_property(gflaglib  gflags IMPORTED_LOCATION_${BUILD_TYPE})
        target_link_libraries(gflags INTERFACE ${gflaglib})
        remove_shared(gflags)
        remove_pthread(gflags)
    else()
        message(STATUS "config_result: ${config_result}")
        message(STATUS "build_result: ${build_result}")
        message(FATAL_ERROR "gflags could not be downloaded.")
    endif()

else()
    message(FATAL_ERROR "Dependency gflags not found and DOWNLOAD_MISSING is OFF")
endif()