include(cmake-modules/filterTarget.cmake)
string(TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE)

find_package(gflags
        HINTS ${DIRECTORY_HINTS}
        PATHS ${CMAKE_INSTALL_PREFIX} $ENV{EBROOTGFLAGS} $ENV{GFLAGS_DIR} $ENV{gflags_DIR}
        PATH_SUFFIXES gflags gflags/lib
        NO_DEFAULT_PATH)

if(TARGET gflags)
    message(STATUS "gflags found")
    #Copy the lib to where it belongs: INTERFACE_LINK_LIBRARIES
    get_target_property(gflags_imported_loc   gflags IMPORTED_LOCATION_${BUILD_TYPE})
    get_target_property(gflags_interface_lib  gflags INTERFACE_LINK_LIBRARIES)
    if(gflags_imported_loc)
        target_link_libraries(gflags INTERFACE ${gflags_imported_loc})
    elseif(gflags_interface_lib)
        message("gflags has an interface link library: ${gflags_interface_lib}")
    else()
        message(FATAL_ERROR "Dependency gflags does not have IMPORTED_LOCATION_${BUILD_TYPE}")
    endif()
    remove_shared(gflags)
    remove_pthread(gflags)

elseif(DOWNLOAD_MISSING)
    message(STATUS "gflags will be installed into ${CMAKE_INSTALL_PREFIX}")
    include(${PROJECT_SOURCE_DIR}/cmake-modules/BuildDependency.cmake)
    build_dependency(gflags "${CMAKE_INSTALL_PREFIX}" "")
    find_package(gflags HINTS ${DIRECTORY_HINTS}
            PATHS ${CMAKE_INSTALL_PREFIX} $ENV{EBROOTGFLAGS} $ENV{GFLAGS_DIR} $ENV{gflags_DIR}
            PATH_SUFFIXES gflags gflags/lib
            NO_DEFAULT_PATH)
    if(TARGET gflags)
        message(STATUS "gflags installed successfully")
        #Copy the lib to where it belongs: INTERFACE_LINK_LIBRARIES
        get_target_property(gflags_imported_loc   gflags IMPORTED_LOCATION_${BUILD_TYPE})
        get_target_property(gflags_interface_lib  gflags INTERFACE_LINK_LIBRARIES)
        if(gflags_imported_loc)
            target_link_libraries(gflags INTERFACE ${gflags_imported_loc})
        elseif(gflags_interface_lib)
            message("gflags has an interface link library: ${gflags_interface_lib}")
        else()
            message(FATAL_ERROR "Dependency gflags does not have IMPORTED_LOCATION_${BUILD_TYPE}")
        endif()
        remove_shared(gflags)
        remove_pthread(gflags)
    else()
        message(FATAL_ERROR "gflags could not be downloaded.")
    endif()

else()
    message(FATAL_ERROR "Dependency gflags not found and DOWNLOAD_MISSING is OFF")
endif()