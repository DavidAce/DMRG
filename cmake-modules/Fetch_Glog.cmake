
include(cmake-modules/filterTarget.cmake)
string(TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE)
find_package(glog 0.4
        HINTS ${DIRECTORY_HINTS}
        PATHS $ENV{EBROOTGLOG} $ENV{GLOG_DIR} $ENV{glog_DIR}
        PATH_SUFFIXES glog glog/lib
        NO_DEFAULT_PATH)


if(TARGET glog::glog)
    message(STATUS "glog found")

    string(TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE)
    get_target_property(gloglib  glog::glog IMPORTED_LOCATION_${BUILD_TYPE})
    if(gloglib)
        target_link_libraries(glog::glog INTERFACE ${gloglib})
        remove_shared(glog::glog)
        remove_pthread(glog::glog)
    else()
        message(FATAL_ERROR "Dependency glog::glog does not have IMPORTED_LOCATION_${BUILD_TYPE}")
    endif()

#    include(cmake-modules/PrintTargetProperties.cmake)
#    print_target_properties(glog::glog)


elseif(DOWNLOAD_MISSING)
    message(STATUS "glog will be installed into ${CMAKE_INSTALL_PREFIX}")
    include(${PROJECT_SOURCE_DIR}/cmake-modules/BuildDependency.cmake)
    list(APPEND GLOG_CMAKE_OPTIONS -Dgflags_DIR:PATH=${CMAKE_INSTALL_PREFIX}/gflags)
    build_dependency(glog "${CMAKE_INSTALL_PREFIX}" "${GLOG_CMAKE_OPTIONS}")
    find_package(glog 0.4
            HINTS ${DIRECTORY_HINTS}
            PATHS ${CMAKE_INSTALL_PREFIX} $ENV{EBROOTGLOG} $ENV{GLOG_DIR} $ENV{glog_DIR}
            PATH_SUFFIXES glog glog/lib
            NO_DEFAULT_PATH)
    if(TARGET glog::glog)
        message(STATUS "glog installed successfully")
        get_target_property(gloglib  glog::glog IMPORTED_LOCATION_${BUILD_TYPE})
        if(gloglib)
            target_link_libraries(glog::glog INTERFACE ${gloglib})
            remove_shared(glog::glog)
            remove_pthread(glog::glog)
        else()
            message(FATAL_ERROR "Dependency glog::glog does not have IMPORTED_LOCATION_${BUILD_TYPE}")
        endif()
#        include(cmake-modules/PrintTargetProperties.cmake)
#        print_target_properties(glog::glog)
    else()
        message(STATUS "config_result: ${config_result}")
        message(STATUS "build_result: ${build_result}")
        message(FATAL_ERROR "glog could not be downloaded.")
    endif()

else()
    message(FATAL_ERROR "Dependency glog not found and DOWNLOAD_MISSING is OFF")
endif()