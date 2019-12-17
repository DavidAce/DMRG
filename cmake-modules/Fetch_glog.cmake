


find_package(glog 0.4
        HINTS $ENV{GLOG_DIR} $ENV{glog_DIR} ${CMAKE_INSTALL_PREFIX}
        PATHS $ENV{EBROOTGLOG}
        PATH_SUFFIXES glog glog/lib
        NO_DEFAULT_PATH)

if(NOT TARGET glog::glog)
    if(BUILD_SHARED_LIBS)
        find_package(glog 0.4
                HINTS ${DIRECTORY_HINTS}
                PATHS $ENV{EBROOTGLOG} $ENV{GLOG_DIR} $ENV{glog_DIR} $ENV{CONDA_PREFIX}
                PATH_SUFFIXES glog glog/lib)
    else()
        message(STATUS "Skipping search through conda libs because this is a static build")
    endif()
endif()

if(TARGET glog::glog)
    message(STATUS "glog found")


elseif(DOWNLOAD_MISSING)
    message(STATUS "glog will be installed into ${CMAKE_INSTALL_PREFIX}")
    include(${PROJECT_SOURCE_DIR}/cmake-modules/BuildDependency.cmake)
    list(APPEND GLOG_CMAKE_OPTIONS -Dgflags_DIR:PATH=${CMAKE_INSTALL_PREFIX}/gflags)
    build_dependency(glog "${CMAKE_INSTALL_PREFIX}" "${GLOG_CMAKE_OPTIONS}")
    find_package(glog 0.4
            HINTS ${CMAKE_INSTALL_PREFIX}
            PATH_SUFFIXES glog glog/lib
            NO_DEFAULT_PATH)
    if(TARGET glog::glog)
        message(STATUS "glog installed successfully")
    else()
        message(FATAL_ERROR "glog could not be downloaded.")
    endif()

else()
    message(FATAL_ERROR "Dependency glog not found and DOWNLOAD_MISSING is OFF")
endif()


if(TARGET glog::glog)
    #Copy the lib to where it belongs: INTERFACE_LINK_LIBRARIES
    include(cmake-modules/filterTarget.cmake)
    string(TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE)
    get_target_property(glog_imported_loc_buildtype   glog::glog IMPORTED_LOCATION_${BUILD_TYPE})
    get_target_property(glog_imported_loc_noconfig    glog::glog IMPORTED_LOCATION_NOCONFIG)
    if(glog_imported_loc_buildtype)
        target_link_libraries(glog::glog INTERFACE ${glog_imported_loc_buildtype})
    elseif(glog_imported_loc_noconfig)
        target_link_libraries(glog::glog INTERFACE ${glog_imported_loc_noconfig})
    else()
        message(STATUS "Dependency glog does not have IMPORTED_LOCATION_${BUILD_TYPE}/_NOCONFIG")
    endif()

    remove_shared(glog::glog)
    remove_pthread(glog::glog)
endif()