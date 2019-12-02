
find_package(gflags
        PATHS ${EXTERNAL_INSTALL_DIR}/gflags $ENV{EBROOTGFLAGS} $ENV{GFLAGS_DIR} $ENV{gflags_DIR}
        NO_DEFAULT_PATH)

if(TARGET gflags)
    message(STATUS "gflags found")
    get_target_property(gflaglib  gflags IMPORTED_LOCATION_RELEASE)
    set_target_properties(gflags PROPERTIES INTERFACE_LINK_LIBRARIES "${gflaglib}" )
#    include(cmake-modules/PrintTargetInfo.cmake)
#    print_target_info(gflags)
elseif(DOWNLOAD_MISSING)
    message(STATUS "gflags will be installed into ${EXTERNAL_INSTALL_DIR}/gflags")
    include(cmake-modules/BuildExternalLibs.cmake)
    build_external_libs(
            "gflags"
            "${EXTERNAL_CONFIG_DIR}"
            "${EXTERNAL_BUILD_DIR}"
            "${EXTERNAL_INSTALL_DIR}"
            ""
    )
    find_package(gflags PATHS ${EXTERNAL_INSTALL_DIR}/gflags NO_DEFAULT_PATH REQUIRED)
    if(TARGET gflags)
        message(STATUS "gflags installed successfully")
        get_target_property(gflaglib  gflags IMPORTED_LOCATION_RELEASE)
        set_target_properties(gflags PROPERTIES INTERFACE_LINK_LIBRARIES "${gflaglib}" )
#        include(cmake-modules/PrintTargetInfo.cmake)
#        print_target_info(gflags)
    else()
        message(STATUS "config_result: ${config_result}")
        message(STATUS "build_result: ${build_result}")
        message(FATAL_ERROR "gflags could not be downloaded.")
    endif()

else()
    message(FATAL_ERROR "Dependency gflags not found and DOWNLOAD_MISSING is OFF")
endif()