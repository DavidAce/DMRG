
function(remove_shared target_name)
    get_target_property(target_libs ${target_name} INTERFACE_LINK_LIBRARIES)
    foreach(lib ${target_libs})
        if(NOT ${lib} MATCHES ".so")
            list(APPEND static_libs ${lib})
        endif()
    endforeach()
    set_target_properties(${target_name} PROPERTIES INTERFACE_LINK_LIBRARIES "${static_libs}")
endfunction()

find_package(glog 0.4 PATHS ${EXTERNAL_INSTALL_DIR}/glog $ENV{EBROOTGLOG} $ENV{GLOG_DIR} $ENV{glog_DIR} NO_DEFAULT_PATH)


if(TARGET glog::glog)
    message(STATUS "glog found")
    remove_shared(glog::glog)
    get_target_property(movelib  glog::glog IMPORTED_LOCATION_RELEASE)
    set_target_properties(glog::glog PROPERTIES INTERFACE_LINK_LIBRARIES "${movelib}")
#    include(cmake-modules/PrintTargetInfo.cmake)
#    print_target_info(glog::glog)


elseif(DOWNLOAD_MISSING)
    message(STATUS "glog will be installed into ${EXTERNAL_INSTALL_DIR}/glog")
    include(cmake-modules/BuildExternalLibs.cmake)
    build_external_libs(
            "glog"
            "${EXTERNAL_CONFIG_DIR}"
            "${EXTERNAL_BUILD_DIR}"
            "${EXTERNAL_INSTALL_DIR}"
            ""
    )
    find_package(glog 0.4 PATHS ${EXTERNAL_INSTALL_DIR}/glog NO_DEFAULT_PATH REQUIRED)
    if(TARGET glog::glog)
        message(STATUS "glog installed successfully")
        remove_shared(glog::glog)
        get_target_property(movelib  glog::glog IMPORTED_LOCATION_RELEASE)
        set_target_properties(glog::glog PROPERTIES INTERFACE_LINK_LIBRARIES "${movelib}")
#        include(cmake-modules/PrintTargetInfo.cmake)
#        print_target_info(glog::glog)
    else()
        message(STATUS "config_result: ${config_result}")
        message(STATUS "build_result: ${build_result}")
        message(FATAL_ERROR "glog could not be downloaded.")
    endif()

else()
    message(FATAL_ERROR "Dependency glog not found and DOWNLOAD_MISSING is OFF")
endif()