include(GNUInstallDirs)
find_package(spdlog 1.3 NO_DEFAULT_PATH PATHS ${INSTALL_DIRECTORY}/spdlog/${CMAKE_INSTALL_LIBDIR}/spdlog/cmake ${spdlog_DIR} )
if(spdlog_FOUND)
    get_target_property(spdlog_LIBRARIES spdlog::spdlog   IMPORTED_LOCATION_RELEASE)
    set_target_properties(spdlog::spdlog PROPERTIES INTERFACE_LINK_LIBRARIES "${spdlog_LIBRARIES}")
    message(STATUS "SPDLOG FOUND IN SYSTEM: ${spdlog_LIBRARIES}")


elseif (DOWNLOAD_SPDLOG OR DOWNLOAD_ALL)
    message(STATUS "Spdlog will be installed into ${INSTALL_DIRECTORY}/spdlog on first build.")
    include(ExternalProject)
    ExternalProject_Add(external_SPDLOG
            GIT_REPOSITORY https://github.com/gabime/spdlog.git
            GIT_TAG v1.x
            GIT_PROGRESS 1
            UPDATE_COMMAND ""
            TEST_COMMAND ""
            PREFIX      ${BUILD_DIRECTORY}/spdlog
            INSTALL_DIR ${INSTALL_DIRECTORY}/spdlog
            CMAKE_ARGS
            -DSPDLOG_BUILD_TESTS:BOOL=ON
            -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
            )


    ExternalProject_Get_Property(external_SPDLOG INSTALL_DIR)
    add_library(spdlog INTERFACE)
    add_library(spdlog::spdlog ALIAS spdlog)
    set(spdlog_DIR ${INSTALL_DIR}/${CMAKE_INSTALL_LIBDIR}/spdlog/cmake)

    add_dependencies(spdlog external_SPDLOG)
    target_include_directories(spdlog INTERFACE ${INSTALL_DIR}/${CMAKE_INSTALL_INCLUDEDIR})
    target_link_libraries(spdlog INTERFACE ${INSTALL_DIR}/${CMAKE_INSTALL_LIBDIR}/spdlog/libspdlog${CMAKE_STATIC_LIBRARY_SUFFIX})
    target_link_libraries (spdlog INTERFACE ${PTHREAD_LIBRARY})
else()
    message("WARNING: Dependency spdlog not found and DOWNLOAD_SPDLOG is OFF. Build will fail.")

endif()