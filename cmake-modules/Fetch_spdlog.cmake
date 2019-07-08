find_package(spdlog 1.3 NO_DEFAULT_PATH PATHS ${H5PP_INSTALL_DIR_THIRD_PARTY}/spdlog/lib/spdlog/cmake ${spdlog_DIR} )
if(spdlog_FOUND)
    get_target_property(spdlog_lib     spdlog::spdlog   INTERFACE_LINK_LIBRARIES)
    message(STATUS "SPDLOG FOUND IN SYSTEM: ${spdlog_lib}")
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
    set(spdlog_DIR ${INSTALL_DIR}/lib/cmake/spdlog)
    add_dependencies(spdlog external_SPDLOG)
    target_include_directories(spdlog INTERFACE ${INSTALL_DIR}/include)
    target_link_libraries(spdlog INTERFACE ${INSTALL_DIR}/lib/spdlog/libspdlog${CUSTOM_SUFFIX})
    target_link_libraries (spdlog INTERFACE ${PTHREAD_LIBRARY})
else()
    message("WARNING: Dependency spdlog not found and DOWNLOAD_SPDLOG is OFF. Build will fail.")

endif()