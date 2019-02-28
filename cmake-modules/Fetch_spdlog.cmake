
message(STATUS "Spdlog will be installed into ${INSTALL_DIRECTORY}/spdlog on first build.")

include(ExternalProject)
ExternalProject_Add(external_SPDLOG
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG v1.x
        GIT_PROGRESS 1
        UPDATE_COMMAND ""
        TEST_COMMAND ""
        PREFIX "${INSTALL_DIRECTORY}/spdlog"
        CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
        )


ExternalProject_Get_Property(external_SPDLOG INSTALL_DIR)
add_library(spdlog INTERFACE)
set(SPDLOG_INCLUDE_DIR ${INSTALL_DIR}/include)
message("SPDLOG INCLUDE DIR: ${SPDLOG_INCLUDE_DIR}" )
add_dependencies(spdlog external_SPDLOG)

set_target_properties(spdlog PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES   "${SPDLOG_INCLUDE_DIR}"
        )
