
if(NOT TARGET glog::glog)
    message(STATUS "Searching for glog 0.4.x")
    find_package(glog 0.4 PATHS ${INSTALL_DIRECTORY}/glog $ENV{EBROOTGLOG} $ENV{GLOG_DIR} NO_DEFAULT_PATH)
    if(TARGET glog::glog)
        get_target_property(GLOG_LIBRARIES   glog::glog INTERFACE_LINK_LIBRARIES)
        message(STATUS "Searching for glog 0.4.x - Success: ${GLOG_LIBRARIES}")
    else()
        message(STATUS "Searching for glog 0.4.x - failed")
    endif()
endif()
if(TARGET glog::glog)
    get_target_property(GLOG_INCLUDE_DIR glog::glog INTERFACE_INCLUDE_DIRECTORIES)
    get_target_property(GLOG_LIBRARIES   glog::glog INTERFACE_LINK_LIBRARIES)
    set_target_properties(glog::glog PROPERTIES INTERFACE_LINK_LIBRARIES "${INSTALL_DIRECTORY}/glog/lib/libglog${CUSTOM_SUFFIX}")
    set_target_properties(glog::glog PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${GLOG_INCLUDE_DIR}")
    message(STATUS "glog FOUND IN SYSTEM: ${GLOG_INCLUDE_DIR}")
else()
    message(STATUS "glog will be installed into ${INSTALL_DIRECTORY}/glog on first build.")

    include(ExternalProject)
    ExternalProject_Add(external_GLOG
            GIT_REPOSITORY https://github.com/google/glog.git
            GIT_TAG master
            GIT_PROGRESS 1
            PREFIX      ${BUILD_DIRECTORY}/glog
            INSTALL_DIR ${INSTALL_DIRECTORY}/glog
            UPDATE_COMMAND ""
            CMAKE_ARGS
            -DCMAKE_BUILD_TYPE=Release
            -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
            -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
            )

    ExternalProject_Get_Property(external_GLOG INSTALL_DIR)
    add_library(glog INTERFACE)
    add_library(glog::glog ALIAS glog)
    set(GLOG_INCLUDE_DIR ${INSTALL_DIR}/include)
    add_dependencies(glog external_GLOG)
    target_link_libraries(glog INTERFACE ${INSTALL_DIR}/lib/libglog${CUSTOM_SUFFIX})
    target_include_directories(glog INTERFACE ${GLOG_INCLUDE_DIR})
endif()
