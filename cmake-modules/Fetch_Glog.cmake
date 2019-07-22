


message(STATUS "Searching for glog 0.4.x")
find_package(glog 0.4 PATHS ${INSTALL_DIRECTORY}/glog $ENV{EBROOTGLOG} $ENV{GLOG_DIR} $ENV{glog_DIR} NO_DEFAULT_PATH)

if(TARGET glog::glog)
    # For some reason glog imports libunwind.so even though we asked for static libraries.
    # In addition, libglog.a hides in the property "LOCATION" instead of its rightful
    # place "INTERFACE_LINK_LIBRARIES".

    get_target_property(GLOG_INCLUDE_DIR glog::glog INTERFACE_INCLUDE_DIRECTORIES)
    get_target_property(GLOG_LIBRARIES   glog::glog LOCATION)
    # In addition, the library may be .so
    get_filename_component(glog_dir ${GLOG_LIBRARIES} DIRECTORY)
    get_filename_component(glog_we  ${GLOG_LIBRARIES} NAME_WE)
    set(GLOG_LIBRARIES "${glog_dir}/${glog_we}${CUSTOM_SUFFIX}")

    set_target_properties(glog::glog PROPERTIES INTERFACE_LINK_LIBRARIES "${GLOG_LIBRARIES}")
    set_target_properties(glog::glog PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${GLOG_INCLUDE_DIR}")
    message(STATUS "Searching for glog 0.4.x - Success: LIB: ${GLOG_LIBRARIES}")
    message(STATUS "Searching for glog 0.4.x - Success: INC: ${GLOG_INCLUDE_DIR}")
else()
    message(STATUS "Searching for glog 0.4.x - failed")
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
            DEPENDS gflags::gflags
            )

    ExternalProject_Get_Property(external_GLOG INSTALL_DIR)
    add_library(glog INTERFACE)
    add_library(glog::glog ALIAS glog)
    set(GLOG_INCLUDE_DIR ${INSTALL_DIR}/include)
    set(glog_DIR ${INSTALL_DIR}/lib/cmake/glog)
    add_dependencies(glog external_GLOG gflags::gflags)
    target_link_libraries(glog INTERFACE ${INSTALL_DIR}/lib/libglog${CUSTOM_SUFFIX} gflags::gflags)
    target_include_directories(glog INTERFACE ${GLOG_INCLUDE_DIR})
endif()
