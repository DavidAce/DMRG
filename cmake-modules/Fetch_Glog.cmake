
#find_package(glog 0.4 PATHS ${INSTALL_DIRECTORY}/glog $ENV{EBROOTGLOG} $ENV{GLOG_DIR} $ENV{glog_DIR} NO_DEFAULT_PATH)
message(STATUS "Searching for glog")
#find_library(
#        GLOG_LIBRARIES
#        NAMES libglog${CUSTOM_SUFFIX}
#        PATH_SUFFIXES lib lib64
#        PATHS  ${INSTALL_DIRECTORY}/glog $ENV{EBROOTGLOG} $ENV{GLOG_DIR} $ENV{glog_DIR}
#        NO_DEFAULT_PATH)
#
#find_path(
#        GLOG_INCLUDE_DIR
#        NAMES glog/logging.h
#        PATH_SUFFIXES include glog/include
#        PATHS  ${INSTALL_DIRECTORY}/glog $ENV{EBROOTGLOG} $ENV{GLOG_DIR} $ENV{glog_DIR}
#        NO_DEFAULT_PATH)


if(GLOG_LIBRARIES AND GLOG_INCLUDE_DIR)
    # For some reason glog imports libunwind.so even though we asked for static libraries.
    # In addition, libglog.a hides in the property "LOCATION" instead of its rightful
    # place "INTERFACE_LINK_LIBRARIES"... so we do this manually
    add_library(glog INTERFACE)
    add_dependencies(glog gflags)
    get_filename_component(GLOG_LIBRARY_DIR ${GLOG_LIBRARIES} DIRECTORY)
    set_target_properties(glog PROPERTIES INTERFACE_LINK_LIBRARIES ${GLOG_LIBRARIES})
    set_target_properties(glog PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${GLOG_INCLUDE_DIR})
    target_link_libraries(glog INTERFACE gflags Threads::Threads)
    message(STATUS "Searching for glog - Success: LIB: ${GLOG_LIBRARIES}")
    message(STATUS "Searching for glog - Success: INC: ${GLOG_INCLUDE_DIR}")
    include(cmake-modules/PrintTargetProperties.cmake)
    print_target_properties(glog)
else()
    message(STATUS "Searching for glog - failed")
    message(STATUS "glog will be installed into ${INSTALL_DIRECTORY}/glog on first build.")
    unset(FLAGS CACHE)
    if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
        set(FLAGS "${FLAGS} -stdlib=libstdc++ ${GCC_TOOLCHAIN}")
    endif()
    include(ExternalProject)
    ExternalProject_Add(external_GLOG
            GIT_REPOSITORY https://github.com/google/glog.git
            GIT_TAG master
            GIT_PROGRESS false
            GIT_SHALLOW true
            PREFIX      ${BUILD_DIRECTORY}/glog
            INSTALL_DIR ${INSTALL_DIRECTORY}/glog
            UPDATE_COMMAND ""
            CMAKE_ARGS
            -DCMAKE_BUILD_TYPE=Release
            -DCMAKE_CXX_FLAGS:STRING=${FLAGS}
            -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
            -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
            DEPENDS gflags
            )

    ExternalProject_Get_Property(external_GLOG INSTALL_DIR)
    add_library(glog INTERFACE)
    include(GNUInstallDirs)
    set(GLOG_INCLUDE_DIR ${INSTALL_DIR}/${CMAKE_INSTALL_INCLUDEDIR})
    set(GLOG_LIBRARY_DIR ${INSTALL_DIR}/${CMAKE_INSTALL_LIBDIR})
    set(GLOG_LIBRARIES   ${GLOG_LIBRARY_DIR}/libglog${CUSTOM_SUFFIX})
    set(glog_DIR ${GLOG_LIBRARY_DIR}/cmake/glog)
    add_dependencies(glog external_GLOG gflags)
    target_link_libraries(glog INTERFACE ${GLOG_LIBRARIES} gflags Threads::Threads)
    target_include_directories(glog INTERFACE ${GLOG_INCLUDE_DIR})
endif()


