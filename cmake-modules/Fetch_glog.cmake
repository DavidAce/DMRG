include(cmake-modules/CheckGlogCompiles.cmake)

    # Glog should only look in conda on shared builds!
set(GLOG_HINTS $ENV{EBROOTGLOG} ${CMAKE_INSTALL_PREFIX} )
if(BUILD_SHARED_LIBS)
    list(APPEND GLOG_HINTS ${CONDA_HINTS})
endif()

find_package(glog 0.4 HINTS ${GLOG_HINTS} PATH_SUFFIXES glog glog/lib)


if(NOT TARGET glog::glog)
    message(STATUS "Looking for glog in system")
    find_library(GLOG_LIBRARIES     glog           HINTS ${GLOG_HINTS})
    find_path   (GLOG_INCLUDE_DIR   glog/logging.h HINTS ${GLOG_HINTS})
    check_glog_compiles("lib_header" "" "" "${GLOG_LIBRARIES}" "${GLOG_INCLUDE_DIR}" "gflags::gflags")
    if (GLOG_COMPILES_lib_header)
        add_library(glog::glog ${LINK_TYPE} IMPORTED)
        set_target_properties(glog::glog PROPERTIES IMPORTED_LOCATION "${GLOG_LIBRARIES}")
        target_include_directories(glog::glog SYSTEM INTERFACE ${GLOG_INCLUDE_DIR})
        target_link_libraries(glog::glog INTERFACE gcc_eh unwind lzma)
        message(STATUS "Found system glog: Don't forget to also install and link to libraries unwind and lzma")
    endif()
endif()


if(TARGET glog::glog)
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
    get_target_property(GLOG_TYPE glog::glog TYPE)
    if(NOT GLOG_TYPE MATCHES "${LINK_TYPE}")
        include(cmake-modules/PrintTargetProperties.cmake)
        print_target_properties(glog::glog)
        message(FATAL_ERROR "Found shared glog library on a static build!")
    endif()

    include(cmake-modules/filterTarget.cmake)
    if(NOT BUILD_SHARED_LIBS)
        remove_pthread_shallow(glog::glog)
        remove_shared(glog::glog)
    endif()

    # Modernize
    get_property(imp_loc_set TARGET glog::glog PROPERTY IMPORTED_LOCATION SET) # Returns a boolean if set
    get_property(loc_set     TARGET glog::glog PROPERTY LOCATION SET) # Returns a boolean if set
    if(loc_set AND NOT imp_loc_set)
        get_target_property(imp_loc glog::glog LOCATION)
        set_target_properties(glog::glog PROPERTIES IMPORTED_LOCATION ${imp_loc})
    endif()
    if(TARGET gflags::gflags)
        target_link_libraries(glog::glog INTERFACE gflags::gflags)
    endif()
    check_glog_compiles("target" "-std=c++17" "" "" "" "glog::glog")
    if(NOT GLOG_COMPILES_target)
        include(cmake-modules/PrintTargetProperties.cmake)
        print_target_properties(glog::glog)
        message(FATAL_ERROR "Could not compile a simple glog program")
    endif()
endif()