


find_package(glog 0.4
        HINTS $ENV{GLOG_DIR} $ENV{glog_DIR} ${CMAKE_INSTALL_PREFIX}
        PATHS $ENV{EBROOTGLOG}
        PATH_SUFFIXES glog glog/lib)

if(NOT TARGET glog::glog)
    find_library(GLOG_LIBRARIES     glog           HINTS $ENV{GLOG_DIR} $ENV{glog_DIR} ${CMAKE_INSTALL_PREFIX} ${DIRECTORY_HINTS})
    find_path   (GLOG_INCLUDE_DIR   glog/logging.h HINTS $ENV{GLOG_DIR} $ENV{glog_DIR} ${CMAKE_INSTALL_PREFIX} ${DIRECTORY_HINTS})
    if (GLOG_LIBRARIES AND GLOG_INCLUDE_DIR)
        add_library(glog::glog UNKNOWN IMPORTED)
        string(TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE)
        set_target_properties(glog::glog PROPERTIES
                IMPORTED_LOCATION_${BUILD_TYPE} "${GLOG_LIBRARIES}"
                IMPORTED_LINK_INTERFACE_LIBRARIES "gcc_eh;unwind;lzma"
                INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${GLOG_INCLUDE_DIR}")
#        target_link_libraries(glog::glog INTERFACE gcc_eh unwind lzma)
        message(STATUS "Found system glog: Don't forget to also install and link to libraries unwind and lzma")
    endif()
endif()



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
    include(cmake-modules/filterTarget.cmake)
    remove_shared(glog::glog)
    remove_pthread_shallow(glog::glog)
endif()