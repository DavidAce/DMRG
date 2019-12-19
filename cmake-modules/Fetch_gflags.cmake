
find_package(gflags COMPONENTS nothreads_static
        HINTS $ENV{GFLAGS_DIR} $ENV{gflags_DIR} ${CMAKE_INSTALL_PREFIX}
        PATHS $ENV{EBROOTGFLAGS}
        PATH_SUFFIXES gflags gflags/lib)

if(NOT TARGET gflags)
    if(BUILD_SHARED_LIBS)
        find_package(gflags
                HINTS ${DIRECTORY_HINTS}
                PATHS ${CMAKE_INSTALL_PREFIX} $ENV{EBROOTGFLAGS} $ENV{GFLAGS_DIR} $ENV{gflags_DIR} $ENV{CONDA_PREFIX}
                PATH_SUFFIXES gflags gflags/lib)
    else()
        message(STATUS "Skipping search through conda libs because this is a static build")
    endif()
endif()
if(TARGET gflags)
    message(STATUS "gflags found")

elseif(DOWNLOAD_MISSING)
    message(STATUS "gflags will be installed into ${CMAKE_INSTALL_PREFIX}")
    include(${PROJECT_SOURCE_DIR}/cmake-modules/BuildDependency.cmake)
    build_dependency(gflags "${CMAKE_INSTALL_PREFIX}" "")
    find_package(gflags
            HINTS ${CMAKE_INSTALL_PREFIX}
            PATH_SUFFIXES gflags gflags/lib
            NO_DEFAULT_PATH)
    if(TARGET gflags)
        message(STATUS "gflags installed successfully")
    else()
        message(FATAL_ERROR "gflags could not be downloaded.")
    endif()

else()
    message(FATAL_ERROR "Dependency gflags not found and DOWNLOAD_MISSING is OFF")
endif()

if(TARGET gflags)
    include(cmake-modules/filterTarget.cmake)
    remove_shared(gflags)
    remove_pthread(gflags)
endif()