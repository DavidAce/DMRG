
# Can't use conda/apt here since they only have shared libraries.
# We also can't use apt-versions since they hardcode usage of Eigen3/Glog/Gflags
# but we need our own patched Eigen3.
find_package(Ceres
        HINTS $ENV{CERES_DIR} $ENV{ceres_DIR} ${CERES_DIR} ${ceres_DIR} ${CMAKE_INSTALL_PREFIX}/ceres
        PATH_SUFFIXES ceres ceres/lib
        NO_DEFAULT_PATH)

if(NOT TARGET ceres)
    if(BUILD_SHARED_LIBS)
    find_package(Ceres
            HINTS  ${DIRECTORY_HINTS}
            PATHS ${CMAKE_INSTALL_PREFIX} $ENV{EBROOTCERES} $ENV{CERES_DIR} $ENV{ceres_DIR} $ENV{CONDA_PREFIX}
            PATH_SUFFIXES ceres ceres/lib)
#            NO_DEFAULT_PATH)
    else()
        message(STATUS "Skipping search through conda libs because this is a static build")
    endif()
endif()

if(TARGET ceres)
    message(STATUS "ceres found")
elseif(DOWNLOAD_MISSING)
    message(STATUS "Ceres will be installed into ${CMAKE_INSTALL_PREFIX} on first build.")
    get_target_property(EIGEN3_INCLUDE_DIR Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
    list (GET EIGEN3_INCLUDE_DIR 0 EIGEN3_INCLUDE_DIR)
    get_target_property(GFLAGS_TYPE    gflags      TYPE)
    get_target_property(GLOG_TYPE      glog::glog  TYPE)
    set(LIB_PROPERTIES INTERFACE_LINK_LIBRARIES IMPORTED_LOCATION_${BUILD_TYPE} IMPORTED_LOCATION_NOCONFIG LOCATION LOCATION_${CMAKE_BUILD_TYPE})
    foreach(prop ${LIB_PROPERTIES})
        if(GLOG_TYPE MATCHES "INTERFACE" AND prop MATCHES "IMPORTED|LOCATION")
            continue()
        endif()
        get_target_property(GLOG_LIBRARIES glog::glog ${prop})
        if(GLOG_LIBRARIES)
            break()
        endif()
    endforeach()
    foreach(prop ${LIB_PROPERTIES})
        if(GFLAGS_TYPE MATCHES "INTERFACE" AND prop MATCHES "IMPORTED|LOCATION")
            continue()
        endif()
        get_target_property(GFLAGS_LIBRARIES gflags ${prop})
        message("${GFLAGS_LIBRARIES}")
        if(GFLAGS_LIBRARIES)
            break()
        endif()
    endforeach()
    get_target_property(GFLAGS_INCLUDE_DIR  gflags      INTERFACE_INCLUDE_DIRECTORIES)
    get_target_property(GLOG_INCLUDE_DIR    glog::glog  INTERFACE_INCLUDE_DIRECTORIES)
    message("GLOG_LIBRARIES   : ${GLOG_LIBRARIES}")
    message("GFLAGS_LIBRARIES : ${GFLAGS_LIBRARIES}")
    unset(CERES_FLAGS CACHE)
    unset(CERES_FLAGS)
    if("${CMAKE_BUILD_TYPE}" MATCHES "Debug")
        set(CERES_FLAGS "${CERES_FLAGS} -O0 -g3 -fstack-protector -D_FORTIFY_SOURCE=2")
    elseif("${CMAKE_BUILD_TYPE}" MATCHES "RelWithDebInfo")
        set(CERES_FLAGS "${CERES_FLAGS} -O1 -g -fstack-protector  -D_FORTIFY_SOURCE=2")
    else()
        set(CERES_FLAGS "${CERES_FLAGS} -O3 -g -DNDEBUG")
    endif()
    # It is VERY IMPORTANT to set this flag when using CERES
    # For some reason the Eigen destructor starts calling std::free instead of
    # the handmade allocator free in calls inside CERES. I do not know the reason.
    # Setting this flag here or globally helps.
#    set(CERES_FLAGS "${CERES_FLAGS} -DEIGEN_MALLOC_ALREADY_ALIGNED=0")

#    set(CERES_FLAGS "${CERES_FLAGS} -I${GLOG_INCLUDE_DIR} -I${GFLAGS_INCLUDE_DIR} -DEIGEN_MALLOC_ALREADY_ALIGNED=0")
    string (REPLACE ";" " " CERES_FLAGS "${CERES_FLAGS}")
    list(APPEND CERES_CMAKE_OPTIONS  -DCMAKE_CXX_FLAGS:STRING=${CERES_FLAGS})
    list(APPEND CERES_CMAKE_OPTIONS  -DEIGEN_INCLUDE_DIR_HINTS:PATH=${CMAKE_INSTALL_PREFIX})
    list(APPEND CERES_CMAKE_OPTIONS  -DEigen3_DIR:PATH=${CMAKE_INSTALL_PREFIX}/Eigen3/share/eigen3/cmake)
    list(APPEND CERES_CMAKE_OPTIONS  -Dgflags_DIR:PATH=${CMAKE_INSTALL_PREFIX}/gflags)
    list(APPEND CERES_CMAKE_OPTIONS  -DGFLAGS_INCLUDE_DIR_HINTS:PATH=${CMAKE_INSTALL_PREFIX})
    list(APPEND CERES_CMAKE_OPTIONS  -DGFLAGS_LIBRARY_DIR_HINTS:PATH=${CMAKE_INSTALL_PREFIX})
    list(APPEND CERES_CMAKE_OPTIONS  -DGFLAGS_LIBRARY:PATH=${GFLAGS_LIBRARIES})
    list(APPEND CERES_CMAKE_OPTIONS  -Dglog_DIR:PATH=${CMAKE_INSTALL_PREFIX}/glog)
    list(APPEND CERES_CMAKE_OPTIONS  -DGLOG_INCLUDE_DIR_HINTS:PATH=${CMAKE_INSTALL_PREFIX})
    list(APPEND CERES_CMAKE_OPTIONS  -DGLOG_LIBRARY_DIR_HINTS:PATH=${CMAKE_INSTALL_PREFIX})
    list(APPEND CERES_CMAKE_OPTIONS  -DGLOG_LIBRARY:PATH=${GLOG_LIBRARIES})

#    message(STATUS "ceres flags  : ${CERES_FLAGS}")
    message(STATUS "ceres options: ${CERES_CMAKE_OPTIONS}")
    include(${PROJECT_SOURCE_DIR}/cmake-modules/BuildDependency.cmake)
    build_dependency(ceres "${CMAKE_INSTALL_PREFIX}" "${CERES_CMAKE_OPTIONS}" )
    find_package(Ceres HINTS ${CMAKE_INSTALL_PREFIX}/ceres NO_DEFAULT_PATH)
    if(TARGET ceres)
        message(STATUS "ceres installed successfully")
    else()
        message(FATAL_ERROR "ceres could not be downloaded.")
    endif()
else()
    message(FATAL_ERROR "Dependency ceres not found and DOWNLOAD_MISSING is OFF")
endif()


if(TARGET ceres)
    include(cmake-modules/filterTarget.cmake)
    remove_shared (ceres)
    remove_pthread_shallow(ceres)
    target_link_libraries(ceres INTERFACE Eigen3::Eigen glog::glog gflags)
endif()