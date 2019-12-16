
# Can't use conda here since it only has shared libraries.
# We also can't use apt-versions since they hardcode usage of Eigen3/Glog/Gflags
# but we need our own patched Eigen3.
find_package(Ceres
        HINTS $ENV{CERES_DIR} $ENV{ceres_DIR} ${CERES_DIR} ${ceres_DIR} ${CMAKE_INSTALL_PREFIX}/ceres
        PATH_SUFFIXES ceres ceres/lib
        NO_DEFAULT_PATH)

if(NOT TARGET ceres AND BUILD_SHARED_LIBS)
    find_package(Ceres
            HINTS  ${DIRECTORY_HINTS}
            PATHS ${CMAKE_INSTALL_PREFIX} $ENV{EBROOTCERES} $ENV{CERES_DIR} $ENV{ceres_DIR} $ENV{CONDA_PREFIX}
            PATH_SUFFIXES ceres ceres/lib
            NO_DEFAULT_PATH)
endif()


if(TARGET ceres)
    message(STATUS "ceres found")
else()


    message(STATUS "Ceres will be installed into ${CMAKE_INSTALL_PREFIX} on first build.")
    get_target_property(EIGEN3_INCLUDE_DIR Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
    list (GET EIGEN3_INCLUDE_DIR 0 EIGEN3_INCLUDE_DIR)
    if("${CMAKE_BUILD_TYPE}" MATCHES "Debug")
        set(CERES_LIBSUFFIX -debug)
    endif()

    get_target_property(GFLAGS_LIBRARIES    gflags      IMPORTED_LOCATION_${BUILD_TYPE})
    get_target_property(GFLAGS_INCLUDE_DIR  gflags      INTERFACE_INCLUDE_DIRECTORIES)
    if(NOT GFLAGS_LIBRARIES)
        get_target_property(GFLAGS_LIBRARIES    gflags      IMPORTED_LOCATION_NOCONFIG)
    endif()
    get_target_property(GLOG_LIBRARIES      glog::glog  IMPORTED_LOCATION_${BUILD_TYPE})
    get_target_property(GLOG_INCLUDE_DIR    glog::glog  INTERFACE_INCLUDE_DIRECTORIES)
    if(NOT GLOG_LIBRARIES)
        get_target_property(GLOG_LIBRARIES    glog::glog      IMPORTED_LOCATION_NOCONFIG)
    endif()
    unset(CERES_FLAGS CACHE)
    unset(CERES_FLAGS)
    if("${CMAKE_BUILD_TYPE}" MATCHES "Debug")
        set(CERES_FLAGS "${CERES_FLAGS} -O0 -g3 -fstack-protector -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -D_FORTIFY_SOURCE=2")
    elseif("${CMAKE_BUILD_TYPE}" MATCHES "RelWithDebInfo")
        set(CERES_FLAGS "${CERES_FLAGS} -O1 -g -fstack-protector -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -D_FORTIFY_SOURCE=2")
    else()
        set(CERES_FLAGS "${CERES_FLAGS} -O3 -g -DNDEBUG")
    endif()
    # It is VERY IMPORTANT to set this flag when using CERES
    # For some reason the Eigen destructor starts calling std::free instead of
    # the handmade allocator free in calls inside CERES. I do not know the reason.
    # Setting this flag here or globally helps.
    set(CERES_FLAGS "${CERES_FLAGS} -DEIGEN_MALLOC_ALREADY_ALIGNED=0")
#    set(CERES_FLAGS "${CERES_FLAGS} -I${GLOG_INCLUDE_DIR} -I${GFLAGS_INCLUDE_DIR} -DEIGEN_MALLOC_ALREADY_ALIGNED=0")
    string (REPLACE ";" " " CERES_FLAGS "${CERES_FLAGS}")
    list(APPEND CERES_CMAKE_OPTIONS  -DCMAKE_CXX_FLAGS:STRING=${CERES_FLAGS})
    list(APPEND CERES_CMAKE_OPTIONS  -DEIGEN_INCLUDE_DIR_HINTS:PATH=${CMAKE_INSTALL_PREFIX})
    list(APPEND CERES_CMAKE_OPTIONS  -DEigen3_DIR:PATH=${CMAKE_INSTALL_PREFIX}/Eigen3)
    list(APPEND CERES_CMAKE_OPTIONS  -Dglog_DIR:PATH=${CMAKE_INSTALL_PREFIX}/glog)
    list(APPEND CERES_CMAKE_OPTIONS  -Dgflags_DIR:PATH=${CMAKE_INSTALL_PREFIX}/gflags)
    list(APPEND CERES_CMAKE_OPTIONS  -DGFLAGS_INCLUDE_DIR_HINTS:PATH=${CMAKE_INSTALL_PREFIX})
    list(APPEND CERES_CMAKE_OPTIONS  -DGFLAGS_LIBRARY_DIR_HINTS:PATH=${CMAKE_INSTALL_PREFIX})
    list(APPEND CERES_CMAKE_OPTIONS  -DGLOG_INCLUDE_DIR_HINTS:PATH=${CMAKE_INSTALL_PREFIX})
    list(APPEND CERES_CMAKE_OPTIONS  -DGLOG_LIBRARY_DIR_HINTS:PATH=${CMAKE_INSTALL_PREFIX})

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

endif()


if(TARGET ceres)
    include(cmake-modules/filterTarget.cmake)
    string(TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE)
    get_target_property(ceres_imported_loc_buildtype  ceres IMPORTED_LOCATION_${BUILD_TYPE})
    get_target_property(ceres_imported_loc_noconfig  ceres IMPORTED_LOCATION_NOCONFIG)
    if(ceres_imported_loc_buildtype)
        target_link_libraries(gflags INTERFACE ${ceres_imported_loc_buildtype})
    elseif(ceres_imported_loc_noconfig)
        target_link_libraries(gflags INTERFACE ${ceres_imported_loc_noconfig})
    else()
        message(STATUS "Dependency ceres does not have IMPORTED_LOCATION_${BUILD_TYPE}/_NOCONFIG")
    endif()
    target_include_directories(ceres INTERFACE ${CERES_INCLUDE_DIR})
    remove_shared(ceres)
    remove_pthread(ceres)
    target_link_libraries(ceres INTERFACE Eigen3::Eigen glog::glog gflags)
endif()