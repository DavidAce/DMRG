include(cmake-modules/filterTarget.cmake)
string(TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE)
find_package(Ceres PATHS ${CMAKE_INSTALL_PREFIX}/ceres)

if(TARGET ceres)
    message(STATUS "ceres found")
    get_target_property(cereslib  ceres IMPORTED_LOCATION_${BUILD_TYPE})
    target_link_libraries(ceres INTERFACE ${cereslib})
    target_include_directories(ceres INTERFACE ${CERES_INCLUDE_DIR})
    remove_shared(ceres)
    remove_pthread(ceres)
    target_link_libraries(ceres INTERFACE Eigen3::Eigen glog::glog gflags Threads::Threads)


else()

    message(STATUS "Ceres will be installed into ${CMAKE_INSTALL_PREFIX} on first build.")
    get_target_property(EIGEN3_INCLUDE_DIR Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
    list (GET EIGEN3_INCLUDE_DIR 0 EIGEN3_INCLUDE_DIR)
    if("${CMAKE_BUILD_TYPE}" MATCHES "Debug")
        set(CERES_LIBSUFFIX -debug)
    endif()

    get_target_property(GFLAGS_LIBRARIES    gflags      INTERFACE_LINK_LIBRARIES)
    get_target_property(GFLAGS_INCLUDE_DIR  gflags      INTERFACE_INCLUDE_DIRECTORIES)
    get_target_property(GLOG_LIBRARIES      glog::glog  INTERFACE_LINK_LIBRARIES)
    get_target_property(GLOG_INCLUDE_DIR    glog::glog  INTERFACE_INCLUDE_DIRECTORIES)
    unset(CERES_FLAGS CACHE)
    unset(CERES_FLAGS)
    if("${CMAKE_BUILD_TYPE}" MATCHES "Debug")
        set(CERES_FLAGS "${CERES_FLAGS} -O0 -g3 -fstack-protector -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -D_FORTIFY_SOURCE=2")
    elseif("${CMAKE_BUILD_TYPE}" MATCHES "RelWithDebInfo")
        set(CERES_FLAGS "${CERES_FLAGS} -O1 -g -fstack-protector -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -D_FORTIFY_SOURCE=2")
    else()
        set(CERES_FLAGS "${CERES_FLAGS} -O3 -DNDEBUG")
    endif()

    set(CERES_FLAGS "${CERES_FLAGS} -I${GLOG_INCLUDE_DIR} -I${GFLAGS_INCLUDE_DIR} ")
    string (REPLACE ";" " " CERES_FLAGS "${CERES_FLAGS}")
    list(APPEND CERES_CMAKE_OPTIONS  -DCERES_FLAGS:STRING=${CERES_FLAGS})
    list(APPEND CERES_CMAKE_OPTIONS  -DCUSTOM_SUFFIX:STRING=${CUSTOM_SUFFIX})
    list(APPEND CERES_CMAKE_OPTIONS  -DEigen3_DIR:PATH=${CMAKE_INSTALL_PREFIX}/Eigen3)
    list(APPEND CERES_CMAKE_OPTIONS  -DEIGEN3_INCLUDE_DIR:PATH=${EIGEN3_INCLUDE_DIR})
    list(APPEND CERES_CMAKE_OPTIONS  -Dgflags_DIR:PATH=${CMAKE_INSTALL_PREFIX}/gflags)
    list(APPEND CERES_CMAKE_OPTIONS  -DGFLAGS_INCLUDE_DIR:PATH=${GFLAGS_INCLUDE_DIR})
    list(APPEND CERES_CMAKE_OPTIONS  -DGFLAGS_LIBRARIES:PATH=${GFLAGS_LIBRARIES})
    list(APPEND CERES_CMAKE_OPTIONS  -Dglog_DIR:PATH=${CMAKE_INSTALL_PREFIX}/glog)
    list(APPEND CERES_CMAKE_OPTIONS  -DGLOG_INCLUDE_DIR:PATH=${GLOG_INCLUDE_DIR})
    list(APPEND CERES_CMAKE_OPTIONS  -DGLOG_LIBRARIES:PATH=${GLOG_LIBRARIES})

    message(STATUS "ceres flags  : ${CERES_FLAGS}")
    message(STATUS "ceres options: ${CERES_CMAKE_OPTIONS}")
    include(${PROJECT_SOURCE_DIR}/cmake-modules/BuildDependency.cmake)
    build_dependency(ceres "${CERES_CMAKE_OPTIONS}")
    find_package(Ceres HINTS ${CMAKE_INSTALL_PREFIX}/ceres)
    if(TARGET ceres)
        message(STATUS "ceres installed successfully")
        get_target_property(cereslib  ceres IMPORTED_LOCATION_${BUILD_TYPE})
        target_link_libraries(ceres INTERFACE ${cereslib})
        target_link_libraries(ceres INTERFACE Eigen3::Eigen glog::glog gflags Threads::Threads)
        remove_shared(ceres)
        remove_pthread(ceres)
        target_include_directories(ceres INTERFACE ${CERES_INCLUDE_DIR})
    else()
        message(STATUS "config_result: ${config_result}")
        message(STATUS "build_result: ${build_result}")
        message(FATAL_ERROR "ceres could not be downloaded.")
    endif()

endif()
