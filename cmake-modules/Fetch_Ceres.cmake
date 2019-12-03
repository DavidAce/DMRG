
#find_package(Ceres PATHS ${EXTERNAL_INSTALL_DIR}/ceres/lib/cmake/Ceres ${ceres_DIR})

#if(Ceres_FOUND)
#    message(STATUS "ceres FOUND IN SYSTEM: ${CERES_ROOT}")
#    #    get_target_property(CERES_LIBRARY ceres  INTERFACE_LINK_LIBRARIES)
#    #    set_target_properties(ceres PROPERTIES INTERFACE_LINK_LIBRARIES "")
#    #    set_target_properties(ceres PROPERTIES INTERFACE_COMPILE_FEATURES "")
#    #    target_link_libraries(ceres INTERFACE  blas lapack lapacke ${CERES_LIBRARY} Threads::Threads )
#else()


    message(STATUS "Ceres will be installed into ${EXTERNAL_INSTALL_DIR}/ceres on first build.")
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
    message(STATUS "ceres flags: ${CERES_FLAGS}")

    include(ExternalProject)
    ExternalProject_Add(external_CERES
            GIT_REPOSITORY https://github.com/ceres-solver/ceres-solver.git
            GIT_TAG "486d81812e83f9274daaca356153d302c5ba58e0"
            GIT_PROGRESS false
            PREFIX      ${EXTERNAL_BUILD_DIR}/ceres
            INSTALL_DIR ${EXTERNAL_INSTALL_DIR}/ceres
#            TEST_COMMAND ${CMAKE_MAKE_PROGRAM} -j test
            CMAKE_ARGS
            -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
            -DCMAKE_CXX_FLAGS:STRING=${CERES_FLAGS}
            -DCMAKE_VERBOSE_MAKEFILE:BOOL=OFF
            -DBUILD_TESTING:BOOL=OFF
            -DBUILD_EXAMPLES:BOOL=OFF
            -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
            -DCMAKE_FIND_LIBRARY_SUFFIXES:STRING=${CUSTOM_SUFFIX}
            -DGFLAGS:BOOL=ON
            -DSUITESPARSE:BOOL=OFF
            -DCXSPARSE:BOOL=OFF
            -DSCHUR_SPECIALIZATIONS:BOOL=OFF
            -DCUSTOM_BLAS:BOOL=OFF
            -DEIGEN_INCLUDE_DIR:PATH=${EIGEN3_INCLUDE_DIR}
            -DEIGEN_INCLUDE_DIR_HINTS:PATH=${EIGEN3_INCLUDE_DIR}
            -DEIGEN_PREFER_EXPORTED_EIGEN_CMAKE_CONFIGURATION:BOOL=FALSE
            -DGFLAGS_PREFER_EXPORTED_GFLAGS_CMAKE_CONFIGURATION:BOOL=TRUE
            -DGFLAGS_INCLUDE_DIR:PATH=${GFLAGS_INCLUDE_DIR}
            -DGFLAGS_LIBRARY:PATH=${GFLAGS_LIBRARIES}
            -DGLOG_PREFER_EXPORTED_GLOG_CMAKE_CONFIGURATION:BOOL=TRUE
            -DGLOG_INCLUDE_DIR:PATH=${GLOG_INCLUDE_DIR}
            -DGLOG_LIBRARY:PATH=${GLOG_LIBRARIES}
            -DLAPACK:BOOL=OFF
            -DEIGENSPARSE:BOOL=ON
            -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
            -DCMAKE_INSTALL_MESSAGE=NEVER #Avoid unnecessary output to console
            DEPENDS Eigen3::Eigen gflags glog::glog
            )

    ExternalProject_Get_Property(external_CERES INSTALL_DIR)
    add_library(ceres INTERFACE)
    include(GNUInstallDirs)
    set(CERES_INCLUDE_DIR ${INSTALL_DIR}/${CMAKE_INSTALL_INCLUDEDIR})
    set(CERES_LIBRARY_DIR ${INSTALL_DIR}/${CMAKE_INSTALL_LIBDIR})
    set(CERES_LIBRARY     ${CERES_LIBRARY_DIR}/libceres${CERES_LIBSUFFIX}${CUSTOM_SUFFIX})
#    set(CERES_LIBRARY     ${CERES_LIBRARY_DIR}/libceres${CUSTOM_SUFFIX})
    add_dependencies(ceres external_CERES)
    add_dependencies(ceres glog::glog gflags Eigen3::Eigen blas lapack lapacke)
    target_link_libraries(ceres INTERFACE  ${CERES_LIBRARY} glog::glog gflags Eigen3::Eigen blas lapack lapacke Threads::Threads )
    target_include_directories(ceres SYSTEM INTERFACE ${CERES_INCLUDE_DIR})
#endif()
