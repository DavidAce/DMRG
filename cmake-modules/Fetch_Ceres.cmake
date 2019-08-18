
#find_package(Ceres PATHS ${INSTALL_DIRECTORY}/ceres/lib/cmake/Ceres ${ceres_DIR})

if(Ceres_FOUND)
    message(STATUS "ceres FOUND IN SYSTEM: ${CERES_ROOT}")
    #    get_target_property(CERES_LIBRARY ceres  INTERFACE_LINK_LIBRARIES)
    #    set_target_properties(ceres PROPERTIES INTERFACE_LINK_LIBRARIES "")
    #    set_target_properties(ceres PROPERTIES INTERFACE_COMPILE_FEATURES "")
    #    target_link_libraries(ceres INTERFACE  blas lapack lapacke ${CERES_LIBRARY} Threads::Threads )
else()
    message(STATUS "ceres will be installed into ${INSTALL_DIRECTORY}/ceres on first build.")
    get_target_property(EIGEN3_INCLUDE_DIR Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
#    get_target_property(GFLAGS_INCLUDE_DIR   gflags    INTERFACE_INCLUDE_DIRECTORIES)
#    get_target_property(GFLAGS_LIBRARIES     gflags    INTERFACE_LINK_LIBRARIES)
#
#    get_target_property(GLOG_INCLUDE_DIR   glog    INTERFACE_INCLUDE_DIRECTORIES)
#    get_target_property(GLOG_LIBRARIES     glog    INTERFACE_LINK_LIBRARIES)

    list (GET EIGEN3_INCLUDE_DIR 0 EIGEN3_INCLUDE_DIR)

#    set(FLAGS "-DEIGEN_MAX_STATIC_ALIGN_BYTES=0 -DNDEBUG -O3 -fstack-protector  -g -fno-omit-frame-pointer -D_GLIBCXX_DEBUG_PEDANTIC -D_GLIBCXX_DEBUG -D_FORTIFY_SOURCE=2")
#    set(FLAGS "-DEIGEN_MAX_STATIC_ALIGN_BYTES=0 -DNDEBUG -O3 -fstack-protector  -g -fno-omit-frame-pointer -D_FORTIFY_SOURCE=2")
    unset(FLAGS CACHE)
    set(FLAGS "${COMMON_OPTIONS}")
    if(CMAKE_BUILD_TYPE MATCHES Debug)
        set(CERES_LIBSUFFIX -debug)
        set(FLAGS "${FLAGS} ${DEBUG_OPTIONS}")
    elseif(CMAKE_BUILD_TYPE MATCHES Release)
        set(FLAGS "${FLAGS} ${RELEASE_OPTIONS}")
    elseif(CMAKE_BUILD_TYPE MATCHES RelWithDebInfo)
        set(FLAGS "${FLAGS} ${RELWITHDEBINFO_OPTIONS}")
    endif()


    ################################
    ### Compiler-dependent flags ###
    ################################


    set(FLAGS "${FLAGS} -I${GLOG_INCLUDE_DIR} -L${GLOG_LIBRARY_DIR} -I${GFLAGS_INCLUDE_DIR}  -L${GFLAGS_LIBRARY_DIR}")
    string (REPLACE ";" " " FLAGS "${FLAGS}")



    message(STATUS "ceres flags: ${FLAGS}")
    message(STATUS "GFLAGS_LIBRARY_DIR  : ${GFLAGS_LIBRARY_DIR}")
    message(STATUS "GFLAGS_INCLUDE_DIR  : ${GFLAGS_INCLUDE_DIR}")
    message(STATUS "GLOG_LIBRARY_DIR    : ${GLOG_LIBRARY_DIR}")
    message(STATUS "GLOG_INCLUDE_DIR    : ${GLOG_INCLUDE_DIR}")

    include(ExternalProject)
    ExternalProject_Add(external_CERES
            GIT_REPOSITORY https://github.com/ceres-solver/ceres-solver.git
            GIT_TAG master
            GIT_PROGRESS false
            GIT_SHALLOW true
            PREFIX      ${BUILD_DIRECTORY}/ceres
            INSTALL_DIR ${INSTALL_DIRECTORY}/ceres
            UPDATE_COMMAND ""
            CMAKE_ARGS
            -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
            -DCMAKE_VERBOSE_MAKEFILE:BOOL=OFF
            -DBUILD_TESTING:BOOL=OFF
            -DBUILD_EXAMPLES:BOOL=OFF
            -DCMAKE_CXX_FLAGS:STRING=${FLAGS}
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
            -DGFLAGS_PREFER_EXPORTED_GFLAGS_CMAKE_CONFIGURATION:BOOL=FALSE
            -DGFLAGS_INCLUDE_DIR:PATH=${GFLAGS_INCLUDE_DIR}
            -DGFLAGS_LIBRARY:PATH=${GFLAGS_LIBRARIES}
            -DGLOG_PREFER_EXPORTED_GLOG_CMAKE_CONFIGURATION:BOOL=FALSE
            -DGLOG_INCLUDE_DIR:PATH=${GLOG_INCLUDE_DIR}
            -DGLOG_LIBRARY:PATH=${GLOG_LIBRARIES}
            -DLAPACK:BOOL=OFF
            -DEIGENSPARSE:BOOL=ON
            -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
            -DCMAKE_INSTALL_MESSAGE=NEVER #Avoid unnecessary output to console
            DEPENDS Eigen3::Eigen gflags glog
            )

    ExternalProject_Get_Property(external_CERES INSTALL_DIR)
    add_library(ceres INTERFACE)
    include(GNUInstallDirs)
    set(CERES_INCLUDE_DIR ${INSTALL_DIR}/${CMAKE_INSTALL_INCLUDEDIR})
    set(CERES_LIBRARY_DIR ${INSTALL_DIR}/${CMAKE_INSTALL_LIBDIR})
    set(CERES_LIBRARY     ${CERES_LIBRARY_DIR}/libceres${CERES_LIBSUFFIX}${CUSTOM_SUFFIX})
    add_dependencies(ceres external_CERES)
    add_dependencies(ceres glog gflags Eigen3::Eigen blas lapack lapacke)
    target_link_libraries(ceres INTERFACE  ${CERES_LIBRARY} glog gflags Eigen3::Eigen blas lapack lapacke Threads::Threads )
    target_include_directories(ceres INTERFACE ${CERES_INCLUDE_DIR})
endif()
