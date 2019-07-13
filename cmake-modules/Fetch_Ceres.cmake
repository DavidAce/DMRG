
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
    list (GET EIGEN3_INCLUDE_DIR 0 EIGEN3_INCLUDE_DIR)

    get_target_property(LAPACK_LIBRARIES lapack INTERFACE_LINK_LIBRARIES)
    string (REPLACE ";" "$<SEMICOLON>" BLAS_LIBRARIES_GENERATOR     "${BLAS_LIBRARIES}")
    string (REPLACE ";" "$<SEMICOLON>" LAPACK_LIBRARIES_GENERATOR   "${LAPACK_LIBRARIES}")

    include(ExternalProject)
    set(FLAGS "-g -DEIGEN_MAX_STATIC_ALIGN_BYTES=0 -Wpedantic -Wextra -O3 -fno-omit-frame-pointer -D_GLIBCXX_DEBUG_PEDANTIC -D_GLIBCXX_DEBUG -D_FORTIFY_SOURCE=2")
    ExternalProject_Add(external_CERES
            GIT_REPOSITORY https://github.com/ceres-solver/ceres-solver.git
            GIT_TAG master
            GIT_PROGRESS 1
            PREFIX      ${BUILD_DIRECTORY}/ceres
            INSTALL_DIR ${INSTALL_DIRECTORY}/ceres
            UPDATE_COMMAND ""
            CMAKE_ARGS
            -DCMAKE_BUILD_TYPE=Release
            -DCMAKE_VERBOSE_MAKEFILE:BOOL=OFF
            -DCMAKE_CXX_FLAGS:STRING=${FLAGS}
            -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
            -DGFLAGS:BOOL=OFF
            -DSUITESPARSE:BOOL=OFF
            -DCXSPARSE:BOOL=OFF
            -DSCHUR_SPECIALIZATIONS:BOOL=OFF
            -DCUSTOM_BLAS:BOOL=OFF

            -DEIGEN_INCLUDE_DIR:PATH=${EIGEN3_INCLUDE_DIR}
            -DEIGEN_INCLUDE_DIR_HINTS:PATH=${EIGEN3_INCLUDE_DIR}
            -DEIGEN_PREFER_EXPORTED_EIGEN_CMAKE_CONFIGURATION:BOOL=OFF
            -DLAPACK:BOOL=OFF
            -DEIGENSPARSE:BOOL=OFF
            -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
            -DCMAKE_INSTALL_MESSAGE=NEVER #Avoid unnecessary output to console
            DEPENDS Eigen3::Eigen glog::glog
            )

    ExternalProject_Get_Property(external_CERES INSTALL_DIR)
    add_library(ceres INTERFACE)
    set(CERES_LIBRARY ${INSTALL_DIR}/lib/libceres${CUSTOM_SUFFIX})
    set(CERES_INCLUDE_DIR ${INSTALL_DIR}/include)
    add_dependencies(ceres external_CERES)
    add_dependencies(ceres glog::glog Eigen3::Eigen blas lapack lapacke)
    target_link_libraries(ceres INTERFACE Threads::Threads glog::glog Eigen3::Eigen ${CERES_LIBRARY} )
    target_include_directories(ceres INTERFACE ${CERES_INCLUDE_DIR})
endif()
