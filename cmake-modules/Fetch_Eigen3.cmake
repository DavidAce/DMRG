
find_package(Eigen3 3.3.4)

if(NOT EIGEN3_FOUND)
    # Try finding arpack as module library
    message(STATUS "SEARCHING FOR EIGEN3 IN LOADED MODULES")
    find_package(Eigen3 3.3.4 PATHS "$ENV{EIGEN3_CMAKE_DIR}")
    if (NOT EIGEN3_FOUND)
    find_path(EIGEN3_INCLUDE_DIR
            NAMES Core
            PATHS "$ENV{EIGEN3_INCLUDE_DIR}/Eigen"
            )
    endif()
endif()


if(EIGEN3_FOUND OR EIGEN3_INCLUDE_DIR)
    message(STATUS "EIGEN FOUND IN SYSTEM: ${EIGEN3_INCLUDE_DIR}")
    add_library(EIGEN3 INTERFACE)
else()
    message(STATUS "Eigen3 will be installed into ${INSTALL_DIRECTORY}/eigen3 on first build.")

    include(ExternalProject)
    ExternalProject_Add(library_EIGEN3
            GIT_REPOSITORY https://github.com/eigenteam/eigen-git-mirror.git
            GIT_TAG 3.3.5
            GIT_PROGRESS 1
            PREFIX "${INSTALL_DIRECTORY}/eigen3"
            UPDATE_COMMAND ""
            TEST_COMMAND ""
            INSTALL_COMMAND ""
            CONFIGURE_COMMAND ""
            BUILD_COMMAND
                ${CMAKE_COMMAND} -E make_directory <INSTALL_DIR>/include/eigen3 && find <INSTALL_DIR>/include/eigen3 -maxdepth 1 -type l -delete &&
                ln -s <SOURCE_DIR>/Eigen/ <SOURCE_DIR>/unsupported/ <INSTALL_DIR>/include/eigen3/
        )


    ExternalProject_Get_Property(library_EIGEN3 INSTALL_DIR)
    add_library(EIGEN3 INTERFACE)
    set(EIGEN3_INCLUDE_DIR ${INSTALL_DIR}/include/eigen3)
    add_dependencies(EIGEN3 library_EIGEN3)
endif()

if(BLAS_LIBRARIES)
    set(EIGEN3_COMPILER_FLAGS  ) # -Wno-parentheses
    if("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU" )
        list(APPEND EIGEN3_COMPILER_FLAGS -Wno-unused-but-set-variable)
    endif()
    if(MKL_FOUND)
        list(APPEND EIGEN3_COMPILER_FLAGS -DEIGEN_USE_MKL_ALL)
        message(STATUS "Eigen3 will use MKL")
    else()
        list(APPEND EIGEN3_COMPILER_FLAGS -DEIGEN_USE_BLAS)
        message(STATUS "Eigen3 will use BLAS")
    endif()
endif()

#if (OpenMP_FOUND)
#    list(APPEND EIGEN_3_COMPILER_FLAGS ${OpenMP_CXX_FLAGS})
#endif()


set_target_properties(EIGEN3 PROPERTIES
        INTERFACE_INCLUDE_DIRECTORY     "${EIGEN3_INCLUDE_DIR}"
        INTERFACE_COMPILE_OPTIONS       "${EIGEN3_COMPILER_FLAGS}"
        )
target_link_libraries(${PROJECT_NAME} PRIVATE EIGEN3)
# Add SYSTEM flag to suppress warnings
target_include_directories(${PROJECT_NAME} SYSTEM PRIVATE ${EIGEN3_INCLUDE_DIR})
target_compile_options(${PROJECT_NAME} PRIVATE ${EIGEN3_COMPILER_FLAGS})
