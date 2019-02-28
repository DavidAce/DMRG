
find_package(Eigen3 3.3.4)

if(NOT EIGEN3_FOUND)
    # Try finding arpack as module library
    message(STATUS "SEARCHING FOR eigen3 IN LOADED MODULES")
    find_package(Eigen3 3.3.4 PATHS "$ENV{EIGEN3_CMAKE_DIR}" NO_DEFAULT_PATH)
    if (NOT EIGEN3_FOUND)
    find_path(EIGEN3_INCLUDE_DIR
            NAMES Core
            PATHS "$ENV{EIGEN3_INCLUDE_DIR}/Eigen"
            )
    endif()
endif()


if(EIGEN3_FOUND OR EIGEN3_INCLUDE_DIR)
    message(STATUS "EIGEN FOUND IN SYSTEM: ${EIGEN3_INCLUDE_DIR}")
    add_library(eigen3 INTERFACE)
else()
    message(STATUS "Eigen3 will be installed into ${INSTALL_DIRECTORY}/eigen3 on first build.")

    include(ExternalProject)
    ExternalProject_Add(external_EIGEN3
            GIT_REPOSITORY https://github.com/eigenteam/eigen-git-mirror.git
            GIT_TAG 3.3.7
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


    ExternalProject_Get_Property(external_EIGEN3 INSTALL_DIR)
    add_library(eigen3 INTERFACE)
    set(EIGEN3_INCLUDE_DIR ${INSTALL_DIR}/include/eigen3)
    add_dependencies(eigen3 external_EIGEN3)
endif()

if(BLAS_LIBRARIES)
    set(EIGEN3_COMPILER_FLAGS  -Wno-parentheses) # -Wno-parentheses
    if("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU" )
        list(APPEND EIGEN3_COMPILER_FLAGS -Wno-unused-but-set-variable)
    endif()
    if(MKL_FOUND)
        list(APPEND EIGEN3_COMPILER_FLAGS -DEIGEN_USE_MKL_ALL)
        list(APPEND EIGEN3_COMPILER_FLAGS -DEIGEN_USE_LAPACKE_STRICT)
        list(APPEND EIGEN3_INCLUDE_DIR ${MKL_INCLUDE_DIR})
        message(STATUS "Eigen3 will use MKL")
    elseif (BLAS_FOUND)
        list(APPEND EIGEN3_COMPILER_FLAGS -DEIGEN_USE_BLAS)
        list(APPEND EIGEN3_COMPILER_FLAGS -DEIGEN_USE_LAPACKE)
        message(STATUS "Eigen3 will use BLAS and LAPACKE")
    endif()
endif()


set_target_properties(eigen3 PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES   "${EIGEN3_INCLUDE_DIR}"
        INTERFACE_COMPILE_OPTIONS       "${EIGEN3_COMPILER_FLAGS}"
        )
#target_link_libraries(${PROJECT_NAME} PRIVATE eigen3)
