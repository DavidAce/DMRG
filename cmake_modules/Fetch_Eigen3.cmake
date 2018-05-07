


find_package(Eigen3 3.3)
if(EIGEN3_FOUND)
    message(STATUS "EIGEN FOUND IN SYSTEM: ${EIGEN3_INCLUDE_DIR}")
    add_library(EIGEN3 INTERFACE)
else()
    message(STATUS "Eigen3 will be installed into ${INSTALL_DIRECTORY}/eigen3 on first build.")

    include(ExternalProject)
    ExternalProject_Add(project_EIGEN3
            GIT_REPOSITORY https://github.com/eigenteam/eigen-git-mirror.git
            GIT_TAG 3.3.4
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


    ExternalProject_Get_Property(project_EIGEN3 INSTALL_DIR)
    add_library(EIGEN3 INTERFACE)
    set(EIGEN3_INCLUDE_DIR ${INSTALL_DIR}/include/eigen3)
#    target_compile_options(${PROJECT_NAME} PRIVATE -DEIGEN_USE_BLAS)
#    add_definitions(-DEIGEN_USE_LAPACK)

    add_dependencies(EIGEN3 project_EIGEN3)
endif()

target_link_libraries(${PROJECT_NAME} EIGEN3)
target_include_directories(${PROJECT_NAME} PRIVATE ${EIGEN3_INCLUDE_DIR})
