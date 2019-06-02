message(STATUS "CppNumSolvers will be installed into ${INSTALL_DIRECTORY}/CppNumSolvers on first build.")

include(ExternalProject)
ExternalProject_Add(external_CppNumSolvers
        GIT_REPOSITORY https://github.com/PatWie/CppNumericalSolvers.git
        GIT_TAG master
        GIT_PROGRESS 1
        PREFIX      ${BUILD_DIRECTORY}/CppNumSolvers
        INSTALL_DIR ${INSTALL_DIRECTORY}/CppNumSolvers
        UPDATE_COMMAND ""
        TEST_COMMAND ""
        INSTALL_COMMAND ""
        CONFIGURE_COMMAND ""
        BUILD_COMMAND
        ${CMAKE_COMMAND} -E copy_directory <SOURCE_DIR>/include <INSTALL_DIR>/include
        DEPENDS Eigen3::Eigen
        )


ExternalProject_Get_Property(external_CppNumSolvers INSTALL_DIR)
add_library(CppNumSolvers INTERFACE)
add_dependencies(CppNumSolvers external_CppNumSolvers)
target_link_libraries(CppNumSolvers INTERFACE Eigen3::Eigen)
target_include_directories(CppNumSolvers INTERFACE ${INSTALL_DIR}/include)
set(CppNumSolvers_INCLUDE_DIR ${INSTALL_DIR}/include)
