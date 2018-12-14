message(STATUS "CppNumSolvers will be installed into ${INSTALL_DIRECTORY}/CppNumSolvers on first build.")

include(ExternalProject)
ExternalProject_Add(library_CppNumSolvers
        GIT_REPOSITORY https://github.com/PatWie/CppNumericalSolvers.git
        GIT_TAG master
        GIT_PROGRESS 1
        PREFIX "${INSTALL_DIRECTORY}/CppNumSolvers"
        UPDATE_COMMAND ""
        TEST_COMMAND ""
        INSTALL_COMMAND ""
        CONFIGURE_COMMAND ""
        BUILD_COMMAND
        find <INSTALL_DIR> -maxdepth 1 -type l -delete &&
        ln -s <SOURCE_DIR>/include <INSTALL_DIR>/include
        )


ExternalProject_Get_Property(library_CppNumSolvers INSTALL_DIR)
add_library(CppNumSolvers INTERFACE)
set(CppNumSolvers_INCLUDE_DIR ${INSTALL_DIR}/include)
add_dependencies(CppNumSolvers library_CppNumSolvers)

set_target_properties(CppNumSolvers PROPERTIES
        INTERFACE_INCLUDE_DIRECTORY     "${CppNumSolvers_INCLUDE_DIR}"
        )
target_link_libraries(${PROJECT_NAME} PRIVATE CppNumSolvers)
target_include_directories(${PROJECT_NAME} PRIVATE ${CppNumSolvers_INCLUDE_DIR})
