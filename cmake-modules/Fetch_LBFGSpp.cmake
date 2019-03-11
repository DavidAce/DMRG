message(STATUS "LBFGS++ will be installed into ${INSTALL_DIRECTORY}/LBFGS++ on first build.")

include(ExternalProject)
ExternalProject_Add(external_LBFGSpp
        GIT_REPOSITORY https://github.com/yixuan/LBFGSpp.git
        GIT_TAG master
        GIT_PROGRESS 1
        PREFIX      ${BUILD_DIRECTORY}/LBFGS++
        INSTALL_DIR ${INSTALL_DIRECTORY}/LBFGS++
        UPDATE_COMMAND ""
        TEST_COMMAND ""
        INSTALL_COMMAND ""
        CONFIGURE_COMMAND ""
        BUILD_COMMAND
        find <INSTALL_DIR> -maxdepth 1 -type l -delete &&
        ln -s <SOURCE_DIR>/include <INSTALL_DIR>/include
        )


ExternalProject_Get_Property(external_LBFGSpp INSTALL_DIR)
add_library(LBFGSpp INTERFACE)
add_dependencies(LBFGSpp external_LBFGSpp)
target_include_directories(LBFGSpp INTERFACE ${INSTALL_DIR}/include)

#target_link_libraries(${PROJECT_NAME} PRIVATE LBFGSpp)
