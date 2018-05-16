


#message(STATUS "SEARCHING FOR ARPACKPP IN SYSTEM...")
#find_library(ARPACKPP_LIBRARIES
#        NAMES "arpack" "libarpack.a"
#        PATH_SUFFIXES "lib" "lib32" "lib64"
#        )

if(ARPACKPP_LIBRARIES)
    message(STATUS "Arpackpp FOUND IN SYSTEM: ${ARPACKPP_LIBRARIES}")
    message(STATUS "FOUND ARPACKPP:   ${ARPACKPP_LIBRARIES}")
    target_link_libraries(${PROJECT_NAME} ${ARPACKPP_LIBRARIES})
else()
    message(STATUS "Arpackpp will be installed into ${INSTALL_DIRECTORY}/arpackpp on first build.")
    include(ExternalProject)
    ExternalProject_Add(library_ARPACKPP
            GIT_REPOSITORY      https://github.com/m-reuter/arpackpp.git
            GIT_TAG             master
            PREFIX              "${INSTALL_DIRECTORY}/arpackpp"
#            UPDATE_COMMAND ""
#            INSTALL_COMMAND ""
            UPDATE_COMMAND ""
            TEST_COMMAND ""
            INSTALL_COMMAND ""
            CONFIGURE_COMMAND ""
            BUILD_COMMAND
            ${CMAKE_COMMAND} -E make_directory <INSTALL_DIR>/include && find <INSTALL_DIR>/include -maxdepth 1 -type l -delete &&
#            ln -sT <SOURCE_DIR>/include <INSTALL_DIR>/include/arpackpp
            ${CMAKE_COMMAND} -E create_symlink <SOURCE_DIR>/include <INSTALL_DIR>/include/arpackpp
            DEPENDS blas lapack arpack
    )

    ExternalProject_Get_Property(library_ARPACKPP INSTALL_DIR)
    add_library(arpackpp INTERFACE)
    set(ARPACKPP_INCLUDE_DIR ${INSTALL_DIR}/include)
    set_target_properties(arpackpp PROPERTIES
            INTERFACE_LINK_LIBRARIES arpack
            INTERFACE_LINK_LIBRARIES blas
            INTERFACE_LINK_LIBRARIES lapack
            INTERFACE_LINK_LIBRARIES EIGEN3
            )
    add_dependencies(arpackpp library_ARPACKPP)
    target_link_libraries(${PROJECT_NAME} PRIVATE arpackpp)
    target_include_directories(${PROJECT_NAME} PRIVATE ${ARPACKPP_INCLUDE_DIR})
endif()



