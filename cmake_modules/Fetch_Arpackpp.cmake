



message(STATUS "SEARCHING FOR ARPACKPP IN SYSTEM...")
find_library(ARPACKPP_LIBRARIES
        NAMES arpackpp arpack++ libarpack2++ libarpack++ libarpackpp
        PATH_SUFFIXES lib lib32 lib64
        )
find_path(ARPACKPP_INCLUDE_DIR
        NAMES ardsnsym.h
        PATHS /usr/include/arpack++ /usr/include /usr/local/
        )
message(STATUS "Arpack++ library found in system: ${ARPACKPP_LIBRARIES}")
message(STATUS "Arpack++ include found in system: ${ARPACKPP_INCLUDE_DIR}")
message(WARNING
"Due to a bug, older versions of Arpack++ require -fpermissive flag to compile. \
Therefore, local installations of Arpack++ will be ignored in favor of the most resent version.")




if (false)
#if(ARPACKPP_LIBRARIES AND ARPACKPP_INCLUDE_DIR)

    add_library(arpack++ UNKNOWN IMPORTED)
    set_target_properties(arpack++ PROPERTIES
            IMPORTED_LOCATION "${ARPACKPP_LIBRARIES}"
            INCLUDE_DIRECTORIES "${ARPACKPP_INCLUDE_DIR}")
    target_link_libraries(${PROJECT_NAME} PRIVATE arpackpp)

else()
    message(STATUS "Arpack++ will be installed into ${INSTALL_DIRECTORY}/arpackpp on first build.")
    enable_language(Fortran)
    include(ExternalProject)
    ExternalProject_Add(library_ARPACK++
            GIT_REPOSITORY      https://github.com/m-reuter/arpackpp.git
            GIT_TAG             master
            PREFIX              "${INSTALL_DIRECTORY}/arpack++"
            UPDATE_COMMAND ""
            TEST_COMMAND ""
            INSTALL_COMMAND ""
            CONFIGURE_COMMAND ""
            BUILD_COMMAND
            ${CMAKE_COMMAND} -E make_directory <INSTALL_DIR>/include && find <INSTALL_DIR>/include -maxdepth 1 -type l -delete &&
            ${CMAKE_COMMAND} -E create_symlink <SOURCE_DIR>/include <INSTALL_DIR>/include/arpack++
            DEPENDS blas lapack arpack
    )

    ExternalProject_Get_Property(library_ARPACK++ INSTALL_DIR)
    add_library(arpack++ INTERFACE)
    set(ARPACKPP_INCLUDE_DIR ${INSTALL_DIR}/include)
    set_target_properties(arpack++ PROPERTIES
            INTERFACE_LINK_LIBRARIES arpack
            INTERFACE_LINK_LIBRARIES blas
            INTERFACE_LINK_LIBRARIES lapack
            INTERFACE_INCLUDE_DIRECTORIES ${ARPACKPP_INCLUDE_DIR}
            )
    add_dependencies(arpack++ library_ARPACK++)
    target_link_libraries(${PROJECT_NAME} PRIVATE arpack++)
endif()



