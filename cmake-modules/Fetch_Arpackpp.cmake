


if(NOT ARPACKPP_LIBRARIES)
message(STATUS "Searching for Arpack++ in system")
    find_library(ARPACKPP_LIBRARIES
        NAMES arpackpp arpack++ libarpack2++ libarpack++ libarpackpp
        PATH_SUFFIXES lib lib32 lib64 x86_64-linux-gnu lib/x86_64-linux-gnu
        PATHS /usr /usr/local
        )
    find_path(ARPACKPP_INCLUDE_DIR
        NAMES ardsnsym.h
        PATHS /usr/include/arpack++ /usr/include /usr/local/
        )
    if(ARPACKPP_LIBRARIES AND ARPACKPP_INCLUDE_DIR)
        message(STATUS "Searching for Arpack++ in system - Success: ${ARPACKPP_LIBRARIES}")
        message(STATUS "Note that old versions of Arpack++ (e.g. the default in Ubuntu Trusty 14.04 LTS) may fail to compile, they require '-fpermissive'.")
    else()
        message(STATUS "Searching for Arpack++ in system - failed")
    endif()
endif()


if (NOT ARPACKPP_LIBRARIES OR NOT ARPACKPP_INCLUDE_DIR)
    # Try finding arpack as module library
    set(ARPACKPP_LIBRARIES "")
    message(STATUS "Searching for Arpack++ in module")
    find_path(ARPACKPP_INCLUDE_DIR
            NAMES ardsnsym.h arpack++
            PATH_SUFFIXES include include/arpack++ include/arpackpp
            PATHS
            $ENV{EBROOTARPACKPLUSPLUS}
            $ENV{ARPACKPP_DIR}
            NO_DEFAULT_PATH
            )
    if(ARPACKPP_INCLUDE_DIR)
        get_filename_component(DIRNAME ${ARPACKPP_INCLUDE_DIR} NAME)
        if("${DIRNAME}" STREQUAL "arpack++")
            get_filename_component(ARPACKPP_INCLUDE_DIR ${ARPACKPP_INCLUDE_DIR} DIRECTORY)
            message(STATUS "Searching for Arpack++ in module - Success: ${ARPACKPP_INCLUDE_DIR}")
        elseif("${DIRNAME}" STREQUAL "arpackpp")
            add_definitions(-DARPACKPP_ALTDIR)
            get_filename_component(ARPACKPP_INCLUDE_DIR ${ARPACKPP_INCLUDE_DIR} DIRECTORY)
            message(STATUS "Searching for Arpack++ in module - Success: ${ARPACKPP_INCLUDE_DIR}")
            message(STATUS "Defined -DARPACKPP_ALTDIR")

        endif()
    else()
        message(STATUS "Searching for Arpack++ in module - failed")
    endif()
endif()


if (ARPACKPP_LIBRARIES OR ARPACKPP_INCLUDE_DIR)
    message(STATUS "Arpack++ library found: ${ARPACKPP_LIBRARIES}")
    message(STATUS "Arpack++ include found: ${ARPACKPP_INCLUDE_DIR}")
    add_library(arpack++ INTERFACE)
    target_link_libraries(arpack++ INTERFACE ${ARPACKPP_LIBRARIES} arpack)
    target_include_directories(arpack++ INTERFACE ${ARPACKPP_INCLUDE_DIR})
else()
    message(STATUS "Arpack++ will be installed into ${INSTALL_DIRECTORY}/arpackpp on first build.")
    include(ExternalProject)
    ExternalProject_Add(external_ARPACK++
            GIT_REPOSITORY      https://github.com/m-reuter/arpackpp.git
            GIT_TAG             master
            PREFIX      ${BUILD_DIRECTORY}/arpack++
            INSTALL_DIR ${INSTALL_DIRECTORY}/arpack++
#            CMAKE_ARGS
#            -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
#            -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
#            -DARPACK_LIB=${ARPACK_LIBRARIES_GENERATOR}
#            -LAPACK_LIBRARIES=${LAPACK_LIBRARIES_GENERATOR}
            UPDATE_COMMAND ""
            TEST_COMMAND ""
            INSTALL_COMMAND ""
            CONFIGURE_COMMAND ""
#            BUILD_COMMAND
#            ${CMAKE_COMMAND} -E make_directory <INSTALL_DIR>/include && find <INSTALL_DIR>/include -maxdepth 1 -type l -delete &&
#            ${CMAKE_COMMAND} -E create_symlink <SOURCE_DIR>/include <INSTALL_DIR>/include/arpack++
            BUILD_COMMAND
            ${CMAKE_COMMAND} -E copy_directory <SOURCE_DIR>/include <INSTALL_DIR>/include/arpack++
            DEPENDS arpack
            )

    ExternalProject_Get_Property(external_ARPACK++ INSTALL_DIR)
    add_library(arpack++ INTERFACE)
    add_dependencies(arpack++ external_ARPACK++ arpack)
    target_link_libraries(arpack++ INTERFACE arpack)
    target_include_directories(arpack++ INTERFACE ${INSTALL_DIR}/include)
endif()



