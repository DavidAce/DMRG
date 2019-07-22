



message(STATUS "SEARCHING FOR ARPACK++ IN SYSTEM...")
find_library(ARPACKPP_LIBRARIES
        NAMES arpackpp arpack++ libarpack2++ libarpack++ libarpackpp
        PATH_SUFFIXES lib lib32 lib64
        PATHS /usr/include/arpack++ /usr/include /usr/local/
        )
find_path(ARPACKPP_INCLUDE_DIR
        NAMES ardsnsym.h
        PATHS /usr/include/arpack++ /usr/include /usr/local/
        )

if (NOT ARPACKPP_LIBRARIES OR NOT ARPACKPP_INCLUDE_DIR)
    # Try finding arpack as module library
    set(ARPACKPP_LIBRARIES "")
    message(STATUS "SEARCHING FOR ARPACK IN LOADED MODULES")
    find_path(ARPACKPP_INCLUDE_DIR
            NAMES ardsnsym.h arpack++
            PATHS
            $ENV{EBROOTARPACKPLUSPLUS}/include
            $ENV{ARPACKPP_DIR}/include
            NO_DEFAULT_PATH
            )
endif()


message(STATUS "Note that old versions of Arpack++ (e.g. the default in Ubuntu Trusty 14.04 LTS) may fail to compile, requiring '-fpermissive'.")
if (ARPACKPP_LIBRARIES OR ARPACKPP_INCLUDE_DIR AND NOT "${OS_PROPERTIES}" MATCHES "trusty" )
    message(STATUS "Arpack++ library found in system: ${ARPACKPP_LIBRARIES}")
    message(STATUS "Arpack++ include found in system: ${ARPACKPP_INCLUDE_DIR}")
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



