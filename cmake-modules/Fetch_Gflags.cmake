message(STATUS "Searching for gflags ")
find_library(
        GFLAGS_LIBRARIES
        NAMES libgflags${CUSTOM_SUFFIX}
        PATH_SUFFIXES lib lib64
        PATHS  ${INSTALL_DIRECTORY}/gflags $ENV{EBROOTGFLAGS} $ENV{GFLAGS_DIR} $ENV{gflags_DIR}
        NO_DEFAULT_PATH)

find_path(
        GFLAGS_INCLUDE_DIR
        NAMES gflags/gflags.h
        PATH_SUFFIXES include gflags/include
        PATHS  ${INSTALL_DIRECTORY}/gflags $ENV{EBROOTGFLAGS} $ENV{GFLAGS_DIR} $ENV{gflags_DIR}
        NO_DEFAULT_PATH)

if(GFLAGS_LIBRARIES AND GFLAGS_INCLUDE_DIR)
    add_library(gflags INTERFACE)
    target_link_libraries(gflags INTERFACE ${GFLAGS_LIBRARIES})
    target_include_directories(gflags INTERFACE ${GFLAGS_INCLUDE_DIR})
    target_compile_definitions(gflags  INTERFACE "GFLAGS_IS_A_DLL=0")
    message(STATUS "Searching for gflags - Success: LIB: ${GFLAGS_LIBRARIES}")
    message(STATUS "Searching for gflags - Success: INC: ${GFLAGS_INCLUDE_DIR}")
    include(cmake-modules/PrintTargetProperties.cmake)
    print_target_properties(gflags)
else()
    message(STATUS "Searching for gflags - failed")
    message(STATUS "gflags will be installed into ${INSTALL_DIRECTORY}/gflags on first build.")

    include(ExternalProject)
    ExternalProject_Add(external_GFLAGS
            GIT_REPOSITORY https://github.com/gflags/gflags.git
            GIT_TAG v2.2.2
            GIT_PROGRESS 1
            PREFIX      ${BUILD_DIRECTORY}/gflags
            INSTALL_DIR ${INSTALL_DIRECTORY}/gflags
            UPDATE_COMMAND ""
            CMAKE_ARGS
            -DCMAKE_BUILD_TYPE=Release
            -DGFLAGS_BUILD_SHARED_LIBS:BOOL=ON
            -DGFLAGS_BUILD_STATIC_LIBS:BOOL=ON
            -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
            )

    ExternalProject_Get_Property(external_GFLAGS INSTALL_DIR)
    add_library(gflags INTERFACE)
    set(GFLAGS_INCLUDE_DIR ${INSTALL_DIR}/include)
    set(gflags_DIR ${INSTALL_DIR}/lib/cmake/gflags)
    add_dependencies(gflags external_GFLAGS)
    target_link_libraries(gflags INTERFACE ${INSTALL_DIR}/lib/libgflags${CUSTOM_SUFFIX})
    target_include_directories(gflags INTERFACE ${GFLAGS_INCLUDE_DIR})
    target_compile_definitions(gflags  INTERFACE "GFLAGS_IS_A_DLL=0")
endif()




