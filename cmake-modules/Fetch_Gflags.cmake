message(STATUS "Searching for gflags ")
find_package(gflags PATHS ${INSTALL_DIRECTORY}/gflags $ENV{EBROOTGFLAGS} $ENV{GFLAGS_DIR} $ENV{gflags_DIR} NO_DEFAULT_PATH)
include(cmake-modules/PrintTargetProperties.cmake)
print_target_properties(gflags)

if(TARGET gflags)
    # For some reason gflags imports libunwind.so even though we asked for static libraries.
    # In addition, libgflags.a hides in the property "LOCATION" instead of its rightful
    # place "INTERFACE_LINK_LIBRARIES".
#    add_library(gflags::gflags ALIAS gflags)

    get_target_property(GFLAGS_INCLUDE_DIR gflags INTERFACE_INCLUDE_DIRECTORIES)
    get_target_property(GFLAGS_LIBRARIES   gflags LOCATION)
    # In addition, the library may be .so
    get_filename_component(gflags_dir ${GFLAGS_LIBRARIES} DIRECTORY)
    get_filename_component(gflags_we  ${GFLAGS_LIBRARIES} NAME_WE)
    set(GFLAGS_LIBRARIES "${gflags_dir}/${gflags_we}${CUSTOM_SUFFIX}")
    set_target_properties(gflags PROPERTIES INTERFACE_LINK_LIBRARIES "${GFLAGS_LIBRARIES}")
    set_target_properties(gflags PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${GFLAGS_INCLUDE_DIR}")
    message(STATUS "Searching for gflags - Success: LIB: ${GFLAGS_LIBRARIES}")
    message(STATUS "Searching for gflags - Success: INC: ${GFLAGS_INCLUDE_DIR}")
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
    add_library(gflags::gflags ALIAS gflags)
    set(GFLAGS_INCLUDE_DIR ${INSTALL_DIR}/include)
    set(gflags_DIR ${INSTALL_DIR}/lib/cmake/gflags)
    add_dependencies(gflags external_GFLAGS)
    target_link_libraries(gflags INTERFACE ${INSTALL_DIR}/lib/libgflags${CUSTOM_SUFFIX})
    target_include_directories(gflags INTERFACE ${GFLAGS_INCLUDE_DIR})
endif()
