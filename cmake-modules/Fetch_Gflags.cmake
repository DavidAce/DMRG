message(STATUS "Searching for gflags ")
#find_library(
#        GFLAGS_LIBRARIES
#        NAMES libgflags${CUSTOM_SUFFIX}
#        PATH_SUFFIXES lib lib64
#        PATHS  ${EXTERNAL_INSTALL_DIR}/gflags $ENV{EBROOTGFLAGS} $ENV{GFLAGS_DIR} $ENV{gflags_DIR}
#        NO_DEFAULT_PATH)
#
#find_path(
#        GFLAGS_INCLUDE_DIR
#        NAMES gflags/gflags.h
#        PATH_SUFFIXES include gflags/include
#        PATHS  ${EXTERNAL_INSTALL_DIR}/gflags $ENV{EBROOTGFLAGS} $ENV{GFLAGS_DIR} $ENV{gflags_DIR}
#        NO_DEFAULT_PATH)

if("${CMAKE_BUILD_TYPE}" MATCHES "Debug")
    set(GFLAGS_LIBSUFFIX _debug)
endif()

if(GFLAGS_LIBRARIES AND GFLAGS_INCLUDE_DIR)
    add_library(gflags INTERFACE)
    get_filename_component(GFLAGS_LIBRARY_DIR ${GFLAGS_LIBRARIES} DIRECTORY)
    target_link_libraries(gflags INTERFACE ${GFLAGS_LIBRARIES})
    target_include_directories(gflags SYSTEM INTERFACE ${GFLAGS_INCLUDE_DIR})
    target_compile_definitions(gflags  INTERFACE "GFLAGS_IS_A_DLL=0")
    message(STATUS "Searching for gflags - Success: LIB: ${GFLAGS_LIBRARIES}")
    message(STATUS "Searching for gflags - Success: INC: ${GFLAGS_INCLUDE_DIR}")
    include(cmake-modules/PrintTargetProperties.cmake)
    print_target_properties(gflags)
else()
    message(STATUS "Searching for gflags - failed")
    message(STATUS "gflags will be installed into ${EXTERNAL_INSTALL_DIR}/gflags on first build.")
    unset(FLAGS CACHE)
    if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
        set(FLAGS "${FLAGS}  -stdlib=libstdc++ ${GCC_TOOLCHAIN}")
    endif()
    include(ExternalProject)
    ExternalProject_Add(external_GFLAGS
            GIT_REPOSITORY https://github.com/gflags/gflags.git
            GIT_TAG v2.2.2
            GIT_PROGRESS false
            GIT_SHALLOW true
            PREFIX      ${EXTERNAL_BUILD_DIR}/gflags
            INSTALL_DIR ${EXTERNAL_INSTALL_DIR}/gflags
            UPDATE_COMMAND ""
            CMAKE_ARGS
            -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
            -DCMAKE_CXX_FLAGS:STRING=${FLAGS}
            -DGFLAGS_BUILD_SHARED_LIBS:BOOL=ON
            -DGFLAGS_BUILD_STATIC_LIBS:BOOL=ON
            -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
            -DCMAKE_INSTALL_MESSAGE=NEVER #Avoid unnecessary output to console
            )

    ExternalProject_Get_Property(external_GFLAGS INSTALL_DIR)
    add_library(gflags INTERFACE)
    set(GFLAGS_INCLUDE_DIR ${INSTALL_DIR}/include)
    set(GFLAGS_LIBRARY_DIR ${INSTALL_DIR}/lib)
    set(GFLAGS_LIBRARIES   ${GFLAGS_LIBRARY_DIR}/libgflags${GFLAGS_LIBSUFFIX}${CUSTOM_SUFFIX})
    set(gflags_DIR ${INSTALL_DIR}/lib/cmake/gflags)
    add_dependencies(gflags external_GFLAGS)
    target_link_libraries(gflags INTERFACE ${GFLAGS_LIBRARIES})
    target_include_directories(gflags  SYSTEM INTERFACE ${GFLAGS_INCLUDE_DIR})
    target_compile_definitions(gflags  INTERFACE "GFLAGS_IS_A_DLL=0")
endif()




