
# Try finding in system. While avoiding Python anaconda's version.
# Usually there is an executable present in system, "h5c++"

find_file(HDF5_C_COMPILER_EXECUTABLE NAMES h5cc PATHS /usr/bin /usr/local/bin NO_DEFAULT_PATH)
find_file(HDF5_CXX_COMPILER_EXECUTABLE NAMES h5c++ h5cc h5fc h5pfc PATHS /usr/bin /usr/local/bin NO_DEFAULT_PATH)
find_file(HDF5_Fortran_COMPILER_EXECUTABLE NAMES h5fc h5pfc PATHS /usr/bin /usr/local/bin NO_DEFAULT_PATH)

message("HDF5 C compiler wrapper             : ${HDF5_C_COMPILER_EXECUTABLE}")
message("HDF5 CXX compiler wrapper           : ${HDF5_CXX_COMPILER_EXECUTABLE}")
message("HDF5 Fortran compiler wrapper       : ${HDF5_Fortran_COMPILER_EXECUTABLE}")
set(HDF5_USE_STATIC_LIBRARIES ON)
set(HDF5_FIND_DEBUG OFF)
find_package(HDF5 COMPONENTS C CXX HL)
if(HDF5_LIBRARIES MATCHES "anaconda")
    message("Found anaconda version. Ignoring...")
    set(HDF5_ANACONDA ON)
    unset(HDF5_FOUND)
    unset(HDF5_LIBRARIES)
    unset(HDF5_INCLUDE_DIR)
    unset(HDF5_DEFINITIONS)
endif()



if(HDF5_FOUND AND HDF5_LIBRARIES AND HDF5_CXX_LIBRARIES AND HDF5_HL_LIBRARIES AND HDF5_CXX_HL_LIBRARIES AND NOT HDF5_ANACONDA)
    message(STATUS "HDF5 FOUND IN SYSTEM: ${HDF5_LIBRARIES}")
    message(STATUS "FOUND HDF5")
    message(STATUS "   In path: ${HDF5_INCLUDE_DIR}")
    message(STATUS "   HDF5 DEFINITIONS         : ${HDF5_DEFINITIONS}")
    message(STATUS "   HDF5_LIBRARIES           : ${HDF5_CXX_LIBRARY_hdf5}")
    message(STATUS "   HDF5_HL_LIBRARIES        : ${HDF5_CXX_LIBRARY_hdf5_hl}")
    message(STATUS "   HDF5_CXX_LIBRARIES       : ${HDF5_CXX_LIBRARY_hdf5_cpp}")
    message(STATUS "   HDF5_HL_CXX_LIBRARIES    : ${HDF5_CXX_LIBRARY_hdf5_hl_cpp}")
    message(STATUS "   HDF5 LDFLAGS             : ${HDF5_LDFLAGS}")
    set(HDF5_LIBRARIES_ORIGINAL         ${HDF5_LIBRARIES})

    set(HDF5_LIBRARIES         ${HDF5_CXX_LIBRARY_hdf5})
    set(HDF5_HL_LIBRARIES      ${HDF5_CXX_LIBRARY_hdf5_hl})
    set(HDF5_CXX_LIBRARIES     ${HDF5_CXX_LIBRARY_hdf5_cpp})
    set(HDF5_HL_CXX_LIBRARIES  ${HDF5_CXX_LIBRARY_hdf5_hl_cpp})


    set(HDF5_LINKER_FLAGS ${HDF5_CXX_LIBRARY_NAMES})
    list(REMOVE_ITEM HDF5_LINKER_FLAGS hdf5_cpp hdf5)

    add_definitions(${HDF5_DEFINITIONS})
    # Add convenience libraries
    add_library(hdf5::hdf5           STATIC IMPORTED)
    add_library(hdf5::hdf5_hl        STATIC IMPORTED)
    add_library(hdf5::hdf5_cpp       STATIC IMPORTED)
    add_library(hdf5::hdf5_hl_cpp    STATIC IMPORTED)
    get_cmake_property(_variableNames VARIABLES)
    foreach (_variableName ${_variableNames})
        if("${_variableName}" MATCHES "HDF5" OR "${_variableName}" MATCHES "hdf5")
            message(STATUS "${_variableName}=${${_variableName}}")
        endif()
    endforeach()

else()
    message(STATUS "HDF5 will be installed into ${INSTALL_DIRECTORY}/hdf5 on first build.")

    include(ExternalProject)
    ExternalProject_Add(library_HDF5
            URL     https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.1/src/hdf5-1.10.1.tar.bz2 # version 1.10.2
            PREFIX              "${INSTALL_DIRECTORY}/hdf5"
            UPDATE_DISCONNECTED 1
            TEST_COMMAND ""
            CMAKE_ARGS
            -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
            -DBUILD_STATIC_LIBS:BOOL=ON
            -DBUILD_SHARED_LIBS:BOOL=OFF
            -DCMAKE_ANSI_CFLAGS:STRING=-fPIC
            -DCMAKE_BUILD_WITH_INSTALL_RPATH:BOOL=OFF
            -DHDF5_ENABLE_PARALLEL=OFF
            -DALLOW_UNSUPPORTED=ON
            -DBUILD_TESTING:BOOL=OFF
            -DHDF5_BUILD_TOOLS:BOOL=ON
            -DHDF5_BUILD_EXAMPLES:BOOL=OFF
            -DHDF5_BUILD_FORTRAN:BOOL=OFF
            -DHDF5_BUILD_JAVA:BOOL=OFF
            -DCMAKE_INSTALL_MESSAGE=NEVER #Avoid unnecessary output to console
            -DCMAKE_C_FLAGS=-w
            )

    ExternalProject_Get_Property(library_HDF5 INSTALL_DIR)
    add_library(hdf5::hdf5           STATIC IMPORTED)
    add_library(hdf5::hdf5_hl        STATIC IMPORTED)
    add_library(hdf5::hdf5_cpp       STATIC IMPORTED)
    add_library(hdf5::hdf5_hl_cpp    STATIC IMPORTED)

    add_dependencies(hdf5::hdf5          library_HDF5)
    add_dependencies(hdf5::hdf5_hl       library_HDF5)
    add_dependencies(hdf5::hdf5_cpp      library_HDF5)
    add_dependencies(hdf5::hdf5_hl_cpp   library_HDF5)
    set(HDF5_LIBRARIES        ${INSTALL_DIR}/lib/libhdf5-static${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(HDF5_HL_LIBRARIES     ${INSTALL_DIR}/lib/libhdf5_hl-static${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(HDF5_CXX_LIBRARIES    ${INSTALL_DIR}/lib/libhdf5_cpp-static${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(HDF5_HL_CXX_LIBRARIES ${INSTALL_DIR}/lib/libhdf5_hl_cpp-static${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(HDF5_INCLUDE_DIR      ${INSTALL_DIR}/include)
endif()

set_target_properties(hdf5::hdf5 PROPERTIES
        IMPORTED_LOCATION ${HDF5_LIBRARIES}
        INTERFACE_LINK_LIBRARIES "${HDF5_LDFLAGS};${HDF5_LINKER_FLAGS}"
        INCLUDE_DIRECTORIES ${HDF5_INCLUDE_DIR})

set_target_properties(hdf5::hdf5_hl PROPERTIES
        IMPORTED_LOCATION ${HDF5_HL_LIBRARIES}
        INTERFACE_LINK_LIBRARIES hdf5::hdf5
        INCLUDE_DIRECTORIES ${HDF5_INCLUDE_DIR})

set_target_properties(hdf5::hdf5_cpp PROPERTIES
        IMPORTED_LOCATION ${HDF5_CXX_LIBRARIES}
        INTERFACE_LINK_LIBRARIES "hdf5::hdf5;hdf5::hdf5_hl"
        INCLUDE_DIRECTORIES ${HDF5_INCLUDE_DIR})

set_target_properties(hdf5::hdf5_hl_cpp PROPERTIES
        IMPORTED_LOCATION ${HDF5_HL_CXX_LIBRARIES}
        INTERFACE_LINK_LIBRARIES "hdf5::hdf5;hdf5::hdf5_hl;hdf5::hdf5_cpp"
        INCLUDE_DIRECTORIES ${HDF5_INCLUDE_DIR})


target_link_libraries(${PROJECT_NAME} PRIVATE hdf5::hdf5)
target_link_libraries(${PROJECT_NAME} PRIVATE hdf5::hdf5_hl)
target_link_libraries(${PROJECT_NAME} PRIVATE hdf5::hdf5_cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE hdf5::hdf5_hl_cpp)
target_include_directories(${PROJECT_NAME} PRIVATE ${HDF5_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} PRIVATE ${HDF5_LDFLAGS} ${HDF5_LINKER_FLAGS})
