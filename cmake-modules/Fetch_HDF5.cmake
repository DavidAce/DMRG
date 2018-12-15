

if (EXISTS "$ENV{HDF5_DIR}")
    # Try finding in system. While avoiding Python anaconda's version.
    # Usually there is an executable present in system, "h5c++"
    cmake_policy(SET CMP0074 NEW)
    set(HDF5_ROOT "$ENV{HDF5_DIR}")
    set(HDF5_USE_STATIC_LIBRARIES ${STATIC_BUILD})
    set(HDF5_FIND_DEBUG OFF)
    find_package(HDF5 COMPONENTS C CXX HL)
    set(HDF5_LINKER_FLAGS      -Wl,--no-as-needed -ldl -lm -lz -lsz -lpthread )
#    if (HDF5_IS_PARALLEL)
#        list(APPEND HDF5_LINKER_FLAGS $ENV{MPI_LIB}/libmpi${CMAKE_STATIC_LIBRARY_SUFFIX})
#        list(APPEND HDF5_INCLUDE_DIR  $ENV{MPI_INCLUDE})
#    endif()
else()
    # Try finding in system. While avoiding Python anaconda's version.
    # Usually there is an executable present in system, "h5c++"

    find_file(HDF5_C_COMPILER_EXECUTABLE NAMES h5cc PATHS /usr/bin /usr/local/bin NO_DEFAULT_PATH)
    find_file(HDF5_CXX_COMPILER_EXECUTABLE NAMES h5c++ h5cc h5fc h5pfc PATHS /usr/bin /usr/local/bin NO_DEFAULT_PATH)
    find_file(HDF5_Fortran_COMPILER_EXECUTABLE NAMES h5fc h5pfc PATHS /usr/bin /usr/local/bin NO_DEFAULT_PATH)
#    message("HDF5 C compiler wrapper             : ${HDF5_C_COMPILER_EXECUTABLE}")
    message("HDF5 CXX compiler wrapper           : ${HDF5_CXX_COMPILER_EXECUTABLE}")

    set(HDF5_USE_STATIC_LIBRARIES ${STATIC_BUILD})
    set(HDF5_FIND_DEBUG OFF)
    find_package(HDF5 COMPONENTS C CXX HL)
    set(HDF5_C_LIBRARY         ${HDF5_CXX_LIBRARY_hdf5})
    set(HDF5_C_HL_LIBRARY      ${HDF5_CXX_LIBRARY_hdf5_hl})
    set(HDF5_CXX_LIBRARY       ${HDF5_CXX_LIBRARY_hdf5_cpp})
    set(HDF5_CXX_HL_LIBRARY    ${HDF5_CXX_LIBRARY_hdf5_hl_cpp})
    set(HDF5_LINKER_FLAGS       -Wl,-Bsymbolic-functions -Wl,-z,relro -lpthread -lz -ldl -lm)


endif()

if(HDF5_LIBRARIES MATCHES "anaconda")
    message("Found anaconda version. Ignoring...")
    set(HDF5_ANACONDA ON)
endif()

if(HDF5_FOUND AND NOT HDF5_ANACONDA)
    # To print all variables, use the code below:
    #
    get_cmake_property(_variableNames VARIABLES)
    foreach (_variableName ${_variableNames})
        if("${_variableName}" MATCHES "HDF5" OR "${_variableName}" MATCHES "hdf5" OR "${_variableName}" MATCHES "h5")
            message(STATUS "${_variableName}=${${_variableName}}")
        endif()
    endforeach()
    include(cmake-modules/Fetch_Szip.cmake)
    list(APPEND HDF5_LINKER_FLAGS ${SZIP_LIBRARY})


    message(STATUS "HDF5 FOUND IN SYSTEM: ${HDF5_LIBRARIES}")
    message(STATUS "FOUND HDF5")
    message(STATUS "   In path: ${HDF5_ROOT}")
    message(STATUS "   HDF5_DEFINITIONS         : ${HDF5_DEFINITIONS}")
    message(STATUS "   HDF5_C_LIBRARY           : ${HDF5_C_LIBRARY}")
    message(STATUS "   HDF5_C_HL_LIBRARY        : ${HDF5_C_HL_LIBRARY}")
    message(STATUS "   HDF5_CXX_LIBRARY         : ${HDF5_CXX_LIBRARY}")
    message(STATUS "   HDF5_CXX_HL_LIBRARY      : ${HDF5_CXX_HL_LIBRARY}")
    message(STATUS "   HDF5_LINKER_FLAGS        : ${HDF5_LINKER_FLAGS}")


    # Add convenience libraries
    add_library(hdf5::hdf5           UNKNOWN IMPORTED)
    add_library(hdf5::hdf5_hl        UNKNOWN IMPORTED)
    add_library(hdf5::hdf5_cpp       UNKNOWN IMPORTED)
    add_library(hdf5::hdf5_hl_cpp    UNKNOWN IMPORTED)
    add_dependencies(hdf5::hdf5             SZIP)
    add_dependencies(hdf5::hdf5_hl          SZIP)
    add_dependencies(hdf5::hdf5_cpp         SZIP)
    add_dependencies(hdf5::hdf5_hl_cpp      SZIP)


else()
    message(STATUS "HDF5 will be installed into ${INSTALL_DIRECTORY}/hdf5 on first build.")

    include(ExternalProject)
    set(HDF5_IS_PARALLEL OFF)
    ExternalProject_Add(library_HDF5
            URL     https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.3/src/hdf5-1.10.3.tar.bz2 # version 1.10.2
            PREFIX              "${INSTALL_DIRECTORY}/hdf5"
            UPDATE_DISCONNECTED 1
            TEST_COMMAND ""
            CMAKE_ARGS
            -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
            -DBUILD_STATIC_LIBS:BOOL=ON
            -DBUILD_SHARED_LIBS:BOOL=ON
            -DCMAKE_ANSI_CFLAGS:STRING=-fPIC
            -DCMAKE_BUILD_WITH_INSTALL_RPATH:BOOL=OFF
            -DHDF5_ENABLE_PARALLEL=${HDF5_IS_PARALLEL}
            -DALLOW_UNSUPPORTED=ON
            -DBUILD_TESTING:BOOL=OFF
            -DHDF5_ENABLE_Z_LIB_SUPPORT:BOOL=ON
            -DHDF5_BUILD_TOOLS:BOOL=ON
            -DHDF5_BUILD_EXAMPLES:BOOL=OFF
            -DHDF5_BUILD_FORTRAN:BOOL=OFF
            -DHDF5_BUILD_JAVA:BOOL=OFF
            -DCMAKE_INSTALL_MESSAGE=NEVER #Avoid unnecessary output to console
            -DCMAKE_C_FLAGS=-w
            )

    ExternalProject_Get_Property(library_HDF5 INSTALL_DIR)
    add_library(hdf5::hdf5           UNKNOWN IMPORTED)
    add_library(hdf5::hdf5_hl        UNKNOWN IMPORTED)
    add_library(hdf5::hdf5_cpp       UNKNOWN IMPORTED)
    add_library(hdf5::hdf5_hl_cpp    UNKNOWN IMPORTED)

    add_dependencies(hdf5::hdf5          library_HDF5)
    add_dependencies(hdf5::hdf5_hl       library_HDF5)
    add_dependencies(hdf5::hdf5_cpp      library_HDF5)
    add_dependencies(hdf5::hdf5_hl_cpp   library_HDF5)
    set(HDF5_C_LIBRARY        ${INSTALL_DIR}/lib/libhdf5${CUSTOM_SUFFIX})
    set(HDF5_C_HL_LIBRARY     ${INSTALL_DIR}/lib/libhdf5_hl${CUSTOM_SUFFIX})
    set(HDF5_CXX_LIBRARY      ${INSTALL_DIR}/lib/libhdf5_cpp${CUSTOM_SUFFIX})
    set(HDF5_CXX_HL_LIBRARY   ${INSTALL_DIR}/lib/libhdf5_hl_cpp${CUSTOM_SUFFIX})
    set(HDF5_LINKER_FLAGS      -Wl,--no-as-needed -ldl -lm -lz -lsz -lpthread )
#    set(HDF5_LINKER_FLAGS     -lpthread -ldl -lm)
    set(HDF5_INCLUDE_DIR      ${INSTALL_DIR}/include)
#    if (HDF5_IS_PARALLEL)
#        list(APPEND HDF5_LINKER_FLAGS ${MPI_LIBRARIES})
#        list(APPEND HDF5_INCLUDE_DIR  ${MPI_INCLUDE_PATH})
#    endif()
endif()


set_target_properties(hdf5::hdf5 PROPERTIES
        IMPORTED_LOCATION "${HDF5_C_LIBRARY}"
        INTERFACE_LINK_LIBRARIES "${HDF5_LINKER_FLAGS}"
        INCLUDE_DIRECTORIES "${HDF5_INCLUDE_DIR}")

set_target_properties(hdf5::hdf5_hl PROPERTIES
        IMPORTED_LOCATION "${HDF5_C_HL_LIBRARY}"
        INTERFACE_LINK_LIBRARIES "hdf5::hdf5;${HDF5_LINKER_FLAGS}"
        INCLUDE_DIRECTORIES "${HDF5_INCLUDE_DIR}")

set_target_properties(hdf5::hdf5_cpp PROPERTIES
        IMPORTED_LOCATION "${HDF5_CXX_LIBRARY}"
        INTERFACE_LINK_LIBRARIES "hdf5::hdf5;hdf5::hdf5_hl;${HDF5_LINKER_FLAGS}"
        INCLUDE_DIRECTORIES "${HDF5_INCLUDE_DIR}")

set_target_properties(hdf5::hdf5_hl_cpp PROPERTIES
        IMPORTED_LOCATION "${HDF5_CXX_HL_LIBRARY}"
        INTERFACE_LINK_LIBRARIES "hdf5::hdf5;hdf5::hdf5_hl;hdf5::hdf5_cpp;${HDF5_LINKER_FLAGS}"
        INCLUDE_DIRECTORIES "${HDF5_INCLUDE_DIR}")


target_link_libraries(${PROJECT_NAME} PRIVATE hdf5::hdf5)
target_link_libraries(${PROJECT_NAME} PRIVATE hdf5::hdf5_hl)
target_link_libraries(${PROJECT_NAME} PRIVATE hdf5::hdf5_cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE hdf5::hdf5_hl_cpp)

target_include_directories(${PROJECT_NAME} PRIVATE ${HDF5_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} PRIVATE ${HDF5_LINKER_FLAGS})