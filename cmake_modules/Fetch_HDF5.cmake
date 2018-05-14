
# Try finding in system. While avoiding Python anaconda's version.
# Usually there is an executable present in system, "h5c++"
# Execute "find /usr -name h5c++ -print -quit" in terminal to find out where it is.
# -quit returns after first match


#execute_process(
#        COMMAND  find /usr -name libhdf5.* -exec dirname {} \; -quit
#        OUTPUT_VARIABLE HDF5_ROOT
#)

#Try to find the compiler wrapper
#execute_process(
#        COMMAND  find /usr/bin -name h5c++ -print -quit
#        OUTPUT_VARIABLE HDF5_CXX_COMPILER_EXECUTABLE
#)
#string(STRIP "${HDF5_CXX_COMPILER_EXECUTABLE}" HDF5_CXX_COMPILER_EXECUTABLE)

set(HDF5_USE_STATIC_LIBRARIES ON)
set(HDF5_FIND_DEBUG OFF)
find_package(HDF5 COMPONENTS CXX HL)
if(HDF5_FOUND AND HDF5_LIBRARIES AND HDF5_CXX_LIBRARIES AND HDF5_HL_LIBRARIES AND HDF5_CXX_HL_LIBRARIES)
    message(STATUS "HDF5 FOUND IN SYSTEM: ${HDF5_LIBRARIES}")
    add_definitions(${HDF5_DEFINITIONS})
    list(APPEND HDF5_LDFLAGS ${HDF5_CXX_LIBRARY_NAMES} ${HDF5_CXX_HL_LIBRARY_NAMES})
    list(APPEND HDF5_LIBRARIES ${HDF5_HL_LIBRARIES})
    list(REMOVE_DUPLICATES HDF5_LIBRARIES)

    message(STATUS "FOUND HDF5")
    message(STATUS "   In path: ${HDF5_INCLUDE_DIR}")
    message(STATUS "   HDF5 DEFINITIONS: ${HDF5_DEFINITIONS}")
    message(STATUS "   HDF5 LIBRARIES  : ${HDF5_LIBRARIES}")

    target_include_directories(${PROJECT_NAME} PRIVATE ${HDF5_INCLUDE_DIR})
    target_link_libraries(${PROJECT_NAME} ${HDF5_LIBRARIES} -ldl)
else()
    message(STATUS "HDF5 will be installed into ${INSTALL_DIRECTORY}/hdf5 on first build.")

    include(ExternalProject)
    ExternalProject_Add(library_HDF5
            URL      https://www.hdfgroup.org/package/source-bzip2-2/?wpdmdl=11811 # version 1.10.2
            PREFIX              "${INSTALL_DIRECTORY}/hdf5"
            UPDATE_DISCONNECTED 1
            TEST_COMMAND ""
            CMAKE_ARGS
            -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
            -DBUILD_SHARED_LIBS:BOOL=OFF
            -DHDF5_ENABLE_PARALLEL=OFF
            -DALLOW_UNSUPPORTED=ON
            -DBUILD_TESTING:BOOL=OFF
            -DHDF5_BUILD_TOOLS:BOOL=OFF
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


    set_target_properties(hdf5::hdf5 PROPERTIES
            IMPORTED_LOCATION ${INSTALL_DIR}/lib/libhdf5${CMAKE_STATIC_LIBRARY_SUFFIX}
            INTERFACE_LINK_LIBRARIES "m;dl;dl"
            INCLUDE_DIRECTORIES ${INSTALL_DIR}/include)

    set_target_properties(hdf5::hdf5_hl PROPERTIES
            IMPORTED_LOCATION ${INSTALL_DIR}/lib/libhdf5_hl${CMAKE_STATIC_LIBRARY_SUFFIX}
            INTERFACE_LINK_LIBRARIES hdf5::hdf5
            INCLUDE_DIRECTORIES ${INSTALL_DIR}/include)

    set_target_properties(hdf5::hdf5_cpp PROPERTIES
            IMPORTED_LOCATION ${INSTALL_DIR}/lib/libhdf5_cpp${CMAKE_STATIC_LIBRARY_SUFFIX}
            INTERFACE_LINK_LIBRARIES hdf5::hdf5
            INCLUDE_DIRECTORIES ${INSTALL_DIR}/include)

    set_target_properties(hdf5::hdf5_hl_cpp PROPERTIES
            IMPORTED_LOCATION ${INSTALL_DIR}/lib/libhdf5_hl_cpp${CMAKE_STATIC_LIBRARY_SUFFIX}
            INTERFACE_LINK_LIBRARIES "hdf5::hdf5_hl;hdf5::hdf5"
            INCLUDE_DIRECTORIES ${INSTALL_DIR}/include)


    add_dependencies(hdf5::hdf5          library_HDF5)
    add_dependencies(hdf5::hdf5_hl       library_HDF5)
    add_dependencies(hdf5::hdf5_cpp      library_HDF5)
    add_dependencies(hdf5::hdf5_hl_cpp   library_HDF5)

    target_link_libraries(${PROJECT_NAME} hdf5::hdf5)
    target_link_libraries(${PROJECT_NAME} hdf5::hdf5_hl)
    target_link_libraries(${PROJECT_NAME} hdf5::hdf5_cpp)
    target_link_libraries(${PROJECT_NAME} hdf5::hdf5_hl_cpp)
    target_link_libraries(${PROJECT_NAME} -ldl)
    target_include_directories(${PROJECT_NAME} PRIVATE ${INSTALL_DIR}/include)
    #For convenience, define these variables
    get_target_property(HDF5_LIBRARIES          hdf5::hdf5        IMPORTED_LOCATION)
    get_target_property(HDF5_HL_LIBRARIES       hdf5::hdf5_hl     IMPORTED_LOCATION)
    get_target_property(HDF5_CXX_LIBRARIES      hdf5::hdf5_cpp    IMPORTED_LOCATION)
    get_target_property(HDF5_CXX_HL_LIBRARIES   hdf5::hdf5_hl_cpp IMPORTED_LOCATION)
endif()



