
# This script attempts to find HDF5 installed from apt or from sources
# In both cases, if HDF5 is found a target "hdf5" is generated so the user can simply do:
#
#  target_link_libraries(mylibrary INTERFACE hdf5)
#
#
#  The user can guide the find pattern with variables:
#       BUILD_SHARED_LIBS           ON/OFF for shared/static libs
#       HDF5_REQUIRED               to require HDF5 to be found, set to ON
#       HDF5_WANT_VERSION           sets the required version (default 1.10)
#




if(BUILD_SHARED_LIBS)
    set(HDF5_TARGET_SUFFIX "shared")
    set(HDF5_LIBRARY_SUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
    set(HDF5_USE_STATIC_LIBRARIES OFF)
else()
    set(HDF5_TARGET_SUFFIX "static")
    set(HDF5_LIBRARY_SUFFIX ${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(HDF5_USE_STATIC_LIBRARIES ON)
endif()
set(HDF5_FIND_DEBUG OFF)
if(NOT HDF5_WANT_VERSION)
    set(HDF5_WANT_VERSION 1.10)
endif()

set(HDF5_ROOT ${HDF5_ROOT} $ENV{HDF5_ROOT} $ENV{HDF5_DIR} $ENV{EBROOTHDF5} $ENV{HOME}/.conda $ENV{HOME}/anaconda3 $ENV{HOME}/miniconda3 /usr /usr/local $ENV{PATH})
set(HDF5_DIR  ${HDF5_ROOT})
find_file(HDF5_CXX_COMPILER_EXECUTABLE      NAMES h5c++ PATH_SUFFIXES bin PATHS ${HDF5_ROOT})


if(HDF5_REQUIRED)
    find_package(HDF5 ${HDF5_WANT_VERSION} COMPONENTS C HL REQUIRED)
else()
    find_package(HDF5 ${HDF5_WANT_VERSION} COMPONENTS C HL)
endif()

# To print all variables, use the code below:
##
#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    if("${_variableName}" MATCHES "HDF5" OR "${_variableName}" MATCHES "hdf5" OR "${_variableName}" MATCHES "h5")
#        message(STATUS "${_variableName}=${${_variableName}}")
#    endif()
#endforeach()


if(HDF5_FOUND)
    # Add convenience libraries to collect all the hdf5 libraries
    add_library(hdf5    INTERFACE)
    add_library(hdf5::hdf5 ALIAS hdf5)
    target_link_libraries(hdf5
            INTERFACE
            ${HDF5_hdf5_hl_LIBRARY}
            ${HDF5_C_LIBRARY_hdf5_hl}
            ${HDF5_C_HL_LIBRARY}
            ${HDF5_hdf5_LIBRARY}
            ${HDF5_C_LIBRARY}
            ${HDF5_C_LIBRARY_hdf5}
            $<LINK_ONLY:-ldl -lm -lz>
            Threads::Threads
            )
    target_include_directories(hdf5 INTERFACE  ${HDF5_INCLUDE_DIR})
#    if(HDF5_ENABLE_Z_LIB_SUPPORT)
#        target_link_libraries(hdf5 INTERFACE $<LINK_ONLY:-lz>  )
#    endif()
    if (_HDF5_LPATH)
        set(HDF5_ROOT ${_HDF5_LPATH})
    endif()
    if(HDF5_C_LIBRARY_sz)
        target_link_libraries(hdf5 INTERFACE $<LINK_ONLY:-lsz>  )
    endif()
    if(HDF5_IS_PARALLEL)
        if(HDF5_USE_STATIC_LIBRARIES)
            message(WARNING "The HDF5 library is parallel and you are linking with -static. You have to link to MPI!")
        endif()
#        find_package(MPI)
#        target_link_libraries(hdf5 INTERFACE -lmpi )
        #        target_link_libraries(hdf5 INTERFACE MPI::MPI_C )
#        get_target_property(MPI_LIB MPI::MPI_C INTERFACE_LINK_LIBRARIES)
#        message("Found MPI Library: ${MPI_LIB}")
#        find_package(MPI REQUIRED)
#        target_link_libraries(hdf5 INTERFACE MPI::MPI_C )
#        get_target_property(MPI_LIB MPI::MPI_C INTERFACE_LINK_LIBRARIES)
#        message("Found MPI Library: ${MPI_LIB}")
#        include(cmake-modules/PrintTargetProperties.cmake)
#        print_target_properties(MPI::MPI_C)
    endif()
#    if(HDF5_C_LIBRARY_z)
#        target_link_libraries(hdf5 INTERFACE $<LINK_ONLY:-lz>  )
#    endif()

#
#    if(TARGET hdf5::hdf5-${HDF5_TARGET_SUFFIX} AND DEFINED HDF5_C_HL_LIBRARY AND DEFINED HDF5_C_LIBRARY)
#        set(HDF5_DIR              ${HDF5_BUILD_DIR}/share/cmake/hdf5)
#        set(HDF5_ROOT             ${HDF5_BUILD_DIR})
#
#        target_link_libraries(hdf5
#                INTERFACE
#                ${HDF5_C_HL_LIBRARY}
#                ${HDF5_C_LIBRARY}
#                $<LINK_ONLY:-ldl -lm>
#                ${PTHREAD_LIBRARY}
#                )
#        if(HDF5_ENABLE_Z_LIB_SUPPORT)
#            target_link_libraries(hdf5 INTERFACE $<LINK_ONLY:-lz>  )
#        endif()
#        target_include_directories(hdf5 INTERFACE  ${HDF5_INCLUDE_DIR})
#    elseif(DEFINED HDF5_hdf5_hl_LIBRARY AND DEFINED HDF5_hdf5_LIBRARY)
#
#        target_link_libraries(hdf5
#                INTERFACE
#                ${HDF5_hdf5_hl_LIBRARY}
#                ${HDF5_hdf5_LIBRARY}
#                $<LINK_ONLY:-ldl -lm -lz>
#                ${PTHREAD_LIBRARY}
#                )
#        target_include_directories(hdf5 INTERFACE  ${HDF5_INCLUDE_DIR})
#
#    elseif(DEFINED HDF5_C_LIBRARY_hdf5 AND DEFINED HDF5_C_LIBRARY_hdf5_hl)
#        #        add_dependencies(hdf5  SZIP)
#        if (_HDF5_LPATH)
#            set(HDF5_ROOT ${_HDF5_LPATH})
#        endif()
#
#        target_link_libraries(hdf5
#                INTERFACE
#                ${HDF5_C_LIBRARY_hdf5_hl}
#                ${HDF5_C_LIBRARY_hdf5}
#                $<LINK_ONLY:-ldl -lm>
#                ${PTHREAD_LIBRARY}
#                )
#        if(HDF5_C_LIBRARY_sz)
#            target_link_libraries(hdf5 INTERFACE $<LINK_ONLY:-lsz>  )
#        endif()
#        if(HDF5_C_LIBRARY_z)
#            target_link_libraries(hdf5 INTERFACE $<LINK_ONLY:-lz>  )
#        endif()
#        target_include_directories(hdf5 INTERFACE ${HDF5_INCLUDE_DIR})
#    endif()
endif()
