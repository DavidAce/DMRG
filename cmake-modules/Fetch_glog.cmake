include(cmake-modules/CheckGlogCompiles.cmake)

# Glog can be compiled with or without libunwind
# unwind can sometimes need gcc_eh and lzma as dependencies.
# Yes, it's a mess. Here we add targets if they are present in the system
# If not, we hope for the best making dummy targets.

# ... Additionally, if only lzma is linked we get this weird error when linking ceres:
#       undefined reference to
#       `ceres::GradientProblemSolver::Summary::Summary()'
# This could also be caused by using openblas from anaconda on Clang builds.
# In any case the above error is more likely on Clang.

find_package(Unwind) # If found defines target unwind::unwind
find_library(LZMA_LIB NAMES lzma)
find_library(GCC_EH_LIB NAMES gcc_eh)
add_library(unwind::full INTERFACE IMPORTED)
if(TARGET unwind::unwind AND LZMA_LIB AND GCC_EH_LIB)
    target_link_libraries(unwind::full INTERFACE gcc_eh unwind::unwind lzma)
endif()

# Glog should only look in conda on shared builds! Conda does not give us the static version
set(GLOG_HINTS $ENV{EBROOTGLOG} ${CMAKE_INSTALL_PREFIX} )
if(BUILD_SHARED_LIBS)
    list(APPEND GLOG_HINTS ${CONDA_HINTS})
endif()
find_package(glog 0.4 HINTS ${GLOG_HINTS} PATH_SUFFIXES glog glog/lib NO_DEFAULT_PATH)


if(NOT TARGET glog::glog)
    message(STATUS "Looking for glog in system")
    find_library(GLOG_LIBRARIES     glog           HINTS ${GLOG_HINTS})
    find_path   (GLOG_INCLUDE_DIR   glog/logging.h HINTS ${GLOG_HINTS})
    if(GLOG_LIBRARIES AND GLOG_INCLUDE_DIR)
        check_glog_compiles("lib_header" "" "" "${GLOG_LIBRARIES};unwind::full;gflags::gflags;pthread" "${GLOG_INCLUDE_DIR}" "")
    endif()
    if (GLOG_COMPILES_lib_header)
        add_library(glog::glog ${LINK_TYPE} IMPORTED)
        set_target_properties(glog::glog PROPERTIES IMPORTED_LOCATION "${GLOG_LIBRARIES}")
        target_include_directories(glog::glog SYSTEM INTERFACE ${GLOG_INCLUDE_DIR})
        message(STATUS "Found system glog: Don't forget to also install and link to libraries unwind and lzma")
    endif()
endif()


if(TARGET glog::glog)

elseif(NOT ${DOWNLOAD_METHOD} MATCHES "none")
    message(STATUS "glog will be installed into ${CMAKE_INSTALL_PREFIX}")
    include(${PROJECT_SOURCE_DIR}/cmake-modules/BuildDependency.cmake)
    list(APPEND GLOG_CMAKE_OPTIONS -Dgflags_DIR:PATH=${CMAKE_INSTALL_PREFIX}/gflags)
    build_dependency(glog "${CMAKE_INSTALL_PREFIX}" "${GLOG_CMAKE_OPTIONS}")
    find_package(glog 0.4
            HINTS ${CMAKE_INSTALL_PREFIX}
            PATH_SUFFIXES glog glog/lib
            NO_DEFAULT_PATH)
    if(TARGET glog::glog)
        message(STATUS "glog installed successfully")
    else()
        message(FATAL_ERROR "glog could not be downloaded.")
    endif()

else()
    message(FATAL_ERROR "Dependency glog not found and DOWNLOAD_METHOD = ${DOWNLOAD_METHOD}")
endif()


if(TARGET glog::glog)
    get_target_property(GLOG_TYPE glog::glog TYPE)
    if(GLOG_TYPE MATCHES "SHARED" AND NOT BUILD_SHARED_LIBS)
        include(cmake-modules/PrintTargetProperties.cmake)
        print_target_properties(glog::glog)
        message(FATAL_ERROR "Found shared glog library on a static build!")
    endif()

    include(cmake-modules/TargetFilters.cmake)
    remove_library_shallow(glog::glog "Threads::Threads|pthread|unwind|gflags")
    target_link_libraries(glog::glog INTERFACE unwind::full gflags::gflags pthread )

    # Modernize
    get_property(imp_loc_set TARGET glog::glog PROPERTY IMPORTED_LOCATION SET) # Returns a boolean if set
    get_property(loc_set     TARGET glog::glog PROPERTY LOCATION SET) # Returns a boolean if set
    if(loc_set AND NOT imp_loc_set)
        get_target_property(imp_loc glog::glog LOCATION)
        set_target_properties(glog::glog PROPERTIES IMPORTED_LOCATION ${imp_loc})
    endif()

    check_glog_compiles("target" "${GLOG_FLAGS}" "" "" "" "glog::glog")
    if(NOT GLOG_COMPILES_target)
        include(cmake-modules/PrintTargetProperties.cmake)
        print_target_properties(glog::glog)
        message(FATAL_ERROR "Could not compile a simple glog program")
    endif()
endif()