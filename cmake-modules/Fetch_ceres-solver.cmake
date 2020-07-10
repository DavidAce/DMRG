if(DMRG_DOWNLOAD_METHOD MATCHES "find|fetch")
    foreach (tgt glog::glog;gflags::gflags;Eigen3::Eigen)
        if(NOT TARGET ${tgt})
            list(APPEND CERES_MISSING_TARGET ${tgt})
            mark_as_advanced(CERES_MISSING_TARGETS)
        endif()
    endforeach()
    if(CERES_MISSING_TARGET)
        message(FATAL_ERROR "Ceres: dependencies missing [${CERES_MISSING_TARGET}]")
    endif()
endif()


if(NOT TARGET Ceres::ceres AND NOT TARGET ceres AND DMRG_DOWNLOAD_METHOD MATCHES "find|fetch")
    include(cmake-modules/CheckCeresCompiles.cmake)

    # Can't use conda here since they only have shared libraries.
    # Can't use config-mode on anaconda either since they have weird requirements
    # on predefined targets for gflags and glog. Therefore -> NO_DEFAULT_PATH
    find_package(Ceres
            HINTS ${CMAKE_INSTALL_PREFIX}
            PATH_SUFFIXES ceres ceres/lib
            NO_DEFAULT_PATH)

    if(NOT TARGET ceres)
        message(STATUS "Looking for ceres in system")
        find_path   (CERES_INCLUDE_DIR        NAMES  ceres/ceres.h                   )
        find_path   (SUITESPARSE_INCLUDE_DIR  NAMES  suitesparse/SuiteSparse_config.h)
        if(CERES_INCLUDE_DIR AND SUITESPARSE_INCLUDE_DIR)
            # We may have a chance at finding Ceres in the system
            find_library(CERES_LIB                NAMES ceres                        )
            find_library(SUITESPARSE_LIB          NAMES suitesparse suitesparseconfig)
            find_library(CXSPARSE_LIB             NAMES cxsparse                     )
            find_library(METIS_LIB                NAMES metis metis_static           )
            find_library(CHOLMOD_LIB              NAMES cholmod                      )
            find_library(COLAMD_LIB               NAMES colamd                       )
            find_library(CCOLAMD_LIB              NAMES ccolamd                      )
            find_library(AMD_LIB                  NAMES amd                          )
            find_library(CAMD_LIB                 NAMES camd                         )
            # Make sure this is the correct linking order
            set(ceres_lib_names ceres cholmod colamd ccolamd amd camd metis cxsparse suitesparse)
            # Check that all libs are present before going forward with defining targets
            foreach(lib ${ceres_lib_names})
                string(TOUPPER ${lib} LIB)
                if(NOT EXISTS "${${LIB}_LIB}")
                    set(CERES_MISSING_LIBS TRUE)
                    message(STATUS "Missing Ceres dependency in system: ${lib}")
                endif()
            endforeach()
            if(NOT CERES_MISSING_LIBS)
                # Now we are confident all libraries are in the system. Note that metis (static) is not available on bionic
                # which is why we do this
                foreach(lib ${ceres_lib_names})
                    string(TOUPPER ${lib} LIB)
                    if(EXISTS "${${LIB}_LIB}")
                        add_library(ceres::${lib} ${LINK_TYPE} IMPORTED)
                        set_target_properties(ceres::${lib} PROPERTIES IMPORTED_LOCATION "${${LIB}_LIB}")
                        if(NOT "${lib}" MATCHES "ceres")
                            target_link_libraries(Ceres::ceres INTERFACE ceres::${lib})
                        endif()
                        if(NOT "${lib}" MATCHES "ceres|suitesparse")
                            target_link_libraries(ceres::${lib} INTERFACE ceres::suitesparse)
                        endif()
                    endif()
                endforeach()
                target_include_directories(Ceres::ceres SYSTEM INTERFACE ${CERES_INCLUDE_DIR})
                target_include_directories(Ceres::ceres SYSTEM INTERFACE ${SUITESPARSE_INCLUDE_DIR})
            endif()
        endif()
    endif()
    if(TARGET Ceres::ceres OR TARGET ceres)
        message(STATUS "Found Ceres")
    endif()
endif()


if(NOT TARGET Ceres::ceres AND NOT TARGET ceres AND DMRG_DOWNLOAD_METHOD MATCHES "fetch")
    message(STATUS "Ceres will be installed into ${CMAKE_INSTALL_PREFIX} on first build.")
    get_target_property(EIGEN3_INCLUDE_DIR Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
    list (GET EIGEN3_INCLUDE_DIR 0 EIGEN3_INCLUDE_DIR)
    get_target_property(GFLAGS_LIBRARIES    gflags::gflags  IMPORTED_LOCATION)
    get_target_property(GFLAGS_INCLUDE_DIR  gflags::gflags  INTERFACE_INCLUDE_DIRECTORIES)
    get_target_property(GLOG_LIBRARIES      glog::glog      IMPORTED_LOCATION)
    get_target_property(GLOG_INCLUDE_DIR    glog::glog      INTERFACE_INCLUDE_DIRECTORIES)
    if(NOT GLOG_LIBRARIES)
        include(cmake-modules/PrintTargetProperties.cmake)
        print_target_properties(glog::glog)
        message(FATAL_ERROR "Could not extract glog libraries: ${GLOG_LIBRARIES}")
    endif()
    if(NOT GFLAGS_LIBRARIES)
        include(cmake-modules/PrintTargetInfo.cmake)
        print_target_info(gflags::gflags)
        message(FATAL_ERROR "Could not extract gflags libraries: ${GFLAGS_LIBRARIES}")
    endif()

    list(APPEND CERES_CMAKE_OPTIONS  -DEIGEN_INCLUDE_DIR_HINTS:PATH=${CMAKE_INSTALL_PREFIX})
    list(APPEND CERES_CMAKE_OPTIONS  -DEigen3_DIR:PATH=${CMAKE_INSTALL_PREFIX}/Eigen3/share/eigen3/cmake)
    list(APPEND CERES_CMAKE_OPTIONS  -Dgflags_DIR:PATH=${CMAKE_INSTALL_PREFIX}/gflags)
    list(APPEND CERES_CMAKE_OPTIONS  -DGFLAGS_INCLUDE_DIR_HINTS:PATH=${CMAKE_INSTALL_PREFIX})
    list(APPEND CERES_CMAKE_OPTIONS  -DGFLAGS_LIBRARY_DIR_HINTS:PATH=${CMAKE_INSTALL_PREFIX})
    list(APPEND CERES_CMAKE_OPTIONS  -DGFLAGS_LIBRARY:PATH=${GFLAGS_LIBRARIES})
    list(APPEND CERES_CMAKE_OPTIONS  -Dglog_DIR:PATH=${CMAKE_INSTALL_PREFIX}/glog)
    list(APPEND CERES_CMAKE_OPTIONS  -DGLOG_INCLUDE_DIR_HINTS:PATH=${CMAKE_INSTALL_PREFIX})
    list(APPEND CERES_CMAKE_OPTIONS  -DGLOG_LIBRARY_DIR_HINTS:PATH=${CMAKE_INSTALL_PREFIX})
    list(APPEND CERES_CMAKE_OPTIONS  -DGLOG_LIBRARY:PATH=${GLOG_LIBRARIES})
    list(APPEND CERES_CMAKE_OPTIONS  -DCMAKE_VERBOSE_MAKEFILE:BOOL=${CMAKE_VERBOSE_MAKEFILE})
    message(STATUS "ceres options: ${CERES_CMAKE_OPTIONS}")
    include(${PROJECT_SOURCE_DIR}/cmake-modules/BuildDependency.cmake)
    build_dependency(ceres "${CMAKE_INSTALL_PREFIX}" "${CERES_CMAKE_OPTIONS}" )
    find_package(Ceres HINTS ${CMAKE_INSTALL_PREFIX} NO_DEFAULT_PATH)
    if(TARGET Ceres::ceres)
        message(STATUS "Ceres installed successfully")
    else()
        message(FATAL_ERROR "Ceres could not be installed")
    endif()
endif()



#if(TARGET ceres AND NOT TARGET Ceres::ceres )
#    # Use this for the ceres targets defined by CONFIG mode find_package
#    # These find_packages have the tendency to do the wrong thing, like
#    #   - injecting shared libraries into static builds
#    #   - using "-lpthread" instead of "pthread"
#
#    get_target_property(CERES_TYPE ceres TYPE)
#    if(CERES_TYPE MATCHES "SHARED" AND NOT BUILD_SHARED_LIBS)
#        include(cmake-modules/PrintTargetProperties.cmake)
#        print_target_properties(ceres)
#        message(FATAL_ERROR "Found shared ceres library on a static build!")
#    endif()
#
#    #Remove any shared libraries like unwind etc which pollute static builds
#    # As a matter of fact... just relink it entirely
#    include(cmake-modules/TargetFilters.cmake)
#    remove_library_shallow(ceres "gcc_eh|unwind|lzma|Threads::Threads|pthread|glog|gflags")
#
##    if(NOT BUILD_SHARED_LIBS)
##    include(cmake-modules/TargetFilters.cmake)
##    remove_library_shallow(ceres "gcc_eh|unwind|lzma|Threads::Threads|pthread|unwind|glog|gflags")
##    target_link_libraries(ceres INTERFACE gcc_eh unwind lzma glog::glog gflags::gflags pthread )
##    endif()
#
#    # Modernize
##    get_property(imp_loc_set TARGET ceres PROPERTY IMPORTED_LOCATION SET) # Returns a boolean if set
##    get_property(loc_set     TARGET ceres PROPERTY LOCATION SET) # Returns a boolean if set
##    if(loc_set AND NOT imp_loc_set)
##        get_target_property(imp_loc ceres LOCATION)
##        set_target_properties(ceres PROPERTIES IMPORTED_LOCATION ${imp_loc})
##    endif()
#
#
#    # Copy ceres to Ceres::ceres to follow proper naming convention
#    include(cmake-modules/CopyTarget.cmake)
#    copy_target(Ceres::ceres ceres)
#endif()


if(TARGET Ceres::ceres)
    include(cmake-modules/CheckCeresCompiles.cmake)
    check_ceres_compiles("Ceres::ceres" "" "" "" "")
    if(NOT CERES_COMPILES)
        include(cmake-modules/PrintTargetProperties.cmake)
        print_target_properties(Ceres::ceres)
        message(FATAL_ERROR "Could not compile simple ceres program")
    endif()
endif()