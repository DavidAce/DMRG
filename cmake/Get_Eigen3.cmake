function(find_eigen3)
    find_package(Eigen3 3.3.7 ${N5} ${N6} ${N7} ${N8} ${REQUIRED}) # Flags ignore system packages. See cmake/SetupPaths.cmake

    if(NOT Eigen3_FOUND AND DMRG_PACKAGE_MANAGER MATCHES "cmake")
        message(STATUS "Eigen3 will be installed into ${DMRG_DEPS_INSTALL_DIR}")
        include(cmake/InstallPackage.cmake)
        install_package(Eigen3 "${DMRG_DEPS_INSTALL_DIR}" "")
        find_package(Eigen3 3.3.7 HINTS ${DMRG_DEPS_INSTALL_DIR} NO_DEFAULT_PATH REQUIRED)
    endif()
    if(Eigen3_FOUND AND TARGET Eigen3::Eigen)
        message(STATUS "Found Eigen3: ${EIGEN3_INCLUDE_DIR}")
        target_include_directories(Eigen3::Eigen SYSTEM INTERFACE ${EIGEN3_INCLUDE_DIR})
    else()
        message(FATAL_ERROR "Eigen3 not found: ${Eigen3_FOUND}")
    endif()
endfunction()

find_eigen3()

#
#if(NOT TARGET Eigen3::Eigen AND DMRG_PACKAGE_MANAGER STREQUAL "find")
#    message(STATUS "DMRG needs a patched version of Eigen 3.3.7 and will downloaded despite DMRG_PACKAGE_MANAGER=find")
#    set(EIGEN_REQUIRED REQUIRED)
#endif()
#
#if(NOT TARGET Eigen3::Eigen AND DMRG_PACKAGE_MANAGER MATCHES "find|cmake")
#    # We want to find our own Eigen3 to make sure we patch it properly
#    find_package(Eigen3
#        HINTS ${DMRG_DEPS_INSTALL_DIR}
#        NO_SYSTEM_PATH # IMPORTANT TO ONLY LOOK IN DMRG'S OWN MODULE DIRECTORY
#        ${EIGEN_REQUIRED})
#    if(Eigen3_FOUND AND TARGET Eigen3::Eigen)
#        message(STATUS "Found Eigen3: ${EIGEN3_INCLUDE_DIR}")
#        target_include_directories(Eigen3::Eigen SYSTEM INTERFACE ${EIGEN3_INCLUDE_DIR})
#    endif()
#endif()
#
#if(NOT Eigen3_FOUND OR NOT TARGET Eigen3::Eigen AND DMRG_PACKAGE_MANAGER MATCHES "cmake")
#    message(STATUS "Eigen3 will be installed into ${DMRG_DEPS_INSTALL_DIR}")
#    include(cmake/InstallPackage.cmake)
#    install_package(Eigen3 "${DMRG_DEPS_INSTALL_DIR}" "")
#    find_package(Eigen3 3.3.7
#            HINTS ${DMRG_DEPS_INSTALL_DIR}
#            NO_SYSTEM_PATH
#            REQUIRED)
#    if(TARGET Eigen3::Eigen)
#        message(STATUS "Eigen3 installed successfully")
#        target_include_directories(Eigen3::Eigen SYSTEM INTERFACE ${EIGEN3_INCLUDE_DIR})
#    else()
#        message(FATAL_ERROR "Eigen3 could not be installed")
#    endif()
#endif()


