
if(NOT TARGET Eigen3::Eigen AND DMRG_PACKAGE_MANAGER STREQUAL "find")
    message(STATUS "DMRG needs a patched version of Eigen 3.3.7 and will downloaded despite DMRG_PACKAGE_MANAGER=find")
endif()

if(NOT TARGET Eigen3::Eigen AND DMRG_PACKAGE_MANAGER MATCHES "find|cmake")
    # We want to find our own Eigen3 to make sure we patch it properly
    find_package(Eigen3
        HINTS ${CMAKE_INSTALL_PREFIX}
        NO_DEFAULT_PATH) # IMPORTANT TO ONLY LOOK IN DMRG'S OWN PLACE
    if(TARGET Eigen3::Eigen)
        message(STATUS "Found Eigen3: ${EIGEN3_INCLUDE_DIR}")
        target_include_directories(Eigen3::Eigen SYSTEM INTERFACE ${EIGEN3_INCLUDE_DIR})
    endif()
endif()

if(NOT TARGET Eigen3::Eigen AND DMRG_PACKAGE_MANAGER MATCHES "find|cmake")
    message(STATUS "Eigen3 will be installed into ${CMAKE_INSTALL_PREFIX}")
    include(${PROJECT_SOURCE_DIR}/cmake-modules/BuildDependency.cmake)
    build_dependency(Eigen3 "${CMAKE_INSTALL_PREFIX}" "")
    find_package(Eigen3 3.3.7
            HINTS ${CMAKE_INSTALL_PREFIX}
            NO_DEFAULT_PATH)
    if(TARGET Eigen3::Eigen)
        message(STATUS "Eigen3 installed successfully")
        target_include_directories(Eigen3::Eigen SYSTEM INTERFACE ${EIGEN3_INCLUDE_DIR})
    else()
        message(FATAL_ERROR "Eigen3 could not be installed")
    endif()
endif()

