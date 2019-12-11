function(find_Arpackpp)
    if (NOT TARGET arpack++)
        message(STATUS "Searching for Arpack++ headers")
        find_path(ARPACKPP_INCLUDE_DIR
                NAMES arpack++/ardsnsym.h arpackpp/ardsnsym.h
                HINTS ${DIRECTORY_HINTS}
                PATHS
                    ${CMAKE_INSTALL_PREFIX}
                    $ENV{EBROOTARPACKPLUSPLUS}
                    $ENV{ARPACKPP_DIR}
                PATH_SUFFIXES
                    arpack++/include/arpack++
                    arpackpp/include/arpackpp
                    arpack++/include arpackpp/include
                    include/arpack++ include/arpackpp
                    arpack++ arpackpp include
                NO_DEFAULT_PATH
                )
        if(ARPACKPP_INCLUDE_DIR)
            message(STATUS "Searching for Arpack++ headers - Success: ${ARPACKPP_INCLUDE_DIR}")
            add_library(arpack++ INTERFACE)
            target_include_directories(arpack++ SYSTEM INTERFACE ${ARPACKPP_INCLUDE_DIR})
        else()
            message(STATUS "Searching for Arpack++ headers - failed")
        endif()
    endif()


    if(NOT TARGET arpack++)
        message(STATUS "Searching for Arpack++ in system")
        find_library(ARPACKPP_LIBRARIES
                NAMES arpackpp arpack++ libarpack2++ libarpack++ libarpackpp
                HINTS ${DIRECTORY_HINTS}
                PATHS
                    ${CMAKE_INSTALL_PREFIX}
                    $ENV{EBROOTARPACKPLUSPLUS}
                    $ENV{ARPACKPP_DIR}
                PATH_SUFFIXES arpack++/lib arpackpp/lib lib lib32 lib64 x86_64-linux-gnu lib/x86_64-linux-gnu
                )
        find_path(ARPACKPP_INCLUDE_DIR
                NAMES ardsnsym.h
                HINTS ${DIRECTORY_HINTS}
                PATHS
                    ${CMAKE_INSTALL_PREFIX}
                    $ENV{EBROOTARPACKPLUSPLUS}
                    $ENV{ARPACKPP_DIR}
                PATH_SUFFIXES
                    arpack++/include/arpack++
                    arpackpp/include/arpackpp
                    arpack++/include arpackpp/include
                    include/arpack++ include/arpackpp
                    arpack++ arpackpp include
                NO_DEFAULT_PATH
                )
        if(ARPACKPP_LIBRARIES AND ARPACKPP_INCLUDE_DIR)
            message(STATUS "Searching for Arpack++ in system - Success: ${ARPACKPP_LIBRARIES}")
            message(STATUS "Note that old versions of Arpack++ (e.g. the default in Ubuntu Trusty 14.04 LTS) may fail to compile, they require '-fpermissive'.")
            add_library(arpack++ INTERFACE IMPORTED)
            target_link_libraries(arpack++ INTERFACE ${ARPACKPP_LIBRARIES})
            target_include_directories(arpack++ SYSTEM INTERFACE ${ARPACKPP_INCLUDE_DIR})
        else()
            message(STATUS "Searching for Arpack++ in system - failed")
        endif()
    endif()
endfunction()