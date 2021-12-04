
function(find_arpackpp)
    if (NOT TARGET ARPACK::ARPACK)
        message(FATAL_ERROR "arpack++ needs to link to target [ARPACK::ARPACK]")
    endif()
    if (NOT TARGET lapacke::lapacke)
        message(FATAL_ERROR "arpack++ needs to link to target [lapacke::lapacke]")
    endif()

    if(ARPACKPP_NO_DEFAULT_PATH)
        set(NO_DEFAULT_PATH NO_DEFAULT_PATH)
    endif()
    if(ARPACKPP_NO_CMAKE_PACKAGE_REGISTRY)
        set(NO_CMAKE_PACKAGE_REGISTRY NO_CMAKE_PACKAGE_REGISTRY)
    endif()

    if (NOT TARGET arpack++::arpack++)
        include(GNUInstallDirs)
        unset(ARPACKPP_LIBRARY)
        unset(ARPACKPP_LIBRARY CACHE)
        unset(ARPACKPP_INCLUDE_DIR)
        unset(ARPACKPP_INCLUDE_DIR CACHE)
        find_library(ARPACKPP_LIBRARY
                NAMES arpackpp arpack++
                ${NO_DEFAULT_PATH}
                ${NO_CMAKE_PACKAGE_REGISTRY}
                )
        find_path(ARPACKPP_INCLUDE_DIR
                NAMES arpack++/arcomp.h arpackpp/arcomp.h
                PATH_SUFFIXES
                    include
                    arpack++ arpackpp
                    arpack++/include arpackpp/include
                    include/arpack++ include/arpackpp
                    arpack++/include/arpack++
                    arpackpp/include/arpackpp
                ${NO_DEFAULT_PATH}
                ${NO_CMAKE_PACKAGE_REGISTRY}
                )

        if(NOT ARPACKPP_LIBRARY)
            set(ARPACKPP_LIBRARY "")
        endif()
        if(NOT ARPACKPP_INCLUDE_DIR)
            set(ARPACKPP_INCLUDE_DIR "")
        endif()
    endif()
endfunction()



find_arpackpp()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(arpack++ DEFAULT_MSG ARPACKPP_INCLUDE_DIR)

if(arpack++_FOUND)
    if(NOT TARGET arpack++::arpack++)
        if(ARPACKPP_LIBRARY)
            add_library(arpack++::arpack++ UNKNOWN IMPORTED)
            set_target_properties(arpack++::arpack++ PROPERTIES IMPORTED_LOCATION ${ARPACKPP_LIBRARY})
        else()
            add_library(arpack++::arpack++ INTERFACE IMPORTED)
        endif()
    endif()
    target_include_directories(arpack++::arpack++ SYSTEM INTERFACE ${ARPACKPP_INCLUDE_DIR})
    target_link_libraries(arpack++::arpack++ INTERFACE lapacke::lapacke ARPACK::ARPACK)
endif()