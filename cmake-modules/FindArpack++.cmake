function(CheckArpackppCompiles TAG REQUIRED_FLAGS REQUIRED_DEFINITIONS REQUIRED_LIBRARIES_UNPARSED REQUIRED_INCLUDES REQUIRED_TARGET)
    #    message(STATUS "Checking if arpack headers are working")
    set(REQUIRED_LIBRARIES)
    include(cmake-modules/getExpandedTarget.cmake)
    foreach(elem ${REQUIRED_LIBRARIES_UNPARSED})
        if(TARGET ${elem})
            expand_target_libs(${elem} expanded_libs)
            list(APPEND REQUIRED_LIBRARIES ${expanded_libs})
        else()
            list(APPEND REQUIRED_LIBRARIES ${elem})
        endif()
    endforeach()

    foreach(elem ${REQUIRED_TARGET})
        if(TARGET ${elem})
            expand_target_libs(${elem} expanded_libs)
            expand_target_incs(${elem} expanded_incs)
            expand_target_opts(${elem} expanded_opts)
            list(APPEND REQUIRED_LIBRARIES ${expanded_libs})
            list(APPEND REQUIRED_INCLUDES  ${expanded_incs})
            list(APPEND REQUIRED_FLAGS     ${expanded_opts})
        endif()
    endforeach()
    string(REPLACE ";" " " REQUIRED_FLAGS      "${REQUIRED_FLAGS}") # Needs to be a space-separated list


    #   Test features
    include(CheckCXXSourceCompiles)
    set(CMAKE_REQUIRED_FLAGS        ${REQUIRED_FLAGS})
    set(CMAKE_REQUIRED_DEFINITIONS  ${REQUIRED_DEFINITIONS})
    set(CMAKE_REQUIRED_LIBRARIES    ${REQUIRED_LIBRARIES})
    set(CMAKE_REQUIRED_INCLUDES     ${REQUIRED_INCLUDES})
    message(STATUS "ARPACK++ TEST COMPILE CMAKE_REQUIRED_FLAGS        ${CMAKE_REQUIRED_FLAGS}")
    message(STATUS "ARPACK++ TEST COMPILE CMAKE_REQUIRED_DEFINITIONS  ${CMAKE_REQUIRED_DEFINITIONS}")
    message(STATUS "ARPACK++ TEST COMPILE CMAKE_REQUIRED_LIBRARIES    ${CMAKE_REQUIRED_LIBRARIES}")
    message(STATUS "ARPACK++ TEST COMPILE CMAKE_REQUIRED_INCLUDES     ${CMAKE_REQUIRED_INCLUDES}")

    check_cxx_source_compiles("
        #include <complex>
        #if __has_include(<arpackpp/arcomp.h>)
        #include <arpackpp/arcomp.h>
        #include <arpackpp/ardscomp.h>
        #include <arpackpp/ardnsmat.h>
        #elif __has_include(<arpack++/arcomp.h>)
        #include <arpack++/arcomp.h>
        #include <arpack++/ardscomp.h>
        #include <arpack++/ardnsmat.h>
        #else
        #error Could not include arpack headers correctly
        #endif
        int main(){

            int                n = 4;
            int                nev = 2;
            int                ncv = 4;
            char ritz[10];
            std::string(\"LR\").copy(ritz,2);
            using Scalar = std::complex<double>;
            std::vector<Scalar> data = {Scalar(1.0,+0.0),Scalar(0.0,+1.0),Scalar(0.0,+1.0),Scalar(0.0,+1.0),
                                        Scalar(0.0,-1.0),Scalar(2.0,+0.0),Scalar(0.0,+2.0),Scalar(0.0,+2.0),
                                        Scalar(0.0,-1.0),Scalar(0.0,+2.0),Scalar(3.0,+3.0),Scalar(0.0,+3.0),
                                        Scalar(0.0,-1.0),Scalar(0.0,+2.0),Scalar(0.0,-3.0),Scalar(4.0,+0.0)};
            ARdsNonSymMatrix<Scalar, double> A(n, data.data());
            ARluCompStdEig<double> solver(nev, A,ritz,ncv);
            solver.FindEigenvectors();
            if (solver.EigenvaluesFound()){return 0;}
            else {return 1;}
        }
        " ARPACKPP_COMPILES_${TAG})
    if(ARPACKPP_COMPILES_${TAG})
        set(ARPACKPP_COMPILES_${TAG} TRUE PARENT_SCOPE)
    else()
        set(ARPACKPP_COMPILES_${TAG} FALSE PARENT_SCOPE)
    endif()
endfunction()



function(find_Arpackpp)
    if (NOT TARGET arpack)
        message(WARNING "Arpack++ needs to link to target [arpack]")
    endif()
    if (NOT TARGET lapacke)
        message(WARNING "Arpack++ needs to link to target [lapacke]")
    endif()

    if (NOT TARGET arpack++)
        message(STATUS "Searching for Arpack++ header-only in system")
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
                )
        if(ARPACKPP_INCLUDE_DIR)
            CheckArpackppCompiles("header_only" "-std=c++17 -g3"  ""  "" "${ARPACKPP_INCLUDE_DIR}" "arpack;lapacke")
            if(ARPACKPP_COMPILES_header_only)
                message(STATUS "Searching for Arpack++ headers - Success: ${ARPACKPP_INCLUDE_DIR}")
                add_library(arpack++ INTERFACE)
                target_link_libraries(arpack++ INTERFACE lapacke arpack blas lapack gfortran)
                target_include_directories(arpack++ SYSTEM INTERFACE ${ARPACKPP_INCLUDE_DIR})
            endif()
        else()
            message(STATUS "Searching for Arpack++ headers - failed - not found")
        endif()
    endif()


    if(NOT TARGET arpack++)
        message(STATUS "Searching for Arpack++ lib in system")
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
                )
        if(ARPACKPP_LIBRARIES AND ARPACKPP_INCLUDE_DIR)
            CheckArpackppCompiles("lib" "-std=c++17 -g3"  ""  "${ARPACKPP_LIBRARIES}" "${ARPACKPP_INCLUDE_DIR}" "arpack;lapacke")
            if(ARPACKPP_COMPILES_lib)
                message(STATUS "Searching for Arpack++ in system - Success: ${ARPACKPP_LIBRARIES}")
                message(STATUS "Note that old versions of Arpack++ (e.g. the default in Ubuntu Trusty 14.04 LTS) may fail to compile, they require '-fpermissive'.")
                add_library(arpack++ INTERFACE IMPORTED)
                target_link_libraries(arpack++ INTERFACE ${ARPACKPP_LIBRARIES} lapacke arpack blas lapack gfortran)
                target_include_directories(arpack++ SYSTEM INTERFACE ${ARPACKPP_INCLUDE_DIR})
            endif()
        else()
            message(STATUS "Searching for Arpack++ lib in system - failed - not found")
        endif()
    endif()
endfunction()