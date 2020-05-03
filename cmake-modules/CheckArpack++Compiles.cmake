function(check_arpackpp_compiles TARGETS LIBS INCS OPTS DEFS)
    if(NOT BUILD_SHARED_LIBS)
        list(APPEND CMAKE_REQUIRED_LIBRARIES -static)
    endif()
    list(APPEND CMAKE_REQUIRED_LIBRARIES     ${LIBS} ${TARGETS})
    list(APPEND CMAKE_REQUIRED_INCLUDES      ${INCS})
    list(APPEND CMAKE_REQUIRED_FLAGS         ${OPTS} -std=c++17)
    list(APPEND CMAKE_REQUIRED_DEFINITIONS   ${DEFS})


    list(TRANSFORM "CMAKE_REQUIRED_DEFINITIONS" PREPEND "-D")  # Definitions should start with "-D"
    string(REPLACE ";" " "  CMAKE_REQUIRED_FLAGS          "${CMAKE_REQUIRED_FLAGS}")        # Needs to be a space-separated list

    include(CheckCXXSourceCompiles)
    if(DMRG_PRINT_CHECKS OR NOT arpack++_FIND_VERBOSE)
        message(STATUS "ARPACK++ COMPILE TEST CMAKE_REQUIRED_LIBRARIES    ${CMAKE_REQUIRED_LIBRARIES}")
        message(STATUS "ARPACK++ COMPILE TEST CMAKE_REQUIRED_INCLUDES     ${CMAKE_REQUIRED_INCLUDES}")
        message(STATUS "ARPACK++ COMPILE TEST CMAKE_REQUIRED_FLAGS        ${CMAKE_REQUIRED_FLAGS}")
        message(STATUS "ARPACK++ COMPILE TEST CMAKE_REQUIRED_DEFINITIONS  ${CMAKE_REQUIRED_DEFINITIONS}")
    endif()
    #   Test features
    check_cxx_source_compiles("
        #include <complex>
        #include <vector>
        #if __has_include(<arpackpp/arcomp.h>)
            #include <arpackpp/arcomp.h>
            #include <arpackpp/ardscomp.h>
            #include <arpackpp/ardnsmat.h>
        #elif __has_include(<arpack++/arcomp.h>)
            #include <arpack++/arcomp.h>
            #include <arpack++/ardscomp.h>
            #include <arpack++/ardnsmat.h>
        #elif __has_include(<arcomp.h>)
            #include <arcomp.h>
            #include <ardscomp.h>
            #include <ardnsmat.h>
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
        " ARPACKPP_COMPILES)
    set(ARPACKPP_COMPILES ${ARPACKPP_COMPILES} PARENT_SCOPE)
    if(NOT ARPACKPP_COMPILES)
        unset(ARPACKPP_COMPILES CACHE)
        unset(ARPACKPP_COMPILES PARENT_SCOPE)
        if(DMRG_PRINT_CHECKS AND EXISTS "${CMAKE_BINARY_DIR}/CMakeFiles/CMakeError.log")
            file(READ "${CMAKE_BINARY_DIR}/CMakeFiles/CMakeError.log" ERROR_LOG)
            message(STATUS "CMakeError.log:\n ${ERROR_LOG}")
        endif()
    endif()
endfunction()
