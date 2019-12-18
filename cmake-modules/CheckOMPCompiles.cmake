function(strip_genex input_string output_string)
    set(result_string)
    foreach(elem ${input_string})
        if(${elem} MATCHES "[$]")
            string(REGEX MATCHALL "([-=]+[a-z0-9]+)" filtered_string "${input_string}")
            list(APPEND result_string ${filtered_string})
        else()
            list(APPEND result_string ${elem})
        endif()
    endforeach()
    set(${output_string} ${result_string} PARENT_SCOPE)
endfunction()


function(check_omp_compiles REQUIRED_FLAGS REQUIRED_LIBRARIES_UNPARSED REQUIRED_INCLUDES)
    include(CheckIncludeFileCXX)
    include(cmake-modules/getExpandedTarget.cmake)
    expand_target_libs("${REQUIRED_LIBRARIES_UNPARSED}" expanded_libs)
    expand_target_incs("${REQUIRED_LIBRARIES_UNPARSED}" expanded_incs)
    expand_target_opts("${REQUIRED_LIBRARIES_UNPARSED}" expanded_opts)

    strip_genex("${expanded_libs}"     expanded_libs)
    strip_genex("${expanded_incs}"     expanded_incs)
    strip_genex("${expanded_opts}"     expanded_opts)
    strip_genex("${REQUIRED_INCLUDES}" REQUIRED_INCLUDES)
    strip_genex("${REQUIRED_FLAGS}"    REQUIRED_FLAGS)

    set(CMAKE_REQUIRED_LIBRARIES "${expanded_libs}") # Can be a ;list
    set(CMAKE_REQUIRED_INCLUDES  "${REQUIRED_INCLUDES};${expanded_incs}") # Can be a ;list
    string(REPLACE ";" " " CMAKE_REQUIRED_FLAGS "${REQUIRED_FLAGS} ${expanded_opts}") # Needs to be a space-separated list

    message(STATUS "OPENMP TEST COMPILE CMAKE_REQUIRED_FLAGS        ${CMAKE_REQUIRED_FLAGS}")
    message(STATUS "OPENMP TEST COMPILE CMAKE_REQUIRED_DEFINITIONS  ${CMAKE_REQUIRED_DEFINITIONS}")
    message(STATUS "OPENMP TEST COMPILE CMAKE_REQUIRED_LIBRARIES    ${CMAKE_REQUIRED_LIBRARIES}")
    message(STATUS "OPENMP TEST COMPILE CMAKE_REQUIRED_INCLUDES     ${CMAKE_REQUIRED_INCLUDES}")

    unset(has_omp_h)
    unset(has_omp_h CACHE)
    unset(OMP_COMPILES CACHE)
    unset(OMP_COMPILES)

    check_include_file_cxx(omp.h    has_omp_h)
    include(CheckCXXSourceCompiles)
    check_cxx_source_compiles("
            #include <omp.h>
            #include <iostream>
            #ifndef _OPENMP
            #error You forgot to define _OPENMP
            #endif
            int main() {
                omp_set_num_threads(4);
                std::cout << \"OMP Threads \" << omp_get_max_threads() << std::endl;
                int counter = 0;
                #pragma omp parallel for shared(counter)
                for(int i = 0; i < 16; i++){
                    #pragma omp atomic
                    counter++;
                }
                std::cout << \"Counter is: \" << counter << std::endl;

                return 0;
            }
            " OMP_COMPILES
            )
    if(NOT OMP_COMPILES)
        set(OMP_COMPILES FALSE PARENT_SCOPE)
    else()
        set(OMP_COMPILES TRUE PARENT_SCOPE)
    endif()
endfunction()