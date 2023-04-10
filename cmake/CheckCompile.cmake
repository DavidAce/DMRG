cmake_minimum_required(VERSION 3.17)

function(check_compile pkg tgt)
    set(options REQUIRED)
    cmake_parse_arguments(PARSE_ARGV 1 CHECK "${options}" "${oneValueArgs}" "${multiValueArgs}")


    if(NOT ${pkg}_compiles_${tgt})
        #list(APPEND CMAKE_REQUIRED_LIBRARIES ${tgt})
        message(CHECK_START "Test compile -- ${pkg} [${tgt}]")
        try_compile(${pkg}_compiles_${tgt}
                ${CMAKE_BINARY_DIR}
                ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/compile/${pkg}.cpp
                OUTPUT_VARIABLE compile_out
                LINK_LIBRARIES ${tgt}
                CXX_STANDARD 17
                CXX_EXTENSIONS OFF
                )
        if(${pkg}_compiles_${tgt})
            message(CHECK_PASS "Success")
            file(APPEND ${CMAKE_BINARY_DIR}/CMakeFiles/CMakeOutput.log "${compile_out}")
            set(${pkg}_compiles_${tgt} "${${pkg}_compiles_${tgt}}" CACHE BOOL "" FORCE)
            mark_as_advanced(${pkg}_compiles_${tgt})
        else()
            message(CHECK_FAIL "Failed")
            file(APPEND ${CMAKE_BINARY_DIR}/CMakeFiles/CMakeError.log "${compile_out}")
            include(${CMAKE_CURRENT_FUNCTION_LIST_DIR}/PrintTargetInfo.cmake)
            print_target_info_recursive(${tgt})
            if(CHECK_REQUIRED)
                message(FATAL_ERROR "Failed to compile ${pkg} with target [${tgt}]")
            endif()
        endif()
    endif()
endfunction()
