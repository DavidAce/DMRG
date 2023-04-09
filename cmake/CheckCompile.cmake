cmake_minimum_required(VERSION 3.17)

function(check_compile pkg tgt file)
    if(NOT ${pkg}_compiles_${tgt})
        list(APPEND CMAKE_REQUIRED_LIBRARIES ${tgt})
        message(CHECK_START "Test compile -- ${pkg} [${tgt}]")
        try_compile(${pkg}_compiles_${tgt}
                ${CMAKE_BINARY_DIR}
                ${file}
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
        endif()
    endif()
endfunction()
