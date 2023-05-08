function(find_quadmath)

    find_library(QUADMATH_LIBRARY NAMES quadmath)
    find_path(QUADMATH_INCLUDE_DIR NAMES quadmath.h PATH_SUFFIXES include)

    include(CheckTypeSize)
    check_type_size(__float128 FLOAT128_EXISTS BUILTIN_TYPES_ONLY LANGUAGE CXX)

    if(NOT QUADMATH_LIBRARY OR NOT QUADMATH_INCLUDE_DIR)
        include(CMakePushCheckState)
        include(CheckCXXSourceCompiles)
        cmake_push_check_state(RESET)
        list(APPEND CMAKE_REQUIRED_LIBRARIES "quadmath")
        check_cxx_source_compiles("
            #include <quadmath.h>
            int main(void){
                __float128 foo = ::sqrtq(123.456);
            }"
            QUADMATH_LINK_ONLY
        )
        cmake_pop_check_state()
        if (QUADMATH_LINK_ONLY)
            set(QUADMATH_INCLUDE_DIR "unused" CACHE PATH "" FORCE)
            set(QUADMATH_LIBRARY "quadmath" CACHE FILEPATH "" FORCE)
        elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
            message(STATUS "Failed to compile a simple quadmath program: Note that Clang does not support quadmath")
        endif()
    endif()



endfunction()

find_quadmath()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(quadmath DEFAULT_MSG QUADMATH_LIBRARY QUADMATH_INCLUDE_DIR FLOAT128_EXISTS)

if(quadmath_FOUND)
    if(NOT TARGET quadmath::quadmath)
        if(QUADMATH_LINK_ONLY)
            add_library(quadmath::quadmath INTERFACE IMPORTED)
            target_link_libraries(quadmath::quadmath INTERFACE ${QUADMATH_LIBRARY})
        else()
            add_library(quadmath::quadmath UNKNOWN IMPORTED)
            set_target_properties(quadmath::quadmath PROPERTIES IMPORTED_LOCATION "${QUADMATH_LIBRARY}")
            if(QUADMATH_INCLUDE_DIR AND NOT QUADMATH_INCLUDE_DIR MATCHES "unused")
                target_include_directories(quadmath::quadmath SYSTEM INTERFACE ${QUADMATH_INCLUDE_DIR} )
            endif()
        endif()
        message(DEBUG "Defined target quadmath::quadmath for library: ${QUADMATH_LIBRARY}")
    endif()
endif()