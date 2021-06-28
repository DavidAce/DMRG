cmake_minimum_required(VERSION 3.15)
function(try_compile_std_optional)
    if(NOT has_std_optional)
        file(WRITE ${CMAKE_BINARY_DIR}/CMakeFiles/CMakeTmp/check-std-optional.cpp
        "
        // Include optional or experimental/optional
        namespace check{}
        #if __has_include(<optional>)
        #include <optional>
        #elif __has_include(<experimental/optional>)
        #include <experimental/optional>
        namespace check{
            namespace std{
                constexpr const std::experimental::nullopt_t &nullopt = std::experimental::nullopt ;
                template<typename T> using std::optional = std::experimental::optional<T>;
            }
        }
        #else
            #error Could not find <optional> or <experimental/optional>
        #endif


        int main(){
            using namespace check;
            std::optional<int> optVar = std::nullopt;
            return 0;
        }
        ")
        message(STATUS "Performing Test has_std_optional")
        try_compile(has_std_optional
                ${CMAKE_BINARY_DIR}
                ${CMAKE_BINARY_DIR}/CMakeFiles/CMakeTmp/check-std-optional.cpp
                OUTPUT_VARIABLE OPTIONAL_COMPILES_OUT
                CXX_STANDARD 17
                CXX_EXTENSIONS OFF
                )
        if(has_std_optional)
            file(APPEND ${CMAKE_BINARY_DIR}/CMakeFiles/CMakeOutput.log "${OPTIONAL_COMPILES_OUT}")
            message(STATUS "Performing Test has_std_optional - Success")
        else()
            file(APPEND ${CMAKE_BINARY_DIR}/CMakeFiles/CMakeError.log "${OPTIONAL_COMPILES_OUT}")
            message(STATUS "Performing Test has_std_optional - Failed")
        endif()
        set(has_std_optional ${has_std_optional} CACHE BOOL "")
        mark_as_advanced(has_std_optional)
    endif()
endfunction()

function(CheckCXXOptional)
    cmake_policy(SET CMP0075 NEW)
    include(CheckIncludeFileCXX)
    try_compile_std_optional()
    if(NOT has_std_optional)
        check_include_file_cxx(optional    has_optional)
        check_include_file_cxx(experimental/optional has_experimental_optional )
        if(NOT has_optional AND NOT has_experimental_optional)
            message(FATAL_ERROR "
                Missing both <optional> and <experimental/optional> headers.\n\
                Consider using a newer compiler (GCC 6 or above, Clang 5 or above),\n\
                or checking the compiler flags. If using Clang, try passing the variable \n\
                GCC_TOOLCHAIN=<path> \n\
                where path is the install directory of a recent GCC to use its standard library.
        ")
        endif()
    endif()
endfunction()