
configure_file(
        ${CMAKE_SOURCE_DIR}/cmake_modules/directory.h.in
        ${CMAKE_BINARY_DIR}/directory/directory.h
)

include_directories(${CMAKE_BINARY_DIR}/directory)

