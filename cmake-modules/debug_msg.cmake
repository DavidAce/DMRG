set(COUT_STRING "")
if(${CMAKE_BUILD_TYPE} MATCHES Debug)
    set(COUT_STRING "std::cout << str << std::endl;")
endif()

configure_file(
        ${CMAKE_SOURCE_DIR}/cmake-modules/debug_msg.h.in
        ${CMAKE_BINARY_DIR}/debug_msg/debug_msg.h
)

include_directories(${CMAKE_BINARY_DIR}/debug_msg)
