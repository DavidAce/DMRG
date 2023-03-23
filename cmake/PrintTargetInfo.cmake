cmake_minimum_required(VERSION 3.18)

if(NOT WIN32)
    if("$ENV{CLICOLOR}" OR "$ENV{CLICOLOR_FORCE}")
        string(ASCII 27 Esc)
        set(ColourReset "${Esc}[m")
        set(ColourBold "${Esc}[1m")
        set(Red "${Esc}[31m")
        set(Green "${Esc}[32m")
        set(Yellow "${Esc}[33m")
        set(Blue "${Esc}[34m")
        set(Magenta "${Esc}[35m")
        set(Cyan "${Esc}[36m")
        set(White "${Esc}[37m")
        set(BoldRed "${Esc}[1;31m")
        set(BoldGreen "${Esc}[1;32m")
        set(BoldYellow "${Esc}[1;33m")
        set(BoldBlue "${Esc}[1;34m")
        set(BoldMagenta "${Esc}[1;35m")
        set(BoldCyan "${Esc}[1;36m")
        set(BoldWhite "${Esc}[1;37m")
    endif()
endif()

function(pad_string OUT_VARIABLE DESIRED_LENGTH FILL_CHAR VALUE)
    string(LENGTH "${VALUE}" VALUE_LENGTH)
    math(EXPR REQUIRED_PADS "${DESIRED_LENGTH} - ${VALUE_LENGTH}")
    set(PAD ${VALUE})
    if(REQUIRED_PADS GREATER 0)
        math(EXPR REQUIRED_MINUS_ONE "${REQUIRED_PADS} - 1")
        foreach(FOO RANGE ${REQUIRED_MINUS_ONE})
            set(PAD "${PAD}${FILL_CHAR}")
        endforeach()
    endif()
    set(${OUT_VARIABLE} "${ColourBold}${PAD}${ColourReset}" PARENT_SCOPE)
endfunction()

function(print_target_info target_name prefix)
    if(TARGET ${target_name})
        get_target_property(PROP_INC ${target_name} INTERFACE_INCLUDE_DIRECTORIES)
        get_target_property(PROP_LIB ${target_name} INTERFACE_LINK_LIBRARIES)
        get_target_property(PROP_OPT ${target_name} INTERFACE_COMPILE_OPTIONS)
        get_target_property(PROP_DEF ${target_name} INTERFACE_COMPILE_DEFINITIONS)
        get_target_property(PROP_FTR ${target_name} INTERFACE_COMPILE_FEATURES)
        get_target_property(PROP_TYP ${target_name} TYPE)
        get_target_property(PROP_IMP ${target_name} IMPORTED)
        get_imported_location(PROP_LOC ${target_name})

        # Remove genexp patterns such as $<$<CONFIG:Release>:> or $<LINK_ONLY:tb-flags>
        list(TRANSFORM PROP_INC REPLACE "[\$<]+[A-Za-z0-9_]+[:]|[A-Za-z0-9_]+[>]+[,:]|>" "")
        list(TRANSFORM PROP_LIB REPLACE "[\$<]+[A-Za-z0-9_]+[:]|[A-Za-z0-9_]+[>]+[,:]|>" "")
        list(TRANSFORM PROP_OPT REPLACE "[\$<]+[A-Za-z0-9_]+[:]|[A-Za-z0-9_]+[>]+[,:]|>" "")
        list(TRANSFORM PROP_DEF REPLACE "[\$<]+[A-Za-z0-9_]+[:]|[A-Za-z0-9_]+[>]+[,:]|>" "")
        list(TRANSFORM PROP_FTR REPLACE "[\$<]+[A-Za-z0-9_]+[:]|[A-Za-z0-9_]+[>]+[,:]|>" "")
        # strip leading semicolons
        string(REGEX REPLACE "^;" "" PROP_INC "${PROP_INC}")
        string(REGEX REPLACE "^;" "" PROP_LIB "${PROP_LIB}")
        string(REGEX REPLACE "^;" "" PROP_OPT "${PROP_OPT}")
        string(REGEX REPLACE "^;" "" PROP_DEF "${PROP_DEF}")
        string(REGEX REPLACE "^;" "" PROP_FTR "${PROP_FTR}")

        list(REMOVE_DUPLICATES PROP_INC)
        list(REMOVE_DUPLICATES PROP_LIB)
        list(REMOVE_DUPLICATES PROP_OPT)
        list(REMOVE_DUPLICATES PROP_DEF)
        list(REMOVE_DUPLICATES PROP_FTR)

        pad_string(padded_target "60" " " "${prefix}[${target_name}]")
        list(APPEND CMAKE_MESSAGE_INDENT "${padded_target}")
        if(PROP_LIB)
            message(STATUS "LIBRARY : ${PROP_LIB}")
        endif()
        if(PROP_INC)
            message(STATUS "INCLUDE : ${PROP_INC}")
        endif()
        if(PROP_OPT)
            message(STATUS "OPTIONS : ${PROP_OPT}")
        endif()
        if(PROP_DEF)
            message(STATUS "DEFINES : ${PROP_DEF}")
        endif()
        if(PROP_FTR)
            message(STATUS "FEATURE : ${PROP_FTR}")
        endif()
        if(PROP_LOC)
            message(STATUS "IMPORTS : ${PROP_LOC}")
        endif()
        list(POP_BACK CMAKE_MESSAGE_INDENT)
    endif()
endfunction()

function(string_target_info target_name prefix strout)
    if(TARGET ${target_name})
        get_target_property(PROP_INC ${target_name} INTERFACE_INCLUDE_DIRECTORIES)
        get_target_property(PROP_LIB ${target_name} INTERFACE_LINK_LIBRARIES)
        get_target_property(PROP_OPT ${target_name} INTERFACE_COMPILE_OPTIONS)
        get_target_property(PROP_DEF ${target_name} INTERFACE_COMPILE_DEFINITIONS)
        get_target_property(PROP_FTR ${target_name} INTERFACE_COMPILE_FEATURES)
        get_target_property(PROP_TYP ${target_name} TYPE)
        get_target_property(PROP_IMP ${target_name} IMPORTED)
        get_imported_location(PROP_LOC ${target_name})

        pad_string(padded_target "60" " " "${prefix}[${target_name}]")
        if(PROP_LIB)
            string(APPEND str "${padded_target} LIBRARY : ${PROP_LIB}\n")
        endif()
        if(PROP_INC)
            string(APPEND str "${padded_target} INCLUDE : ${PROP_INC}\n")
        endif()
        if(PROP_OPT)
            string(APPEND str "${padded_target} OPTIONS : ${PROP_OPT}\n")
        endif()
        if(PROP_DEF)
            string(APPEND str "${padded_target} DEFINES : ${PROP_DEF}\n")
        endif()
        if(PROP_FTR)
            string(APPEND str"${padded_target} FEATURE : ${PROP_FTR}\n")
        endif()
        if(PROP_LOC)
            string(APPEND str "${padded_target} IMPORTS : ${PROP_LOC}\n")
        endif()
    endif()
    set(${strout} "${${strout}}${str}" PARENT_SCOPE)
endfunction()

function(print_compiler_info prefix)
    pad_string(padded_string "60" " " "${prefix}[${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}]")
    list(APPEND CMAKE_MESSAGE_INDENT "${padded_string}")
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
        message(STATUS "OPTIONS : ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
    endif()
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        message(STATUS "OPTIONS : ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
    endif()
    if(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
        message(STATUS "OPTIONS : ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
    endif()
    if(CMAKE_BUILD_TYPE STREQUAL "MinSizeRel")
        message(STATUS "OPTIONS : ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_MINSIZEREL}")
    endif()
    list(POP_BACK CMAKE_MESSAGE_INDENT)
endfunction()

function(string_compiler_info prefix strout)
    pad_string(padded_string "60" " " "${prefix}[${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}]")
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
        set(padded_string "${padded_string} OPTIONS : ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}\n")
    endif()
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        set(padded_string "${padded_string} OPTIONS : ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}\n")
    endif()
    if(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
        set(APPEND padded_string "${padded_string} OPTIONS : ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}\n")
    endif()
    if(CMAKE_BUILD_TYPE STREQUAL "MinSizeRel")
        set(APPEND padded_string "${padded_string} OPTIONS : ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_MINSIZEREL}\n")
    endif()
    set(${strout} "${${strout}}${padded_string}" PARENT_SCOPE)
endfunction()

# Print summary of project targets
function(print_project_summary prj)
    include(${CMAKE_CURRENT_FUNCTION_LIST_DIR}/getExpandedTarget.cmake)
    if(NOT TARGET ${prj})
        message(FATAL_ERROR "${prj} is not a valid target")
    endif()
    message(STATUS "| PROJECT SUMMARY [${prj}]")
    message(STATUS "|--------------------------")
    print_compiler_info("| ")
    expand_target_all_targets(${prj} TARGET_EXPANDED "")
    foreach(t ${TARGET_EXPANDED})
        print_target_info(${t} "| ")
    endforeach()
    message(STATUS "|--------------------------")
endfunction()

# Print summary of project targets
function(write_project_summary prj)
    include(cmake/getExpandedTarget.cmake)
    if(NOT TARGET ${prj})
        message(FATAL_ERROR "${prj} is not a valid target")
    endif()
    string_compiler_info("| " target_strinfo)
    expand_target_all_targets(${prj} TARGET_EXPANDED)
    foreach(t ${TARGET_EXPANDED})
        string_target_info(${t} "| " target_strinfo)
    endforeach()
    string(REPLACE "<LINK_ONLY:" "<1:" target_strinfo "${target_strinfo}")
    file(GENERATE OUTPUT ${CMAKE_BINARY_DIR}/${PROJECT_NAME}_$<CONFIG>_$<COMPILE_LANGUAGE>_target_info.txt
         CONTENT "${target_strinfo}"
         TARGET ${PROJECT_NAME})
endfunction()

# Print summary of project targets
function(print_and_write_project_summary prj)
    include(cmake/getExpandedTarget.cmake)
    if(NOT TARGET ${prj})
        message(FATAL_ERROR "${prj} is not a valid target")
    endif()
    message(STATUS "| PROJECT SUMMARY [${prj}]")
    message(STATUS "|--------------------------")
    print_compiler_info("| ")
    string_compiler_info("| " target_strinfo)
    expand_target_all_targets(${prj} TARGET_EXPANDED "")
    foreach(t ${TARGET_EXPANDED})
        print_target_info(${t} "| ")
        string_target_info(${t} "| " target_strinfo)
    endforeach()
    message(STATUS "|--------------------------")
    string(REPLACE "<LINK_ONLY:" "<1:" target_strinfo "${target_strinfo}")
    file(GENERATE OUTPUT ${CMAKE_BINARY_DIR}/${PROJECT_NAME}_$<CONFIG>_$<COMPILE_LANGUAGE>_target_info.txt
         CONTENT "${target_strinfo}"
         TARGET ${PROJECT_NAME})
endfunction()