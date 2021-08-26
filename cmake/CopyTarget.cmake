function(copy_target target_old target_new)
    if(NOT TARGET ${target_old})
        message(DEBUG "There is no target named '${target_old}'")
        return()
    endif()
    if(TARGET ${target_new})
        message(DEBUG "There is already a target named '${target_new}'")
        return()
    endif()


    execute_process(COMMAND cmake --help-property-list OUTPUT_VARIABLE CMAKE_PROPERTY_LIST)
    list(APPEND CMAKE_PROPERTY_LIST IMPORTED_LOCATION IMPORTED_LOCATION_RELEASE IMPORTED_LOCATION_DEBUG)
    # Convert command output into a CMake list
    STRING(REGEX REPLACE ";" "\\\\;" CMAKE_PROPERTY_LIST "${CMAKE_PROPERTY_LIST}")
    STRING(REGEX REPLACE "\n" ";" CMAKE_PROPERTY_LIST "${CMAKE_PROPERTY_LIST}")
    # For some reason, "TYPE" shows up twice - others might too?
    list(REMOVE_DUPLICATES CMAKE_PROPERTY_LIST)

    # build whitelist by filtering down from CMAKE_PROPERTY_LIST in case cmake is
    # a different version, and one of our hardcoded whitelisted properties
    # doesn't exist!
    unset(CMAKE_WHITELISTED_PROPERTY_LIST)
    foreach(prop ${CMAKE_PROPERTY_LIST})
        if(prop MATCHES "^(INTERFACE|[_a-z]|IMPORTED_LIBNAME_|MAP_IMPORTED_CONFIG_)|^(COMPATIBLE_INTERFACE_(BOOL|NUMBER_MAX|NUMBER_MIN|STRING)|EXPORT_NAME|IMPORTED(_GLOBAL|_CONFIGURATIONS|_LIBNAME)?|NAME|TYPE|NO_SYSTEM_FROM_IMPORTED)$")
            list(APPEND CMAKE_WHITELISTED_PROPERTY_LIST ${prop})
        endif()
    endforeach()

    get_target_property(target_type     ${target_old} TYPE)
    get_target_property(target_imported ${target_old} IMPORTED)
    get_target_property(target_global   ${target_old} IMPORTED_GLOBAL)
    if(target_global)
        set(GLOBAL GLOBAL)
    endif()

    if(target_type MATCHES "INTERFACE" AND target_imported)
        add_library(${target_new} INTERFACE IMPORTED ${GLOBAL})
    elseif(target_type MATCHES "INTERFACE" AND NOT target_imported)
        add_library(${target_new} INTERFACE ${GLOBAL})
    elseif(target_type MATCHES "STATIC")
        add_library(${target_new} STATIC IMPORTED ${GLOBAL})
    elseif(target_type MATCHES "SHARED")
        add_library(${target_new} SHARED IMPORTED ${GLOBAL})
    else()
        message(FATAL_ERROR "Copy target not supported for type: ${target_type}" )
    endif()

    if(target_type MATCHES "INTERFACE")
        set(PROP_LIST ${CMAKE_WHITELISTED_PROPERTY_LIST})
    else()
        set(PROP_LIST ${CMAKE_PROPERTY_LIST})
    endif()

    foreach (prop ${PROP_LIST})
        # Skip read-only properties
        if("${prop}" MATCHES "IMPORTED_GLOBAL|NAME|TYPE")
            continue()
        endif()
        string(REPLACE "<CONFIG>" "${CMAKE_BUILD_TYPE}" prop ${prop})
        # message ("Checking ${prop}")
        get_property(propval TARGET ${target_old} PROPERTY ${prop} SET) # Returns a boolean if set
        if (propval)
            get_target_property(propval ${target_old} ${prop})
#            message ("${target_old} ${prop} = ${propval}")
            set_target_properties(${target_new} PROPERTIES "${prop}" "${propval}")
            # Modernize
            if("${prop}" MATCHES "^LOCATION")
                set_target_properties(${target_new} PROPERTIES IMPORTED_LOCATION "${propval}")
            endif()
        endif()


    endforeach()

endfunction()