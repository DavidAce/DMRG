function(get_num_threads num_threads)
    cmake_host_system_information(RESULT _host_name   QUERY HOSTNAME)
    if(${_host_name} MATCHES "travis|TRAVIS|Travis")
        message("Using make -j 2 (travis-ci)")
        set(${num_threads} 2 PARENT_SCOPE)
        return()
    endif()
    if(CMAKE_BUILD_PARALLEL_LEVEL)
        message("Using make -j ${CMAKE_BUILD_PARALLEL_LEVEL}")
        set(${num_threads} ${CMAKE_BUILD_PARALLEL_LEVEL} PARENT_SCOPE)
        set(ENV{CMAKE_BUILD_PARALLEL_LEVEL} ${CMAKE_BUILD_PARALLEL_LEVEL})
        return()
    endif()

    if(${_host_name} MATCHES "ThinkStation")
        message("Using make -j 8")
        set(${num_threads} 8 PARENT_SCOPE)
        return()
    endif()
    if(${_host_name} MATCHES "raken")
        message("NUM_THREADS : 32")
        set(${num_threads} 32 PARENT_SCOPE)
        return()
    endif()
    if(DEFINED ENV{MAKEFLAGS})
        string(REGEX MATCH "-j[ ]*[0-9]+|--jobs[ ]*[0-9]+" REGMATCH "$ENV{MAKEFLAGS}")
        string(REGEX MATCH "[0-9]+" NUM "${REGMATCH}")
        if(NUM)
            message("NUM_THREADS : ${NUM}")
            set(${num_threads} ${NUM} PARENT_SCOPE)
            return()
        endif()
    endif()
    if(MAKE_THREADS)
        string(REGEX MATCH "-j[ ]*[0-9]+|--jobs[ ]*[0-9]+" REGMATCH "${MAKE_THREADS}")
        string(REGEX MATCH "[0-9]+" NUM "${REGMATCH}")
        if(NUM)
            message("Using make -j ${NUM}")
            set(${num_threads} ${NUM} PARENT_SCOPE)
            return()
        endif()
    endif()
    include(ProcessorCount)
    ProcessorCount(NUM)
    if(NUM)
        message("Using make -j ${NUM}")
        set(${num_threads} ${NUM} PARENT_SCOPE)
        return()
    endif()


    set(${num_threads} 1 PARENT_SCOPE)

endfunction()