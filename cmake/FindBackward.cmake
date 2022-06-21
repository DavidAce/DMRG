unset(BACKWARD_LIBRARY)
unset(BACKWARD_LIBRARY CACHE)

find_path(BACKWARD_INCLUDE_DIR
        backward.hpp
        HINTS ${DMRG_DEPS_INSTALL_DIR}
        PATH_SUFFIXES include backward/include
        NO_CMAKE_ENVIRONMENT_PATH
        NO_SYSTEM_ENVIRONMENT_PATH
        NO_CMAKE_SYSTEM_PATH
        )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Backward
        DEFAULT_MSG
        BACKWARD_INCLUDE_DIR)


if (Backward_FOUND AND NOT TARGET Backward::Backward)
    add_library(Backward::Backward INTERFACE IMPORTED)
    target_include_directories(Backward::Backward SYSTEM INTERFACE ${BACKWARD_INCLUDE_DIR})

    # Configure
    find_package(libunwind)


    if (NOT TARGET Backward::libunwind)
        if (TARGET libunwind::libunwind)
            add_library(Backward::libunwind IMPORTED INTERFACE)
            target_link_libraries(Backward::libunwind INTERFACE libunwind::libunwind)
        else ()
            find_library(LIBUNWIND_LIBRARY unwind)
            find_path(LIBUNWIND_INCLUDE_DIR libunwind.h)
            if (LIBUNWIND_LIBRARY AND LIBUNWIND_INCLUDE_DIR)
                add_library(Backward::libunwind IMPORTED INTERFACE)
                set_target_properties(Backward::libunwind PROPERTIES IMPORTED_LOCATION "${LIBUNWIND_LIBRARY}")
                target_include_directories(Backward::libunwind SYSTEM INTERFACE ${LIBUNWIND_INCLUDE_DIR})
            endif ()
        endif ()
    endif ()
    if (TARGET Backward::libunwind)
        target_link_libraries(Backward::Backward INTERFACE Backward::libunwind)
        target_compile_definitions(Backward::Backward INTERFACE
                BACKWARD_HAS_LIBUNWIND=1
                BACKWARD_HAS_UNWIND=0
                BACKWARD_HAS_BACKTRACE=0 # Not considered
                )
    else ()
        target_compile_definitions(Backward::Backward INTERFACE
                BACKWARD_HAS_LIBUNWIND=0
                BACKWARD_HAS_UNWIND=1
                BACKWARD_HAS_BACKTRACE=0 # Not considered
                )
    endif ()


    if (NOT TARGET Backward::libdw)
        # find libdw
        find_library(LIBDW_LIBRARY dw)
        find_path(LIBDW_INCLUDE_DIR NAMES "elfutils/libdw.h" "elfutils/libdwfl.h")
        if (LIBDW_LIBRARY AND LIBDW_INCLUDE_DIR)
            add_library(Backward::libdw UNKNOWN IMPORTED)
            set_target_properties(Backward::libdw PROPERTIES IMPORTED_LOCATION "${LIBDW_LIBRARY}")
            target_include_directories(Backward::libdw SYSTEM INTERFACE ${LIBDW_INCLUDE_DIR})
            target_link_libraries(Backward::Backward INTERFACE Backward::libdw)
            target_compile_definitions(Backward::Backward INTERFACE
                    BACKWARD_HAS_DW=1
                    BACKWARD_HAS_DWARF=0
                    BACKWARD_HAS_BFD=0 # Not considered
                    BACKWARD_HAS_BACKTRACE_SYMBOL=0
                    )
        endif ()
    endif ()
    if (NOT TARGET Backward::libdw AND NOT TARGET Backward::libdwarf)
        # find libdwarf
        find_path(LIBDWARF_INCLUDE_DIR NAMES "libdwarf.h" PATH_SUFFIXES libdwarf)
        find_path(LIBELF_INCLUDE_DIR NAMES "libelf.h")
        find_path(LIBDL_INCLUDE_DIR NAMES "dlfcn.h")
        find_library(LIBDWARF_LIBRARY dwarf)
        find_library(LIBELF_LIBRARY elf)
        find_library(LIBDL_LIBRARY dl)
        if (LIBDWARF_LIBRARY AND LIBDWARF_INCLUDE_DIR AND LIBELF_LIBRARY AND LIBELF_INCLUDE_DIR AND LIBDL_LIBRARY AND LIBDL_INCLUDE_DIR)
            add_library(Backward::libdwarf INTERFACE IMPORTED)
            add_library(Backward::libelf INTERFACE IMPORTED)
            add_library(Backward::libdl INTERFACE IMPORTED)
            set_target_properties(Backward::libdwarf PROPERTIES IMPORTED_LOCATION "${LIBDWARF_LIBRARY}")
            set_target_properties(Backward::libelf PROPERTIES IMPORTED_LOCATION "${LIBELF_LIBRARY}")
            set_target_properties(Backward::libdl PROPERTIES IMPORTED_LOCATION "${LIBDL_LIBRARY}")
            target_include_directories(Backward::libdwarf SYSTEM INTERFACE ${LIBDWARF_INCLUDE_DIR})
            target_include_directories(Backward::libelf SYSTEM INTERFACE ${LIBELF_INCLUDE_DIR})
            target_include_directories(Backward::libdl SYSTEM INTERFACE ${LIBDL_INCLUDE_DIR})
            target_link_libraries(Backward::libelf INTERFACE Backward::libdl)
            target_link_libraries(Backward::libdwarf INTERFACE Backward::libelf Backward::libdl)
            target_link_libraries(Backward::Backward INTERFACE Backward::libdwarf)
            target_compile_definitions(Backward::Backward INTERFACE
                    BACKWARD_HAS_DW=0
                    BACKWARD_HAS_DWARF=1
                    BACKWARD_HAS_BFD=0 # Not considered
                    BACKWARD_HAS_BACKTRACE_SYMBOL=0
                    )
        endif ()
    endif ()
    if (NOT TARGET Backward::libdw AND NOT TARGET Backward::libdwarf)
        target_compile_definitions(Backward::Backward INTERFACE
                BACKWARD_HAS_DW=0
                BACKWARD_HAS_DWARF=0
                BACKWARD_HAS_BFD=0 # Not considered
                BACKWARD_HAS_BACKTRACE_SYMBOL=1
                )
    endif ()

endif ()