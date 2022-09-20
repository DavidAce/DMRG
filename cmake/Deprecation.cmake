
if(DMRG_DOWNLOAD_METHOD)
    message(FATAL_ERROR "The variable [DMRG_DOWNLOAD_METHOD] has been deprecated. Replace by:\n"
            "DMRG_PACKAGE_MANAGER:STRING=[find|cmake|fetch|find-or-cmake|find-or-fetch|conan]")
endif()

if(DMRG_DEPS_IN_SUBDIR)
    message(DEPRECATION "The option [DMRG_DEPS_IN_SUBDIR] will be ignored.")
endif()

if(DMRG_PRINT_INFO)
    message(DEPRECATION "The option [DMRG_PRINT_INFO] has been deprecated. Replace by:\n"
            "CMake CLI option --loglevel=[TRACE|DEBUG|VERBOSE|STATUS...] or\n"
            "set CMAKE_MESSAGE_LOG_LEVEL=[TRACE|DEBUG|VERBOSE|STATUS...]"
            )
endif()

if(DMRG_PROFILE_BUILD)
    message(FATAL_ERROR "The option [DMRG_PROFILE_BUILD] has been deprecated. Replace by:\n"
            "COMPILER_PROFILE_BUILD:BOOL=[TRUE|FALSE]")
endif()

if(DMRG_ENABLE_ASAN)
    message(FATAL_ERROR "The option [DMRG_ENABLE_ASAN] has been deprecated. Replace by:\n"
            "COMPILER_ENABLE_ASAN:BOOL=[TRUE|FALSE]")
endif()

if(DMRG_ENABLE_USAN)
    message(FATAL_ERROR "The option [DMRG_ENABLE_USAN] has been deprecated. Replace by:\n"
            "COMPILER_ENABLE_USAN:BOOL=[TRUE|FALSE]")
endif()

if(DMRG_ENABLE_LTO)
    message(FATAL_ERROR "The option [DMRG_ENABLE_LTO] has been deprecated. Replace by:\n"
            "COMPILER_ENABLE_LTOCMAKE_INTERPROCEDURAL_OPTIMIZATION:BOOL=[TRUE|FALSE]")
endif()

if(DMRG_ENABLE_PCH)
    message(FATAL_ERROR "The option [DMRG_ENABLE_PCH] has been deprecated. Replace by:\n"
            "COMPILER_ENABLE_PCH:BOOL=[TRUE|FALSE]")
endif()

if(DMRG_ENABLE_CCACHE)
    message(FATAL_ERROR "The option [DMRG_ENABLE_CCACHE] has been deprecated. Replace by:\n"
            "COMPILER_ENABLE_CCACHE:BOOL=[TRUE|FALSE]")
endif()

if(COMPILER_MARCH)
    message(FATAL_ERROR "The option [COMPILER_MARCH] has been deprecated.")
endif()
if(COMPILER_MTUNE)
    message(FATAL_ERROR "The option [COMPILER_MTUNE] has been deprecated.")
endif()
if(COMPILER_ENABLE_LTO)
    message(FATAL_ERROR "The option [COMPILER_ENABLE_LTO] has been deprecated.. Replace by:\n"
            "CMAKE_INTERPROCEDURAL_OPTIMIZATION:BOOL=[TRUE|FALSE]")
endif()

if(DMRG_MARCH)
    message(FATAL_ERROR "The option [DMRG_MARCH] has been deprecated.")
endif()
if(DMRG_MTUNE)
    message(FATAL_ERROR "The option [DMRG_MTUNE] has been deprecated.")
endif()