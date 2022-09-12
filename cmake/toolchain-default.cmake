foreach(lang C CXX CUDA)
    set(CMAKE_${lang}_STANDARD 17)
    set(CMAKE_${lang}_STANDARD_REQUIRED ON)
    set(CMAKE_${lang}_EXTENSIONS OFF)

endforeach()
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE STRING "") ### Write compile commands to file
set(CMAKE_CXX_FLAGS_INIT "-g -fno-strict-aliasing -fdiagnostics-color=always" CACHE STRING "")
set(CMAKE_CXX_FLAGS_RELEASE_INIT "-ffp-contract=fast" CACHE STRING "")
set(CMAKE_CXX_FLAGS_DEBUG_INIT "-fno-omit-frame-pointer -fstack-protector-strong -D_FORTIFY_SOURCE=2" CACHE STRING "")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "-fno-omit-frame-pointer -fstack-protector-strong -D_FORTIFY_SOURCE=2" CACHE STRING "")

set(CMAKE_CUDA_FLAGS_INIT "-Xcompiler -fno-strict-aliasing -Xcompiler -fdiagnostics-color=always --expt-relaxed-constexpr -Xcudafe --diag_suppress=esa_on_defaulted_function_ignored" CACHE STRING "")
set(CMAKE_CUDA_ARCHITECTURES "75;80;86" CACHE STRING "")
set(CMAKE_CUDA_SEPARABLE_COMPILATION "OFF" CACHE STRING "")
set(CMAKE_CUDA_HOST_COMPILER "g++-10" CACHE STRING "") # CUDA 11 does not support versions > 10

set(CMAKE_FIND_LIBRARY_USE_LIB64_PATHS ON)
set(CMAKE_FIND_LIBRARY_USE_LIB32_PATHS ON)

###############################
# Settings for shared builds
# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already (but later on when installing)
# Note: Since TB++ is often run from the build folder we want to keep the build-folder RPATH in the executable.
#       Therefore itt makes sense to keep this setting "FALSE" here but "TRUE" for dependencies that are
#       installed with in "fetch" mode with externalproject_add
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


