
if(DMRG_DOWNLOAD_METHOD MATCHES "conan")
    ##################################################################
    ### Install conan-modules/conanfile.txt dependencies          ###
    ### This uses conan to get spdlog,eigen3,h5pp,ceres-solver    ###
    ###    ceres-solver/2.0.0@davidace/development                ###
    ###    h5pp/1.5.1@davidace/stable                             ###
    ###    eigen/3.3.7@davidace/patched                           ###
    ##################################################################
    if(DMRG_DOWNLOAD_METHOD MATCHES "find|fetch|native")
        if(NOT TARGET Eigen3::Eigen)
            list(APPEND CONAN_REQUIRES REQUIRES eigen/3.3.7@davidace/patched)
        endif()
        if(NOT TARGET h5pp::h5pp)
            list(APPEND CONAN_REQUIRES REQUIRES h5pp/1.7.0@davidace/stable)
        endif()
        if(NOT TARGET ceres::ceres)
            list(APPEND CONAN_REQUIRES REQUIRES ceres-solver/2.0.0@davidace/development)
        endif()
    else()
        list(APPEND CONAN_REQUIRES CONANFILE conanfile.txt)
    endif()


    find_program (
            CONAN_COMMAND
            conan
            HINTS ${CONAN_PREFIX} ${CONDA_PREFIX} $ENV{CONAN_PREFIX} $ENV{CONDA_PREFIX}
            PATHS $ENV{HOME}/anaconda3 $ENV{HOME}/miniconda3 $ENV{HOME}/.conda $ENV{HOME}/.local
            PATH_SUFFIXES bin envs/dmrg/bin
    )
    message(STATUS "Found conan: ${CONAN_COMMAND}")

    # Download automatically, you can also just copy the conan.cmake file
    if(NOT EXISTS "${CMAKE_BINARY_DIR}/conan.cmake")
        message(STATUS "Downloading conan.cmake from https://github.com/conan-io/cmake-conan")
        file(DOWNLOAD "https://github.com/conan-io/cmake-conan/raw/v0.15/conan.cmake"
                "${CMAKE_BINARY_DIR}/conan.cmake")
    endif()

    include(${CMAKE_BINARY_DIR}/conan.cmake)
    conan_add_remote(NAME bincrafters        URL https://api.bintray.com/conan/bincrafters/public-conan)
    conan_add_remote(NAME conan-community    URL https://api.bintray.com/conan/conan-community/conan)
    conan_add_remote(NAME conan-dmrg INDEX 1 URL https://api.bintray.com/conan/davidace/conan-dmrg)

    if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang")
        # Let it autodetect libcxx
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
        # There is no libcxx
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        list(APPEND conan_libcxx compiler.libcxx=libstdc++11)
    endif()
    conan_cmake_run(
            ${CONAN_REQUIRES}
            CONAN_COMMAND ${CONAN_COMMAND}
            SETTINGS compiler.cppstd=17
            SETTINGS "${conan_libcxx}"
            BUILD_TYPE ${CMAKE_BUILD_TYPE}
            BASIC_SETUP CMAKE_TARGETS
            BUILD missing
    )
endif()
