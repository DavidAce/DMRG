


message(STATUS "SEARCHING FOR ARPACK IN SYSTEM...")
find_library(ARPACK_LIBRARIES
        NAMES arpack libarpack.a libarpack2.a libarpack.so libarpack2.so
        PATH_SUFFIXES lib lib32 lib64
        )

if(ARPACK_LIBRARIES)
    message(STATUS "ARPACK found in system:   ${ARPACK_LIBRARIES}")
    add_library(arpack UNKNOWN IMPORTED)
    set_target_properties(arpack PROPERTIES
            IMPORTED_LOCATION "${ARPACK_LIBRARIES}")
    target_link_libraries(${PROJECT_NAME} PUBLIC arpack)
    return()
else()
    message(STATUS "Arpack-ng will be installed into ${INSTALL_DIRECTORY}/arpack-ng on first build.")
    include(ExternalProject)
    ExternalProject_Add(library_ARPACK
            GIT_REPOSITORY      https://github.com/opencollab/arpack-ng.git
#            GIT_TAG             master # Latest version has problems with fortran linking. so stick with this version instead.
            GIT_TAG             3.5.0 # Latest version has problems with fortran linking. so stick with this version instead.
            PREFIX              "${INSTALL_DIRECTORY}/arpack-ng"
            UPDATE_COMMAND ""
#            BUILD_IN_SOURCE 1
#            CONFIGURE_COMMAND
#                ./bootstrap &&
#                ./configure --prefix=<INSTALL_DIR>
#            BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} && ${CMAKE_MAKE_PROGRAM} check
#            INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install

            CMAKE_ARGS
            -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
            -DCMAKE_INSTALL_MESSAGE=NEVER #Avoid unnecessary output to console
            -DCMAKE_C_FLAGS=-w
            -DEXAMPLES=OFF
            -DCMAKE_BUILD_TYPE=Release
            -DMPI=OFF
            -DBUILD_SHARED_LIBS=OFF
            -DBLAS_LIBRARIES:PATH=${BLAS_LIBRARIES}
            -DLAPACK_LIBRARIES:PATH=${LAPACK_LIBRARIES}
            DEPENDS blas lapack
            )

    ExternalProject_Get_Property(library_ARPACK INSTALL_DIR)
    set(ARPACK_INCLUDE_DIRS ${INSTALL_DIR}/include)
    add_library(arpack UNKNOWN IMPORTED)
    set_target_properties(arpack PROPERTIES
            IMPORTED_LOCATION ${INSTALL_DIR}/lib/libarpack${CMAKE_STATIC_LIBRARY_SUFFIX}
            INTERFACE_LINK_LIBRARIES "${GFORTRAN_LIB};${BLAS_LIBRARIES};${LAPACK_LIBRARIES}"
            INCLUDE_DIRECTORIES ${INSTALL_DIR}/include)
    add_dependencies(arpack library_ARPACK)

    target_link_libraries(${PROJECT_NAME} PUBLIC arpack)
    target_include_directories(${PROJECT_NAME} PUBLIC ${ARPACK_INCLUDE_DIRS})
    #For convenience, define these variables
    get_target_property(ARPACK_LIBRARIES arpack IMPORTED_LOCATION)
endif()