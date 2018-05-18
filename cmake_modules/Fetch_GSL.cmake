

enable_language(Fortran)
include(cmake_modules/FindGFortran.cmake)

find_package(GSL)
if(GSL_FOUND)
    message(STATUS "GSL FOUND IN SYSTEM: ${GSL_LIBRARIES}")
else()
    message(STATUS "GSL will be installed into ${INSTALL_DIRECTORY}/gsl on first build.")
    include(ExternalProject)
    ExternalProject_Add(library_GSL
            URL      http://ftp.acc.umu.se/mirror/gnu.org/gnu/gsl/gsl-2.4.tar.gz
            PREFIX              "${INSTALL_DIRECTORY}/gsl"
            CONFIGURE_COMMAND
                cd <SOURCE_DIR> &&
                pwd &&
                ./configure --enable-silent-rules CFLAGS= --enable-shared=no --prefix=<INSTALL_DIR>
            BUILD_COMMAND
                cd <SOURCE_DIR> &&
                pwd &&
                ${CMAKE_MAKE_PROGRAM} --quiet --silent
            INSTALL_COMMAND
                cd <SOURCE_DIR> &&
                pwd &&
                ${CMAKE_MAKE_PROGRAM} --quiet install
            )

    ExternalProject_Get_Property(library_GSL INSTALL_DIR)
    set(GSL_INCLUDE_DIRS ${INSTALL_DIR}/include)
    add_library(GSL::gsl UNKNOWN IMPORTED)
    set_target_properties(GSL::gsl PROPERTIES
            IMPORTED_LOCATION ${INSTALL_DIR}/lib/libgsl${CMAKE_STATIC_LIBRARY_SUFFIX}
            INCLUDE_DIRECTORIES ${GSL_INCLUDE_DIRS})
    add_dependencies(GSL::gsl library_GSL)

    add_library(GSL::gslcblas UNKNOWN IMPORTED)
    set_target_properties(GSL::gslcblas PROPERTIES
            IMPORTED_LOCATION ${INSTALL_DIR}/lib/libgslcblas${CMAKE_STATIC_LIBRARY_SUFFIX}
            INCLUDE_DIRECTORIES ${GSL_INCLUDE_DIRS})
    add_dependencies(GSL::gslcblas library_GSL)
endif()


target_link_libraries(${PROJECT_NAME} PRIVATE GSL::gsl)
target_link_libraries(${PROJECT_NAME} PRIVATE GSL::gslcblas)
target_include_directories(${PROJECT_NAME} PRIVATE ${GSL_INCLUDE_DIRS})
get_target_property(GSL_LIBRARIES        GSL::gsl       IMPORTED_LOCATION)
get_target_property(GSLCBLAS_LIBRARIES   GSL::gslcblas  IMPORTED_LOCATION)
