
function(target_include_mkl_directories lpk tgt)
    message(CHECK_START "Looking for Intel MKL include directory in target ${tgt}")
    get_target_property(MKL_INCLUDE_DIR ${tgt} INTERFACE_INCLUDE_DIRECTORIES)
    if(NOT MKL_INCLUDE_DIR)
        message(CHECK_FAIL "not found")
        message(CHECK_START "Looking for Intel MKL include directory in system")
        set(MKL_ROOT_SEARCH_PATHS
            $ENV{MKLROOT} ${MKLROOT}
            $ENV{EBROOTIMKL} ${EBROOTIMKL}
            $ENV{BLASROOT} ${BLASROOT}
            /opt
            /opt/intel
            /opt/intel/oneapi
            )

        set(MKL_PATH_SUFFIXES
            intel/oneapi/mkl/latest
            oneapi/mkl/latest
            intel/oneapi/mkl
            oneapi/mkl
            intel/mkl
            mkl
            )

        find_path(MKL_ROOT_DIR
                  include/mkl.h
                  PATHS ${MKL_ROOT_SEARCH_PATHS}
                  PATH_SUFFIXES ${MKL_PATH_SUFFIXES}
                  )
        if(MKL_ROOT_DIR)
            find_path(MKL_INCLUDE_DIR
                      mkl_lapacke.h
                      HINTS ${MKL_ROOT_DIR}/include
                      )
        endif()
    endif()
    if(MKL_INCLUDE_DIR)
            message(CHECK_PASS "found: ${MKL_INCLUDE_DIR}")
            target_include_directories(${lpk} INTERFACE ${MKL_INCLUDE_DIR})
            target_compile_definitions(${lpk} INTERFACE MKL_AVAILABLE)
            set(LAPACKE_INCLUDE_DIR "${MKL_INCLUDE_DIR}" PARENT_SCOPE)
    else()
            message(CHECK_FAIL "not found")
            message(FATAL_ERROR "Could not find Intel MKL directory:\n"
                    "MKLROOT                : ${MKLROOT}\n"
                    "MKL_ROOT               : ${MKL_ROOT}\n"
                    "EBROOTIMKL             : ${EBROOTIMKL}\n"
                    "BLASROOT               : ${BLASROOT}\n"
                    "BLAS_ROOT              : ${BLAS_ROOT}\n"
                    "ENV MKLROOT            : $ENV{MKLROOT}\n"
                    "ENV MKL_ROOT           : $ENV{MKL_ROOT}\n"
                    "ENV EBROOTIMKL         : $ENV{EBROOTIMKL}\n"
                    "ENV BLASROOT           : $ENV{BLASROOT}\n"
                    "ENV BLAS_ROOT          : $ENV{BLAS_ROOT}\n"
                    "MKL_ROOT_SEARCH_PATHS  : ${MKL_ROOT_SEARCH_PATHS}\n"
                    "MKL_ROOT_DIR           : ${MKL_ROOT_DIR}\n"
                    "MKL_INCLUDE_DIR        : ${MKL_INCLUDE_DIR}\n")

    endif()
endfunction()

function(target_include_openblas_directories lpk tgt)
    message(CHECK_START "Looking for OpenBLAS include directory in target ${tgt}")
    get_target_property(OpenBLAS_INCLUDE_DIR ${tgt} INTERFACE_INCLUDE_DIRECTORIES)
    if(NOT OpenBLAS_INCLUDE_DIR)
        message(CHECK_FAIL "not found")
        message(CHECK_START "Looking for OpenBLAS include directory in system")
        set(OpenBLAS_ROOT_SEARCH_PATHS
          $ENV{OpenBLASROOT} ${OpenBLASROOT}
          $ENV{OpenBLAS_ROOT} ${OpenBLAS_ROOT}
          $ENV{OpenBLAS_HOME} ${OpenBLAS_HOME}
          $ENV{BLASROOT} ${BLASROOT}
          $ENV{EBROOTOPENBLAS} ${EBROOTOPENBLAS}
        )
        find_path(OpenBLAS_INCLUDE_DIR NAMES openblas/lapacke.h
            PATHS ${OpenBLAS_ROOT_SEARCH_PATHS}
            PATH_SUFFIXES include
        )
        find_path(LAPACKE_INCLUDE_DIR NAMES lapacke.h
            HINTS ${OpenBLAS_ROOT_SEARCH_PATHS}
            PATH_SUFFIXES include openblas include/openblas
        )
    endif()
    if(OpenBLAS_INCLUDE_DIR)
        message(CHECK_PASS "found: ${OpenBLAS_INCLUDE_DIR}")
        set(LAPACKE_INCLUDE_DIR "${OpenBLAS_INCLUDE_DIR}" PARENT_SCOPE)
        target_include_directories(${lpk} INTERFACE ${OpenBLAS_INCLUDE_DIR})
        target_compile_definitions(${lpk} INTERFACE OPENBLAS_AVAILABLE)
    elseif(LAPACKE_INCLUDE_DIR)
        message(CHECK_PASS "found: ${LAPACKE_INCLUDE_DIR}")
        set(LAPACKE_INCLUDE_DIR "${LAPACKE_INCLUDE_DIR}" PARENT_SCOPE)
        target_include_directories(${lpk} INTERFACE ${LAPACKE_INCLUDE_DIR})
        target_compile_definitions(${lpk} INTERFACE LAPACKE_AVAILABLE)
    else()
        message(CHECK_FAIL "not found")
        message(FATAL_ERROR "Could not find OpenBLAS include directory:\n"
                    "OpenBLASROOT               : ${OpenBLASROOT}\n"
                    "OpenBLAS_ROOT              : ${OpenBLAS_ROOT}\n"
                    "OpenBLAS_HOME              : ${OpenBLAS_HOME}\n"
                    "BLASROOT                   : ${BLASROOT}\n"
                    "EBROOTOPENBLAS             : ${EBROOTOPENBLAS}\n"
                    "ENV OpenBLASROOT           : ${OpenBLASROOT}\n"
                    "ENV OpenBLAS_ROOT          : ${OpenBLAS_ROOT}\n"
                    "ENV OpenBLAS_HOME          : ${OpenBLAS_HOME}\n"
                    "ENV BLASROOT               : ${BLASROOT}\n"
                    "ENV EBROOTOPENBLAS         : ${EBROOTOPENBLAS}\n"
                    "OpenBLAS_ROOT_SEARCH_PATHS : ${OpenBLAS_ROOT_SEARCH_PATHS}\n"
                    "OpenBLAS_INCLUDE_DIR       : ${OpenBLAS_INCLUDE_DIR}\n"
                    "LAPACKE_INCLUDE_DIR        : ${LAPACKE_INCLUDE_DIR}\n")

    endif()
endfunction()

function(target_include_lapacke_directories lpk src)
    message(DEBUG "FindLapacke: Inspecting BLA_VENDOR: ${BLA_VENDOR}")
    if(BLA_VENDOR MATCHES Intel OR $ENV{BLA_VENDOR} MATCHES Intel)
        target_include_mkl_directories(${lpk} ${src})
    elseif(BLA_VENDOR MATCHES OpenBLAS OR $ENV{BLA_VENDOR} MATCHES OpenBLAS)
        target_include_openblas_directories(${lpk} ${src})
    endif()
    if(LAPACKE_INCLUDE_DIR)
        message(DEBUG "Lapacke include directory: ${LAPACKE_INCLUDE_DIR}")
    endif()
endfunction()


function(find_Lapacke)
    if(BLA_VENDOR MATCHES FlexiBLAS OR ENV{BLA_VENDOR} MATCHES FlexiBLAS)
        find_package(lapacke REQUIRED CONFIG)
    elseif(BLA_VENDOR MATCHES Intel OR ENV{BLA_VENDOR} MATCHES Intel)
        set(ENABLE_BLAS95 ON)
        set(ENABLE_LAPACK95 ON)
        find_package(MKL REQUIRED BYPASS_PROVIDER)
    else()
        find_package(BLAS REQUIRED)
        find_package(LAPACK REQUIRED)
    endif()

    include(cmake/CheckCompile.cmake)

    set(lapacke_tgt_candidates lapacke Lapacke::Lapacke MKL::MKL LAPACK::LAPACK BLAS::BLAS mkl::mkl OpenBLAS::OpenBLAS)
    foreach(tgt ${lapacke_tgt_candidates})
        if (TARGET ${tgt} AND NOT TARGET lapacke::lapacke)
            message(DEBUG "Looking for Lapacke in ${tgt}")
            if (EXISTS ${PROJECT_SOURCE_DIR}/cmake/compile/Lapacke.cpp)
                if(NOT Lapacke_compiles)
                    check_compile(Lapacke ${tgt} ${PROJECT_SOURCE_DIR}/cmake/compile/Lapacke.cpp)
                endif()
                if (Lapacke_compiles)
                    add_library(lapacke::lapacke INTERFACE IMPORTED)
                    target_link_libraries(lapacke::lapacke INTERFACE ${tgt})
                    target_include_lapacke_directories(lapacke::lapacke ${tgt})
                    message(DEBUG "Looking for Lapacke in ${tgt} - found")
                    break()
                endif()
            else()
                add_library(lapacke::lapacke INTERFACE IMPORTED)
                target_link_libraries(lapacke::lapacke INTERFACE ${tgt})
                message(DEBUG "Assuming Lapacke is in ${tgt}")
                break()
            endif()
        endif()
    endforeach()
    if(NOT TARGET lapacke::lapacke)
        # Lapacke not found in the usual targets, perhaps OpenBLAS was built without lapacke, for instance
        message(DEBUG "Looking for lapacke directly")
        find_path(LAPACKE_INCLUDE_DIR NAMES lapacke.h
                    PATHS ${OpenBLAS_ROOT_SEARCH_PATHS}
                    PATH_SUFFIXES include
                    )

        find_library(LAPACKE_LIBRARY NAMES lapacke)

        if(LAPACKE_INCLUDE_DIR)
            message(DEBUG "Found LAPACKE_INCLUDE_DIR: ${LAPACKE_INCLUDE_DIR}")
            add_library(lapacke::lib INTERFACE IMPORTED )
            target_include_directories(lapacke::lib SYSTEM INTERFACE ${LAPACKE_INCLUDE_DIR})
            if(LAPACKE_LIBRARY)
                message(DEBUG "Found LAPACKE_LIBRARY: ${LAPACKE_LIBRARY}")
                target_link_libraries(lapacke::lib INTERFACE ${LAPACKE_LIBRARY})
                set_target_properties(lapacke::lib PROPERTIES IMPORTED_LOCATION ${LAPACKE_LIBRARY})
            endif()
            target_link_libraries(lapacke::lib INTERFACE BLAS::BLAS LAPACK::LAPACK)
            if(NOT Lapacke_compiles)
                check_compile(Lapacke lapacke::lib ${PROJECT_SOURCE_DIR}/cmake/compile/Lapacke.cpp)
            endif()
            if (Lapacke_compiles)
                add_library(lapacke::lapacke INTERFACE IMPORTED)
                target_link_libraries(lapacke::lapacke INTERFACE lapacke::lib)
                target_compile_definitions(lapacke::lapacke INTERFACE LAPACKE_AVAILABLE)
                message(DEBUG "Looking for Lapacke in lapacke::lib - found")
            endif()
        endif()
    endif()

endfunction()

find_Lapacke()
if(TARGET lapacke::lapacke)
    set(LAPACKE_TARGET lapacke::lapacke)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Lapacke DEFAULT_MSG LAPACKE_TARGET)
