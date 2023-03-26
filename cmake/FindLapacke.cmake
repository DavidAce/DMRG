
function(target_include_mkl_directories tgt)
    set(MKL_ROOT_SEARCH_PATHS
        ${MKL_ROOT_DIR}
        $ENV{MKL_DIR} ${MKL_DIR}
        $ENV{MKLDIR} ${MKLDIR}
        $ENV{MKLROOT} ${MKLROOT}
        $ENV{MKL_ROOT} ${MKL_ROOT}
        $ENV{mkl_root} ${mkl_root}
        $ENV{EBROOTIMKL} ${EBROOTIMKL}
        $ENV{HOME}
        /opt
        /opt/intel
        /opt/intel/oneapi
        /usr/lib/x86_64-linux-gnu
        /usr
        /Library/Frameworks/Intel_MKL.framework/Versions/Current/lib/universal
        "Program Files (x86)/Intel/ComposerXE-2011/mkl"
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
              HINTS ${MKL_ROOT_SEARCH_PATHS}
              PATH_SUFFIXES ${MKL_PATH_SUFFIXES}
              )
    if(MKL_ROOT_DIR)
        find_path(MKL_INCLUDE_DIR
                  mkl_lapacke.h
                  HINTS ${MKL_ROOT_DIR}/include
                  )

        if(MKL_INCLUDE_DIR)
            target_include_directories(${tgt} INTERFACE ${MKL_INCLUDE_DIR})
            target_compile_definitions(${tgt} INTERFACE MKL_AVAILABLE)
            set(LAPACKE_INCLUDE_DIR "${MKL_INCLUDE_DIR}" PARENT_SCOPE)
        else()
            message(WARNING "Found MKL_ROOT_DIR but not MKL_INCLUDE_DIR:"
                    "MKL_ROOT_DIR   : ${MKL_ROOT_DIR})"
                    "MKL_INCLUDE_DIR: ${MKL_INCLUDE_DIR})")
        endif()
    endif()
endfunction()

function(target_include_openblas_directories tgt)
    set(OpenBLAS_ROOT_SEARCH_PATHS
      $ENV{OpenBLAS_ROOT} ${OpenBLAS_ROOT}
      $ENV{OpenBLAS_HOME} ${OpenBLAS_HOME}
      $ENV{BLAS_ROOT} ${BLAS_ROOT}
      $ENV{BLASROOT} ${BLASROOT}
      $ENV{EBROOTOPENBLAS} ${EBROOTOPENBLAS}
    )
    find_path(OpenBLAS_INCLUDE_DIR NAMES openblas/lapacke.h
        HINTS ${OpenBLAS_ROOT_SEARCH_PATHS}
        PATH_SUFFIXES include
    )
    find_path(LAPACKE_INCLUDE_DIR NAMES lapacke.h
        HINTS ${OpenBLAS_ROOT_SEARCH_PATHS}
        PATH_SUFFIXES include openblas include/openblas
    )

    if(OpenBLAS_INCLUDE_DIR)
        set(LAPACKE_INCLUDE_DIR "${OpenBLAS_INCLUDE_DIR}" PARENT_SCOPE)
        target_include_directories(${tgt} INTERFACE ${OpenBLAS_INCLUDE_DIR})
        target_compile_definitions(${tgt} INTERFACE OPENBLAS_AVAILABLE)
    elseif(LAPACKE_INCLUDE_DIR)
        set(LAPACKE_INCLUDE_DIR "${LAPACKE_INCLUDE_DIR}" PARENT_SCOPE)
        target_include_directories(${tgt} INTERFACE ${LAPACKE_INCLUDE_DIR})
        target_compile_definitions(${tgt} INTERFACE LAPACKE_AVAILABLE)
    endif()
endfunction()

function(target_include_lapacke_directories tgt)
    message(DEBUG "FindLapacke: Inspecting BLA_VENDOR: ${BLA_VENDOR}")
    if(BLA_VENDOR MATCHES Intel OR $ENV{BLA_VENDOR} MATCHES Intel)
        target_include_mkl_directories(${tgt})
    elseif(BLA_VENDOR MATCHES OpenBLAS OR $ENV{BLA_VENDOR} MATCHES OpenBLAS)
        target_include_openblas_directories(${tgt})
    endif()
    if(LAPACKE_INCLUDE_DIR)
        message(DEBUG "Lapacke include directory: ${LAPACKE_INCLUDE_DIR}")
    endif()
endfunction()


function(find_Lapacke)
    if(NOT DEFINED BLA_VENDOR AND NOT DEFINED $ENV{BLA_VENDOR})
        message(FATAL_ERROR "BLA_VENDOR is unset" )
    endif()
    if(NOT DEFINED BLA_VENDOR AND DEFINED $ENV{BLA_VENDOR})
        set(BLA_VENDOR $ENV{BLA_VENDOR})
    endif()

    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
    include(cmake/CheckCompile.cmake)

    set(lapacke_tgt_candidates LAPACK::LAPACK mkl::mkl OpenBLAS::OpenBLAS)
    foreach(tgt ${lapacke_tgt_candidates})
        if (TARGET ${tgt} AND NOT TARGET lapacke::lapacke)
            message(DEBUG "Looking for Lapacke in ${tgt}")
            if (EXISTS ${PROJECT_SOURCE_DIR}/cmake/compile/Lapacke.cpp)
                unset(check_compile_Lapacke CACHE)
                check_compile(Lapacke ${tgt} ${PROJECT_SOURCE_DIR}/cmake/compile/Lapacke.cpp)
                if (check_compile_Lapacke)
                    add_library(lapacke::lapacke INTERFACE IMPORTED)
                    target_link_libraries(lapacke::lapacke INTERFACE ${tgt})
                    target_include_lapacke_directories(lapacke::lapacke)
                    message(DEBUG "Looking for Lapacke in ${tgt} - found")
                    break()
                else ()
                    unset(check_compile_Lapacke CACHE)
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

        find_library(LAPACK_LIBRARY NAMES lapack)
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
            if(LAPACKE_LIBRARY)
                message(DEBUG "Found LAPACK_LIBRARY: ${LAPACK_LIBRARY}")
                target_link_libraries(lapacke::lib INTERFACE ${LAPACK_LIBRARY})
            endif()
            target_link_libraries(lapacke::lib INTERFACE BLAS::BLAS LAPACK::LAPACK)
            unset(check_compile_Lapacke CACHE)
            check_compile(Lapacke lapacke::lib ${PROJECT_SOURCE_DIR}/cmake/compile/Lapacke.cpp)
            if (check_compile_Lapacke)
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
