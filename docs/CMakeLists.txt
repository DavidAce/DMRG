
list(APPEND CONDA_ROOT
        ${CONDA_PREFIX}
        $ENV{HOME}/anaconda3
        $ENV{HOME}/miniconda3
        $ENV{HOME}/anaconda3/envs/dmrg
        $ENV{HOME}/miniconda3/envs/dmrg)

list(APPEND Doxygen_ROOT ${CONDA_ROOT})
list(APPEND Sphinx_ROOT ${CONDA_ROOT})

set(DOXYGEN_GENERATE_XML YES)
set(DOXYGEN_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

find_package(Doxygen REQUIRED)
doxygen_add_docs(dmrg-doxygen
        ALL
        ${PROJECT_SOURCE_DIR}/source
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMENT "Generating Doxygen XML")


find_package(Sphinx REQUIRED)

set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/source)
set(SPHINX_BUILD ${CMAKE_CURRENT_SOURCE_DIR}/sphinx)

add_custom_target(dmrg-sphinx ALL
        COMMAND
        ${SPHINX_EXECUTABLE} -b html
        -Dbreathe_projects.DMRG++=${CMAKE_CURRENT_SOURCE_DIR}/xml # Tell Breathe where to find the Doxygen XML output
        ${SPHINX_SOURCE} ${SPHINX_BUILD}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/source
        COMMENT "Generating documentation with Sphinx")
