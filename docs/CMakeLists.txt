set(Doxygen_ROOT
        ${CONDA_PREFIX}
        $ENV{HOME}/anaconda3
        $ENV{HOME}/miniconda3
        $ENV{HOME}/anaconda3/envs/dmrg
        $ENV{HOME}/miniconda3/envs/dmrg)

set(Sphinx_ROOT
        ${CONDA_PREFIX}
        $ENV{HOME}/anaconda3
        $ENV{HOME}/miniconda3
        $ENV{HOME}/anaconda3/envs/dmrg
        $ENV{HOME}/miniconda3/envs/dmrg
        )

set(DOXYGEN_GENERATE_XML YES)

find_package(Doxygen REQUIRED)
doxygen_add_docs(dmrg-doxygen
        ${PROJECT_SOURCE_DIR}/source
        COMMENT "Generating docs")


find_package(Sphinx REQUIRED)

set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR})
set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/docs/sphinx)

add_custom_target(dmrg-sphinx ALL
        COMMAND
        ${SPHINX_EXECUTABLE} -b html
        -Dbreathe_projects.DMRG++=${CMAKE_CURRENT_BINARY_DIR}/xml # Tell Breathe where to find the Doxygen output
        ${SPHINX_SOURCE} ${SPHINX_BUILD}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating documentation with Sphinx")



# Find all the public headers
#get_target_property(DMRG_INCLUDE_DIR DMRG++ INTERFACE_INCLUDE_DIRECTORIES)
#file(GLOB_RECURSE DMRG_HEADERS ${DMRG_INCLUDE_DIR}/*.h)
#
#
#set(DOXYGEN_INPUT_DIR ${PROJECT_SOURCE_DIR}/source)
#set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/docs/doxygen)
#set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIR}/html/index.html)
#set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
#set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
#
##Replace variables inside @@ with the current values
#configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)
#
#file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR}) #Doxygen won't create this for us
#add_custom_command(OUTPUT ${DOXYGEN_INDEX_FILE}
#                   DEPENDS ${DMRG_HEADERS}
#                   COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
#                   MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN}
#                   COMMENT "Generating docs")
#
#
#add_custom_target(Doxygen ALL DEPENDS ${DOXYGEN_INDEX_FILE})



