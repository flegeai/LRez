################################################################################
#  TOOLS
################################################################################

MACRO(SUBDIRLIST result curdir)
    FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
    SET (dirlist "")
    FOREACH(child ${children})
      IF(IS_DIRECTORY ${curdir}/${child})
          LIST(APPEND dirlist ${child})
      ENDIF()
    ENDFOREACH()
    SET(${result} ${dirlist})
ENDMACRO()

# We get the list of subdirectories names
SUBDIRLIST(PROGRAMS ${PROJECT_SOURCE_DIR}/tools)

find_package(HTSlib REQUIRED)

# We create a target for each 
FOREACH (program ${PROGRAMS})
 MESSAGE("Program : " ${program})
  set (PROGRAM_SOURCE_DIR ${PROJECT_SOURCE_DIR}/tools/${program}/src)
  include_directories (${PROGRAM_SOURCE_DIR} ${HTSlib_INCLUDE_DIRS})
  file (GLOB_RECURSE  ProjectFiles  ${PROGRAM_SOURCE_DIR}/*.cpp)
  MESSAGE("Prfojectfiles : " ${ProjectFiles})
  add_executable(${program} ${ProjectFiles})
  target_link_libraries(${program} ${HTSlib_LIBRARIES})
ENDFOREACH (program)

################################################################################
#  INSTALLATION 
################################################################################
install (TARGETS ${PROGRAMS} DESTINATION bin)
