project(LRez)

cmake_minimum_required(VERSION 3.1.0)
set (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
MESSAGE("CMAKE MODULE PATH : "  ${CMAKE_MODULE_PATH} )

################################################################################
# TOOLS
################################################################################

ADD_SUBDIRECTORY(tools)

################################################################################
# Packaging
################################################################################
# We set the version number
SET (CPACK_PACKAGE_DESCRIPTION_SUMMARY  "LRez")
SET (CPACK_PACKAGE_VENDOR               "INRAE")
SET (CPACK_PACKAGE_VERSION_MAJOR        1)
SET (CPACK_PACKAGE_VERSION_MINOR        0)
SET (CPACK_PACKAGE_VERSION_PATCH        0)
SET (CPACK_PACKAGE_VERSION              ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${PACK_PACKAGE_VERSION_PATCH})

# We set the kind of archive
SET (CPACK_GENERATOR                    "TGZ")
SET (CPACK_SOURCE_GENERATOR             "TGZ")

# Packaging the source ; we ignore unwanted files 
SET (CPACK_SOURCE_IGNORE_FILES          
    "^${PROJECT_SOURCE_DIR}/build/"  
)

# Packaging the binary ; we want to include some additional files
INSTALL (FILES   ${CMAKE_CURRENT_SOURCE_DIR}/LICENCE                     DESTINATION .)

include (CPack)

