# FindBam
# --------
# You can set BAM_ROOT variable to specify path
#
# This module defines the following variables:
#
#   BAMTOOLS_INCLUDE_DIRS - include directories for BAM
#   BAMTOOLS_LIBRARIES - libraries to link against BAM
#   BAMTOOLS_FOUND - true if BAM has been found and can be used

find_path(BAMTOOLS_INCLUDE_DIR api/BamReader.h PATHS ${BAMTOOLS_ROOT}/include PATH_SUFFIXES bamtools ../include ../include/bamtools)
find_library(BAMTOOLS_LIBRARY bamtools PATHS ${BAMTOOLS_ROOT}/lib64 ${BAMTOOLS_ROOT}/lib PATH_SUFFIXES bamtools ../lib)

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(BAMTOOLS
        REQUIRED_VARS BAMTOOLS_LIBRARY BAMTOOLS_INCLUDE_DIR)

if (BAMTOOLS_FOUND)
    set(BAMTOOLS_INCLUDE_DIRS ${BAMTOOLS_INCLUDE_DIR})
    set(BAMTOOLS_LIBRARIES ${BAMTOOLS_LIBRARY})
endif()

mark_as_advanced(BAMTOOLS_INCLUDE_DIR BAMTOOLS_LIBRARY)
