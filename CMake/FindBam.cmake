# FindBam
# --------
# You can set BAM_ROOT variable to specify path
# Result Variables
# ^^^^^^^^^^^^^^^^
#
# This module defines the following variables:
#
#   BAM_INCLUDE_DIRS - include directories for BAM
#   BAM_LIBRARIES - libraries to link against BAM
#   BAM_FOUND - true if BAM has been found and can be used

find_path(BAM_INCLUDE_DIR bam.h PATHS ${BAM_ROOT}/include PATH_SUFFIXES samtools)
find_library(BAM_LIBRARY bam PATHS ${BAM_ROOT}/lib)

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(BAM
        REQUIRED_VARS BAM_LIBRARY BAM_INCLUDE_DIR)

if (BAM_FOUND)
    set(BAM_INCLUDE_DIRS ${BAM_INCLUDE_DIR})
    set(BAM_LIBRARIES ${BAM_LIBRARY})
endif()


mark_as_advanced(BAM_INCLUDE_DIR BAM_LIBRARY)