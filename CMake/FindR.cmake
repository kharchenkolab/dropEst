# FindR
# --------
#
# Find R, RCpp and RInside
#
# You can set R_ROOT to specify R paths
#
# Result Variables
# ^^^^^^^^^^^^^^^^
#
# This module defines the following variables:
#
#   R_INCLUDE_DIRS - include directories for R
#   R_LIBRARIES - libraries to link against R
#   R_FOUND - true if R has been found and can be used

find_path(R_INCLUDE_DIR R.h PATHS ${R_ROOT}/include /usr/share/R/include)
find_library(R_LIBRARY NAME R PATHS ${R_ROOT}/lib PATH_SUFFIXES R/lib lib)

get_filename_component(R_LIBRARY ${R_LIBRARY} REALPATH)
get_filename_component(R_LIB_DIR ${R_LIBRARY} PATH)
set(R_LIBS_PATHS ${R_ROOT}/library ${R_LOCAL_LIBS} ${R_ROOT}/site-library ${R_LIB_DIR}/../ ${R_LIB_DIR}/../../ /usr/local/lib /usr/lib)


find_path(RCPP_INCLUDE_DIR Rcpp.h PATHS ${R_LIBS_PATHS} PATH_SUFFIXES R/site-library/Rcpp/include/ Rcpp/include/)
find_path(RINSIDE_INCLUDE_DIR RInside.h PATHS ${R_LIBS_PATHS} PATH_SUFFIXES R/site-library/RInside/include RInside/include)

find_library(RINSIDE_LIBRARY RInside PATH ${RINSIDE_INCLUDE_DIR}/../lib)

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(R
        REQUIRED_VARS R_LIBRARY R_INCLUDE_DIR RCPP_INCLUDE_DIR  RINSIDE_INCLUDE_DIR RINSIDE_LIBRARY)

if (R_FOUND)
    set(R_INCLUDE_DIRS ${R_INCLUDE_DIR} ${RCPP_INCLUDE_DIR} ${RINSIDE_INCLUDE_DIR})
    set(R_LIBRARIES ${R_LIBRARY} ${RINSIDE_LIBRARY})
endif()


mark_as_advanced(R_INCLUDE_DIR R_LIBRARY RCPP_INCLUDE_DIR RINSIDE_INCLUDE_DIR RINSIDE_LIBRARY)
