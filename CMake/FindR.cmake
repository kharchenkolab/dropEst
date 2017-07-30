# FindR
# --------
#
# Find R, RCpp, RCpArmadillo and RInside
#
# You can set R_ROOT and R_PACKAGES to specify R paths
#
# Result Variables
# ^^^^^^^^^^^^^^^^
#
# This module defines the following variables:
#
#   R_INCLUDE_DIRS: include directories for R
#   R_LIBRARIES: libraries to link against R
#   R_PACKAGES_DIRS: directories with installed R packages
#   R_FOUND: true if R has been found and can be used

# detection for OSX
if(APPLE)
    find_library(R_LIBRARIES R)
    if(R_LIBRARIES)
        set(R_ROOT "${R_LIBRARIES}/Resources" CACHE PATH "R root directory")
        set(R_INCLUDE_DIRS "${R_ROOT}/include" CACHE PATH "R include directory")
        set(R_EXECUTABLE "${R_ROOT}/R" CACHE PATH "R executable")
    endif()
# Find R executable and paths (UNIX)
else()
    # find executable
    if(NOT R_ROOT)
        find_program(R_EXECUTABLE R)
    else()
        find_program(R_EXECUTABLE R
                HINTS ${LIBRARY_ARCH_HINT_PATH} ${R_ROOT}/bin)
    endif()

    if(R_EXECUTABLE-NOTFOUND)
        message(STATUS "Unable to locate R executable")
    endif()

    # ask R for the home path
    if(NOT R_ROOT)
        execute_process(
                COMMAND ${R_EXECUTABLE} "--slave" "--vanilla" "-e" "cat(R.home())"
                OUTPUT_VARIABLE R_ROOT
        )
        if(R_ROOT)
            set(R_ROOT ${R_ROOT} CACHE PATH "R root directory")
        endif()
    endif()

    # ask R for the include dir
    if(NOT R_INCLUDE_DIRS)
        execute_process(
                COMMAND ${R_EXECUTABLE} "--slave" "--no-save" "-e" "cat(R.home('include'))"
                OUTPUT_VARIABLE R_INCLUDE_DIRS
        )
        if(R_INCLUDE_DIRS)
            set(R_INCLUDE_DIRS ${R_INCLUDE_DIRS} CACHE PATH "R include directory")
        endif()
    endif()

    # ask R for the lib dir
    if(NOT R_LIB_DIR)
        execute_process(
                COMMAND ${R_EXECUTABLE} "--slave" "--no-save" "-e" "cat(R.home('lib'))"
                OUTPUT_VARIABLE R_LIB_DIR
        )
    endif()

    # look for the core R library
    find_library(R_CORE_LIBRARY NAMES R
            HINTS ${R_LIB_DIR} ${R_ROOT}/bin)
    if(R_CORE_LIBRARY)
        set(R_LIBRARIES ${R_CORE_LIBRARY})
    else()
        message(STATUS "Could not find libR shared library.")
    endif()

    # cache R_LIBRARIES
    if(R_LIBRARIES)
        set(R_LIBRARIES ${R_LIBRARIES} CACHE PATH "R runtime libraries")
    endif()
endif()

# look for necessary R packages
if (NOT R_PACKAGES_DIRS)
    execute_process(
            COMMAND ${R_EXECUTABLE} "--slave" "--no-save" "-e" "cat(.libPaths(), sep=';')"
            OUTPUT_VARIABLE R_PACKAGES_DIRS
    )
endif()

find_path(RCPP_ARM_INCLUDE_DIR RcppArmadillo.h
        HINTS ${R_PACKAGES_DIRS}
        PATH_SUFFIXES /RcppArmadillo/include/)
find_path(RCPP_INCLUDE_DIR Rcpp.h
        HINTS ${R_PACKAGES_DIRS}
        PATH_SUFFIXES /Rcpp/include)
find_path(RINSIDE_INCLUDE_DIR RInside.h
        HINTS ${R_PACKAGES_DIRS}
        PATH_SUFFIXES /RInside/include)

find_library(RINSIDE_LIBRARY RInside PATH ${RINSIDE_INCLUDE_DIR}/../lib)

# Additional packages
execute_process(
        COMMAND ${R_EXECUTABLE} "--slave" "--no-save" "-e" "cat('Matrix' %in% rownames(installed.packages()))"
        OUTPUT_VARIABLE R_MATRIX_PKG
)

if(${R_MATRIX_PKG} STREQUAL FALSE)
    message(STATUS "Unable to locate R package 'Matrix'")
    set(R_MATRIX_PKG )
endif()

# define find requirements
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(R DEFAULT_MSG
        R_LIBRARIES
        R_INCLUDE_DIRS
        RCPP_INCLUDE_DIR
        RCPP_ARM_INCLUDE_DIR
        RINSIDE_INCLUDE_DIR
        RINSIDE_LIBRARY
        R_MATRIX_PKG
)

if (R_FOUND)
    set(R_INCLUDE_DIRS ${R_INCLUDE_DIRS} ${RCPP_INCLUDE_DIR} ${RCPP_ARM_INCLUDE_DIR} ${RINSIDE_INCLUDE_DIR})
    set(R_LIBRARIES ${R_LIBRARIES} ${RINSIDE_LIBRARY})
    message(STATUS "Found R: ${R_ROOT}")
endif()

# mark low-level variables from FIND_* calls as advanced
mark_as_advanced(
        R_CORE_LIBRARY
        RCPP_ARM_INCLUDE_DIR
        RCPP_INCLUDE_DIR
        RINSIDE_INCLUDE_DIR
        R_MATRIX_PKG
)