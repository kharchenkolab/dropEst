if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.8) # TODO: test it
        message(FATAL_ERROR "DropEst requires gcc version 4.8 or later")
    endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.2)
        message(FATAL_ERROR "DropEst requires clang version 3.2 or later")
    endif()
else()
    message(WARNING "Unsupported compiler is detected. DropEst compilation was not tested on it and may fail")
endif()

set(Boost_USE_STATIC_LIBS        OFF)
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME     OFF)
add_definitions(-DBOOST_ALL_DYN_LINK)

# set(BOOST_ROOT "/home/vp76/local/")
find_package(Boost 1.54.0 COMPONENTS filesystem iostreams log system thread unit_test_framework REQUIRED) # Thread is required for log on some systems

if(Boost_FOUND)
    message(STATUS "** Boost Include: ${Boost_INCLUDE_DIR}")
    message(STATUS "** Boost Libraries: ${Boost_LIBRARY_DIRS}")
    message(STATUS "** Boost Link-Libs: ${Boost_LIBRARIES}")
endif(Boost_FOUND)

find_package (ZLIB REQUIRED)

#set(BAMTOOLS_ROOT "/opt/bamtools-2.2.3/")
find_package (BamTools REQUIRED)
find_package (Threads REQUIRED)

# set(R_ROOT "/opt/R-3.3.1/lib64/R/")
find_package (R REQUIRED)