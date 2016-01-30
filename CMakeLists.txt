cmake_minimum_required(VERSION 3.2)
project(cp)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/CMake")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wno-deprecated -m64 -DBOOST_LOG_DYN_LINK")

set(BOOST_ROOT "/home/victor/libs/boost_1_58_0/")
find_package(Boost COMPONENTS unit_test_framework serialization thread system log iostreams REQUIRED)
find_package (Threads REQUIRED)
find_package (ZLIB REQUIRED)

set(RLIB_DIR "/usr/local/lib/R/site-library")

include_directories(${Boost_INCLUDE_DIR} ${PROJECT_SOURCE_DIR} ${ZLIB_INCLUDE_DIRS}
        /usr/include/samtools ${RLIB_DIR}/Rcpp/include/ ${RLIB_DIR}/RInside/include /usr/share/R/include/)
set(LIBRARIES ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

FILE(GLOB TagsSearchSources TagsSearch/*.cpp TagsSearch/Counters/*.cpp)
add_library(TagsSearch ${TagsSearchSources})

add_executable(indroptag indroptag.cpp)
target_link_libraries(indroptag TagsSearch ${LIBRARIES})



FILE(GLOB EstimationSources Estimation/*.cpp)
add_library(Estimation ${EstimationSources})

set(EST_LIBRARIES /usr/lib/libbam.a /usr/lib/R/lib/libR.so ${RLIB_DIR}/RInside/lib/libRInside.so ${ZLIB_LIBRARIES} ${LIBRARIES})

add_executable(indropest indropest.cpp)
target_link_libraries(indropest Estimation ${EST_LIBRARIES})



set (TEST TestTagsSearch)
add_executable (${TEST} ${PROJECT_SOURCE_DIR}/Tests/TestTagsSearch.cpp)
target_link_libraries (${TEST} TagsSearch ${LIBRARIES})
enable_testing ()
add_test (${TEST} ${TEST})