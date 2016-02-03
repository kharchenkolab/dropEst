cmake_minimum_required(VERSION 2.8)
project(cp)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/CMake")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wno-deprecated -m64 -DBOOST_LOG_DYN_LINK")

set(BOOST_ROOT "/home/victor/libs/boost_1_58_0/")
find_package(Boost 1.54.0 COMPONENTS unit_test_framework serialization thread system log iostreams REQUIRED)
find_package (Threads REQUIRED)
find_package (ZLIB REQUIRED)
find_package (Bam REQUIRED)
find_package (R)

set(INCLUDE_DIRS ${Boost_INCLUDE_DIR} ${PROJECT_SOURCE_DIR} ${ZLIB_INCLUDE_DIRS} ${BAM_INCLUDE_DIRS})
set(LIBRARIES ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
set(EST_LIBRARIES ${BAM_LIBRARIES} ${ZLIB_LIBRARIES} ${LIBRARIES})

if (R_FOUND)
    add_definitions(-DR_LIBS)
    set(INCLUDE_DIRS ${INCLUDE_DIRS} ${R_INCLUDE_DIRS})
    set(LIBRARIES ${LIBRARIES} ${R_LIBRARIES})
endif()

include_directories(${INCLUDE_DIRS})

FILE(GLOB TagsSearchSources TagsSearch/*.cpp TagsSearch/Counters/*.cpp)
add_library(TagsSearch ${TagsSearchSources})

add_executable(indroptag indroptag.cpp)
target_link_libraries(indroptag TagsSearch ${LIBRARIES})


FILE(GLOB EstimationSources Estimation/*.cpp)
add_library(Estimation ${EstimationSources})

add_executable(indropest indropest.cpp)
target_link_libraries(indropest Estimation ${EST_LIBRARIES})



set (TEST TestTagsSearch)
add_executable (${TEST} ${PROJECT_SOURCE_DIR}/Tests/TestTagsSearch.cpp)
target_link_libraries (${TEST} TagsSearch ${LIBRARIES})
enable_testing ()
add_test (${TEST} ${TEST})