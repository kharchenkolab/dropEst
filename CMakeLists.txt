cmake_minimum_required(VERSION 2.8)

#set(CMAKE_CXX_COMPILER /opt/gcc/4.8.5/bin/g++)

project(cp)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/CMake")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wno-deprecated -m64 -DBOOST_LOG_DYN_LINK")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wno-deprecated -m64 -DBOOST_LOG_DYN_LINK -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb -gdwarf-2 -g3 -Wno-deprecated -m64 -DBOOST_LOG_DYN_LINK -std=c++11")
add_definitions(-DPROJ_DATA_PATH=\"${PROJECT_SOURCE_DIR}/data\")

#set(BOOST_ROOT "/home/vp76/libs/boost_1_58_0/")
find_package(Boost 1.54.0 COMPONENTS unit_test_framework serialization thread system log iostreams REQUIRED)
find_package (Threads REQUIRED)
find_package (ZLIB REQUIRED)

#set(BAMTOOLS_ROOT "/opt/bamtools-2.2.3/")
find_package (BamTools REQUIRED)

#set(R_ROOT "/home/vp76/local/lib64/R")
find_package (R)

if (R_FOUND)
    add_definitions(-DR_LIBS)
    set(INCLUDE_DIRS ${INCLUDE_DIRS} ${R_INCLUDE_DIRS})
    set(EST_LIBRARIES ${R_LIBRARIES})
endif()

set(INCLUDE_DIRS ${Boost_INCLUDE_DIR} ${PROJECT_SOURCE_DIR} ${ZLIB_INCLUDE_DIRS} ${BAMTOOLS_INCLUDE_DIRS} ${INCLUDE_DIRS})
include_directories(${INCLUDE_DIRS})

FILE(GLOB DropToolsSources Tools/*.cpp)
add_library(DropTools ${DropToolsSources})

FILE(GLOB TagsSearchSources TagsSearch/*.cpp TagsSearch/Counters/*.cpp)
add_library(TagsSearch ${TagsSearchSources})

FILE(GLOB EstimationSources Estimation/*.cpp Estimation/Results/*.cpp Estimation/Merge/*.cpp)
add_library(Estimation ${EstimationSources})
target_link_libraries(Estimation DropTools)

set(LIBRARIES DropTools ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
set(EST_LIBRARIES ${EST_LIBRARIES} ${BAMTOOLS_LIBRARIES} ${ZLIB_LIBRARIES} ${LIBRARIES})

add_executable(indroptag indroptag.cpp)
target_link_libraries(indroptag TagsSearch ${LIBRARIES})

add_executable(indropest indropest.cpp)
target_link_libraries(indropest Estimation ${EST_LIBRARIES})


set (TEST_TS TestTagsSearch)
add_executable (${TEST_TS} ${PROJECT_SOURCE_DIR}/Tests/TestTagsSearch.cpp)
target_link_libraries (${TEST_TS} TagsSearch ${LIBRARIES})
enable_testing ()
add_test (${TEST_TS} ${TEST_TS})


set (TEST_EST TestEstimation)
add_executable (${TEST_EST} ${PROJECT_SOURCE_DIR}/Tests/TestEstimation.cpp)
target_link_libraries (${TEST_EST} Estimation ${EST_LIBRARIES})
add_test (${TEST_EST} ${TEST_EST})


set (TEST_TOOLS TestTools)
add_executable (${TEST_TOOLS} ${PROJECT_SOURCE_DIR}/Tests/TestTools.cpp)
target_link_libraries (${TEST_TOOLS} DropTools ${LIBRARIES})
add_test (${TEST_TOOLS} ${TEST_TOOLS})
