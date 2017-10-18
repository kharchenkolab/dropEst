cmake_minimum_required(VERSION 2.8)

project(dropEst)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/CMake")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -std=c++11 -Wno-attributes -Wno-deprecated -Wno-deprecated-declarations")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -ggdb -gdwarf-2 -g3")

if (NOT "${CMAKE_BUILD_TYPE}" STREQUAL Debug)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")
endif (NOT "${CMAKE_BUILD_TYPE}" STREQUAL Debug)

add_definitions(-DPROJ_DATA_PATH=\"${PROJECT_SOURCE_DIR}/data\")
add_definitions(-DPROJ_BIN_PATH=\"${CMAKE_CURRENT_BINARY_DIR}\")
add_definitions(-DVERSION=\"0.7.1\")

include(deps)

set(INCLUDE_DIRS ${Boost_INCLUDE_DIR} ${PROJECT_SOURCE_DIR} ${ZLIB_INCLUDE_DIRS} ${BAMTOOLS_INCLUDE_DIRS} ${R_INCLUDE_DIRS})
include_directories(${INCLUDE_DIRS})

set(BASE_LIBRARIES ${R_LIBRARIES} ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

FILE(GLOB DropToolsSources Tools/*.cpp)
add_library(DropTools ${DropToolsSources})
target_link_libraries(DropTools ${BASE_LIBRARIES})

FILE(GLOB TagsSearchSources TagsSearch/*.cpp TagsSearch/Counters/*.cpp)
add_library(TagsSearch ${TagsSearchSources})

FILE(GLOB EstimationSources Estimation/*.cpp  Estimation/Merge/*.cpp Estimation/Merge/UMIs/*.cpp
        Estimation/Merge/BarcodesParsing/*.cpp Estimation/BamProcessing/*.cpp)
add_library(Estimation ${EstimationSources})
target_link_libraries(Estimation DropTools)


set(LIBRARIES DropTools ${BASE_LIBRARIES})
set(EST_LIBRARIES ${EST_LIBRARIES} ${BAMTOOLS_LIBRARIES} ${ZLIB_LIBRARIES} ${LIBRARIES})

add_executable(droptag droptag.cpp)
target_link_libraries(droptag TagsSearch ${LIBRARIES})

add_executable(dropest dropest.cpp)
target_link_libraries(dropest Estimation ${EST_LIBRARIES})

#TODO: check it
add_custom_command(TARGET dropest POST_BUILD
        COMMAND cp
        ${CMAKE_SOURCE_DIR}/scripts/Report/*.R ${CMAKE_CURRENT_BINARY_DIR}/)

set (TEST_TS TestTagsSearch)
add_executable (${TEST_TS} ${PROJECT_SOURCE_DIR}/Tests/TestTagsSearch.cpp)
target_link_libraries (${TEST_TS} TagsSearch ${LIBRARIES})
enable_testing ()
add_test (${TEST_TS} ${TEST_TS})

set (TEST_EST TestEstimation)
add_executable (${TEST_EST} ${PROJECT_SOURCE_DIR}/Tests/TestEstimation.cpp)
target_link_libraries (${TEST_EST} Estimation ${EST_LIBRARIES})
add_test (${TEST_EST} ${TEST_EST})

set (TEST_EST_MERGE_PROBS TestEstimationMergeProbs)
add_executable (${TEST_EST_MERGE_PROBS} ${PROJECT_SOURCE_DIR}/Tests/TestEstimationMergeProbs.cpp)
target_link_libraries (${TEST_EST_MERGE_PROBS} Estimation ${EST_LIBRARIES})
add_test (${TEST_EST_MERGE_PROBS} ${TEST_EST_MERGE_PROBS})

set (TEST_TOOLS TestTools)
add_executable (${TEST_TOOLS} ${PROJECT_SOURCE_DIR}/Tests/TestTools.cpp)
target_link_libraries (${TEST_TOOLS} DropTools ${LIBRARIES})
add_test (${TEST_TOOLS} ${TEST_TOOLS})
