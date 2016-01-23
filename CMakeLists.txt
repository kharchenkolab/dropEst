cmake_minimum_required(VERSION 3.2)
project(cp)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/CMake")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wno-deprecated -m64 -DBOOST_LOG_DYN_LINK")

set(BOOST_ROOT "/home/victor/libs/boost_1_58_0/")
find_package(Boost COMPONENTS unit_test_framework thread system log iostreams REQUIRED)
find_package (Threads REQUIRED)

include_directories(${Boost_INCLUDE_DIR} ${PROJECT_SOURCE_DIR})

FILE(GLOB TagsSearchSources TagsSearch/*.cpp TagsSearch/Counters/*.cpp)
add_library(TagsSearch ${TagsSearchSources})

set(LIBRARIES ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARIES})

#set_source_files_properties(indropest.cpp PROPERTIES COMPILE_FLAGS
#        "-lz -lm -lbam -lboost_iostreams -lboost_serialization -std=c++11")
#set_source_files_properties(indropest.cpp PROPERTIES COMPILE_FLAGS   "-lz -lm -lbam -lboost_iostreams -lboost_serialization -std=c++11 -I /home/victor/libs/samtools/trunk")

add_executable(indroptag indroptag.cpp)
target_link_libraries(indroptag  TagsSearch ${LIBRARIES})

set (TEST TestTagsSearch)
add_executable (${TEST} ${PROJECT_SOURCE_DIR}/Tests/TestTagsSearch.cpp)
target_link_libraries (${TEST} TagsSearch ${LIBRARIES})
enable_testing ()
add_test (${TEST} ${TEST})