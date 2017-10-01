cmake_minimum_required(VERSION 2.8)

project(dropEstR)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/../CMake")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -std=c++11 -Wno-attributes -Wno-deprecated -Wno-deprecated-declarations")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -ggdb -gdwarf-2 -g3")

if (NOT "${CMAKE_BUILD_TYPE}" STREQUAL Debug)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")
endif (NOT "${CMAKE_BUILD_TYPE}" STREQUAL Debug)

#set(FIND_RCPP_PROGRESS TRUE)
find_package (R REQUIRED)

set(INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/inst/include ${R_INCLUDE_DIRS})
include_directories(${INCLUDE_DIRS})

set(BASE_LIBRARIES ${R_LIBRARIES})

FILE(GLOB DropEstRSources src/*.cpp)
add_library(dropEstR ${DropEstRSources})
target_link_libraries(dropEstR ${BASE_LIBRARIES})