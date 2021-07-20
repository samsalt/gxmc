
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

project(gxc)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g2 -ggdb -fopenmp")  
SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall -fopenmp")
SET(CMAKE_BUILD_TYPE "Debug")
SET(CMAKE_BUILD_TYPE "Release")

if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
    message("Debug mode")
    add_definitions(-DDEBUG)
else()
    message("Release mode")
endif ()


SET(EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/bin") 
SET(LOCAL_LIB_PATH "${PROJECT_SOURCE_DIR}/lib") 

INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/include") 

AUX_SOURCE_DIRECTORY(src DIR_SRCS) 

# link_directories(~/lib_local/seacas/lib)
# add_library(${LOCAL_LIB_PATH}/libexodusgxc.a STATIC)

add_executable(gxc MACOSX_BUNDLE ${DIR_SRCS})
# target_link_libraries(${PROJECT_NAME} ${LOCAL_LIB_PATH}/libexodus.so.2)
target_link_libraries(${PROJECT_NAME} /usr/lib/x86_64-linux-gnu/libexodus.so.2)
