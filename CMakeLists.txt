cmake_minimum_required(VERSION 3.27)
project(Math)
include(FetchContent)
include(GoogleTest)
enable_testing()
include(CTest)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
FetchContent_declare(
    googletest
    GIT_REPOSITORY  https://github.com/google/googletest.git
    GIT_TAG         main
)
FetchContent_makeAvailable(googletest)
file(GLOB SOURCES "src/*.cpp")
file(GLOB TESTS "test/*.cpp")
set(INCLUDE "include")
include_directories(${INCLUDE})
add_library(Math STATIC ${SOURCES})
add_executable(tests ${TESTS} ${SOURCES})
target_include_directories(tests PUBLIC ${googletest_SOURCE_DIR}/googletest/include  ${gtest_SOURCE_DIR})
target_link_libraries(tests gtest_main Math)
gtest_add_tests(TARGET tests)
add_test(NAME tests COMMAND tests)


