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
file(GLOB VECTOR_TEST "test/VectorTest.cpp")
file(GLOB MATRIX_TEST "test/MatrixTest.cpp")
set(INCLUDE "include")
include_directories(${INCLUDE})

add_executable(vector_test ${VECTOR_TEST})
target_include_directories(vector_test PUBLIC ${googletest_SOURCE_DIR}/googletest/include ${gtest_SOURCE_DIR})
target_link_libraries(vector_test gtest_main Math)
gtest_add_tests(TARGET vector_test)
add_test(NAME vector_test COMMAND vector_test)

add_executable(matrix_test ${MATRIX_TEST})
target_include_directories(matrix_test PUBLIC ${googletest_SOURCE_DIR}/googletest/include ${gtest_SOURCE_DIR})
target_link_libraries(matrix_test gtest_main Math)
gtest_add_tests(TARGET matrix_test)
add_test(NAME matrix_test COMMAND matrix_test)


