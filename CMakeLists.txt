cmake_minimum_required(VERSION 3.20)

project(Fermi VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

file(GLOB SOURCE_FILES src/*.cc)
list(FILTER SOURCE_FILES EXCLUDE REGEX ".*_test\\.cc")

add_library(fermi STATIC ${SOURCE_FILES})
target_link_libraries(fermi PRIVATE sparsex)

set(SPARSEX_DISABLE_TESTS ON)
add_subdirectory(lib/sparsex)

target_include_directories(fermi PUBLIC src lib)

file(GLOB EXEC_FILES "src/main/main_*.cc")
foreach(EXEC_FILE ${EXEC_FILES})
    get_filename_component(EXEC_NAME ${EXEC_FILE} NAME_WE)
    add_executable(${EXEC_NAME} ${EXEC_FILE})
    target_link_libraries(${EXEC_NAME} PRIVATE fermi sparsex)
endforeach()

file(GLOB TEST_FILES "src/*_test.cc")
add_subdirectory(lib/googletest)
include(GoogleTest)

foreach(TEST_FILE ${TEST_FILES})
    get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)
    add_executable(${TEST_NAME} ${TEST_FILE})
    target_link_libraries(${TEST_NAME} fermi gtest gtest_main pthread)
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
	gtest_discover_tests(${TEST_NAME})
endforeach()

enable_testing()
