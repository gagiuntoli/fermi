cmake_minimum_required(VERSION 3.20)

project(Fermi VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB SOURCE_FILES src/*.cc)

add_library(FermiLib STATIC ${SOURCE_FILES})

target_include_directories(FermiLib PUBLIC src)

add_executable(main_keff src/main/main_keff.cc)
target_link_libraries(main_keff PRIVATE FermiLib)

file(GLOB TEST_FILES "src/*_test.cc")
include(GoogleTest)

foreach(TEST_FILE ${TEST_FILES})
    get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)
    add_executable(${TEST_NAME} ${TEST_FILE})
    target_link_libraries(${TEST_NAME} FermiLib gtest gtest_main pthread)
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
	gtest_discover_tests(${TEST_NAME})
endforeach()

enable_testing()