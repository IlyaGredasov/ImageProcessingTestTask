cmake_minimum_required(VERSION 3.10)
project(image_processing_test_task)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)

option(BUILD_TESTS "Build tests" OFF)

find_package(Qt6 COMPONENTS Core Widgets Gui REQUIRED)
qt_standard_project_setup()

file(GLOB LIB_SOURCES src/*.cpp)
file(GLOB LIB_HEADERS include/*.hpp)

add_library(image_processing_lib ${LIB_SOURCES} ${LIB_HEADERS})

target_include_directories(image_processing_lib PUBLIC
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_BINARY_DIR}/image_processing_test_task/include
)

target_link_libraries(image_processing_lib Qt::Core Qt::Widgets Qt::Gui)

add_executable(image_processing_test_task src/main.cpp)
target_link_libraries(image_processing_test_task image_processing_lib)

if (WIN32)
    get_target_property(windeployqt_exe Qt6::windeployqt LOCATION)
    add_custom_command(TARGET image_processing_test_task POST_BUILD
        COMMAND "${windeployqt_exe}"
        $<TARGET_FILE:image_processing_test_task>
    )
endif()

if (BUILD_TESTS)
    find_package(GTest REQUIRED)
    file(GLOB TEST_SOURCES tests/*.cpp)

    add_executable(tests ${TEST_SOURCES})
    target_include_directories(tests PRIVATE ${CMAKE_SOURCE_DIR}/include)
    target_link_libraries(tests GTest::gtest_main)

    enable_testing()
    gtest_discover_tests(tests)
endif()
