cmake_minimum_required(VERSION 3.25)
project(ImageProcessingTestTask LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(BUILD_TESTS "Build tests" OFF)

find_package(Qt6 COMPONENTS Core Widgets Gui REQUIRED)
qt_standard_project_setup()

file(GLOB LIB_SOURCES src/*.cpp)
file(GLOB LIB_HEADERS include/*.hpp)

add_library(ImageProcessingLib ${LIB_SOURCES} ${LIB_HEADERS})

target_include_directories(ImageProcessingLib PUBLIC
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_BINARY_DIR}/ImageProcessingTestTask_autogen/include
)

target_link_libraries(ImageProcessingLib
    Qt::Core Qt::Widgets Qt::Gui
)

add_executable(ImageProcessingTestTask src/main.cpp)
target_link_libraries(ImageProcessingTestTask ImageProcessingLib)

if (WIN32)
    get_target_property(windeployqt_exe Qt6::windeployqt LOCATION)
    add_custom_command(TARGET ImageProcessingTestTask POST_BUILD
        COMMAND "${windeployqt_exe}"
        $<TARGET_FILE:ImageProcessingTestTask>
    )
endif()

if (BUILD_TESTS)
    find_package(GTest REQUIRED)
    file(GLOB TEST_SOURCES tests/*.cpp)

    add_executable(Tests ${TEST_SOURCES})
    target_include_directories(Tests PRIVATE
        ${CMAKE_SOURCE_DIR}/include
    )
    target_link_libraries(Tests
        GTest::gtest_main
    )

    enable_testing()
    gtest_discover_tests(Tests)
endif()
