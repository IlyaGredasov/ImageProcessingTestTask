cmake_minimum_required(VERSION 3.25)
project(ImageProcessingTestTask)

set(CMAKE_CXX_STANDARD 17)

find_package(Qt6 COMPONENTS Core Widgets Gui REQUIRED)
qt_standard_project_setup()
file(GLOB SOURCES "${CMAKE_SOURCE_DIR}/src/*.cpp")
file(GLOB HEADERS "${CMAKE_SOURCE_DIR}/include/*.hpp")
file(GLOB UIS "${CMAKE_SOURCE_DIR}/src/*.ui")
add_executable(ImageProcessingTestTask ${SOURCES} ${HEADERS} ${UIS})
target_include_directories(ImageProcessingTestTask PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_BINARY_DIR}/ImageProcessingTestTask_autogen/include)
target_link_libraries(ImageProcessingTestTask Qt::Core Qt::Widgets Qt::Gui)
if (WIN32)
    get_target_property(windeployqt_exe Qt6::windeployqt LOCATION)
    add_custom_command(TARGET ImageProcessingTestTask POST_BUILD
            COMMAND "${windeployqt_exe}" --qmldir ${CMAKE_SOURCE_DIR}
            $<TARGET_FILE:ImageProcessingTestTask>
            COMMENT "Running windeployqt")
endif ()

find_package(GTest REQUIRED)
file(GLOB TESTS "${CMAKE_SOURCE_DIR}/tests/*.cpp")
add_executable(Tests ${TESTS})
enable_testing()
target_link_libraries(Tests GTest::gtest_main)
target_include_directories(Tests PRIVATE ${CMAKE_SOURCE_DIR}/include)
gtest_discover_tests(Tests)