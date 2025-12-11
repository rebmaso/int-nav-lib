cmake_minimum_required(VERSION 3.7)

project(
    intnavlib
    VERSION 1.0
    DESCRIPTION "integrated navigation c++ library"
)

# ===================== Global Settings =====================

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

option(BUILD_COVERAGE "Enable coverage reporting" OFF)
option(BUILD_TESTS "Build tests" ON)
option(BUILD_VISION "Build vision module" ON)

# ===================== Dependencies ========================

find_package(Eigen3 REQUIRED)

# ===================== Library Sources =====================

file(GLOB PRIVATE_SOURCE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")

add_library(intnavlib STATIC ${PRIVATE_SOURCE_FILES})

target_include_directories(
    intnavlib
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/third_party>
        $<INSTALL_INTERFACE:include>
)

set_target_properties(
    intnavlib
    PROPERTIES
        VERSION ${PROJECT_VERSION}
        SOVERSION ${PROJECT_VERSION_MAJOR}
        POSITION_INDEPENDENT_CODE ON
        LINKER_LANGUAGE CXX
)

target_link_libraries(intnavlib PUBLIC 
                                    Eigen3::Eigen
                                    ${OPENTHREADS_LIBRARY}
                                )

# ============= Optional Vision Build ==========

if(BUILD_VISION)

target_compile_definitions(intnavlib PUBLIC BUILD_VISION)

find_package(glog REQUIRED)
find_package(OpenCV 4 REQUIRED)
find_package(OsgEarth REQUIRED)
set(OpenGL_GL_PREFERENCE LEGACY)
find_package(OpenGL REQUIRED)
find_package(EGL REQUIRED)
find_package(OSG REQUIRED COMPONENTS osg osgDB osgGA osgUtil osgViewer)
find_package(Fontconfig)

target_link_libraries(intnavlib PUBLIC
                                    Eigen3::Eigen
                                    ${OPENTHREADS_LIBRARY}
                                    glog::glog
                                    ${OpenCV_LIBS}
                                    ${EGL_LIBRARY}
                                    ${OSG_LIBRARY}
                                    ${OSGVIEWER_LIBRARY}
                                    ${OSGDB_LIBRARY}
                                    ${OSGGA_LIBRARY}
                                    ${OSGUTIL_LIBRARY}
                                    ${OSGEARTH_LIBRARY}
                                )

endif()

# ===================== Testing =================

if(BUILD_TESTS)

    # Enable compiler coverage on test targets
    function(enable_coverage target)
        if(BUILD_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
            target_compile_options(${target} PRIVATE -O0 -g --coverage)
            target_link_options(${target} PRIVATE --coverage)
        endif()
    endfunction()

    find_package(glog REQUIRED)

    include(CTest)
    enable_testing()

    # GTest: find or fetch
    find_package(GTest)
    if(NOT GTest_FOUND)
        include(FetchContent)
        FetchContent_Declare(
            googletest
            URL https://github.com/google/googletest/archive/5376968f6948923e2411081fd9372e71a59d8e77.zip
        )
        set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
        FetchContent_MakeAvailable(googletest)
    endif()

    # Collect test sources
    file(GLOB TEST_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp")

    # Create each test
    foreach(test_file ${TEST_SOURCES})
        get_filename_component(test_name ${test_file} NAME_WE)

        add_executable(${test_name} ${test_file})

        target_link_libraries(
            ${test_name}
            PRIVATE
                gtest
                gtest_main
                glog::glog
                intnavlib
        )

        add_test(NAME ${test_name} COMMAND ${test_name})
        enable_coverage(${test_name})
    endforeach()

endif()

# ===================== Install ==============================

include(GNUInstallDirs)

install(
    TARGETS intnavlib
    EXPORT intnavlibTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Install public headers
install(DIRECTORY include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/intnavlib
        FILES_MATCHING PATTERN "*.h")

install(DIRECTORY third_party/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/intnavlib
        FILES_MATCHING PATTERN "*.h")

# Install the export file
install(
    EXPORT intnavlibTargets
    FILE intnavlibTargets.cmake
    NAMESPACE intnavlib::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/intnavlib
)

# CMake package configuration 
include(CMakePackageConfigHelpers)

# Write version and config
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/intnavlibConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/intnavlibConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/intnavlibConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/intnavlib
)

install(
    FILES
        "${CMAKE_CURRENT_BINARY_DIR}/intnavlibConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/intnavlibConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/intnavlib
)
