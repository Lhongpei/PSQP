cmake_minimum_required(VERSION 3.20)
project(PSQP LANGUAGES C CXX) # QP Presolver, forked from PSLP

set(PSQP_VERSION_MAJOR 0)
set(PSQP_VERSION_MINOR 1)
set(PSQP_VERSION_PATCH 0)
set(PSQP_VERSION "${PSQP_VERSION_MAJOR}.${PSQP_VERSION_MINOR}.${PSQP_VERSION_PATCH}")
add_compile_definitions(PSQP_VERSION="${PSQP_VERSION}")

message(STATUS "Configuring PSQP (version ${PSQP_VERSION}) - QP Presolver forked from PSLP")

include(GNUInstallDirs)

# Build position - independent code
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Disable C compiler extensions
set(CMAKE_C_EXTENSIONS OFF)

# This is needed for clock_gettime on Linux without compiler extensions
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    add_compile_definitions(_POSIX_C_SOURCE=200809L)
endif()

# Append custom CMake modules
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Option for shared / static libraries
option(BUILD_SHARED_LIBS "Build libraries as shared as opposed to static" ON)

# Enable RPATH support
include(AddInstallRPATHSupport)
add_install_rpath_support(
  BIN_DIRS "${CMAKE_INSTALL_FULL_BINDIR}"
  LIB_DIRS "${CMAKE_INSTALL_FULL_LIBDIR}"
  INSTALL_NAME_DIR "${CMAKE_INSTALL_FULL_LIBDIR}"
  USE_LINK_PATH
)

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Add uninstall target
include(AddUninstallTarget)

# Specify the C standard
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

#Collect source files
file(GLOB CORE_SOURCES CONFIGURE_DEPENDS
    src/core/*.c
    src/explorers/*.c
)

# distinguish between public and private headers
set(PUBLIC_HEADERS
  include/PSQP/PSQP_API.h
 
)

file(GLOB PRIVATE_HEADERS CONFIGURE_DEPENDS
    include/core/*.h
    include/data_structures/*.h
    include/explorers/*.h
)


# Handle debugger conditionally
option(ENABLE_DEBUGGER "Enable debugger functionality" OFF)

if(NOT (ENABLE_DEBUGGER OR CMAKE_BUILD_TYPE STREQUAL "Debug"))
    list(REMOVE_ITEM CORE_SOURCES 
        "${CMAKE_CURRENT_SOURCE_DIR}/src/core/Debugger.c"
    )
    message(STATUS "Debugger disabled")
else()
    message(STATUS "Debugger enabled")
endif()

# Create PSQP library
add_library(PSQP ${CORE_SOURCES})
target_include_directories(PSQP
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/PSQP>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/PSQP>
  PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/core>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/data_structures>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/explorers>
)

# Common compile options
if (MSVC)
  target_compile_options(PSQP PRIVATE /W4)
  set(CMAKE_DEBUG_POSTFIX "d")
  set_target_properties(PSQP PROPERTIES DEBUG_POSTFIX "d")
  set_target_properties(PSQP PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
  target_compile_definitions(PSQP PRIVATE PSLP_EXPORTS)
  target_compile_options(PSQP PRIVATE
    $<$<CONFIG:Debug>:/Od>
    $<$<CONFIG:Release>:/O2 /DNDEBUG>
    $<$<CONFIG:RelWithDebInfo>:/O2 /DNDEBUG>
    $<$<CONFIG:MinSizeRel>:/O1 /DNDEBUG>
  )
else()
  target_compile_options(PSQP PRIVATE -Wall -Wextra -Wpedantic -Wno-sign-compare -Wconversion)
  target_compile_options(PSQP PRIVATE
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>:-O3 -DNDEBUG>
    $<$<CONFIG:RelWithDebInfo>:-O2 -g -DNDEBUG>
    $<$<CONFIG:MinSizeRel>:-Os -DNDEBUG>
  )

endif()

if(ENABLE_DEBUGGER OR CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_definitions(PSQP PRIVATE DEBUGGER_ENABLED)
else()
    target_compile_definitions(PSQP PRIVATE NDEBUG)
endif()

# Link pthread library or compatibility layer
if(NOT WIN32)
    find_package(Threads REQUIRED)
    target_link_libraries(PSQP PRIVATE Threads::Threads)
endif()

# Testing
option(PSQP_BUILD_TESTING "Build PSLP tests" OFF)
if(PSQP_BUILD_TESTING)
  include(CTest)
  enable_testing()
  message(STATUS "PSQP testing enabled (this should be disabled for production builds!!!)")

  add_executable(all_tests tests/all_tests.c tests/test_macros.c)

  # Link PSQP library
  target_link_libraries(all_tests PRIVATE PSQP)
  target_compile_definitions(PSQP PRIVATE TESTING)

    # Add include directories explicitly for the test target
  target_include_directories(all_tests PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/include/core
      ${CMAKE_CURRENT_SOURCE_DIR}/include/data_structures
      ${CMAKE_CURRENT_SOURCE_DIR}/include/explorers
      ${CMAKE_CURRENT_SOURCE_DIR}/include/PSQP
  )

  target_compile_definitions(all_tests PRIVATE TESTING)

  add_test(NAME all_tests COMMAND all_tests)
else()
  message(STATUS "Testing disabled")
endif()

# Link math library on UNIX
if(UNIX)
  target_link_libraries(PSQP PRIVATE m)
endif()

option(COVERAGE "Enable coverage reporting" OFF)

if(COVERAGE)
    message(STATUS "Building with coverage flags")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -g --coverage")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} --coverage")
endif()

add_executable(test_cpp EXCLUDE_FROM_ALL tests/test_cpp.cpp)
target_link_libraries(test_cpp PRIVATE PSQP)

# QR example executable
add_executable(qr_example EXCLUDE_FROM_ALL tests/qr_example.c)
target_link_libraries(qr_example PRIVATE PSQP)
target_include_directories(qr_example PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include/core
    ${CMAKE_CURRENT_SOURCE_DIR}/include/data_structures
    ${CMAKE_CURRENT_SOURCE_DIR}/include/explorers
    ${CMAKE_CURRENT_SOURCE_DIR}/include/PSLP
)

# Test for advanced QR optimizations
add_executable(test_qr_advanced tests/test_qr_advanced.c)
target_link_libraries(test_qr_advanced PRIVATE PSQP)
target_include_directories(test_qr_advanced PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include/core
    ${CMAKE_CURRENT_SOURCE_DIR}/include/data_structures
    ${CMAKE_CURRENT_SOURCE_DIR}/include/explorers
    ${CMAKE_CURRENT_SOURCE_DIR}/include/PSLP
)

# Test for DtonsEq with pure Q diagonal
add_executable(test_qr_dton tests/test_qr_dton.c)
target_link_libraries(test_qr_dton PRIVATE PSQP)
target_include_directories(test_qr_dton PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include/core
    ${CMAKE_CURRENT_SOURCE_DIR}/include/data_structures
    ${CMAKE_CURRENT_SOURCE_DIR}/include/explorers
    ${CMAKE_CURRENT_SOURCE_DIR}/include/PSLP
)

# QP Debug Test
add_executable(test_qp_debug tests/test_psqp_all/test_qp_debug.c)
target_link_libraries(test_qp_debug PRIVATE PSQP)
target_include_directories(test_qp_debug PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include/core
    ${CMAKE_CURRENT_SOURCE_DIR}/include/data_structures
    ${CMAKE_CURRENT_SOURCE_DIR}/include/explorers
    ${CMAKE_CURRENT_SOURCE_DIR}/include/PSLP
)

# Comprehensive PSQP Test Suite
add_executable(test_psqp_all
    tests/test_psqp_all/test_main.c
    tests/test_psqp_all/test_qr_core.c
    tests/test_psqp_all/test_problem.c
    tests/test_psqp_all/test_presolve.c
    tests/test_psqp_all/test_qp_specific.c
    tests/test_psqp_all/test_qp_step4.c
    tests/test_psqp_all/test_qp_step5.c
    tests/test_psqp_all/test_qp_step6.c
    tests/test_psqp_all/test_qp_step7.c
    tests/test_psqp_all/test_qp_step8.c
    tests/test_psqp_all/test_qp_full.c
    tests/test_psqp_all/test_edge_cases.c
    tests/test_psqp_all/test_edge_full.c
    tests/test_psqp_all/test_qp_postsolve.c
    tests/test_psqp_all/test_qp_kkt.c
    tests/test_psqp_all/test_qp_partial.c
)
target_link_libraries(test_psqp_all PRIVATE PSQP m)
target_include_directories(test_psqp_all PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include/core
    ${CMAKE_CURRENT_SOURCE_DIR}/include/data_structures
    ${CMAKE_CURRENT_SOURCE_DIR}/include/explorers
    ${CMAKE_CURRENT_SOURCE_DIR}/include/PSLP
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_psqp_all
)

# Installation and export
install(TARGETS PSQP
  EXPORT PSQPTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/PSQP
)

install(FILES ${PUBLIC_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/PSQP)

# Export the targets file
install(EXPORT PSQPTargets
  FILE PSQPTargets.cmake
  NAMESPACE PSQP::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PSQP
)

# Create and install Config file
include(CMakePackageConfigHelpers)
configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/PSQPConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/PSQPConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PSQP
)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/PSQPConfigVersion.cmake
  VERSION ${PSQP_VERSION}
  COMPATIBILITY AnyNewerVersion
)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/PSQPConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/PSQPConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PSQP
)