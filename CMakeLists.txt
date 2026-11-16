# kinesis: dynamical systems as objects of study and of inference — a language in which they are declared, what a
# declared law does over time, and laws recovered from records.
#
# Every module is a header-only target kinesis::<module>. kinesis uses the targets of calx, morphe, physis and tyche:
# inside the infrastructure, the siblings under foundations/ (each added here when no larger project has added it).
# Built on its own, kinesis builds its checks (KINESIS_CHECKS); added to a larger project, only when asked.

cmake_minimum_required(VERSION 3.20)
project(kinesis VERSION 0.1.0 LANGUAGES CXX)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(KINESIS_TOP_LEVEL ON)
else()
  set(KINESIS_TOP_LEVEL OFF)
endif()
option(KINESIS_CHECKS "Build kinesis's checks" ${KINESIS_TOP_LEVEL})

if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 20)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
endif()
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

# kinesis_sibling(<directory> <a target it defines>): the sibling library, added once.
function(kinesis_sibling directory target)
  if(NOT TARGET ${target})
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../${directory}/CMakeLists.txt)
      add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../${directory} ${CMAKE_BINARY_DIR}/${directory})
    else()
      message(FATAL_ERROR "kinesis needs ${directory}: build it inside the infrastructure (foundations/${directory}).")
    endif()
  endif()
endfunction()

kinesis_sibling(calx calx::linalg)
kinesis_sibling(morphe morphe::topology)
kinesis_sibling(physis physis::dynamics)
kinesis_sibling(tyche tyche::probability)
if(KINESIS_CHECKS)
  enable_testing()
endif()

# kinesis_module(<name> [LIBRARIES ...] [CHECK_LIBRARIES ...]): the target kinesis::<name> and, with the checks,
# check_kinesis_<name>, which may use more libraries than the module does.
function(kinesis_module name)
  cmake_parse_arguments(ARG "" "" "LIBRARIES;CHECK_LIBRARIES" ${ARGN})
  add_library(kinesis_${name} INTERFACE)
  add_library(kinesis::${name} ALIAS kinesis_${name})
  target_include_directories(kinesis_${name} INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/${name}/include)
  target_compile_features(kinesis_${name} INTERFACE cxx_std_20)
  target_link_libraries(kinesis_${name} INTERFACE ${ARG_LIBRARIES})
  if(KINESIS_CHECKS AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${name}/tests/check_${name}.cpp)
    calx_add_check(check_kinesis_${name} SOURCES ${name}/tests/check_${name}.cpp LIBRARIES kinesis::${name} ${ARG_CHECK_LIBRARIES})
  endif()
endfunction()

# The modules in their order (decision 0003).
kinesis_module(kernel)
kinesis_module(terms LIBRARIES kinesis::kernel)
kinesis_module(language LIBRARIES kinesis::terms morphe::formal)
kinesis_module(systems LIBRARIES kinesis::language physis::dynamics physis::quantum calx::linalg calx::random calx::integrators morphe::formal)
kinesis_module(analysis LIBRARIES kinesis::systems calx::linalg calx::random tyche::descriptive tyche::time_series)
kinesis_module(braids LIBRARIES kinesis::kernel morphe::algebra morphe::arithmetic calx::linalg)
kinesis_module(induction LIBRARIES kinesis::systems calx::linalg calx::random)
kinesis_module(coupling LIBRARIES kinesis::systems CHECK_LIBRARIES kinesis::analysis)
kinesis_module(catalogue LIBRARIES kinesis::language CHECK_LIBRARIES kinesis::systems kinesis::analysis kinesis::induction)
kinesis_module(balance LIBRARIES kinesis::induction calx::linalg CHECK_LIBRARIES calx::random)
