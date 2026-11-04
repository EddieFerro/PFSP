cmake_minimum_required(VERSION 3.20)
project(pfsp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

option(PFSP_NATIVE   "Compile with -march=native (only if you run where you build)" OFF)
option(PFSP_SANITIZE "Build with AddressSanitizer + UndefinedBehaviorSanitizer"       OFF)
# option(PFSP_TESTS    "Build the test suite (run with ctest)"                          ON)
# option(PFSP_BENCH    "Build the benchmark input simulator"                            ON)

find_package(ZLIB REQUIRED)

# -------------------------------------------------------------------------
# Third-party sources: pinned versions
# -------------------------------------------------------------------------
include(FetchContent)

# CLI11 (single-header release asset)
FetchContent_Declare(cli11
  URL      https://github.com/CLIUtils/CLI11/releases/download/v2.4.2/CLI11.hpp
  DOWNLOAD_NO_EXTRACT TRUE
)

# kseq.h (from seqtk: 64-bit lengths, memchr line scan)
FetchContent_Declare(kseq
  URL      https://raw.githubusercontent.com/lh3/seqtk/94e707082d39b0a038f234df676e32d9802c0dc7/kseq.h
  DOWNLOAD_NO_EXTRACT TRUE
)

# MurmurHash3 (smhasher)
set(SMHASHER_RAW https://raw.githubusercontent.com/aappleby/smhasher/07bb4de10a63e8cc2e1724865454eba635742383/src)
FetchContent_Declare(murmur3_h
  URL      ${SMHASHER_RAW}/MurmurHash3.h
  DOWNLOAD_NO_EXTRACT TRUE
)
FetchContent_Declare(murmur3_cpp
  URL      ${SMHASHER_RAW}/MurmurHash3.cpp
  DOWNLOAD_NO_EXTRACT TRUE
)

FetchContent_MakeAvailable(cli11 kseq murmur3_h murmur3_cpp)

# Gather the single files into one include directory in the build tree.
set(PFSP_THIRD_PARTY_DIR ${CMAKE_BINARY_DIR}/third_party)
file(MAKE_DIRECTORY ${PFSP_THIRD_PARTY_DIR})
file(COPY ${cli11_SOURCE_DIR}/CLI11.hpp
          ${kseq_SOURCE_DIR}/kseq.h
          ${murmur3_h_SOURCE_DIR}/MurmurHash3.h
          ${murmur3_cpp_SOURCE_DIR}/MurmurHash3.cpp
     DESTINATION ${PFSP_THIRD_PARTY_DIR})

# -------------------------------------------------------------------------
# Common compile settings
# -------------------------------------------------------------------------
add_library(pfsp_options INTERFACE)
target_include_directories(pfsp_options INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${PFSP_THIRD_PARTY_DIR}
)
target_compile_options(pfsp_options INTERFACE -Wall)
if(PFSP_NATIVE)
  target_compile_options(pfsp_options INTERFACE -march=native)
endif()
if(PFSP_SANITIZE)
  target_compile_options(pfsp_options INTERFACE
    -O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined -fno-sanitize-recover=all)
  target_link_options(pfsp_options INTERFACE -fsanitize=address,undefined)
endif()

# MurmurHash3: third-party, compiled without warnings
add_library(murmurhash3 STATIC ${PFSP_THIRD_PARTY_DIR}/MurmurHash3.cpp)
target_include_directories(murmurhash3 PUBLIC ${PFSP_THIRD_PARTY_DIR})
target_compile_options(murmurhash3 PRIVATE -w)
if(PFSP_SANITIZE)
  target_compile_options(murmurhash3 PRIVATE -fsanitize=address,undefined)
endif()

# -------------------------------------------------------------------------
# pfsp
# -------------------------------------------------------------------------
add_executable(pfsp src/main.cpp)
target_link_libraries(pfsp PRIVATE pfsp_options murmurhash3 ZLIB::ZLIB)

# -------------------------------------------------------------------------
# Tests: cmake --build build && ctest --test-dir build --output-on-failure
# -------------------------------------------------------------------------
# if(PFSP_TESTS)
#   enable_testing()
#   set(PFSP_TEST_NAMES
#     test_alphabet test_scan test_sort test_parse test_phrase test_match test_e2e)
#   foreach(t IN LISTS PFSP_TEST_NAMES)
#     add_executable(${t} tests/${t}.cpp)
#     target_link_libraries(${t} PRIVATE pfsp_options murmurhash3)
#     add_test(NAME ${t} COMMAND ${t} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
#   endforeach()
#   set_tests_properties(test_e2e PROPERTIES ENVIRONMENT "PFSP_BIN=$<TARGET_FILE:pfsp>")
# endif()

# -------------------------------------------------------------------------
# Benchmarks: python3 bench/run_bench.py --pfsp build/pfsp --simulate build/simulate
# -------------------------------------------------------------------------
# if(PFSP_BENCH)
#   add_executable(simulate bench/simulate.cpp)
#   target_compile_options(simulate PRIVATE -Wall)
# endif()
