# Check support for native CPU architecture optimizations
#
# This module will define the following variables:
#
# NATIVE_FLAGS - compiler flags to enable supported native optimizations
# HAVE_SSE2    - SSE2 supported
# HAVE_SSE4_2  - SSE4.2 supported
# HAVE_AVX     - AVX supported
# HAVE_AVX2    - AVX2 supported
# HAVE_NEON    - ARM NEON supported
#
# Copyright (C) 2018 by Matt Swain <m.swain@me.com>
# Redistribution and use is allowed according to the terms of the BSD license.

include(CheckCXXSourceRuns)

if(MSVC)
  set(SSE2_FLAG "/arch:SSE2")
  set(SSE4_2_FLAG "/arch:SSE2")  # /arch:SSE4 flag doesn't exist?
  set(AVX_FLAG "/arch:AVX")
  set(AVX2_FLAG "/arch:AVX2")
else()
  set(SSE2_FLAG "-msse2")
  set(SSE4_2_FLAG "-msse4.2")
  set(AVX_FLAG "-mavx")
  set(AVX2_FLAG "-mavx2")
endif()

# Store current CMAKE_REQUIRED_FLAGS so we can reset to this afterwards
set(PREVIOUS_CMAKE_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS})

set(CMAKE_REQUIRED_FLAGS ${SSE2_FLAG})
check_cxx_source_runs("
  #include <emmintrin.h>
  int main() {
    __m128d a, b;
    double vals[2] = {0};
    a = _mm_loadu_pd(vals);
    b = _mm_add_pd(a, a);
    _mm_storeu_pd(vals, b);
    return 0;
  }"
  HAVE_SSE2)

set(CMAKE_REQUIRED_FLAGS ${SSE4_2_FLAG})
check_cxx_source_runs("
  #include <emmintrin.h>
  #include <nmmintrin.h>
  int main() {
    long long a[2] = { 1, 2};
    long long b[2] = {-1, 3};
    long long c[2];
    __m128i va = _mm_loadu_si128((__m128i*)a);
    __m128i vb = _mm_loadu_si128((__m128i*)b);
    __m128i vc = _mm_cmpgt_epi64(va, vb);
    _mm_storeu_si128((__m128i*)c, vc);
    if (c[0] == -1LL && c[1] == 0LL)
      return 0;
    else
      return 1;
  }"
  HAVE_SSE4_2)

set(CMAKE_REQUIRED_FLAGS ${AVX_FLAG})
check_cxx_source_runs("
  #include <immintrin.h>
  int main()
  {
    __m256 a;
    a = _mm256_set1_ps(0);
    return 0;
  }"
  HAVE_AVX)

set(CMAKE_REQUIRED_FLAGS ${AVX2_FLAG})
check_cxx_source_runs("
  #include <immintrin.h>
  int main()
  {
    __m256i a = {0};
    a = _mm256_abs_epi16(a);
    return 0;
  }"
  HAVE_AVX2)

# Reset CMAKE_REQUIRED_FLAGS
set(CMAKE_REQUIRED_FLAGS ${PREVIOUS_CMAKE_REQUIRED_FLAGS})

# Check if __ARM_NEON is defined  TODO: Not sure this works
CHECK_SYMBOL_EXISTS(__ARM_NEON "" HAVE_NEON) 

# Set native flags to all supported
if(NOT MSVC)
  # For gcc, clang, intel, just use -march=native to automatically enable all supported flags
  set(NATIVE_FLAGS "-march=native")
  # -msse2 -msse4.2 -mavx -mavx2
else()
  # For MSVC, enable flags individually (Don't set /arch:SSE2 if CMAKE_CL_64 - it is on by default)
  if(HAVE_AVX2)
    set(NATIVE_FLAGS "${AVX2_FLAG}")
  elseif(HAVE_AVX)
    set(NATIVE_FLAGS "${AVX_FLAG}")
  elseif(HAVE_SSE4_2 AND NOT CMAKE_CL_64)
    set(NATIVE_FLAGS "${SSE4_2_FLAG}")
  elseif(HAVE_SSE2 AND NOT CMAKE_CL_64)
    SET(NATIVE_FLAGS "${SSE2_FLAG}")
  endif()
endif()

mark_as_advanced(NATIVE_FLAGS HAVE_SSE2 HAVE_SSE4_2 HAVE_AVX HAVE_AVX2 HAVE_NEON)
