# Find RapidJSON
#
# This module makes use of the following variable:
#
# RAPIDJSON_INCLUDE_DIR   - Path to RapidJSON include directory
#
# This module will define the following variables:
#
# RAPIDJSON_FOUND         - True if RapidJSON was found
# RAPIDJSON_INCLUDE_DIRS  - Path to RapidJSON include directory
# RAPIDJSON_VERSION       - The version of RapidJSON that was found
#
# Copyright (C) 2018 by Matt Swain <m.swain@me.com>
# Redistribution and use is allowed according to the terms of the BSD license.

if(RAPIDJSON_INCLUDE_DIRS)

  # Already in cache
  set(RAPIDJSON_FOUND TRUE)

else()

  # Use pkg-config to get the include dir and use it as a hint in find_path
  find_package(PkgConfig)
  if(PKG_CONFIG_FOUND)
    pkg_check_modules(PC_RAPIDJSON RapidJSON>=1.1.0 QUIET)
  endif()

  # Find include path with hints from manually specified RAPIDJSON_INCLUDE_DIR and pkg-config
  find_path(
    RAPIDJSON_INCLUDE_DIRS
    NAMES rapidjson/rapidjson.h
    HINTS ${RAPIDJSON_INCLUDE_DIR} ${PC_RAPIDJSON_INCLUDEDIR}
    DOC "Include directory for the RapidJSON library."
  )
  mark_as_advanced(RAPIDJSON_INCLUDE_DIRS)

  # Get the version from pkg-config otherwise parse from the header file
  if(PC_RAPIDJSON_VERSION)
    set(RAPIDJSON_VERSION ${PC_RAPIDJSON_VERSION})
  elseif(EXISTS ${RAPIDJSON_INCLUDE_DIRS}/rapidjson/rapidjson.h)
    file(STRINGS ${RAPIDJSON_INCLUDE_DIRS}/rapidjson/rapidjson.h DEFINES REGEX "#define RAPIDJSON_(MAJOR|MINOR|PATCH)_VERSION ([0-9]+)")
    string(REGEX REPLACE ".+([0-9]+).+([0-9]+).+([0-9]+)" "\\1.\\2.\\3" RAPIDJSON_VERSION "${DEFINES}")
  endif()
  mark_as_advanced(RAPIDJSON_VERSION)

  find_package_handle_standard_args(RapidJSON REQUIRED_VARS RAPIDJSON_INCLUDE_DIRS VERSION_VAR RAPIDJSON_VERSION)
endif()
