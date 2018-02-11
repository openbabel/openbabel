/**********************************************************************
json.h -

Copyright (C) 2018 by Matt Swain

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#ifndef OPENBABEL_JSON_H
#define OPENBABEL_JSON_H

#include <rapidjson/document.h>
#include <rapidjson/error/en.h>
#include <rapidjson/istreamwrapper.h>
#include <rapidjson/ostreamwrapper.h>
#include <rapidjson/writer.h>
#include <rapidjson/prettywriter.h>


// These should already be defined by CMake
#if defined(OPTIMIZE_NATIVE)
  #if !defined(RAPIDJSON_SSE42) && defined(__SSE4_2__)
    #define RAPIDJSON_SSE42
  #elif !defined(RAPIDJSON_SSE2) && defined(__SSE2__)
    #define RAPIDJSON_SSE2
  #elif!defined(RAPIDJSON_NEON) &&  defined(__ARM_NEON)
    #define RAPIDJSON_NEON
  #endif
#endif

#endif //OPENBABEL_JSON_H
