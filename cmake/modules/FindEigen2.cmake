# - Try to find Eigen2 lib
# Once done this will define
#
#  EIGEN2_FOUND - system has eigen lib
#  EIGEN2_INCLUDE_DIR - the eigen include directory

# Copyright (c) 2006, 2007 Montel Laurent, <montel@kde.org>
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.

if (EIGEN2_INCLUDE_DIR)

  # in cache already
  set(EIGEN2_FOUND TRUE)
  
  message(STATUS "Eigen2 from chache: ${EIGEN2_INCLUDE_DIR}")

else (EIGEN2_INCLUDE_DIR)

find_path(EIGEN2_INCLUDE_DIR NAMES Eigen/Core
     PATHS
     ${INCLUDE_INSTALL_DIR}/eigen2
     /usr/local/include/eigen2
   )

if(EIGEN2_INCLUDE_DIR)
  set(EIGEN2_FOUND TRUE)
endif(EIGEN2_INCLUDE_DIR)

if(EIGEN2_FOUND)
  message(STATUS "Found Eigen2: ${EIGEN2_INCLUDE_DIR}")
else(EIGEN2_FOUND)
  message(FATAL_ERROR "Could NOT find Eigen2")
endif(EIGEN2_FOUND)

mark_as_advanced(EIGEN2_INCLUDE_DIR)

endif(EIGEN2_INCLUDE_DIR)

