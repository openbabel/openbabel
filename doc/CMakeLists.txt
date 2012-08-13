file(GLOB man_1 *.1)
install(FILES ${man_1} DESTINATION share/man/man1)
install(FILES splash.png DESTINATION share/openbabel/${BABEL_VERSION})

OPTION(BUILD_DOCS "Build Open Babel documentation" OFF)
IF(BUILD_DOCS)
    find_package(Doxygen)
    if(NOT DOXYGEN_EXECUTABLE)
      return()
    endif(NOT DOXYGEN_EXECUTABLE)
    IF(DOT)
        SET(HAVE_DOT YES)
    ELSE(DOT)
        SET(HAVE_DOT NO)
    ENDIF(DOT)
    # This processes our Doxyfile.in and substitutes paths to generate
    # a final Doxyfile
    CONFIGURE_FILE(${openbabel_SOURCE_DIR}/doc/Doxyfile.in
        ${openbabel_SOURCE_DIR}/Doxyfile
    )
    # This creates a new target to build documentation.
    # It runs ${DOXYGEN} which is the full path and executable to
    # Doxygen on your system, set by the FindDoxygen.cmake module
    # It runs the final generated Doxyfile against it.
    # The DOT_PATH is substituted into the Doxyfile.
    ADD_CUSTOM_TARGET(docs
        COMMAND ${DOXYGEN_EXECUTABLE} ${openbabel_SOURCE_DIR}/Doxyfile
        WORKING_DIRECTORY ${openbabel_SOURCE_DIR}
    )
ENDIF(BUILD_DOCS)

OPTION(BUILD_EXAMPLES "Build Open Babel examples" OFF)
IF(BUILD_EXAMPLES)
   add_subdirectory(examples)
ENDIF(BUILD_EXAMPLES)
