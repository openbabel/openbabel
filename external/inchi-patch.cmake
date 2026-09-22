# InChI does not have install targets
file(APPEND "${ORIG}" "
### OPENBABEL PATCH ###
install(FILES
  \"\${PROJECT_SOURCE_DIR}/../../../INCHI_BASE/src/bcf_s.h\"
  \"\${PROJECT_SOURCE_DIR}/../../../INCHI_BASE/src/inchi_api.h\"
  \"\${PROJECT_SOURCE_DIR}/../../../INCHI_BASE/src/ixa.h\"
  DESTINATION include/inchi
)
install(TARGETS libinchi LIBRARY DESTINATION ${LIB_INSTALL_DIR})")
