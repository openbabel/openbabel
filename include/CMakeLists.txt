
file(GLOB openbabel_headers openbabel/*.h)
file(GLOB math_headers      openbabel/math/*.h)
file(GLOB stereo_headers      openbabel/stereo/*.h)

install(FILES inchi_api.h          DESTINATION include/inchi)
install(FILES chemdrawcdx.h        DESTINATION include/openbabel-2.0/openbabel)
install(FILES ${openbabel_headers} DESTINATION include/openbabel-2.0/openbabel)
install(FILES ${math_headers}      DESTINATION include/openbabel-2.0/openbabel/math)
install(FILES ${stereo_headers}      DESTINATION include/openbabel-2.0/openbabel/stereo)
