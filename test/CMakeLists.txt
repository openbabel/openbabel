# define TESTDATADIR for tests that need input files
cmake_policy(SET CMP0005 NEW) # needed for line below (escapes value)
add_definitions( -DTESTDATADIR="${CMAKE_SOURCE_DIR}/test/files/" )

############################################################
#  neew tests using obtest.h
############################################################

############################################################
#  old tests using "not ok"
############################################################

set(old_tests
    aromatest
    atom
    bond
    cansmi
    cmlreadfile
    conversion
    data
    ffghemical
    ffmmff94
    ffuff
    formalcharge
    format
    formula
#    inchiwrite
    internalcoord
    invalidsmarts
    invalidsmiles
    iterators
    logp_psa
    math
    mol
    phmodel
    residue
    ringtest
    smartstest
    smilesmatch
    unitcell
   )
                

foreach(test ${old_tests})
  add_executable(${test} ${test}.cpp)
  target_link_libraries(${test} ${BABEL_LIBRARY})
  add_test(${test}_Test ${CMAKE_BINARY_DIR}/test/${test})
  set_tests_properties(${test}_Test PROPERTIES FAIL_REGULAR_EXPRESSION "not ok")
endforeach(test ${old_tests})

############################################################  
#  roundtrip 
############################################################  
  
add_executable(roundtrip roundtrip.cpp)
target_link_libraries(roundtrip ${BABEL_LIBRARY})
install(TARGETS roundtrip
                RUNTIME DESTINATION bin
                LIBRARY DESTINATION lib
                ARCHIVE DESTINATION lib
)
