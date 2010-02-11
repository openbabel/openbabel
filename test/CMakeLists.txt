# define TESTDATADIR for tests that need input files
add_definitions(-DTESTDATADIR="\\"${CMAKE_SOURCE_DIR}/test/files/\\"")
############################################################
#  new tests using obtest.h
############################################################

set(tests
    uniqueid
    smiles
    stereo
    tetrahedral
    tetranonplanar
    tetraplanar
#    squareplanar
    cistrans
    stereoperception
   )

set(TEST_PATH ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

foreach(test ${tests})
  add_executable(${test}test ${test}test.cpp)
  target_link_libraries(${test}test ${BABEL_LIBRARY})
  if(MSVC)
    add_test(${test}_Test ${CMAKE_BINARY_DIR}/${test}test)
  else(MSVC)
    add_test(${test}_Test ${TEST_PATH}/${test}test)
  endif(MSVC)
  set_tests_properties(${test}_Test PROPERTIES
    FAIL_REGULAR_EXPRESSION "ERROR;FAIL;Test failed")
  if(MSVC)
    set_target_properties(${test}test PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
  endif(MSVC)
endforeach(test ${tests})

 
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
    datatest
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
  if(MSVC)
    add_test(${test}_Test ${CMAKE_BINARY_DIR}/${test})
  else(MSVC)
    add_test(${test}_Test ${TEST_PATH}/${test})
  endif(MSVC)
  set_tests_properties(${test}_Test PROPERTIES
    FAIL_REGULAR_EXPRESSION "not ok")
  if(MSVC)
    set_target_properties(${test} PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
  endif(MSVC)
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
