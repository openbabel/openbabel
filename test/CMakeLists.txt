# define TESTDATADIR for tests that need input files
add_definitions(-DTESTDATADIR="\\"${CMAKE_SOURCE_DIR}/test/files/\\"")

# define FORMATDIR for location of format plugin binaries
add_definitions(-DFORMATDIR="\\"${openbabel_BINARY_DIR}/lib/\\"")

###########################################################
#  new tests using obtest.h
############################################################

set(tests
    automorphism
    builder
    canonconsistent
    canonstable
    cistrans
    graphsym
    implicitH
    lssr
    isomorphism
    rotor
    shuffle
    smiles
    spectrophore
    squareplanar
    stereo
    stereoperception
    tetrahedral
    tetranonplanar
    tetraplanar
    uniqueid
  )

if (EIGEN2_FOUND)
set(tests
        align
        ${tests}
)
endif (EIGEN2_FOUND)

set(TEST_PATH ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

foreach(test ${tests})
  add_executable(test_${test} ${test}test.cpp)
  target_link_libraries(test_${test} ${BABEL_LIBRARY})
  add_test(${test}_Test ${TEST_PATH}/test_${test})
  set_tests_properties(${test}_Test PROPERTIES
    FAIL_REGULAR_EXPRESSION "ERROR;FAIL;Test failed"
    ENVIRONMENT "BABEL_LIBDIR=${CMAKE_BINARY_DIR}/lib"
    ENVIRONMENT "BABEL_DATADIR=${CMAKE_SOURCE_DIR}/data")
endforeach(test ${tests})

############################################################
#  old-style tests using "not ok"
############################################################

set(old_tests
    aromatest
    atom
    bond
    cansmi
    charge-mmff94
    charge-gasteiger
    conversion
    datatest
    ffgaff
    ffghemical
    ffmmff94
    ffuff
    formalcharge
    format
    formula
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
    smartsparse
    smilesmatch
    unitcell
   )
if(LIBXML2_FOUND)
  set(oldtests ${oldtests} cmlreadfile)
endif(LIBXML2_FOUND)

if(WITH_INCHI)
  add_executable(test_inchiwrite inchiwrite.cpp ../src/formats/getinchi.cpp)
  target_link_libraries(test_inchiwrite ${BABEL_LIBRARY})
  # files in test/inchi -- both .sdf and .txt
  set(inchitests Samples SamplesTechMan)
  set(inchidata ${CMAKE_SOURCE_DIR}/test/inchi)
  foreach(test ${inchitests})
    add_test(inchi${test}_Test 
             ${TEST_PATH}/test_inchiwrite ${inchidata}/${test}.sdf ${inchidata}/${test}.txt)
    set_tests_properties(inchi${test}_Test PROPERTIES
                         FAIL_REGULAR_EXPRESSION "Not ok"
                         ENVIRONMENT "BABEL_LIBDIR=${CMAKE_BINARY_DIR}/lib"
                         ENVIRONMENT "BABEL_DATADIR=${CMAKE_SOURCE_DIR}/data"
	 	        )
  endforeach(test ${inchitests})
endif(WITH_INCHI)

foreach(test ${old_tests})
  add_executable(test_${test} ${test}.cpp)
  target_link_libraries(test_${test} ${BABEL_LIBRARY})
  add_test(${test}_Test ${TEST_PATH}/test_${test})
  set_tests_properties(${test}_Test PROPERTIES
    FAIL_REGULAR_EXPRESSION "not ok"
    ENVIRONMENT "BABEL_LIBDIR=${CMAKE_BINARY_DIR}/lib"
    ENVIRONMENT "BABEL_DATADIR=${CMAKE_SOURCE_DIR}/data"
 )
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

###########################
# Tests wrapped in Python #
###########################

# The following tests should work okay in MinGW, it's just that the calling
# script needs some work. On Cygwin, there's some weird DLL problems when
# calling from Python.
if(NOT MINGW AND NOT CYGWIN)
  include(UsePythonTest)
  if(PYTHON_EXECUTABLE)
    set(pytests
        babel sym fastsearch unique kekule)
    foreach(pytest ${pytests})
      ADD_PYTHON_TEST(pytest_${pytest} test${pytest}.py)
      set_tests_properties(pytest_${pytest} PROPERTIES
           FAIL_REGULAR_EXPRESSION "ERROR;FAIL;Test failed"
           ENVIRONMENT "BABEL_LIBDIR=${CMAKE_BINARY_DIR}/lib"
           ENVIRONMENT "BABEL_DATADIR=${CMAKE_SOURCE_DIR}/data"
      )
    endforeach(pytest ${pytests})
  endif(PYTHON_EXECUTABLE)
endif(NOT MINGW AND NOT CYGWIN)

###############################
# Tests using Python bindings #
###############################

if (PYTHON_BINDINGS)
  set(pybindtests
      bindings _pybel example)
  foreach(pybindtest ${pybindtests})
    ADD_PYTHON_TEST(pybindtest_${pybindtest} test${pybindtest}.py)
    set_tests_properties(pybindtest_${pybindtest} PROPERTIES
         FAIL_REGULAR_EXPRESSION "ERROR;FAIL;Test failed"
         ENVIRONMENT "BABEL_LIBDIR=${CMAKE_BINARY_DIR}/lib"
         ENVIRONMENT "BABEL_DATADIR=${CMAKE_SOURCE_DIR}/data"
    )
  endforeach(pybindtest ${pybindtests})
endif (PYTHON_BINDINGS)

