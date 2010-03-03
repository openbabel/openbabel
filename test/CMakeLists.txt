# define TESTDATADIR for tests that need input files
add_definitions(-DTESTDATADIR="\\"${CMAKE_SOURCE_DIR}/test/files/\\"")

# define FORMATDIR for location of format plugin binaries
add_definitions(-DFORMATDIR="\\"${openbabel_BINARY_DIR}/lib/\\"")

###########################################################
#  new tests using obtest.h
############################################################

set(tests
    uniqueid
    smiles
    stereo
    tetrahedral
    tetranonplanar
    tetraplanar
    squareplanar
    cistrans
    stereoperception
    graphsym
   )

set(TEST_PATH ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

foreach(test ${tests})
  add_executable(test_${test} ${test}test.cpp)
  target_link_libraries(test_${test} ${BABEL_LIBRARY})
  if(MSVC)
    add_test(${test}_Test ${CMAKE_BINARY_DIR}/test_${test})
  else(MSVC)
    add_test(${test}_Test ${TEST_PATH}/test_${test})
  endif(MSVC)
  set_tests_properties(${test}_Test PROPERTIES
    FAIL_REGULAR_EXPRESSION "ERROR;FAIL;Test failed")
  if(MSVC)
    set_target_properties(test_${test} PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
  endif(MSVC)
endforeach(test ${tests})

 
############################################################
#  old-style tests using "not ok"
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

# TODO: forgot that inchiwrite requires arguments -- it compiles now, but
#   doesn't run a test
if(WITH_INCHI)
  add_executable(test_inchiwrite inchiwrite.cpp ../src/formats/getinchi.cpp)
  target_link_libraries(test_inchiwrite ${BABEL_LIBRARY})
  if(MSVC)
    set_target_properties(test_inchiwrite PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
  endif(MSVC)
  # files in test/inchi -- both .sdf and .txt
  set(inchitests Samples SamplesTechMan)
  set(inchidata ${CMAKE_SOURCE_DIR}/test/inchi)
  foreach(test ${inchitests})
	  if(MSVC)
	    add_test(NAME inchi${test}_Test 
							COMMAND ${CMAKE_BINARY_DIR}/test_inchiwrite ${inchidata}/${test}.sdf ${inchidata}/${test}.txt)
	  else(MSVC)
	    add_test(NAME inchi${test}_Test 
							COMMAND ${TEST_PATH}/test_inchiwrite ${inchidata}/${test}.sdf ${inchidata}/${test}.txt)
    endif(MSVC)
	  set_tests_properties(inchi${test}_Test PROPERTIES
		  FAIL_REGULAR_EXPRESSION "Not ok"
			ENVIRONMENT "BABEL_DATADIR=${CMAKE_SOURCE_DIR}/data"
	 	)
  endforeach(test ${inchitests})
endif(WITH_INCHI)

foreach(test ${old_tests})
  add_executable(test_${test} ${test}.cpp)
  target_link_libraries(test_${test} ${BABEL_LIBRARY})
  if(MSVC)
    add_test(${test}_Test ${CMAKE_BINARY_DIR}/test_${test})
  else(MSVC)
    add_test(${test}_Test ${TEST_PATH}/test_${test})
  endif(MSVC)
  set_tests_properties(${test}_Test PROPERTIES
    FAIL_REGULAR_EXPRESSION "not ok"
    ENVIRONMENT "BABEL_DATADIR=${CMAKE_SOURCE_DIR}/data"
 )
  if(MSVC)
    set_target_properties(test_${test} PROPERTIES
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

###########################
# Tests wrapped in Python #
###########################

include(UsePythonTest)
ADD_PYTHON_TEST(pytest_babel testbabel.py)
ADD_PYTHON_TEST(pytest_sym testsym.py)
ADD_PYTHON_TEST(pytest_fastsearch testfastsearch.py)
ADD_PYTHON_TEST(pytest_unique testunique.py)
ADD_PYTHON_TEST(pytest_kekule testkekule.py)
foreach(pytest babel sym fastsearch unique)

  set_tests_properties(pytest_${pytest} PROPERTIES
       FAIL_REGULAR_EXPRESSION "ERROR;FAIL;Test failed"
       ENVIRONMENT "BABEL_LIBDIR=${CMAKE_BINARY_DIR}/lib"
  )
endforeach(pytest babel sym fastsearch unique)

