# define TESTDATADIR for tests that need input files
add_definitions(-DTESTDATADIR="${CMAKE_SOURCE_DIR}/test/files/")

# define FORMATDIR for location of format plugin binaries
if(BINDINGS_ONLY)
  set(FORMATDIR "${OB_MODULE_PATH}/")
else()
  set(FORMATDIR "${openbabel_BINARY_DIR}/lib${LIB_SUFFIX}/")
endif()

add_definitions(-DFORMATDIR="${FORMATDIR}/")

# TRUE for build types with optimization; slow tests are skipped otherwise
if(CMAKE_BUILD_TYPE MATCHES "^(Release|RelWithDebInfo|MinSizeRel)$")
  set(OB_OPTIMIZED_BUILD TRUE)
else()
  set(OB_OPTIMIZED_BUILD FALSE)
endif()

# ##########################################################
# new tests using obtest.h
# ###########################################################

# ############### Add new tests here
set(cpptests
  alias automorphism builder canonconsistent canonfragment canonstable carspacegroup cifspacegroup
  cistrans conversion fuzzregress graphsym gzip addh
  implicitH ketformat lssr isomorphism multicml periodic regressions rotor shuffle smiles spectrophore
  squareplanar stereo stereoperception tautomer tetrahedral
  tetranonplanar tetraplanar uniqueid
)
set(alias_parts 1)
set(automorphism_parts 1 2 3 4 5 6 7 8 9 10)
set(builder_parts 1 2 3 4 5 6 7 8 9 10)
set(canonconsistent_parts 1 2 3)
set(canonfragment_parts 1)
set(canonstable_parts 1)
set(carspacegroup_parts 1 2 3 4)
set(cifspacegroup_parts 1 2 3 4 5 6 7 8 9 10 11 12 13)
set(cistrans_parts 1 2 3 4 5 6 7 8 9)
set(conversion_parts 1)
set(fuzzregress_parts 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27)
set(graphsym_parts 1 2 3 4 5)
set(gzip_parts 1)
set(addh_parts 1)
set(implicitH_parts 1)
set(ketformat_parts 1 2 3 4 5 6 7 8 9 10 11 12 13)
set(lssr_parts 1 2 3 4 5)
set(isomorphism_parts 1 2 3 4 5 6 7 8 9)
set(multicml_parts 1)
set(periodic_parts 1 2 3 4)
set(regressions_parts 1 2 221 222 223 224 225 226 227 228 229 240 241 242 1794 2111 2428 2646 2677)
set(rotor_parts 1 2 3 4)
set(shuffle_parts 1 2 3 4 5)
set(smiles_parts 1 2 3)
set(spectrophore_parts 1 2 3 4 5)
set(squareplanar_parts 1 2 3 4 5)
set(stereo_parts 1 2 3 4 5 6)
set(stereoperception_parts 1 2 3 4)
set(tautomer_parts 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29)
set(tetrahedral_parts 1 2 3 4 5)
set(tetranonplanar_parts 1)
set(tetraplanar_parts 1)
set(uniqueid_parts 1 2)

if(TARGET Eigen3::Eigen)
  set(cpptests
    align distgeom gen3d ${cpptests})
  set(align_parts 1 2 3 4 5)
  # Parts 1-5 are fast; parts 6-7 are large molecules (~10-60s)
  # Parts 8-19 are the pr2317 hard cases (new distgeom fallback tests).
  # Excluded (fused/bridged ring systems that exhaust L-BFGS on every trial):
  #   13 (quinine/quinidine), 15 (polycyclic alkaloids), 16 (aminoglycoside)
  # Excluded (too large for 30 s wall-clock limit):
  #   19 (disaccharide, 38+ atoms)
  # All excluded cases are covered by gen3dtest via the builder path.
  if(NOT OB_OPTIMIZED_BUILD)
    # Skip slow tests in unoptimized builds (Debug, sanitizer builds, etc.)
    # Part 6/20 (stereoisomer coverage) uses small molecules and is fast.
    set(distgeom_parts 1 2 3 4 5 20)
    set(gen3d_parts 1 2 3 4 6)
  else()
    # In optimized builds, include comprehensive test coverage
    set(distgeom_parts 1 2 3 4 5 6 7 8 9 10 11 12 14 17 18 20)
    set(gen3d_parts 1 2 3 4 5 6)
  endif()
endif()

if(WITH_MAEPARSER)
  set(cpptests ${cpptests}
    maereader)
  set(maereader_parts 1 2)
endif()

set(origtests
  aromatest atom bond cansmi charge_mmff94 charge_gasteiger conversion
  datatest ffgaff ffghemical ffmmff94 ffuff formalcharge format formula
  internalcoord invalidsmarts invalidsmiles iterators logp_psa math mol
  pdbreadfile phmodel residue ringtest smartstest smartsparse smilesmatch
  unitcell
)
set(atom_parts 1 2 3 4)
set(ffmmff94_parts 1 2 3 4 5 6)
set(math_parts 1 2 3 4)
set(pdbreadfile_parts 1 2 3 4)

if(BUILD_SHARED)
  if(LIBXML2_FOUND)
    set(origtests ${origtests} cmlreadfile)
  endif(LIBXML2_FOUND)
else()
  if(WITH_STATIC_LIBXML)
    set(origtests ${origtests} cmlreadfile)
  endif()
endif()

# #################################
foreach(cpptest ${cpptests})
  set(cpptestsrc ${cpptestsrc} ${cpptest}test.cpp)
endforeach()

foreach(origtest ${origtests})
  set(cpptestsrc ${cpptestsrc} ${origtest}.cpp)
endforeach()

foreach(origtest ${origtests})
  if(NOT DEFINED "${origtest}_parts")
    set(${origtest}_parts "1")
  endif()
endforeach()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
create_test_sourcelist(srclist test_runner.cpp
  ${cpptestsrc})

set(TEST_PATH ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

if(NOT BUILD_SHARED)
  set(libs openbabel)

  if(WITH_STATIC_LIBXML)
    set(libs ${libs} ${LIBXML2_LIBRARIES})
  endif()

  if(WITH_STATIC_INCHI)
    set(libs ${libs} ${INCHI_LIBRARY})
  endif()

  if(NOT MSVC)
    if(BUILD_MIXED)
      if(CMAKE_COMPILER_IS_GNUCXX)
        if(CMAKE_SYSTEM_NAME MATCHES Linux)
          #
          # Relocatable binaries on linux using -static-libgcc. (GCC only)
          #
          set(CMAKE_CXX_FLAGS "-static-libgcc ${CMAKE_CXX_FLAGS}")

          # make sure the linker will find the static libstdc++
          file(REMOVE "${CMAKE_BINARY_DIR}/libstdc++.a")
          execute_process(COMMAND "g++" "-print-file-name=libstdc++.a"
            OUTPUT_VARIABLE libstdcpp
            OUTPUT_STRIP_TRAILING_WHITESPACE)
          execute_process(COMMAND "ln" "-s" "${libstdcpp}"
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

          link_directories(${CMAKE_BINARY_DIR})
        endif(CMAKE_SYSTEM_NAME MATCHES Linux)
      endif(CMAKE_COMPILER_IS_GNUCXX)
    else(BUILD_MIXED)
      if(NOT APPLE)
        set(CMAKE_CXX_FLAGS "-static ${CMAKE_CXX_FLAGS}")
      endif()
    endif(BUILD_MIXED)
  endif()
else()
  set(libs openbabel)
endif()

add_executable(test_runner ${srclist} obtest.cpp)
target_link_libraries(test_runner ${libs})
set_target_properties(test_runner PROPERTIES ENABLE_EXPORTS TRUE)

if(NOT BUILD_SHARED AND NOT BUILD_MIXED)
  set_target_properties(test_runner PROPERTIES LINK_SEARCH_END_STATIC TRUE)
endif()

foreach(cpptest ${cpptests})
  foreach(part ${${cpptest}_parts})
    add_test(test_${cpptest}_${part}
      ${TEST_PATH}/test_runner ${cpptest}test ${part})
    set_tests_properties(test_${cpptest}_${part} PROPERTIES
      FAIL_REGULAR_EXPRESSION "ERROR;FAIL;Test failed"
      ENVIRONMENT "BABEL_DATADIR=${CMAKE_SOURCE_DIR}/data;LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib${LIB_SUFFIX}:$ENV{LD_LIBRARY_PATH}")
  endforeach()
endforeach()

foreach(cpptest ${origtests})
  foreach(part ${${cpptest}_parts})
    add_test(test_${cpptest}_${part}
      ${TEST_PATH}/test_runner ${cpptest} ${part})
    set_tests_properties(test_${cpptest}_${part} PROPERTIES
      FAIL_REGULAR_EXPRESSION "not ok"
      ENVIRONMENT "BABEL_DATADIR=${CMAKE_SOURCE_DIR}/data;LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib${LIB_SUFFIX}:$ENV{LD_LIBRARY_PATH}")
  endforeach()
endforeach()

# ###########################################################
# old-style tests using "not ok"
# ###########################################################
if(NOT BUILD_SHARED)
  set(WITH_INCHI ${WITH_STATIC_INCHI})
endif()

if(WITH_INCHI)
  add_executable(test_inchiwrite inchiwrite.cpp ../src/formats/getinchi.cpp)
  target_link_libraries(test_inchiwrite ${libs})

  if(NOT BUILD_SHARED AND NOT BUILD_MIXED)
    set_target_properties(test_inchiwrite PROPERTIES LINK_SEARCH_END_STATIC TRUE)
  endif()

  # files in test/inchi -- both .sdf and .txt
  set(inchitests Samples.sdf SamplesTechMan.sdf Steffen_PubChem.smi)
  set(inchidata ${CMAKE_SOURCE_DIR}/test/inchi)

  foreach(test ${inchitests})
    add_test(inchi${test}_Test
      ${TEST_PATH}/test_inchiwrite ${inchidata}/${test} ${inchidata}/${test}.txt)
    set_tests_properties(inchi${test}_Test PROPERTIES
      FAIL_REGULAR_EXPRESSION "Not ok"
      ENVIRONMENT "BABEL_DATADIR=${CMAKE_SOURCE_DIR}/data;LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib${LIB_SUFFIX}:$ENV{LD_LIBRARY_PATH}"
    )
  endforeach(test ${inchitests})
endif(WITH_INCHI)

# ###########################################################
# roundtrip
# ###########################################################
add_executable(roundtrip roundtrip.cpp)
target_link_libraries(roundtrip ${libs})

if(NOT BUILD_SHARED AND NOT BUILD_MIXED)
  set_target_properties(roundtrip PROPERTIES LINK_SEARCH_END_STATIC TRUE)
endif()

install(TARGETS roundtrip
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib${LIB_SUFFIX}
  ARCHIVE DESTINATION lib${LIB_SUFFIX}
)

# ##########################
# Tests wrapped in Python #
# ##########################

# The following tests should work okay in MinGW, it's just that the calling
# script needs some work. On Cygwin, there's some weird DLL problems when
# calling from Python.
if(NOT MINGW AND NOT CYGWIN)
  include(UsePythonTest)

  if(Python_EXECUTABLE)
    # Exclude distgeom from unoptimized builds (Debug, sanitizers) — it times out
    if(NOT OB_OPTIMIZED_BUILD)
      set(pytests
        babel gauss sym smartssym fastsearch unique kekule pdbformat RInChI)
    else()
      set(pytests
        babel gauss sym smartssym fastsearch distgeom unique kekule pdbformat RInChI)
    endif()

    foreach(pytest ${pytests})
      ADD_TEST(NAME pytest_${pytest} COMMAND ${Python_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/test${pytest}.py")
      set_tests_properties(pytest_${pytest} PROPERTIES
        FAIL_REGULAR_EXPRESSION "ERROR;FAIL;Test failed"
        ENVIRONMENT "PYTHONPATH=${CMAKE_SOURCE_DIR}/scripts/python:${CMAKE_BINARY_DIR}/lib${LIB_SUFFIX};BABEL_LIBDIR=${FORMATDIR};BABEL_DATADIR=${CMAKE_SOURCE_DIR}/data;LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib${LIB_SUFFIX}:$ENV{LD_LIBRARY_PATH}"
      )
    endforeach(pytest ${pytests})
  endif(Python_EXECUTABLE)
endif(NOT MINGW AND NOT CYGWIN)

# ##############################
# Tests using Python bindings #
# ##############################
if(PYTHON_BINDINGS)
  include(UsePythonTest)
  set(pybindtests
    bindings
    _pybel
    example
    obconv_writers
    cdjsonformat
    ketformat
    pcjsonformat
    roundtrip
  )

  foreach(pybindtest ${pybindtests})
    if(NOT MSVC)
      set(PYTHONPATH "${CMAKE_SOURCE_DIR}/scripts/python:${CMAKE_BINARY_DIR}/lib${LIB_SUFFIX}")
    else()
      set(PYTHONPATH "${CMAKE_SOURCE_DIR}/scripts/python\;${CMAKE_BINARY_DIR}/bin/Release")
    endif()

    ADD_TEST(NAME pybindtest_${pybindtest} COMMAND ${Python_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/test${pybindtest}.py")
    set_tests_properties(pybindtest_${pybindtest} PROPERTIES
      FAIL_REGULAR_EXPRESSION "ERROR;FAIL;Test failed"
      ENVIRONMENT "PYTHONPATH=${PYTHONPATH};BABEL_LIBDIR=${FORMATDIR};BABEL_DATADIR=${CMAKE_SOURCE_DIR}/data;LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib${LIB_SUFFIX}:$ENV{LD_LIBRARY_PATH}"
    )
  endforeach(pybindtest ${pybindtests})
endif(PYTHON_BINDINGS)

# ###########################################################
# Test 2D depiction using ImageMagick utilities
# ###########################################################

# Example: add_depict_test(${CMAKE_SOURCE_DIR}/test/files/aromatics.smi aromatics ${CMAKE_SOURCE_DIR}/test/files/alias.png)
function(add_depict_test testname obfile image_file ref_image)
  set(depict_dir ${CMAKE_BINARY_DIR}/depict_test)
  add_test(NAME test_${testname} COMMAND ${CMAKE_SOURCE_DIR}/test/test_depiction.sh
    ${CMAKE_BINARY_DIR}/bin/obabel ${ImageMagick_convert_EXECUTABLE}
    ${ImageMagick_compare_EXECUTABLE} ${obfile} ${depict_dir}
    ${image_file} ${ref_image})
  set_tests_properties(test_${testname} PROPERTIES
    PASS_REGULAR_EXPRESSION "0 @ 0,0")
endfunction()

option(TEST_DEPICTION "Test depiction using ImageMagick utilities" OFF)

if(TEST_DEPICTION AND UNIX)
  find_package(ImageMagick COMPONENTS convert compare)

  if(ImageMagick_convert_FOUND AND ImageMagick_compare_FOUND)
    set(depict_dir ${CMAKE_BINARY_DIR}/depict_test)
    file(MAKE_DIRECTORY ${depict_dir})
    add_depict_test(depict_aromatics ${CMAKE_SOURCE_DIR}/test/files/aromatics.smi
      aromatics ${CMAKE_SOURCE_DIR}/test/files/aromatics.png)
    add_depict_test(depict_alias ${CMAKE_SOURCE_DIR}/test/files/alias.mol
      alias ${CMAKE_SOURCE_DIR}/test/files/alias.png)
  endif()
endif()

# Fuzz harnesses. When LIB_FUZZING_ENGINE is set (OSS-Fuzz build), link
# against the engine. Otherwise link against a small standalone main so
# saved corpus files can be replayed as ordinary CTest regression cases.
if(DEFINED ENV{LIB_FUZZING_ENGINE})
  set(_fuzz_main_lib $ENV{LIB_FUZZING_ENGINE})
  set(_fuzz_extra_src "")
else()
  set(_fuzz_main_lib "")
  set(_fuzz_extra_src fuzz/standalone_main.cc)
endif()
  
add_executable(fuzz_obconversion_smiles fuzz/fuzz_obconversion.cpp ${_fuzz_extra_src})
target_compile_definitions(fuzz_obconversion_smiles PRIVATE -DFUZZ_INPUT_FORMAT="SMILES")
target_link_libraries(fuzz_obconversion_smiles ${libs} ${_fuzz_main_lib})

add_executable(fuzz_obconversion_sdf fuzz/fuzz_obconversion.cpp ${_fuzz_extra_src})
target_compile_definitions(fuzz_obconversion_sdf PRIVATE -DFUZZ_INPUT_FORMAT="SDF")
target_link_libraries(fuzz_obconversion_sdf ${libs} ${_fuzz_main_lib})

add_executable(fuzz_convert fuzz/fuzz_convert.cpp ${_fuzz_extra_src})
target_link_libraries(fuzz_convert ${libs} ${_fuzz_main_lib})
# fuzz_convert.cpp pulls in FuzzedDataProvider.h, which requires C++17
# (e.g. std::is_enum_v). The top-level project default is still C++11.
set_target_properties(fuzz_convert PROPERTIES CXX_STANDARD 17)

add_executable(fuzz_molecule fuzz/fuzz_molecule.cpp ${_fuzz_extra_src})
target_link_libraries(fuzz_molecule ${libs} ${_fuzz_main_lib})
set_target_properties(fuzz_molecule PROPERTIES CXX_STANDARD 17)

add_executable(fuzz_smart fuzz/fuzz_smart.cpp ${_fuzz_extra_src})
target_link_libraries(fuzz_smart ${libs} ${_fuzz_main_lib})
set_target_properties(fuzz_smart PROPERTIES CXX_STANDARD 17)

add_executable(fuzz_empty_write fuzz/fuzz_empty_write.cpp ${_fuzz_extra_src})
target_link_libraries(fuzz_empty_write ${libs} ${_fuzz_main_lib})
set_target_properties(fuzz_empty_write PROPERTIES CXX_STANDARD 17)

#add_executable(fuzz_depict fuzz/fuzz_depict.cpp ${_fuzz_extra_src})
#target_link_libraries(fuzz_depict ${libs} ${_fuzz_main_lib})
#set_target_properties(fuzz_depict PROPERTIES CXX_STANDARD 17)

# Replay each saved corpus file as a regression test. CMake re-evaluates
# the glob at configure time, so dropping a new file into the corpus
# directory and re-running cmake is enough to register a new test.
if(NOT DEFINED ENV{LIB_FUZZING_ENGINE})
  set(_fuzz_targets fuzz_convert fuzz_obconversion_sdf fuzz_obconversion_smiles fuzz_molecule fuzz_smart fuzz_empty_write) #fuzz_depict
  foreach(_target ${_fuzz_targets})
    file(GLOB _corpus_files CONFIGURE_DEPENDS
      "${CMAKE_CURRENT_SOURCE_DIR}/fuzz/corpus/${_target}/*")
    foreach(_corpus_file ${_corpus_files})
      get_filename_component(_corpus_name "${_corpus_file}" NAME)
      add_test(NAME ${_target}_${_corpus_name}
        COMMAND ${TEST_PATH}/${_target} "${_corpus_file}")
      set_tests_properties(${_target}_${_corpus_name} PROPERTIES
        ENVIRONMENT "BABEL_DATADIR=${CMAKE_SOURCE_DIR}/data;BABEL_LIBDIR=${FORMATDIR};LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib${LIB_SUFFIX}:$ENV{LD_LIBRARY_PATH}")
    endforeach()
  endforeach()
endif()
