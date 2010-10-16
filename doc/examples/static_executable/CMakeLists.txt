#
# This script can be used to create static executables linking to the static
# OpenBabel2 library.
#
# This script requires OpenBabel to be build and installed. For example:
#
#  cd openbabel-2.3
#  mkdir build
#  cd build
#  cmake -DBUILD_SHARED=OFF -DCMAKE_INSTALL_PREFIX=/home/me/some/path ..
#  make
#  make install
#
# To compile your static executable:
#
#  cd myproject
#  mkdir build
#  cd build
#  cmake -DOpenBabel2_DIR=/home/me/some/path/lib/openbabel ..
#  make
#
# All plugins are inside the static libopenbabel.a but the symbols for the
# plugin classes have to be undefined. Plugins can be disabled by removing
# the class names from the format_classes, descriptor_classes, ... lists below.
#
# Compilers:
# - GNU GCC: tested on linux
# - Intel ICC: not tested, should be the same as gcc
# - MSVC: not tested, uses different name mangling
# - MinGW: not tested, should be the same as gcc
#

# This line is required for cmake backwards compatibility.
cmake_minimum_required(VERSION 2.6)

# Name of your project
project(myproject)
# Create a list of source files (easier to maintain)
set(sources myexe.cpp)
# Set the name for the executable
set(executable_target myexe)

#
# The plugin class names to include in the executable.
#
set(format_classes
    ACRFormat
    ADFOutputFormat
    ADFInputFormat
    OBT41Format
    AlchemyFormat
    AmberPrepFormat
    OBAPIInterface
    BallStickFormat
    BGFFormat
    BoxFormat
    CacaoFormat
    CacheFormat
    CARFormat
    CCCFormat
    CHEM3D1Format
    CHEM3D2Format
    ChemDrawBinaryFormat
    ChemDrawFormat
    ChemKinFormat
    CHTFormat
    CIFFormat
    CopyFormat
    CRK2DFormat
    CRK3DFormat
    CSRFormat
    CSSRFormat
    DlpolyConfigFormat
    DlpolyHISTORYFormat
    DMolFormat
    XXXFormat
    FASTAFormat
    FastSearchFormat
    FCHKFormat
    FEATFormat
    FenskeZmatFormat
    FHIainsFormat
    FingerprintFormat
    FreeFormFractionalFormat
    GAMESSOutputFormat
    GAMESSInputFormat
    GAMESSUKFormat
    GAMESSUKInputFormat
    GAMESSUKOutputFormat
    OBGaussianCubeFormat
    GaussianOutputFormat
    GaussianInputFormat
    GaussianZMatrixInputFormat
    GenBankFormat
    GhemicalFormat
    GROMOS96Format
    GULPFormat
    HINFormat
    JaguarOutputFormat
    JaguarInputFormat
    MCDLFormat
    MDLFormat
    MOLFormat
    SDFFormat
    mmCIFFormat
    MacroModFormat
    MNAFormat
    MOL2Format
    OBMoldenFormat
    MolproOutputFormat
    MolproInputFormat
    MolReportFormat
    MOPACFormat
    MOPACCARTFormat
    MOPACINTFormat
    MPDFormat
    MPQCFormat
    MPQCInputFormat
    MSIFormat
    OBMSMSFormat
    NulFormat
    NWChemOutputFormat
    NWChemInputFormat
    OBOpenDXCubeFormat
    OutputFormat
    PCModelFormat
    PDBFormat
    PDBQTFormat
    PNGFormat
    PovrayFormat
    PQRFormat
    PQSFormat
    PWscfFormat
    QChemOutputFormat
    QChemInputFormat
    ReportFormat
    SmiReactFormat
    RXNFormat
    ShelXFormat
    SMIFormat
    CANSMIFormat
    FIXFormat
    SVGFormat
    TextFormat
    ThermoFormat
    TinkerFormat
    TitleFormat
    TurbomoleFormat
    UniChemFormat
    VASPFormat
    ViewMolFormat
    XEDFormat
    XTCFormat
    XYZFormat
    YOBFormat
    ZINDOFormat
   )

set(descriptor_classes
    CanSmiles
    CompoundFilter
    MWFilter
    SmartsFilter
    TitleFilter
    FormulaDescriptor
    FPCount
    InChIFilter
    SmartsDescriptor
   )

set(fingerprint_classes
    fingerprint2
    PatternFP
   )

set(forcefield_classes
    OBForceFieldGaff
    OBForceFieldGhemical
    OBForceFieldMMFF94
    OBForceFieldUFF
   )

set(op_classes
    OpAddInIndex
    OpAddPolarH
    OpCanonical
    OpConformer
    OpFillUC
    OpEnergy
    OpMinimize
    OpGen2D
    OpGen3D
    Op2D
    OpNewS
    OpPartialCharge
    OpReadConformers
    OpSort
    OpUnique
    OpExtraOut
   )

set(charge_classes
    GasteigerCharges
    MMFF94Charges
    QEqCharges
    QTPIECharges
   )

################################################################################

#
# Set compile flags for various compilers.
#
if(MSVC)
  # Set cl flags for static compiling
  set(CMAKE_CXX_FLAGS_DEBUG "/MTd")
  set(CMAKE_EXE_LINKER_FLAGS_DEBUG "/INCREMENTAL:NO /NODEFAULTLIB:MSVCRT")
  set(CMAKE_CXX_FLAGS_RELEASE "/MT /O2 /Ob2 /D NDEBUG")
  set(CMAKE_EXE_LINKER_FLAGS_RELEASE "/INCREMENTAL:NO /NODEFAULTLIB:MSVCRT")
  # Note: static libraries are specified when running cmake

  #
  # Perform name mangling for plugin symbols
  #
  # FIXME...
else()
  # Use -static flag to create static executable
  set(CMAKE_CXX_FLAGS "-static ${CMAKE_CXX_FLAGS}")
  # Make sure we find static libraries
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")

  #
  # Perform name mangling for plugin symbols
  #
  foreach(plugin_classes ${format_classes} ${descriptor_classes} ${fingerprint_classes} ${forcefield_classes} ${op_classes} ${charge_classes})
    foreach(plugin_class ${plugin_classes})
      # Perform simple name mangling for symbol
      string(LENGTH ${plugin_class} length)
      set(plugin_symbols_link_flags "${plugin_symbols_link_flags} -u _ZTVN9OpenBabel${length}${plugin_class}E")
    endforeach()
  endforeach()
  set(CMAKE_CXX_FLAGS "${plugin_symbols_link_flags} ${CMAKE_CXX_FLAGS}")
endif()

# Set the path containing OpenBabel2Config.cmake, needed for find_package below.
find_path(OpenBabel2_DIR OpenBabel2Config.cmake PATHS
          ${OpenBabel2_DIR}
          "/usr/lib/openbabel"
         "/usr/local/lib/openbabel")

#
# Find and setup OpenBabel2.
#
find_package(OpenBabel2 REQUIRED)
include_directories(${OpenBabel2_INCLUDE_DIRS})

# The executable
add_executable(${executable_target} ${sources})
# Link against imported openbabel target
target_link_libraries(${executable_target} openbabel)
# Prevent -Wl,-Bdynamic from being added to the end of the link line.
set_target_properties(${executable_target} PROPERTIES
                      LINK_SEARCH_END_STATIC TRUE)
install(TARGETS ${executable_target} DESTINATION bin)
