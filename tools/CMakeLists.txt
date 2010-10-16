if(MINIMAL_BUILD)
  set(tools
        babel
  )
else(MINIMAL_BUILD)
  set(tools
        babel
        obabel
        obchiral
        obconformer
        obenergy
        obfit
        obgen
        obminimize
        obprobe
        obprop
        obrotamer
        obrotate
  )
endif(MINIMAL_BUILD)

INCLUDE(CheckFunctionExists)
SET(GETOPT getopt)
CHECK_FUNCTION_EXISTS(${GETOPT} GETOPT_FOUND)

if(BUILD_SHARED)
  foreach(tool ${tools})
    add_executable(${tool} ${tool}.cpp)
    target_link_libraries(${tool} ${BABEL_LIBRARY})
    install(TARGETS ${tool}
                    RUNTIME DESTINATION bin
                    LIBRARY DESTINATION lib
                    ARCHIVE DESTINATION lib
    )
  endforeach(tool)

  if(NOT MINIMAL_BUILD)
    # obspectrophore -- requires getopt
    if (NOT GETOPT_FOUND)
      add_executable(obspectrophore obspectrophore.cpp getopt.c)
    else(NOT GETOPT_FOUND)
      add_executable(obspectrophore obspectrophore.cpp)
    endif(NOT GETOPT_FOUND)
    target_link_libraries(obspectrophore ${BABEL_LIBRARY})
    install(TARGETS obspectrophore
                    RUNTIME DESTINATION bin
                    LIBRARY DESTINATION lib
                    ARCHIVE DESTINATION lib
    )

    # obgrep -- requires getopt
    if (NOT GETOPT_FOUND)
      add_executable(obgrep obgrep.cpp getopt.c)
    else(NOT GETOPT_FOUND)
      add_executable(obgrep obgrep.cpp)
    endif(NOT GETOPT_FOUND)
      target_link_libraries(obgrep ${BABEL_LIBRARY})
      install(TARGETS obgrep
                      RUNTIME DESTINATION bin
                      LIBRARY DESTINATION lib
                      ARCHIVE DESTINATION lib
      )
  endif(NOT MINIMAL_BUILD)

else(BUILD_SHARED)
###############################################################################
#
# Build STATIC babel with all plugins taken from the static libopenbabel.a
#
# To ensure the plugins are not removed from the final executable, the relevant
# plugin class (the classes inheriting OBPlugin) symbols have to be undefined
# while linking to ensure the symbols are taken from the static libopenabbel.a.
#
# GNU GCC & Intel ICC
#
# name mangling: _ZTVN9OpenBabel${length}${plugin_class}E
# undefine symbol compiler option: -u symbol
#
# MSVC:
#
# name mangling: FIXME (not the same as gcc,icc)
# undefine symbol compiler option: /U symbol
#
###############################################################################

  # These are the class names for the formats.
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


  if(MSVC)
    #
    # Microsoft Visual C++ compiler
    #
  else(MSVC)
    #
    # GNU GCC or Intel ICC compiler.
    #
    # NOTE: not tested for icc but the options and name mangling should be the same.
    foreach(plugin_classes ${format_classes} ${descriptor_classes} ${fingerprint_classes} ${forcefield_classes} ${op_classes} ${charge_classes})
      foreach(plugin_class ${plugin_classes})
        # Perform simple name mangling for symbol
        string(LENGTH ${plugin_class} length)
        set(plugin_symbols_link_flags "${plugin_symbols_link_flags} -u _ZTVN9OpenBabel${length}${plugin_class}E")
      endforeach(plugin_class ${plugin_classes})
    endforeach(plugin_classes ${format_classes} ${descriptor_classes} ${fingerprint_classes} ${forcefield_classes} ${op_classes} ${charge_classes})
    set(CMAKE_CXX_FLAGS "${plugin_symbols_link_flags} ${CMAKE_CXX_FLAGS}")

    if (BUILD_MIXED)
      #
      # Relocatable binaries on linux using -static-libgcc. (GCC only)
      #
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc")
      # make sure the linker will find the static libstdc++
      file(REMOVE "${CMAKE_BINARY_DIR}/libstdc++.a")
      execute_process(COMMAND "g++" "-print-file-name=libstdc++.a"
                      OUTPUT_VARIABLE libstdcpp
                      OUTPUT_STRIP_TRAILING_WHITESPACE)
      execute_process(COMMAND "ln" "-s" "${libstdcpp}"
                      WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
      link_directories(${CMAKE_BINARY_DIR})
    else(BUILD_MIXED)
      #
      # Default static build using GNU-like compiler.
      #
      set(CMAKE_CXX_FLAGS "-static ${CMAKE_CXX_FLAGS}")
    endif(BUILD_MIXED)
  endif(MSVC)

  # build a static babel only
  add_executable(babel babel.cpp)
  target_link_libraries(babel openbabel)
  # Prevent -Wl,-Bdynamic from being added to the end of the link line.
  set_target_properties(babel PROPERTIES
                        LINK_SEARCH_END_STATIC TRUE)
  install(TARGETS babel
          RUNTIME DESTINATION bin
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)

endif(BUILD_SHARED)
