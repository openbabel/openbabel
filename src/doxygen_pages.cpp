// This file contains doxygen documentation only
//
// pages:
//  cmake_project
//  generic_data

namespace OpenBabel {

/**
 * @page cmake_project Creating your own projects using CMake
 *
 * CMake provides modules to find a large number of modules to find
 * common dependencies. OpenBabel provides configuration files for CMake
 * so it can automatically found using the find_package command. Calling
 * find_package will try to locate the configuration files of OpenBabel,
 * setting 3 variables:
 *
 * @li OpenBabel3_FOUND
 * @li OpenBabel3_INCLUDE_DIRS
 * @li OpenBabel3_LIBRARIES
 *
 * The find_package command allows you to specify the package is required
 * and cmake will handle this further. If openbabel is optional, the
 * first variable can be used in your cmake logic to optionally build the
 * additional code. Since find_package only sets variables, you still
 * need to call include_directories with OpenBabel3_INCLUDE_DIRS in the
 * argument list. The OpenBabel3_LIBRARIES variable can be used directly
 * in your target_link_libraries command.
 *
 * Below is a minimal but working example of a project. For simplicity,
 * only one CMakeLists.txt file and one source file is used. This can be
 * used as a template to get started. The <a href="http://www.cmake.org/cmake/help/documentation.html">
 * cmake documentation</a> can be consultated as your project becomes more
 * complex.
 *
 * @b CMakeLists.txt
 * @code
# this line is required for cmake backwards compatibility
cmake_minimum_required(VERSION 3.1)

# name of your project
project(myproject)

# find and setup openbabel
find_package(OpenBabel3 REQUIRED)
include_directories(${OpenBabel3_INCLUDE_DIRS})

# create a list of source files (easier to maintain)
set(sources main.cpp)

# the executable
add_executable(myexe ${sources})
target_link_libraries(myexe ${OpenBabel3_LIBRARIES})
install(TARGETS myexe DESTINATION bin)
   @endcode
   @b main.cpp
   @code
#include <openbabel/mol.h>

using namespace OpenBabel;

int main()
{
  OBMol mol;

  // ... see examples for code ...

 return 0;
}
   @endcode
 *
 */


 /**
  * @page generic_data Working with Generic Data
  * @since 2.3
  *
  * @section generic_data_intro Introduction
  * Generic data is a concept used in OpenBabel to store additional information in objects. The objects are
  * usually molecules, atoms or bonds (OBMol, OBAtom and OBBond). The data can literally be anything since
  * OBPairTemplate allows any datatype (classes have to be copyable though) to be stored. For example, a file
  * format contains some strings or numbers (e.g. QM energy, biological activity, chemical supplier & price,
  * ...) for each molecule and these can be stored in the OBMol object. When the file format is used to read
  * a file, the program can access and use this data. A concrete example is the PDB file format which
  * specifies a large number of protein specific data types. All data which cannot be stored using the API is
  * stored as strings in the OBMol object. The program (e.g. a 3D molecular viewer) can retrieve the data
  * (e.g. secondary structure) and use it. It would not be possible to add API methods for all this.
  *
  * @section generic_data_design Design
  * There are two abstract classes defining the interfaces. The OBGenericData interface makes it possible
  * to work with derived classes without knowing anything about the data itself. It contains methods
  * (@ref OBGenericData::SetAttribute and @ref OBGenericData::GetAttribute) for associating the data with a name.
  * To use std::map<std::string, T> analogy, the attribute is the key for the data T. GetValue always returns
  * a std::string and derived classes should convert their data to a string when possible. Returning an empty
  * string is acceptable though.
  * The second OBBase class defines an interface to store/retrieve/remove OBGenericData objects by attribute,
  * type or source. To use std::map analogy again, classes derived from OBBase are the map.
  *
  * @section generic_data_str_num Storing strings and numbers
  * In many cases storing strings and numbers is all you need. Strings can be stored using the OBPairData
  * class. For numbers there is OBPairInteger and OBPairFloatingPoint. Although the interface is almost the
  * same for these classes multiple examples are given to make it easier to copy/paste.
  *
  * Storing and retrieving a string:
  * @code
  * // storing a string
  * OBPairData *supplier = new OBPairData;
  * supplier->SetAttribute("supplier"); // the name or key for the data
  * supplier->SetValue("some supplier name/id"); // reading from a file for example
  * mol.SetData(supplier);
  *
  * // retrieve the string by attribute
  * if (mol.HasData("supplier")) {
  *   OBPairData *supplier = dynamic_cast<OBPairData*>(mol.GetData("supplier"));
  *   cout << "supplier: " << supplier->GetValue() << endl;
  * }
  * @endcode
  * Storing and retrieving an integer:
  * @code
  * // storing an integer
  * OBPairInteger *data = new OBPairInteger;
  * data->SetAttribute("numAromRings"); // the name or key for the data
  * data->SetValue(numAromRings); // computed before
  * mol.SetData(data);
  *
  * // retrieve the integer by attribute
  * if (mol.HasData("numAromRings")) {
  *   OBPairInteger *data = dynamic_cast<OBPairInteger*>(mol.GetData("numAromRings"));
  *   cout << "number of aromatic rings: " << data->GetGenericValue() << endl;
  * }
  * @endcode
  *
  * There is a small difference between strings and numbers. The main reason is that
  * GetValue always returns a string. OBPairInteger and OBPairFloatingPoint
  * are actually typedefs for OBPairTemplate which defines the appropriate GetGenericValue
  * method to return the numeric data type.
  *
  * Storing and retrieving a floating point value:
  * @code
  * // storing an integer
  * OBPairFloatingPoint *data = new OBPairFloatingPoint;
  * data->SetAttribute("activity"); // the name or key for the data
  * data->SetValue(8.3); // computed before
  * mol.SetData(data);
  *
  * // retrieve the integer by attribute
  * if (mol.HasData("activity")) {
  *   OBPairFloatingPoint *data = dynamic_cast<OBPairFloatingPoint*>(mol.GetData("activity"));
  *   cout << "biological activity: " << data->GetGenericValue() << endl;
  * }
  * @endcode
  *
  * @section generic_data_template Truly generic data using OBPairTemplate
  * Although there are a number of classes for specific data types, using OBPairTemplate
  * the same can be accomplished with less code. The second example illustrates this but
  * a simpler example is given first.
  *
  * Storing a list of suppliers in an OBMol object:
  * @code
  * typedef OBPairTemplate< std::vector<std::string> > SupplierData;
  * // storing the supplier list
  * SupplierData *data = new SupplierData;
  * data->SetAttribute("suppliers");
  * data->SetValue(suppliers);
  * mol.SetData(data);
  *
  * // retrieve the supplier list
  * if (mol.HasData("suppliers")) {
  *   SupplierData *data = dynamic_cast<SupplierData*>(mol.GetData("suppliers"));
  *   std::vector<std::string> &suppliers = data->GetGenericData();
  *   for (unsigned int i = 0; i < suppliers.size(); ++i)
  *     cout << suppliers[i] << endl;
  * }
  * @endcode
  *
  * Storing complex data in an OBMol object:
  * @code
  * // data representation struct
  * struct MyDataRepr {
  *   double value, error;
  *   string unit;
  * };
  * typedef OBPairTemplate< MyDataRepr > MyData;
  *
  * // storing the supplier list
  * MyData *data = new MyData;
  * data->SetAttribute("mydata");
  * MyDataRepr repr;
  * repr.value = 5.3;
  * repr.error = 0.3;
  * repr.unit = "kJ/mol";
  * data->SetValue(repr);
  * mol.SetData(data);
  *
  * // retrieve the supplier list
  * if (mol.HasData("mydata")) {
  *   MyData *data = dynamic_cast<MyData*>(mol.GetData("mydata"));
  *   MyDataRepr &repr = data->GetGenericData();
  *   cout << repr.value << " +/- " << repr.error << " " << repr.unit << endl;
  * }
  * @endcode
  *
  * @section generic_data_specific Specific data types
  * A number of specific OBGenericData subclasses are provided for frequently used
  * data types: AliasData, OBAngleData, OBAtomClassData, OBChiralData, OBCommentData,
  * OBConformerData, OBDOSData, OBElectronicTransitionData, OBExternalBondData,
  * OBGridData, OBMatrixData, OBNasaThermoData, OBOrbitalEnergyData, OBPairData,
  * OBRateData, OBRingData, OBRotamerList, OBRotationData, OBSerialNums, OBSetData,
  * OBStereoBase, OBSymmetryData, OBTorsionData, OBUnitCell, OBVectorData,
  * OBVibrationData, OBVirtualBond. Consult the documentation for these classes for
  * more information.
  *
  * @section generic_data_formats Generic data & file formats
  * Various file formats read and write generic data. This section contains an overview
  * of the data used by file formats. When adding or extending a file format it is highly
  * recommended to update this section.
  *
  * @subsection generic_data_specific_by_data Read data ordered by data type
  * This section only contains information on data types used in a similar way by at
  * least two file formats.
  *
  * OBUnitCell:
  *   @li @b cacaoformat: *.caccrt
  *   @li @b carformat: *.car *.arc
  *
  * @subsection generic_data_specific_by_format Read data ordered by format
  *
  * @b adfformat: *.adfout
  *   @li Attribute: "Dipole Moment", Type: OBVectorData, Value: dipile moment vector
  *   @li Attribute: "PartialCharges", Type: OBPairData, Value: "Mulliken"
  *   @li Attribute: "GridData", Type: OBGridData, Value: ??, Multiple
  *
  * @b cacaoformat: *.caccrt
  *   @li Attribute: "", Type: OBUnitCell, Value: the unit cell data
  *
  * @b carformat: *.car *.arc
  *   @li Attribute: "", Type: OBUnitCell, Value: the unit cell data
  *
  * @b chemkinformat: *.ck
  *   @li Attribute: "Rate data", Type: OBRateData, Value: the reaction rate data
  *
  *
  *
  *
  *
  *
  */


}

/// @file doxygen_pages.cpp
/// @brief Additional doxygen ocumentation.
