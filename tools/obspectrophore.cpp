#include <iostream>

#ifdef WIN32
#include "getopt.h"
#else
#include <unistd.h>  // For the getopt() function
#endif

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/spectrophore.h>
#include <cstdlib>



void showParameters(const OpenBabel::OBSpectrophore&, const std::string&);
void showImplementationDetails(const char*);
void showError(const std::string&);
void showHelp(const char*);
bool isValidValue(const char, char*);
OpenBabel::OBSpectrophore::NormalizationOption stringToNormalizationOption(char*);
OpenBabel::OBSpectrophore::AccuracyOption stringToAccuracyOption(char*);
OpenBabel::OBSpectrophore::StereoOption stringToStereoOption(char*);



int main(int argc,char **argv)
{
	std::string ifile = "";
	OpenBabel::OBSpectrophore::AccuracyOption accuracy = OpenBabel::OBSpectrophore::AngStepSize20;
	OpenBabel::OBSpectrophore::StereoOption stereo = OpenBabel::OBSpectrophore::NoStereoSpecificProbes;
	OpenBabel::OBSpectrophore::NormalizationOption normalization = OpenBabel::OBSpectrophore::NoNormalization;
	double resolution = 3.0;
	int c;
	
	opterr = 0;
	std::string msg;
	
	while ((c = getopt(argc, argv, "ui:n:a:s:r:h")) != -1)
	{
		switch (c)
		{
			case 'u':
            showImplementationDetails(argv[0]);
            exit(1);
				break;
				
			case 'i':
				if (!isValidValue('i', optarg))
				{
					msg = "Option -i is followed by an invalid argument: ";
					msg += optarg;
					showError(msg);
					exit(1);
				}
				else
				{
					ifile = optarg;
				}
				break;
				
			case 'n':
				if (!isValidValue('n', optarg))
				{
					msg = "Option -n is followed by an invalid argument: ";
					msg += optarg;
					showError(msg);
					exit(1);
				}
				else normalization = stringToNormalizationOption(optarg);
				break;
				
			case 'a':
				if (!isValidValue('a', optarg))
				{
					msg = "Option -a is followed by an invalid argument: ";
					msg += optarg;
					showError(msg);
					exit(1);
				}
				else accuracy = stringToAccuracyOption(optarg);
				break;
				
			case 's':
				if (!isValidValue('s', optarg))
				{
					msg = "Option -s is followed by an invalid argument: ";
					msg += optarg;
					showError(msg);
					exit(1);
				}
				else stereo = stringToStereoOption(optarg);
				break;
				
			case 'r':
				if (!isValidValue('r', optarg))
				{
					msg = "Option -r is followed by an invalid argument: ";
					msg += optarg;
					showError(msg);
					exit(1);
				}
				else
				{
					resolution = atof(optarg);
					if (resolution <= 0)
					{
						msg = "Resolution -r should be larger than 0.0: ";
						msg += optarg;
						showError(msg);
						exit(1);
					}
				}
				break;
				
			case 'h':
				showHelp(argv[0]);
				exit(0);
				break;
				
			case '?':
				if ((optopt == 'i') || 
						(optopt == 'n') || 
						(optopt == 'a') || 
						(optopt == 's') || 
						(optopt == 'r'))
				{
					msg = "Option -";
					msg += optopt;
					msg += " requires an argument.";
					showError(msg);
					exit(1);
				}
				else
				{
					msg = "Unknown option -";
					msg += optopt;
					msg += ".";
					showError(msg);
					exit(1);
				}
				break;
				
			default:
				showError("Unknown option");
				exit(1);
				break;
		}
	}
	
	// The input file (-i) is the only required option
	if (ifile.empty())
	{
		msg = "Input file specification is required (option -i).";
		showError(msg);
		exit(1);
	}
	OpenBabel::OBConversion obconversion;
	OpenBabel::OBFormat *format = obconversion.FormatFromExt(ifile.c_str());
	if (!format)
	{
		msg = "Could not find file format for ";
		msg += ifile;
		showError(msg);
		exit(1);
	}
	obconversion.SetInFormat(format);
	std::ifstream ifs;
	ifs.open(ifile.c_str());
	obconversion.SetInStream(&ifs);
	
	// Start calculations
	OpenBabel::OBMol mol;
	OpenBabel::OBSpectrophore spec;
	spec.SetAccuracy(accuracy);
	spec.SetNormalization(normalization);
	spec.SetStereo(stereo);
	spec.SetResolution(resolution);
	showParameters(spec, ifile);
	unsigned int count(0);
	while (obconversion.Read(&mol))
	{
		std::vector<double> result = spec.GetSpectrophore(&mol);
		if (result.empty()) {
			std::cerr << "Error calculating Spectrophore from molecule number ";
			std::cerr << count;
			std::cerr << " (counting starts at 0)!";
			std::cerr << std::endl;
			exit(1);
		}
		else
		{
			std::cout << mol.GetTitle() << "\t";
			for (unsigned int i(0); i < result.size(); ++i)
			{
				std::cout << result[i] << "\t";
			}
			std::cout << std::endl;
		}
		mol.Clear();
		++count;
	}
	return 0;
}



bool
isValidValue(const char c, char* v)
{
	std::string o = v;
	for (unsigned int i = 0; i < o.length(); i++) o[i] = toupper(o[i]);
	switch (c)
	{
		case 'i':
			if (o.empty())
         {
            return false;
         }
         else
         {
            return true;
         }
		case 'n':
			if ((!o.compare("NO")) ||
					(!o.compare("ZEROMEAN")) ||
					(!o.compare("UNITSTD")) ||
					(!o.compare("ZEROMEANANDUNITSTD")))
         {
            return true;
         }
			else
         {
            return false;
         }
		case 'a':
			if ((!o.compare("1")) ||
					(!o.compare("2")) ||
					(!o.compare("5")) ||
					(!o.compare("10")) ||
					(!o.compare("15")) ||
					(!o.compare("20")) ||
					(!o.compare("30")) ||
					(!o.compare("36")) ||
					(!o.compare("45")) ||
					(!o.compare("60")))
         {
            return true;
         }
			else
         {
            return false;
         }
		case 's':
			if ((!o.compare("NO")) ||
					(!o.compare("UNIQUE")) ||
					(!o.compare("MIRROR")) ||
					(!o.compare("ALL")))
         {
            return true;
         }
         else
         {
            return false;
         }
		case 'r':
			// Check if real positive number
			if (o.empty() || (o[0] == '-'))
         {
            return false;
         }
			else
			{
				// Remove '+'as first character
				if (o[0] == '+') o.erase(0,1);
				
				// Remove maximal one '.'
				for (unsigned int i = 0; i < o.length(); i++) 
				{
					if (o[i] == '.')
					{
						o.erase(i,1);
						break;
					}
				}
				
				// Check if remaining characters are digits
				for (unsigned int i = 0; i < o.length(); i++)
            {
               if (!isdigit(o[i]))
               {
                  return false;
               }
            }
			}
	}
	return true;
}



OpenBabel::OBSpectrophore::NormalizationOption
stringToNormalizationOption(char* v)
{
	std::string o = v;
	for (unsigned int i = 0; i < o.length(); i++) o[i] = toupper(o[i]);
	if (!o.compare("NO"))                 return OpenBabel::OBSpectrophore::NoNormalization;
	if (!o.compare("ZEROMEAN"))           return OpenBabel::OBSpectrophore::NormalizationTowardsZeroMean;
	if (!o.compare("UNITSTD"))            return OpenBabel::OBSpectrophore::NormalizationTowardsUnitStd;
	if (!o.compare("ZEROMEANANDUNITSTD")) return OpenBabel::OBSpectrophore::NormalizationTowardsZeroMeanAndUnitStd;
	return OpenBabel::OBSpectrophore::NoNormalization;
}



OpenBabel::OBSpectrophore::AccuracyOption
stringToAccuracyOption(char* v)
{
	std::string o = v;
	if (!o.compare("1"))  return OpenBabel::OBSpectrophore::AngStepSize1;
	if (!o.compare("2"))  return OpenBabel::OBSpectrophore::AngStepSize2;
	if (!o.compare("5"))  return OpenBabel::OBSpectrophore::AngStepSize5;
	if (!o.compare("10")) return OpenBabel::OBSpectrophore::AngStepSize10;
	if (!o.compare("15")) return OpenBabel::OBSpectrophore::AngStepSize15;
	if (!o.compare("20")) return OpenBabel::OBSpectrophore::AngStepSize20;
	if (!o.compare("30")) return OpenBabel::OBSpectrophore::AngStepSize30;
	if (!o.compare("36")) return OpenBabel::OBSpectrophore::AngStepSize36;
	if (!o.compare("45")) return OpenBabel::OBSpectrophore::AngStepSize45;
	if (!o.compare("60")) return OpenBabel::OBSpectrophore::AngStepSize60;
	return OpenBabel::OBSpectrophore::AngStepSize20;
}



OpenBabel::OBSpectrophore::StereoOption
stringToStereoOption(char* v)
{
	std::string o = v;
	for (unsigned int i = 0; i < o.length(); i++) o[i] = toupper(o[i]);
	if (!o.compare("NO"))     return OpenBabel::OBSpectrophore::NoStereoSpecificProbes;
	if (!o.compare("UNIQUE")) return OpenBabel::OBSpectrophore::UniqueStereoSpecificProbes;
	if (!o.compare("MIRROR")) return OpenBabel::OBSpectrophore::MirrorStereoSpecificProbes;
	if (!o.compare("ALL"))    return OpenBabel::OBSpectrophore::AllStereoSpecificProbes;
	return OpenBabel::OBSpectrophore::NoStereoSpecificProbes;
}



void
showParameters(const OpenBabel::OBSpectrophore& spec, const std::string& ifile)
{
	std::string msg;
	
	std::cout << std::endl;
	std::cout << "*******************************************" << std::endl;
	std::cout << "SPECTROPHORE(TM) CALCULATOR: OBSPECTROPHORE" << std::endl;
	std::cout << "*******************************************" << std::endl;
	std::cout << std::endl;
	
	std::cout << "Input file:       " << ifile << std::endl;
	
	OpenBabel::OBSpectrophore::NormalizationOption n = spec.GetNormalization();
	switch (n)
	{
		case OpenBabel::OBSpectrophore::NoNormalization:
			msg = "No";
			break;
		case OpenBabel::OBSpectrophore::NormalizationTowardsZeroMean:
			msg = "ZeroMean";
			break;
		case OpenBabel::OBSpectrophore::NormalizationTowardsUnitStd:
			msg = "UnitStd";
			break;
		case OpenBabel::OBSpectrophore::NormalizationTowardsZeroMeanAndUnitStd:
			msg = "ZeroMeanAndUnitStd";
			break;
	}
	std::cout << "Normalization:    " << msg << std::endl;
	
	OpenBabel::OBSpectrophore::AccuracyOption a = spec.GetAccuracy();
	switch (a)
	{
		case OpenBabel::OBSpectrophore::AngStepSize1:
			msg = "1";
			break;
		case OpenBabel::OBSpectrophore::AngStepSize2:
			msg = "2";
			break;
		case OpenBabel::OBSpectrophore::AngStepSize5:
			msg = "5";
			break;
		case OpenBabel::OBSpectrophore::AngStepSize10:
			msg = "10";
			break;
		case OpenBabel::OBSpectrophore::AngStepSize15:
			msg = "15";
			break;
		case OpenBabel::OBSpectrophore::AngStepSize20:
			msg = "20";
			break;
		case OpenBabel::OBSpectrophore::AngStepSize30:
			msg = "30";
			break;
		case OpenBabel::OBSpectrophore::AngStepSize36:
			msg = "36";
			break;
		case OpenBabel::OBSpectrophore::AngStepSize45:
			msg = "45";
			break;
		case OpenBabel::OBSpectrophore::AngStepSize60:
			msg = "60";
			break;
	}
	std::cout << "Accuracy:         " << msg << " degrees" << std::endl;
	
	OpenBabel::OBSpectrophore::StereoOption s = spec.GetStereo();
	switch (s)
	{
		case OpenBabel::OBSpectrophore::NoStereoSpecificProbes:
			msg = "No";
			break;
		case OpenBabel::OBSpectrophore::UniqueStereoSpecificProbes:
			msg = "Unique";
			break;
		case OpenBabel::OBSpectrophore::MirrorStereoSpecificProbes:
			msg = "Mirror";
			break;
		case OpenBabel::OBSpectrophore::AllStereoSpecificProbes:
			msg = "All";
			break;
	}
	std::cout << "Stereo treatment: " << msg << std::endl;
	
	double r = spec.GetResolution();
	std::cout << "Resolution:       " << r << " Angstrom" << std::endl;
	
	std::cout << std::endl;
}



void
showError(const std::string& msg)
{
	std::cerr << "\n*******************************************\n";
	std::cerr << "SPECTROPHORE(TM) CALCULATOR: OBSPECTROPHORE\n";
	std::cerr << "*******************************************";
	std::cerr << std::endl;
	std::cerr << std::endl;
   std::cerr << "ERROR: " << msg;
	std::cerr << std::endl;
}



void
showHelp(const char* cmd)
{
	std::cout << "\n*******************************************\n";
	std::cout << "SPECTROPHORE(TM) CALCULATOR: OBSPECTROPHORE\n";
	std::cout << "*******************************************";
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << "DISCLAIMER:" << std::endl;
	std::cout << std::endl;
	std::cout << "Copyright (C) 2005-2010 by Silicos NV" << std::endl;
	std::cout << std::endl;
	std::cout << "This program and file are part of the Open Babel project." << std::endl;
	std::cout << "For more information, see <http://openbabel.org/>" << std::endl;
	std::cout << std::endl;
	std::cout << "This program is free software; you can redistribute it and/or modify" << std::endl;
	std::cout << "it under the terms of the GNU General Public License as published by" << std::endl;
	std::cout << "the Free Software Foundation version 2 of the License." << std::endl;
	std::cout << std::endl;
	std::cout << "The algorithm in this software has been covered by patent WO2009146735." << std::endl;
	std::cout << "However, Silicos NV and the inventors of the above mentioned patent assure " << std::endl;
	std::cout << "that no patent infringment claims will be issued against individuals or" << std::endl;
	std::cout << "institutions that use this software under the GNU General Public License." << std::endl;
	std::cout << std::endl;
	std::cout << "This program is distributed in the hope that it will be useful," << std::endl;
	std::cout << "but WITHOUT ANY WARRANTY; without even the implied warranty of" << std::endl;
	std::cout << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the" << std::endl;
	std::cout << "GNU General Public License for more details." << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << "USAGE:" << std::endl;
	std::cout << std::endl;
	std::cout << cmd << " [parameters]" << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << "PARAMETERS:" << std::endl;
	std::cout << std::endl;
	std::cout << "-i [required]: input file" << std::endl;
	std::cout << "-n [optional]: normalization" << std::endl;
	std::cout << "-a [optional]: accuracy" << std::endl;
	std::cout << "-s [optional]: stereo treatment" << std::endl;
	std::cout << "-r [optional]: resolution" << std::endl;
	std::cout << "-u [optional]: implementation details" << std::endl;
	std::cout << "-h [optional]: help" << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << "PARAMETER DETAILS:" << std::endl;
	std::cout << std::endl;
	std::cout << "-i : Specifies the molecular input file from which Spectrophores(TM)" << std::endl;
	std::cout << "     are to be calculated. The filetype is automatically detected" << std::endl;
	std::cout << "     from the file extension." << std::endl;
	std::cout << std::endl;
	std::cout << "-n : Specifies the kind of normalization that should be performed." << std::endl;
	std::cout << "     Valid values are (without quotes):" << std::endl;
	std::cout << "         No (default)" << std::endl;
	std::cout << "         ZeroMean" << std::endl;
	std::cout << "         UnitStd" << std::endl;
	std::cout << "         ZeroMeanAndUnitStd" << std::endl;
	std::cout << std::endl;
	std::cout << "-a : Specifies the required accuracy expressed as the angular stepsize." << std::endl;
	std::cout << "     Only the following discrete values are allowed:" << std::endl;
	std::cout << "         1, 2, 5, 10, 15, 20 (default), 30, 36, 45, 60" << std::endl;
	std::cout << std::endl;
	std::cout << "-s : Specifies the kind of cages that should be used in terms of" << std::endl;
	std::cout << "     the underlying pointgroup: P1 or P-1. Valid values are (without quotes):" << std::endl;
	std::cout << "         No (default)" << std::endl;
	std::cout << "         Unique" << std::endl;
	std::cout << "         Mirror" << std::endl;
	std::cout << "         All" << std::endl;
	std::cout << std::endl;
	std::cout << "-r : Specifies the required resolution expressed as a real positive number." << std::endl;
	std::cout << "     The default value is 3.0 Angstrom. Negative values or a value equal" << std::endl;
	std::cout << "     to 0 generate an error message" << std::endl;
	std::cout << std::endl;
	std::cout << "-u : Displays implementation details." << std::endl;
	std::cout << std::endl;
	std::cout << "-h : Displays this help file." << std::endl;
	std::cout << std::endl;
}



void
showImplementationDetails(const char* cmd)
{
	std::cout << "\n*******************************************\n";
	std::cout << "SPECTROPHORE(TM) CALCULATOR: OBSPECTROPHORE\n";
	std::cout << "*******************************************";
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << "DISCLAIMER:" << std::endl;
	std::cout << std::endl;
	std::cout << "Copyright (C) 2005-2010 by Silicos NV" << std::endl;
	std::cout << std::endl;
	std::cout << "This program and file are part of the Open Babel project." << std::endl;
	std::cout << "For more information, see <http://openbabel.org/>" << std::endl;
	std::cout << std::endl;
	std::cout << "This program is free software; you can redistribute it and/or modify" << std::endl;
	std::cout << "it under the terms of the GNU General Public License as published by" << std::endl;
	std::cout << "the Free Software Foundation version 2 of the License." << std::endl;
	std::cout << std::endl;
	std::cout << "The algorithm in this software has been covered by patent WO2009146735." << std::endl;
	std::cout << "However, Silicos NV and the inventors of the above mentioned patent assure " << std::endl;
	std::cout << "that no patent infringment claims will be issued against individuals or" << std::endl;
	std::cout << "institutions that use this software under the GNU General Public License." << std::endl;
	std::cout << std::endl;
	std::cout << "This program is distributed in the hope that it will be useful," << std::endl;
	std::cout << "but WITHOUT ANY WARRANTY; without even the implied warranty of" << std::endl;
	std::cout << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the" << std::endl;
	std::cout << "GNU General Public License for more details." << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << "IMPLEMENTATION DETAILS:" << std::endl;
	std::cout << std::endl;
	std::cout << "Spectrophores(TM) are one-dimensional descriptors generated from the property" << std::endl;
	std::cout << "fields surrounding the molecules. The technology allows the accurate description" << std::endl;
	std::cout << "of molecules in terms of their surface properties or fields. Comparison of" << std::endl;
	std::cout << "molecules’ property fields provides a robust structure-independent method of" << std::endl;
	std::cout << "aligning actives from different chemical classes. When applied to molecules such" << std::endl;
	std::cout << "as ligands and drugs, Spectrophores(TM) can be used as powerful molecular " << std::endl;
	std::cout << "descriptors in the fields of chemoinformatics, virtual screening, and QSAR" << std::endl;
	std::cout << "modeling." << std::endl;
	std::cout << "The computation of Spectrophores(TM) is independent of the position and" << std::endl;
	std::cout << "orientation of the molecule and this enables easy and fast comparison of" << std::endl;
	std::cout << "Spectrophores(TM) between different molecules. Molecules having similar" << std::endl;
	std::cout << "three-dimensional properties and shapes always yield similar Spectrophores(TM)." << std::endl;
	std::cout << "A Spectrophore(TM) is calculated by surrounding the three-dimensional" << std::endl;
	std::cout << "conformation of the molecule by a three-dimensional arrangement of points," << std::endl;
	std::cout << "followed by calculating the interaction between each of the atom properties and" << std::endl;
	std::cout << "the surrounding the points. The three-dimensional arrangement of the points" << std::endl;
	std::cout << "surrounding the molecule can be regarded as an ‘artificial’ cage or receptor," << std::endl;
	std::cout << "and the interaction calculated between the molecule and the cage can be regarded" << std::endl;
	std::cout << "as an artificial representation of an affinity value between molecule and cage." << std::endl;
	std::cout << "Because the calculated interaction is dependent on the relative orientation of" << std::endl;
	std::cout << "the molecule within the cage, the molecule is rotated in discrete angles and the" << std::endl;
	std::cout << "most favorable interaction value is kept as final result. The angular stepsize" << std::endl;
	std::cout << "at which the molecule is rotated along its three axis can be specified by the" << std::endl;
	std::cout << "user and influences the accuracy of the method." << std::endl;
	std::cout << std::endl;
	std::cout << "Atomic properties" << std::endl;
	std::cout << std::endl;
	std::cout << "The calculation of a Spectrophore(TM) starts by calculating the atomic" << std::endl;
	std::cout << "contributions of each property from which one wants to calculate a" << std::endl;
	std::cout << "Spectrophore(TM) from. In the current implementation, four atomic properties are" << std::endl;
	std::cout << "converted into a Spectrophore(TM); these four properties include the atomic" << std::endl;
	std::cout << "partial charges, the atomic lipohilicities, the atomic shape deviations and the" << std::endl;
	std::cout << "atomic electrophilicities. The atomic partial charges and atomic electrophilicity" << std::endl;
	std::cout << "properties are calculated using the electronegativity equalisation method (EEM)" << std::endl;
	std::cout << "as described by Bultinck and coworkers (J. Phys. Chem. 2002, A106, 7895-7901;" << std::endl;
	std::cout << "J. Chem. Inf. Comput. Sci. 2003, 43, 422-428). Atomic lipophilic potential" << std::endl;
	std::cout << "parameters are calculated using a rule-based method. Finally, the atomic shape" << std::endl;
	std::cout << "deviation is generated by calculating, for each atom, the atom’s deviation from" << std::endl;
	std::cout << "the average molecular radius. This is done in a four steps process:" << std::endl;
	std::cout << "  - The molecular center of geometry (COG) is calculated;" << std::endl;
	std::cout << "  - The distances between each atom and the molecular COG are calculated;" << std::endl;
	std::cout << "  - The average molecular radius is calculated by averaging all the atomic " << std::endl;
	std::cout << "    distances." << std::endl;
	std::cout << "  - The distances between each atom and the COG are then divided by the average" << std::endl;
	std::cout << "    molecular radius and centered on zero." << std::endl;
	std::cout << std::endl;
	std::cout << "Interaction between the atoms and cage points" << std::endl;
	std::cout << std::endl;
	std::cout << "Following the calculation of all required atomic properties, the next step in" << std::endl;
	std::cout << "the calculation of a Spectrophore(TM) consists of determining the total" << std::endl;
	std::cout << "interaction value V(c,p) between each of the atomic contributions of" << std::endl;
	std::cout << "property p with a set of interaction points on an artificial cage c" << std::endl;
	std::cout << "surrounding the molecular conformation. For this purpose, each of these" << std::endl;
	std::cout << "interaction points i on cage c is assigned a value P(c,i)" << std::endl;
	std::cout << "which is either +1 or -1, with the constraint that the sum of all interaction" << std::endl;
	std::cout << "points on a particular cage should be zero. In a typical Spectrophore(TM)" << std::endl;
	std::cout << "calculation, a cage is represented as a rectangular box encompassing the" << std::endl;
	std::cout << "molecular conformation in all three dimensions, with the centers of the box" << std::endl;
	std::cout << "edges being the interaction points. Such a configuration gives twelve" << std::endl;
	std::cout << "interaction points per cage, and, in the case of a non-stereospecific" << std::endl;
	std::cout << "distribution of the interaction points, leads to 12 different cages. Although" << std::endl;
	std::cout << "there are no particular requirements as to the dimensions of the rectangular" << std::endl;
	std::cout << "cage, the distance between the interaction points and the geometrical extremes" << std::endl;
	std::cout << "of the molecule should be such that a meaningful interaction value between each" << std::endl;
	std::cout << "cage point and the molecular entity can be calculated. In this respect, the" << std::endl;
	std::cout << "default dimensions of the cage are constantly adjusted to enclose the molecule" << std::endl;
	std::cout << "at a minimum distance of 3 A along all dimensions. This cage size can be" << std::endl;
	std::cout << "modified by the user and influences the resolution of the Spectrophore(TM)." << std::endl;
	std::cout << "The total interaction value V(c,p) between the atomic contribution values" << std::endl;
	std::cout << "A(j,p) of property p for a given molecular conformation and the" << std::endl;
	std::cout << "cage interaction values P(c,i) for a given cage c is calculated" << std::endl;
	std::cout << "according a standard interaction energy equation. It takes into account the" << std::endl;
	std::cout << "Euclidean distance between each atom and each cage point. This total interaction" << std::endl;
	std::cout << "V(c,p) for a given property p and cage c for a given" << std::endl;
	std::cout << "molecular conformation is minimized by sampling the molecular orientation along" << std::endl;
	std::cout << "the three axis in angular steps and the calculation of the interaction value for" << std::endl;
	std::cout << "each orientation within the cage. The final total interaction V(c,p) for" << std::endl;
	std::cout << "a given cage c and property p corresponds to the lowest" << std::endl;
	std::cout << "interaction value obtained this way, and corresponds to the c’th value in" << std::endl;
	std::cout << "the one-dimensional Spectrophore(TM) vector calculated for molecular property" << std::endl;
	std::cout << "p. As a result, a Spectrophore(TM) is organized as a vector of minimized" << std::endl;
	std::cout << "interaction values V, each of these organized in order of cages and" << std::endl;
	std::cout << "property values. Since for a typical Spectrophore(TM) implementation twelve" << std::endl;
	std::cout << "different cages are used, the total length of a Spectrophore(TM) vector equals" << std::endl;
	std::cout << "to 12 times the number of properties. Since four different properties are used" << std::endl;
	std::cout << "in the current implementation (electrostatic, lipophilic, electrophilic" << std::endl;
	std::cout << "potentials, and an additional shape index as described before), this leads to a" << std::endl;
	std::cout << "total Spectrophore(TM) length of 48 real values per molecular conformation." << std::endl;
	std::cout << "Since Spectrophore(TM) descriptors are dependent on the actual" << std::endl;
	std::cout << "three-dimensional conformation of the molecule, a typical analysis includes the" << std::endl;
	std::cout << "calculation of Spectrophores(TM) from a reasonable set of different" << std::endl;
	std::cout << "conformations. It is then up to the user to decide on the most optimal strategy" << std::endl;
	std::cout << "for processing the different Spectrophore(TM) vectors. In a typical virtual" << std::endl;
	std::cout << "screening application, calculating the average Spectrophore(TM) vector from all" << std::endl;
	std::cout << "conformations of a single molecule may be a good strategy; other applications" << std::endl;
	std::cout << "have benefit from calculating a weighted average or the minimal values." << std::endl;
	std::cout << "For each molecule in the input file, a Spectrophore(TM) is calculated and" << std::endl;
	std::cout << "printed to standard output as a vector of 48 numbers (in the case of a" << std::endl;
	std::cout << "non-stereospecific Spectrophore(TM). The 48 doubles are organised into 4 sets" << std::endl;
	std::cout << "of 12 doubles each:" << std::endl;
	std::cout << "   - numbers 01-11: Spectrophore(TM) values calculated from the atomic partial charges;" << std::endl;
	std::cout << "   - numbers 13-24: Spectrophore(TM) values calculated from the atomic lipophilicity properties;" << std::endl;
	std::cout << "   - numbers 25-36: Spectrophore(TM) values calculated from the atomic shape deviations;" << std::endl;
	std::cout << "   - numbers 37-48: Spectrophore(TM) values calculated from the atomic electrophilicity properties;" << std::endl;
	std::cout << std::endl;
	std::cout << "Accuracy" << std::endl;
	std::cout << std::endl;
	std::cout << "As already mentioned, the total interaction between cage and molecule for a" << std::endl;
	std::cout << "given property is minimized by sampling the molecular orientation in angular" << std::endl;
	std::cout << "steps of a certain magnitude. As a typical angular step size, 20 degrees was found to" << std::endl;
	std::cout << "be the best compromise between accuracy and computer speed. Larger steps sizes" << std::endl;
	std::cout << "are faster to calculate but have the risk of missing the global interaction" << std::endl;
	std::cout << "energy minimum, while smaller angular steps sizes do sample the rotational space" << std::endl;
	std::cout << "more thoroughly but at a significant computational cost. The accuracy can be" << std::endl;
	std::cout << "specified by the user using the -a option." << std::endl;
	std::cout << std::endl;
	std::cout << "Resolution" << std::endl;
	std::cout << std::endl;
	std::cout << "Spectrophores(TM) capture information about the property fields surrounding the" << std::endl;
	std::cout << "molecule, and the amount of detail that needs to be captured can be regulated by" << std::endl;
	std::cout << "the user. This is done by altering the minimal distance between the molecule and" << std::endl;
	std::cout << "the surrounding cage. The resolution can be specified by the user with the" << std::endl;
	std::cout << "-r option. The default distance along all dimensions is 3.0 Angstrom." << std::endl;
	std::cout << "The larger the distance, the lower the resolution. With a higher resolution," << std::endl;
	std::cout << "more details of the property fields surrounding the molecule are contained by" << std::endl;
	std::cout << "the Spectrophore(TM). On the contrary, low resolution settings may lead to a more" << std::endl;
	std::cout << "general representation of the property fields, with little or no emphasis on" << std::endl;
	std::cout << "small local variations within the fields. Using a low resolution can be the" << std::endl;
	std::cout << "method of choice during the initial virtual screening experiments in order to get" << std::endl;
	std::cout << "an initial, but not so discriminative, first selection. This initial selection" << std::endl;
	std::cout << "can then further be refined during subsequent virtual screening steps using a" << std::endl;
	std::cout << "higher resolution. In this setting, small local differences in the fields between" << std::endl;
	std::cout << "pairs of molecules will be picked up much more easily." << std::endl;
	std::cout << "The absolute values of the individual Spectrophore(TM) data points are dependent" << std::endl;
	std::cout << "on the used resolution. Low resolution values lead to small values of the" << std::endl;
	std::cout << "calculated individual Spectrophore(TM) data points, while high resolutions will" << std::endl;
	std::cout << "lead to larger data values. It is therefore only meaningful to compare only" << std::endl;
	std::cout << "Spectrophores(TM) that have been generated using the same resolution settings or" << std::endl;
	std::cout << "after some kind of normalization is performed." << std::endl;
	std::cout << "Computation time is not influenced by the specified resolution, hence the" << std::endl;
	std::cout << "computation time is identical for all different resolution settings." << std::endl;
	std::cout << std::endl;
	std::cout << "Stereospecificity" << std::endl;
	std::cout << std::endl;
	std::cout << "Some of the cages that are used to calculated Spectrophores(TM) have a" << std::endl;
	std::cout << "stereospecific distribution of the interaction points. The resulting" << std::endl;
	std::cout << "interaction valus resulting from these cages are therefore sensitive to the" << std::endl;
	std::cout << "enantiomeric configuration of the molecule within the cage. The fact that both" << std::endl;
	std::cout << "stereoselective as well as stereo non-selective cages can be used makes it" << std::endl;
	std::cout << "possible to include or exclude stereospecificity in the virtual screening" << std::endl;
	std::cout << "search. Depending on the desired output, the stereospecificity of" << std::endl;
	std::cout << "Spectrophores(TM) can be specified by the user using the -s option:" << std::endl;
	std::cout << "  - No stereospecificity (default). Spectrophores(TM) are generated using cages" << std::endl;
	std::cout << "    that are not stereospecific. For most applications, these Spectrophores(TM)" << std::endl;
	std::cout << "    will suffice." << std::endl;
	std::cout << "  - Unique stereospecificity. Spectrophores(TM) are generated using unique" << std::endl;
	std::cout << "    stereospecific cages." << std::endl;
	std::cout << "  - Mirror stereospecificity. Mirror stereospecific Spectrophores(TM) are" << std::endl;
	std::cout << "    Spectrophores(TM) resulting from the mirror enantiomeric form of the input" << std::endl;
	std::cout << "    molecules." << std::endl;
	std::cout << "The differences between the corresponding data points of unique and mirror" << std::endl;
	std::cout << "stereospecific Spectrophores(TM) are very small and require very long" << std::endl;
	std::cout << "calculation times to obtain a sufficiently high quality level. This increased" << std::endl;
	std::cout << "quality level is triggered by the accuracy setting and will result in" << std::endl;
	std::cout << "calculation times being increased by at least a factor 100. As a consequence, it" << std::endl;
	std::cout << "is recommended to apply this increased accuracy only in combination with a" << std::endl;
	std::cout << "limited number of molecules, and when the small differences between the" << std::endl;
	std::cout << "stereospecific Spectrophores(TM) are really critical. However, for the vast" << std::endl;
	std::cout << "majority of virtual screening applications, this increased accuracy is not" << std::endl;
	std::cout << "required as long as it is not the intention to draw conclusions about" << std::endl;
	std::cout << "differences in the underlying molecular stereoselectivity. Non-stereospecific" << std::endl;
	std::cout << "Spectrophores(TM) will therefore suffice for most applications." << std::endl;
	std::cout << std::endl;
	std::cout << "Normalisation" << std::endl;
	std::cout << std::endl;
	std::cout << "It may sometimes be desired to focus on the relative differences between the" << std::endl;
	std::cout << "Spectrophore(TM) data points rather than focussing on the absolute differences." << std::endl;
	std::cout << "In these cases, normalization of Spectrophores(TM) may be required. The current" << std::endl;
	std::cout << "implementation offers with the -n option the possibility to normalize in four" << std::endl;
	std::cout << "different ways:" << std::endl;
	std::cout << "  - No normalization (default);" << std::endl;
	std::cout << "  - Normalization towards zero mean;" << std::endl;
	std::cout << "  - Normalization towards standard deviation;" << std::endl;
	std::cout << "  - Normalization towards zero mean and unit standard deviation." << std::endl;
	std::cout << "In all these cases, normalization is performed on a ‘per-property’ basis, which" << std::endl;
	std::cout << "means that the data points belonging to the same property set are treated as a" << std::endl;
	std::cout << "single set and that normalization is only performed on the data points within" << std::endl;
	std::cout << "each of these sets and not across all data points." << std::endl;
	std::cout << "Normalization may be important when comparing the Spectrophores(TM) of charged" << std::endl;
	std::cout << "molecules with those of neutral molecules. For molecules carrying a global" << std::endl;
	std::cout << "positive charge, the resulting Spectrophore(TM) data points of the charge and" << std::endl;
	std::cout << "electrophilicity properties will both be shifted in absolute value compared to" << std::endl;
	std::cout << "the corresponding data points of the respective neutral species. Normalization" << std::endl;
	std::cout << "of the Spectrophores(TM) removes the original magnitude differences for the data" << std::endl;
	std::cout << "points corresponding to the charge and electrophilicity properties of charged" << std::endl;
	std::cout << "and neutral species. Therefore, if the emphasis of the virtual screening" << std::endl;
	std::cout << "consists of the identification of molecules with similar property fields without" << std::endl;
	std::cout << "taking into account differences in absolute charge, then Spectrophores(TM)" << std::endl;
	std::cout << "should be normalized towards zero mean. However, if absolute charge differences" << std::endl;
	std::cout << "should be taken into account to differentiate between molecules, unnormalized" << std::endl;
	std::cout << "Spectrophores(TM) are recommended." << std::endl;
	std::cout << std::endl;
}
