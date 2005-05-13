/**********************************************************************
Copyright (C) 2004 by Chris Morley

This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#ifndef DLHANDLER_H
#define DLHANDLER_H

#if HAVE_CONFIG_H
#include "config.h"
#endif
#include <string>
#include <vector>

// These macros are used in DLL builds. If they have not
// been set in babelconfig.h, define them as nothing.
#ifndef OBCONV
	#define OBCONV
#endif
#ifndef OBDLL
	#define OBDLL
#endif

/**
	* Interface for dynamic libraries. This class defines an interface for
	* finding and opening dynamic libraries on different platforms 
	* (the cpp files are different).
	* It has only what is needed for OpenBabel and is not intended to be 
	* general purpose.
	*/
class OBCONV DLHandler
{
public:

	/** Provides the path from which the conversion dynamic library,
	  * (OBConv.dll in Windows) was loaded.
	  * This is the default directory for the format files (*.obf in Windows)
	  */ 
	static bool getConvDirectory(std::string& convPath);

	/** Searches a directory specified by path for files whose name matches
	  * a pattern which can include * as a wildcard.
		* The path name should include a final separator (\ or /).
	  * The routine fills a vector of strings with the matching file names (including path).
	  * Note that this is not recursive: it only matches files in the specified path.
	  * For example, if path = e:\\path\\to\\ and pattern = *.obf it will return
	  * vector entries lik e:\\path\\to\\cmlformat.obf
	  * Returns the number of valid files.
	*/
	static int findFiles (std::vector <std::string>& file_list, 
			const std::string& pattern, const std::string& path);

	/** Searches for files which match a full filename (including the path) which
	  * contains a wildcard.
	  * The routine adds matching file names (including path) to a vector of strings .
	  * Returns the number of matching files.
	  * If no wildcard in name adds name to vector and returns -1.
	*/
	static int findFiles (std::vector<std::string>& file_list,const std::string &filename);

	/** Opens a dynamic library */
	static bool openLib(const std::string& lib_name);

	//To select OB format files
	static const char* getFormatFilePattern();

	static char DLHandler::getSeparator();
	static void DLHandler::Sleep(int n);

};

#endif	/* DLHANDLER_H*/
