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

#include "babelconfig.h"
#include <string>
#include <vector>


/**
	* Interface for dynamic libraries. This class defines an interface for
	* finding and opening dynamic libraries on different platforms 
	* (the cpp files are different).
	* It has only what is needed for OpenBabel and is not intended to be 
	* general purpose.
	*/
class DLHandler
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
	  * Returns the number of valid files.  */
	static int findFiles (std::vector <std::string>& file_list, 
			const std::string& pattern, const std::string& path);

	/** Opens a dynamic library */
	static bool openLib(const std::string& lib_name);

	//To select OB format files
	static const char* getFormatFilePattern();

	static char DLHandler::getSeparator();

};

#endif	/* DLHANDLER_H*/
