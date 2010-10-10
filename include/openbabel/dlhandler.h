/**********************************************************************
dlhandler.h - Dynamic loader for file format modules.

Copyright (C) 2004-2005 by Chris Morley

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#ifndef OB_DLHANDLER_H
#define OB_DLHANDLER_H

#include <openbabel/babelconfig.h>

#include <string>
#include <vector>

// These macros are used in DLL builds. If they have not
// been set in babelconfig.h, define them as nothing.
#ifndef OBERROR
	#define OBERROR
#endif
#ifndef OBDLL
	#define OBDLL
#endif

/** \class DLHandler dlhandler.h <openbabel/dlhandler.h>
    \brief Interface for dynamic libraries.

    This class defines an interface for finding and opening dynamic
    loadable libraries on different platforms (e.g., modular plugins)
    via different source code files.
    It has only what is needed for OpenBabel and is not intended to be
    general purpose. Internally, it is used for dynamic loading and unloading
    OBFormat file translation modules.
**/
class OBERROR DLHandler
{
public:

	/**
          * Get the directory containing the dynamic library plugins.
          *
          * The result is stored in @p convPath and depends on the operating
          * system and build tools.
          *
          * Linux: always returns OB_MODULE_PATH (defined in src/config.h.cmake).
          *
          * Windows cygwin: Uses windows' GetModuleHandle to get a handle to the
          * application module. From this, the applications path can be extracted.
          * The appname.exe is removed and ..\\lib\\openbabel\\BABEL_VERSION\\ is appended.
          *
          * Windows MSVC: Uses windows' GetModuleHandle to get a handle to the
          * OBError dll module. From this, the dll's path can be extracted.
          *
          * @sa findFiles
	  */
	static bool getConvDirectory(std::string& convPath);

	/** Searches a directory specified by path for files whose name matches
	  * a pattern which can include * as a wildcard.
          * If the BABEL_LIBDIR environment variable is set, this will override
          * the @p path parameter.
          * The path name should include a final separator (\ or /).
	  * The routine fills a vector of strings with the matching file names (including path).
	  * Note that this is not recursive: it only matches files in the specified path.
	  * For example, if path = e:\\path\\to\\ and pattern = *.obf it will return
	  * vector entries lik e:\\path\\to\\cmlformat.obf
	  * \return the number of valid files.
	  */
	static int findFiles (std::vector <std::string>& file_list,
			const std::string& pattern, const std::string& path);

	/** Searches for files which match a full filename (including the path) which
	  * contains a wildcard.
	  * The routine adds matching file names (including path) to a vector of strings .
	  * \return the number of matching files.
	  * If no wildcard in name adds name to vector and returns -1.
	*/
	static int findFiles (std::vector<std::string>& file_list,const std::string &filename);

	//! Open a dynamic library with path @p lib_name
	static bool openLib(const std::string& lib_name);

	//! \return The file extension pattern for Open Babel plugin modules (e.g. *.obf on Windows)
	static const char* getFormatFilePattern();

  //! \return The system directory separator (i.e. "\" on Windows, "/" on UNIX)
	static char getSeparator();

  //! Call the system routine to wait (sleep) this process for @p n seconds
	static void Sleep(int n);

};

#endif	/* DLHANDLER_H*/

//! \file dlhandler.h
//! \brief Dynamic loader for file format modules.
