//*** Completely untested ***

#include "dlhandler.h"

#include <unistd.h>
#include <dirent.h>
#include <sys/stat.h>

//Globals for scandir()
string targetPattern;
int matchFiles (const struct dirent *entry_p)
{
	return !strcmp(entry_p->d_name, targetPattern)
};

bool DLHandler::getConvDirectory(std::string& convPath)
{
	//Need to provide the directory from which this shared library was loaded.
	//This is the default directory for format shared library files.
	return false;
}

int DLHandler::findFiles (std::vector <std::string>& file_list,const std::string& pattern, const std::string& path);
{
	if(path.empty()) path="./"; //use current dir ectory if path is empty

	targetPattern=pattern; //make accessible to global function

	struct dirent **entries_pp;

	int count = scandir (path.c_str(), &entries_pp, matchFiles, NULL);

	for(int i=0;i<count;i++)
		file_list.push_back(path + (*entries_pp)->d_name);
	return count;
}

bool DLHandler::openLib(const std::string& lib_name)
{
	return dlopen(lib_name, RTLD_LAZY) != 0;
}

const char* DLHandler::getFormatFilePattern()
{
	return "lib*.so.*";
}

char DLHandler::getSeparator()
{ return '/';}