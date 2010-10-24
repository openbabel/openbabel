#include <openbabel/plugin.h>

using namespace OpenBabel;

int main()
{
  //
  // Print out all plugins sorted by plugin type.
  //
  std::cout << "This executable contains the following plugins:" << std::endl;
  // Get the plugin types.
  std::vector<std::string> pluginTypes;
  OBPlugin::ListAsVector("plugins", 0, pluginTypes);

  for (std::size_t i = 0; i < pluginTypes.size(); ++i) {
    std::cout << std::endl << "  " << pluginTypes[i] << ":" << std::endl << std::endl;

    // Get the plugins for the current plugin types
    std::vector<std::string> plugins;
    OBPlugin::ListAsVector(pluginTypes[i].c_str(), 0, plugins);
    for (std::size_t j = 0; j < plugins.size(); ++j)
      std::cout << "    " << plugins[j] << std::endl;

  }

  return 0;
}

