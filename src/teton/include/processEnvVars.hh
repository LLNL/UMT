#include "conduit/conduit_node.hpp"

namespace Teton
{

void logOverride(const std::string &envVar,
                 const char *value,
                 bool verbose = false,
                 const std::string &extraMessage = "");
void processEnvVars(conduit::Node &options, bool verbose);

} // namespace Teton
