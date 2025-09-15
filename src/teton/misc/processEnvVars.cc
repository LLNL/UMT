#include "processEnvVars.hh"
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

namespace Teton
{

void logOverride(const std::string &envVar, const char *value, bool verbose, const std::string &extraMessage)
{
   if (verbose)
   {
      std::cerr << "Teton: Overriding " << envVar << " to " << value << std::endl;
      if (!extraMessage.empty())
      {
         std::cout << extraMessage << std::endl;
      }
   }
}

void processEnvVars(conduit::Node &options, bool verbose)
{
   // Enable use of environment variables to override Teton behavior.
   const char *enableEnvVars = getenv("TETON_ENABLE_ENV_VARS");
   if (enableEnvVars && atoi(enableEnvVars) > 0)
   {
      if (verbose)
      {
         std::cout << "Teton: Enabling use of environment variables to override runtime behavior. "
                   << "These may change or be deprecated at any time. For developer use only.\n"
                   << "Supported env vars are:\n"
                   << "TETON_VERBOSE\nTETON_DUMP_INPUT_AT_CYCLE\nTETON_DUMP_VIZ\nTETON_DUMP_METRICS\n"
                   << "TETON_NUM_CPU_THREADS\nTETON_DEVICE_NUM_PROCESSORS\nTETON_SWEEP_MAX_HYPERDOMAINS\n"
                   << "TETON_GTA_MAX_HYPERDOMAINS\nTETON_MIN_GSET_SIZE\nTETON_SWEEP_VERSION\n"
#if !defined(TETON_ENABLE_MINIAPP_BUILD)
                   << "TETON_OPERATOR_SPLIT_PDV\n"
#endif
                   << std::endl;
      }

      // Map of environment variables to their corresponding options keys.
      const std::map<std::string, std::string> envVarMap = {
         {"TETON_VERBOSE", "verbose"},
         {"TETON_DUMP_INPUT_AT_CYCLE", "dump_input_at_cycle"},
         {"TETON_DUMP_VIZ", "dump_viz"},
         {"TETON_DUMP_METRICS", "dump_metrics"},
         {"TETON_NUM_CPU_THREADS", "concurrency/omp_cpu_max_threads"},
         {"TETON_DEVICE_NUM_PROCESSORS", "concurrency/omp_device_num_processors"},
         {"TETON_MIN_GSET_SIZE", "sweep/min_group_set_size"},
         {"TETON_SWEEP_VERSION", "sweep/kernel/version"},
         {"TETON_SWEEP_MAX_HYPERDOMAINS", "sweep/max_hyperdomains"},
         {"TETON_GTA_MAX_HYPERDOMAINS", "gta/max_hyperdomains"}};

      for (const auto &entry : envVarMap)
      {
         const char *envValue = getenv(entry.first.c_str());
         if (envValue)
         {
            options[entry.second] = atoi(envValue);
            logOverride(entry.first, envValue, verbose);
         }
      }

#if !defined(TETON_ENABLE_MINIAPP_BUILD)
      const char *operatorSplitPdv = getenv("TETON_OPERATOR_SPLIT_PDV");
      if (operatorSplitPdv)
      {
         options["operator_split_pdv"] = atoi(operatorSplitPdv);
         logOverride(
            "TETON_OPERATOR_SPLIT_PDV",
            operatorSplitPdv,
            verbose,
            "Teton will now call teton_applypdv to apply the PdV work in an operator-split explicit step before each time step.");
      }
#endif
   }
}

} // namespace Teton
