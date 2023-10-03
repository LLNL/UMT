
#ifndef TetonSurfaceTallies_HH__
#define TetonSurfaceTallies_HH__

#include "conduit/conduit.hpp"
#include <string>

namespace TetonSurfaceTallies
{

void dumpTallyToJson(const conduit::Node &blueprint, const conduit::Node &options, int mpi_rank = 0);

}

#endif // TetonSurfaceTallies_HH__
