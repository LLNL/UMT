#include "utilities/DBC.hh"
namespace dbc {
static bool lock = false;

bool 
assertionLock (void) {
    if (lock) return false;
    lock = true;
    return true;
}

void 
assertionUnLock (void) {
    lock = false;
}

} // namespace dbc
