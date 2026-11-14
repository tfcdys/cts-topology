// File-based FLUTE LUT build: the base64 LUT variables are not used, but the
// generic initLUT() helper still references them at link time.
#include <string>
namespace Flute {
std::string powv9;
std::string post9;
}  // namespace Flute
