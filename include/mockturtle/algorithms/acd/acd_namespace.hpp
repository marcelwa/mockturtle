/* Namespace shim for the vendored Ashenhurst-Curtis decomposition sources.
 *
 * The files ac_decomposition.hpp, acd66.hpp, acdXX.hpp and kitty_*.hpp in this
 * directory are copied verbatim from ABC (src/map/if/acd/, author Alessandro
 * Tempia Calvino, EPFL, MIT licence).  They are self-contained C++17 apart from
 * two ABC macros and the ABC scalar type `word`.  Defining the macros to open
 * our own namespace keeps the vendored `kitty` -- a reduced private copy -- from
 * colliding with mockturtle's real `kitty` dependency: everything below lands
 * in `abcacd::kitty` and `abcacd::acd`.
 *
 * This header must be included before any of the vendored files.
 */
#pragma once

#include <cstdint>

#ifndef ABC_NAMESPACE_CXX_HEADER_START
#define ABC_NAMESPACE_CXX_HEADER_START namespace abcacd {
#define ABC_NAMESPACE_CXX_HEADER_END }
#endif

namespace abcacd
{
using word = uint64_t;
}
