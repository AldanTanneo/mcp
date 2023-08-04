#pragma once

#include "mcp-common.hpp"
#include "mcp-matrix+formula.hpp"
#include <set>
#include <string>
#include <unordered_map>
#include <utility>

namespace bucket {

struct Pattern {
  Sign lsign;    // left sign
  Sign rsign;    // right sign
  size_t lcoord; // left coordinate
  size_t rcoord; // right coordinate

  Pattern() = default;
  ~Pattern() = default;
  constexpr Pattern(Sign ls, size_t lc, Sign rs, size_t rc) noexcept
      : lsign(ls), rsign(rs), lcoord(lc), rcoord(rc) {}

  constexpr bool operator==(const Pattern &other) const noexcept {
    return lsign == other.lsign && rsign == other.rsign &&
           lcoord == other.lcoord && rcoord == other.rcoord;
  }
};

typedef std::array<integer, 2> Point;
constexpr size_t X = 0, Y = 1;

struct Clause {
  Sign lsign, rsign;
  integer lval, rval;
  size_t lcoord, rcoord;
};

constexpr bool valid(const Clause &c) noexcept {
  return c.lsign != lnone && c.rsign != lnone;
}

//------------------------------------------------------------------------------

typedef std::unordered_map<Pattern, std::set<Point>> Bucket;

void insert(const Clause &, Bucket &);
Formula get_formula(const Bucket &);
bool sat_bucket(const Row &, const Bucket &);
void print_bucket(const Bucket &);

} // namespace bucket

template <> struct std::hash<bucket::Pattern> {
  using argument_type = bucket::Pattern;
  using result_type = std::size_t;

  result_type operator()(const argument_type &pattern) const {
    size_t left_hash = hash_combine(std::hash<Sign>{}(pattern.lsign),
                                    std::hash<size_t>{}(pattern.lcoord));
    size_t right_hash = hash_combine(std::hash<Sign>{}(pattern.rsign),
                                     std::hash<size_t>{}(pattern.rcoord));
    return hash_combine(left_hash, right_hash);
  }
};
