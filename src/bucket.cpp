
#include "mcp-common.hpp"
#include "mcp-matrix+formula.hpp"

using GeneralClause = Clause;

#include "bucket.hpp"

#include <algorithm>
#include <iostream>

using namespace std;

//------------------------------------------------------------------------------

namespace bucket {

void insert(const Clause &c, Bucket &B) {

  Pattern pattern(c.lsign, c.lcoord, c.rsign, c.rcoord);
  Point point{c.lval, c.rval};

  if (B.count(pattern) > 0) {
    auto it = B[pattern].begin();
    while (it != B[pattern].end()) {
      if (*it == point)
        return;                                 // already there
      if (c.lsign == lpos && c.rsign == lpos) { // SW
        if ((*it)[X] <= point[X] && (*it)[Y] <= point[Y]) {
          it = B[pattern].erase(it);
          continue;
        } else if (point[X] <= (*it)[X] && point[Y] <= (*it)[Y])
          return;
      } else if (c.lsign == lneg && c.rsign == lpos) { // SE
        if ((*it)[X] >= point[X] && (*it)[Y] <= point[Y]) {
          it = B[pattern].erase(it);
          continue;
        } else if (point[X] >= (*it)[X] && point[Y] <= (*it)[Y])
          return;
      } else if (c.lsign == lpos && c.rsign == lneg) { // NW
        if ((*it)[X] <= point[X] && (*it)[Y] >= point[Y]) {
          it = B[pattern].erase(it);
          continue;
        } else if (point[X] <= (*it)[X] && point[Y] >= (*it)[Y])
          return;
      } else if (c.lsign == lneg && c.rsign == lneg) { // NE
        if ((*it)[X] >= point[X] && (*it)[Y] >= point[Y]) {
          it = B[pattern].erase(it);
          continue;
        } else if (point[X] >= (*it)[X] && point[Y] >= (*it)[Y])
          return;
      } else {
        cerr << "+++ bucket insert: you should not be here" << endl;
        exit(1);
      }
      it++;
    }
    if (X == Y) {
      it = B[pattern].begin();
      while (it != B[pattern].end()) {
        if ((*it)[X] < point[X] && point[X] < (*it)[Y] && (*it)[Y] < point[Y])
          point[X] = (*it)[X]; // stretch to left
        else if (point[X] < (*it)[X] && (*it)[X] < point[Y] &&
                 point[Y] < (*it)[Y])
          point[Y] = (*it)[Y]; // stretch to right
        else {
          it++;
          continue;
        }
        it = B[pattern].erase(it);
      }
    }
  }
  B[pattern].insert(point);
}

Formula get_formula(size_t arity, const Bucket &B) {
  Formula f;
  for (const auto &b : B) {
    Pattern pattern = b.first;

    for (const Point &point : b.second) {
      GeneralClause c(arity, Literal(lnone, 0, 0));

      c[pattern.rcoord].sign = Sign(c[pattern.rcoord].sign | pattern.rsign);
      c[pattern.lcoord].sign = Sign(c[pattern.lcoord].sign | pattern.lsign);

      if (pattern.rsign & lpos)
        c[pattern.rcoord].pval = point[0];
      else
        c[pattern.rcoord].nval = point[0];

      if (pattern.lsign & lpos)
        c[pattern.lcoord].pval = point[1];
      else
        c[pattern.lcoord].nval = point[1];

      f.push_back(std::move(c));
    }
  }
  return f;
}

// Bucket is a formula encoded differently.
// Test the satisfiability of each clause.
bool sat_bucket(const Row &t, const Bucket &B) {
  for (const auto &b : B) {
    Pattern pattern = b.first;
    for (const Point &point : b.second) {
      if ((pattern.lsign == lpos && t[pattern.lcoord] < point[0]) &&
          (pattern.lsign == lneg && t[pattern.lcoord] > point[0]) &&
          (pattern.rsign == lpos && t[pattern.rcoord] < point[1]) &&
          (pattern.rsign == lneg && t[pattern.rcoord] > point[1]))
        return false;
    }
  }
  return true;
}

} // namespace bucket

string to_string(const bucket::Point &p) {
  return "[" + to_string(p[0]) + "," + to_string(p[1]) + "]";
}

string to_string(const bucket::Pattern &p) {
  return "[x_" + to_string(p.lcoord) + sign_string[p.lsign] + " , " + "x_" +
         to_string(p.rcoord) + sign_string[p.rsign] + "]";
}

void bucket::print_bucket(const Bucket &bucket) {
  cout << "*** Buckets:" << endl;
  for (const auto &b : bucket) {
    Pattern pattern = b.first;
    cout << "... pattern " << to_string(pattern) << ":";
    for (const Point &point : b.second)
      cout << " (" << point[0] << ", " << point[1] << ")";
    cout << endl;
  }
  cout << endl;
}
