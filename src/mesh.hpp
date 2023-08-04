#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <unordered_set>
#include <vector>

#include "bucket.hpp"

namespace mesh {
// ordinal directions
enum Direction : char { NW = 1, NE = 2, SE = 3, SW = 4 };
const std::string dir_string[5] = {"?", "NW", "NE", "SE", "SW"};

class ColNode {
public:
  integer value_j;
  std::shared_ptr<ColNode> north;
  std::shared_ptr<ColNode> south;

  ColNode() = default;
  ~ColNode() = default;
  inline ColNode(const bucket::Point &p) {
    value_j = p[bucket::Y];
    north = nullptr;
    south = nullptr;
  }

  bool isolated(const bucket::Point &, const Direction) const;
};

class ColTree {
public:
  ColTree() { root = nullptr; }
  ~ColTree() = default;

  bool empty() const { return root == nullptr; }
  std::shared_ptr<ColNode> tree() const { return root; }
  void insert(const bucket::Point &);
  bool isolated(const bucket::Point &, const Direction) const;

private:
  std::shared_ptr<ColNode> root;
};

class RowNode {
public:
  integer value_i;
  std::shared_ptr<RowNode> west;
  std::shared_ptr<RowNode> east;
  std::shared_ptr<ColTree> column;

  RowNode() = default;
  ~RowNode() = default;
  RowNode(const bucket::Point &point) {
    value_i = point[bucket::X];
    west = nullptr;
    east = nullptr;
    column = std::make_shared<ColTree>();
    column->insert(point);
  }

  bool isolated(const bucket::Point &, const Direction) const;
};

class RowTree {
public:
  RowTree() { root = nullptr; }
  ~RowTree() = default;

  bool empty() const { return root == nullptr; }
  std::shared_ptr<RowNode> tree() const { return root; }
  void insert(const bucket::Point &);
  void insert(const std::vector<bucket::Point> &);
  bool isolated(const bucket::Point &, const Direction) const;

private:
  std::shared_ptr<RowNode> root;
  void build(const std::vector<bucket::Point> &, const int, const int);
};
typedef std::vector<std::vector<RowTree>> Mesh;
typedef std::vector<std::unordered_set<integer>> Strip;

bucket::Clause isolation(const bucket::Point &p,
                         const std::array<size_t, 2> &xy,
                         const Direction quadrant);

void init(Mesh &);
void init(Strip &);

} // namespace mesh

std::ostream &operator<<(std::ostream &output, const bucket::Point &p);
std::ostream &operator<<(std::ostream &output, const mesh::ColNode &cn);
std::ostream &operator<<(std::ostream &output, const mesh::ColTree &ct);
std::ostream &operator<<(std::ostream &output, const mesh::RowNode &rn);
std::ostream &operator<<(std::ostream &output, const mesh::RowTree &rt);
