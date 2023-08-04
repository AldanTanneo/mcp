/**************************************************************************
 *                                                                        *
 *                                                                        *
 *	           Multiple Characterization Problem (MCP)                    *
 *                                                                        *
 *	Author:   Miki Hermann                                                *
 *	e-mail:   hermann@lix.polytechnique.fr                                *
 *	Address:  LIX (CNRS UMR 7161), Ecole Polytechnique, France            *
 *                                                                        *
 *	Author: Gernot Salzer                                                 *
 *	e-mail: gernot.salzer@tuwien.ac.at                                    *
 *	Address: Technische Universitaet Wien, Vienna, Austria                *
 *                                                                        *
 *	Version: all                                                          *
 *      File:    mcp-matrix+formula.hpp                                   *
 *                                                                        *
 *      Copyright (c) 2019 - 2023                                         *
 *                                                                        *
 * Data structures for row, matrices, literals, clauses, and formula.     *
 *                                                                        *
 * This software has been created within the ACCA Project.                *
 *                                                                        *
 *                                                                        *
 **************************************************************************/

#pragma once

#include <deque>
#include <functional>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#define GLOBAL_VERSION "1.05d-mekong-"

//------------------------------------------------------------------------------

extern std::string version;
extern const std::string print_strg[];
extern const std::string display_strg[];

// domain value type
using integer = uint16_t;
// maximum value of the domain (cardinality - 1)
constexpr integer DCARD = std::numeric_limits<integer>::max();

class RowView;

// mask into a row
using Mask = std::vector<bool>;

// data structure representing a matrix row
class Row {
public:
  using container = std::vector<integer>;

private:
  container data;

public:
  inline explicit Row() = default;
  inline explicit Row(size_t size) : data(container(size)) {}

  // constructs a row from rvalue data
  inline explicit Row(container &&data) : data(data) {}

  Row(const Row &) = delete;
  Row(Row &&) = default;
  Row &operator=(Row &&) = default;
  // explicitely clone the row
  inline Row clone() const & { return Row(container(data)); }

  // get the size of the row
  constexpr size_t size() const noexcept { return data.size(); }
  // constant index operation
  constexpr const integer &operator[](size_t index) const & noexcept {
    return data[index];
  }
  // mutable index operation
  constexpr integer &operator[](size_t index) & noexcept { return data[index]; }
  // add a value to the end of the row
  inline void push_back(integer value) & { data.push_back(value); }
  // preallocate enough memory for the specified number of values
  inline void reserve(size_t size) & { data.reserve(size); }
  // resize the row to the given size
  inline void resize(size_t size) & { data.resize(size); }

  // set the row to the (element-wise) minimum between itself and other
  void inplace_minimum(const Row &other) &;
  // set the row to the (element-wise) minimum between itself and other
  void inplace_minimum(const RowView &other) &;

  // total alphabetical order on rows.
  // - if self < other, return -1
  // - if self == other, return 0
  // - if self > other, return 1
  inline int total_order(const Row &other) const & {
    size_t i = 0;
    for (; i < size() && i < other.size(); ++i) {
      if (data[i] < other[i])
        return -1;
      if (data[i] > other[i])
        return 1;
    }
    if (i < size())
      return 1;
    if (i < other.size())
      return -1;
    return 0;
  }

  // restricts the row to the given set of columns
  void restrict(const Mask &m);

  bool operator>=(const Row &) const;
  bool operator>(const Row &) const;
  bool operator==(const Row &) const;
};

// basic matrix class
class Matrix {
public:
  using container = std::vector<Row>;

private:
  // the actual data contained in the matrix
  container data;

public:
  // empty matrix
  explicit Matrix() = default;
  // constructs a matrix from rvalue data
  inline explicit Matrix(container &&data) : data(data) {}

  Matrix(const Matrix &) = delete;
  Matrix(Matrix &&) = default;
  Matrix &operator=(Matrix &&) = default;
  // explicitely clone the matrix
  inline Matrix clone() const & { return Matrix(container(data)); }

  // checks wether the matrix is empty
  constexpr bool empty() const noexcept { return data.empty(); }
  // reserves space for the rows
  inline void reserve(size_t size) { data.reserve(size); }
  // returns the number of rows
  constexpr size_t num_rows() const noexcept { return data.size(); }
  // returns the number of columns
  constexpr size_t num_cols() const noexcept {
    return data.size() > 0 ? data[0].size() : 0;
  }

  // equivalent to M[row][col]
  constexpr integer get(size_t row, size_t col) const noexcept {
    return data[row][col];
  }
  // returns a const view into a row
  inline const Row &operator[](size_t index) const { return data[index]; }

  // returns a mutable view into a row
  inline Row &operator[](size_t index) { return data[index]; }

  // add a new row to the matrix
  inline void add_row(Row new_row) { data.push_back(new_row); }

  // deletes a row without preserving row order
  void delete_row(size_t index);

  // restricts the matrix to the given set of columns
  void restrict(const Mask &m);
  // sort the matrix according to the total order on rows
  void sort();
  // remove duplicate rows in a sorted matrix. useful after a restriction.
  void remove_duplicates();

  friend class MatrixMask;
};

// immutable view into a masked row: behaves as if only a subset of columns were
// present
class RowView {
public:
  using permutation = std::vector<size_t>;

private:
  const Row &data;
  const permutation &cols;

public:
  // construct a new masked row from the original data and the column mask
  inline explicit RowView(const Row &data, const permutation &cols)
      : data(data), cols(cols) {}

  // get the size of the masked row
  inline size_t size() const { return cols.size(); }

  // equivalent to R[index] with the mask applied.
  inline integer operator[](size_t index) const { return data[cols[index]]; }

  bool operator>=(const RowView &) const;
  bool operator>(const RowView &) const;
  bool operator==(const RowView &) const;
  bool operator==(const Row &) const;

  // copies the masked row to a new one.
  Row to_row() const;
};

// immutable view into a masked matrix: behave as if only a subset of columns
// were present
class MatrixMask {
public:
  using permutation = RowView::permutation;

private:
  permutation cols;
  std::reference_wrapper<const Matrix> matrix;

public:
  inline explicit MatrixMask(const Matrix &mat, const Mask &column_mask)
      : cols(std::vector<size_t>()), matrix(mat) {
    for (size_t i = 0; i < column_mask.size(); ++i) {
      if (column_mask[i]) {
        cols.push_back(i);
      }
    }
  }
  inline explicit MatrixMask(const Matrix &mat)
      : cols(permutation()), matrix(mat) {
    cols.reserve(mat.num_cols());
    for (size_t i = 0; i < mat.num_cols(); ++i) {
      cols.push_back(i);
    }
  }
  inline explicit MatrixMask(const MatrixMask &other)
      : cols(other.cols), matrix(other.matrix) {}

  inline MatrixMask &operator=(const MatrixMask &&other) {
    MatrixMask mat(other);
    std::swap(mat, *this);
    return (MatrixMask &)(*this);
  }

  // hides a column in the mask without erasing it from the original matrix
  inline void hide_column(size_t index) {
    cols.erase(std::begin(cols) + index);
  }

  // checks if the underlying matrix is empty
  constexpr bool empty() const { return matrix.get().empty(); }
  // returns the number of rows
  constexpr size_t num_rows() const { return matrix.get().num_rows(); }
  // returns the number of columns
  constexpr size_t num_cols() const { return cols.size(); }

  // equivalent to M[row][col], with the mask applied
  constexpr integer get(size_t row, size_t col) const {
    return matrix.get().data[row][cols[col]];
  }
  // returns a const view into a masked row
  inline RowView operator[](size_t index) const {
    return RowView(matrix.get().data[index], cols);
  }

  // copies the masked matrix to a new one.
  Matrix to_matrix() const &;
};

#define CARDlimit 50

// sign enumeration
// - lneg is `x <= nval`
// - lpos is `x >= pval`
// - lboth is `x <= nval || x >= pval`
// - lnone means no literal.
enum Sign : char { lnone = 0b00, lneg = 0b01, lpos = 0b10, lboth = 0b11 };
// string representations of Sign
extern const std::string sign_string[];

// variable <= nval or variable >= pval (or both)
struct Literal {
  Sign sign;
  integer pval, nval;

  constexpr Literal() noexcept : sign(lnone), pval(0), nval(0) {}
  ~Literal() = default;

  constexpr explicit Literal(Sign s, integer pv, integer nv)
      : sign(s), pval(pv), nval(nv) {}

  static constexpr Literal none() { return Literal(lnone, 0, 0); }
  static constexpr Literal neg(integer nval) { return Literal(lneg, 0, nval); }
  static constexpr Literal pos(integer pval) { return Literal(lpos, pval, 0); }
  static constexpr Literal both(integer pval, integer nval) {
    return Literal(lboth, pval, nval);
  }

  // inverts the literal.
  // `x <= d` becomes `x >= d + 1`.
  // `x >= d` becomes `x <= d - 1`.
  inline Literal swap() const noexcept {
    Literal res;
    if (sign & lneg && nval < DCARD) {
      res.sign = lpos;
      res.pval = nval + 1;
    }
    if (sign & lpos && pval > 0) {
      res.sign = Sign(res.sign | lneg);
      res.nval = pval - 1;
    }
    return res;
  }

  constexpr bool operator==(const Literal &other) const {
    return this->sign == other.sign && this->pval == other.pval &&
           this->nval == other.nval;
  }

  constexpr bool operator!=(const Literal &other) const {
    return this->sign != other.sign || this->pval != other.pval ||
           this->nval != other.nval;
  }

  constexpr bool operator<(const Literal &other) const {
    return this->sign < other.sign ||
           (this->sign == other.sign && this->nval < other.nval) ||
           (this->sign == other.sign && this->pval < other.pval);
  }

  // check wether a value satisfies a literal
  constexpr bool sat(integer val) const {
    return (sign & lneg && val <= nval) || (sign & lpos && val >= pval);
  }
};

// clause type, a disjunction of literals.
// the literal at index `i` is related to variable `x_i`
using Clause = std::vector<Literal>;
// formula type, a conjunction of clauses.
using Formula = std::deque<Clause>;

typedef std::map<std::string, Matrix> Group_of_Matrix;
extern Group_of_Matrix group_of_matrix;
extern std::vector<std::string> grps;

enum Action : char { aONE = 0, aALL = 1, aNOSECT = 2, aSELECTED = 3 };
enum Print : char { pVOID = 0, pCLAUSE = 1, pIMPL = 2, pMIX = 3, pDIMACS = 4 };
enum Display : char {
  yUNDEF = 0,
  yHIDE = 1,
  yPEEK = 2,
  ySECTION = 3,
  ySHOW = 4
};

const std::string print_string[] = {"void", "clause", "implication", "mix",
                                    "dimacs"};

extern std::string varid;
extern bool varswitch;
extern std::vector<std::string> varnames;
enum NAME : char { nOWN = 0, nPOSITIVE = 1, nNEGATIVE = 2 };

extern const int SENTINEL;
extern const double RSNTNL;
extern const int MTXLIMIT;

extern Action action;
extern std::string selected;
extern std::string suffix;
extern int arity;
extern int offset;
extern Print print;
extern Display display;

//------------------------------------------------------------------------------

// read a matrix from a CSV input
void read_matrix(Group_of_Matrix &matrix);
// print a matrix
void print_matrix(const Group_of_Matrix &matrix);
// read a formula from its Extended DIMACS representation
void read_formula(std::vector<int> &names, Formula &formula);
// split a string along the specified delimiters
std::vector<std::string> split(std::string strg, std::string delimiters);
// get the Extended DIMACS representation of a formula
std::string formula2dimacs(const std::vector<int> &names,
                           const Formula &formula);
// get the string representation of a formula
std::string formula2string(const std::vector<int> &names,
                           const Formula &formula);
// get the latex representation of a formula
std::string formula2latex(const std::vector<int> &names,
                          const Formula &formula);
// checks that a row satisfies a clause
bool sat_clause(const RowView &tuple, const Clause &clause);
// checks that a row satisfies a formula
bool sat_formula(const RowView &tuple, const Formula &formula);

// checks that a row satisfies a clause
bool sat_clause(const Row &tuple, const Clause &clause);
// checks that a row satisfies a formula
bool sat_formula(const Row &tuple, const Formula &formula);

// checks that all rows in a matrix satisfy a clause
bool sat_clause(const Matrix &matrix, const Clause &clause);
// checks that all rows in a matrix satisfy a formula
bool sat_formula(const Matrix &matrix, const Formula &formula);

// display a row
std::ostream &operator<<(std::ostream &output, const Row &row);
// display a row view
std::ostream &operator<<(std::ostream &output, const RowView &row);
// display a matrix
std::ostream &operator<<(std::ostream &output, const Matrix &M);

// combine two hash values
constexpr size_t hash_combine(size_t a, size_t b) {
  if (sizeof(size_t) >= 8)
    a ^= b + 0x517cc1b727220a95 + (a << 6) + (a >> 2);
  else
    a ^= b + 0x9e3779b9 + (a << 6) + (a >> 2);
  return a;
}

// specialization of hash for Row
template <> class std::hash<Row> {
public:
  size_t operator()(const Row &r) const noexcept {
    size_t res = 0;
    for (size_t i = 0; i < r.size(); ++i) {
      hash_combine(res, std::hash<integer>{}(r[i]));
    }
    return res;
  }
};

// specialization of hash for RowView
template <> class std::hash<RowView> {
public:
  size_t operator()(const RowView &r) const noexcept {
    size_t res = 0;
    for (size_t i = 0; i < r.size(); ++i) {
      hash_combine(res, std::hash<integer>{}(r[i]));
    }
    return res;
  }
};
