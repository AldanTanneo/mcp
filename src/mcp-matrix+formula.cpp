/**************************************************************************
 *                                                                        *
 *                                                                        *
 *        Multiple Characterization Problem (MCP)                         *
 *                                                                        *
 * Author:   Miki Hermann                                                 *
 * e-mail:   hermann@lix.polytechnique.fr                                 *
 * Address:  LIX (CNRS UMR 7161), Ecole Polytechnique, France             *
 *                                                                        *
 * Author:   Gernot Salzer                                                *
 * e-mail:   gernot.salzer@tuwien.ac.at                                   *
 * Address:  Technische Universitaet Wien, Vienna, Austria                *
 *                                                                        *
 * Author:   César Sagaert                                                *
 * e-mail:   cesar.sagaert@ensta-paris.fr                                 *
 * Address:  ENSTA Paris, Palaiseau, France                               *
 *                                                                        *
 * Version: all                                                           *
 *     File:    src/mcp-matrix+formula.cpp                                *
 *                                                                        *
 *      Copyright (c) 2019 - 2023                                         *
 *                                                                        *
 **************************************************************************/

#include <cstdint>
#include <iomanip>
#include <iostream>
// #include <sstream>
#include "mcp-matrix+formula.hpp"
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

//------------------------------------------------------------------------------

string version = GLOBAL_VERSION;
const int SENTINEL = -1;
const double RSNTNL = -1.0;
const size_t MTXLIMIT = 4000;

const string print_strg[] = {"void", "clause", "implication", "mixed",
                             "DIMACS"};
const string display_strg[] = {"undefined", "hide", "peek", "section", "show"};
const string sign_string[] = {/*lnone*/ "?",
                              /*lneg*/ "<=", /*lpos*/ ">=", /*lboth*/ "<=>"};

integer DMAX = 0;
Group_of_Matrix group_of_matrix;
vector<string> grps;

string varid = "x";
bool varswitch = false;
vector<string> varnames;

Action action = aALL;
string selected = "";
Print print = pVOID;
Display display = yUNDEF;
string suffix;
size_t arity = 0;
size_t offset = 0;

//------------------------------------------------------------------------------

void Row::inplace_minimum(const Row &other) & {
  for (size_t i = 0; i < other.size(); ++i) {
    if (operator[](i) > other[i]) {
      operator[](i) = other[i];
    }
  }
}

void Row::inplace_minimum(const RowView &other) & {
  for (size_t i = 0; i < other.size(); ++i) {
    if (operator[](i) > other[i]) {
      operator[](i) = other[i];
    }
  }
}

void Row::restrict(const Mask &m) {
  size_t col = 0;
  for (size_t j = 0; j < m.size(); ++j) {
    if (m[j]) {
      // invariant: col <= j
      data[col] = data[j];
      col++;
    }
  }
  data.resize(col);
}

bool Row::operator>=(const Row &rhs) const {
  if (rhs.size() != size())
    return false;

  for (size_t i = 0; i < rhs.size(); ++i) {
    if (rhs[i] > operator[](i))
      return false;
  }
  return true;
}

bool Row::operator>(const Row &rhs) const {
  if (rhs.size() != size())
    return false;

  for (size_t i = 0; i < rhs.size(); ++i) {
    if (rhs[i] >= operator[](i))
      return false;
  }
  return true;
}

bool Row::operator==(const Row &rhs) const {
  if (rhs.size() != size())
    return false;

  for (size_t i = 0; i < rhs.size(); ++i) {
    if (rhs[i] != operator[](i))
      return false;
  }
  return true;
}

bool RowView::operator>=(const RowView &rhs) const {
  if (rhs.size() != size())
    return false;

  for (size_t i = 0; i < rhs.size(); ++i) {
    if (rhs[i] > operator[](i))
      return false;
  }
  return true;
}

bool RowView::operator>=(const Row &rhs) const {
  if (rhs.size() != size())
    return false;

  for (size_t i = 0; i < rhs.size(); ++i) {
    if (rhs[i] > operator[](i))
      return false;
  }
  return true;
}

bool RowView::operator>(const RowView &rhs) const {
  if (rhs.size() != size())
    return false;

  for (size_t i = 0; i < rhs.size(); ++i) {
    if (rhs[i] >= operator[](i))
      return false;
  }
  return true;
}

bool RowView::operator>(const Row &rhs) const {
  if (rhs.size() != size())
    return false;

  for (size_t i = 0; i < rhs.size(); ++i) {
    if (rhs[i] >= operator[](i))
      return false;
  }
  return true;
}

bool RowView::operator==(const RowView &rhs) const {
  if (rhs.size() != size())
    return false;

  for (size_t i = 0; i < rhs.size(); ++i) {
    if (rhs[i] != operator[](i))
      return false;
  }
  return true;
}

bool RowView::operator==(const Row &rhs) const {
  if (rhs.size() != size())
    return false;

  for (size_t i = 0; i < rhs.size(); ++i) {
    if (rhs[i] != operator[](i))
      return false;
  }
  return true;
}

Row RowView::to_row() const {
  Row res(cols.size());
  for (size_t i = 0; i < cols.size(); ++i) {
    res[i] = data[cols[i]];
  }
  return res;
}

void Matrix::delete_row(size_t index) {
  // swap and forget last row
  size_t last = num_rows() - 1;

  if (index != last)
    std::swap(data[index], data[last]);
  data.resize(last);
}

Matrix MatrixMask::to_matrix() const & {
  vector<Row> data;
  data.reserve(num_rows());

  for (size_t i = 0; i < num_rows(); ++i) {
    data.push_back(operator[](i).to_row());
  }

  return Matrix(std::move(data));
}

void Matrix::restrict(const Mask &m) {
  // can be parallelised easily
  // #pragma omp parallel for
  for (size_t i = 0; i < num_rows(); ++i) {
    data[i].restrict(m);
  }
}

size_t partition_matrix(Matrix &mtx, size_t low, size_t high) {
  const Row &pivot = mtx[high];
  size_t p_index = low;

  for (size_t i = low; i < high; i++) {
    if (mtx[i].total_order(pivot) <= 0) {
      std::swap(mtx[i], mtx[p_index]);
      p_index++;
    }
  }
  std::swap(mtx[high], mtx[p_index]);

  return p_index;
}

void sort_matrix(Matrix &mtx, size_t low, size_t high) {
  if (low < high) {
    size_t p_index = partition_matrix(mtx, low, high);
    sort_matrix(mtx, low, p_index - 1);
    sort_matrix(mtx, p_index + 1, high);
  }
}

void Matrix::sort() { sort_matrix((Matrix &)(*this), 0, num_rows() - 1); }

void Matrix::remove_duplicates() {
  auto ip = unique(data.begin(), data.end());
  data.resize(size_t(distance(data.begin(), ip)));
}

//------------------------------------------------------------------------------

// template <typename T>
// bool sat_clause(const T &tuple, const Clause &clause) {
//   // does the tuple satisfy the clause?
//   for (size_t i = 0; i < tuple.size(); ++i) {
//     if (clause[i].sat(tuple[i])) {
//       return true;
//     }
//   }
//   return false;
// }

// template <typename T>
// bool sat_formula(const T &tuple, const Formula &formula) {
//   // does the tuple satisfy the formula?
//   for (Clause cl : formula) {
//     if (!sat_clause(tuple, cl)) {
//       return false;
//     }
//   }
//   return true;
// }

bool sat_clause(const RowView &tuple, const Clause &clause) {
  // does the tuple satisfy the clause?
  for (size_t i = 0; i < tuple.size(); ++i) {
    if (clause[i].sat(tuple[i])) {
      return true;
    }
  }
  return false;
}

bool sat_formula(const RowView &tuple, const Formula &formula) {
  // does the tuple satisfy the formula?
  for (Clause cl : formula) {
    if (!sat_clause(tuple, cl)) {
      return false;
    }
  }
  return true;
}

bool sat_clause(const Row &tuple, const Clause &clause) {
  // does the tuple satisfy the clause?
  for (size_t i = 0; i < tuple.size(); ++i) {
    if (clause[i].sat(tuple[i])) {
      return true;
    }
  }
  return false;
}

bool sat_formula(const Row &tuple, const Formula &formula) {
  // does the tuple satisfy the formula?
  for (Clause cl : formula) {
    if (!sat_clause(tuple, cl)) {
      return false;
    }
  }
  return true;
}

bool sat_clause(const Matrix &matrix, const Clause &clause) {
  for (size_t i = 0; i < matrix.num_rows(); ++i) {
    if (!sat_clause(matrix[i], clause)) {
      return false;
    }
  }
  return true;
}

bool sat_formula(const Matrix &matrix, const Formula &formula) {
  for (size_t i = 0; i < matrix.num_rows(); ++i) {
    if (!sat_formula(matrix[i], formula)) {
      return false;
    }
  }
  return true;
}

vector<string> split(string strg, string delimiters) {
  // splits a string into chunks separated by delimiters (split in perl)
  vector<string> chunks;

  for (size_t i = 0; i < strg.length(); ++i)
    if (!isprint(strg[i])) {
      cerr << "+++ string on input has a non-printable character on position "
           << i << endl;
      exit(2);
    }
  // get rid of non-printable characters at the end of the string without warning
  // because Linux, Windows, and MacOS all terminate the string differently
  while (! strg.empty() && ! isprint(strg.back()))
    strg.pop_back();

  while (!strg.empty()) {
    size_t found = strg.find_first_not_of(delimiters);
    if (found == string::npos)
      break;
    strg.erase(0, found);
    found = strg.find_first_of(delimiters);
    if (found == string::npos) { // strg can be nonempty
      if (strg.length() > 0)
        chunks.push_back(strg);
      break;
    }
    chunks.push_back(strg.substr(0, found));
    strg.erase(0, found);
  }
  return chunks;
}

string clause2dimacs(const vector<size_t> &names, const Clause &clause) {
  // transforms clause into readable clausal form in (extended) DIMACS format to
  // print
  string output = "\t";
  for (size_t lit = 0; lit < clause.size(); ++lit) {
    string var = to_string(offset + names[lit]);
    if (clause[lit].sign & lpos) {
      output += var + ":" + to_string(clause[lit].pval) + " ";
    }
    if (clause[lit].sign & lneg) {
      output += "-" + var + ":" + to_string(clause[lit].nval) + " ";
    }
  }
  output += "0";
  return output;
}

string formula2dimacs(const vector<size_t> &names, const Formula &formula) {
  // transforms formula into readable clausal form in 'extended) DIMACS format
  // to print
  if (formula.empty())
    return " ";

  string output;
  for (Clause clause : formula)
    output += clause2dimacs(names, clause) + "\n";
  return output;
}

// TODO:
// add "boolean" flag / env var, for pretty printing boolean formulas?
// else switch on lit.sign and print var_name >= pval + var_name <= nval, etc
string literal2string(const size_t &litname, const Literal lit) {
  string output;
  if (varswitch) {
    // TODO: fix later
    /* vector<string> new_names = split(varnames[litname], ":");
    if (new_names.size() > 1) // positive or negative (boolean case)
      output +=
          (lit.sign == lneg) ? new_names[nNEGATIVE] : new_names[nPOSITIVE];
    else // own name only
      output += new_names[nOWN] + (lit.sign == lneg ? "<=" : ">="); */
  } else { // variable without name
    string var_name = varid + to_string(offset + litname);
    if (lit.sign & lpos) {
      output += var_name + ">=" + to_string(lit.pval);
    }
    if (lit.sign == lboth) {
      output += " + ";
    }
    if (lit.sign & lneg) {
      output += var_name + "<=" + to_string(lit.nval);
    }
  }
  return output;
}

string rlcl2string(const vector<size_t> &names, const Clause &clause) {
  string output;
  bool plus = false;
  for (size_t lit = 0; lit < clause.size(); ++lit) {
    if (clause[lit].sign != lnone) {
      if (plus == true)
        output += " + ";
      else
        plus = true;
      output += literal2string(names[lit], clause[lit]);
    }
  }
  return output;
}

string impl2string(const vector<size_t> &names, const Clause &clause) {
  string output;
  for (size_t lit = 0; lit < clause.size(); ++lit)
    if (clause[lit].sign & lneg) {
      Literal l = clause[lit];
      l.sign = lpos;
      l.pval = l.nval + 1;
      output += literal2string(names[lit], l) + " ";
    }
  output += "->";
  for (size_t lit = 0; lit < clause.size(); ++lit)
    if (clause[lit].sign & lpos) {
      Literal l = clause[lit];
      l.sign = lpos;
      output += " " + literal2string(names[lit], l);
    }
  return output;
}

string clause2string(const vector<size_t> &names, const Clause &clause) {
  string output = "(";
  if (print == pCLAUSE) {
    output += rlcl2string(names, clause);
  } else if (print == pIMPL) {
    output += impl2string(names, clause);
  } else if (print == pMIX) {
    int pneg = 0;
    int ppos = 0;
    for (size_t lit = 0; lit < clause.size(); ++lit)
      if (clause[lit].sign & lneg)
        pneg++;
      else if (clause[lit].sign & lpos)
        ppos++;
    output += (pneg == 0 || ppos == 0) ? rlcl2string(names, clause)
                                       : impl2string(names, clause);
  }
  output += ')';
  return output;
}

string formula2string(const vector<size_t> &names, const Formula &formula) {
  // transforms formula into readable clausal, implication or mixed form to
  // print
  if (formula.empty())
    return " ";

  if (print == pDIMACS)
    return formula2dimacs(names, formula);

  string output;
  bool times = false;
  for (Clause clause : formula) {
    output += (times == true) ? "\n\t* " : "\t  ";
    times = true;
    output += clause2string(names, clause);
  }
  return output;
}

string literal2latex(const size_t &litname, const Literal lit) {
  string output;
  if (varswitch) {
    vector<string> new_names = split(varnames[litname], ":");
    // TODO: fix later
    /* if (new_names.size() > 1) // positive or negative
      output += lit == lneg && print == pCLAUSE ? new_names[nNEGATIVE]
                                                : new_names[nPOSITIVE];
    else // own name only
      output += (lit == lneg && print == pCLAUSE ? "\\neg " + new_names[nOWN]
                                                 : new_names[nOWN]); */
  } else { // variable without name
    string var_name = varid + to_string(offset + litname);
    if (lit.sign & lpos) {
      output += var_name + "\\geq" + to_string(lit.pval);
    }
    if (lit.sign == lboth) {
      output += " \\lor ";
    }
    if (lit.sign & lneg) {
      output += var_name + "\\leq" + to_string(lit.nval);
    }
  }
  return output;
}

string rlcl2latex(const vector<size_t> &names, const Clause &clause) {
  string output;
  bool lor = false;
  for (size_t lit = 0; lit < clause.size(); ++lit) {
    if (clause[lit].sign != lnone) {
      if (lor == true)
        output += " \\lor ";
      else
        lor = true;
      output += literal2latex(names[lit], clause[lit]);
    }
  }
  return output;
}

string impl2latex(const vector<size_t> &names, const Clause &clause) {
  string output;
  for (size_t lit = 0; lit < clause.size(); ++lit)
    if (clause[lit].sign & lneg) {
      Literal l = clause[lit];
      l.sign = lpos;
      l.pval = l.nval + 1;
      output += literal2latex(names[lit], l) + " ";
    }
  output += "\\to";
  for (size_t lit = 0; lit < clause.size(); ++lit)
    if (clause[lit].sign & lpos) {
      Literal l = clause[lit];
      l.sign = lpos;
      output += " " + literal2latex(names[lit], l);
    }
  return output;
}

string clause2latex(const vector<size_t> &names, const Clause &clause) {
  string output = "\\left(";
  if (print == pCLAUSE) {
    output += rlcl2latex(names, clause);
  } else if (print == pIMPL) {
    output += impl2latex(names, clause);
  } else if (print == pMIX) {
    int pneg = 0;
    int ppos = 0;
    for (size_t lit = 0; lit < clause.size(); ++lit)
      if (clause[lit].sign & lneg)
        pneg++;
      else if (clause[lit].sign & lpos)
        ppos++;
    output += (pneg == 0 || ppos == 0) ? rlcl2latex(names, clause)
                                       : impl2latex(names, clause);
  }
  output += "\\right)";
  return output;
}

string formula2latex(const vector<size_t> &names, const Formula &formula) {
  // transforms formula into readable clausal form in LaTeX format to print
  if (formula.empty())
    return " ";

  string output;
  bool land = false;
  for (Clause clause : formula) {
    output += (land == true) ? "\n\t\\land " : "\t  ";
    land = true;
    output += clause2latex(names, clause);
  }
  return output;
}

void read_formula(vector<size_t> &names, Formula &formula) {
  // formula read instructions

  int nvars;
  cin >> suffix >> arity >> nvars >> offset;

  // cerr << "*** suffix = " << suffix
  //      << ", arity = " << arity
  //      << ", #vars = " << nvars
  //      << ", offset = " << offset
  //      << endl;

  vector<int> validID;
  int dummy;
  for (int i = 0; i < nvars; ++i) {
    cin >> dummy;
    validID.push_back(dummy);
  }

  for (size_t i = 0; i < arity; ++i)
    names.push_back(i);

  string lit;
  Clause clause(arity, Literal::none());
  while (cin >> lit) {
    if (lit == "0") { // end of clause in DIMACS
      formula.push_back(clause);
      for (size_t i = 0; i < arity; ++i)
        clause[i] = Literal::none();
    } else {
      vector<string> parts = split(lit, ":");
      if (parts.size() == 0 || parts.size() > 2) {
        cerr << "+++ " << quoted(lit)
             << " is not a valid Extended DIMACS literal" << endl;
        exit(2);
      }
      try {
        long long var = stoll(parts[0]);

        if (find(cbegin(validID), cend(validID), abs(var)) == cend(validID)) {
          cerr << "+++ " << abs(var) << " outside allowed variable names"
               << endl;
          exit(2);
        }

        unsigned long long val =
            parts.size() == 2 ? stoull(parts[1]) : (var < 0 ? 0 : 1);
        if (val > (unsigned long long)(std::numeric_limits<integer>::max())) {
          throw out_of_range(parts[1]);
        }

        Literal l = clause.at(size_t(abs(var)) - 1 - offset);
        if (var < 0) {
          l.sign = Sign(l.sign | lneg);
          l.nval = integer(val);
        } else {
          l.sign = Sign(l.sign | lpos);
          l.pval = integer(val);
        }
        clause.at(size_t(abs(var)) - 1 - offset) = l;

      } catch (invalid_argument const &ex) {
        cerr << "+++ The Extended DIMACS literal " << quoted(lit)
             << " contains invalid integer values: " << ex.what() << endl;
        exit(2);
      } catch (out_of_range const &ex) {
        cerr << "+++ The Extended DIMACS literal " << quoted(lit)
             << " contains an integer value that falls out of range: "
             << ex.what() << endl;
        exit(2);
      }
    }
  }
}

ostream &operator<<(ostream &output, const Row &row) {
  // overloading ostream to print a row
  // transforms a tuple (row) to a printable form
  // for (bool bit : row)
  for (size_t i = 0; i < row.size(); ++i)
    // output << to_string(bit); // bit == true ? 1 : 0;
    // output << to_string(row[i]);
    output << row[i] << " ";
  return output;
}

ostream &operator<<(ostream &output, const RowView &row) {
  // overloading ostream to print a row
  // transforms a tuple (row) to a printable form
  // for (bool bit : row)
  for (size_t i = 0; i < row.size(); ++i)
    // output << to_string(bit); // bit == true ? 1 : 0;
    // output << to_string(row[i]);
    output << row[i] << " ";
  return output;
}

template <typename M>
ostream &internal_matrix_display(ostream &output, const M &m) {
  // overloading ostream to print a matrix
  // transforms a matrix to a printable form
  const size_t n = m.num_rows();
  const size_t max_to_show = (display == yPEEK && n >= 5) ? 5 : n;
  for (size_t i = 0; i < max_to_show; ++i) {
    output << "\t";
    for (size_t j = 0; j < m.num_cols(); ++j) {
      output << m[i][j] << " ";
    }
    output << "\n";
  }
  if (display == yPEEK && n >= 5)
    output << "\t...\n" << endl;
  return output;
}

ostream &operator<<(ostream &output, const Matrix &M) {
  return internal_matrix_display<Matrix>(output, M);
}
ostream &operator<<(ostream &output, const MatrixMask &M) {
  return internal_matrix_display<MatrixMask>(output, M);
}
