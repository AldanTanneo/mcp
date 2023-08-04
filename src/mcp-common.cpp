/**************************************************************************
 *                                                                        *
 *                                                                        *
 *	       Multiple Characterization Problem (MCP)                    *
 *                                                                        *
 *	Author:   Miki Hermann                                            *
 *	e-mail:   hermann@lix.polytechnique.fr                            *
 *	Address:  LIX (CNRS UMR 7161), Ecole Polytechnique, France        *
 *                                                                        *
 *	Author: Gernot Salzer                                             *
 *	e-mail: gernot.salzer@tuwien.ac.at                                *
 *	Address: Technische Universitaet Wien, Vienna, Austria            *
 *                                                                        *
 *	Version: common to all                                            *
 *      File:    mcp-common.cpp                                           *
 *                                                                        *
 *      Copyright (c) 2019 - 2023                                         *
 *                                                                        *
 * Given several  Boolean matrices  representing the  presence/absence of *
 * attributes in  observations, this software generates  Horn, dual Horn, *
 * or bijunctive formulas from positive and negative examples represented *
 * by these matrices.                                                     *
 *                                                                        *
 * This software has been created within the ACCA Project.                *
 *                                                                        *
 *                                                                        *
 **************************************************************************/

#include "mcp-common.hpp"
#include "mcp-matrix+formula.hpp"
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <unordered_set>

using namespace std;

// string version = GLOBAL_VERSION;
bool debug = false;
// string varid = "x";
// bool varswitch = false;
// vector<string> varnames;

map<Row, int> pred;        // predecessor function for Zanuttini's algorithm
map<Row, int> succ;        // successor function for Zanuttini's algorithm
map<Row, vector<int>> sim; // sim table for Zanuttini's algorithm

// const int SENTINEL     = -1;
const string STDIN = "STDIN";
const string STDOUT = "STDOUT";
// const int MTXLIMIT     = 4000;

// Action action       = aALL;
Closure closure = clHORN;
Cooking cooking = ckWELLDONE;
Direction direction = dBEGIN;
// Print print         = pVOID;
bool setcover = true;
Strategy strategy = sLARGE;
// Display display     = yUNDEF;
string input = STDIN;
string output = STDOUT;
bool disjoint = true;
// int arity = 0;

// int offset          = 0;
string tpath = "/tmp/"; // directory where the temporary files will be stored
bool np_fit = false;
int chunkLIMIT = 4096; // heavily hardware dependent; must be optimized
string latex = "";     // file to store latex output

ifstream infile;
ofstream outfile;
ofstream latexfile;
string formula_output; // prefix of files, where formulas will be stored

const string action_strg[] = {"One to One", "One to All Others",
                              "One to All Others, Nosection",
                              "Selected to All Others"};
const string closure_strg[] = {"Horn", "dual Horn", "bijunctive", "affine",
                               "CNF"};
const string cooking_strg[] = {"raw", "bleu", "medium", "well done"};
const string direction_strg[] = {
    "begin", "end", "optimum", "random", "low cardinality", "high cardinality"};
const string pcl_strg[] = {"Horn", "Horn", "bijunctive", "affine", "cnf"};
const string strategy_strg[] = {"large", "exact"};
// const string print_strg[]     = {"void",       "clause",     "implication",
// "mixed",   "DIMACS"}; const string display_strg[]   = {"undefined",  "hide",
// "peek",        "section", "show"};
const string arch_strg[] = {"seq", "mpi", "pthread", "hybrid"};

//--------------------------------------------------------------------------------------------------

void read_arg(int argc, char *argv[]) { // reads the input parameters
  int argument = 1;
  while (argument < argc) {
    string arg = argv[argument];
    if (arg == "--action" || arg == "-a") {
      string act = argv[++argument];
      if (act == "one" || act == "1") {
        action = aONE;
      } else if (act == "all" || act == "a") {
        action = aALL;
      } else if (act == "nosection" || act == "nosect" || act == "nos" ||
                 act == "no" || act == "ns" || act == "n") {
        action = aNOSECT;
      } else if (act == "selected" || act == "select" || act == "sel" ||
                 act == "s") {
        if (argument < argc - 1) {
          action = aSELECTED;
          selected = argv[++argument];
        } else
          cerr << "+++ no group selected, revert to default" << endl;
      } else
        cerr << "+++ unknown action " << act << endl;
    } else if (arg == "-d" || arg == "--direction") {
      string dir = argv[++argument];
      if (dir == "begin" || dir == "b") {
        direction = dBEGIN;
      } else if (dir == "end" || dir == "e") {
        direction = dEND;
      } else if (dir == "lowcard" || dir == "lcard" || dir == "lc") {
        direction = dLOWCARD;
      } else if (dir == "highcard" || dir == "hcard" || dir == "hc") {
        direction = dHIGHCARD;
      } else if (dir == "random" || dir == "rand" || dir == "r") {
        direction = dRAND;
        // } else if (dir == "optimum"
        // 		 || dir == "optimal"
        // 		 || dir == "opt"
        // 		 || dir ==  "x") {
        // 	direction = dOPT;
      } else
        cerr << "+++ unknown direction option " << dir << endl;
    } else if (arg == "--closure" || arg == "-clo") {
      string cl = argv[++argument];
      if (cl == "horn" || cl == "Horn" || cl == "HORN" || cl == "h") {
        closure = clHORN;
      } else if (cl == "dhorn" || cl == "dHorn" || cl == "dualHorn" ||
                 cl == "dual-Horn" || cl == "dual_Horn" || cl == "dh") {
        closure = clDHORN;
      } else if (cl == "bij" || cl == "bijunctive" || cl == "b") {
        closure = clBIJUNCTIVE;
      } else if (cl == "general" || cl == "gen" || cl == "cnf" || cl == "CNF") {
        closure = clCNF;
      } else
        cerr << "+++ unknown closure option " << cl << endl;
    } else if (arg == "-pr" || arg == "--print") {
      string prt = argv[++argument];
      if (prt == "clause" || prt == "clausal" || prt == "cl" || prt == "c") {
        print = pCLAUSE;
      } else if (prt == "implication" || prt == "impl" || prt == "imp" ||
                 prt == "im" || prt == "i") {
        print = pIMPL;
      } else if (prt == "mix" || prt == "mixed" || prt == "m") {
        print = pMIX;
      } else if (prt == "dimacs" || prt == "DIMACS") {
        print = pDIMACS;
      } else
        cerr << "+++ unknown print option " << prt << endl;
    } else if (arg == "-s" || arg == "--strategy") {
      string strtgy = argv[++argument];
      if (strtgy == "e" || strtgy == "ex" || strtgy == "exact") {
        strategy = sEXACT;
      } else if (strtgy == "l" || strtgy == "lg" || strtgy == "large") {
        strategy = sLARGE;
      } else
        cerr << "+++ unknown strategy option " << strtgy << endl;
    } else if (arg == "-sc" || arg == "--setcover" || arg == "--SetCover") {
      string sc = argv[++argument];
      if (sc == "y" || sc == "Y" || sc == "yes" || sc == "YES") {
        setcover = true;
      } else if (sc == "n" || sc == "N" || sc == "no" || sc == "NO") {
        setcover = false;
      } else
        cerr << "+++ unknown set cover option " << sc << endl;
    } else if ((arch == archMPI || arch == archHYBRID) &&
               (arg == "-f" || arg == "--fit")) {
      string fnp = argv[++argument];
      if (fnp == "y" || fnp == "Y" || fnp == "yes" || fnp == "YES") {
        np_fit = true;
      } else if (fnp == "n" || fnp == "N" || fnp == "no" || fnp == "NO") {
        np_fit = false;
      } else
        cerr << "+++ unknown fit number of processes option " << fnp << endl;
    } else if (arg == "-ck" || arg == "--cook" || arg == "--cooking") {
      string ck = argv[++argument];
      if (ck == "r" || ck == "raw") {
        cooking = ckRAW;
      } else if (ck == "b" || ck == "bl" || ck == "bleu") {
        cooking = ckBLEU;
      } else if (ck == "m" || ck == "med" || ck == "medium") {
        cooking = ckMEDIUM;
      } else if (ck == "wd" || ck == "done" || ck == "well" ||
                 ck == "well_done" || ck == "welldone" || ck == "all") {
        cooking = ckWELLDONE;
      } else
        cerr << "+++ unknown cooking option " << ck << endl;
    } else if (arg == "--input" || arg == "-i") {
      input = argv[++argument];
    } else if (arg == "--output" || arg == "-o") {
      output = argv[++argument];
    } else if (arg == "--formula" || arg == "--logic" || arg == "-l") {
      formula_output = argv[++argument];
    } else if (arg == "--matrix" || arg == "--mtx" || arg == "-m") {
      string mtx = argv[++argument];
      if (mtx == "yes" || mtx == "y" || mtx == "show") {
        display = ySHOW;
      } else if (mtx == "peek") {
        display = yPEEK;
      } else if (mtx == "section") {
        display = ySECTION;
      } else if (mtx == "no" || mtx == "n" || mtx == "hide") {
        display = yHIDE;
      } else if (mtx == "undefined" || mtx == "undef" || mtx == "u") {
        display = yUNDEF;
      } else
        cerr << "+++ unknown matrix print option " << mtx << endl;
    } else if (arg == "--offset" || arg == "-of" || arg == "--shift" ||
               arg == "-sh") {
      offset = stoi(argv[++argument]);
    } else if (arch > archMPI && (arg == "--chunk" || arg == "-ch")) {
      chunkLIMIT = stoi(argv[++argument]);
    } else if (arch != archSEQ && (arg == "--tpath" || arg == "-tp")) {
      tpath = argv[++argument];
    } else if (arg == "--latex") {
      latex = argv[++argument];
    } else if (arg == "--debug") {
      debug = true;
    } else
      cerr << "+++ unknown option " << arg << endl;
    ++argument;
  }
}

int hamming_distance(const Row &u, const Row &v) {
  // Hamming distance between two tuples
  if (u.size() != v.size())
    return SENTINEL;

  int sum = 0;
  for (int i = 0; i < u.size(); ++i)
    sum += abs((long int)u[i] - (long int)v[i]);
  return sum;
}

// is the tuple row in the Horn closure of matrix M?
// TODO: REWRITE in a zero-copy way in
// [x] mpi
// [ ] posix
// [x] seq
/*
bool InHornClosure(const RowView &row, const Matrix &M) {
  Matrix P = ObsGeq(row, M);

  if (P.empty()) {
    return false;
  } else if (row == MIN(P)) {
    return true;
  } else {
    return false;
  }
}
*/

// is the intersection of F and of the Horn closure of T empty?
//
// T = MinimizeObs(T);	// Optional, may not be worth the effort
bool SHCPsolvable(const MatrixMask &T, const MatrixMask &F) {
  for (size_t i = 0; i < F.num_rows(); ++i) {
    if (InHornClosure(F[i], T))
      return false;
  }
  return true;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// checks if the intersection of T and F is non empty
bool isect_nonempty(const MatrixMask &T, const MatrixMask &F) {
  unordered_set<RowView> orig{};

  // insert rowviews into the hashset...
  for (size_t i = 0; i < T.num_rows(); ++i) {
    orig.insert(T[i]);
  }

  // check for the presence of any row from F in the hashset
  for (size_t i = 0; i < F.num_rows(); ++i) {
    if (orig.find(F[i]) != orig.end())
      return true;
  }

  return false;
}

bool inadmissible(const MatrixMask &T, const MatrixMask &F) {
  if (closure < clBIJUNCTIVE)
    return !SHCPsolvable(T, F);
  else
    return isect_nonempty(T, F);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// TODO: statistics heuristic? -> Z score, Median Absolute Deviation...
// These heuristics need global data to be calculated, however.

int hamming_weight(const vector<bool> &row) { // Hamming weight of a tuple
  int sum = accumulate(cbegin(row), cend(row), 0);
  return sum;
}

// computes the minimal section for Horn or dual Horn closures
Mask minsect(Matrix &T, Matrix &F) {
  const int n = T.num_cols();

  if (inadmissible(T, F)) {
    disjoint = false;
    Mask emptymask(n, false);
    return emptymask;
  }

  // boolean mask
  Mask mask(n, true);
  MatrixMask Tm(T), Fm(F);

  vector<size_t> coords(n);

  switch (direction) {
  case dBEGIN:
    for (size_t i = 0; i < n; ++i) {
      coords[i] = n - 1 - i;
    }
    break;
  case dEND:
    for (size_t i = 0; i < n; ++i) {
      coords[i] = i;
    }
    break;
  case dRAND:
    for (size_t i = 0; i < n; ++i) {
      coords[i] = i;
    }
    std::random_shuffle(coords.begin(), coords.end());
    break;
  case dLOWCARD:
  case dHIGHCARD:
    std::cerr << "Lowcard / Highcard sorting need to be replaced with "
                 "statistic heuristics (Z score, etc)"
              << std::endl;
    std::cerr << "Terminating." << std::endl;
    exit(2);
    break;
  case dOPT:
    std::cerr << "Unsupported direction: dOPT" << std::endl;
    exit(2);
    break;
  }

  for (int i : coords) {
    mask[i] = false;
    Tm.hide_column(i);
    Fm.hide_column(i);
    if (inadmissible(Tm, Fm)) {
      mask[i] = true;
      Tm = MatrixMask(T, mask);
      Fm = MatrixMask(F, mask);
    }
  }

  return mask;
}

void w_f(const string &filename, const string suffix, const vector<int> &names,
         const Formula &formula) {
  // open the file and write the formula in it
  ofstream formfile;
  formfile.open(filename);
  if (!formfile.is_open()) {
    cerr << "+++ Cannot open formula output file " << filename << endl;
    cerr << "+++ Formula not written" << endl;
  } else {
    formfile << suffix << " " << arity << " " << formula.cbegin()->size() << " "
             << offset << endl;
    int old_offset = offset;
    offset = 1;
    for (int n : names)
      formfile << " " << n + offset;
    formfile << endl;
    formfile << formula2dimacs(names, formula) << endl;
    offset = old_offset;
    formfile.close();
  }
}

void write_formula(const string &suffix1, const string &suffix2,
                   const vector<int> &names, const Formula &formula) {
  // write formula to a file in DIMACS format
  // offset begins at 1, if not set otherwise
  w_f(formula_output + "_" + suffix1 + "_" + suffix2 + ".log", suffix1, names,
      formula);
}

void write_formula(const string &suffix, const vector<int> &names,
                   const Formula &formula) {
  // write formula to a file in DIMACS format
  // offset begins at 1, if not set otherwise
  w_f(formula_output + "_" + suffix + ".log", suffix, names, formula);
}

bool satisfied_by(const Clause &clause, const Matrix &T) {
  // is the clause satified by all tuples in T?
  for (size_t i = 0; i < T.num_rows(); ++i) {
    bool satisfied = false;
    for (int j = 0; j < clause.size(); ++j) {
      if (clause[j].sat(T.get(i, j))) {
        satisfied = true;
        break;
      }
    }
    if (!satisfied)
      return false;
  }
  return true;
}

int numlit(const Clause &clause) {
  // number of literals in a clause
  int i = 0;
  for (Literal lit : clause)
    switch (lit.sign) {
    case lneg:
    case lpos:
      i++;
      break;
    case lboth:
      i += 2;
      break;
    case lnone:
      break;
    }
  return i;
}

int firstlit(const Clause &clause) {
  // coordinate of first literal
  int i = 0;
  while (i < clause.size() && clause[i].sign == lnone)
    i++;
  return i;
}

bool clauseLT(const Clause &a, const Clause &b) {
  // is clause a < clause b in literals ?
  for (int i = 0; i < a.size(); ++i)
    if (a[i] < b[i])
      return true;
  return false;
}

// This overloading is necessay because deque implements >= differently
bool clauseGE(const Clause &a, const Clause &b) {
  // overloading >=
  // is clause a >= clause b?
  // order on clauses:
  // 1. number of literals
  // 2. coordinate of the first literal
  // 3. order on positive / negative literals
  if (a == b)
    return true;
  if (numlit(a) < numlit(b))
    return false;
  if (numlit(a) == numlit(b) && firstlit(a) < firstlit(b))
    return false;
  if (numlit(a) == numlit(b) && firstlit(a) == firstlit(b))
    return clauseLT(b, a);
  return true;
}

int partition_formula(Formula &formula, int low, int high) {
  Clause pivot = formula[high];
  int p_index = low;

  for (int i = low; i < high; i++)
    if (pivot >= formula[i]) {
      Clause t = formula[i];
      formula[i] = formula[p_index];
      formula[p_index] = t;
      p_index++;
    }
  Clause t = formula[high];
  formula[high] = formula[p_index];
  formula[p_index] = t;

  return p_index;
}

void sort_formula(Formula &formula, int low, int high) {
  if (low < high) {
    int p_index = partition_formula(formula, low, high);
    sort_formula(formula, low, p_index - 1);
    sort_formula(formula, p_index + 1, high);
  }
}

// restricts matrix A to columns determined by the bitvector sect
void restrict(const Mask &sect, Matrix &A) {
  A.restrict(sect);
  A.sort();
  A.remove_duplicates();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Row minHorn(const Matrix &M) {
  // computes the minimal tuple of a matrix coordinate wise
  Row mh = M[0].clone();
  for (int i = 1; i < M.num_rows(); ++i)
    mh.inplace_minimum(M[i]);
  return mh;
}

Formula unitres(const Formula &formula) { // unit resolution
  Formula units;
  Formula clauses;
  /*
  for (Clause cl : formula)
    if (numlit(cl) == 1)
      units.push_back(cl);
    else
      clauses.push_back(cl);

  Formula resUnits = units;
  while (!units.empty() && !clauses.empty()) {
    Clause &unit = units.front();
    units.pop_front();
    int index = 0;
    while (unit[index] == lnone)
      index++;
    for (int j = 0; j < clauses.size(); j++) {
      Clause clause = clauses[j];
      if ((unit[index] == lpos && clause[index] == lneg) ||
          (unit[index] == lneg && clause[index] == lpos))
        clauses[j][index] = lnone;
    }

    auto it = clauses.begin();
    while (it != clauses.end())
      if (numlit(*it) == 1) {
        units.push_back(*it);
        resUnits.push_back(*it);
        it = clauses.erase(it);
      } else
        ++it;
  }

  // sort(resUnits.begin(), resUnits.end());
  sort_formula(resUnits, 0, resUnits.size() - 1);
  auto last1 = unique(resUnits.begin(), resUnits.end());
  resUnits.erase(last1, resUnits.end());

  // test if there are no two unit clauses having literals of opposite parity
  int ru_bound = resUnits.size();
  for (int i = 0; i < ru_bound - 1; ++i) {
    Clause unit_i = resUnits[i];
    int index = 0;
    while (index < unit_i.size() && unit_i[index] == lnone)
      index++;
    if (index < ru_bound) {
      for (int j = i + 1; j < ru_bound; ++j) {
        Clause unit_j = resUnits[j];
        if (unit_i[index] == lpos && unit_j[index] == lneg ||
            unit_i[index] == lneg && unit_j[index] == lpos) {
          Clause emptyClause(unit_i.size(), lnone);
          Formula emptyFormula;
          emptyFormula.push_back(emptyClause);
          return emptyFormula;
        }
      }
    }

  }


  clauses.insert(clauses.end(), resUnits.begin(), resUnits.end());
  // sort(clauses.begin(), clauses.end()), cmp_numlit;
  sort_formula(clauses, 0, clauses.size() - 1);
  auto last2 = unique(clauses.begin(), clauses.end());
  clauses.erase(last2, clauses.end());
  */

  return clauses;
}

bool subsumes(const Clause &cla, const Clause &clb) {
  // does clause cla subsume clause clb ?
  // cla must be smaller than clb
  for (int i = 0; i < cla.size(); ++i)
    if ((cla[i].sign & clb[i].sign & lneg && cla[i].nval < clb[i].nval) ||
        (cla[i].sign & clb[i].sign & lpos && cla[i].pval > clb[i].pval))
      return false;
  return true;
}

// perform subsumption on clauses of a formula
// clauses must be sorted by length --- IS GUARANTEED
Formula subsumption(Formula &formula) {
  Formula res;
  // sort_formula not necessary, clauses are GUARANTEED SORTED
  // sort(formula.begin(), formula.end(), cmp_numlit);
  // sort_formula(formula, 0, formula.size()-1);
  while (!formula.empty()) {
    Clause clause = formula.front();
    formula.pop_front();
    auto it = formula.begin();
    while (it != formula.end())
      if (subsumes(clause, *it))
        it = formula.erase(it);
      else
        ++it;
    res.push_back(clause);
  }
  return res;
}

bool empty_clause(const Clause &clause) { // is the clause empty ?
  for (Literal lit : clause)
    if (lit.sign != lnone)
      return false;
  return true;
}

Formula redundant(const Formula &formula) { // eliminating redundant clauses
                                            // clauses must be sorted by length
                                            // --- IS GUARANTEED
  const int lngt = formula[0].size();

  Formula prefix, suffix;
  // prefix.insert(prefix.end(), formula.begin(), formula.end());
  prefix = formula;
  // sort_formula not necessary, clauses are GUARANTEED SORTED
  // sort(prefix.begin(), prefix.end(), cmp_numlit);
  // sort_formula(prefix, 0, prefix.size()-1);

  int left = 0;
  while (left < prefix.size() && numlit(prefix[left]) == 1)
    left++;

  while (left < prefix.size()) {
    Clause pivot = prefix.back();
    prefix.pop_back();
    Formula newUnits;
    for (int i = 0; i < pivot.size(); ++i) {
      if (pivot[i].sign != lnone) {
        Clause newclause(lngt, Literal::none());
        newclause[i] = pivot[i].swap();
        newUnits.push_back(newclause);
      }
    }
    newUnits.insert(newUnits.end(), prefix.begin(), prefix.end());
    newUnits.insert(newUnits.end(), suffix.begin(), suffix.end());
    // sort(newUnits.begin(), newUnits.end(), cmp_numlit);
    sort_formula(newUnits, 0, newUnits.size() - 1);
    auto last3 = unique(newUnits.begin(), newUnits.end());
    newUnits.erase(last3, newUnits.end());
    Formula bogus = unitres(newUnits);
    bool keep = true;
    for (Clause bgcl : bogus)
      if (empty_clause(bgcl)) {
        keep = false;
        break;
      }
    if (keep)
      suffix.push_front(pivot);
  }
  prefix.insert(prefix.end(), suffix.begin(), suffix.end());
  return prefix;
}

Formula SetCover(const Matrix &Universe, const Formula &SubSets) {
  // perform set cover optimizing the clauses as subsets falsified by tuples as
  // universe Universe = tuples in F SubSets  = clauses of a formula
  enum Presence { NONE = 0, ABSENT = 1, PRESENT = 2 };
  typedef pair<reference_wrapper<const Row>, reference_wrapper<const Clause>>
      Falsification;
  // Incidence matrix indicating which tuple (row) falsifies which clause
  map<Falsification, Presence> incidence;
  vector<size_t> R; // tuples still active == not yet falsified
  Formula selected; // selected clauses for falsification

  for (size_t i = 0; i < Universe.num_rows(); ++i) {
    const Row &tuple = Universe[i];
    R.push_back(i);
    for (const Clause &clause : SubSets)
      incidence[make_pair(cref(tuple), cref(clause))] =
          sat_clause(tuple, clause) ? ABSENT : PRESENT;
  }

  // perform set cover
  while (!R.empty()) {
    // How many tuples does the clause falsify (intersect)?
    map<reference_wrapper<const Clause>, int> intersect;
    for (const Clause &clause : SubSets)
      intersect[clause] = 0;
    for (size_t i : R) {
      const Row &tuple = Universe[i];
      for (const Clause &clause : SubSets)
        if (incidence[make_pair(cref(tuple), cref(clause))] == PRESENT)
          ++intersect[clause];
    }

    int maxw = 0;
    Clause maxset;
    for (const Clause &clause : SubSets) {
      if (intersect[clause] > maxw) {
        maxw = intersect[clause];
        maxset = clause;
      }
    }

    if (maxw == 0)
      break;
    selected.push_back(maxset);

    auto it = R.begin();
    while (it != R.end()) {
      const Row &tuple = Universe[*it];
      if (incidence[make_pair(cref(tuple), cref(maxset))] == PRESENT) {
        for (const Clause &clause : SubSets)
          incidence[make_pair(cref(tuple), cref(clause))] = ABSENT;
        it = R.erase(it);
      } else
        ++it;
    }
  }
  // sort(selected.begin(), selected.end(), cmp_numlit);
  sort_formula(selected, 0, selected.size() - 1);
  return selected;
}

void cook(Formula &formula) {
  // perform redundancy elimination on formula according to cooking
  if (!formula.empty()) {
    if (cooking == ckRAW)
      sort_formula(formula, 0, formula.size() - 1);
    // sort(formula.begin(), formula.end(), cmp_numlit);
    if (cooking >= ckBLEU) {
      formula = unitres(formula);
      // formula = binres(formula);
    }
    if (cooking >= ckMEDIUM)
      formula = subsumption(formula);
    if (cooking == ckWELLDONE)
      formula = redundant(formula);
  }
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void predecessor(const Matrix &R) {
  // predecessor function of Zanuttini's algorithm
  // R must be lexicographically sorted
  pred.clear();
  pred[R[0]] = SENTINEL;
  auto p = cref(R[0]);
  for (int i = 1; i < R.num_rows(); ++i) {
    const Row &m = R[i];
    int j = 0;
    for (int k = 0; k < m.size(); ++k)
      if (m[k] == p.get()[k])
        j++;
      else
        break;
    pred[R[i]] = j;
    p = m;
  }
}

void successor(const Matrix &R) {
  // sucessor function of Zanuttini's algorithm
  // R must be lexicographically sorted
  succ.clear();
  auto m = cref(R[0]);
  for (int i = 1; i < R.num_rows(); ++i) {
    const Row &s = R[i];
    int j = 0;
    for (int k = 0; k < m.get().size(); ++k)
      if (m.get()[k] == s[k])
        j++;
      else
        break;
    succ[R[i - 1]] = j;
    m = s;
  }
  succ[R[R.num_rows() - 1]] = SENTINEL;
}

void simsim(const Matrix &R) { // sim array of Zanuttini's algorithm
  sim.clear();
  for (size_t i = 0; i < R.num_rows(); ++i) {
    const Row &mm = R[i];
    const vector<int> dummy(mm.size(), SENTINEL);
    sim[mm] = dummy;
    for (size_t k = 0; k < R.num_rows(); ++k) {
      const Row &m1m = R[k];
      if (mm == m1m)
        continue;
      int j0 = 0;
      while (mm[j0] == m1m[j0])
        ++j0;
      int j = j0;
      while (j < m1m.size() && (m1m[j] == true || mm[j] == false)) {
        if (j > succ[mm] && mm[j] == false && m1m[j] == true)
          sim[mm][j] = max(sim[mm][j], j0);
        j++;
      }
    }
  }
}

// TODO: what the frick
Clause hext(const Row &m, const int &j) {
  // generate clauses with Zanuttini's algorithm
  Clause clause(m.size(), Literal::none());
  /*
  for (int i = 0; i < j; ++i)
    if (m[i] == true)
      clause[i] = lneg;
  if (j > pred[m] && m[j] == true)
    clause[j] = lpos;
  else if (j > succ[m] && m[j] == false && sim[m][j] == SENTINEL)
    clause[j] = lneg;
  else if (j > succ[m] && m[j] == false && sim[m][j] != SENTINEL) {
    clause[j] = lneg;
    clause[sim[m][j]] = lpos;
  }
  */
  return clause;
}

// TODO: What the frick
// phi is a CNF formula and M is a set of tuples, such that sol(phi) = M
// constructs a reduced prime formula phiPrime, such that sol(phi) =
// sol(phiPrime)
Formula primality(const Formula &phi, const Matrix &M) {
  Formula phiPrime;
  /*
  const int card = M.num_rows();
  const int lngt = M.num_cols();
  auto last = make_unique<int[]>(card);

  for (Clause clause : phi) {
    if (clause.size() != lngt) {
      cerr << "+++ Clause size and vector length do not match" << endl;
      exit(2);
    }
    for (int k = 0; k < card; ++k) {
      last[k] = SENTINEL;
      for (int j = 0; j < clause.size(); ++j)
        if (M[k][j] == true && clause[j] == lpos ||
            M[k][j] == false && clause[j] == lneg)
          last[k] = j;
    }
    Clause cPrime(clause.size(), lnone);
    for (int j = 0; j < clause.size(); ++j) {
      bool d;
      if (clause[j] == lpos) {
        d = true;
        for (int k = 0; k < card; ++k)
          if (last[k] == j)
            d = d && M[k][j];
        if (d == true)
          cPrime[j] = lpos;
      } else if (clause[j] == lneg) {
        d = false;
        for (int k = 0; k < card; ++k)
          if (last[k] == j)
            d = d || M[k][j];
        if (d == false)
          cPrime[j] = lneg;
      }
      if (cPrime[j] != lnone)
        for (int k = 0; k < card; ++k)
          for (int i = 0; i < lngt; ++i)
            if (M[k][i] == true && cPrime[j] == lpos ||
                M[k][i] == false && cPrime[j] == lneg) {
              last[k] = SENTINEL;
              break;
            }
    }
    phiPrime.push_back(cPrime);
  }
  */
  return phiPrime;
}

// TODO: What the frick
// learn the exact Horn clause from positive examples T
// uses Zanuttini's algorithm
Formula learnHornExact(const Matrix &T) {
  Formula H;

  const int lngt = T.num_cols();
  if (T.num_rows() == 1) { // T has only one row / tuple
    const Row &t = T[0];
    for (int i = 0; i < lngt; ++i) {
      Clause clause(lngt, Literal::none());
      clause[i] = Literal::pos(t[i]);
      H.push_back(clause);
      clause[i] = Literal::neg(t[i]);
      H.push_back(clause);
    }
    return H;
  }

  /*
  // What
  // The
  // Frick
  sort_matrix(T, 0, T.num_rows() - 1);
  successor(T);
  predecessor(T);
  simsim(T);

  for (size_t i = 0; i < T.num_rows(); ++i) {
    const Row &m = T[i];
    for (int j = 0; j < lngt; ++j)
      if ((j > pred[m] && m[j] == true) || (j > succ[m] && m[j] == false))
        H.push_back(hext(m, j));
  }

  */
  H = primality(H, T);
  cook(H);
  return H;
}

// learn general CNF formula with large strategy
// from negative examples F
Formula learnCNFlarge(const Matrix &F) {
  Formula formula;
  for (size_t j = 0; j < F.num_rows(); ++j) {
    const Row &row = F[j];
    Clause clause;
    // for (bool bit : row)
    for (int i = 0; i < row.size(); ++i) {
      // isolate row[i]
      Literal lit;
      if (row[i] > 0) {
        lit.sign = lneg;
        lit.nval = row[i] - 1;
      }
      if (row[i] < DCARD) {
        lit.sign = Sign(lit.sign | lpos);
        lit.pval = row[i] + 1;
      }
      clause.push_back(lit);
    }
    formula.push_back(clause);
  }
  cook(formula);
  return formula;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// UNSAFE if m1 == m2
int fork(const Row &m1, const Row &m2) {
  int i = 0;
  // we can guarantee that always m1 != m2, therefore
  // we can drop i < m1.size()
  while (m1[i] == m2[i])
    ++i;
  return i;
}

deque<int> fork(const Matrix &M) {
  deque<int> frk(M.num_rows() - 1, 0);
  for (int ell = 0; ell < M.num_rows() - 1; ++ell)
    frk[ell] = fork(M[ell], M[ell + 1]);
  frk.push_front(SENTINEL);
  frk.push_back(SENTINEL);
  return frk;
}

/*
Clause negTerm(const Row &m, const int &i) {
  Clause clause(m.size(), Literal::none());
  for (int j = 0; j < i; ++j)
    clause[j] = m[j] == false ? lpos : lneg;
  return clause;
}

Clause negLeft(const Row &m, const int &i) {
  Clause clause = negTerm(m, i);
  clause[i] = lpos; // negLT
  return clause;
}

Clause negRight(const Row &m, const int &i) {
  Clause clause = negTerm(m, i);
  clause[i] = lneg; // negGT
  return clause;
}
*/

// TODO: What the frick
// learn general CNF formula with exact strategy
// from positive examples T
Formula learnCNFexact(const Matrix &T) {
  Formula formula;
  /*
  // sort(T.begin(), T.end());
  sort_matrix(T, 0, T.size() - 1);
  auto ip = unique(T.begin(), T.end());
  T.resize(distance(T.begin(), ip));

  deque<int> frk = fork(T); // frk.size() == T.size()+1
  for (int ell = 0; ell < T.size(); ++ell) {
    for (int i = frk[ell] + 1; i < arity; ++i)
      if (T[ell][i] == true)
        formula.push_back(negLeft(T[ell], i));
    for (int i = frk[ell + 1] + 1; i < arity; ++i)
      if (T[ell][i] == false)
        formula.push_back(negRight(T[ell], i));
  }
  formula = primality(formula, T);
  */
  cook(formula);
  return formula;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// swap the polarity of values in a tuple
void polswap_row(Row &row) {
  for (size_t i = 0; i < row.size(); ++i) {
    row[i] = DCARD - row[i];
  }
}

// swap polarity of every tuple in a matrix
void polswap_matrix(Matrix &A) {
  for (size_t i = 0; i < A.num_rows(); ++i) {
    polswap_row(A[i]);
  }
}

// swap polarity of literals in a clause
void polswap_clause(Clause &clause) {
  for (Literal &literal : clause)
    literal = literal.swap();
}

// swap polarity of every clause of a formula
void polswap_formula(Formula &formula) {
  for (Clause &clause : formula)
    polswap_clause(clause);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

string time2string(int seconds) {
  enum TimeUnit { second = 0, minute = 1, hour = 2, day = 3 };
  const string tu_name[] = {" second(s) ", " minute(s) ", " hour(s) ",
                            " day(s) "};
  const int timevalue[] = {60, 60, 24};

  int timeunit[4] = {0, 0, 0, 0};

  if (seconds == 0)
    return "0 seconds";

  for (int t = second; t < day; t++) {
    timeunit[t] = seconds % timevalue[t];
    seconds /= timevalue[t];
  }
  if (seconds > 0)
    timeunit[day] = seconds;

  string output;
  for (int t = day; t >= second; t--) {
    if (timeunit[t] > 0)
      output += to_string(timeunit[t]) + tu_name[t];
  }
  return output;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
