/**************************************************************************
 *                                                                        *
 *                                                                        *
 *	       Multiple Characterization Problem (MCP)                        *
 *                                                                        *
 *	Author:   Miki Hermann                                                *
 *	e-mail:   hermann@lix.polytechnique.fr                                *
 *	Address:  LIX (CNRS UMR 7161), Ecole Polytechnique, France            *
 *                                                                        *
 *	Author:   Gernot Salzer                                               *
 *	e-mail:   gernot.salzer@tuwien.ac.at                                  *
 *	Address:  Technische Universitaet Wien, Vienna, Austria               *
 *                                                                        *
 * Author:   César Sagaert                                                *
 * e-mail:   cesar.sagaert@ensta-paris.fr                                 *
 * Address:  ENSTA Paris, Palaiseau, France                               *
 *                                                                        *
 *	Version: all                                                          *
 *     File:    src/mcp-seq.cpp                                           *
 *                                                                        *
 *      Copyright (c) 2019 - 2023                                         *
 *                                                                        *
 **************************************************************************/

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <queue>
#include <sstream>

#include "mcp-bucket.hpp"
#include "mcp-common.hpp"
#include "mcp-matrix+formula.hpp"
#include "mcp-mesh.hpp"

using namespace std;

Arch arch = archSEQ;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void adjust() { // adjusts the input parameters
  if (input != STDIN) {
    infile.open(input);
    if (infile.is_open())
      cin.rdbuf(infile.rdbuf());
    else {
      cerr << "+++ Cannot open input file " << input << endl;
      exit(2);
    }
  }

  if (output != STDOUT) {
    outfile.open(output);
    if (outfile.is_open())
      cout.rdbuf(outfile.rdbuf());
    else {
      cerr << "+++ Cannot open output file " << output << endl;
      exit(2);
    }
  }

  if (latex.length() > 0) {
    if (latex.find(".") == string::npos || latex.length() < 4 ||
        latex.substr(latex.length() - 4) != ".tex")
      latex += ".tex";
    latexfile.open(latex);
    if (!latexfile.is_open()) {
      cout << "+++ Cannot open latex file " << latex << endl;
      exit(2);
    }
  }

  if (offset < 0 && print != pDIMACS) {
    cout << "+++ WARNING: offset reset to 0" << endl;
    offset = 0;
  } else if (offset <= 0 && print == pDIMACS) {
    cout << "+++ WARNING: offset reset to 1" << endl;
    offset = 1;
  }

  if (print == pVOID)
    print = (closure <= clDHORN) ? pMIX : pCLAUSE;
  if (strategy == sEXACT && setcover) {
    cout << "+++ WARNING: set cover abandoned" << endl;
    setcover = false;
  }
}

void print_arg() {
  cout << "@@@ Parameters:" << endl;
  cout << "@@@ ===========" << endl;
  cout << "@@@ version       = " << version << endl;
  cout << "@@@ input         = " << input << endl;
  cout << "@@@ output        = " << output << endl;
  cout << "@@@ latex output  = " << (latex.length() > 0 ? latex : "no") << endl;
  cout << "@@@ action        = " << action_strg[action] << endl;
  cout << "@@@ closure       = " << closure_strg[closure] << endl;
  cout << "@@@ direction     = " << direction_strg[direction] << endl;
  if (direction == dRAND)
    cout << "@@@ random seed   = " << random_seed << endl;
  cout << "@@@ strategy      = " << strategy_strg[strategy] << endl;
  cout << "@@@ cooking       = " << cooking_strg[cooking] << endl;
  cout << "@@@ set cover     = " << (setcover ? "yes" : "no") << endl;
  cout << "@@@ var. offset   = " << offset << endl;
  cout << "@@@ print matrix  = " << display_strg[display]
       << (display == yUNDEF ? " (will be changed)" : "") << endl;
  cout << "@@@ print formula = " << print_strg[print] << endl;
  cout << "@@@ out   formula = "
       << (formula_output.empty() ? "none" : formula_output) << endl;
  cout << "@@@ debug         = " << (debug ? "yes" : "no") << endl << endl;
}

void read_matrix(Group_of_Matrix &matrix) {
  // reads the input matrices
  int ind_a, ind_b;
  string line;

  getline(cin, line);
  istringstream inds(line);
  inds >> ind_a >> ind_b;
  cout << "+++ Indication line: " << ind_a << " " << ind_b << endl;

  if (ind_a == 1) {
    cout << "+++ Own names for variables" << endl;
    varswitch = true;

    getline(cin, line);
    istringstream vars(line);
    string vname;
    while (vars >> vname)
      varnames.push_back(vname);
    arity = varnames.size();
  }
  if (ind_b == 1)
    getline(cin, line);

  string group;
  size_t numline = 0;
  while (getline(cin, line)) {
    numline++;
    /*
    istringstream nums(line);
    nums >> group;
    Row temp;
    integer number;
    while (nums >> number)
      temp.push_back(number);
    */
    const auto nums = split(line, ",");
    group = nums.at(0);
    Row temp;
    for (size_t i = 1; i < nums.size(); ++i) {
      integer x = integer(stoull(nums.at(i)));
      DMAX = max(x, DMAX);
      temp.push_back(x);
    }
    if (arity == 0)
      arity = temp.size();
    else if (arity != temp.size())
      cout << "*** arity discrepancy on line " << numline << endl;

    if (matrix.find(group) == matrix.end()) {
      matrix.insert({group, Matrix()});
    }
    matrix[group].add_row(std::move(temp));
  }

  if (input != STDIN)
    infile.close();

  if (display == yUNDEF) {
    display = (numline * arity > MTXLIMIT) ? yHIDE : yPEEK;
    cout << "@@@ print matrix  = " << display_strg[display] << " (redefined)"
         << endl;
  }
}

void print_matrix(const Group_of_Matrix &matrix) {
  // prints the matrices
  cout << "+++ Arity = " << arity << endl;
  for (const auto &group : matrix) {
    cout << "+++ Group " << group.first;
    grps.push_back(group.first);
    const Matrix &gmtx = group.second;
    cout << " [" << gmtx.num_rows() << "]:" << endl;
    if (display == yPEEK || display == ySHOW)
      cout << gmtx << endl;
  }
  sort(grps.begin(), grps.end());
  cout << "+++ Number of groups = " << grps.size() << endl;

  if (action == aSELECTED) {
    if (find(grps.begin(), grps.end(), selected) == grps.end()) {
      cout << "+++ Selected group " << selected << " does not exist" << endl;
      exit(2);
    } else
      cout << "@@@ selected group   = " << selected << endl;
  }

  cout << endl;
}

// selects rows of matrix above the input row, computes their minimum
// then checks wether the row is equal to the minimum
template <typename R, typename M>
bool InHornClosure(const R &row, const M &matrix) {
  if (matrix.empty()) {
    return false;
  }
  Row MIN;
  bool any_above = false;
  for (size_t i = 0; i < matrix.num_rows(); ++i) {
    const R &r = matrix[i];
    if (r >= row) {
      if (any_above)
        MIN.inplace_minimum(r);
      else {
        MIN = r.to_row();
        any_above = true;
      }
    }
  }
  return any_above && row == MIN;
}

bool InHornClosure(const RowView &r, const MatrixMask &m) {
  return InHornClosure<RowView, MatrixMask>(r, m);
}
bool InHornClosure(const Row &r, const Matrix &m) {
  return InHornClosure<Row, Matrix>(r, m);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// learn a Horn clause from positive examples T and negative examples F
// with the large strategy
Formula learnHornLarge(const Matrix &positiveT, const Matrix &negativeF) {
  Formula varphi;
  size_t arity = positiveT.num_cols();

  for (size_t i = 0; i < negativeF.num_rows(); ++i) {
    const Row &f = negativeF[i];
    if (!sat_formula(f, varphi)) {
      continue;
    } else if (InHornClosure(f, positiveT)) {
      cerr << "+++ negative example present in Horn closure of T" << endl;
      cerr << "+++ the negative culprit is '" << f << "'" << endl;
      exit(2);
    }

    Clause c;
    c.reserve(arity);
    for (size_t i = 0; i < arity; ++i) {
      Sign sign = f[i] > 0 ? lneg : lnone;
      integer val = f[i] > 0 ? f[i] - 1 : 0;
      Literal lit(sign, 0, val);
      c.push_back(lit);
    }
    if (sat_clause(positiveT, c)) {
      varphi.push_back(std::move(c));
    } else {
      bool eliminated = false;
      size_t i = 0;
      Literal old;
      while (!eliminated && i < arity) {
        if (f[i] < DMAX) {
          old = c[i];
          c[i].sign = (Sign)(c[i].sign | lpos);
          c[i].pval = f[i] + 1;
          if (sat_clause(positiveT, c)) {
            varphi.push_back(std::move(c));
            eliminated = true;
          } else {
            c[i] = old;
          }
        }
        i++;
      }
      if (!eliminated) {
        cerr << "+++ FAILURE: could not eliminate..." << endl;
      }
    }
  }

  if (debug)
    cerr << "*** before cook in learnHornLarge" << endl;
  cook(varphi);
  if (debug)
    cerr << "*** after cook" << endl;
  return varphi;
}

// learn a bijunctive clause from positive examples T and negative examples F
Formula learnBijunctive(const Matrix &positiveT, const Matrix &negativeF) {
  bucket::Bucket bucket;
  mesh::Strip strip;
  mesh::Mesh mesh;

  size_t arity = positiveT.num_cols();
  mesh::init(mesh, arity);
  mesh::init(strip, arity);

  for (size_t k = 0; k < negativeF.num_rows(); ++k) {
    const Row &f = negativeF[k];

    if (!(bucket::sat_bucket(f, bucket))) {
      continue;
    }

    bucket::Clause c;
    bool eliminated = false;
    size_t i = 0;

    while (!eliminated && i < arity) {
      c = mesh::isolation({f[i], f[i]}, {i, i}, mesh::SE);
      // Test if f[i] is in strip[i] and valid
      if (strip[i].count(f[i]) == 0 && bucket::valid(c))
        eliminated = true;

      size_t j = i + 1;
      while (!eliminated && j < arity) {
        const bucket::Point f2 = {f[i], f[j]};
        const std::array<size_t, 2> ij = {i, j};
        // Test if f2 can be isolated by trying all directions

        mesh::Direction q = mesh::NW;
        while (!eliminated) {
          if (mesh[i][j].isolated(f2, q)) {
            c = mesh::isolation(f2, ij, q);
            eliminated = bucket::valid(c);
          }

          if (q == mesh::SW)
            break;
          // undefined behavior when q == mesh::SW
          q = mesh::Direction(q + 1);
        }

        j++;
      }

      if (eliminated)
        bucket::insert(c, bucket);
      i++;
    }
    if (!eliminated) {
      cerr << "+++ 2SAT formula not possible:    No isolated point" << endl;
      exit(2);
    }
  }

  Formula B = get_formula(bucket, arity);
  cook(B);
  return B;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Formula post_prod(const vector<size_t> &A, const Matrix &F,
                  const Formula &formula) {
  Formula schf;
  if (setcover == true) {
    cout << "+++ " << pcl_strg[closure] << " formula before set cover ["
         << formula.size() << "] =" << endl;
    cout << formula2string(A, formula) << endl;
    if (latex.length() > 0) {
      latexfile << "% " << pcl_strg[closure] << " formula before set cover ["
                << formula.size() << "] =" << endl;
      latexfile << "\\varphi = " << endl;
      latexfile << formula2latex(A, formula) << endl << endl;
    }

    schf = SetCover(F, formula);
    cout << "+++ " << pcl_strg[closure] << " formula after set cover ["
         << schf.size() << "] =" << endl;
    cout << formula2string(A, schf) << endl;
    if (latex.length() > 0) {
      latexfile << "% " << pcl_strg[closure] << " formula after set cover ["
                << schf.size() << "] =" << endl;
      latexfile << "\\varphi = " << endl;
      latexfile << formula2latex(A, schf) << endl << endl;
    }
  } else {
    cout << "+++ " << pcl_strg[closure] << " formula [" << formula.size()
         << "] =" << endl;
    cout << formula2string(A, formula) << endl;
    if (latex.length() > 0) {
      latexfile << "% " << pcl_strg[closure] << " formula [" << formula.size()
                << "] =" << endl;
      latexfile << "\\varphi = " << endl;
      latexfile << formula2latex(A, formula) << endl << endl;
    }
    schf = formula;
  }

  if (closure == clDHORN) {
    cout << "+++ swapping the formula back to dual Horn" << endl;
    // Formula dschf = polswap_formula(schf);
    polswap_formula(schf);
    cout << "+++ final dual Horn formula [" << schf.size() << "] =" << endl;
    cout << formula2string(A, schf) << endl;
    if (latex.length() > 0) {
      latexfile << "% swapping the formula back to dual Horn" << endl;
      latexfile << "% final dual Horn formula [" << schf.size()
                << "] =" << endl;
      latexfile << "\\varphi = " << endl;
      latexfile << formula2latex(A, schf) << endl << endl;
    }
  }
  return schf;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// one group as positive agains one group of negative examples
void OneToOne() {

  for (size_t i = 0; i < grps.size(); ++i) {
    Matrix T = group_of_matrix[grps[i]].clone();
    for (size_t j = 0; j < grps.size(); ++j) {
      if (j == i)
        continue;
      Matrix F = group_of_matrix[grps[j]].clone();

      if (closure == clDHORN) {
        cout << "+++ swapping polarity of vectors and treating swapped vectors "
                "as Horn"
             << endl;
        polswap_matrix(T);
        polswap_matrix(F);
      }

      Mask sect = minsect(T, F);
      cout << "+++ Section of groups T=" << grps[i] << " and F=" << grps[j]
           << ":" << endl;
      if (!disjoint) {
        cout << "+++ Matrices <T> and F are not disjoint, therefore I cannot "
                "infer a formula"
             << endl
             << endl;
      } else {
        int hw = std::accumulate(sect.cbegin(), sect.cend(), 0);
        cout << "+++ Relevant variables [" << hw << "]: ";
        vector<size_t> A;
        for (size_t k = 0; k < sect.size(); ++k)
          if (sect[k] == true) {
            A.push_back(k);
            if (varswitch) {
              vector<string> new_names = split(varnames[k], ":");
              cout << new_names[nOWN];
            } else
              cout << varid + to_string(offset + k);
            cout << " ";
          }
        cout << endl;
        cout << "+++ A [" << hw << "] = { ";
        for (size_t coord : A)
          cout << offset + coord << " ";
        cout << "}" << endl;
        T.restrict(sect);
        F.restrict(sect);

        cout << "+++ T|_A [" << T.num_rows() << "]";
        if (display >= ySECTION) {
          cout << " = { " << endl;
          cout << T;
          cout << "+++ }";
        }
        cout << endl;

        cout << "+++ F|_A [" << F.num_rows() << "]";
        if (display >= ySECTION) {
          cout << " = { " << endl;
          cout << F;
          cout << "+++ }";
        }
        cout << endl;

        // Matrix HC;
        // if (closure == clHORN || closure == clDHORN) {
        //   HC = HornClosure(TsectA);
        //   cout << "+++ " << pcl_strg[closure] << "Closure(T|_A) [" <<
        //   HC.size() << "]"; if (display >= ySECTION) {
        //     cout << " = { " << endl;
        //     cout << HC;
        //     cout << "+++ }";
        //   }
        //   cout << endl;
        // } else if (closure == clBIJUNCTIVE)
        //   cout << "+++ " << pcl_strg[closure] << "Closure not computed" <<
        //   endl;

        Formula formula;
        if (closure == clHORN || closure == clDHORN)
          formula =
              strategy == sEXACT ? learnHornExact(T) : learnHornLarge(T, F);
        else if (closure == clBIJUNCTIVE) {
          formula = learnBijunctive(T, F);
          if (formula.empty()) {
            cout << "+++ bijunctive formula not possible for this configuration"
                 << endl
                 << endl;
            continue;
          }
        } else if (closure == clCNF)
          formula = strategy == sLARGE ? learnCNFlarge(F) : learnCNFexact(T);

        Formula schf = post_prod(A, F, formula);
        if (!formula_output.empty())
          write_formula(grps[i], grps[j], A, schf);
      }
      disjoint = true;
      cout << endl;
    }
  }
}

// selected group of positive examples against all other groups together as
// negative examples
void SelectedToAll(const string &grp) {
  Matrix T = group_of_matrix.at(grp).clone();
  Matrix F;
  vector<string> index;
  for (size_t j = 0; j < grps.size(); ++j) {
    if (grps[j] == grp)
      continue;

    F.reserve(F.num_rows() + group_of_matrix[grps[j]].num_rows());
    for (size_t i = 0; i < group_of_matrix[grps[j]].num_rows(); ++i) {
      F.add_row(group_of_matrix[grps[j]][i].clone());
    }

    index.push_back(grps[j]);
  }
  sort(index.begin(), index.end());

  if (closure == clDHORN) {
    cout << "+++ swapping polarity of vectors and treating swapped vectors as "
            "Horn"
         << endl;
    polswap_matrix(T);
    polswap_matrix(F);
  }

  Mask sect = minsect(T, F);
  cout << "+++ Section of groups T=" << grp << " and F=( ";
  for (string coord : index)
    cout << coord << " ";
  cout << "):" << endl;
  if (!disjoint)
    cout << "+++ Matrices <T> and F are not disjoint, therefore I cannot infer "
            "a formula"
         << endl
         << endl;
  else {
    int hw = std::accumulate(sect.cbegin(), sect.cend(), 0);
    cout << "+++ Relevant variables [" << hw << "]: ";
    vector<size_t> A;
    for (size_t k = 0; k < sect.size(); ++k)
      if (sect[k]) {
        A.push_back(k);
        if (varswitch) {
          vector<string> new_names = split(varnames[k], ":");
          cout << new_names[nOWN];
        } else
          cout << varid << to_string(offset + k);
        cout << " ";
      }
    cout << endl;
    cout << "+++ A [" << hw << "] = { ";
    for (size_t var : A)
      cout << offset + var << " ";
    cout << "}" << endl;
    T.restrict(sect);
    F.restrict(sect);

    cout << "+++ T|_A [" << T.num_rows() << "," << T.num_cols() << "]";
    if (display >= ySECTION) {
      cout << " = { " << endl;
      cout << T;
      cout << "+++ }";
    }
    cout << endl;

    cout << "+++ F|_A [" << F.num_rows() << "," << F.num_cols() << "]";
    if (display >= ySECTION) {
      cout << " = { " << endl;
      cout << F;
      cout << "+++ }";
    }
    cout << endl;

    Formula formula;
    if (closure == clHORN || closure == clDHORN)
      formula = strategy == sEXACT ? learnHornExact(T) : learnHornLarge(T, F);
    else if (closure == clBIJUNCTIVE) {
      formula = learnBijunctive(T, F);
      if (formula.empty()) {
        cout << "+++ bijunctive formula not possible for this configuration"
             << endl
             << endl;
        return;
      }
    } else if (closure == clCNF)
      formula = strategy == sLARGE ? learnCNFlarge(F) : learnCNFexact(T);

    Formula schf = post_prod(A, F, formula);
    if (!formula_output.empty())
      write_formula(grp, A, schf);
  }
  disjoint = true;
  cout << endl;
}

// one group of positive examples against all other groups together as negative
// examples
void OneToAll() {
  for (auto &grp : grps)
    SelectedToAll(grp);
}

// one group of positive examples against all other groups together as negative
// examples no section is done
void OneToAllNosection() {

  vector<size_t> names(arity);
  for (size_t nms = 0; nms < arity; ++nms)
    names[nms] = nms;

  for (size_t i = 0; i < grps.size(); ++i) {
    cout << "group " << i << ": " << grps[i] << endl;
    Matrix T = group_of_matrix[grps[i]].clone();
    Matrix F;
    vector<string> index;
    for (size_t j = 0; j < grps.size(); ++j) {
      if (j == i)
        continue;

      F.reserve(F.num_rows() + group_of_matrix[grps[j]].num_rows());
      for (size_t k = 0; k < group_of_matrix[grps[j]].num_rows(); ++k) {
        F.add_row(group_of_matrix[grps[j]][k].clone());
      }

      index.push_back(grps[j]);
    }
    sort(index.begin(), index.end());

    if (closure == clDHORN) {
      cout << "+++ swapping polarity of vectors and treating swapped vectors "
              "as Horn"
           << endl;
      polswap_matrix(T);
      polswap_matrix(F);
    }

    cout << "+++ Groups T=" << grps[i] << " and F=( ";
    for (string coord : index)
      cout << coord << " ";
    cout << "):" << endl;

    if (display == ySHOW) {
      cout << "+++ T [" << T.num_rows() << "] = { " << endl;
      cout << T;
      cout << "+++ }" << endl;

      cout << "+++ F [" << F.num_rows() << "] = { " << endl;
      cout << F;
      cout << "+++ }" << endl;
    } else {
      cout << "+++ T [" << T.num_rows() << "]" << endl;
      cout << "+++ F [" << F.num_rows() << "]" << endl;
    }

    // Matrix HC;
    // if (closure == clHORN || closure == clDHORN) {
    //   HC = HornClosure(T);
    //   cout << "+++ " << pcl_strg[closure] << "Closure(T) [" << HC.size() <<
    //   "]"; if (display == ySHOW) {
    // 	cout << " = { " << endl;
    // 	cout << HC;
    // 	cout << "+++ }";
    //   }
    //   cout << endl;
    // } else if (closure == clBIJUNCTIVE)
    //   cout << "+++ " << pcl_strg[closure] << "Closure not computed" << endl;

    if (inadmissible(T, F)) {
      cout << "*** F not disjoint from <T>" << endl;
      disjoint = false;
      cout << "+++ Matrices <T> and F are not disjoint, therefore I cannot "
              "infer a formula"
           << endl
           << endl;
    }

    if (disjoint) {
      Formula formula;
      if (closure == clHORN || closure == clDHORN)
        formula = strategy == sEXACT ? learnHornExact(T) : learnHornLarge(T, F);
      else if (closure == clBIJUNCTIVE) {
        formula = learnBijunctive(T, F);
        if (formula.empty()) {
          cout << "+++ bijunctive formula not possible for this configuration"
               << endl
               << endl;
          continue;
        }
      } else if (closure == clCNF)
        formula = strategy == sLARGE ? learnCNFlarge(F) : learnCNFexact(T);

      Formula schf = post_prod(names, F, formula);
      if (!formula_output.empty())
        write_formula(grps[i], names, schf);
    }
    disjoint = true;
    cout << endl;
  }
}

//////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
  time_t start_time = time(nullptr);

  version += arch_strg[arch];

  read_arg(argc, argv);
  adjust();
  print_arg();
  read_matrix(group_of_matrix);
  print_matrix(group_of_matrix);

  switch (action) {
  case aONE:
    OneToOne();
    break;
  case aALL:
    OneToAll();
    break;
  case aNOSECT:
    OneToAllNosection();
    break;
  case aSELECTED:
    SelectedToAll(selected);
    break;
  }

  time_t finish_time = time(nullptr);
  cout << "+++ time = " << time2string(int(difftime(finish_time, start_time)))
       << endl;

  cout << "+++ end of run +++" << endl;
  return 0;
}
////////////////////////////////////////////////////////////////////////////////
