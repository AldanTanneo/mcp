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
 *	Version: common for all                                           *
 *      File:    mcp-common.hpp                                           *
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

#pragma once

#include "mcp-matrix+formula.hpp"
#include <deque>
#include <map>
#include <string>
#include <vector>

// #define GLOBAL_VERSION "1.04-c++-"

//--------------------------------------------------------------------------------------------------

extern std::string version;
extern bool debug;

// extern string varid;
// extern bool varswitch;
// extern vector<string> varnames;
// enum NAME {nOWN = 0, nPOSITIVE = 1, nNEGATIVE = 2};

extern std::map<Row, int>
    pred; // predecessor function for Zanuttini's algorithm
extern std::map<Row, int> succ; // successor function for Zanuttini's algorithm
extern std::map<Row, std::vector<int>>
    sim; // sim table for Zanuttini's algorithm

// enum Action    {aONE    = 0, aALL    = 1, aNOSECT      = 2};
enum Closure : char {
  clHORN = 0,
  clDHORN = 1,
  clBIJUNCTIVE = 2,
  clAFFINE = 3,
  clCNF = 4
};
enum Cooking : char { ckRAW = 0, ckBLEU = 1, ckMEDIUM = 2, ckWELLDONE = 3 };
enum Direction : char {
  dBEGIN = 0,
  dEND = 1,
  dOPT = 2,
  dRAND = 3,
  dLOWCARD = 4,
  dHIGHCARD = 5
};
// enum Print     {pVOID   = 0, pCLAUSE = 1, pIMPL        = 2, pMIX       = 3,
// pDIMACS   = 4};
enum Strategy : char { sLARGE = 0, sEXACT = 1 };
// enum Display   {yUNDEF  = 0, yHIDE   = 1, yPEEK        = 2, ySECTION   = 3,
// ySHOW     = 4};
enum Arch : char { archSEQ = 0, archMPI = 1, archPTHREAD = 2, archHYBRID = 3 };

// extern const int SENTINEL;
extern const std::string STDIN;
extern const std::string STDOUT;
// extern const int MTXLIMIT;

// extern Action action;
extern Closure closure;
extern Cooking cooking;
extern Direction direction;
// extern Print print;
extern bool setcover;
extern Strategy strategy;
// extern Display display;
extern std::string input;
extern std::string output;
extern bool disjoint;
// extern int arity;

// extern int offset;
extern std::string tpath; // directory where the temporary files will be stored
extern bool np_fit;
extern int chunkLIMIT; // heavily hardware dependent; must be optimized
extern Arch arch;
extern std::string latex; // file to store latex output

extern std::ifstream infile;
extern std::ofstream outfile;
extern std::ofstream latexfile;
extern std::string formula_output;

extern const std::string action_strg[];
extern const std::string closure_strg[];
extern const std::string cooking_strg[];
extern const std::string direction_strg[];
extern const std::string pcl_strg[];
extern const std::string strategy_strg[];
// extern const string print_strg[];
// extern const string display_strg[];
extern const std::string arch_strg[];

//--------------------------------------------------------------------------------------------------

void read_arg(int argc, char *argv[]);
int hamming_distance(const Row &u, const Row &v);

// ostream& operator<< (ostream &output, const Row &row);
// ostream& operator<< (ostream &output, const Matrix &M);

bool InHornClosure(const RowView &a, const MatrixMask &M);
bool InHornClosure(const Row &a, const Matrix &M);

bool inadmissible(const Matrix &T, const Matrix &F);

int hamming_weight(const Row &row);

std::vector<bool> minsect(const Matrix &T, const Matrix &F);

bool satisfied_by(const Clause &clause, const Matrix &T);

Matrix restrict(const std::vector<bool> &sect, const Matrix &A);

Matrix HornClosure(const Matrix &M);

Row minHorn(const Matrix &M);

// void sort_formula (Formula &formula, int low, int high);
// Formula unitres (const Formula &formula);
// Formula binres (const Formula &formula);
// Formula subsumption (Formula formula);
// Formula redundant (Formula formula);
Formula primality(const Formula &phi, const Matrix &M);
Formula SetCover(const Matrix &Universe, const Formula &SubSets);
void cook(Formula &formula);
Formula learnHornExact(const Matrix &T);
Formula learnCNFlarge(const Matrix &F);
Formula learnCNFexact(const Matrix &T);
void write_formula(const std::string &suffix1, const std::string &suffix2,
                   const std::vector<int> &names, const Formula &formula);
void write_formula(const std::string &suffix, const std::vector<int> &names,
                   const Formula &formula);

void polswap_formula(Formula &formula);
void polswap_matrix(Matrix &);

std::string time2string(int seconds);

//==================================================================================================
