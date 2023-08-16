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
 *     File:    src/mcp-parallel.hpp                                      *
 *                                                                        *
 *      Copyright (c) 2019 - 2023                                         *
 *                                                                        *
 **************************************************************************/

#pragma once

#include "mcp-matrix+formula.hpp"
#include <string>

using namespace std;

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void adjust();
void print_arg();
Formula learnHornLarge(ofstream &process_outfile, const Matrix &T,
                       const Matrix &F);
Formula learnBijunctive(ofstream &process_outfile, const Matrix &T,
                        const Matrix &F);
// void OneToOne (ofstream &process_outfile, const int &i);
// void OneToAll (ofstream &process_outfile, const int &i);
// void OneToAllNosection (ofstream &process_outfile, const int &i);
void split_action(ofstream &popr, ofstream &latpr, const int &process_id);
void crash();

//==================================================================================================
