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
 *     File:    src/mcp-mpi.cpp                                           *
 *                                                                        *
 *      Copyright (c) 2019 - 2023                                         *
 *                                                                        *
 **************************************************************************/

#include <fstream>
#include <iostream>
#include <memory>
#include <mpi.h>
// #include <mutex>
#include "mcp-common.hpp"
#include "mcp-matrix+formula.hpp"
#include "mcp-parallel.hpp"

using namespace std;

Arch arch = archMPI;

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void _row_inplace_minimum(Row &MIN, const RowView &other) {
  for (size_t i = 0; i < other.size(); ++i) {
    if (MIN[i] > other[i]) {
      MIN[i] = other[i];
    }
  }
}

bool InHornClosure(const RowView &row, const MatrixMask &M) {
  // selects rows of M above the input row, computes their minimum
  // then checks wether the row is equal to the minimum
  Row MIN;
  bool any_above = false;
  for (size_t i = 0; i < M.num_rows(); ++i) {
    if (row >= M[i]) {
      if (any_above) {
        _row_inplace_minimum(MIN, M[i]);
      } else {
        any_above = true;
        MIN = Row(M[i].size());
        for (size_t j = 0; j < M[i].size(); ++j) {
          MIN[j] = M[i][j];
        }
      }
    }
  }
  return any_above && row == MIN;
}
////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
  int num_procs;    // for MPI
  int process_rank; // for MPI

  time_t start_time = time(nullptr);

  version += arch_strg[arch];
  set_terminate(crash);

  read_arg(argc, argv);
  adjust();
  print_arg();
  read_matrix(group_of_matrix);
  print_matrix(group_of_matrix);

  if (action == aSELECTED) {
    cerr << endl;
    cerr << "*** no need to run the selected action in parallel" << endl;
    cerr << "*** use mcp-seq to run this example" << endl;
    exit(2);
  }

  const string temp_prefix = tpath + "mcp-tmp-";
  const string basename = temp_prefix + to_string(start_time);
  // ofstream *process_outfile = new ofstream[grps.size()];	// replaced with
  // smart pointer ofstream *latex_outfile = new ofstream[grps.size()];	//
  // replaced with smart pointer
  auto process_outfile = make_unique<ofstream[]>(grps.size());
  auto latex_outfile = make_unique<ofstream[]>(grps.size());

  int ierr = MPI_Init(&argc, &argv);
  if (ierr != 0) {
    cerr << endl << "*** MCP - fatal error" << endl;
    cerr << "*** MPI_Init returned error = " << ierr << endl;
    exit(1);
  }

  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);

  if (process_rank == 0) {
    outfile << "+++ Number of processes = " << num_procs << endl << endl;
    if (np_fit && grps.size() < num_procs) {
      cerr << "+++ " << num_procs - grps.size() << " unused processes" << endl;
      outfile << "+++ " << num_procs - grps.size() << " unused processes"
              << endl;
      cerr << "+++ launch 'mpirun -np " << grps.size() << " " << argv[0]
           << " --fit yes ...'" << endl;
      outfile << "+++ launch 'mpirun -np " << grps.size() << " " << argv[0]
              << " --fit yes ...'" << endl;
      MPI_Abort(MPI_COMM_WORLD, 99);
    }
  }

  string filename = basename + "-" + to_string(process_rank) + ".txt";
  string latexname = basename + "-" + to_string(process_rank) + ".tex";

  if (process_rank < grps.size()) {
    // semaphore.lock();
    process_outfile[process_rank].open(filename);
    // semaphore.unlock();
    if (!process_outfile[process_rank].is_open()) {
      cerr << "+++ Cannot open output file " << filename << endl;
      exit(2);
    }
    process_outfile[process_rank] << "+++ Start output of process "
                                  << process_rank << endl
                                  << endl;
    if (latex.length() > 0) {
      latex_outfile[process_rank].open(latexname);
      if (!latex_outfile[process_rank].is_open()) {
        cerr << "+++ Cannot open latex file " << latexname << endl;
        exit(2);
      }
      latex_outfile[process_rank] << "% Start output of process "
                                  << process_rank << endl
                                  << endl;
    }

    for (int process_id = process_rank; process_id < grps.size();
         process_id += num_procs)
      split_action(process_outfile[process_rank], latex_outfile[process_rank],
                   process_id);

    process_outfile[process_rank] << "+++ End output of process "
                                  << process_rank << endl;
    // semaphore.lock();
    process_outfile[process_rank].close();
    // semaphore.unlock();
    if (latex.length() > 0) {
      latex_outfile[process_rank] << "% End output of process " << process_rank
                                  << endl;
      latex_outfile[process_rank].close();
    }
  }

  // MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  if (process_rank == 0) {
    // delete [] process_outfile;	// useless, since we use smart pointers
    // delete [] latex_outfile;		// useless, since we use smart pointers

    int gsize = grps.size();
    for (int proc_id = 0; proc_id < min(num_procs, gsize); ++proc_id) {
      string filename = basename + "-" + to_string(proc_id) + ".txt";

      ifstream in(filename);
      string line;

      while (getline(in, line))
        outfile << line << endl;
      outfile << endl;
      // semaphore.lock();
      in.close();
      remove(filename.c_str());
      // semaphore.unlock();

      if (latex.length() > 0) {
        string latexname = basename + "-" + to_string(proc_id) + ".tex";
        ifstream in(latexname);
        while (getline(in, line))
          latexfile << line << endl;
        latexfile << endl;
        in.close();
        remove(latexname.c_str());
      }
    }

    time_t finish_time = time(nullptr);
    outfile << "+++ time = " << time2string(difftime(finish_time, start_time))
            << endl;

    outfile << "+++ end of run +++" << endl;
    outfile.close();
    if (latex.length() > 0)
      latexfile.close();
  }
  return 0;
}
////////////////////////////////////////////////////////////////////////////////
