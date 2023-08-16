##**************************************************************************
##*                                                                        *
##*                                                                        *
##*	       Multiple Characterization Problem (MCP)                         *
##*                                                                        *
##*	Author:   Miki Hermann                                                 *
##*	e-mail:   hermann@lix.polytechnique.fr                                 *
##*	Address:  LIX (CNRS UMR 7161), Ecole Polytechnique, France             *
##*                                                                        *
##*	Author:   Gernot Salzer                                                *
##*	e-mail:   gernot.salzer@tuwien.ac.at                                   *
##*	Address:  Technische Universitaet Wien, Vienna, Austria                *
##*                                                                        *
##* Author:   César Sagaert                                                *
##* e-mail:   cesar.sagaert@ensta-paris.fr                                 *
##* Address:  ENSTA Paris, Palaiseau, France                               *
##*                                                                        *
##*	Version: all                                                           *
##*     File:    src/Makefile                                              *
##*                                                                        *
##*      Copyright (c) 2019 - 2023                                         *
##*                                                                        *
##**************************************************************************

BIN   := ../bin
FLAGS := -std=c++14 -O4 
# -Wall -Werror -Wfatal-errors -Wconversion -Os -fsanitize=address -fsanitize=undefined
CXX   := g++ $(FLAGS)

.PHONY: compile compile-no-mpi seq export clean scratch uniq

compile: seq mpi pthread hybrid trans guess cnf split check sparse \
	predict chk2tst uniq

compile-no-mpi: seq pthread trans guess cnf split check sparse predict \
	chk2tst uniq


#---------------------------------------------------------------------------------------------------

mcp-matrix+formula-seq.o: mcp-matrix+formula.cpp mcp-matrix+formula.hpp
	$(CXX) -c -o $@ mcp-matrix+formula.cpp

mcp-common-seq.o: mcp-common.cpp mcp-common.hpp mcp-matrix+formula.hpp
	$(CXX) -c -o $@ mcp-common.cpp

mcp-bucket-seq.o: mcp-bucket.cpp mcp-bucket.hpp mcp-matrix+formula.hpp
	$(CXX) -c -o $@ mcp-bucket.cpp

mcp-mesh-seq.o: mcp-mesh.cpp mcp-mesh.hpp mcp-bucket.hpp mcp-matrix+formula.hpp
	$(CXX) -c -o $@ mcp-mesh.cpp

#---------------------------------------------------------------------------------------------------

seq: mcp-matrix+formula-seq.o mcp-common-seq.o mcp-seq.o mcp-mesh-seq.o mcp-bucket-seq.o
	$(CXX) -o $(BIN)/mcp-seq mcp-seq.o mcp-common-seq.o mcp-matrix+formula-seq.o mcp-bucket-seq.o mcp-mesh-seq.o


mcp-seq.o: mcp-seq.cpp mcp-bucket.hpp mcp-matrix+formula.hpp mcp-common.hpp mcp-mesh.hpp
	$(CXX) -c -o mcp-seq.o mcp-seq.cpp

#---------------------------------------------------------------------------------------------------

mpi: matrix+formula-mpi.o common-mpi.o parallel-mpi.o mpi.o
	mpic++ -O4 -o $(BIN)/mcp-mpi mcp-mpi.o \
				     mcp-common-mpi.o \
				     mcp-matrix+formula-mpi.o \
				     mcp-parallel-mpi.o

mpi.o: mcp-mpi.cpp mcp-common.hpp mcp-parallel.hpp
	mpic++ -O4 -c -o mcp-mpi.o mcp-mpi.cpp

parallel-mpi.o: mcp-parallel.cpp mcp-parallel.hpp mcp-common.hpp
	mpic++ -O4 -c -o mcp-parallel-mpi.o mcp-parallel.cpp

common-mpi.o: mcp-common.cpp mcp-common.hpp
	mpic++ -O4 -c -o mcp-common-mpi.o mcp-common.cpp

matrix+formula-mpi.o: mcp-matrix+formula.cpp mcp-matrix+formula.hpp
	mpic++ -O4 -c -o mcp-matrix+formula-mpi.o mcp-matrix+formula.cpp

#---------------------------------------------------------------------------------------------------

pthread: matrix+formula-pthread.o common-pthread.o parallel-pthread.o posix-pthread.o pthread.o
	$(CXX) -pthread -O4 -o $(BIN)/mcp-pthread mcp-pthread.o \
					       mcp-common-pthread.o \
					       mcp-matrix+formula-pthread.o \
					       mcp-parallel-pthread.o \
					       mcp-posix-pthread.o

pthread.o: mcp-pthread.cpp mcp-common.hpp mcp-parallel.hpp mcp-posix.hpp
	$(CXX) -pthread -c -O4 -o mcp-pthread.o mcp-pthread.cpp

posix-pthread.o: mcp-posix.cpp mcp-posix.hpp mcp-parallel.hpp mcp-common.hpp
	$(CXX) -pthread -c -O4 -o mcp-posix-pthread.o mcp-posix.cpp

parallel-pthread.o: mcp-parallel.cpp mcp-parallel.hpp mcp-common.hpp
	$(CXX) -pthread -c -O4 -o mcp-parallel-pthread.o mcp-parallel.cpp

common-pthread.o: mcp-common.cpp mcp-common.hpp
	$(CXX) -pthread -c -O4 -o mcp-common-pthread.o mcp-common.cpp

matrix+formula-pthread.o: mcp-matrix+formula.cpp mcp-matrix+formula.hpp
	$(CXX) -pthread -c -O4 -o mcp-matrix+formula-pthread.o mcp-matrix+formula.cpp

#---------------------------------------------------------------------------------------------------

hybrid: matrix+formula-hybrid.o common-hybrid.o parallel-hybrid.o posix-hybrid.o hybrid.o
	mpic++ -pthread -O4 -o $(BIN)/mcp-hybrid mcp-hybrid.o \
						 mcp-common-hybrid.o \
					         mcp-matrix+formula-hybrid.o \
						 mcp-parallel-hybrid.o \
						 mcp-posix-hybrid.o

hybrid.o: mcp-hybrid.cpp mcp-common.hpp mcp-parallel.hpp mcp-posix.hpp
	mpic++ -pthread -O4 -c -o mcp-hybrid.o mcp-hybrid.cpp


posix-hybrid.o: mcp-posix.cpp mcp-posix.hpp mcp-parallel.hpp mcp-common.hpp
	mpic++ -pthread -c -O4 -o mcp-posix-hybrid.o mcp-posix.cpp

parallel-hybrid.o: mcp-parallel.cpp mcp-parallel.hpp mcp-common.hpp
	mpic++ -pthread -c -O4 -o mcp-parallel-hybrid.o mcp-parallel.cpp

common-hybrid.o: mcp-common.cpp mcp-common.hpp
	mpic++ -pthread -c -O4 -o mcp-common-hybrid.o mcp-common.cpp

matrix+formula-hybrid.o: mcp-matrix+formula.cpp mcp-matrix+formula.hpp
	mpic++ -pthread -c -O4 -o mcp-matrix+formula-hybrid.o mcp-matrix+formula.cpp

#---------------------------------------------------------------------------------------------------

trans: mcp-trans.cpp matrix+formula-trans.o
	$(CXX) -o $(BIN)/mcp-trans mcp-trans.cpp mcp-matrix+formula-trans.o
matrix+formula-trans.o: mcp-matrix+formula.cpp mcp-matrix+formula.hpp
	$(CXX) -c -O4 -o mcp-matrix+formula-trans.o mcp-matrix+formula.cpp

#---------------------------------------------------------------------------------------------------

guess: mcp-guess.cpp
	$(CXX) -o $(BIN)/mcp-guess mcp-guess.cpp

#---------------------------------------------------------------------------------------------------

cnf: mcp-cnf
	cp -f mcp-cnf $(BIN)/mcp-cnf

#---------------------------------------------------------------------------------------------------

split: mcp-split.cpp
	$(CXX) -o $(BIN)/mcp-split mcp-split.cpp

#---------------------------------------------------------------------------------------------------

check:  matrix+formula-check.o mcp-check.o
	$(CXX) -o $(BIN)/mcp-check mcp-check.o \
				    mcp-matrix+formula-check.o

mcp-check.o: mcp-check.cpp mcp-matrix+formula.hpp
	$(CXX) -c -o mcp-check.o mcp-check.cpp

matrix+formula-check.o: mcp-matrix+formula.cpp mcp-matrix+formula.hpp
	$(CXX) -c -o mcp-matrix+formula-check.o mcp-matrix+formula.cpp

#---------------------------------------------------------------------------------------------------

sparse: matrix+formula-sparse.o mcp-sparse.o
	$(CXX) -o $(BIN)/mcp-sparse mcp-sparse.o \
				   mcp-matrix+formula-sparse.o

mcp-sparse.o: mcp-sparse.cpp mcp-matrix+formula.hpp
	$(CXX) -c -o mcp-sparse.o mcp-sparse.cpp

matrix+formula-sparse.o: mcp-matrix+formula.cpp mcp-matrix+formula.hpp
	$(CXX) -c -o mcp-matrix+formula-sparse.o mcp-matrix+formula.cpp

#---------------------------------------------------------------------------------------------------

predict: matrix+formula-predict.o mcp-predict.o
	$(CXX) -o $(BIN)/mcp-predict mcp-predict.o \
				    mcp-matrix+formula-predict.o

mcp-predict.o: mcp-predict.cpp mcp-matrix+formula.hpp
	$(CXX) -c -o mcp-predict.o mcp-predict.cpp

matrix+formula-predict.o: mcp-matrix+formula.cpp mcp-matrix+formula.hpp
	$(CXX) -c -o mcp-matrix+formula-predict.o mcp-matrix+formula.cpp

#---------------------------------------------------------------------------------------------------

chk2tst: matrix+formula-chk2tst.o mcp-chk2tst.o
	$(CXX) -o $(BIN)/mcp-chk2tst mcp-chk2tst.o \
				    mcp-matrix+formula-chk2tst.o

mcp-chk2tst.o: mcp-chk2tst.cpp mcp-matrix+formula.hpp
	$(CXX) -c -o mcp-chk2tst.o mcp-chk2tst.cpp

matrix+formula-chk2tst.o: mcp-matrix+formula.cpp mcp-matrix+formula.hpp
	$(CXX) -c -o mcp-matrix+formula-chk2tst.o mcp-matrix+formula.cpp

#---------------------------------------------------------------------------------------------------

uniq: matrix+formula-uniq.o mcp-uniq.o
	$(CXX) -o $(BIN)/mcp-uniq mcp-uniq.o \
			mcp-matrix+formula-uniq.o

mcp-uniq.o: mcp-uniq.cpp mcp-matrix+formula.hpp
	$(CXX) -c -o mcp-uniq.o mcp-uniq.cpp

matrix+formula-uniq.o: mcp-matrix+formula.cpp mcp-matrix+formula.hpp
	$(CXX) -c -o mcp-matrix+formula-uniq.o mcp-matrix+formula.cpp

#---------------------------------------------------------------------------------------------------

export:
	cp -f $(BIN)/mcp-* ..

clean:
	rm -f *.d *.o $(BIN)/mcp-*

scratch: clean
	rm -f *~

#---------------------------------------------------------------------------------------------------
# EOF
