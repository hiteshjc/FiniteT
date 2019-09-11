#ifndef QDDM_HEADER
#define QDDM_HEADER

#include"lanczos_more.h"
#include"ed.h"
#include"printing_functions.h"
#include"density_matrix.h"
#include"hamiltonian_spin_functions.h"

void calc_qddm_fermions(Ham &h,
                        int iterations, 
                        std::vector<double> &eigs);
void calc_qddm_spins(Ham &h,
                        int iterations, 
                        std::vector<double> &eigs);
void calc_qddm_spin1(int sz,Ham &h,
               int iterations, 
	       std::vector<int> kept_indices,
               std::vector<double> &eigs, int nlow);
void calc_qddm_spin1b2(int sz,Ham &h,
               int iterations, 
	       std::vector<int> kept_indices,
               std::vector<double> &eigs, int nlow);


#endif
