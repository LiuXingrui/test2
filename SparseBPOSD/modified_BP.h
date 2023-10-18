#include <iostream>
#include<cstring>
#include<string>
#include<vector>
#include <fstream>
#include <math.h> 
#include <random>
#include <sstream>

#include"BP.h"
#include "Sparse.h"

#ifndef MBP
#define MBP
/** @brief Belief-propagation functions
   
    @author Xingrui Liu <xliu@ucr.edu>
    @date August 2022
    */

using namespace std;
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
using namespace itpp;

/** @brief Modified Belief-propagation functions
   
    @author Xingrui Liu <xliu@ucr.edu>
    @date February 2023
    */

/*! \relates nodes
 *Perform the matrix reduction for H and D, for stabilizer codes, if H=Hx, D=Hz, I call it D because it represents the degeneracy

\param A stores the columns of D which are added.
\param num_of_row_reduction is the number of how many rows are reduced.
 */
//
#ifndef modified_BP
#define modified_BP
extern int protocol_B_suc;
extern int protocol_B_max_fail;
extern int protocol_B_syn_fail;
extern int max_fail;
extern int syn_fail;
extern int OSD_suc;
extern int  fix_suc;
extern int fix_fail;
extern int converge_to_trapping_set;
extern int did_not_converge_for_large_iterations;
extern int OSD_fail;
extern int syn_tilde_suc;
extern int  num_no_error;
extern int OSD_0_suc;
extern int OSD_1_suc;
extern int OSD_2_suc;
extern int OSD_higher_order_suc;

extern ivec global_sort_index;
extern ivec global_perm2;
extern GF2mat global_e_S;
extern int global_new_wt;

extern int  Ran_order_0_suc;
extern int  Ran_order_1_suc;
extern int  Ran_order_2_suc;
extern int  Ran_order_higher_order_suc;

extern int wt_1_residual_error;
extern int wt_2_residual_error;
extern int large_wt_residual_error;
extern int avg_red_qubits;
extern int BP_suc;
extern int avg_BP_suc;
extern int avg_OSD_suc;
extern int freeze_rows_suc;
#endif


void Mat_Trans(const GF2mat& H, const GF2mat& D, vector<double> &K, vector<double>& K_tilde, GF2mat& H_tilde, GF2mat& D_tilde, vector<vector<int>>& A, int num_of_row_reduction,int debug=0);


void row_reduction(const int which_row, GF2mat& H, GF2mat& D, GF2mat& H2, GF2mat& D2, vector<double> &input_K,vector<double>&output_K,vector<vector<int>>& A, int debug=0);

void construct_B_tau(vector<double> &B_tau,int w, vector<double>& input_K, vector<int> Ai,vector<vector<int>>&  Tau);

void Add_cols(GF2mat& H,GF2mat& D,vector<int>& Ai,const int w,int debug=0);

void K_trans(vector<double>& input_K, vector<double> &output_K, vector<int>& Ai, int w,vector<vector<int>>& b);

void Add_cols_and_rows(GF2mat& H,GF2mat& D,vector<int>& Ai,int w,GF2mat& H2,GF2mat& D2,int debug,vector<vector<int>>& b);

void construct_BF(GF2mat& D,vector<int> &Ai, int w,int c, GF2mat& DB, GF2mat& F,vector<vector<int>>& b,int n,int debug=0);


/*! \relates nodes
 *useless decoder
 */
//bool new_decoder1(GF2mat& H2,GF2mat &H, GF2mat &G,const nodes checks[],const nodes errors[],const vec &pv,const vec& pv_dec,double& num_iter, int lmax,int& max_fail,int&syn_fail,int debug, vec &LR,int rankH,int options);


/*! \relates nodes
 *useless decoder
 */
bool modified_decoder(GF2mat& H,GF2mat& G,GF2mat &H_tilde, const vector<vector<int>>& A,const nodes checks[],const nodes errors[],const vec &pv,const vec& pv_dec,const int&  num_row_red,double& num_iter, int lmax,int debug, vec &LR,int rankH);
void inverse_trans(GF2mat& output_e,const GF2mat& output_e_tilde,const vector<vector<int>>& A);

int Weight(const bvec &cw);

double energy (bool empty_G,const GF2mat& G, const GF2mat& eT, const vector<double> &K);
double ML_suc_rate (vec p,const GF2mat& D,const GF2mat& H, const GF2mat& H_tilde_star, vector<double>&K,int d,int num_large_wt_error);

void ML_decoder_verify(vec p,const GF2mat& H, const GF2mat D,const GF2mat& H_star, const GF2mat D_star,const GF2mat H_tilde, const GF2mat H_tilde_star,vector<double>&K,int max_num_cws,double& after_trans_suc_rate,double&suc_rate,double &after_trans_theoric_suc_rate,double& min_wt_suc_count,int debug=0);

void ML_decoder(const GF2mat& input_e, const GF2mat& D,const GF2mat& L,vector<double>& K,GF2mat& output_e);
void algebraic_decoder(const GF2mat& H,const GF2mat& syndrome,GF2mat &output_e);



bool   alternative_decoder(const GF2mat &Hx, const GF2mat &Hz,const GF2mat &Gx,GF2mat G,const GF2mat &real_e,GF2mat &final_output_e,const nodes checks1[],const nodes errors1[],const vec &pv, vec&pv_dec,const double& pmin,const double& pmax,double& num_iter, int lmax,int &wt, int debug, vec &LR2,int rankH,double M,int&LLR_fail,int& LLR_suc,int&num_fixed_LLR,bool& start_row_reduction,double alpha,int other_decoder,int channel=0, int max_row_freeze=0,int max_rep=1);
void del_col (GF2mat & H, int col_ind);
void restore_e(GF2mat& output_e,  const GF2mat& original_output_e, const vector<vector<vector<int>>> &pos_deleted_entries,const vector<vector<int>>& value_deleted_entries);
void insert_col(GF2mat& output_e,int which_col,const GF2mat& B);
void insert_row(GF2mat& output_e,int which_row,const GF2mat& B);
void small_row_reduction(GF2mat &H,GF2mat& G, GF2mat& output_e, int col_ind,  vector<vector<int>>& Ai, mat& mcv, mat& mvc);
void del_row (GF2mat & H, int row_ind);

/*
mat merge_mat_hori(const mat &left,const mat &right);
void del_col (mat & H, int col_ind);
*/

void right_small_row_reduction(bool& empty_G,GF2mat& syndrome,GF2mat &H,GF2mat& G, GF2mat& output_e, int col_ind,  vector<vector<int>>& Ai, mat& mcv, mat& mvc, vector<GF2mat> & rows_del,vector<int>& B,int& hardSelectValue,int&LLR_fail,int& LLR_suc,int debug);
void right_restore_e(GF2mat& output_e,  const GF2mat& original_output_e, const vector<vector<int>>&pos_deleted_entries, const vector<int>& value_deleted_entries,vector<GF2mat>& rows_del_list);
GF2mat row_gaussian(const GF2mat H, itpp::GF2mat& Perm) ;
//perform column gaussian elimination on H, return the matrix P such that HP is the column gaussian elimination of H,  and also permute all all-zero columns 
//to the right side of H
GF2mat column_gaussian(itpp::GF2mat H,int rankH) ;
GF2mat row_gaussian(itpp::GF2mat H,int rankH) ;
bool  Fast_OSD(const GF2mat &Hx, const GF2mat &Gx,const nodes checks[],const nodes errors[],const vec &pv, vec&pv_dec,const double& pmin,
const double& pmax,double& num_iter, int lmax,int &wt,  int debug, vec &LR2,int rankH,int& max_rep);
ivec swap_zero_rows_to_bottom(GF2mat& H,int rankH);

ivec reliability_sort_index (const GF2mat & H, const vec& LR);
void reliability_perm_H (const GF2mat & H, int which_row,GF2mat &output_H,GF2mat& output_col_perm,GF2mat& output_row_perm, GF2mat& H2,GF2mat& A,GF2mat& output_H1);
bool reliability_solve(const GF2mat&H1,const GF2mat&A, const GF2mat &s1, const GF2mat &e2,GF2mat & e1);
GF2mat reliability_restored_e(const GF2mat col_perm, const GF2mat &e1,const GF2mat &e2);
bool   circuit_decoder(const Sparse_GF2 &sparseH,const Sparse_GF2 &sparseL,const GF2mat &Hx, const GF2mat Gx, GF2mat G,const GF2mat &real_e,GF2mat &final_output_e,const nodes checks1[],const nodes errors1[],const vec&pv_dec,double& num_iter, int lmax,int &wt,  int debug, vec &LR2,int rankH,double M,int&LLR_fail,int&LLR_suc,int& num_fixed_LLR,bool& start_row_reduction,double alpha,int other_decoder,int channel,int max_row_freeze,int OSD_order=100);

bool MBP_del_cols(const Sparse_GF2 &sparse_MH,const Sparse_GF2 &sparse_ML,const Sparse_GF2 &sparseH,
const Sparse_GF2 &sparseL,const Sparse_GF2&Sparse_real_Me,
const GF2mat &real_e,const std::vector<vector<int>>& Positions,const nodes checks[],const nodes errors[],const vec&pv_dec,double& num_iter, 
const int& lmax,const int &wt,  const int& debug, const int& schedule,const int& OSD_order,Sparse_GF2& Sparse_zero_rvec, GF2mat& zero_mat1);
#endif