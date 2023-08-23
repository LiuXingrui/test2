#include <iostream>
#include<cstring>
#include<string>
#include<vector>
#include <fstream>
#include <math.h> 
#include <random>
#include <sstream>
#include <chrono>
#include <algorithm>
#include <itpp/itbase.h>
#include <itpp/itcomm.h>

#include"BP.h"

using namespace std;
using namespace itpp;

#ifndef SPARSE_GF2_H
#define SPARSE_GF2_H

class Sparse_GF2 {
private:
    std::vector<std::vector<int>> row_mat;
    std::vector<std::vector<int>> col_mat;
    int num_rows;
    int num_cols;

public:
    Sparse_GF2();
    Sparse_GF2(const itpp::GF2mat& mat);
    Sparse_GF2(int a, int b);
    
    // Declare all other member functions here, but don't define them.
    int rows() const;
    int cols() const;
    std::vector<int> get_row(const int& i) const;
    std::vector<int> get_col(const int& j) const;
   Sparse_GF2 operator*(const Sparse_GF2& other) const;
    bool operator==(const Sparse_GF2& other) const;
	bool operator!=(const Sparse_GF2& other) const;
     Sparse_GF2 operator+(const Sparse_GF2& other) const;
     void row_to_col();
      void col_to_row();
       void swap_rows(const int& i, const int& j,const bool& update_col=true);
       void swap_cols(const int& i, const int& j,const bool &update_row=true);
        itpp::GF2mat to_GF2mat() const;
        void add_rows(const int& i, const int& j,const bool& update_col=true);
        void add_cols(const int& i, const int& j,const bool &update_row=true);
        Sparse_GF2 row_gaussian(int& rankH) const;
        Sparse_GF2 inverse() const;
		
		
        Sparse_GF2 col_gaussian(int& rankH) const;	
		Sparse_GF2 permute_cols(const itpp::ivec& Perm) ;
        static Sparse_GF2 col_permutation_matrix(const itpp::ivec& perm) ;
        Sparse_GF2 transpose() const;
        void del_row(const int& row_idx);
         void del_col(const int&  col_idx);
         Sparse_GF2 get_submatrix(const int&  begin_row, const int&  begin_col,const int&  end_row, const int&  end_col) ;
          int col_weight(const int&  i) const;
           int row_weight(const int&  i) const;
           void set(const int&  i, const int&  j, const int&  value);  
};
void Sparse_OSD(const Sparse_GF2& H,vec& LR,const GF2mat& denseH,const GF2mat& syndrome,Sparse_GF2& Sparse_output_eOSD,int& suc_order,int OSD_order=100);
bool Fcircuit_decoder(const Sparse_GF2 &sparseH,const Sparse_GF2 &sparseL,const GF2mat &Hx, const GF2mat& Gx, const GF2mat& G,const GF2mat &real_e,const nodes checks[],const nodes errors[],const vec&pv_dec,double& num_iter, const int& lmax,const int &wt,  const int& debug, const int& OSD_order,Sparse_GF2& Sparse_zero_rvec, GF2mat& zero_mat1);
 void quan_s_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const GF2mat& syndrome,const vec &pv,int&c, int& v,  Sparse_GF2& output_e, vec &LR,const int &debug,double alpha=1);
 void quan_p_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const GF2mat& syndrome,const vec &pv,int c, int v,  Sparse_GF2& output_e, vec &LR,int debug,double alpha=1);
#endif // SPARSE_GF2_H
