#include <iostream>
#include<cstring>
#include<string>
#include<vector>
#include <fstream>
#include <math.h> 
#include <random>
#include <sstream>
#include <bits/stdc++.h>


#include"BP.h"
#include"modified_BP.h"
#include "Sparse.h"
using namespace std;
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
using namespace itpp;

int protocol_B_suc=0;
int protocol_B_max_fail=0;
int protocol_B_syn_fail=0;
int max_fail=0;
int syn_fail=0;
int OSD_suc=0;
int  fix_suc=0;
int fix_fail=0;
int converge_to_trapping_set=0;
int did_not_converge_for_large_iterations=0;
int OSD_fail=0;
int syn_tilde_suc=0;
int num_no_error=0;
int OSD_0_suc=0;
int OSD_1_suc=0;
int OSD_2_suc=0;
int OSD_higher_order_suc=0;

ivec global_sort_index;
ivec global_perm2;
GF2mat global_e_S;
int global_new_wt=0;

int  Ran_order_0_suc=0;
int  Ran_order_1_suc=0;
int  Ran_order_2_suc=0;
int  Ran_order_higher_order_suc=0;

int wt_1_residual_error=0;
int wt_2_residual_error=0;
int large_wt_residual_error=0;
int avg_red_qubits=0;
int BP_suc=0;
int avg_BP_suc=0;
int avg_OSD_suc=0;
int freeze_rows_suc=0;
	
int Weight(const bvec &cw)
{
  int n=cw.size();
  int wt=0;
  for (int i=0;i<n;i++)
    {
      if(cw(i)==1){wt++;}
    }
  return wt;
}






//row echelon form by row permutation
GF2mat row_gaussian(const GF2mat H, itpp::GF2mat& Perm) 
{
	
	//cout<<"meet an empty col of G:"<<endl;
    GF2mat H_tilde = H; // initialize H_tilde as H
	bmat P_bmat;
	P_bmat=eye_b(H.rows());
	GF2mat P(P_bmat);
    int n_rows = H_tilde.rows();
   // int n_cols = H_tilde.cols();

    // Perform row echelon form
    int pivot_col=0;
	int pivot_row=0;
   while (pivot_row<H_tilde.rows() and pivot_col<H_tilde.cols())
	{
        int temp=-1;
		 
	   for (int i=pivot_row;i<n_rows;i++)
	   {
			if(H_tilde.get(i, pivot_col) == 1) {temp=i;break;}
	   }
	   
	  // cout<<"pivot is "<<pivot_row<<endl;
	    
	    if (temp==-1){pivot_col++;continue;}
	
		H_tilde.swap_rows(temp,pivot_row);		
		P.swap_rows(temp,pivot_row);
		
		for (int r=0;r<n_rows;r++)
		{
			if (H_tilde.get(r,pivot_col)==1 and r!=pivot_row) {H_tilde.add_rows(r,pivot_row); P.add_rows(r,pivot_row);}
		}
		pivot_col++;
		pivot_row++;
		
	}
	
	Perm=P;
	//if (P*H==H_tilde){}
	//else{}
//	cout<<"row_echelon_form: H_tilde is\n"<<H_tilde<<"\n P*H is\n"<<P*H<<endl;
    return H_tilde;
}

//perform row gaussian elimination on H, return the matrix P such that PH is the column gaussian elimination of H,  and also permute all all-zero rows
//to the bottom of H
GF2mat row_gaussian(itpp::GF2mat H,int rankH) {
  int r = H.rows();
  int c = H.cols();
  
 // if (c>r){}
  //else{cout<<"permute_GF2mat_columns error:H.cols()<H.rows()"<<endl;}
  // Initialize the permutation vector as an identity permutation
  GF2mat permutation(r,r);
  for (int i = 0; i < r; ++i) {
    permutation.set(i,i,1);
  }

int pivot_col=0;
int pivot_row=0;
while (pivot_row<r)
{
	pivot_col=-1;
	for (int i=0;i<c;i++)
	{
		if (H(pivot_row,i)==1){pivot_col=i;break;}
	}
 
	if (pivot_col!=-1)
	{
		for (int i=0;i<r;i++)
		{
			if (H(i,pivot_col)==1 &&i!=pivot_row)
			{
				for (int j=0;j<c;j++) 
				{
					H.set(i,j,H(i,j)+H(pivot_row,j));
					
				}
					for (int j=0;j<r;j++) 
				{
					permutation.set(i,j,permutation(i,j)+permutation(pivot_row,j));			
				}
			}
		}
	}
	pivot_row++;
}

	
		int current_row=0;

		int swaped_row=rankH;
		while (current_row<rankH and swaped_row<r)
		{
			bool empty_row=true;
			for (int i=0;i<c;i++)
			{
				if (H(current_row,i)==1){empty_row=false;break;}
			}
			if (empty_row)
			{
			H.swap_rows(current_row, swaped_row);
			permutation.swap_rows(current_row, swaped_row);
			
			swaped_row++;
			}
			else 
			{
			current_row++;
			}
		}
	
	
  return permutation;
}

//perform column gaussian elimination on H, return the matrix P such that HP is the column gaussian elimination of H,  and also permute all all-zero columns 
//to the right side of H
GF2mat column_gaussian(itpp::GF2mat H,int rankH) {
  int r = H.rows();
  int c = H.cols();
  
 // if (c>r){}
  //else{cout<<"permute_GF2mat_columns error:H.cols()<H.rows()"<<endl;}
  // Initialize the permutation vector as an identity permutation
  GF2mat permutation(c,c);
  for (int i = 0; i < c; ++i) {
    permutation.set(i,i,1);
  }
  
int pivot_col=0;
int pivot_row=0;
while (pivot_col<c)
{
	pivot_row=-1;
	for (int i=0;i<r;i++)
	{
		if (H(i,pivot_col)==1){pivot_row=i;break;}
	}
 
	if (pivot_row!=-1)
	{
		for (int i=0;i<c;i++)
		{
			if (H(pivot_row,i)==1 &&i!=pivot_col)
			{
				for (int j=0;j<r;j++) 
				{
					H.set(j,i,H(j,i)+H(j,pivot_col));
					
				}
					for (int j=0;j<c;j++) 
				{
					permutation.set(j,i,permutation(j,i)+permutation(j,pivot_col));			
				}
			}
		}
	}
	pivot_col++;
}

		int current_col=0;

		int swaped_col=rankH;
		while (current_col<rankH and swaped_col<c)
		{
			bool empty_col=true;
			for (int i=0;i<r;i++)
			{
				if (H(i,current_col)==1){empty_col=false;break;}
			}
			if (empty_col)
			{
			H.swap_cols(current_col, swaped_col);
			permutation.swap_cols(current_col, swaped_col);
			
			swaped_col++;
			}
			else 
			{
			current_col++;
			}
		}
	
  return permutation;
}



