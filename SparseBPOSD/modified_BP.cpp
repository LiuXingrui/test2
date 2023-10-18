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


#include <itpp/itbase.h>  // For the `vec` type from the IT++ library

// ... (Other includes and previous code)

void Sparse_GF2::del_identical_cols(Sparse_GF2& B,int max_cols, std::vector<vector<int>>& Positions, itpp::vec& p) {
    int num_identical_pairs = 0;
    
    for (int i = 0; i < num_cols - 1 && num_identical_pairs < max_cols; i++) {
        for (int j = i + 1; j < num_cols && num_identical_pairs < max_cols; j++) {
           
		   if (col_mat[i] == col_mat[j]) {  // If columns are identical
                // Delete j-th column from A
				vector<int> identical_pair;
				identical_pair.push_back(i);
				identical_pair.push_back(j);
				
                col_mat.erase(col_mat.begin() + j);
                num_cols--;
				B.col_mat.erase(B.col_mat.begin() + j);
                B.num_cols--;

                // Update the row representation
				col_to_row();
				B.col_to_row();

                // Update vector p
                p[i] = p[i] + p[j] - p[i] * p[j];
                p.del(j);

                // Store the position
                Positions.push_back(identical_pair);
                
                num_identical_pairs++;
            }
        }
    }
}

Sparse_GF2 Sparse_GF2::res_deleted_cols(const std::vector<vector<int>>& Positions) const{
    if (num_cols != 1) {
        throw std::runtime_error("res_deleted_cols:  should have only one column.");
    }
    Sparse_GF2 E = *this;
    for (int i = Positions.size() - 1; i >= 0; i--) {
        std::vector<int> zero_row;
        E.row_mat.insert(E.row_mat.begin() + Positions[i][1], zero_row);
        E.num_rows++;
    }

    // Update the column representation
    E.row_to_col();
	return E;
}

//error is a column vector
void Sparse_GF2::error_del_cols(const std::vector<vector<int>>& Positions) {
    if (num_cols != 1) {
        throw std::runtime_error("error_del_cols:  should have only one column.");
    }

    for (int i = 0; i <Positions.size() ; i++) {
	//d_rows(Positions[i][1],Positions[i][0]);
		row_mat.erase(row_mat.begin() +  Positions[i][1]);
        num_rows--;
		
    }
row_to_col();
    // Update the column representation
}

bool MBP_del_cols(const Sparse_GF2 &sparse_MH,const Sparse_GF2 &sparse_ML,const Sparse_GF2 &sparseH,
const Sparse_GF2 &sparseL,const Sparse_GF2&Sparse_real_Me,
const GF2mat &real_e,const std::vector<vector<int>>& Positions,const nodes checks[],const nodes errors[],const vec&pv_dec,double& num_iter, 
const int& lmax,const int &wt,  const int& debug, const int& schedule,const int& OSD_order,Sparse_GF2& Sparse_zero_rvec, GF2mat& zero_mat1)
{
	int c=sparseH.rows();
 int v=sparseH.cols();
 int Mv=sparse_MH.cols();

 
 //if no error, break

  Sparse_GF2 Sparse_real_e(real_e);
  Sparse_GF2 Sparse_syndrome=sparseH*Sparse_real_e;

  Sparse_GF2 Sparse_Msyndrome=sparse_MH*Sparse_real_Me;

  
  GF2mat syndrome=Sparse_syndrome.to_GF2mat();
  
  //GF2mat Msyndrome=Sparse_Msyndrome.to_GF2mat();
	
  //is the syndrome a zero vector?
  if (syndrome==zero_mat1)
	{
	  // is e a stablizer?
	  if (sparseL*Sparse_real_e==Sparse_zero_rvec)
	    {
			num_no_error++;
			return true;
	    }        
	  else  
	    {
	      syn_fail++;
	      return false;	     
	    }  
	}
	
	//this part can be optimized:
      mat mcv(c,Mv);   // the messages from check nodes to variable nodes, which are p0/p1
      mat mvc(c,Mv);   // the messages from variable nodes to check nodes
      mcv.zeros();
      mvc.zeros();
	  
	  GF2mat MH=sparse_MH.to_GF2mat();
      initialize_massages( mcv,mvc, MH); //initialize to all-1 matrix	
	  
	  //GF2mat output_e(v,1);	
	  vec LR(Mv);	
	
      for (int l=1;l<=lmax;l++)
	{		 
			Sparse_GF2 Sparse_output_Me(Mv,1);
			Sparse_GF2 Sparse_output_e(v,1);
			if (schedule==0)
			{
				quan_p_update(checks,errors, mcv,mvc,syndrome,pv_dec, c, Mv,Sparse_output_Me,LR,debug);				
			}
	      else if (schedule==1)
			{		
				quan_s_update(checks,errors, mcv,mvc,syndrome,pv_dec, c, Mv,Sparse_output_Me,LR,debug);
			}
			 else if (schedule==2)
			 {
				 quan_s_C_update(checks,errors, mcv,mvc,syndrome,pv_dec, c, Mv,Sparse_output_Me,LR,debug);
			 }
			 else if (schedule==3)
			{		
				quan_Ran_s_update(checks,errors, mcv,mvc,syndrome,pv_dec, c, Mv,Sparse_output_Me,LR,debug);
			}
	    
			//Sparse_GF2 Sparse_output_e(output_e);
			
			 //sparseH*Sparse_output_e==Sparse_syndrome
	
	
			 
			if ( sparse_MH.checkProductEquality(Sparse_output_Me, Sparse_syndrome) )		   
			{				
					
					//sparseL*(Sparse_output_e+Sparse_real_e)==Sparse_zero_rvec
				  if(sparse_ML.checkProductEquality(Sparse_output_Me+Sparse_real_Me, Sparse_zero_rvec))
					{
					
						Sparse_output_e=Sparse_output_Me.res_deleted_cols(Positions);
						 if(sparseL.checkProductEquality(Sparse_output_e+Sparse_real_e, Sparse_zero_rvec))
						 {
							
							num_iter= num_iter+l;		
							BP_suc++;
							return true;
						 }
						 else
							{ 
						         if (debug &1)
								 {
									 cout<<"BP givies an error with the same syndrome but failed: output_e is \n"<<endl;
									 print_surf_lattice(Sparse_output_e); 
									  cout<<"real_e is \n"<<endl;
									 print_surf_lattice(Sparse_real_e);
									 cout<<"\nor "<<endl;
									Sparse_real_e.print_cols();
									
									 cout<<"transformed real_e is "<<endl;
									Sparse_real_Me.print_cols();
									
									 
									 cout<<"\noriginal syndrome is \n "<<endl;
									 Sparse_syndrome.print_cols();
									 /*
									 cout<<"transformed syndrome is \n "<<endl;
									 Sparse_Msyndrome.print_cols();
									 */
								 }
								syn_fail++;
								return false; 
							}
					}	    
				else
				{ 
					if (debug &1)
						{
									 cout<<"BP givies an error with the same syndrome but failed: output_e is \n"<<endl;
									 print_surf_lattice(Sparse_output_e); 
									  cout<<"real_e is \n"<<endl;
									 print_surf_lattice(Sparse_real_e);
									 cout<<"\nor "<<endl;
									Sparse_real_e.print_cols();
									
									 cout<<"transformed real_e is "<<endl;
									Sparse_real_Me.print_cols();
									
									 
									 cout<<"\noriginal syndrome is \n "<<endl;
									 Sparse_syndrome.print_cols();
									 /*
									 cout<<"transformed syndrome is \n "<<endl;
									 Sparse_Msyndrome.print_cols();
									 */
						}
					syn_fail++;
					return false;
				}	    	  
		}
	}
      						
      if (debug&8)
	{
		int suc_order=0;
			  Sparse_GF2 Sparse_output_MeOSD(Mv,1);
			  Sparse_GF2 Sparse_output_eOSD(v,1);
			  
	  Sparse_OSD(sparse_MH,LR,MH,syndrome,Sparse_output_MeOSD,suc_order,OSD_order);

	Sparse_output_eOSD=Sparse_output_MeOSD.res_deleted_cols(Positions);
	
	  if(sparseL.checkProductEquality(Sparse_output_eOSD+Sparse_real_e, Sparse_zero_rvec))
	    {
	      OSD_suc++;
		  if (suc_order==0){OSD_0_suc++;}
		else if(suc_order==1){OSD_1_suc++;}
		else if(suc_order==2){OSD_2_suc++;}
		else {OSD_higher_order_suc++;}
	      return true;
	    }	 
		else 
		{								         
			if (debug &1)
				 {
					 cout<<"original syndrome is \n "<<endl;
					Sparse_syndrome.print_cols();
					cout<<"transformed syndrome is \n "<<endl;
					Sparse_Msyndrome.print_cols();
									 
					cout<<"OSD fails: output_e is \n"<<endl;
					print_surf_lattice(Sparse_output_eOSD);
					cout<<"real_e is \n"<<endl;
					print_surf_lattice(Sparse_real_e);
			}
			OSD_fail++;
			syn_fail++;			
			return false;
		}
	}
      max_fail++;
      return false;
 }
 