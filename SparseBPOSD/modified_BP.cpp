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
	


//Perform the matrix reduction for H and D, for stabilizer codes, if H=Hx, D=Hz, I call it D because it represents the degeneracy
void Mat_Trans(const GF2mat& H, const GF2mat& D, vector<double> &K, vector<double>& K_tilde, GF2mat& H_tilde, GF2mat& D_tilde, vector<vector<int>>& A, int num_of_row_reduction,int debug){

  GF2mat input_H=H;
  GF2mat input_D=D;
  GF2mat output_H, output_D;
  vector<double> input_K=K;
  vector<double> output_K=K;
  // cout<<"start"<<endl;
  
  int real_num_of_row_reduction;
  if (num_of_row_reduction>=0){real_num_of_row_reduction=num_of_row_reduction;}
  else{real_num_of_row_reduction=-num_of_row_reduction;} //negative num_of_row_reduction for combine e at the bottom of D


    if (debug&4)
	{
	  GF2mat zero_mat1(output_H.rows(),output_D.rows());
	  cout<<"before row reduction H is"<<H<<" \n D: \n"<<D<<endl;	  
	  cout<<"K_tilde is \n:"<<endl;
	  for(auto ii:input_K)
	    {
	      cout<<ii<<" ";
	    }
	  cout<<endl;
	}
      
  for (int i=0;i<real_num_of_row_reduction;i++)
    {
      // cout<<i<<"th row reduction"<<endl;
      int min_wt=H.cols();
      int this_row=0;

      if (num_of_row_reduction>=0)
		{
			for (int j=0;j<input_D.rows();j++)
			{  
				bvec temp_vec=input_D.get_row(j);
				int temp_wt=Weight(temp_vec);
				if (temp_wt<min_wt){min_wt=temp_wt;this_row=j;}
			}
		}
      
      else
		{
			for (int j=0;j<input_D.rows()-1;j++) //donot touch the last row- the error vector
			{	  
				bvec temp_vec=input_D.get_row(j);
				int temp_wt=Weight(temp_vec);
				if (temp_wt<min_wt){min_wt=temp_wt;this_row=j;}
			}
		}
      
   
      if (debug&4)
		{
			cout<<i+1<<"th row reduction, delete "<<this_row<<"th row:"<<endl;
		}
      //perform row reduction for the column with minimum weight of D
      // cout<<121<<endl;
      if (min_wt>12){cout<<"for "<<i<<"th row reduction (first index is 0) weight of minimum-weight row of D is "<<min_wt<<", so stop here"<<endl;break;}
      row_reduction(this_row,input_H, input_D, output_H,output_D,input_K,output_K,A,debug);//I havn't delete the rows and columns
      
      //  cout<<22<<endl;
      //  check if the prog is right:
    
			GF2mat zero_mat1(output_H.rows(),output_D.rows());
			//cout<<"after "<<i<<"th row reduction before deleting rows/cols \n H:\n"<<output_H<<" \n D: \n"<<output_D<<endl;
			if((output_H*output_D.transpose())==zero_mat1) {}
			else{cout<<"before deleting: output_H*output_D^T!=0 for the row "<<i<<endl;}
	
      

      if (A[i][0]==0) //if delete the 0th column from H and D
		{
			output_H=output_H.get_submatrix(0,1,output_H.rows()-1,output_H.cols()-1);
			if (i!=D.rows()-1) { output_D=output_D.get_submatrix(0,1,output_D.rows()-1,output_D.cols()-1);      }
		}
      
      else if(A[i][0]==output_H.cols()-1)//if delete the last column
		{
			output_H=output_H.get_submatrix(0,0,output_H.rows()-1,output_H.cols()-2);
			if (i!=D.rows()-1) { output_D=output_D.get_submatrix(0,0,output_D.rows()-1,output_D.cols()-2);}
		}

      
      else //delete A[i][0]th column from H and D
		{
			GF2mat H1=output_H.get_submatrix(0,0,output_H.rows()-1,A[i][0]-1);
			GF2mat H2=output_H.get_submatrix(0,A[i][0]+1,output_H.rows()-1,output_H.cols()-1);
			output_H=merge_mat_hori(H1,H2);
			if (i!=D.rows()-1)
			{
				GF2mat D1=output_D.get_submatrix(0,0,output_D.rows()-1,A[i][0]-1);
				GF2mat D2=output_D.get_submatrix(0,A[i][0]+1,output_D.rows()-1,output_D.cols()-1);
				output_D=merge_mat_hori(D1,D2);
			}	  
		}
      if (i!=D.rows()-1)
		{
			if (this_row==0) //if delete 0th row from D
			{
				output_D=output_D.get_submatrix(1,0,output_D.rows()-1,output_D.cols()-1);
			}
			else if(this_row==output_D.rows()-1)
			{
				output_D=output_D.get_submatrix(0,0,output_D.rows()-2,output_D.cols()-1);
			}
			else 
			{
				GF2mat D3=output_D.get_submatrix(0,0,this_row-1,output_D.cols()-1);
				// cout<<D3.cols()<<endl;
				GF2mat D4=output_D.get_submatrix(this_row+1,0,output_D.rows()-1,output_D.cols()-1);
				// cout<<D4.cols()<<endl;
				output_D=merge_mat_vert(D3,D4);
			}
     
		}
      output_K.erase(output_K.begin()+this_row);

      //need delete an element from K
      
      if (debug&4)
		{
			GF2mat zero_mat1(output_H.rows(),output_D.rows());
			cout<<"after deleting "<<i+1<<" rows (delete row: "<<this_row<<" and column: "<<A[i][0]<<" )\n H:\n"<<output_H<<" \n D: \n"<<output_D<<endl;
			if (i!=D.rows()-1)
			{
				if((output_H*output_D.transpose())==zero_mat1) {}
				else{cout<<"after deleting rows and cols: output_H*output_D^T!=0 for the row the result is"<<output_H*output_D.transpose()<<endl;return;}
			}
	  
			cout<<"K_tilde is \n:"<<endl;
			for(auto ii:output_K){cout<<ii<<" ";}
			cout<<endl;
		}
          
      input_H=output_H;
      input_D=output_D;
      input_K=output_K;
    }

  H_tilde=input_H;
  D_tilde=input_D;
  K_tilde=input_K;
  // cout<<"end"<<endl;
}


void row_reduction(const int which_row, GF2mat& H, GF2mat& D, GF2mat& H2, GF2mat& D2, vector<double> &input_K,vector<double>&output_K,vector<vector<int>>& A, int debug){
  vector<int> Ai;
  vector<vector<int>> b;
  int n=D.cols();
  int k=D.rows();
  int w=0;

  for (int j=0;j<n;j++)
    {
      if (D(which_row,j)==1){Ai.push_back(j);w++;}  // find where are the 1s
    }

  //cout<<1<<"wt is "<<w<<endl;
  // if (w==72){cout<<D<<endl;}
  Ai.push_back(w);
  A.push_back(Ai);
  if (w>1)
    {
      Add_cols(H,D,Ai,w,debug);// now H -> H*A^T^(-1), D-> DA
      //  cout<<2<<endl;
      if (w>2)
	{
	  // cout<<"before add cols and rows"<<endl;
	  Add_cols_and_rows(H,D,Ai,w,H2,D2,debug,b);  // now H-> \frac{H*A^T^(-1)}{F}, D-> DA | DAB
	  //  cout<<3<<endl;
	  // cout<<"after add cols and rows"<<endl;
	}
      else{H2=H;D2=D;}
      K_trans(input_K,output_K,Ai,w,b);
      //  cout<<4<<endl;
      // cout<<"after K_trans"<<endl;
    }
  else{H2=H;D2=D;}
}

  // b stores the binary form of 1,2...2^(w-1)-1
void K_trans(vector<double>& input_K, vector<double> &output_K, vector<int>& Ai, int w,vector<vector<int>>& b){
  int num_change_cols=pow(2,w-1)-1;// w=2: 1; w=3: 3; w=4:7
  int num_B_tau=pow(2,w);
  int change_col=0;
  vector<double> B_tau;
 
  vector<vector<int>>  Tau;

  int temp_K_changed=0;
  
  // get the binray form of i, store in the vector called "tau",notice the ith element in "tau" is the ith digit of i, so if i=1010, tau=[0,1,0,1]
  for (int i=0;i<num_B_tau;i++)
    {
      int tempi=i;
      vector<int> tau;  
      for (int j=0;j<w;j++)
	{
	  int tempi_red=tempi%2;
	  tempi=tempi/2;
	  
	  // if jth digit is one, set corresponding entries in co and ro to be one
	  if (tempi_red==1){tau.push_back(1);}
	  else {tau.push_back(0);};
	}
      Tau.push_back(tau);   
    }
  //now construct all B_tau
  construct_B_tau(B_tau,w,input_K,Ai,Tau);
  /*
  cout<<"w is "<<w<<"size of B_tau is "<<B_tau.size()<<endl;
  cout<<"B_tau is"<<endl;
  
  for (auto ii:B_tau)
    {
      cout<<ii<<" ";
    }
  cout<<endl;
  */
  
  for (int i=1;i<=num_change_cols;i++) // every changed col gives a changed K
    {
      double tempK=1;
      vector<int> i_binary_pos_plus_1; //so that is the position of corresponding tau for K
      int wt_i=0;
      int tempi=i;
      //get the binary form of i,store the position of 1 in "i_binary_pos"
      for (int ii=0;ii<=w-2;ii++)
	{
	  int tempi_red=tempi%2;
	  tempi=tempi/2;
	  if (tempi_red==1){i_binary_pos_plus_1.push_back(ii+1);wt_i++;}
	  
	}
      //get the tau for K
      vector<int> tau_for_K_pos=i_binary_pos_plus_1;
      if (wt_i%2==1) {tau_for_K_pos.push_back(0);}
  
      //because the vector i has weight w-1, but the vector tau has weight w, every even weight tau gives a i, 
      // when convert i back to tau,just add a digit at the end of i, if i has even weight, the first digit of tau=0, if bi has odd weight, the first digit of tau=1.
      double temp=1;
      for (int j=0;j<num_B_tau;j++)
	{
	  int wt_tau_bi=0;
	  for (int ii=0;ii<tau_for_K_pos.size();ii++)
	    {
	      if (Tau[j][tau_for_K_pos[ii]]==1){wt_tau_bi++;}//Tau[j] is the binary form of j, check the number of i such that Tau[j][i]=tau[i]=1
	    }
	  if (wt_tau_bi%2==0){temp=temp*B_tau[j];}
	  else{temp=temp/B_tau[j];}
	}      
      if (wt_i==1){output_K[Ai[temp_K_changed]]=1/pow(2,w)*log(temp);temp_K_changed++;}  //change the value of corresponding K
	else {output_K.push_back(1/pow(2,w)*log(temp));}
      }
      
}
      
void construct_B_tau(vector<double> &B_tau,int w, vector<double>& input_K, vector<int> Ai,vector<vector<int>>&  Tau){
   int num_B_tau=pow(2,w);
   double temp=0;
   for (int i=0;i<num_B_tau;i++)
     {
       temp=0;
       for (int j=0;j<w;j++)
	 {
	   temp=temp+input_K[Ai[j]]*pow(-1,Tau[i][j]);
	 }
       temp=cosh(temp);
       B_tau.push_back(temp);
     }	 	 
 }
 
void Add_cols(GF2mat& H,GF2mat& D,vector<int>& Ai,const int w,int debug){

  int k=D.rows();
  int r=H.rows();

  for (int i=1;i<w;i++)  // choose 1st, 2nd to w-1 cols, and 0th column to them,
    {
      for (int j=0;j<k;j++)
	{
	  D.set(j,Ai[i],D(j,Ai[i])+D(j,Ai[0]));
	}
    }
  for (int j=0;j<r;j++)
    {
      H.set(j,Ai[0],0);
    }
  
  /*
  if (debug==1)
    {
      GF2mat zero_mat1(H.rows(),D.rows());
      cout<<"after add cols \n H:\n"<<H<<" \n D: \n"<<D<<endl;
      if((H*D.transpose())==zero_mat1) {}
      else{cout<<"after add_cols, H*D^T!=0,"<<endl;return;}
    }
      */
}


void Add_cols_and_rows(GF2mat& H,GF2mat& D,vector<int>& Ai,int w,GF2mat& H2,GF2mat& D2,int debug,vector<vector<int>>& b){

  int num_add_cols=pow(2,w-1)-w; //number of added cols, also is the num of cols of B,cols of B are binray form of 0,1,...2^(w-1)-1,
                             //but delete w-1 weight=1 cols, and 1 weight=0 col, so num_add_col=2^(w-1)-(w-1)-1=2^(w-1)-w
  int c=D.cols();
  int r=D.rows();
  int n=H.cols();
  GF2mat  DB(c,num_add_cols);
  GF2mat  F(num_add_cols,H.cols()+num_add_cols);
  
  construct_BF(D,Ai,w,c,DB,F,b,n,debug);
  
  GF2mat zero_mat(H.rows(),num_add_cols);
  GF2mat H1=merge_mat_hori(H,zero_mat);


  D2=merge_mat_hori(D,DB);
  H2=merge_mat_vert(H1,F);
  
  /*
  if (debug==1)
    {
       GF2mat zero_mat1(H2.rows(),D2.rows());
      if((H2*D2.transpose())==zero_mat1) {}
       else{cout<<"after add_cols_and_rows, H*D^T!=0, that is \n   H\n"<<H2<<"\n D: \n"<<D2<<"\n  result:\n"<<H2*D2.transpose()<<endl;return;}
    }
	*/
      
}

// b stores the binary form of 0...2^(w-1)-1, and the last element is w
void construct_BF(GF2mat& D,vector<int> &Ai, int w,int c, GF2mat& DB, GF2mat& F,vector<vector<int>>& b,int n,int debug){

  int added_c=0;
  int r=D.rows();
  int num_add_cols=pow(2,w-1)-w;
  int num_change_cols=pow(2,w-1)-1; // number of changed cols,which equals to num_add_cols_+w-1
  GF2mat B(c,num_add_cols);
  bvec co(c); // the column added to the right of D
  bvec ro(F.cols());//the row added to the bottom of H
  ro.zeros();
  co.zeros();
  int col_ind=0;  // a column index
  int temp=0;

  //for B:
  // i starts with 1, because i=0 gives the all-zero col.
  for (int i=1;i<=num_change_cols;i++)
    {
      int wt=0;
      int ii=i;
      co.zeros();
      ro.zeros();
      vector<int> bi;  //bi is the binary form of i
      //calculate bi, and where are the corresponding 1-entries of co
      for (int j=1;j<=w-1;j++)
	{
	  //temp is jth digit of i
	  temp=ii%2;
	  ii=ii/2;
	  // cout<<"j is "<<j<<endl;
	  //   cout<<"temp is "<<temp<<endl;
	  //  cout<<"ii is "<<ii<<endl;
	
	  // if jth digit is one, set corresponding entries in co and ro to be one
	  if (temp==1){ co(Ai[j])=1;ro(Ai[j])=1;bi.push_back(1);wt++;}
	  else {bi.push_back(0);};
	}

      //if wt=1, we don't need to add this col
      if (wt>1)
	{
	  //cout<<"col_ind is"<<col_ind<<endl;
	  B.set_col(col_ind,co);
	  if (wt==2){F.set_row(col_ind,ro);}

	  // if wt>2, we can set F(i,i)=F(i,i-1)=1 first and to find  which entry should equals to 1 that keeps the orthogonality.
	  if (wt>2)
	    {
	     
	      //cout<<"B is\n"<<B<<endl;
	      // cout<<"F is\n"<<F<<endl;
	      F.set(col_ind,n+col_ind-1,1);
	      
	      bvec temp_vec1=B.get_col(col_ind);
	      bvec temp_vec2=B.get_col(col_ind-1);
	  
	      // substract the corrsponding columns of B to find where is the third col.
	      bvec temp_vec21=temp_vec1+temp_vec2;
	      int temp_wt=0;
	      int temp_one;
	      //if wt of temp_vec21 =1, then temp_vec1 is the sum of temp_vec2 and a column from D
	      //if wt >1, is the sum of temp_vec2 and a column from DB
	      for (int iii=1;iii<w;iii++)
		{
		  if (temp_vec21(Ai[iii])==1){temp_wt++;temp_one=iii;}
		}
	      if (temp_wt==1){F.set(col_ind,Ai[temp_one],1);}
	      else
		{
		  for (int k=0;k<col_ind-1;k++)
		    {
		      if (temp_vec21==B.get_col(k)){F.set(col_ind,k+n,1);break;}
		      if (k==col_ind-2){cout<<"k==col_ind-2, something went wrong"<<endl;return;}
		    }
		}
	    }
	  //cout<<3333<<endl;
	  F.set(col_ind,n+col_ind,1);
	  col_ind++;
	}
      bi.push_back(w);// store the weight in the last position of bi
      b.push_back(bi);
    }
  DB=D*B;
 
  //if (debug==1){
  // cout<<"B is \n"<<B<<endl;
  // cout<<"F is\n"<<F<<endl;
  // }
  
}


  

/*

bool new_decoder1(GF2mat& H2,GF2mat &H, GF2mat &G,const nodes checks[],const nodes errors[],const vec &pv,const vec& pv_dec,double& num_iter, int lmax,int& max_fail,int&syn_fail,int debug, vec &LR,int rankH,int options)
{
  int v=H.cols();
  int c=H.rows();

  int wt_real_e=0;
  GF2mat real_e(v,1);
  
  wt_real_e=error_channel(real_e, pv);
  if (wt_real_e==0)
    {	  
      return true;
    }
    
  //if no error, break

  GF2mat zero_rvec(1,v);
   
  LR.zeros();
  //vec LR_avg=LR;
  GF2mat zero_mat1(c,1);
  GF2mat zero_mat2(v,1);
  GF2mat zero_rvec2(v-rankH,1);
  
  GF2mat syndrome=H*real_e;

  //is the syndrome a zero vector?
  if (syndrome==zero_mat1)
	{
	  // is e a stablizer?
	  if (G*real_e==zero_rvec2)
	    {
	      return true;
	    }        
	  else  
	    {
	      syn_fail++;
	      return false;	     
	    }  
	}
  
      mat mcv(c,v);   // the messages from check nodes to variable nodes, which are p0/p1
      mat mvc(c,v);   // the messages from variable nodes to check nodes
      mcv.zeros();
      mvc.zeros();
    
      initialize_massages( mcv,mvc, H); //initialize to all-1 matrix

      GF2mat output_e(v,1);
      
      for (int l=1;l<=lmax;l++)
	{
	  //  cout<<l<<endl;
   	  quan_s_update(checks,errors, mcv,mvc,syndrome,pv_dec, c, v,output_e,LR,1);	
	    
	   if (H*output_e==syndrome)
		{		  
		  if(G*(output_e+real_e)==zero_rvec2)
		    {
		      num_iter= num_iter+l;		  
		      return true;
		    }
		  else
		    {
		      // cout<<"\n syndrome is ok, but decoding fails:"<<endl;
		      syn_fail++;		
		      return false;      	        
		    }	    	  
		}  
	}
      
      cout<<"output e is \n"<<endl;
      err_pos2(output_e);
      cout<<"real e is \n"<<endl;
      err_pos2(real_e);      
       
      max_fail++;
      return false;
 }

*/

//if H is Hx, then G is dual of Hz 
bool modified_decoder(GF2mat& H,GF2mat& G,GF2mat &H_tilde, const vector<vector<int>>& A,const nodes checks[],const nodes errors[],const vec &pv,const vec& pv_dec,const int&  num_row_red,double& num_iter, int lmax,int debug, vec &LR,int rankH)
{
  int v=H.cols();
  int c=H.rows();
  int vt=H_tilde.cols();
  int ct=H_tilde.rows();
  if (num_row_red==0){if(H==H_tilde){}else {cout<<"wrong H_tilde"<<endl;}}

  int wt_real_e=0;
  GF2mat real_e(v,1);
  
  wt_real_e=error_channel(real_e, pv);
  if (wt_real_e==0)
    {
      syn_tilde_suc++;
      num_no_error++;
      return true;
    }
    
  //if no error, break

  GF2mat zero_rvec(1,v); 
  LR.zeros();
  
  //vec LR_avg=LR;
  GF2mat zero_mat1(c,1);
  GF2mat zero_mat2(v,1);
  GF2mat zero_mat3(G.rows(),1);
  GF2mat zero_rvec2(v-rankH,1);
  // cout<<11<<endl;
  GF2mat zero_vec_for_syndrome(ct-c,1);
  // cout<<12<<endl;
  
  GF2mat syndrome=H*real_e;

  //is the syndrome a zero vector?
  if (syndrome==zero_mat1)
	{
	  // is e a stablizer?
	  if (G*real_e==zero_rvec2)
	    {
	      syn_tilde_suc++;
	      return true;
	    }        
	  else  
	    {
	      syn_fail++;
	      return false;	     
	    }  
	}
  GF2mat syndrome_tilde(ct,1);
  // cout<<12.125<<endl;
  // cout<<ct-c<<endl;
  if (ct-c!=0){ syndrome_tilde=merge_mat_vert(syndrome,zero_vec_for_syndrome);}
  else{syndrome_tilde=syndrome;}
  // cout<<12.25<<endl;
  
  mat mcv(ct,vt);   // the messages from check nodes to variable nodes, which are p0/p1
  mat mvc(ct,vt);   // the messages from variable nodes to check nodes
  mcv.zeros();
  mvc.zeros();
  
  //  cout<<12.5<<endl;
  initialize_massages( mcv,mvc, H_tilde); //initialize to all-1 matrix
  
  //  cout<<13<<endl;
  
  GF2mat output_e(v,1);
  GF2mat output_e_tilde(vt,1);
      
  for (int l=1;l<=lmax;l++)
    {
      //  cout<<l<<endl;
      quan_s_update(checks,errors, mcv,mvc,syndrome_tilde,pv_dec, ct, vt,output_e_tilde,LR,1);
      
      //  cout<<14<<endl;
      
      if (H_tilde*output_e_tilde==syndrome_tilde)
		{
			inverse_trans(output_e,output_e_tilde,A);
	  
			/*
			if(num_row_red==0){if (output_e==output_e_tilde){}else{cout<<"wrong output_e"<<endl;}}
			if(num_row_red==0){if (syndrome==syndrome_tilde){}else{cout<<"wrong syndrome"<<endl;}}
			if(num_row_red==0){cout<<"H*output_e is \n"<<H*output_e<<"\n H*real_e is \n"<<H*real_e<<endl;}
			*/
			
			syn_tilde_suc++;
	 
			if(G*(output_e+real_e)==zero_mat3)
			{
				num_iter= num_iter+l;		  
				return true;
			}
			else
			{
				syn_fail++;		
				return false;      	        
			}	    	  
		}  
    }
	            
      max_fail++;
      return false;
 }


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

void inverse_trans(GF2mat& output_e,const GF2mat& output_e_tilde,const vector<vector<int>>& A)
{
  // cout<<1<<endl;
  int num_row_red=A.size();  //
  int v=output_e.rows();
  int temp=0;
  vector<int> temp_vec;
  for (int j=0;j<output_e_tilde.rows();j++){temp_vec.push_back(output_e_tilde.get(j,0));}
  //  cout<<2<<endl;

  //here I give an example, assume G has a row (1111), we have an error e=(1011) or (0100) (they are degenerate, since H must have even number of 1s in this row),
  //after a row reduction, e-> (100) no matter which one it is. To convert it back, just add a zero at the position a column was deleted.
  
  for (int i=num_row_red-1;i>=0;i--)
    {
      temp=A[i][0];
      temp_vec.insert(temp_vec.begin()+temp,0);      
    }
  // cout<<3<<endl;
  for (int i=0;i<v;i++){output_e.set(i,0,temp_vec[i]);}

  // cout<<"inverse_trans: num_row_red="<<num_row_red<<endl;
  // if (num_row_red==0){cout<<"output_e is \n"<<output_e<<"\n e_tilde is \n"<<output_e_tilde<<endl;}
}

//here eT is a row vector, because in other functions usually error is a column vector
// something wrong here, so empty_G should always be true
//and in fact, this function is the negative energy
double energy (bool empty_G,const GF2mat& G, const GF2mat& eT, const vector<double> &K)
{

  if (eT.rows()!=1){cout<<"eT.rows!=1, some thing wrong with the energy function"<<endl;return 0;}
  if (empty_G==true)
    {
      // cout<<"e1"<<endl;
      double temp=0;
      for (int i=0;i<eT.cols();i++)
	{
	  int temp2=eT.get(0,i);//becaues get return a bin, convert it to int
	  temp=K[i]*pow(-1,temp2)+temp;
	}
      //cout<<"e2"<<endl;
      return temp;
    }

  
  if(eT.cols()!=G.cols()){cout<<"error in energy func:sizes of e and G do not match"<<endl;return 0;}

  double temp=0;
  GF2mat tempmat=eT*G;
  for (int i=0;i<eT.cols();i++)
    {
      int temp2=tempmat.get(0,i);
      temp=K[i]*pow(-1,temp2)+temp;
    }
  return temp;
}

double ML_suc_rate (vec p,const GF2mat& D,const GF2mat& H, const GF2mat& H_tilde_star, vector<double>&K,int d,int num_large_wt_error)
{
  cout<<"ML_suc_rate start"<<endl;
  
  double suc_rate=0;
  int v=H.cols();
  int c=H.rows();
  
  int vt=H_tilde_star.cols();

  int wt_real_e=0;
  GF2mat real_e(v,1);
  int large_wt_error=0;
  int num_error=0;
  vector<double> K_tilde;

  while (large_wt_error< num_large_wt_error)
    {
      K_tilde.clear();
      num_error++;
      
      int wt_real_e=0;
      GF2mat real_e(v,1); 
      wt_real_e=error_channel(real_e, p);

      // cout<<0<<endl;
      
      if (wt_real_e<d/2){suc_rate++;}   
      else
	{
	  large_wt_error++;
	  GF2mat real_eT=real_e.transpose();
	  
	  GF2mat De=merge_mat_vert(D,real_eT);
	  vector<vector<int>> A;
	  GF2mat H_tilde,e_tilde,empty_tempmat;
	  
	  
	  //  cout<<1<<endl;
	  // here e_tilde is the result transformed matrix De, so its transpose is the error vector transformed
	  Mat_Trans(H,De,K,  K_tilde, H_tilde, e_tilde, A,-D.rows());

	  //	  cout<<2<<endl;
	  if (e_tilde.cols()==H_tilde_star.cols() && e_tilde.rows()==H_tilde_star.rows())
	    {
	      double exp_energy_ratio=exp(energy (true,empty_tempmat, e_tilde, K_tilde))/(exp(energy (true,empty_tempmat, e_tilde, K_tilde))+exp(energy (true,empty_tempmat, e_tilde+H_tilde_star, K_tilde)));
	      //   double exp_e_plus_c_energy=0;
	      
	      if (exp_energy_ratio>0.5){suc_rate++;}
	    }
	  else {cout<<"error: e_tilde.cols()="<<e_tilde.cols()<<" \n H_tilde_star.cols()="<<H_tilde_star.cols()<<"\n  e_tilde.rows()="<<e_tilde.rows()<<"\n H_tilde_star.rows()="<<H_tilde_star.rows()<<endl;}
	}
    }
  cout<<" number of codewords tested: "<<num_error<<endl;

  /*
  cout<<"K_tilde in ML calculating is \n"<<endl;
  for (int ii=0;ii<K_tilde.size();ii++)
    {
      cout<<K_tilde[ii]<<endl;
    }
  */
  
  suc_rate=suc_rate/num_error;
  return suc_rate;

  
}

void ML_decoder_verify(vec p,const GF2mat& H, const GF2mat D,const GF2mat& H_star, const GF2mat D_star,const GF2mat H_tilde, const GF2mat H_tilde_star,vector<double>&K,int max_num_cws,double& after_trans_suc_rate,double&suc_rate,double &after_trans_theoric_suc_rate,double& min_wt_suc_rate,int debug)
{
  // get logical operator:
  GF2mat L;
  // cout<<"ML_suc_rate start"<<endl;
  get_logical(H, D,H_star, L);

  int min_wt_suc_count=0;
  int after_trans_suc_count=0;
  int suc_count=0;
  int after_trans_theoric_suc_count=0;

  int num_cws_count=0;
  
   // bool after_trans_suc=true;
   // bool suc=true;
  
  int v=H.cols();
  int c=H.rows();
  
  int vt=H_tilde_star.cols();

  vector<double> K_tilde;

  //generate a random error
  while (num_cws_count<max_num_cws)
    {
      num_cws_count++;
      K_tilde.clear();
      
      int wt_real_e=0;
      GF2mat real_e(v,1); 
      wt_real_e=error_channel(real_e, p);

      if (debug==1){ cout<<"real_e is \n"<<real_e<<endl;}
 

      // cout<<0<<endl;
      
      GF2mat real_eT=real_e.transpose();
      
      GF2mat De=merge_mat_vert(D,real_eT);
      vector<vector<int>> A;
      GF2mat H_tilde,e_tilde,empty_tempmat;
	    
      // cout<<3<<endl;
      //here e_tilde is a row vector so it is in fact the transposed one
      Mat_Trans(H,De,K,  K_tilde, H_tilde, e_tilde, A,-D.rows());
           if (debug==1)
	     {
	       cout<<"after trans real_e is \n"<<e_tilde<<endl;
	       cout<<"e_tilde+logical_operattor is\n"<<e_tilde+H_tilde_star<<endl;
	     }

      //  cout<<4<<endl;

      //theoric ML suc rate:
      if (e_tilde.cols()==H_tilde_star.cols() && e_tilde.rows()==H_tilde_star.rows())
	{
	  double exp_energy_ratio=exp(energy (true,empty_tempmat, e_tilde, K_tilde))/(exp(energy (true,empty_tempmat, e_tilde, K_tilde))+exp(energy (true,empty_tempmat, e_tilde+H_tilde_star, K_tilde)));
	  //   double exp_e_plus_c_energy=0;
	      
	  if (exp_energy_ratio>=0.5)
	    {
	      after_trans_theoric_suc_count++;
	           if (debug==1)
		     {
		       cout<<"in theory: suc"<<endl;
		     }
	    }
	}
	  else {cout<<"error: e_tilde.cols()="<<e_tilde.cols()<<" \n H_tilde_star.cols()="<<H_tilde_star.cols()<<"\n  e_tilde.rows()="<<e_tilde.rows()<<"\n H_tilde_star.rows()="<<H_tilde_star.rows()<<endl;}


      //ML decoder:
      GF2mat syndrome=H*real_e;
      GF2mat algebraic_e(H.cols(),1);
      GF2mat ML_output_e(H.cols(),1);
      //  cout<<5<<endl;
      
      algebraic_decoder(H,syndrome,algebraic_e);
     
      
      // cout<<6<<endl;
      //  cout<<"H is \n"<<H<<" \n H_star is \n"<<H_star
      if (debug==1)  {cout<<"ML decoder:\n algebraic_e is\n"<<algebraic_e<<endl; }
      // cout<<"ML decoder:\n algebraic_e is\n"<<algebraic_e<<endl;
 
      // cout<<"algebraic_e +L is\n"<<algebraic_e+L.transpose()<<endl;

      int alge_e_wt=s_weight(algebraic_e);
      int alge_e_L_wt=s_weight(algebraic_e+L.transpose());
      GF2mat min_wt_e;
      
      if (alge_e_wt<=alge_e_L_wt) {min_wt_e=algebraic_e;}
      
      else {min_wt_e=algebraic_e+L.transpose();}
      
      GF2mat tempzeromat(D_star.rows(),1);
      if (H*min_wt_e==syndrome)
	{
	  if (D_star*(min_wt_e+real_e)==tempzeromat)
	    {
	      min_wt_suc_count++;
	    }
	}
      
	  ML_decoder(algebraic_e, D,L,K,ML_output_e);
      
           if (debug==1)
	     {
	       cout<<"ML decoder: output_e is\n"<<ML_output_e<<endl;
	     }

      //  cout<<7<<endl;
    
      if (H*ML_output_e==syndrome)
	{
	  //  cout<<"D_star is\n"<<D_star<<endl;
	  // cout<<"sum is \n"<<ML_output_e+real_e<<endl;
	  //	  cout<<"D_star*(ML_output_e+real_e)=\n"<<D_star*(ML_output_e+real_e)<<endl;
	  //	  cout<<"D_star*(ML_output_e+real_e+Logical)=\n"<<D_star*(ML_output_e+real_e+L.transpose())<<endl;
	  //cout<<"zeromat is \n"<<one_x_n_zeromat<<endl;
	  if (D_star*(ML_output_e+real_e)==tempzeromat)
	    {
	      suc_count++;
	      if (debug==1){ cout<<"ML suc"<<endl;}
	    }
	}
      else {cout<<"sth wrong with the ML decoder, H*output_e!=syndrome"<<endl;}
	
    

      //after tansformation ML decoder:
      GF2mat syndrome_tilde(H_tilde.rows(),1);
      GF2mat zero_vec_for_syndrome(H_tilde.rows()-H.rows(),1);
      syndrome_tilde=merge_mat_vert(syndrome,zero_vec_for_syndrome);
      
      GF2mat algebraic_e_tilde(H_tilde.cols(),1);
      GF2mat	ML_output_e_tilde(H_tilde.cols(),1);
      GF2mat wrong_e(H_tilde.cols(),1);
      //  cout<<8<<endl;
      algebraic_decoder(H_tilde,syndrome_tilde,algebraic_e_tilde);
           if (debug==1)
	     {
	       cout<<"after trans algebraic_e: \n"<<algebraic_e_tilde<<endl;
	       cout<<"after trans algebraic_e+logical operator: \n"<<algebraic_e_tilde+H_tilde_star.transpose()<<endl;
	     }
  
      
      //  cout<<9<<endl;
      double exp_energy_ratio_trans=exp(energy (true,empty_tempmat, algebraic_e_tilde.transpose(), K_tilde))/(exp(energy (true,empty_tempmat, algebraic_e_tilde.transpose(), K_tilde))+exp(energy (true,empty_tempmat, algebraic_e_tilde.transpose()+H_tilde_star, K_tilde)));
      //   double exp_e_plus_c_energy=0;
      //  cout<<10<<endl;
      if (exp_energy_ratio_trans>=0.5){ML_output_e_tilde=algebraic_e_tilde;wrong_e=algebraic_e_tilde+H_tilde_star.transpose();} else{ML_output_e_tilde=algebraic_e_tilde+H_tilde_star.transpose();wrong_e=algebraic_e_tilde;}
      //x cout<<"after trans ML_output_e_tilde=\n"<<ML_output_e_tilde<<endl;
      //  cout<<11<<endl;
      GF2mat trans_ML_output_e(H.cols(),1);
      GF2mat trans_ML_wrong_output_e(H.cols(),1);
      
      inverse_trans(trans_ML_output_e,ML_output_e_tilde,A);
      inverse_trans(trans_ML_wrong_output_e,wrong_e,A);
           if (debug==1)
	     {
	       cout<<"after inverse trans :\n"<<trans_ML_output_e<<endl;
	       cout<<"after inverse trans wrong e:\n"<<trans_ML_wrong_output_e<<endl;
	     }
      //cout<<12<<endl;
      //  cout<<H<<endl;
      //  cout<<trans_ML_output_e<<endl;
      if (H*trans_ML_output_e==syndrome)
	{
	  // cout<<13<<endl;
	  if (D_star*(trans_ML_output_e+real_e)==tempzeromat){after_trans_suc_count++;}
	  // else  {cout<<"(D_star*(trans_ML_wrong_output_e+real_e)=\n"<<D_star*(trans_ML_wrong_output_e+real_e)<<endl;}
	  // cout<<14<<endl;
	}
      else {cout<<"sth wrong with the ML decoder, H*output_e!=syndrome"<<endl;}
      
    }

  after_trans_theoric_suc_rate=1.0*after_trans_theoric_suc_count/max_num_cws;
  after_trans_suc_rate=1.0*after_trans_suc_count/max_num_cws;
  suc_rate=1.0*suc_count/max_num_cws;
  min_wt_suc_rate=1.0*min_wt_suc_count/max_num_cws;
}



//it is better to use L instead of H_star, but I havenot wrotten the function to calculate L yet.
void ML_decoder(const GF2mat& input_e, const GF2mat& D,const GF2mat& L,vector<double>& K,GF2mat& output_e){
  
  if(D.rows()>10){cout<<" D>rows()>10, it is better not to use ML decoder"<<endl;return;}
  
  double unnormalized_p=0;
 
  
  GF2mat input_e_T=input_e.transpose();
  GF2mat temp_output_e_T=input_e_T;
  double max_unnormalized_p=exp(energy (true,D, input_e_T, K));

  //for input_e:
  //start with i=1 because alpha=0 is already added
  for (int i=1;i<pow(2,D.rows()-1);i++)  
    {
      int tempi=i;
      GF2mat alpha(1,D.rows());

      // get an alpha
      for (int l=0;l<D.rows();l++)
	{
	  int tempi_red=tempi%2;
	  tempi=tempi/2;
	  alpha.set(0,l,tempi_red);
	}
      
      GF2mat tempmat1= input_e_T+alpha*D;

      max_unnormalized_p=max_unnormalized_p+exp(energy (true,D, tempmat1, K));
    }

  //calculate sum of p(e+alpha*D) over alpha for all possible e
  //go over possible e
  for (int j=0;j<L.rows();j++)
    {
      //  cout<<"j="<<j<<endl;
      //  cout<<"rows="<<H_star.rows()<<endl;
      GF2mat tempmat=L.get_row(j);
      tempmat=tempmat.transpose();

      GF2mat temp_e_T=input_e_T+tempmat;
      unnormalized_p=exp(energy (true,D, temp_e_T, K));
   
      //go over alpha
      for (int i=1;i<pow(2,D.rows()-1);i++)  
	{
	  //  cout<<"i="<<i<<endl;
	  int tempi=i;
	  GF2mat alpha(1,D.rows());

	  // get an alpha
	  for (int l=0;l<D.rows();l++)
	    {
	      int tempi_red=tempi%2;
	      tempi=tempi/2;
	      //  cout<<62<<endl;
	      alpha.set(0,l,tempi_red);
	    }
	  //	  cout<<63<<endl;
	  // cout<<temp_e_T<<endl;
	  //  cout<<D<<endl;
	  //  cout<<alpha<<endl;
	  GF2mat tempmat1= temp_e_T+alpha*D;
	  //	  cout<<64<<endl;
	  unnormalized_p=unnormalized_p+exp(energy (true,D, tempmat1, K));
	  //	  cout<<65<<endl;
	}
      if (unnormalized_p>max_unnormalized_p){max_unnormalized_p=unnormalized_p;temp_output_e_T=temp_e_T;}
      //  cout<<66<<endl;
    }
  //  cout<<67<<endl;
  output_e=temp_output_e_T.transpose();
  //  cout<<68<<endl;

}
  
//GF2mat Logical_O()

//give a output_e satisfy H*output_e=s
void algebraic_decoder(const GF2mat& H,const GF2mat& syndrome,GF2mat &output_e) 
{ 
 
 //cout<<"as"<<endl;
 
  int n=H.cols();
  int r=H.rows();

  GF2mat H1=H;
   
  GF2mat T;
  ivec perm2;
  GF2mat U;
  // now we have T*H1*perm2*perm2_inv*e=T*s, def s1=T*s, U=T*H1*perm2, e1=perm2_inv*e: U*e1=s1

  int rankH=H1.T_fact(T,U,perm2);
  if (rankH==0){return;}
  GF2mat perm2_mat=col_permutation_matrix(perm2);
  GF2mat perm2_mat_inv=perm2_mat.inverse();
 
  GF2mat syndrome1=T*syndrome;
 
    //U may be not a full rank matrix, need to delete some rows. So also need to delete some rows of syndrome1
	//cout<<H2<<endl;
	//cout<<"rank is"<<rankH<<"H1 is\n"<<H1<<endl;
  GF2mat H2=U.get_submatrix(0,0,rankH-1,n-1);


 //  cout<<"H2 is\n"<<H2<<endl;
  GF2mat syndrome2=syndrome1.get_submatrix(0,0,rankH-1,0);
  // cout<<"syndrome2 is\n"<<syndrome2<<endl;
  
  
  // OSD-0: e1=(HS_inv*s2
  //                  0   )
  // cout<<51<<endl;
  GF2mat H_S=H2.get_submatrix(0,0,rankH-1,rankH-1);
 //  cout<<"H_S is \n"<<H_S<<endl;
   GF2mat H_T(1,1);
   
  if (H2.cols()>rankH){H_T=H2.get_submatrix(0,rankH,rankH-1,n-1);}
  // cout<<52<<endl;
  GF2mat HS_inv=H_S.inverse();
  // cout<<"HS_inv is\n"<<HS_inv<<endl;
  
  GF2mat e_S=HS_inv*syndrome2;
    
	//cout<<"e_S is\n"<<e_S<<endl;
 
 GF2mat e_T(1,1);
  if(n-rankH>0){ GF2mat e_T_temp(n-rankH,1); e_T=e_T_temp;}

//cout<<3<<endl;

  // cout<<53<<endl;
  for (int i=0;i<rankH;i++){output_e.set(i,0,e_S(i,0));}
  //  cout<<54<<endl;
  if(n-rankH>0) {for (int i=rankH;i<n;i++){output_e.set(i,0,e_T(i-rankH,0));}}
  // cout<<55<<endl;
  output_e=perm2_mat*output_e;
 // cout<<"output_e is \n"<<output_e<<endl;
   	  	  
  if (H*output_e==syndrome){}
  else
    {
		cout<<"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"<<endl;
		cout<<"algebraic_decoder error! \n rankH is "<<rankH<<endl;
      cout<<"H*output_e!=syndrom H is\n"<<H<<"\n output_e is\n"<<output_e.transpose()<<"\n syndrome is \n"<<syndrome<<endl;
	  cout<<"H2 is \n"<<H2<<"\n HS is\n"<<H_S<<"\n H_t is \n"<<H_T<<"\n e_S is \n"<<e_S.transpose()<<"\n e_T is\n "<<e_T.transpose()<<"\n syndrome2 is\n"<<syndrome2<<endl;
    }  
}

//G is dual of Hx, Gx is dual of Hz.
bool   alternative_decoder(const GF2mat &Hx, const GF2mat &Hz,const GF2mat &Gx,GF2mat G,const GF2mat &real_e,GF2mat &final_output_e,const nodes checks1[],const nodes errors1[],const vec &pv, vec&pv_dec,const double& pmin,const double& pmax,double& num_iter, int lmax,int &wt,  int debug, vec &LR2,int rankH,double M,int&LLR_fail,int&LLR_suc,int& num_fixed_LLR,bool& start_row_reduction,double alpha,int other_decoder,int channel,int max_row_freeze, int max_rep)
{
	// variables with prefix const_ are the variables before any row reduction

 int c0=Hx.rows();
 int v0=Hx.cols();
  GF2mat H=Hx;
  		int v=H.cols();
		int c=H.rows();
		

if (other_decoder==1 and M>=0){start_row_reduction=true;}
else{start_row_reduction=false;}
  
  
  GF2mat real_e_copy=real_e;
  //if no error, break

  GF2mat zero_rvec(1,v);
   
  //LR.zeros();
  //vec LR_avg=LR;
  GF2mat zero_mat1(c,1);
  GF2mat zero_mat2(v,1);
  GF2mat tempmat333=Gx*real_e;
  GF2mat zero_rvec2(tempmat333.rows(),1);
   GF2mat const_zero_rvec2(tempmat333.rows(),1);
  
  GF2mat syndrome=Hx*real_e;
  GF2mat original_syndrome=syndrome;
  

  //is the syndrome a zero vector?
  if (syndrome==zero_mat1)
	{
	  // is e a stablizer?
	  if (Gx*real_e==zero_rvec2)
	    {
			num_no_error++;
			final_output_e=real_e;
	      return true;
	    }        
	  else  
	    {
	      syn_fail++;
	      return false;	     
	    }  
	}
	
      mat mcv(c,v);   // the messages from check nodes to variable nodes, which are p0/p1
      mat mvc(c,v);   // the messages from variable nodes to check nodes
      mcv.zeros();
      mvc.zeros();
    
      initialize_massages( mcv,mvc, H); //initialize to all-1 matrix

      
	  GF2mat res_output_e(v,1);
	  GF2mat const_output_e(v,1);
    
      vector<GF2mat>rows_del;// each element stores a row of G deleted 
	  vector<int>value_deleted_entries;// values deleted 
	  
	  vector<vector<int>> pos_deleted_entries;  // each element stores the columns and rows deleted/changed in a round of row reduction: 
	  //example of an element: [[4,8,-1,2,3],[5,10,-2,7,9]] , means in this round, 4th and 5th columns of H are deleted, 
	  //8th row of G is added to 2h, 3th rows, 10th column is added to 7th,9th rows.
	  //-1 means this row is also added to e, -2 means not.
      int E1=0;
	  
	vec old_LR(v);
	
	vec LR_avg(v);
	vec LR_prod(v);
	for (int i=0;i<v;i++){LR_prod(i)=1;}
	
	
	vec LR_final;
	bool empty_G=false;
	bool converge_again=false; //is BP converged again after row reduction?
	
	//if ((debug/4)%2==1) {cout<<"real_e is\n"<<endl;err_pos2(real_e);}
	//bool start_row_reduction=false;
	//int ltemp=0;
	//int ltemp2=0;
	//int ltemp3=0;
	//int ltemp4=0;
	int l2=0;
	
	if (M>0){max_rep=1;}
	if (debug>8){max_rep=1;}
	GF2mat zero_e(v,1);
	final_output_e=zero_e;
	
	while (l2<max_rep)
	{
		//pro_dist(pmin,pmax,pv_dec);
			
		mcv.zeros();
		mvc.zeros();
    
		initialize_massages( mcv,mvc, H); //initialize to all-1 matrix
	
		GF2mat output_e(v,1);
		
		
      for (int l=1;l<=lmax;l++)
	{
		v=H.cols();
		c=H.rows();
	  //cout<<l<<"th iteration starts"<<endl;
	   nodes* checks = new nodes[H.rows()];
	  nodes* errors = new nodes[H.cols()];  
	  
	   //std::unique_ptr<nodes[]> checks(new nodes[H.rows()]);
   // std::unique_ptr<nodes[]> errors(new nodes[H.cols()]);
	
	  initialize_checks (H, checks,  E1);
	  initialize_errors(H, errors);
	 vec LR(v);
	 
			if (debug&2)
			{
				quan_p_update(checks,errors, mcv,mvc,syndrome,pv_dec, c, v,output_e,LR,debug,alpha);				
			}
	      else
			{
				//cout<<"serial V starts"<<endl;			
				quan_s_update(checks,errors, mcv,mvc,syndrome,pv_dec, c, v,output_e,LR,debug,alpha);
			}
	    
	    
		//cout<<"l="<<l<<" before check after updating H=\n"<<H<<"\noutput_e=\n"<<output_e<<"\nsyndrome=\n"<<syndrome<<endl;
		if (debug&4 and start_row_reduction==true)
		{
			cout<<"l="<<l<<"\n LLR=\n"<<endl;
			for (int aaa=0;aaa<LR.size();aaa++){cout<<log(LR(aaa))<<"  ";}
			cout<<endl;
			cout<<"output_e is\n"<<output_e.transpose()<<endl;
			right_restore_e(res_output_e,output_e,  pos_deleted_entries,value_deleted_entries,rows_del);
			cout<<"after restore"<<endl;
			err_pos2(res_output_e);
			cout<<"H is\n"<<H<<"\n G is\n"<<G<<"\nsyndrome is\n"<<syndrome.transpose()<<endl;
		}
	
			if (debug&32)
			{
					for (int iii=0;iii<v;iii++)
				{
					LR_prod(iii)=pow(LR_prod(iii),0.9)*LR(iii);
					LR_avg(iii)=pow(LR_prod(iii),1.0/(1-pow(0.9,l)));
					output_e.set(iii,0,LR_avg(iii)<1? 1:0);   
				}
			}
			
	   if (H*output_e==syndrome)
		   
		{		  final_output_e=final_output_e+output_e;
				//cout<<"check1 is ok"<<endl;
				if  (H.cols()!=v0)
				{
				 right_restore_e(res_output_e,output_e,  pos_deleted_entries,value_deleted_entries,rows_del);
				}
				else{res_output_e=final_output_e;}
				
				  if(Gx*(res_output_e+real_e)==const_zero_rvec2)
					{
						num_iter= num_iter+l;		
						if(start_row_reduction==true){fix_suc++;}			
							if (checks !=NULL){delete[] checks;}
							if(errors!=NULL) {delete[] errors;}
							BP_suc++;
						
							if (debug&32) {avg_BP_suc++;}
						final_output_e=res_output_e;
						return true;
					}
		    
		  else
		    {
		      // cout<<"\n syndrome is ok, but decoding fails:"<<endl;
			  if(start_row_reduction==true){fix_fail++;}			
		      syn_fail++;		
			  	if (checks !=NULL){delete[] checks;}
		if(errors!=NULL) {delete[] errors;}
		      return false;      	        
		    }	    	  
		}
		

	
	//cout<<"start row reduction0"<<endl;
		//where to perform reduction?
		//cout<<"G is \n"<<G<<"\n H is \n"<<H<<endl;
		if ( start_row_reduction==false and (LR==old_LR or l>0.5*lmax) and M>=0)
		
		{
			//cout<<"M="<<M<<endl;
			start_row_reduction=true;
			if (l<0.5*lmax) {converge_to_trapping_set++;}
			else {did_not_converge_for_large_iterations++;}
			
			if (debug&4 )
			{
			cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
			if (l<0.5*lmax) {cout<<" \n l="<<l<<"  LLR converged, start row reduction: "<<endl;}
			else {cout<<" \n l="<<l<<"  but LLR did not converge, perhaps meets an oscillation or maybe it's a matter of calculation accuracy, start row reduction: "<<endl;}
			cout<<"real_e is"<<real_e.transpose()<<endl;
			err_pos2(real_e);
			cout<<"output_e is\n"<<output_e.transpose()<<endl;		
			err_pos2(output_e);
			
			cout<<"\n LLR=\n"<<endl;
			for (int aaa=0;aaa<LR.size();aaa++){cout<<log(LR(aaa))<<"  ";}
			cout<<endl;
			cout<<"H is\n"<<H<<"\n G is\n"<<G<<"\nsyndrome is\n"<<syndrome.transpose()<<endl;
			}
			
		}
		
		if (  (not empty_G ) and  start_row_reduction and M>=0 )
		{
			//cout<<"del some cols of H"<<endl;
			int num_del_col=0;
			int hardSelectValue=-1;
			bool exists_larger_than_M=true;
			vec abs_LLR(LR.size());
			for (int jjj=0;jjj<LR.size();jjj++)
			{
				abs_LLR(jjj)=abs(log(LR(jjj)));
			}
		    
			while( exists_larger_than_M)
			{
				ivec sortindex=sort_index(abs_LLR);// ascending order
				int max_ind=sortindex(sortindex.size()-1);	
						
				double max_abs_LLR=abs_LLR(max_ind);
				
					if (max_abs_LLR>M)
						{			 
							num_fixed_LLR++;
							if (not empty_G)
							{
								right_small_row_reduction(empty_G,syndrome,H,G,output_e,max_ind, pos_deleted_entries,mcv,mvc,rows_del,value_deleted_entries, hardSelectValue,LLR_fail,LLR_suc,debug);	
								abs_LLR.del(max_ind);
							}
							else{break;}
					
							num_del_col++;
						}
					else {exists_larger_than_M=false;}
			}
					if (debug&4){cout<<"delete "<<num_del_col<<" columns in this round"<<endl;}			
					
					if (not empty_G)
					{
						if (GF2mat_rank(H)+GF2mat_rank(G)!=H.cols()){cout<<"error: after row reduction GF2mat_rank(H)+GF2mat_rank(G)!=H.cols(): H is\n"<<H<<"\n G is \n"<<G<<endl;}
					}
		 }
			
		  //protocol B:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111
			//if ( debug&16 and l==lmax and (start_row_reduction or M==-50))
				//if (l==lmax) {cout<<"should goto pb"<<endl;}
			//ltemp2=l;
			if ( (debug&16) and l==lmax and (start_row_reduction or M==-10))
			{
				    //cout<<"P_B"<<endl;
					algebraic_decoder(H, syndrome,output_e);
					
					right_restore_e(res_output_e,output_e,  pos_deleted_entries,value_deleted_entries,rows_del);
					
					if (debug&4 )
					{
						cout<<"meet max iterations, start protocol B,  which output the error (after restoring):\n"<<endl;
						err_pos2(res_output_e);
					}
					
					 if (Hx*res_output_e==original_syndrome)
					 {
						if (Gx*(res_output_e+real_e)==const_zero_rvec2)
						{				 
							//num_iter= num_iter+l;		
							fix_suc++;
							protocol_B_suc++;
								if (checks !=NULL){delete[] checks;}
		if(errors!=NULL) {delete[] errors;}
		final_output_e=res_output_e;
							return true;					
						}
						else
						{
						syn_fail++;
						protocol_B_syn_fail++;
						fix_fail++;
							if (checks !=NULL){delete[] checks;}
		if(errors!=NULL) {delete[] errors;}
						return false;				 
						}
					 }
					 else
					 {
						// cout<<"PB max fail"<<endl;
						protocol_B_max_fail++;
						fix_fail++;
						max_fail++;
							if (checks !=NULL){delete[] checks;}
		if(errors!=NULL) {delete[] errors;}
						return false;					 
					 }
						
			}		 
		
		
		//cout<<"l="<<l<<"  over , go to next iteration "<<"G is\n"<<G<<"\nH is\n"<<H<<endl;
		if (old_LR.size()==LR.size()){old_LR=LR;}

		if (checks !=NULL){delete[] checks;}
		if(errors!=NULL) {delete[] errors;}
		//ltemp=l;
		LR_final=LR;
	}
	l2++;
	final_output_e=output_e+final_output_e;
	syndrome=syndrome+H*final_output_e;
	}
	syndrome=original_syndrome;
      
		if (max_row_freeze>0)
		{
			if (max_row_freeze>Hz.rows()){max_row_freeze=Hz.rows();}
			bool rel_suc=false;
			ivec R_index=reliability_sort_index ( Hz, LR_final);      
			for (int freeze_ind=0;freeze_ind<max_row_freeze;freeze_ind++)
			{
				GF2mat output_H,output_col_perm,output_row_perm,H2,A,output_e1, output_H1,output_e;
				//cout<<1<<endl;
				reliability_perm_H (Hx, R_index(freeze_ind),output_H,output_col_perm, output_row_perm, H2, A,output_H1);
				GF2mat output_e2(H2.cols(),1);
				//cout<<2<<endl;
				GF2mat rel_syndrome =output_row_perm*syndrome;
				GF2mat syndrome_2=rel_syndrome.get_submatrix(A.rows(),0,syndrome.rows()-1,0);
				GF2mat syndrome_1=rel_syndrome.get_submatrix(0,0,A.rows()-1,0);
				//cout<<3<<endl;
				
				     int  v2=H2.cols();
					int	c2=H2.rows();
	  
		                nodes* checks2 = new nodes[H2.rows()];
	                     nodes* errors2 = new nodes[H2.cols()];  
	  
						initialize_checks (H2, checks2,  E1);
						initialize_errors(H2, errors2);
						vec LR2(v2);
				 mat mcv2(c2,v2);   
					mat mvc2(c2,v2);   
						mcv2.zeros();
						mvc2.zeros();
                
				vec pv_dec2(v2);
				for (int temp_ind=0;temp_ind<v2;temp_ind++){pv_dec2(temp_ind)=pv_dec(temp_ind);}
				
			    initialize_massages( mcv2,mvc2, H2); //initialize to all-1 matrix
				//cout<<4<<endl;
				for (int l=1;l<=lmax;l++)
					{
							 
						if (debug&2)
							{
									quan_p_update(checks2,errors2, mcv2,mvc2,rel_syndrome,pv_dec2, c2, v2,output_e2,LR2,debug,alpha);				
							}
						else
							{	
								quan_s_update(checks2,errors2, mcv2,mvc2,rel_syndrome,pv_dec2, c2, v2,output_e2,LR2,debug,alpha);
							}
	    
			
							if (H2*output_e2==rel_syndrome){break;}
					}
					//cout<<5<<endl;
					rel_suc=reliability_solve(output_H1,A, syndrome_1, output_e2,output_e1);
					//cout<<output_H1<<endl;
					//cout<<6<<endl;
					if (rel_suc==true)
					{
						//cout<<7<<endl;
						//cout<<H<<endl;
						//cout<<output_e1<<endl;
						//cout<<output_e2<<endl;
						//cout<<output_col_perm<<endl;
						output_e=reliability_restored_e( output_col_perm, output_e1,output_e2);
						//cout<<8<<endl;
						  if(Gx*(output_e+real_e)==const_zero_rvec2)
						  {
							  freeze_rows_suc++;
							  return true;
						  }
					}
					
					if (checks2 !=NULL){delete[] checks2;}
		if(errors2!=NULL) {delete[] errors2;}
			}
		}
		   
								
		
		
	
      if (debug&8)
	{
		int suc_order=0;
		avg_red_qubits=avg_red_qubits+H.cols();
		
		if (debug&32){LR_final=LR_avg;}
	  //cout<<"OSD:M="<<M<<endl;
	  //full_OSD(LR_final,H,syndrome,output_e,2, suc_order);
	  GF2mat output_e(v,1);
	  OSD(LR_final,H,syndrome,output_e,suc_order);

	  if (start_row_reduction==true)
	  {
		right_restore_e(res_output_e,output_e,  pos_deleted_entries,value_deleted_entries,rows_del);
	  }
	  else {res_output_e=output_e;}
	  /*
	  	if (debug&4 )
		{
			cout<<"OSD output the error (after restoring):\n"<<endl;
			err_pos2(res_output_e);
		}
		*/
	  if(Gx*(res_output_e+real_e)==const_zero_rvec2)
	    {
	      // cout<<"OSD suc2.1"<<endl;
	      OSD_suc++;
		  if(start_row_reduction==true){fix_suc++;}	
		  if (suc_order==0){OSD_0_suc++;}
		else if(suc_order==1){OSD_1_suc++;}
		else if(suc_order==2){OSD_2_suc++;}
		else {OSD_higher_order_suc++;}
		final_output_e=res_output_e;
	      return true;
	    }	 
		else 
		{
			 
	
	  
			
			if (debug&4)
	 {
		 cout<<"\n\n\nxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
		 cout<<"OSD fails:"<<endl;
	   cout<<"output e is \n"<<endl;
	   err_pos2(res_output_e);
	   
	   cout<<"OSD-0 e_S="<<global_e_S.transpose()<<endl;
	  // cout<<"new_wt="<<global_new_wt<<endl;
	   
	   cout<<"real e is \n"<<endl;
	   err_pos2(real_e);
	   GF2mat sume=real_e+res_output_e;
	   cout<<"residual e is \n"<<endl;
	   err_pos2(sume);
	  // cout<<"H*residual_e="<<(Hx*sume).transpose()<<endl;
	   //  cout<<"mcv is \n"<<mcv<<endl;
	   //  cout<<"\n mvc is \n"<<mvc<<endl;
	   cout<<"LR is \n"<<endl;
	   int dist=sqrt(H.cols());
	   for (int ii=0;ii<LR_final.length();ii++)
	   {
		   cout<<LR_final(ii)<<"   ";
		   if (ii%dist==dist-1){cout<<endl;}
	   }
	   
	   //cout<<"sort index is (ascending)  "<<global_sort_index<<endl;
	   
	  // cout<<"H is\n"<<H<<endl;
	  // GF2mat H1=H;
     //  H1.permute_cols(global_sort_index,0);
	 //  cout<<"after perm1\n"<<H1<<endl;
	   cout<<"perm2 is\n"<<global_perm2<<endl;
	//    GF2mat T;

  //GF2mat U;


 // int xxx=H1.T_fact(T,U,global_perm2);
	  // cout<<"after perm2\n"<<U<<endl;
	 }
			OSD_fail++;
			if(start_row_reduction==true){fix_fail++;}	
			syn_fail++;
			
			return false;
		}
	}
 
	 if(start_row_reduction){fix_fail++;}		
//cout<<"end max fail xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;	 
//cout<<"ltemp="<<ltemp<<endl;

  if (debug&4)
	 {
	   cout<<"output e is \n"<<endl;
	   err_pos2(res_output_e);
	   cout<<"real e is \n"<<endl;
	   err_pos2(real_e);
	   GF2mat sume=real_e+res_output_e;
	   cout<<"residual e is \n"<<endl;
	   err_pos2(sume);
	   //  cout<<"mcv is \n"<<mcv<<endl;
	   //  cout<<"\n mvc is \n"<<mvc<<endl;
	   cout<<"LR is \n"<<endl;
	   int dist=sqrt(H.cols());
	   for (int ii=0;ii<LR_final.length();ii++)
	   {
		   cout<<LR_final(ii)<<"     ";
		   if (ii%dist==dist-1){cout<<endl;}
	   }
	   
	 }
      max_fail++;
      return false;
 }
	
void del_col (GF2mat & H, int col_ind)
{
	if (col_ind<0 or col_ind>H.cols()-1) {cout<<"del_col: wrong col_ind H.cols()="<<H.cols()<<"  col_ind="<<col_ind<<endl;return;}
	if(H.cols()==1){cout<<"del_col the matrix only have one col, cannot delete it (GF2mat cannot creat a matrix with no entires"<<endl;return;}
	
	if (col_ind==H.cols()-1){H=H.get_submatrix(0,0,H.rows()-1,H.cols()-2);}
	
	else if (col_ind==0){ H=H.get_submatrix(0,1,H.rows()-1,H.cols()-1);}
	else
	{
	
		GF2mat H1=H.get_submatrix(0,0,H.rows()-1,col_ind-1);
		GF2mat H2=H.get_submatrix(0,col_ind+1,H.rows()-1,H.cols()-1);
		H=merge_mat_hori(H1,H2);
	}

}

	/*
void del_col (mat & H, int col_ind)
{
	if (col_ind<0 or col_ind>H.cols()-1) {cout<<"del_col: wrong col_ind H.cols()="<<H.cols()<<"  col_ind="<<col_ind<<endl;return;}
	
	if (col_ind==H.cols()-1){H=H.get_submatrix(0,0,H.rows()-1,H.cols()-2);}
	else if (col_ind==0){H=H.get_submatrix(0,1,H.rows()-1,H.cols()-1);}
	else
	{
		mat H1=H.get_submatrix(0,0,H.rows()-1,col_ind-1);
		mat H2=H.get_submatrix(0,col_ind+1,H.rows()-1,H.cols()-1);
		H=merge_mat_hori(H1,H2);
	}
}
*/

//A store the position where output_e's entry are deleted (output_e are changed by each deletion, so its size is changed by each deletion), B stores the elements deleted.
//note: output_e is a column vector
void restore_e(GF2mat& output_e,  const GF2mat& original_output_e, const vector<vector<vector<int>>> &pos_deleted_entries,const vector<vector<int>>& value_deleted_entries)
{
	
	if (pos_deleted_entries.size()!=value_deleted_entries.size()){cout<<"restore_e error: pos_deleted_entries.size="<<pos_deleted_entries.size()<<"  value_deleted_entries.size="<<value_deleted_entries.size()<<endl;return;}
	output_e=original_output_e;
	if (output_e.cols()!=1){cout<<"restore_e error: output_e.cols()!=1"<<endl;return;}	
//	cout<<"r1"<<endl;
	for (int ii=value_deleted_entries.size()-1;ii>=0;ii--)
	{
		vector<int> B=value_deleted_entries[ii];
		vector<vector<int>> Ai=pos_deleted_entries[ii];
		if (Ai.size()!=B.size())
		{
			//cout<<"restore_e error: Ai.size!=B.size"<<endl;
			//cout<<"Ai is \n"<<Ai.size()<<"\n B is\n"<<B.size()<<endl;
			return;
			
			}
//	cout<<"r2"<<endl;
	
		for (int j=Ai.size()-1;j>=0;j--)
		{
			GF2mat value_to_be_res(1,1);
			value_to_be_res.set(0,0,B[j]);
		
			//the first element of Ai[j] is the postion deleted in jth round
			int row_to_be_res=Ai[j][0];
		
		
			//restore the col deleted
		//	cout<<"insert row st"<<endl;
			insert_row(output_e,row_to_be_res,value_to_be_res);	
		//	cout<<"insert row suc"<<endl;
		
			// restore the other cols 
			for (int i=Ai[j].size()-1;i>=1;i--)
			{
				output_e.set(Ai[j][i], 0,output_e(Ai[j][i],0)+output_e(row_to_be_res,0));
			}
		//	cout<<"r3"<<endl;
		}
	}
	
	
}
//insert matrix B at column which_col, after the function, first col of B should be at which_col
void insert_col(GF2mat& output_e,int which_col,const GF2mat& B)
{
	if (which_col<0 or which_col>output_e.cols()) {cout<<"insetr_col error: which_col is out of range"<<endl;return;}
	
	if (which_col==0){output_e=merge_mat_hori(B,output_e);}
	else if (which_col==output_e.cols()){output_e=merge_mat_hori(output_e,B);}
	else
	{
		GF2mat temp1=output_e.get_submatrix(0,0,output_e.rows()-1,which_col-1);
		GF2mat temp2=output_e.get_submatrix(0,which_col,output_e.rows()-1,output_e.cols()-1);
		GF2mat temp3=merge_mat_hori(temp1, B);
	    output_e=merge_mat_hori(temp3,temp2);
	}
}	

//insert matrix B at row which_row, after the function, first row of B should be at which_row
void insert_row(GF2mat& output_e,int which_row,const GF2mat& B)
{
	if (which_row<0 or which_row>output_e.rows()) {cout<<"insert_row error: which_row="<<which_row<<"  output_e.rows()="<<output_e.rows()<<endl;return;}
	if (B.cols()!=output_e.cols()){cout<<"insert_row error: output_e.cols()="<<output_e.cols()<<"  B.cols()="<<B.cols()<<endl;return;}
	
	if (which_row==0){output_e=merge_mat_vert(B,output_e);}
	else if (which_row==output_e.rows()){output_e=merge_mat_vert(output_e,B);}
	else
	{
		GF2mat temp1=output_e.get_submatrix(0,0,which_row-1,output_e.cols()-1);
		GF2mat temp2=output_e.get_submatrix(which_row,0,output_e.rows()-1,output_e.cols()-1);
		GF2mat temp3=merge_mat_vert(temp1, B);
	    output_e=merge_mat_vert(temp3,temp2);
	}
}	
//note:output_e is a column vector
void small_row_reduction(GF2mat &H,GF2mat& G, GF2mat& output_e, int col_ind,  vector<vector<int>>& Ai, mat& mcv, mat& mvc)
{ 		
		if (H.cols()!=G.cols()){cout<<" small_row_reduction error:  H.cols()="<<H.cols()<<"  G.cols()="<<G.cols()<<endl;return;}
	    if (G.rows()==0) {return;}
		if (col_ind<0 or col_ind>H.cols()-1) {cout<<"small_row_reduction error: H.cols()="<<H.cols()<<"  col_ind="<<col_ind<<endl;return;}
		
		//cout<<"s1"<<endl;
		vector<int> A;
		int which_row=-1;
		//the first element of A is the index for the deleted col:
		
		A.push_back(col_ind);
		// find a row to perform row reduction
	//	cout<<"s2"<<"G.rows="<<"G.rows()"<<"  col_ind="<<col_ind<<endl;
		for (int i=0;i<G.rows();i++)
		{
			if (G(i,col_ind)==1) {which_row=i;A.push_back(col_ind);break;}
		}
		
		//find where are the 1s in this row:
		//cout<<2<<endl;
	//	cout<<"  205 G.cols()="<<G.cols()<<"  which_row="<<which_row<<endl;
		//cout<<1<<endl;
		if (which_row!=-1)
		{
			for(int i=0;i<G.cols() ;i++)
			{
			   //cout<<"G.cols()="<<G.cols()<<endl;
			
				if (G(which_row,i)==1 and i!=col_ind)
					//perform row reduction, add col_indth col to ith col
				{
					//cout<<21<<endl;
					A.push_back(i);
					//cout<<3<<endl;
					for (int j=0;j<G.rows();j++)
					{
						G.set(j,i, G(j,i)+G(j,col_ind));
					}
					output_e.set(i,0,output_e(i,0)+output_e(col_ind,0));
					//cout<<4<<"i="<<i<<endl;
				}
			}
			//cout<<2<<endl;
			//cout<<41<<endl;
			if (G.rows()>1){del_row(G,which_row);}
			else 
			{
				GF2mat tempG(1,1);
				G=tempG;
			}
			//cout<<42<<endl;
		}
		/*
		 for (int j=0;j<H.rows();j++)
		{
			H.set(j,col_ind,0);
		}
		
		GF2mat zeromat2(H.rows(),G.rows());
		if (H* G.transpose()==zeromat2){}
		else{cout<<"small_row_reduction error before deleting rows/cols: H*G^T="<<H* G.transpose()<<"\n H is\n"<<H<<"\n G is\n"<<G<<endl;}
		
		*/
		//cout<<"before delete G=\n"<<G<<"\n H=\n"<<H<<endl;
		//cout<<25<<endl;
		Ai.push_back(A);
		if (G.cols()>1) {del_col(G,col_ind);}
	   // cout<<5<<endl;
		del_col (H, col_ind);
		// cout<<6<<endl;
		mcv.del_col (col_ind);
		mvc.del_col (col_ind);
		// cout<<7<<endl;
		del_row(output_e,col_ind);
	//	 cout<<8<<endl;
			//cout<<"after delete G=\n"<<G<<"\n H=\n"<<H<<endl;
		// check if this is right:
		//cout<<3<<endl;
		if (not (G.rows()==1 and G.cols()==1))
		{  
			//cout<<"before check G=\n"<<G<<"\nH=\n"<<H<<endl;
			GF2mat zeromat(H.rows(),G.rows());
		//	cout<<9<<endl;
			if (H* G.transpose()==zeromat){}
			else{cout<<"small_row_reduction error: H*G^T="<<H* G.transpose()<<"\n H is\n"<<H<<"\n G is\n"<<G<<endl;}
		}
}
			
			
	void del_row (GF2mat & H, int row_ind)
{
	if (row_ind<0 or row_ind>H.rows()-1) {cout<<"del_col: wrong col_ind H.cols()="<<H.cols()<<"  row_ind="<<row_ind<<endl;return;}
	
	 if(H.rows()==1){cout<<"del_row: the matrix only have one row, cannot delete it (GF2mat cannot creat a matrix with no entires"<<endl;return;}
	if (row_ind==H.rows()-1){H=H.get_submatrix(0,0,H.rows()-2,H.cols()-1);}
	
	else if (row_ind==0){H=H.get_submatrix(1,0,H.rows()-1,H.cols()-1);}
	else
	{
	
		GF2mat H1=H.get_submatrix(0,0,row_ind-1,H.cols()-1);
		GF2mat H2=H.get_submatrix(row_ind+1,0,H.rows()-1,H.cols()-1);
		H=merge_mat_vert(H1,H2);
		
	}

}
/*
mat merge_mat_hori(const mat &left,const mat &right)
{
  if (left.rows()!=right.rows())
    {
      cout<<"the 2 matrices cannot merge horizontally"<<endl;
      mat error(1,1); 
      return error;
    }

  else
    {
      int r=left.rows();
      int c=left.cols()+right.cols();
      int c1=left.cols();
      int c2=right.cols();
      mat m(r,c);

      for (int i=0;i<r;i++)
	{
	  for (int j1=0;j1<c1;j1++){
	    m.set(i,j1,left(i,j1));
	  }
	  for (int j2=0;j2<c2;j2++){
	    m.set(i,c1+j2,right(i,j2));
	  }	  
	}
       return m;
    }
}
*/

void right_small_row_reduction(bool& empty_G,GF2mat &syndrome,GF2mat &H,GF2mat& G, GF2mat& output_e, int col_ind,  vector<vector<int>>& Ai ,mat& mcv, mat& mvc, vector<GF2mat> & row_del, vector<int>& B,int& hardSelectValue,int&LLR_fail,int&LLR_suc, int debug)
{

		if (H.cols()!=G.cols() or G.cols()!=output_e.rows()){cout<<" small_row_reduction error:  H.cols()="<<H.cols()<<"  G.cols()="<<G.cols()<<"  output_e.rows()="<<output_e.rows()<<endl;return;}
	    if (G.rows()==0) {return;}
		if (col_ind<0 or col_ind>H.cols()-1) {cout<<"small_row_reduction error: H.cols()="<<H.cols()<<"  col_ind="<<col_ind<<endl;return;}
		
		//cout<<"s1"<<endl;
		//vector<int> A;
		int which_row=-1;
		//the first element of A is the index for the deleted col:
		vector<int> A;
		A.push_back(col_ind);

		for (int i=0;i<G.rows();i++)
		{
			if (G(i,col_ind)==1) {which_row=i;A.push_back(col_ind);break;}
		}
		
        A.push_back(which_row);
		
		//if G contains an one in this column:
		if (which_row!=-1)
		{
			
			// add this row to output_e, although this may not affect BP, so I guess I don't need this
			if (output_e.get(col_ind,0)==1)
			{
				//cout<<"add a row of G to output_e"<<endl;
				A.push_back(-1);//-1 for add this row to e
				for (int j=0;j<G.cols();j++)
					{
						output_e.set(j,0,output_e(j,0)+G(which_row,j));
					}
			}
			else{A.push_back(-2);}
			row_del.push_back(G.get_submatrix(which_row,0,which_row,G.cols()-1));  //store  the row added to e
			for(int i=0;i<G.rows() ;i++)
			{
			   //cout<<"G.cols()="<<G.cols()<<endl;
			
				if (G(i,col_ind)==1 and i!=which_row)
					//perform row reduction, add which_rowth row to ith row
				{
					
					//cout<<21<<endl;
					A.push_back(i);
					//cout<<3<<endl;
					for (int j=0;j<G.cols();j++)
					{
						G.set(i,j, G(i,j)+G(which_row,j));
						//output_e.set(j,0,output_e(j,0)+G(which_row,j));
					}
					
					//cout<<4<<"i="<<i<<endl;
				}
			}
			
			if (G.rows()>1){del_row(G,which_row);}
			else 
			{
				empty_G=true;
			}
		
		}
	
	else if(debug%2==1)
	{
		A.push_back(-2);
		
		GF2mat Perm;
		row_del.push_back(Perm);
		GF2mat U=row_gaussian( H, Perm) ;
		//cout<<" after gaussian elimination: H is\n"<<U<<"\nG is (unchanged)\n"<<G<<endl;
		GF2mat new_syndrome=Perm*syndrome;
		int wt_1_row=-1;	
		
		int hard_select=-1;
		bin bin_one=1;
		
		for (int iii=U.rows()-1;iii>=0;iii--)
		{
	
			if(U.get(iii,col_ind)==bin_one)  // if this row has an 1 hit the right column, calculate wt of this row:
			{
				int row_wt=0;
				for (int ii=0;ii<U.cols();ii++)
				{
					//cout<<1<<endl;
					if (U.get(iii,ii)==1){row_wt++;}
					//cout<<2<<endl;
				}
				if (row_wt==1){wt_1_row=iii;break;} //if row_wt==1, this is ok
			}
		}
		if (wt_1_row==-1) 
		{
			//cout<<"error: "<<col_ind<<"th col of H does not hit a wt_1 row: after row gaussian elimination H is \n"<<U<<"\n G is\n"<<G<<endl;
			//an empty column of H gives no information, just del them later.
		}
		else
		{   
			hard_select=new_syndrome.get(wt_1_row,0)==0? 0:1;
			if (output_e(col_ind,0)!=hard_select){output_e.set(col_ind,0,hard_select);LLR_fail++;}
			else{LLR_suc++;}
			// do not need to delete a row from H because only U has a wt=1 row rather than H but need to modify s
			//if (H.rows()>1){del_row(H,wt_1_row);mcv.del_row (wt_1_row);mvc.del_row (wt_1_row);}
			if (new_syndrome.get(wt_1_row,0)==1) {new_syndrome.set(wt_1_row,0,0); syndrome=Perm.inverse()* new_syndrome;}
		}
	}
	
	else{
		A.push_back(-2);
		
		GF2mat Perm;
		row_del.push_back(Perm);
	}
		Ai.push_back(A);
		if (G.cols()>1) {del_col(G,col_ind);}
		else 
			{
				empty_G=true;
			}
	
		del_col (H, col_ind);

		mcv.del_col (col_ind);
		mvc.del_col (col_ind);
		B.push_back(output_e.get(col_ind,0));
		del_row(output_e,col_ind);
	
		if (not 	empty_G)
		{  

			GF2mat zeromat(H.rows(),G.rows());

			if (H* G.transpose()==zeromat){}
			else{cout<<"R_small_row_reduction error: H*G^T="<<H* G.transpose()<<"\n H is\n"<<H<<"\n G is\n"<<G<<endl;}
		}
}
			
void right_restore_e(GF2mat& output_e,  const GF2mat& original_output_e, const vector<vector<int>> &pos_deleted_entries, const vector<int>& value_deleted_entries,vector<GF2mat>& rows_del_list)
{
	
	if (pos_deleted_entries.size()!=value_deleted_entries.size()){cout<<"restore_e error: pos_deleted_entries.size="<<pos_deleted_entries.size()<<"  value_deleted_entries.size="<<value_deleted_entries.size()<<endl;return;}
	
	if (output_e.cols()!=1){cout<<"restore_e error: output_e.cols()!=1"<<endl;return;}	
	output_e=original_output_e;
//	cout<<"r1"<<endl;
   //int counter1=rows_del_list.size()-1;
 //  cout<<"add those values to e:"<<endl;
//  for (int i=0;i<value_deleted_entries.size();i++){cout<<value_deleted_entries[i]<<"  ";}
  //cout<<endl;

	for (int ii=value_deleted_entries.size()-1;ii>=0;ii--)
	{
		
			GF2mat value_to_be_res(1,1);
			value_to_be_res.set(0,0,value_deleted_entries[ii]);
		
			//the first element of Ai[j] is the postion deleted in jth round
			int row_to_be_res=pos_deleted_entries[ii][0];	
			insert_row(output_e,row_to_be_res,value_to_be_res);
			
			if ( pos_deleted_entries[ii][2]==-1)
			{
				//cout<<"s"<<endl;
				output_e=output_e+rows_del_list[ii].transpose();
				//cout<<"e"<<endl;
				//counter1--;
			}
			//else if(counter1<-1) {cout<<"restore error: counter1="<<counter1<<"   rows_del_list.size()="<<rows_del_list.size()<<endl;return;}	
	}
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


/*
ivec permute_GF2mat_columns(itpp::GF2mat H) {
  int r = H.rows();
  int c = H.cols();
  
  if (c>r){}
  else{cout<<"permute_GF2mat_columns error:H.cols()<H.rows()"<<endl;}
  // Initialize the permutation vector as an identity permutation
  ivec permutation(c);
  for (int i = 0; i < c; ++i) {
    permutation[i] = i;
  }

    GF2mat T;
  GF2mat U;
  ivec perm;
  
    int current_col=1;
	int swaped_col=r;
	while (current_col<r)
	{
		GF2mat H_S=H.get_submatrix(0,0,r-1,current_col);
		int rankHS=H_S.T_fact(T,U,perm);
		if (rankHS==current_col+1){current_col++;}
		else 
		{
		H.swap_cols(current_col, swaped_col);
		 std::swap(permutation[current_col], permutation[swaped_col]);
		swaped_col++;
		}
	}


  return permutation;
}
*/

bool  Fast_OSD(const GF2mat &Hx, const GF2mat &Gx,const nodes checks[],const nodes errors[],const vec &pv, vec&pv_dec,
const double& pmin,const double& pmax,double& num_iter, int lmax,int &wt,  int debug, vec &LR2,int rankH,int& max_rep)
{
	// variables with prefix const_ are the variables before any row reduction
  int v=Hx.cols();
  int c=Hx.rows();

  int wt_real_e=0;
  GF2mat real_e(v,1); 

  if (wt==0)  
    {
      wt_real_e=error_channel(real_e, pv);
      if (wt_real_e==0)
		{
			num_no_error++;
			return true;
		}
    }
  else
    {
      error_channel2(real_e,wt);
    }


  GF2mat zero_rvec(1,v);
   
  //LR.zeros();
  //vec LR_avg=LR;
  GF2mat zero_mat1(c,1);
  GF2mat zero_mat2(v,1);
  GF2mat zero_rvec2(Gx.rows(),1);
  
  GF2mat syndrome=Hx*real_e;
  
  //is the syndrome a zero vector?
  if (syndrome==zero_mat1)
	{
	  // is e a stablizer?
	  if (Gx*real_e==zero_rvec2)
	    {
	      return true;
	    }        
	  else  
	    {
	      syn_fail++;
	      return false;	     
	    }  
	}
	
      mat mcv(c,v);   // the messages from check nodes to variable nodes, which are p0/p1
      mat mvc(c,v);   // the messages from variable nodes to check nodes
            
      vec LR(v);

    
	int l2=0;
	GF2mat output_e(v,1);
	GF2mat output_e1(v,1);
	GF2mat new_syndrome=syndrome;
	
	while (l2<max_rep)
	{
				
		mcv.zeros();
		mvc.zeros();
		initialize_massages( mcv,mvc, Hx); //initialize to all-1 matrix		
		LR.zeros();
				
		for (int l=1;l<=lmax;l++)
		{	   				
			quan_s_update(checks,errors, mcv,mvc,new_syndrome,pv_dec, c, v,output_e1,LR,debug);	
					
			if (Hx*(output_e+output_e1)==syndrome)		   
			{		  				
				  if(Gx*(output_e1+output_e+real_e)==zero_rvec2)
					{
						num_iter= num_iter+l;							
						return true;
					}		    
					else
					{					
						syn_fail++;		
						return false;      	        
					}	    
                          
			}					
		}
		
		new_syndrome=new_syndrome+Hx*output_e1;
		output_e=output_e+output_e1;	
		/*
		int wt_new=s_weight(new_syndrome);
		int d=sqrt(v);
		if (wt_new==1)
		{
			for (int i3=0;i3<v;i3++)
			{ 
				
					GF2mat output_e2(v,1);
					output_e2.set(i3,0,1);
					GF2mat new_syn2=Hx*output_e2+new_syndrome;
					int wt_new2=s_weight(new_syn2);
					if (wt_new2==0)
					{										
						if(Gx*(output_e2+output_e+real_e)==zero_rvec2)
						{
							//num_iter= num_iter+l;							
							return true;
						}		    
						else
						{					
							syn_fail++;		
							return false;      	        
						}	    
					}
					else if (wt_new2<wt_new)
					{
						
						cout<<"output_e2  is"<<endl;
						err_pos2(output_e2);
						cout<<"real_e is"<<endl;
						err_pos2(real_e);
						cout<<"new_syn2 is"<<new_syn2.transpose()<<endl;
						cout<<"new_syn is"<<new_syndrome.transpose()<<endl;
					
						new_syndrome=new_syndrome+new_syn2;
						output_e=output_e+output_e2;	
					}
				
			}
			
		}
		*/
		l2++;
		pro_dist(pmin,pmax,pv_dec);	
	}
	
	
      	           
      if (debug&8)
	{
		int suc_order=0;
	 	  
	  OSD(LR,Hx,syndrome,output_e,suc_order);

	 
	  if(Gx*(output_e+real_e)==zero_rvec2)
	    {
	      // cout<<"OSD suc2.1"<<endl;
	      OSD_suc++;		 
		  if (suc_order==0){OSD_0_suc++;}
		else if(suc_order==1){OSD_1_suc++;}
		else if(suc_order==2){OSD_2_suc++;}
		else {OSD_higher_order_suc++;}		
	      return true;
	    }	 
		else 
		{
			if (debug&4)
	 {
		 cout<<"\n\n\nxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
		 cout<<"OSD fails:"<<endl;
	   cout<<"output e is \n"<<output_e.transpose()<<endl;
	   err_pos2(output_e);
	   
	  // cout<<"OSD-0 e_S="<<global_e_S.transpose()<<endl;
	  // cout<<"new_wt="<<global_new_wt<<endl;
	   
	   cout<<"real e is \n"<<real_e.transpose()<<endl;
	   err_pos2(real_e);
	   GF2mat sume=real_e+output_e;
	   cout<<"residual e is \n"<<endl;
	   err_pos2(sume);
	   cout<<"syndrome is"<<syndrome.transpose()<<endl;
	  // cout<<"H*residual_e="<<(Hx*sume).transpose()<<endl;
	   //  cout<<"mcv is \n"<<mcv<<endl;
	   //  cout<<"\n mvc is \n"<<mvc<<endl;
	   cout<<"LR is \n"<<endl;
	   int dist=sqrt(Hx.cols());
	   for (int ii=0;ii<LR.length();ii++)
	   {
		   cout<<LR(ii)<<"   ";
		   if (ii%dist==dist-1){cout<<endl;}
	   }
	   
	 }

			OSD_fail++;		
			syn_fail++;		
			return false;
		}
	}
 
  if (debug&4)
	 {
		 cout<<"\n\n aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"<<endl;
		 GF2mat out_syn=Hx*output_e;
		 GF2mat sume=real_e+output_e;
	   GF2mat res_syn=Hx*sume;	 
		 
		  cout<<"wt of real_e is "<<s_weight(real_e)<<endl;
		  cout<<"wt of syndrome is "<<s_weight(syndrome)<<endl;
		  cout<<"wt  of output_e is "<<s_weight(output_e)<<endl;
		   cout<<"wt of output_e syndrome is "<<s_weight(out_syn)<<endl;
		    cout<<"wt of real_e+output_e is "<<s_weight(real_e+output_e)<<endl;
           cout<<"wt of output_e + real_e syndrome is "<<s_weight(res_syn)<<endl;
		 
		  cout<<"real e is \n"<<real_e.transpose()<<endl;	 
	      err_pos2(real_e);   
	    cout<<"syndrome is"<<syndrome.transpose()<<endl;
			   
	   cout<<"output e is \n"<<output_e.transpose()<<endl;	   
	   err_pos2(output_e);
	   
	   cout<<"output_e syndrome is"<<out_syn.transpose()<<endl;
	   
	   	   
	   cout<<"output_e + real_e is \n"<<sume.transpose()<<endl;
	    err_pos2(sume);
	   cout<<"syndrome of output_e + real_e is"<<res_syn.transpose()<<endl;
		   
	   cout<<"LR is \n"<<endl;
	   int dist=sqrt(Hx.cols());
	   for (int ii=0;ii<LR.length();ii++)
	   {
		   cout<<LR(ii)<<"     ";
		   if (ii%dist==dist-1){cout<<endl;}
	   }
	   
	 }
	 int wt_red=s_weight(output_e+real_e);
	 if (wt_red==1){wt_1_residual_error++;}
	 else if(wt_red==2){wt_2_residual_error++;}
else { large_wt_residual_error++;}
      max_fail++;
      return false;
 }
	ivec reliability_sort_index (const GF2mat & H, const vec& LR)
	{
		 int n=LR.length();
  if (n!=H.cols()){cout<<"reliability error: LR.length!=H.cols"<<endl;}
  int r=H.rows();
  vec abs_LLR(n);
  vec reliability(r);
  for (int i=0;i<n;i++)
    {
      if (LR(i)>=0) 
	{
	  abs_LLR(i)=abs(log(LR(i)));
	}
      else
	{
	  cout<<"reliability error: negative LR!"<<endl;
	  cout<<LR(i)<<endl;
	}
    }
  
  for (int i=0;i<r;i++)
  {
	  reliability(i)=0;
	  for (int j=0;j<n;j++)
	  {
		  if (H.get(i,j)==1){reliability(i)=reliability(i)+abs_LLR(j);}
	  }
  }
	//the sort function gives a ascending  order
	ivec perm=sort_index(reliability);
	return perm;
	}
	
	void reliability_perm_H (const GF2mat & H, int which_row,GF2mat &output_H,GF2mat& output_col_perm,GF2mat& output_row_perm, GF2mat& H2,GF2mat& A,GF2mat& output_H1)
	{
		GF2mat H1=H;
		int n=H.cols();
		int r=H.rows();
		int swapped_col=0;
		bmat I_n_bmat=eye_b(n);
		GF2mat col_perm(I_n_bmat);
		bmat I_r_bmat=eye_b(r);
		GF2mat row_perm(I_r_bmat);
		
		for (int i=0;i<n;i++)
		{
			if (H1.get(which_row,i)==1)
			{
				H1.swap_cols(i,swapped_col);
				col_perm.swap_cols(i,swapped_col);
				swapped_col++;				
			}
		}
		
		int swapped_row=0;
		
		for (int i=0;i<r;i++)
		{
			bool has_one=false;
			for (int j=0;j<swapped_col;j++)
			{
				if (H1.get(i,j)==1){has_one=true;break;}
			}
			
			if (has_one==true)
			{
				H1.swap_rows(i,swapped_row);
				row_perm.swap_rows(i,swapped_row);
				swapped_row++;	
			}
		}
		output_col_perm=col_perm;
		output_row_perm=row_perm;
		output_H=H1;
		H2=H1.get_submatrix(swapped_row,swapped_col,H1.rows()-1,H1.cols()-1);
		
		A=H1.get_submatrix(0,swapped_col,swapped_row-1,H1.cols()-1);
		output_H1=H1.get_submatrix(0,0,swapped_row-1,swapped_col-1);
	/*	
		if (output_H1.cols()+A.cols()!=H1.cols()){cout<<"output_H1.cols()+A.cols()!=H1.cols()"<<endl;}
		if (H2.rows()+A.rows()!=H1.rows())
		{
			cout<<"H2.rows()+A.rows()!=H1.rows()"<<endl;
			cout<<"H is"<<H1<<endl;
			cout<<"H2 is "<<H2<<endl;
			cout<<"A is "<<A<<endl;
	}
	*/
	}
	
	bool reliability_solve(const GF2mat&H1,const GF2mat&A, const GF2mat &s1, const GF2mat &e2,GF2mat & e1)
	{

		
		GF2mat news=s1+A*e2;
		
		int rankH1=GF2mat_rank(H1);
		GF2mat H1news=merge_mat_hori(H1,news);
		int rankH1news=GF2mat_rank(H1news);
		vec LR(H1.cols());
		GF2mat e(H1.cols(),1);
		for (int i=0;i<H1.cols();i++){LR(i)=1;}
		
		int temp;
		if (rankH1<rankH1news){return false;}
		else 
		{
			
			OSD(LR, H1, news,e,temp);
			e1=e;
			return true;
		}
			
		
	}
	
	GF2mat reliability_restored_e(const GF2mat col_perm, const GF2mat &e1,const GF2mat &e2)
	{
		GF2mat e=merge_mat_vert(e1,e2);
		e=col_perm*e;
		return e;
	}
	


bool   circuit_decoder(const Sparse_GF2 &sparseH,const Sparse_GF2 &sparseL,const GF2mat &Hx, const GF2mat Gx, GF2mat G,const GF2mat &real_e,GF2mat &final_output_e,const nodes checks1[],const nodes errors1[],const vec&pv_dec,double& num_iter, int lmax,int &wt,  int debug, vec &LR2,int rankH,double M,int&LLR_fail,int&LLR_suc,int& num_fixed_LLR,bool& start_row_reduction,double alpha,int other_decoder,int channel,int max_row_freeze,int OSD_order)
{
	// variables with prefix const_ are the variables before any row reduction

 int c0=Hx.rows();
 int v0=Hx.cols();
  GF2mat H=Hx;
  		int v=H.cols();
		int c=H.rows();
		GF2mat Hz=G;

if (other_decoder==1 and M>=0){start_row_reduction=true;}
else{start_row_reduction=false;}
  
  
  
  //if no error, break

  GF2mat zero_rvec(1,v);
   
  //LR.zeros();
  //vec LR_avg=LR;
  GF2mat zero_mat1(c,1);
  GF2mat zero_mat2(v,1);
  //GF2mat tempmat333=Gx*real_e;
  GF2mat zero_rvec2(Gx.rows(),1);
   GF2mat const_zero_rvec2(Gx.rows(),1);
  
  GF2mat syndrome=Hx*real_e;
  GF2mat original_syndrome=syndrome;
  

  //is the syndrome a zero vector?
  if (syndrome==zero_mat1)
	{
	  // is e a stablizer?
	  if (Gx*real_e==zero_rvec2)
	    {
			num_no_error++;
			final_output_e=real_e;
	      return true;
	    }        
	  else  
	    {
	      syn_fail++;
	      return false;	     
	    }  
	}
	Sparse_GF2 Sparse_real_e(real_e);
	Sparse_GF2 Sparse_syndrome(syndrome);
	Sparse_GF2 Sparse_zero_rvec2(zero_rvec2);
	
      mat mcv(c,v);   // the messages from check nodes to variable nodes, which are p0/p1
      mat mvc(c,v);   // the messages from variable nodes to check nodes
      mcv.zeros();
      mvc.zeros();
    
      initialize_massages( mcv,mvc, H); //initialize to all-1 matrix

      
	  GF2mat res_output_e(v,1);
	  GF2mat const_output_e(v,1);
    
      vector<GF2mat>rows_del;// each element stores a row of G deleted 
	  vector<int>value_deleted_entries;// values deleted 
	  
	  vector<vector<int>> pos_deleted_entries;  // each element stores the columns and rows deleted/changed in a round of row reduction: 
	  //example of an element: [[4,8,-1,2,3],[5,10,-2,7,9]] , means in this round, 4th and 5th columns of H are deleted, 
	  //8th row of G is added to 2h, 3th rows, 10th column is added to 7th,9th rows.
	  //-1 means this row is also added to e, -2 means not.
      int E1=0;
	  
	vec old_LR(v);
	
	vec LR_avg(v);
	vec LR_prod(v);
	for (int i=0;i<v;i++){LR_prod(i)=1;}
	
	
	vec LR_final;
	bool empty_G=false;
	bool converge_again=false; //is BP converged again after row reduction?
	

	int l2=0;
	
	GF2mat zero_e(v,1);
	final_output_e=zero_e;
	
		//pro_dist(pmin,pmax,pv_dec);
	
		GF2mat output_e(v,1);
		
		
      for (int l=1;l<=lmax;l++)
	{
		v=H.cols();
		c=H.rows();
	  //cout<<l<<"th iteration starts"<<endl;
	   nodes* checks = new nodes[H.rows()];
	  nodes* errors = new nodes[H.cols()];  
	  
	   //std::unique_ptr<nodes[]> checks(new nodes[H.rows()]);
   // std::unique_ptr<nodes[]> errors(new nodes[H.cols()]);
	
	  initialize_checks (H, checks,  E1);
	  initialize_errors(H, errors);
	 vec LR(v);
	 
			if (debug&2)
			{
				quan_p_update(checks,errors, mcv,mvc,syndrome,pv_dec, c, v,output_e,LR,debug,alpha);				
			}
	      else
			{
				//cout<<"serial V starts"<<endl;			
				quan_s_update(checks,errors, mcv,mvc,syndrome,pv_dec, c, v,output_e,LR,debug,alpha);
			}
	    
	    
		//cout<<"l="<<l<<" before check after updating H=\n"<<H<<"\noutput_e=\n"<<output_e<<"\nsyndrome=\n"<<syndrome<<endl;
		if (debug&4 and start_row_reduction==true)
		{
			cout<<"l="<<l<<"\n LLR=\n"<<endl;
			for (int aaa=0;aaa<LR.size();aaa++){cout<<log(LR(aaa))<<"  ";}
			cout<<endl;
			cout<<"output_e is\n"<<output_e.transpose()<<endl;
			right_restore_e(res_output_e,output_e,  pos_deleted_entries,value_deleted_entries,rows_del);
			cout<<"after restore"<<endl;
			err_pos2(res_output_e);
			cout<<"H is\n"<<H<<"\n G is\n"<<G<<"\nsyndrome is\n"<<syndrome.transpose()<<endl;
		}
	
			if (debug&32)
			{
					for (int iii=0;iii<v;iii++)
				{
					LR_prod(iii)=pow(LR_prod(iii),0.9)*LR(iii);
					LR_avg(iii)=pow(LR_prod(iii),1.0/(1-pow(0.9,l)));
					output_e.set(iii,0,LR_avg(iii)<1? 1:0);   
				}
			}
			Sparse_GF2 Sparse_output_e(output_e);
			
	   if (sparseH*Sparse_output_e==Sparse_syndrome)
		   
		{		  final_output_e=final_output_e+output_e;
				//cout<<"check1 is ok"<<endl;
				if  (H.cols()!=v0)
				{
				 right_restore_e(res_output_e,output_e,  pos_deleted_entries,value_deleted_entries,rows_del);
				}
				else{res_output_e=final_output_e;}
				Sparse_GF2 Sparse_res_output_e(res_output_e);
				
				  if(sparseL*(Sparse_res_output_e+Sparse_real_e)==Sparse_zero_rvec2)
					{
						num_iter= num_iter+l;		
						if(start_row_reduction==true){fix_suc++;}			
							if (checks !=NULL){delete[] checks;}
							if(errors!=NULL) {delete[] errors;}
							BP_suc++;
						
							if (debug&32) {avg_BP_suc++;}
						final_output_e=res_output_e;
						return true;
					}
		    
		  else
		    {
		      // cout<<"\n syndrome is ok, but decoding fails:"<<endl;
			  if(start_row_reduction==true){fix_fail++;}			
		      syn_fail++;		
			  	if (checks !=NULL){delete[] checks;}
		if(errors!=NULL) {delete[] errors;}
		      return false;      	        
		    }	    	  
		}
		

	
	//cout<<"start row reduction0"<<endl;
		//where to perform reduction?
		//cout<<"G is \n"<<G<<"\n H is \n"<<H<<endl;
		if ( start_row_reduction==false and (LR==old_LR or l>0.5*lmax) and M>=0)
		
		{
			//cout<<"M="<<M<<endl;
			start_row_reduction=true;
			if (l<0.5*lmax) {converge_to_trapping_set++;}
			else {did_not_converge_for_large_iterations++;}
			
			if (debug&4 )
			{
			cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
			if (l<0.5*lmax) {cout<<" \n l="<<l<<"  LLR converged, start row reduction: "<<endl;}
			else {cout<<" \n l="<<l<<"  but LLR did not converge, perhaps meets an oscillation or maybe it's a matter of calculation accuracy, start row reduction: "<<endl;}
			cout<<"real_e is"<<real_e.transpose()<<endl;
			err_pos2(real_e);
			cout<<"output_e is\n"<<output_e.transpose()<<endl;		
			err_pos2(output_e);
			
			cout<<"\n LLR=\n"<<endl;
			for (int aaa=0;aaa<LR.size();aaa++){cout<<log(LR(aaa))<<"  ";}
			cout<<endl;
			cout<<"H is\n"<<H<<"\n G is\n"<<G<<"\nsyndrome is\n"<<syndrome.transpose()<<endl;
			}
			
		}
		
		if (  (not empty_G ) and  start_row_reduction and M>=0 )
		{
			//cout<<"del some cols of H"<<endl;
			int num_del_col=0;
			int hardSelectValue=-1;
			bool exists_larger_than_M=true;
			vec abs_LLR(LR.size());
			for (int jjj=0;jjj<LR.size();jjj++)
			{
				abs_LLR(jjj)=abs(log(LR(jjj)));
			}
		    
			while( exists_larger_than_M)
			{
				ivec sortindex=sort_index(abs_LLR);// ascending order
				int max_ind=sortindex(sortindex.size()-1);	
						
				double max_abs_LLR=abs_LLR(max_ind);
				
					if (max_abs_LLR>M)
						{			 
							num_fixed_LLR++;
							if (not empty_G)
							{
								right_small_row_reduction(empty_G,syndrome,H,G,output_e,max_ind, pos_deleted_entries,mcv,mvc,rows_del,value_deleted_entries, hardSelectValue,LLR_fail,LLR_suc,debug);	
								abs_LLR.del(max_ind);
							}
							else{break;}
					
							num_del_col++;
						}
					else {exists_larger_than_M=false;}
			}
					if (debug&4){cout<<"delete "<<num_del_col<<" columns in this round"<<endl;}			
					
					if (not empty_G)
					{
						if (GF2mat_rank(H)+GF2mat_rank(G)!=H.cols()){cout<<"error: after row reduction GF2mat_rank(H)+GF2mat_rank(G)!=H.cols(): H is\n"<<H<<"\n G is \n"<<G<<endl;}
					}
		 }
			
		  //protocol B:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111
			//if ( debug&16 and l==lmax and (start_row_reduction or M==-50))
				//if (l==lmax) {cout<<"should goto pb"<<endl;}
			//ltemp2=l;
			if ( (debug&16) and l==lmax and (start_row_reduction or M==-10))
			{
				    //cout<<"P_B"<<endl;
					algebraic_decoder(H, syndrome,output_e);
					
					right_restore_e(res_output_e,output_e,  pos_deleted_entries,value_deleted_entries,rows_del);
					
					if (debug&4 )
					{
						cout<<"meet max iterations, start protocol B,  which output the error (after restoring):\n"<<endl;
						err_pos2(res_output_e);
					}
					
					 if (Hx*res_output_e==original_syndrome)
					 {
						if (Gx*(res_output_e+real_e)==const_zero_rvec2)
						{				 
							//num_iter= num_iter+l;		
							fix_suc++;
							protocol_B_suc++;
								if (checks !=NULL){delete[] checks;}
		if(errors!=NULL) {delete[] errors;}
		final_output_e=res_output_e;
							return true;					
						}
						else
						{
						syn_fail++;
						protocol_B_syn_fail++;
						fix_fail++;
							if (checks !=NULL){delete[] checks;}
		if(errors!=NULL) {delete[] errors;}
						return false;				 
						}
					 }
					 else
					 {
						// cout<<"PB max fail"<<endl;
						protocol_B_max_fail++;
						fix_fail++;
						max_fail++;
							if (checks !=NULL){delete[] checks;}
		if(errors!=NULL) {delete[] errors;}
						return false;					 
					 }
						
			}		 
		
		
		//cout<<"l="<<l<<"  over , go to next iteration "<<"G is\n"<<G<<"\nH is\n"<<H<<endl;
		if (old_LR.size()==LR.size()){old_LR=LR;}

		if (checks !=NULL){delete[] checks;}
		if(errors!=NULL) {delete[] errors;}
		//ltemp=l;
		LR_final=LR;
	}
	
	if (debug&4 )
		{
			cout<<"BP failed:"<<endl;
	
			cout<<"output_e is\n"<<endl;
			v=H.cols();
		for (int index5=0;index5<v;index5++)
		{
			if (output_e(index5,0)==1){cout<<index5<<" ";}
					
		}
		cout<<endl;
		cout<<"real_e is\n"<<endl;
		for (int index5=0;index5<v;index5++)
		{
			if (real_e(index5,0)==1){cout<<index5<<" ";}
					
		}
		cout<<endl;
		}
	final_output_e=output_e+final_output_e;
	syndrome=syndrome+H*final_output_e;
	
	syndrome=original_syndrome;
      
		if (max_row_freeze>0)
		{
			if (max_row_freeze>Hz.rows()){max_row_freeze=Hz.rows();}
			bool rel_suc=false;
			ivec R_index=reliability_sort_index ( Hz, LR_final);      
			for (int freeze_ind=0;freeze_ind<max_row_freeze;freeze_ind++)
			{
				GF2mat output_H,output_col_perm,output_row_perm,H2,A,output_e1, output_H1,output_e;
				//cout<<1<<endl;
				reliability_perm_H (Hx, R_index(freeze_ind),output_H,output_col_perm, output_row_perm, H2, A,output_H1);
				GF2mat output_e2(H2.cols(),1);
				//cout<<2<<endl;
				GF2mat rel_syndrome =output_row_perm*syndrome;
				GF2mat syndrome_2=rel_syndrome.get_submatrix(A.rows(),0,syndrome.rows()-1,0);
				GF2mat syndrome_1=rel_syndrome.get_submatrix(0,0,A.rows()-1,0);
				//cout<<3<<endl;
				
				     int  v2=H2.cols();
					int	c2=H2.rows();
	  
		                nodes* checks2 = new nodes[H2.rows()];
	                     nodes* errors2 = new nodes[H2.cols()];  
	  
						initialize_checks (H2, checks2,  E1);
						initialize_errors(H2, errors2);
						vec LR2(v2);
				 mat mcv2(c2,v2);   
					mat mvc2(c2,v2);   
						mcv2.zeros();
						mvc2.zeros();
                
				vec pv_dec2(v2);
				for (int temp_ind=0;temp_ind<v2;temp_ind++){pv_dec2(temp_ind)=pv_dec(temp_ind);}
				
			    initialize_massages( mcv2,mvc2, H2); //initialize to all-1 matrix
				//cout<<4<<endl;
				for (int l=1;l<=lmax;l++)
					{
							 
						if (debug&2)
							{
									quan_p_update(checks2,errors2, mcv2,mvc2,rel_syndrome,pv_dec2, c2, v2,output_e2,LR2,debug,alpha);				
							}
						else
							{	
								quan_s_update(checks2,errors2, mcv2,mvc2,rel_syndrome,pv_dec2, c2, v2,output_e2,LR2,debug,alpha);
							}
	    
			
							if (H2*output_e2==rel_syndrome){break;}
					}
					//cout<<5<<endl;
					rel_suc=reliability_solve(output_H1,A, syndrome_1, output_e2,output_e1);
					//cout<<output_H1<<endl;
					//cout<<6<<endl;
					if (rel_suc==true)
					{
						//cout<<7<<endl;
						//cout<<H<<endl;
						//cout<<output_e1<<endl;
						//cout<<output_e2<<endl;
						//cout<<output_col_perm<<endl;
						output_e=reliability_restored_e( output_col_perm, output_e1,output_e2);
						//cout<<8<<endl;
						  if(Gx*(output_e+real_e)==const_zero_rvec2)
						  {
							  freeze_rows_suc++;
							  return true;
						  }
					}
					
					if (checks2 !=NULL){delete[] checks2;}
		if(errors2!=NULL) {delete[] errors2;}
			}
		}
		   
								
		
		
	
      if (debug&8)
	{
		int suc_order=0;
		avg_red_qubits=avg_red_qubits+H.cols();
		
		if (debug&32){LR_final=LR_avg;}
	  //cout<<"OSD:M="<<M<<endl;
	  //full_OSD(LR_final,H,syndrome,output_e,2, suc_order);
	  GF2mat output_e(v,1);
	 // Sparse_OSD(sparseH,LR_final,H,syndrome,output_e,suc_order,OSD_order);

	  if (start_row_reduction==true)
	  {
		right_restore_e(res_output_e,output_e,  pos_deleted_entries,value_deleted_entries,rows_del);
	  }
	  else {res_output_e=output_e;}
	  /*
	  	if (debug&4 )
		{
			cout<<"OSD output the error (after restoring):\n"<<endl;
			err_pos2(res_output_e);
		}
		*/
	  if(Gx*(res_output_e+real_e)==const_zero_rvec2)
	    {
	      // cout<<"OSD suc2.1"<<endl;
	      OSD_suc++;
		  if(start_row_reduction==true){fix_suc++;}	
		  if (suc_order==0){OSD_0_suc++;}
		else if(suc_order==1){OSD_1_suc++;}
		else if(suc_order==2){OSD_2_suc++;}
		else {OSD_higher_order_suc++;}
		final_output_e=res_output_e;
	      return true;
	    }	 
		else 
		{
			 
	
	  
			
			if (debug&4)
	 {
		 cout<<"\n\n\nxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
		 cout<<"OSD fails:"<<endl;
	   cout<<"output e is \n"<<endl;
	   
	   
		for (int index5=0;index5<v;index5++)
		{
			if (output_e(index5,0)==1){cout<<index5<<" ";}
					
		}
		cout<<endl;
	   cout<<res_output_e<<endl;
	   err_pos2(res_output_e);
	   
	   cout<<"OSD-0 e_S="<<global_e_S.transpose()<<endl;
	  // cout<<"new_wt="<<global_new_wt<<endl;
	   
	   cout<<"real e is \n"<<endl;
	   err_pos2(real_e);
	   GF2mat sume=real_e+res_output_e;
	   cout<<"residual e is \n"<<endl;
	   err_pos2(sume);
	  // cout<<"H*residual_e="<<(Hx*sume).transpose()<<endl;
	   //  cout<<"mcv is \n"<<mcv<<endl;
	   //  cout<<"\n mvc is \n"<<mvc<<endl;
	   cout<<"LR is \n"<<endl;
	   int dist=sqrt(H.cols());
	   for (int ii=0;ii<LR_final.length();ii++)
	   {
		   cout<<LR_final(ii)<<"   ";
		   if (ii%dist==dist-1){cout<<endl;}
	   }
	   
	   //cout<<"sort index is (ascending)  "<<global_sort_index<<endl;
	   
	  // cout<<"H is\n"<<H<<endl;
	  // GF2mat H1=H;
     //  H1.permute_cols(global_sort_index,0);
	 //  cout<<"after perm1\n"<<H1<<endl;
	   cout<<"perm2 is\n"<<global_perm2<<endl;
	//    GF2mat T;

  //GF2mat U;


 // int xxx=H1.T_fact(T,U,global_perm2);
	  // cout<<"after perm2\n"<<U<<endl;
	 }
			OSD_fail++;
			if(start_row_reduction==true){fix_fail++;}	
			syn_fail++;
			
			return false;
		}
	}
 
	 if(start_row_reduction){fix_fail++;}		
//cout<<"end max fail xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;	 
//cout<<"ltemp="<<ltemp<<endl;

  if (debug&4)
	 {
	   cout<<"output e is \n"<<endl;
	   err_pos2(res_output_e);
	   cout<<"real e is \n"<<endl;
	   err_pos2(real_e);
	   GF2mat sume=real_e+res_output_e;
	   cout<<"residual e is \n"<<endl;
	   err_pos2(sume);
	   //  cout<<"mcv is \n"<<mcv<<endl;
	   //  cout<<"\n mvc is \n"<<mvc<<endl;
	   cout<<"LR is \n"<<endl;
	   int dist=sqrt(H.cols());
	   for (int ii=0;ii<LR_final.length();ii++)
	   {
		   cout<<LR_final(ii)<<"     ";
		   if (ii%dist==dist-1){cout<<endl;}
	   }
	   
	 }
      max_fail++;
      return false;
 }