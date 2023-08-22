#include <iostream>
#include<cstring>
#include<string>
#include<vector>
#include <fstream>
#include <math.h> 
#include <random>
#include <sstream>


/** @brief Belief-propagation functions

   
    @author Xingrui Liu <xliu@ucr.edu>
    @date August 2022
    */

using namespace std;
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
using namespace itpp;

#ifndef BP
#define BP
class  nodes{

public:
  ivec neighbors; /**< store the neighbors of this node */
  int degree; /**< degree of this node */
};

/*! \relates nodes
 * Initilize the check nodes, for each check, get its neighbors and degree. 

\param E stores the number of edges in this graph.
\param H is the parity check matrix.
 */
void initialize_checks (const GF2mat &H, nodes checks[], int & E);


/*! \relates nodes
 *Initilize the variable nodes.
 */
void initialize_errors(const GF2mat &H, nodes errors[]); 

/*! \relates nodes
 *Initilize the messages , set   mcv(i,j) and   mvc(i,j) to 1 if H(i,j)=1. 

\param mcv stores the messages from check nodes to variable nodes. 
\param mvc stores the messages fom variable nodes to check nodes.
 */
void initialize_massages(mat &mcv,mat& mvc,const GF2mat &H);

/*! \relates nodes

 *Do a parallel iteration for classical BP.
\param c is the number of check nodes.
\param v is the number of varaible nodes.
\param output_e is the result error vector.
 */
void p_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const GF2mat& syndrome,double p,int c, int v,  GF2mat& output_e, int debug);

void s_C_update(const nodes checks[],const nodes errors[],mat &mcv,mat &mvc,const GF2mat& syndrome,double p,int c, int v,  GF2mat& output_e,int debug);

/*! \relates nodes

 *Do a serial iteration for classical BP.
\param c is the number of check nodes.
\param v is the number of varaible nodes.
\param output_e is the result error vector.
 */
void s_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const GF2mat& syndrome,double p,int c, int v,  GF2mat& output_e,int debug);
  
/*! \relates nodes
 *  Update the message c_i to v_j, that is mcv(i,j).

\param s is syndrome(i).
 */
 //void update_ci_to_vj(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,int i,int j,bin s);
 //void update_ci_to_vj_MS(const nodes checks[],const nodes errors[],mat& mcv,mat& mvc,int i,int j,bin s);

/*! \relates nodes
 * Update the message v_j to c_i, that is mvc(i,j).

\param ipr =(1-p)/p.
\param alpha is for larger step size.
 */
 //void update_vj_to_ci(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,int j,int i,double ipr,double alpha=1);

inline void update_ci_to_vj(const nodes checks[],const nodes errors[],mat& mcv,mat& mvc,const int& i,const int& j,const bin& s)
{
   int c_degree=checks[i].degree;
   
   // if (c_degree==2) {cout<<"wt=2 check here"<<endl;}  
   double temp=1.0;
   for (int kk=0;kk<c_degree;kk++)
     {
      int vk=(checks[i].neighbors)(kk);
      if (vk!=j)
	{
		double mvc_value = mvc(i, vk);
	temp *= (mvc_value - 1) / (mvc_value + 1);
	}
     }
  
     mcv.set(i,j,s==0? (1+temp)/(1-temp):(1-temp)/(1+temp));
     if (isinf(mcv(i,j)))
       {
	 //cout<<"inf"<<endl;
	 mcv.set(i,j,1);
       }    
}

inline void update_ci_to_vj_MS(const nodes checks[],const nodes errors[],mat& mcv,mat& mvc,const int& i,const int& j,const bin& s)
{
   int c_degree=checks[i].degree;
   
   // if (c_degree==2) {cout<<"wt=2 check here"<<endl;}
   int sign_of_log_mvc=1;
   double min_log_mvc=std::numeric_limits<double>::infinity();
   
   for (int kk=0;kk<c_degree;kk++)
     {
      int vk=(checks[i].neighbors)(kk);
      double log_mvc=0;
      
      if (vk!=j)
	{
	  log_mvc=log(mvc(i,vk));
	  if (abs(log_mvc)<min_log_mvc) {min_log_mvc=log_mvc;}
	  sign_of_log_mvc=sign_of_log_mvc * ((log_mvc>0)-(log_mvc<0)); // notice log_mvc is impossible to be 0, so the expression is safe.	    
	 
	}
     }

   int sint=s;
   double log_mcv=pow(-1,sint)*min_log_mvc*sign_of_log_mvc;
   mcv.set(i,j,exp(log_mcv));
     // if (c_degree==2) {cout<<"mcv is "<<mcv(i,j)<<endl;}
   if (isinf(mcv(i,j)))
     {
       //cout<<"inf"<<endl;
       mcv.set(i,j,1);
     }    
}

//update vj to ci massage:
inline void update_vj_to_ci(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const int& j,const int& i,const double& ipr,double alpha=1)
{ 
  int v_degree=errors[j].degree;
  mvc.set(i,j,ipr);

  for (int kk=0;kk<v_degree;kk++)
    {
      int ck=(errors[j].neighbors)(kk);
      if (ck!=i)
		{
	  mvc.set(i,j,mvc(i,j)*mcv(ck,j));
		}
    }   
}

inline void update_vj_to_ci_larger_step(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const int& j,const int& i,const double& ipr,double alpha=1)
{ 
  int v_degree=errors[j].degree;
  mvc.set(i,j,ipr);

  for (int kk=0;kk<v_degree;kk++)
    {
      int ck=(errors[j].neighbors)(kk);
      if (ck!=i)
	{
	  mvc.set(i,j,mvc(i,j)*pow(mcv(ck,j),1.0/alpha));
	}
      else
	{
	  mvc.set(i,j,mvc(i,j)*pow(mcv(ck,j),1.0-1.0/alpha));
	}   
     }   
}
/*! \relates nodes
 *Generate a error vector, decode it by classical BP, output 1 if secceeds, output 0 if fails.

 *\param v is the number of variable nodes.
 *\param c is the number of check nodes.
\param H is the parity check matrix.
 *\param num_iter stores the total iterations for successful decoding. num_iter/num_of_suc_decoding=averge iterations for seccessful decoding.
 * \param lmax is the maximum number of iterations.
 * \param er stores the number of bit errors after decoding, er/num_of_cws= bit error rate afer decoding.
 *\param debug debug=1 for parallel schedule, debug=0 for serial schedule.
 */
int  cla_decode(int v,int c,const GF2mat &H,const nodes checks[],const nodes errors[],double& num_iter, int lmax,int & er,double p,int debug); 

/*! \relates nodes
 * Output the distance between output_e and real_e.

\param n is the length of the vectors.
 */
int distance(const GF2mat &output_e, const GF2mat &real_e,int n);

/*! \relates nodes
 *  Read a parity check matrix from file_name, return it. The first element in the stored file is 1.

\param  n stores the columns of the martrix
\param  r stores the rows of the matrix.
 */
GF2mat read_matrix (int& n,int &r, string & file_name);

/*! \relates nodes
 *  Read a parity check matrix from file_name, return it. The first element in the stored file is 0.

\param  n stores the columns of the martrix
\param  r stores the rows of the matrix.
 */
GF2mat read_matrix_from_0 (int& n,int &r, string & file_name);
  

/*! \relates nodes
 *Write a parity check matrix H in file_name.
 */
void write_matrix(string file_name, GF2mat &H);

/*! \relates nodes
 * Merge 2 matrices horizentolly, return a matrix=(left,right).
 */
GF2mat merge_mat_hori(const GF2mat &left,const GF2mat &right);

/*! \relates nodes
 * Classical BSC channel with error rate p, input codeword cw, return  weight of the error.
 */
int cla_error_channel(GF2mat &cw, double p);

/*! \relates nodes
 *Quantum BSC channel, e(i)=1 with probalility p(i). Return  weight of the error.
 */
int error_channel(GF2mat &cw, const vec &p);
int depolarizing(GF2mat &cw, double p);
int depolarizing2(GF2mat &ex, GF2mat &ez, double p);
/*! \relates nodes
 *A quantum error channel which generates errors with weight wt.
 */
void error_channel2(GF2mat &error, int wt);

/*! \relates nodes
 *Return the weight of a column vector.
 */
int s_weight(const GF2mat &s);

/*! \relates nodes
 *Return the weight of a row vector.
 */
int weight(const GF2mat &cw);

/*! \relates nodes
 *Get the probability distribution pv, pv(i) is randomely choosen betweem pmin and pmax.
 */
void pro_dist(double pmin,double pmax, vec& pv);
void pro_dist2(double pmin,double pmax, vec& pv);

/*! \relates nodes
 *Generate a error vector, decode it by quantum BP, return true if secceeds, return false if fails.

 *\param H is the parity check matrix.
\param G is the matrix for corresponding logical operators.
 *\param pv is the error rate distribution for error channel.
\param pv_dec is the error rate distribution for decoding.
 *\param num_iter is for calculating average number of iterations.
 *\param lmax is maximum number of iterations.
 *\param wt if wt=0, use error_channel, if wt>=1, use error_channel2.
 *\param max_fail stores the number of failures such that reaches the maximum iterations.
 \param syn_fail stores the number of failures such that output_e+real_e=a logical operator.
 \param debug is for bitwise options, see readme.md
 \param LR is the likelyhood-ratios.
 *\param rankH is the rank of H.
 *\param OSD_suc stores the number of successful decoding by OSD.
 *\param alpha is for larger step size.
 * \param lambda has not been used yet, just set it be 1 is ok. this param is deleted.
 *\param other_decoder if it equals to 1 , use serial C schedule
 */
bool  quan_decode(GF2mat &H, GF2mat &G,const nodes checks[],const nodes errors[],const vec &pv,const vec&pv_dec,const double& pmin,const double& pmax,double& num_iter, int lmax,int &wt,int& max_fail, int&syn_fail,  int debug, vec &LR,int rankH,int &OSD_suc,int& num_no_error,double alpha=1,int other_decoder=0);

/*! \relates nodes
 *Ordered statistical decoder, input LR, get output_e.
 */
void full_OSD(vec& LR,const GF2mat& H,const GF2mat& syndrome,GF2mat& output_e,int OSD_order,int &  suc_order);
void OSD(vec& LR,const GF2mat& H,const GF2mat& syndrome,GF2mat& output_e,int&suc_order);
void GPT_dfs_For_OSD(int index, int n, int OSDorder, int non_zero_count, int& suc_order,GF2mat& new_e_T, GF2mat& e_T,const GF2mat& e_S, int& min_wt, const GF2mat& HS_inv, const GF2mat&H_T,const GF2mat& syndrome2) ;
bool  OSD_decode(GF2mat &H, GF2mat &G,int &wt, vec& pv, int debug,vec &LR,int rankH,int OSD_order);
/*! \relates nodes
 *Generate a column permutation matrix by permutation vector perm.
 */
GF2mat col_permutation_matrix(ivec & perm);
GF2mat col_permutation_matrix_s(ivec & perm);

/*! \relates nodes
 *Do a quantum parallel iteration.
 */
void quan_p_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const GF2mat& syndrome,const vec &pv,int c, int v,  GF2mat& output_e, vec &LR,int debug,double alpha=1);

/*! \relates nodes
 *Do a quantum serial iteration.
 */
void quan_s_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const GF2mat& syndrome,const vec &pv,int c, int v,  GF2mat& output_e, vec &LR,int debug,double alpha=1);

void quan_s_C_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const GF2mat& syndrome,const vec &pv,int c, int v,  GF2mat& output_e, vec &LR,int debug,double alpha=1);

/*! \relates nodes
 *Show the patterns of error for checkerboard codes.
 */
void err_pos2(const GF2mat &error);

/*! \relates nodes
 * Return the rank of a GF2mat matrix
 */
int GF2mat_rank(const GF2mat& H_s);

/*! \relates nodes
 *Convert dense GF2mat matrix to sparse one.
 */
void dense_to_sparse(GF2mat &G,GF2mat_sparse& G_s);

/*! \relates nodes
 *cout i if error(i)==1.
 */
void err_pos1(const GF2mat &error);

/*! \relates nodes
 *Return the logical operators matrix for H.
 */
GF2mat get_gen(const GF2mat &H);
GF2mat get_gen(const GF2mat &Hx, const GF2mat &Hz);

bmat GF2mat_to_bmat(const GF2mat& H);

GF2mat merge_mat_vert(const GF2mat &up,const GF2mat &bottom);
void get_logical(const GF2mat& H, const GF2mat& D,const GF2mat& H_star, GF2mat& L);
bool  decoder_T(GF2mat &H, GF2mat &HT, GF2mat &G, GF2mat &new_H,const nodes checks[],const nodes errors[],const vec &pv,const vec&pv_dec,double& num_iter, int lmax,int &wt,int debug, vec &LR,int rankH);
#endif